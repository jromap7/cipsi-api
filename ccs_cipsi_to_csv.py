__author__ = "Jordi Roma & Maria-Jose Falaguera"
__date__ = "11 Mar 2022"

"""
fill_cipsi.py: Fill CIPSI database with processed file Get patent vs molecule associations from the following sources (see Southan 2015, Drug Dicov. Today):
 - SureChEMBL: https://ftp.ebi.ac.uk/pub/databases/chembl/SureChEMBL/data/map/
 - GooglePatents (from PubChem): https://ftp.ncbi.nlm.nih.gov/pubchem/Other/GooglePatents/
 - IBM (from PubChem): https://ftp.ncbi.nlm.nih.gov/pubchem/Other/IBM/
""" 
import json
import numpy as np
import itertools
import gzip
import rdkit
import psycopg2
import psycopg2.extras

from timeit import default_timer as timer
from multiprocessing import Pool, cpu_count
from cipsi import getPatents
from datetime import datetime
from rdkit import Chem
from rdkit import DataStructs # calculate similarity
from rdkit.Chem import rdFMCS # calculate MCS
from rdkit.Chem import Descriptors # to calculate MW
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds.MurckoScaffold import MakeScaffoldGeneric
from rdkit import RDLogger
# Skips console warning from RDKit
RDLogger.DisableLog('rdApp.*')

# Read config parameters
config_fname = './config/parameters.json'

def read_json(json_fname):
    with open(json_fname) as fh:
        result = json.load(fh)
    return(result)

params = read_json(config_fname)
version = 'v'+str(params['year_version'])
TABLES_VERSION_FOLDER = params['output_tables_folder']+'/'+version

# Database connection parameters
DB_HOST = 'blackadder'
DB_USERNAME = 'jroma'
DB_PASSWORD = params['db_password']
DB_NAME = 'cipsi_db'

def openConnection():
    conn = psycopg2.connect(
        database=DB_NAME,
        user=DB_USERNAME,
        host=DB_HOST,
        password=DB_PASSWORD,
        keepalives=1,
        keepalives_idle=30,
        keepalives_interval=10,
        keepalives_count=5
        #options='-c search_path={}'.format(schema),
    )
    return(conn)


def write_mols_to_csv(values):
    fname = TABLES_VERSION_FOLDER+'/'+'ccs.csv.gz'
    with gzip.open(fname,'at') as fh:
        for value in values:
            fh.write(value+'\n')

def create_table_ccs(values):

     #sql = "INSERT INTO {table} (mol_inchikey, mol_name) VALUES (%s, %s) \
    sql = "INSERT INTO {table} (mol_inchikey, ccs_mol) \
            SELECT mol_inchikey, mol_mol as ccs_mol \
            FROM molecule \
            WHERE mol_inchikey IN ('{WHERE}')\
        ON CONFLICT (mol_inchikey) DO NOTHING".format(table='ccs', WHERE= "','".join(list(values)))

    conn = openConnection()
    cur = conn.cursor()
    cur.execute(sql)
    conn.commit()
    cur.close()
    conn.close()

def getMoleculesForPatent(patent):
    # get molecules associated to patent in SureChEMBL

    SELECT = ['m.mol_inchikey','mol_to_smiles(m.mol_mol)','mol_amw(m.mol_mol)']
    FROM = ['patent p']
    INNER_JOIN = ['pat2mol ptm ON p.patent_id = ptm.patent_id',
                    'molecule m ON m.mol_inchikey = ptm.mol_inchikey']
    WHERE = ["p.patent_id = '{}'".format(patent)]

    query = 'SELECT {} FROM {} INNER JOIN {} WHERE {};'.format(','.join(SELECT), ','.join(FROM), 
                                                            ' LEFT JOIN '.join(INNER_JOIN), ' AND '.join(WHERE)) 
        
    conn = openConnection()
    with conn.cursor() as cur:
        cur.execute(query)
        
        molecules = {}
        for i in cur:
            mw = i[2]
            if mw > float(600): continue
            inchikey = i[0]
            smiles = i[1]
            mol = Chem.MolFromSmiles(smiles)
            fp = AllChem.GetMorganFingerprint(mol, 2)
            
            
            molecules[inchikey] = {'mol': mol, 'fp': fp, 'smiles': smiles, 'mw': mw}    

    conn.close()
    return(molecules)

def calculate_ccs(molecules):

    # calculate molecular pairwise similarity
    similarities = {}
    for mol1, mol2 in itertools.combinations(list(molecules), 2):
        sim = DataStructs.DiceSimilarity(molecules[mol1]['fp'], molecules[mol2]['fp'])
        similarities[tuple(sorted([mol1, mol2]))] = round(sim, 3)

    if len(similarities) == 0:
        return
    # calculate MCS for molecular pairs with sim >= 0.4
    MCS = {}
    for mol1, mol2 in similarities:
        if similarities[(mol1, mol2)] >= 0.4:
            mcs = rdFMCS.FindMCS([molecules[mol1]['mol'],
                                molecules[mol2]['mol']],
                                ringMatchesRingOnly=True, completeRingsOnly=True, timeout=1).smartsString
            if mcs == '': continue # no MCS found
                    
            try:
                _ = MCS[mcs]
            except KeyError:
                MCS[mcs] = {'pairs': set()}
            MCS[mcs]['pairs'].update([(mol1, mol2)])
    
    # asses MCS coverage of molecules
    for mcs in MCS:
        m = Chem.MolFromSmarts(mcs)
        MCS[mcs]['mol'] = m
        MCS[mcs]['mols'] = set()
        
        # save mols in pairs as covered by MCS
        for mol1, mol2 in MCS[mcs]['pairs']:
            MCS[mcs]['mols'].update([mol1, mol2])
        
        # save mols that although do not generated this MCS, do contain the MCS as a substr.
        for mol in set(molecules).difference(MCS[mcs]['mols']):
            if molecules[mol]['mol'].HasSubstructMatch(m):
                MCS[mcs]['mols'].update([mol])
            
        MCS[mcs]['n_mols'] = len(MCS[mcs]['mols'])
        MCS[mcs]['cov'] = round(float(MCS[mcs]['n_mols'])/float(len(molecules)), 3)

    # asses MCS internal homogeneity
    for mcs in MCS:
        hom = []
        for mol1, mol2 in itertools.combinations(MCS[mcs]['mols'], 2):
            sim = similarities[tuple(sorted([mol1, mol2]))]
            hom.append(sim)
        MCS[mcs]['hom'] = np.median(hom)

    # assess MCS substructure inclusion and congeneric degree
    for mcs in MCS:
        MCS[mcs]['congeneric'] = set()
        MCS[mcs]['included'] = set()
    for mcs1, mcs2 in itertools.combinations(MCS, 2):
        mol1 = MCS[mcs1]['mol']
        mol2 = MCS[mcs2]['mol']

        if mol1.HasSubstructMatch(mol2, useQueryQueryMatches=True):
            MCS[mcs1]['included'].update([mcs2]) # MCS1 contains MCS2 as a substructure (MCS1 includes MCS2)
            MCS[mcs2]['congeneric'].update([mcs1]) # mol2 contains mol1 as a substructure (mol1 covers mol2)
        if mol2.HasSubstructMatch(mol1, useQueryQueryMatches=True):
            MCS[mcs2]['included'].update([mcs1])
            MCS[mcs1]['congeneric'].update([mcs2])  

    n_mcs = len(MCS)
    for mcs in MCS:
        MCS[mcs]['inc'] = round(float(len(MCS[mcs]['included']))/float(n_mcs), 3)

    # score MCS
    for mcs in MCS:
        score = round(MCS[mcs]['cov']*MCS[mcs]['hom']*MCS[mcs]['inc'], 4)
        MCS[mcs]['score'] = score
    scores = [ MCS[mcs]['score'] for mcs in MCS ]

    # select top MCS
    if len(scores) > 1:
        cutoff = round(np.quantile(scores, 0.99), 4)
    else:
        return()
    #print(cutoff)
    topMCS = { mcs: data for mcs, data in MCS.items() if data['score'] >= cutoff }
    #print(len(topMCS))

    mols = []
    for mcs in sorted(topMCS, key=lambda x: topMCS[x]['score'], reverse=True):
        mols.append(topMCS[mcs]['mol'])

    # remove redundant MCS (collapse onto better scored ones)
    tmp = set()
    for mcs1 in sorted(topMCS, key=lambda x: topMCS[x]['score'], reverse=True):
        new = True
        for mcs2 in tmp:
            if mcs1 in topMCS[mcs2]['congeneric']: # mcs1 contains mcs2 as a substructure
                new = False
        if new:
            tmp.update([mcs1])
            
    topMCS = { mcs: topMCS[mcs] for mcs in tmp }
    len(topMCS)

    # view CCS
    mols = []
    for mcs in sorted(topMCS, key=lambda x: topMCS[x]['score'], reverse=True):
        mols.append(topMCS[mcs]['mol'])
    
    # recover mols not having the MCS as substructure but that are highly similar to MCS covered molecules
    covered = set([ mol for mcs in topMCS for mol in topMCS[mcs]['mols'] ])
    to_recover = set(molecules).difference(covered)
    recovered = set([])
    for mol1 in covered:
        for mol2 in to_recover:
            sim = similarities[tuple(sorted([mol1, mol2]))]
            if sim >= 0.8:
                recovered.update([mol2])
        to_recover = to_recover.difference(recovered)
    covered = covered.union(recovered)

    # calculate patent global molecular homogeneity
    homogeneity = []
    for mol1, mol2 in itertools.combinations(covered, 2):
        sim = similarities[tuple(sorted([mol1, mol2]))]
        homogeneity.append(sim)
    
    # view SureChEMBLccs molecules
    mols = []
    for mol in covered:
        mols.append(molecules[mol]['mol'])
    
    values = set()  # TODO add patent_id, MCS smarts to values and mol smiles to values
    for mol in mols:
        values.update([Chem.inchi.InchiToInchiKey(Chem.inchi.MolToInchi(mol))])

    write_mols_to_csv(values)

def get_ccs(patent):
    molecules = getMoleculesForPatent(patent)
    if len(molecules) > 1:
        print('Len molecules:', len(molecules))
        values = calculate_ccs(molecules)
        return(values)
    else:
        return()


def worker(patent):
    #logfile = open('logfile_css_cipsi.txt', 'at', buffering=1)
    get_ccs(patent)
            #create_table_ccs(values)

def main():
    start = timer()

    print(f'starting computations on {cpu_count()} cores')

    patents = set()
    fname = TABLES_VERSION_FOLDER+'/'+'patent_test.csv.gz'
    print('Reading patents CSV...')
    with gzip.open(fname,'rt') as fh:
        for line in fh:
            line = line.strip('\n').split('\t')
            patents.update([line[0]])
    print('Done.')

    with Pool() as pool:
        pool.map(worker, patents)
    end = timer()
    print(f'elapsed time: {end - start}')


if __name__ == '__main__':
    main()