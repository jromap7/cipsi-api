#!/usr/bin/env python

__author__ = "María Martínez, Maria-José Falaguera"
__date__ = "26 July 2022"

"""
ccs.py: Script for calculating the CCS for patents
"""

import itertools
from multiprocessing import Pool
import os

import numpy as np
import psycopg2
import psycopg2.extras
from rdkit import Chem
from rdkit import DataStructs # calculate similarity
from rdkit import RDLogger
from rdkit.Chem import rdFMCS # calculate MCS

# Skips console warning from RDKit
RDLogger.DisableLog('rdApp.*')


def calculate_ccs(patent):
    """
    Calculate CCS for the molecules of a patent.

    Args:
        patent (str):   ID of the patent from which the CCS will be extracted.

    Returns:
        CSV file: 
    """

    print('In calculate_ccs')
    print(patent)

    # Define the database connection
    print('Connecting to database')

    DB_HOST = 'erebor.prib.upf.edu'
    DB_USERNAME = 'jroma'
    DB_PASSWORD = 'c1ps1'
    DB_NAME = 'cipsi_db'

    conn = get_db_connection(DB_HOST, DB_USERNAME, DB_PASSWORD, DB_NAME)

    # Get molecules associated to the patent
    print('Getting molecules associated to the patent')

    SELECT = ['distinct(m.mol_inchikey)','mol_send(m.mol_mol)', 'bfp_to_binary_text(m.mol_ecfp4)']
    FROM = ['molecule m']
    INNER_JOIN = ['pat2mol pm ON pm.mol_inchikey = m.mol_inchikey']
    WHERE = ["pm.patent_id = '{}' and m.mol_amw < 600.0".format(patent)]

    query = 'SELECT {} FROM {} INNER JOIN {} WHERE {};'.format(','.join(SELECT), ','.join(FROM),
                                                        ' LEFT JOIN '.join(INNER_JOIN), ' AND '.join(WHERE))

    with conn.cursor() as cur:
        cur.execute(query)
        molecules = {}
        for i in cur:
            inchikey = i[0]
            mol = Chem.Mol(i[1].tobytes())
            fp_expl = DataStructs.CreateFromBinaryText(i[2].tobytes())
            fp_np = np.frombuffer(fp_expl.ToBitString().encode(), 'u1') - ord('0')
            fp = DataStructs.cDataStructs.UIntSparseIntVect(len(fp_np))
            for ix, value in enumerate(fp_np):
                fp[ix] = int(value)
            molecules[inchikey] = {'mol': mol, 'fp': fp}
        cur.close()
    conn.close()

    # Connection to save the new registers
    DB_HOST = 'erebor.prib.upf.edu'
    DB_USERNAME = 'maria'
    DB_PASSWORD = 'c1ps1d3v'
    DB_NAME = 'cipsi_dev'

    # End function if less than 2 molecules were found
    if len(molecules) <= 1:
        print('Less than 2 molecules were found')
        conn = get_db_connection(DB_HOST, DB_USERNAME, DB_PASSWORD, DB_NAME)
        query = '''insert into pat2ccs2mol (patent_id)
                    values (\'''' + patent + '''\')                                 
        '''
        with conn.cursor() as cur:
            try:
                cur.execute(query)
                conn.commit()
            except psycopg2.Error as ex:
                print(query + ":\n" + str(ex),"ERROR")
                query = '''insert into tmp_ccs_errors (patent_id, error)
                            values(\'''' + patent + '''\', 
                            \'''' + query + '''\': \'''' + str(ex) + '''\')
                '''
                cur.execute(query)
                conn.commit()
            finally:
                cur.close()
                conn.close()
                return()

    # 1. Extraction of MCSs
    ## Calculate molecular pairwise similarity
    print('Calculating molecules pairwise similarity')

    similarities = {}
    for mol1, mol2 in itertools.combinations(list(molecules), 2):
        sim = DataStructs.DiceSimilarity(molecules[mol1]['fp'], molecules[mol2]['fp'])
        similarities[tuple(sorted([mol1, mol2]))] = round(sim, 3)

    ## Calculate MCS for molecular pairs with sim >= 0.4
    print('Calculating MCSs')

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

    # 2. Selection of candidate MCSs
    print('Calculating MCS score')

    ## Assess MCS coverage of molecules
    for mcs in MCS:
        m = Chem.MolFromSmarts(mcs)
        MCS[mcs]['mol'] = m
        MCS[mcs]['mols'] = set()

        # save mols in pairs as covered by MCS
        for mol1, mol2 in MCS[mcs]['pairs']:
            MCS[mcs]['mols'].update([mol1, mol2])

        # save mols that although do not generate this MCS, do contain the MCS as a substr.
        for mol in set(molecules).difference(MCS[mcs]['mols']):
            if molecules[mol]['mol'].HasSubstructMatch(m):
                MCS[mcs]['mols'].update([mol])

        MCS[mcs]['n_mols'] = len(MCS[mcs]['mols'])
        MCS[mcs]['cov'] = round(float(MCS[mcs]['n_mols'])/float(len(molecules)), 3)

    ## Assess MCS internal homogeneity
    for mcs in MCS:
        hom = []
        for mol1, mol2 in itertools.combinations(MCS[mcs]['mols'], 2):
            sim = similarities[tuple(sorted([mol1, mol2]))]
            hom.append(sim)
        MCS[mcs]['hom'] = np.median(hom)

    ## Assess MCS substructure inclusion and congeneric degree
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

    ## Score MCS
    for mcs in MCS:
        score = round(MCS[mcs]['cov']*MCS[mcs]['hom']*MCS[mcs]['inc'], 4)
        MCS[mcs]['score'] = score
    scores = [ MCS[mcs]['score'] for mcs in MCS ]

    ## Select top MCS
    print('Selecting top MCSs')

    if len(scores) == 0:
        print('No MCSs were found')
        conn = get_db_connection(DB_HOST, DB_USERNAME, DB_PASSWORD, DB_NAME)
        query = '''insert into pat2ccs2mol (patent_id)
                    values (\'''' + patent + '''\')                                 
        '''
        with conn.cursor() as cur:
            try:
                cur.execute(query)
                conn.commit()
            except psycopg2.Error as ex:
                print(query + ":\n" + str(ex),"ERROR")
                query = '''insert into tmp_ccs_errors (patent_id, error)
                            values(\'''' + patent + '''\', 
                            \'''' + query + '''\': \'''' + str(ex) + '''\')
                '''
                cur.execute(query)
                conn.commit()
            finally:
                cur.close()
                conn.close()
                return()

    cutoff = round(np.quantile(scores, 0.99), 4)
    topMCS = { mcs: data for mcs, data in MCS.items() if data['score'] >= cutoff }

    ## Remove redundant MCS (collapse onto better scored ones)
    print('Removing redundant MCSs')

    tmp = set()
    for mcs1 in sorted(topMCS, key=lambda x: topMCS[x]['score'], reverse=True):
        new = True
        for mcs2 in tmp:
            if mcs1 in topMCS[mcs2]['congeneric']: # mcs1 contains mcs2 as a substructure
                new = False
        if new:
            tmp.update([mcs1])

    topMCS = { mcs: topMCS[mcs] for mcs in tmp }

    ## View CCS
    mols = []
    for mcs in sorted(topMCS, key=lambda x: topMCS[x]['score'], reverse=True):
        mols.append(topMCS[mcs]['mol'])

    # 3. Recover mols not having the MCS as substructure but that are highly similar to MCS covered molecules
    print('Recovering molecules')

    recovered_MCS = {}
    covered = set([ mol for mcs in topMCS for mol in topMCS[mcs]['mols'] ])
    to_recover = set(molecules).difference(covered)
    recovered = set([])
    for mol1 in covered:
        for mol2 in to_recover:
            sim = similarities[tuple(sorted([mol1, mol2]))]
            if sim >= 0.8:
                recovered.update([mol2])
                for mcs in topMCS:
                    for mol in topMCS[mcs]['mols']:
                        if mol == mol1:
                            try:
                                _ = recovered_MCS[mcs]
                            except KeyError:
                                recovered_MCS[mcs] = set()
                            recovered_MCS[mcs].update([mol2])
        to_recover = to_recover.difference(recovered)
    covered = covered.union(recovered)

    # 4. Merge topMCS and recovered_MCS into a dictionary
    print('Calculating final MCSs')

    finalMCS = { mcs: topMCS[mcs]['mols'] for mcs in topMCS }

    for mcs in finalMCS:
        try:
           _ = recovered_MCS[mcs]
           for mol in recovered_MCS[mcs]:
               finalMCS[mcs].update([mol])
        except KeyError:
           pass

    print('Final MCS for patent ' + patent)
    print(finalMCS)

    # 5. Write results into the database
    conn = get_db_connection(DB_HOST, DB_USERNAME, DB_PASSWORD, DB_NAME)
    for mcs in finalMCS:
        values = ''
        print('Generating query')

        for mol in finalMCS[mcs]:
            value = '(\'{patent}\', \'{mcs}\', \'{mol}\')'.format(patent=patent, mcs=mcs, mol=mol)
            values = values + ', ' + value
        values = values[2:]

        query = 'insert into pat2ccs2mol values ' + values

        print('Inserting data into ddbb')
        with conn.cursor() as cur:
                try:
                    cur.execute(query)
                    conn.commit()
                except psycopg2.Error as ex:
                    print(query + ":\n" + str(ex),"ERROR")
                    query = '''insert into tmp_ccs_errors (patent_id, error)
                                values(\'''' + patent + '''\', 
                                \'''' + query + '''\': \'''' + str(ex) + '''\')
                    '''
                    cur.execute(query)
                    conn.commit()

    cur.close()
    conn.close()
    return()


def run_calculate_ccs(patentList, nCpus=os.cpu_count()-4):
    """
    Execute the function calculate_ccs for all the patents in a list.

    Args:
        patentList (str[]):     List with patent IDs.
        nCpus (int):            Number of worker processes to use in the Pool.
    """

    print('In run_calculate_ccs')
    print('Running on', os.cpu_count()-4)
    print(patentList)

    if nCpus == 1:
        for patent in patentList:
            calculate_ccs(patent)
    else:
        with Pool(processes=nCpus) as pl:
            pl.map(calculate_ccs, patentList)


def get_db_connection(host, username, password, dbname):
    print('Connecting to database')

    conn = psycopg2.connect(
        database=dbname,
        user=username,
        host=host,
        password=password,
        keepalives=1,
        keepalives_idle=30,
        keepalives_interval=10,
        )

    return(conn)


if __name__ == "__main__":

    # Get patent list from database
    DB_HOST = 'erebor.prib.upf.edu'
    DB_USERNAME = 'jroma'
    DB_PASSWORD = 'c1ps1'
    DB_NAME = 'cipsi_db'

    conn = get_db_connection(DB_HOST, DB_USERNAME, DB_PASSWORD, DB_NAME)

    SELECT = ['patent_id']
    FROM = ['pat2mol']
    GROUP_BY = ['patent_id']

    query = 'SELECT {} FROM {} GROUP BY {} HAVING COUNT(*) > 1 ORDER BY COUNT(*), patent_id ASC;'.format(','.join(SELECT),
    ','.join(FROM), ','.join(GROUP_BY))

    with conn.cursor() as cur:
        cur.execute(query)
        patentList = []
        for i in cur:
            patent = i[0]
            patentList.append(patent)
        cur.close()
    conn.close()

    # See if the table exists and if not create new one
    DB_HOST = 'erebor.prib.upf.edu'
    DB_USERNAME = 'maria'
    DB_PASSWORD = 'c1ps1d3v'
    DB_NAME = 'cipsi_dev'

    conn = get_db_connection(DB_HOST, DB_USERNAME, DB_PASSWORD, DB_NAME)

    tableName = 'pat2ccs2mol'

    with conn.cursor() as cur:
        cur.execute('select exists(select * from information_schema.tables where table_name=%s)', (tableName,))
        exists = cur.fetchone()[0]
        cur.close()

    with conn.cursor() as cur:
        if exists:
            query = 'delete from ' + tableName
        else:
            query = '''create table ''' + tableName + ''' (
                            patent_id VARCHAR(50),
                            ccs_smarts TEXT,
                            mol_inchikey VARCHAR(50)
                            )
            '''
        cur.execute(query)
        conn.commit()
        cur.close()
    conn.close()

    run_calculate_ccs(patentList)
