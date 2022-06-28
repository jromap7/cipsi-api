#!/usr/bin/env python

"""
tcrd.py: Query Target Central Resource Database (TCRD).
"""

__author__ = "Maria J. Falaguera"
__date__ = "31 Mar 2022"


import json

import psycopg2
import psycopg2.extras

def read_json(json_fname):
    with open(json_fname) as fh:
        result = json.load(fh)
    return(result)

# Read config parameters
config_fname = './config/parameters.json'
params = read_json(config_fname)

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
        #options='-c search_path={}'.format(schema),
    )
    return(conn)

def getPatents(nationality=None):
    """
    Get patents.

    Args:
        nationality (str/list/set/dict):    nationality short 'XX' e.g. 'US'

    Returns:
        set:   {patent_id_1, 2, 3 ...}
    """

    SELECT = ['patent.patent_id']
    FROM = ['patent']
    WHERE = ['1 = 1']

    if nationality is not None:
        if isinstance(nationality, str):
            WHERE.append("patent.nationality = ('{}')".format("','".join(list(nationality))))
        elif isinstance(nationality, list) or isinstance(nationality, set) or isinstance(nationality, dict):
            WHERE.append("patent.nationality IN ('{}')".format("','".join(list(nationality))))

    sql = 'SELECT DISTINCT {} FROM {} WHERE {};'.format(','.join(SELECT), ','.join(FROM), ' AND '.join(WHERE))
    conn = openConnection()
    cur = conn.cursor()
    cur.execute(sql)
    results = cur.fetchall()
    results = set([result[0] for result in results])

    cur.close()
    conn.close()

    return(results)

def countPatents(nationality=None):
    """
    Get patents count.

    Args:
        nationality (str/list/set/dict):    nationality short 'XX' e.g. 'US'

    Returns:
        int: number of patents
    """
    return(len(getPatents(nationality=nationality)))

def getInfoForPatent(patent=None):
    """
    Get info for patent.

    Args:
        patent (str/list/set/dict):    patent id (wO '-')

    Returns:
        dict: {'Source_1':{patent_1, patent_2,...}, 'Source_2':{...}, ...}
    """
    SELECT = ['patent.patent_id', 'patent.nationality']
    FROM = ['patent']
    WHERE = ['1 = 1']

    if patent is not None:
        if isinstance(patent, str):
            WHERE.append("patent.patent_id = ('{}')".format("','".join(list(patent))))
        elif isinstance(patent, list) or isinstance(patent, set) or isinstance(patent, dict):
            WHERE.append("patent.patent_id IN ('{}')".format("','".join(list(patent)))) 

    sql = 'SELECT DISTINCT {} FROM {} WHERE {};'.format(','.join(SELECT), ','.join(FROM),' AND '.join(WHERE))
    results = {}
    conn = openConnection()
    cur = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    cur.execute(sql)

    for row in cur:
        for r, v in row.items():
            if r != 'patent_id':
                aux = {r:v}
                
        results[row['patent_id']] = aux

    cur.close()
    conn.close()

    return(results)

def getPatentForSource(source=None):
    """
    Get patents for each source.

    Args:
        source (str/list/set/dict):    source name

    Returns:
        dict: {'Source_1':{patent_1, patent_2,...}, 'Source_2':{...}, ...}
    """
    SELECT = ['p2s.patent_id','s.source_name']
    FROM = ['source s']
    INNER_JOIN = ['pat2source p2s ON p2s.source_name = s.source_name']
    WHERE = ['1 = 1']

    if source is not None:
        if isinstance(source, str):
            WHERE.append("p2s.source_name = ('{}')".format("','".join(list(source))))
        elif isinstance(source, list) or isinstance(source, set) or isinstance(source, dict):
            WHERE.append("p2s.source_name IN ('{}')".format("','".join(list(source)))) 

    sql = 'SELECT {} FROM {} INNER JOIN {} WHERE {};'.format(','.join(SELECT), ','.join(FROM),' INNER JOIN '.join(INNER_JOIN) ,' AND '.join(WHERE))             


    results = {}
    conn = openConnection()
    cur = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    cur.execute(sql)
    for row in cur:
        try:
            _ = results[row['source_name']]
        except:
            results[row['source_name']] = set()
        results[row['source_name']].update([(row['patent_id'])])
        

    cur.close()
    conn.close()

    return(results)

def countPatentForSource(source=None):#none
    """
    Get patent count associated with source.

    Args:
        source (str/list/set/dict):    source name

    Returns:
        dict: {'Source_1': number of patents for source_1, 'Source_2': number of patents for source_2, ...}
    """
    return({k:len(v) for k,v in getPatentForSource(source=source).items()})

def getSourceForPatent(patent):
    """
    Get sources associated with patent.

    Args:
        patent (str/list/set/dict):    patent id (wO '-')

    Returns:
        dict: {patent_1:{source_1,source_2, ...}, patent_2:{..}, ...}
    """
    SELECT = ['p2s.source_name','p.patent_id']
    FROM = ['patent p']
    LEFT_JOIN = ['pat2source p2s ON p2s.patent_id = p.patent_id']
    WHERE = ['1 = 1']

    if patent is not None:
        if isinstance(patent, str):
            WHERE.append("p2s.patent_id = ('{}')".format(patent))
        elif isinstance(patent, list) or isinstance(patent, set) or isinstance(patent, dict):
            WHERE.append("p2s.patent_id IN ('{}')".format("','".join(list(patent)))) 

    sql = 'SELECT {} FROM {} LEFT JOIN {} WHERE {};'.format(','.join(SELECT), ','.join(FROM),' INNER JOIN '.join(LEFT_JOIN) ,' AND '.join(WHERE))             

    results = {}
    conn = openConnection()
    cur = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    cur.execute(sql)
    for row in cur:
        try:
            _ = results[row['patent_id']]
        except:
            results[row['patent_id']] = set()
        results[row['patent_id']].update([(row['source_name'])])
        

    cur.close()
    conn.close()

    return(results)

def countSourceForPatent(patent=None):
    """
    Get patent count associated with source.

    Returns:
        dict: {'Patent_1': number of patents for source_1, 'Source_2': number of patents for source_2, ...}
    """
    return({k:len(v) for k,v in getSourceForPatent(patent=patent).items()})

def getPatentForMolecule(molecule=None):
    """
    Get patents associated with molecule.

    Returns:
        dict: {'inchikey_1':{patent_1, patent_2,...}, 'inchikey_2':{...}, ...}
    """
    SELECT = ['m.mol_inchikey','p2m.patent_id']
    FROM = ['molecule m']
    INNER_JOIN = ['pat2mol p2m ON p2m.mol_inchikey = m.mol_inchikey']
    WHERE = ['1 = 1']

    if molecule is not None:
        if isinstance(molecule, str):
            WHERE.append("p2m.mol_inchikey = ('{}')".format((molecule)))
        elif isinstance(molecule, list) or isinstance(molecule, set) or isinstance(molecule, dict):
            WHERE.append("p2m.mol_inchikey IN ('{}')".format("','".join(list(molecule)))) 
    print(WHERE)
    sql = 'SELECT {} FROM {} INNER JOIN {} WHERE {};'.format(','.join(SELECT), ','.join(FROM),' INNER JOIN '.join(INNER_JOIN) ,' AND '.join(WHERE))             


    results = {}
    conn = openConnection()
    cur = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    cur.execute(sql)
    for row in cur:
        try:
            _ = results[row['mol_inchikey']]
        except:
            results[row['mol_inchikey']] = set()
        results[row['mol_inchikey']].update([(row['patent_id'])])
        

    cur.close()
    conn.close()

    return(results)

def countPatentForMolecule(molecule=None):
    """
    Get patents associated with molecule.

    Returns:
        dict: {'inchikey_1':number of patents for inchikey_1, 'inchikey_2': ...}
    """
    return({k:len(v) for k,v in getPatentForMolecule(molecule=molecule).items()})

def getMoleculeForPatent(patent=None):
    """
    Get molecules associated with patent.

    Returns:
        dict: {patent_1:{inchikey_1,inchikey_2, ...}, patent_2:{..}, ...}
    """
    SELECT = ['p.patent_id','p2m.mol_inchikey']
    FROM = ['patent p']
    INNER_JOIN = ['pat2mol p2m ON p2m.patent_id = m.patent_id']
    WHERE = ['1 = 1']

    if patent is not None:
        if isinstance(patent, str):
            WHERE.append("p2m.patent_id = ('{}')".format("','".join(list(patent))))
        elif isinstance(patent, list) or isinstance(patent, set) or isinstance(patent, dict):
            WHERE.append("p2m.patent_id IN ('{}')".format("','".join(list(patent)))) 

    sql = 'SELECT {} FROM {} INNER JOIN {} WHERE {};'.format(','.join(SELECT), ','.join(FROM),' INNER JOIN '.join(INNER_JOIN) ,' AND '.join(WHERE))             

    results = {}
    conn = openConnection()
    cur = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    cur.execute(sql)
    for row in cur:
        try:
            _ = results[row['mol_inchikey']]
        except:
            results[row['mol_inchikey']] = set()
        results[row['mol_inchikey']].update([(row['patent_id'])])
        

    cur.close()
    conn.close()

    return(results)

def countMoleculeForPatent(patent=None):
    """
    Get patents associated with molecule.

    Returns:
        dict: {'inchikey_1':number of patents for inchikey_1, 'inchikey_2': ...}
    """
    return({k:len(v) for k,v in getMoleculeForPatent(patent=patent).items()})

def getMolecules():
    """
    Get molecules.

    Returns:
        set:   [inchikey_1, 2, 3 ...]}    
    """
    SELECT = ['molecule.mol_inchikey']
    FROM = ['molecule']
    WHERE = ['1 = 1']

    sql = 'SELECT DISTINCT {} FROM {} WHERE {};'.format(','.join(SELECT), ','.join(FROM), ' AND '.join(WHERE))
    
    conn = openConnection()
    cur = conn.cursor()
    cur.execute(sql)
    results = cur.fetchall()
    results = set([it[0] for it in results]) # r

    cur.close()
    conn.close()

    return(results)

def countMolecules():
    """
    Get molecules count

    Returns:
        int: number of molecules
    """
    return(len(getMolecules()))

def getMoleculeForSource(source=None):
    """
    Get molecules associated with patent.

    Returns:
        dict: {source_1:{inchikey_1,inchikey_2, ...}, source_2:{..}, ...}
    """
    SELECT = ['s.source_name','p2m.mol_inchikey']
    FROM = ['source s']
    INNER_JOIN = ['pat2mol p2m ON p2m.source_name = s.source_name']
    WHERE = ['1 = 1']

    if source is not None:
        if isinstance(source, str):
            WHERE.append("p2m.source_name = ('{}')".format("','".join(list(source))))
        elif isinstance(source, list) or isinstance(source, set) or isinstance(source, dict):
            WHERE.append("p2m.source_name IN ('{}')".format("','".join(list(source)))) 

    sql = 'SELECT {} FROM {} INNER JOIN {} WHERE {};'.format(','.join(SELECT), ','.join(FROM),' INNER JOIN '.join(INNER_JOIN) ,' AND '.join(WHERE))             

    results = {}
    conn = openConnection()
    cur = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    cur.execute(sql)
    for row in cur:
        try:
            _ = results[row['source_name']]
        except:
            results[row['source_name']] = set()
        results[row['source_name']].update([(row['mol_inchikey'])])
        

    cur.close()
    conn.close()

    return(results)

def countMoleculeForSource(source=None):
    """
    Get molecules associated with patent.

    Returns:
        dict: {source_1:{number of molecules for source_1}, source__2:{..}, ...}
    """
    return({k: len(v) for k, v in getMoleculeForSource(source=source).items()})

def getSourceForMolecule(molecule):
    """
    Get molecules associated with patent.

    Returns:
        dict: {inchikey_1:{source_1,source_2, ...}, inchikey_2:{..}, ...}
    """
    SELECT = ['m.mol_inchikey','p2m.source_name']
    FROM = ['molecule m']
    INNER_JOIN = ['pat2mol p2m ON p2m.mol_inchikey = m.mol_inchikey']
    WHERE = ['1 = 1']

    if molecule is not None:
        if isinstance(molecule, str):
            WHERE.append("p2m.mol_inchikey = ('{}')".format("','".join(list(molecule))))
        elif isinstance(molecule, list) or isinstance(molecule, set) or isinstance(molecule, dict):
            WHERE.append("p2m.mol_inchikey IN ('{}')".format("','".join(list(molecule)))) 

    sql = 'SELECT {} FROM {} INNER JOIN {} WHERE {};'.format(','.join(SELECT), ','.join(FROM),' INNER JOIN '.join(INNER_JOIN) ,' AND '.join(WHERE))             

    results = {}
    conn = openConnection()
    cur = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    cur.execute(sql)
    for row in cur:
        try:
            _ = results[row['mol_inchikey']]
        except:
            results[row['mol_inchikey']] = set()
        results[row['mol_inchikey']].update([(row['source_name'])])
        

    cur.close()
    conn.close()

    return(results)

def countSourceForMolecule(molecule=None):
    """
    Get molecules associated with patent.

    Returns:
        dict: {source_1:{inchikey_1,inchikey_2, ...}, source__2:{..}, ...}
    """
    return({k:len(v) for k,v in getMoleculeForSource(molecule=molecule).items()})

def getSources():
    """
    Get sources

    Returns:
        set: {source_1, source_2, ...}
    """
    SELECT = ['source.source_name']
    FROM = ['source']
    WHERE = ['1 = 1']

    sql = 'SELECT DISTINCT {} FROM {} WHERE {};'.format(','.join(SELECT), ','.join(FROM), ' AND '.join(WHERE))
    
    conn = openConnection()
    cur = conn.cursor()
    cur.execute(sql)
    results = cur.fetchall()
    results = set([result[0] for result in results])

    cur.close()
    conn.close()

    return(results)

def countSources():
    """
    Get sources count

    Returns:
        int: sources count 
    """
    return(len(getSources()))