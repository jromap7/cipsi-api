#!/usr/bin/env python

'''
cipsi.py: Chemical Intellectual Property Service @ IMIM. Connections with
          the CIPSI db.
'''

__author__ = 'María Martínez'
__date__ = '23 Jun 2022'

import pandas as pd
import psycopg2
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import IPythonConsole

def get_cipsi_connection():
    '''
    Get CIPSI Connection

    Returns:
        conn:   connection to CIPSI db to be used in the queries.
    '''
    
    conn = psycopg2.connect(
        host="blackadder",
        database="cipsi_test",
        user='syspharm',
        password='whatever')

    return(conn)

def identicalquerybyinchiKey(conn, ID, currentPage, pageSize, orderColumn, orderDirection):
    '''
    Perform identical search when the user enters an inchiKey.

    Args:
        conn:           connection to the CIPSI database.
        ID (string):    molecule inchiKey.
        currentPage:    the page number of the table the user is in.
        pageSize:       the number of rows the user wants to display in the table.

    Returns:
        result (pandas dataframe): {}
    '''

    limit = int(pageSize)
    offset = int(pageSize)*int(currentPage)

    if orderColumn == 'source_name':
        orderTable = 'pm'
    elif orderColumn == 'mol_inchikey' or orderColumn == 'mol_name' or orderColumn == 'mol_mol':
        orderTable = 'm'
    else:
        orderTable = 'p'

    query = '''select m.mol_inchikey, m.mol_name, m.mol_mol, p.patent_id, p.patent_title, pm.source_name, p.geographic_location, p.pub_date, p.priority_date, p.grant_date, p.filing_date, p.abstract_text, p.inventor_name, p.assignee_name_orig, p.url
            from molecule m
            inner join pat2mol pm on pm.mol_inchikey = m.mol_inchikey
            inner join patent p on p.patent_id = pm.patent_id
            where m.mol_inchikey = \'''' + ID + '''\'
            order by ''' + orderTable + '''.''' + orderColumn + ''' ''' + orderDirection + '''
            limit ''' + str(limit) + '''
            offset ''' + str(offset) + '''
    '''

    cur = conn.cursor()
    cur.execute(query)
    result = pd.DataFrame(cur.fetchall())
    cur.close()

    if not result.empty:
        result.columns = ["mol_inchikey", "mol_name", "mol_mol", "patent_id", "patent_title", "source_name", "geographic_location", "pub_date", "priority_date", "grant_date", "filing_date", "abstract_text", "inventor_name", "assignee_name_orig", "url"]
    
    return(result)

def identicalqueryrowsbyinchiKey(conn, ID):
    '''
    Compute total rows of the identical search by inchiKey.

    Args:
        conn:           connection to the CIPSI database.
        ID (string):    molecule inchiKey.

    Returns:
        result (int):   number of total registers found.
    '''
    
    query = '''select count(*)
            from molecule m
            inner join pat2mol pm on pm.mol_inchikey = m.mol_inchikey
            inner join patent p on p.patent_id = pm.patent_id
            where m.mol_inchikey = \'''' + ID + '''\'
    '''

    cur = conn.cursor()
    cur.execute(query)
    result = cur.fetchall()[0][0]
    cur.close()

    return(result)

def getSmiles(conn, ID):
    '''
    Get the smiles of an inchiKey

    Args:
        conn:           connection to the CIPSI database.
        ID (string):    molecule inchiKey.

    Returns:
        result (pandas dataframe): {}
    '''

    query = '''SELECT mol_mol
            from molecule
            where mol_inchikey = \'''' + ID + '''\''''

    cur = conn.cursor()
    cur.execute(query)
    result = cur.fetchall()
    cur.close()

    print(result)
    
    if result:
        result = result[0][0]

    return(result)

def similarityquerybysmiles(conn, ID, currentPage, pageSize, orderColumn, orderDirection, similarityThrs):
    '''
    Perform similarity search when the user enters a smiles.

    Args:
        conn:           connection to the CIPSI database.
        ID (string):    molecule smiles.
        currentPage:    the page number of the table the user is in.
        pageSize:       the number of rows the user wants to display in the table.

    Returns:
        result (pandas dataframe): {}
    '''
    
    limit = int(pageSize)
    offset = int(pageSize)*int(currentPage)

    if orderColumn == 'source_name':
        orderTable = 'pm.'
    elif orderColumn == 'mol_inchikey' or orderColumn == 'mol_name' or orderColumn == 'mol_mol':
        orderTable = 'm.'
    elif orderColumn == 'sml':
        orderTable = ''
    else:
        orderTable = 'p.'

    query = '''select m.mol_inchikey, m.mol_name, tanimoto_sml(morganbv_fp('C1COCCO1'), mol_ecfp4) as sml, m.mol_mol, p.patent_id, p.patent_title, pm.source_name, p.geographic_location, p.pub_date, p.priority_date, p.grant_date, p.filing_date, p.abstract_text, p.inventor_name, p.assignee_name_orig, p.url
            from molecule m
            inner join pat2mol pm on pm.mol_inchikey = m.mol_inchikey
            inner join patent p on p.patent_id = pm.patent_id
            where mol_ecfp4 % morganbv_fp(\'''' + ID + '''\') and tanimoto_sml(morganbv_fp('C1COCCO1'), mol_ecfp4) > ''' + similarityThrs + '''
            order by ''' + orderTable + orderColumn + ''' ''' + orderDirection + '''
            limit ''' + str(limit) + '''
            offset ''' + str(offset) + '''
    '''

    cur = conn.cursor()
    cur.execute(query)
    result = pd.DataFrame(cur.fetchall())
    cur.close()

    if not result.empty:
        result.columns = ["mol_inchikey", "mol_name", "sml", "mol_mol", "patent_id", "patent_title", "source_name", "geographic_location", "pub_date", "priority_date", "grant_date", "filing_date", "abstract_text", "inventor_name", "assignee_name_orig", "url"]

    return(result)

def similarityqueryrowsbysmiles(conn, ID, similarityThrs):
    '''
    Compute total rows of the similarity search by smiles.

    Args:
        conn:                       connection to the CIPSI database.
        ID (string):                molecule inchiKey.
        similarityThrs (string):    threshold to filter similarity.

    Returns:
        result (int):   number of total registers found.
    '''
    
    query = '''select count(*)
            from molecule m
            inner join pat2mol pm on pm.mol_inchikey = m.mol_inchikey
            inner join patent p on p.patent_id = pm.patent_id
            where mol_ecfp4 % morganbv_fp(\'''' + ID + '''\') and tanimoto_sml(morganbv_fp('C1COCCO1'), mol_ecfp4) > ''' + similarityThrs + '''
    '''

    cur = conn.cursor()
    cur.execute(query)
    result = cur.fetchall()[0][0]
    cur.close()

    return(result)

def substructurequerybysmiles(conn, ID, currentPage, pageSize, orderColumn, orderDirection):
    '''
    Perform substructure search when the user enters a smiles.

    Args:
        conn:           connection to the CIPSI database.
        ID (string):    molecule smiles.
        currentPage:    the page number of the table the user is in.
        pageSize:       the number of rows the user wants to display in the table.

    Returns:
        result (pandas dataframe): {}
    '''

    limit = int(pageSize)
    offset = int(pageSize)*int(currentPage)

    if orderColumn == 'source_name':
        orderTable = 'pm.'
    elif orderColumn == 'mol_inchikey' or orderColumn == 'mol_name' or orderColumn == 'mol_mol':
        orderTable = 'm.'
    elif orderColumn == 'sml':
        orderTable = ''
    else:
        orderTable = 'p.'

    query = '''select m.mol_inchikey, m.mol_name, m.mol_mol, p.patent_id, p.patent_title, pm.source_name, p.geographic_location, p.pub_date, p.priority_date, p.grant_date, p.filing_date, p.abstract_text, p.inventor_name, p.assignee_name_orig, p.url
            from molecule m
            inner join pat2mol pm on pm.mol_inchikey = m.mol_inchikey
            inner join patent p on p.patent_id = pm.patent_id
            where m.mol_mol@>\'''' + ID + '''\'
            order by ''' + orderTable + orderColumn + ''' ''' + orderDirection + '''
            limit ''' + str(limit) + '''
            offset ''' + str(offset) + '''
    '''

    cur = conn.cursor()
    cur.execute(query)
    result = pd.DataFrame(cur.fetchall())
    cur.close()

    if not result.empty:
        result.columns = ["mol_inchikey", "mol_name", "mol_mol", "patent_id", "patent_title", "source_name", "geographic_location", "pub_date", "priority_date", "grant_date", "filing_date", "abstract_text", "inventor_name", "assignee_name_orig", "url"]


    return(result)

def substructurequeryrowsbysmiles(conn, ID):
    '''
    Compute total rows of the substructure search by smiles.

    Args:
        conn:                       connection to the CIPSI database.
        ID (string):                molecule inchiKey.
        similarityThrs (string):    threshold to filter similarity.

    Returns:
        result (int):   number of total registers found.
    '''

    query = '''select count(*)
            from molecule m
            inner join pat2mol pm on pm.mol_inchikey = m.mol_inchikey
            inner join patent p on p.patent_id = pm.patent_id
            where m.mol_mol@>\'''' + ID + '''\'
    '''
    print(query)

    cur = conn.cursor()
    cur.execute(query)
    result = cur.fetchall()[0][0]
    cur.close()

    return(result)