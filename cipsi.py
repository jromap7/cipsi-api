#!/usr/bin/env python

'''
cipsi.py: Chemical Intellectual Property Service @ IMIM. Connections with
          the CIPSI db.
'''

__author__ = 'María Martínez'
__date__ = '23 Jun 2022'

import pandas as pd
import psycopg2

def get_cipsi_connection():
    '''
    Get CIPSI Connection

    Returns:
        conn:   connection to CIPSI db to be used in the queries.
    '''
    
    conn = psycopg2.connect(
        host="blackadder",
        database="cipsi_db",
        user='syspharm',
        password='whatever')

    return(conn)

def identicalquerybyinchiKey(conn, ID):
    '''
    Perform identical search when the user enters an inchiKey.

    Args:
        conn:           connection to the CIPSI database.
        ID (string):    molecule inchiKey.

    Returns:
        result (pandas dataframe): {}
    '''
    
    SELECT = ['m.mol_inchikey', 'm.mol_name', 'm.mol_mol', 'pm.patent_id']
    FROM = ['molecule m']
    INNER_JOIN = ['pat2mol pm ON pm.mol_inchikey = m.mol_inchikey']
    WHERE = ['m.mol_inchikey = \'' + ID + '\' limit 50']

    query = 'SELECT {} FROM {} INNER JOIN {} WHERE {};'.format(','.join(SELECT), ','.join(FROM), 
                                                             ' INNER JOIN '.join(INNER_JOIN), ' AND '.join(WHERE))

    cur = conn.cursor()
    cur.execute(query)
    print('query executed')
    result = pd.DataFrame(cur.fetchall())
    print('results fetched')
    cur.close()
    conn.close()

    result.columns = ["inchiKey", "schemblID", "smiles", "patentID"]

    return(result)