#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 10:30:11 2020

@author: vivekmodi
"""

from Rascore import db, Clusters
import sys
import pandas as pd

# df['H1 Cluster'].fillna('NA',inplace=True);df['H2 Cluster'].fillna('NA',inplace=True);df['H3 Cluster'].fillna('NA',inplace=True)
def update_database(df):
    db.drop_all()
    db.create_all()
    for i in df.index:
        ras_clus=Clusters(pdb=df.at[i,'PDB ID'],\
        protein_name=df.at[i,'Protein Name'],\
        bound_protein=df.at[i,'Bound Protein'],\
        mutation_status=df.at[i,'Mutation Status'],\
        nucleotide_class=df.at[i,'Nucleotide Class'],\
        drug_class=df.at[i,'Drug Class'],\
        bound_protein_class=df.at[i,'Bound Protein Class'],\
        homodimer_status=df.at[i,'Homodimer Status'],\
        switch1_cluster=df.at[i,'Switch 1 Cluster'],\
        switch2_cluster=df.at[i,'Switch 2 Cluster'],\
        conformational_state=df.at[i,'Conformational State'])

        db.session.add(ras_clus)
        db.session.commit()

if __name__ == '__main__':
    filename=sys.argv[1]
    df=pd.read_csv(filename,sep='\t',header='infer')
    update_database(df)
