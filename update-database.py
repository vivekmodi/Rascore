#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 10:30:11 2020

@author: vivekmodi
"""

from Rascore import db, Clusters
import sys
import pandas as pd
import pickle


# df['H1 Cluster'].fillna('NA',inplace=True);df['H2 Cluster'].fillna('NA',inplace=True);df['H3 Cluster'].fillna('NA',inplace=True)
def update_database(df):
    db.drop_all()
    db.create_all()
    for i in df.index:
        ras_clus=Clusters(pdb=df.at[i,'PDB ID'],\
        gene_name=df.at[i,'Gene Name'],\
        protein_name=df.at[i,'Protein Name'],\
        bound_protein=df.at[i,'Bound Protein'],\
        bound_peptide=df.at[i,'Bound Peptide'],\
        mutation_status=df.at[i,'Mutation Status'],\
        nucleotide_class=df.at[i,'Nucleotide Class'],\
        drug_class=df.at[i,'Drug Class'],\
        bound_protein_class=df.at[i,'Bound Protein Class'],\
        homodimer_status=df.at[i,'Homodimer Status'],\
        spatial_y32=df.at[i,'Spatial Y32'],\
        spatial_y71=df.at[i,'Spatial Y71'],\
        dihedral_switch1=df.at[i,'Dihedral Switch 1'],\
        dihedral_switch2=df.at[i,'Dihedral Switch 2'],\
        residue_range=df.at[i,'Residue Range'],\
        nucleotide=df.at[i,'Nucleotide'],\
        drug=df.at[i,'Drug'],\
        other_ligands=df.at[i,'Other Ligands'],\
        experiment_type=df.at[i,'Experiment Type'],\
        resolution=df.at[i,'Resolution'],\
        rfactor=df.at[i,'R-Factor'],\
        crystal_form=df.at[i,'Crystal Form'],\
        deposit_date=df.at[i,'Deposit Date'],\
        pmid=df.at[i,'PMID'])

        db.session.add(ras_clus)
        db.session.commit()

# def create_search_dict(Clusters):
#     pdb_set=set();conformation_set=set();bound_protein_set=set();bound_protein_class_set=set();nucleotide_class_set=set();drug_class_set=set()
#     mutation_set=set()

#     for items in Clusters.query.all():
#         pdb_set.add(items.pdb[:4])
#         conformation_set.add(items.conformational_state.strip())
#         nucleotide_class_set.add(items.nucleotide_class.strip())
#         drug_class_set.add(items.drug_class.strip())
#         bound_protein_class_set.add(items.bound_protein_class.strip())   #Check if these variables can also have ',' in them

#         if ',' in items.bound_protein:
#             for protein in items.bound_protein.split(','):
#                 bound_protein_set.add(protein.strip())
#         else:
#             bound_protein_set.add(items.bound_protein.strip())

#         if ',' in items.mutation_status:
#             for mutation in items.mutation_status.split(','):
#                 mutation_set.add(mutation.strip())
#         else:
#             mutation_set.add(items.mutation_status.strip())

#     count_entire_set=dict()
#     for conf in conformation_set:
#         for mut in mutation_set:
#             for bp in bound_protein_set:
#                 for bpc in bound_protein_class_set:
#                     for nucl in nucleotide_class_set:
#                         for drug in drug_class_set:
#                            if Clusters.query.filter(Clusters.conformational_state.contains(conf),Clusters.mutation_status.contains(mut),Clusters.bound_protein.contains(bp),Clusters.bound_protein_class.contains(bpc),Clusters.nucleotide_class.contains(nucl),Clusters.drug_class.contains(drug)).count()>0:
#                                count_entire_set[conf,mut,bp,bpc,nucl,drug]=Clusters.query.filter(Clusters.conformational_state.contains(conf),Clusters.mutation_status.contains(mut),Clusters.bound_protein.contains(bp),Clusters.bound_protein_class.contains(bpc),Clusters.nucleotide_class.contains(nucl),Clusters.drug_class.contains(drug)).count()

#     fhandle_pickle=open('search_dict.pkl','w')
#     pickle.dump(count_entire_set,fhandle_pickle)
#     fhandle_pickle.close()

if __name__ == '__main__':
    filename=sys.argv[1]
    df=pd.read_csv(filename,sep='\t',header='infer')
    update_database(df)
