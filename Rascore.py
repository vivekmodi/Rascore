#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 10:30:11 2020

@author: vivekmodi
"""

import os, subprocess, sys, glob, random, numpy as np
from datetime import datetime
from flask import Flask, render_template, url_for, redirect, request
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate
#from flask_wtf import FlaskForm
#from wtforms import StringField, SubmitField
#from wtforms.validators import DataRequired
#from werkzeug.utils import secure_filename
#from sqlalchemy import func
#from collections import defaultdict
#from flask import Markup
#from sqlalchemy import desc,asc
import pandas as pd
from natsort import natsorted

pwd=os.getcwd()
sys.path.append(pwd+'/scripts')

UPLOAD_FOLDER = (pwd+'/server/uploads')
ALLOWED_EXTENSIONS = {'txt', 'pdb','cif'}

app=Flask(__name__)
app.config['SECRET_KEY'] = 'mysecretkey'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['ALLOWED_EXTENSIONS'] = ALLOWED_EXTENSIONS
app.jinja_env.add_extension('jinja2.ext.loopcontrols')

############SQL Databse and Models############################
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///ras_data.sqlite'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

db = SQLAlchemy(app)
Migrate(app,db)

class Clusters(db.Model):
    pdb=db.Column(db.Text,primary_key=True)
    protein_name=db.Column(db.Text)
    bound_protein=db.Column(db.Text)
    mutation_status=db.Column(db.Text)
    nucleotide_class=db.Column(db.Text)
    drug_class=db.Column(db.Text)
    bound_protein_class=db.Column(db.Text)
    homodimer_status=db.Column(db.Text)
    switch1_cluster=db.Column(db.Text)
    switch2_cluster=db.Column(db.Text)
    conformational_state=db.Column(db.Text)

    
    def __init__(self,pdb=pdb,protein_name=protein_name,bound_protein=bound_protein,mutation_status=mutation_status,nucleotide_class=nucleotide_class,drug_class=drug_class,\
                 bound_protein_class=bound_protein_class,homodimer_status=homodimer_status,switch1_cluster=switch1_cluster,\
                 switch2_cluster=switch2_cluster,conformational_state=conformational_state):

        self.pdb=pdb;self.protein_name=protein_name;self.bound_protein=bound_protein;self.mutation_status=mutation_status;self.nucleotide_class=nucleotide_class;
        self.drug_class=drug_class;self.bound_protein_class=bound_protein_class;self.homodimer_status=homodimer_status;
        self.switch1_cluster=switch1_cluster;self.switch2_cluster=switch2_cluster;self.conformational_state=conformational_state;

    def __repr__(self):
        return f'{self.pdb} {self.protein_name} {self.bound_protein} {self.mutation_status} {self.nucleotide_class} {self.drug_class} {self.bound_protein_class}\
            {self.homodimer_status} {self.switch1_cluster} {self.switch2_cluster} {self.conformational_state}'

def create_lists():
    pdbListDb=list()
    for pdbs in Clusters.query.with_entities(Clusters.pdb):
        #pdbid,chainid,cdr=pdbs[0].split('_')
        pdbListDb.append(pdbs[0][:4])
#     for cluster in perCDR.query.with_entities(perCDR.cluster):
#         clusterListDb.append(cluster[0])
#         cluster_cdr_list.append(cluster[0])
#         cdr_length='-'.join(cluster[0].split('-')[0:2])  #Get the first two element after the split like H1-10 from H1-10-1
#         cluster_cdr_list.append(cdr_length)     #This list contains both clusters and cdr-lengths

#     clusterListDb=list(natsorted(set(sorted(clusterListDb))))
#     cluster_cdr_list=list(natsorted(set(sorted(cluster_cdr_list))))
#     pdbListDb=list(natsorted(set(sorted(pdbListDb))))

#     for cluster in cluster_cdr_list:
#             clusterChainCount[cluster]=perCDR.query.filter(perCDR.cluster.contains(cluster)).count()  #Incorrect count for cdr-length due to substring mismatch

#     clusterChainCount['H1-All']=perCDR.query.filter(perCDR.cluster.contains('H1')).count()
#     clusterChainCount['H2-All']=perCDR.query.filter(perCDR.cluster.contains('H2')).count()
#     clusterChainCount['H3-All']=perCDR.query.filter(perCDR.cluster.contains('H3')).count()
#     clusterChainCount['L1-All']=perCDR.query.filter(perCDR.cluster.contains('L1')).count()
#     clusterChainCount['L2-All']=perCDR.query.filter(perCDR.cluster.contains('L2')).count()
#     clusterChainCount['L3-All']=perCDR.query.filter(perCDR.cluster.contains('L3')).count()
    return pdbListDb

def count_entries(Clusters):
    pdb_list=dict();entry_count=dict()
    pdb_list['All']=list();pdb_list['KRAS']=list();pdb_list['HRAS']=list();pdb_list['NRAS']=list();
    pdb_list['All']=[pdbs.pdb[:4] for pdbs in Clusters.query.all()]
    pdb_list['HRAS']=[pdbs.pdb[:4] for pdbs in Clusters.query.filter(Clusters.protein_name.contains('HRAS'))]
    pdb_list['KRAS']=[pdbs.pdb[:4] for pdbs in Clusters.query.filter(Clusters.protein_name.contains('KRAS'))]
    pdb_list['NRAS']=[pdbs.pdb[:4] for pdbs in Clusters.query.filter(Clusters.protein_name.contains('NRAS'))]
    entry_count['All']=len(set(pdb_list['All']));entry_count['HRAS']=len(set(pdb_list['HRAS']))
    entry_count['KRAS']=len(set(pdb_list['KRAS']));entry_count['NRAS']=len(set(pdb_list['NRAS']))
    
    return entry_count
    
    
def write_text_file(sublist,tsvFile,cdr_format):
    print(cdr_format)
    fhandle_textFile=open(f'{pwd}/static/{tsvFile}','w')
    if cdr_format:     #Header for perCDR files
        fhandle_textFile.write('CDR\tCDR Length\tPDB\tChain\tResolution\tAHO Resnum\tAuthor Resnum\tSequence\tGermline Sequence\tGene\tPDB Species\tFrame Germline\tCluster\tDistance\tCDR Germline\tCDR Seqid\tRama4\tBeta Turns\tMinimum EDIA\n')
        for item in sublist:
            fhandle_textFile.write(f'{item.cdr}\t{item.cdr_length}\t{item.pdb}\t{item.chain}\t{item.resolution}\t{item.aho_resnum}\t\
                                   {item.author_resnum}\t{item.sequence}\t{item.germline_sequence}\t{item.gene}\t{item.pdb_species}\t\
                                   {item.frame_germline}\t{item.cluster}\t{item.distance}\t{item.cdr_germline}\t{item.cdr_seqid}\t\
                                   {item.rama4}\t{item.beta_turns}\t{item.minimum_edia}\n')

    else:    #Header for perVRegion files
        fhandle_textFile.write('PDB\tVH Chain\tVH Framework\tVH Framework Seqid\tH1 Cluster\tH1 Cluster Distance\tH1 Seqid\tH2 Cluster\tH2 Cluster Distance\tH2 Seqid\tH3 Cluster\tH3 Cluster Distance\tVL Chain\tVL Framework\tVL Framework Seqid\tL1 Cluster\tL1 Cluster Distance\tL1 Seqid\tL2 Cluster\tL2 Cluster Distance\tL2 Seqid\tL3 Cluster\tL3 Cluster Distance\n')
        for item in sublist:
            fhandle_textFile.write(f'{item.pdb}\t{item.vh_chain}\t{item.vh_framework}\t{item.vh_framework_seqid}\t{item.h1_cluster}\t{item.h1_cluster_distance}\t{item.h1_seqid}\t\
                 {item.h2_cluster}\t{item.h2_cluster_distance}\t{item.h2_seqid}\t{item.h3_cluster}\t{item.h3_cluster_distance}\t{item.vl_chain}\t{item.vl_framework}\t{item.vl_framework_seqid}\t\
                 {item.l1_cluster}\t{item.l1_cluster_distance}\t{item.l1_seqid}\t{item.l2_cluster}\t{item.l2_cluster_distance}\t{item.l2_seqid}\t{item.l3_cluster}\t{item.l3_cluster_distance}\n')
    fhandle_textFile.close()

@app.route('/index')
@app.route('/home')
@app.route('/')
def home():
    return render_template('home.html')

@app.route('/about')
def about():
    return render_template('about.html')

@app.route('/browse')
def browse():
    return redirect(url_for('uniqueQuery',settings='All',queryname='All'))

    


@app.route('/statistics')
def statistics():
#     pdb_species_unique=dict();gene_unique=dict();entries=dict()
#     (pdbListDb,chainListDb,clusterListDb,clusterChainCount,cluster_cdr_list)=create_lists()
#     for cluster in clusterListDb:
#         #cluster_list=perCDR.query.filter(perCDR.cluster.contains(queryname)).all()
#         entries[cluster]=perCDR.query.filter(perCDR.cluster.contains(cluster)).count()
#         pdb_species_list=perCDR.query.filter(perCDR.cluster.contains(cluster)).with_entities('pdb_species')
#         pdb_species_unique[cluster]=pd.Series(names[0] for names in pdb_species_list).unique()
#         gene_list=perCDR.query.filter(perCDR.cluster.contains(cluster)).with_entities('gene')
#         gene_unique[cluster]=pd.Series(names[0] for names in gene_list).unique()

    return('home.html')
#     return render_template('statistics.html',clusterListDb=clusterListDb,entries=entries,pdb_species_unique=pdb_species_unique,gene_unique=gene_unique)

@app.route('/formSearch', methods=['GET','POST'])
def formSearch():
    (pdbListDb)=create_lists()

    if request.method=='POST':
        inputString=request.form['inputString'].lower()
        print(pdbListDb,inputString)
        if inputString in pdbListDb:    #match without chain
            return redirect(url_for('uniqueQuery',queryname=inputString,settings='PDB'))
       
        else:
            return render_template('nomatch.html')

    return render_template('search.html')

# @app.route('/formSearchMultiple',methods=['GET','POST'])
# def formSearchMultiple():
#     (pdbListDb,chainListDb,clusterListDb,clusterChainCount,cluster_cdr_list)=create_lists()
#     if request.method=='POST':
#         cdr_select=request.form['cdr_select']

#         if cdr_select=='All':
#             return redirect(url_for('browse'))
#         elif 'All' in cdr_select:
#             cdr_select=cdr_select[0:-4]
#             return redirect(url_for('uniqueQuery',settings='cdr',queryname=cdr_select))
#         elif cdr_select in clusterListDb:
#             return redirect(url_for('uniqueQuery',settings='cluster',queryname=cdr_select))
#         else:
#             return redirect(url_for('uniqueQuery',settings='cdr_length',queryname=cdr_select))

#     return render_template ('search.html', pdbListDb=pdbListDb, clusterListDb=clusterListDb, clusterChainCount=clusterChainCount, cluster_cdr_list=cluster_cdr_list)


@app.route('/webserver', methods=['GET','POST'])
def webserver():
    return render_template('webserver.html')

@app.route('/download')
def download():
    return render_template('download.html')

@app.route('/help')
def help():
    return render_template('help.html')

@app.route('/contact')
def contact():
    return render_template('contact.html')

@app.route('/dunbrackLab')
def dunbrackLab():
    return redirect("http://dunbrack.fccc.edu")


@app.route('/<settings>/<queryname>')
def uniqueQuery(settings,queryname):
    if settings=='PDB':
        queryname=queryname.upper()
        if len(queryname)==5:
            queryname=queryname[0:4]           #Remove chain so that only PDB id is always searched
        pdb_list=Clusters.query.filter(Clusters.pdb.contains(queryname)).all()
        pdb_count=Clusters.query.filter(Clusters.pdb.contains(queryname)).count()
        #pdb_resolution=perCDR.query.filter(perCDR.pdb.contains(queryname)).with_entities('resolution').first()[0]

        #vl_species=perCDR.query.filter(perCDR.pdb.contains(queryname),perCDR.cdr.contains('L')).with_entities('pdb_species').first()[0]
        #tsvFile=f'downloads/text-files/{queryname}.tab'
        #write_text_file(pdb_list,tsvFile,True)

        return render_template('pdbs.html',queryname=queryname,pdb_list=pdb_list,pdb_count=pdb_count)

    if settings=='All':
        retrieve_str=dict();chain_count=dict();entry_count=dict()
        retrieve_str['All']=Clusters.query.all();chain_count['All']=Clusters.query.count()
        retrieve_str['HRAS']=Clusters.query.filter(Clusters.protein_name.contains('HRAS'));chain_count['HRAS']=Clusters.query.filter(Clusters.protein_name.contains('HRAS')).count()
        retrieve_str['KRAS']=Clusters.query.filter(Clusters.protein_name.contains('KRAS'));chain_count['KRAS']=Clusters.query.filter(Clusters.protein_name.contains('KRAS')).count()
        retrieve_str['NRAS']=Clusters.query.filter(Clusters.protein_name.contains('NRAS'));chain_count['NRAS']=Clusters.query.filter(Clusters.protein_name.contains('NRAS')).count()
        entry_count=count_entries(Clusters)
    
        return render_template('browse.html',retrieve_str=retrieve_str,chain_count=chain_count,entry_count=entry_count)

        
#     if settings=='cluster':
#         cluster_list=perCDR.query.filter(perCDR.cluster.contains(queryname)).all()
#         cluster_cdr_length=perCDR.query.filter(perCDR.cluster.contains(queryname)).with_entities('cdr_length').first()[0]
#         cluster_cdr=perCDR.query.filter(perCDR.cluster.contains(queryname)).with_entities('cdr').first()[0]
#         pdb_count=perCDR.query.filter(perCDR.cluster.contains(queryname)).with_entities('pdb')
#         pdb_unique_count=len(pd.Series(names[0] for names in pdb_count).unique())
#         chain_count=perCDR.query.filter(perCDR.cluster.contains(queryname)).count()    #This is number of chains with a given cluster
#         per_loop=percent_loop_length(chain_count,queryname)
#         seq_in_cluster=perCDR.query.filter(perCDR.cluster.contains(queryname)).with_entities('sequence')
#         seq_unique_count=len(pd.Series(names[0] for names in seq_in_cluster).unique())
#         pdb_species_list=perCDR.query.filter(perCDR.cluster.contains(queryname)).with_entities('pdb_species')
#         pdb_species_unique=pd.Series(names[0] for names in pdb_species_list).sort_values().unique()
#         gene_list=perCDR.query.filter(perCDR.cluster.contains(queryname)).with_entities('gene')
#         gene_unique=pd.Series(names[0] for names in gene_list).unique()
#         rama_list=perCDR.query.filter(perCDR.cluster.contains(queryname)).with_entities('rama4')
#         (most_common_rama_count,most_common_rama_string )=most_common_rama(rama_list)
#         tsvFile=f'downloads/text-files/{queryname}.tab'
#         write_text_file(cluster_list,tsvFile,True)
#         return render_template('cluster.html',queryname=queryname,cluster_list=cluster_list,cluster_cdr_length=cluster_cdr_length,\
#                                cluster_cdr=cluster_cdr,pdb_unique_count=pdb_unique_count,chain_count=chain_count,per_loop=per_loop,\
#                                    seq_unique_count=seq_unique_count,pdb_species_unique=pdb_species_unique,gene_unique=gene_unique,\
#                                        most_common_rama_count=most_common_rama_count,most_common_rama_string=most_common_rama_string,tsvFile=tsvFile)

#     if settings=='cdr':
#         cdr_list=perCDR.query.filter(perCDR.cdr.contains(queryname)).all()
#         pdb_count=perCDR.query.filter(perCDR.cdr.contains(queryname)).with_entities('pdb')
#         pdb_unique_count=len(pd.Series(names[0] for names in pdb_count).unique())
#         chain_count=perCDR.query.filter(perCDR.cdr.contains(queryname)).count()    #This is number of chains with a given cluster
#         pdb_species_list=perCDR.query.filter(perCDR.cdr.contains(queryname)).with_entities('pdb_species')
#         pdb_species_unique=pd.Series(names[0] for names in pdb_species_list).unique()
#         gene_list=perCDR.query.filter(perCDR.cdr.contains(queryname)).with_entities('gene')
#         gene_unique=pd.Series(names[0] for names in gene_list).unique()
#         tsvFile=f'downloads/text-files/{queryname}.tab'
#         write_text_file(cdr_list,tsvFile,True)
#         return render_template('cdr.html',queryname=queryname,cdr_list=cdr_list,pdb_unique_count=pdb_unique_count,chain_count=chain_count,\
#                                pdb_species_unique=pdb_species_unique,gene_unique=gene_unique,tsvFile=tsvFile)

#     if settings=='cdr_length':
#         cdr_length_list=perCDR.query.filter(perCDR.cdr_length.contains(queryname)).all()
#         cdr_cdr_length=perCDR.query.filter(perCDR.cdr_length.contains(queryname)).with_entities('cdr').first()[0]
#         pdb_count=perCDR.query.filter(perCDR.cdr_length.contains(queryname)).with_entities('pdb')
#         pdb_unique_count=len(pd.Series(names[0] for names in pdb_count).unique())
#         chain_count=perCDR.query.filter(perCDR.cdr_length.contains(queryname)).count()    #This is number of chains with a given cluster
#         pdb_species_list=perCDR.query.filter(perCDR.cdr_length.contains(queryname)).with_entities('pdb_species')
#         pdb_species_unique=pd.Series(names[0] for names in pdb_species_list).unique()
#         gene_list=perCDR.query.filter(perCDR.cdr_length.contains(queryname)).with_entities('gene')
#         gene_unique=pd.Series(names[0] for names in gene_list).unique()
#         tsvFile=f'downloads/text-files/{queryname}.tab'
#         write_text_file(cdr_length_list,tsvFile,True)
#         return render_template('cdr_length.html',queryname=queryname,cdr_cdr_length=cdr_cdr_length,cdr_length_list=cdr_length_list,\
#                                pdb_unique_count=pdb_unique_count,chain_count=chain_count,pdb_species_unique=pdb_species_unique,\
#                                    gene_unique=gene_unique,tsvFile=tsvFile)

#     if settings=='germline':
#         queryname=queryname.split('*')[0]
#         germline_list=perVRegion.query.filter(perVRegion.vh_framework.contains(queryname)).all()

#         if germline_list:    #First check for VH germline and then for VL
#             tsvFile=f'downloads/text-files/{queryname}.tab'
#             write_text_file(germline_list,tsvFile,False)
#             return render_template('germline.html',queryname=queryname,germline_list=germline_list,tsvFile=tsvFile)
#         else:
#             germline_list=perVRegion.query.filter(perVRegion.vl_framework.contains(queryname)).all()
#             tsvFile=f'downloads/text-files/{queryname}.tab'
#             write_text_file(germline_list,tsvFile,False)
#             return render_template('germline.html',queryname=queryname,germline_list=germline_list,tsvFile=tsvFile)

#     if settings=='frame_germline':
#         queryname=queryname.split('*')[0]
#         germline_list=perCDR.query.filter(perCDR.frame_germline.contains(queryname)).all()
#         tsvFile=f'downloads/text-files/{queryname}.tab'
#         write_text_file(germline_list,tsvFile,True)

#         return render_template('germline.html',queryname=queryname,germline_list=germline_list,tsvFile=tsvFile)

#     if settings=='cdr_germline':
#         queryname=queryname.split('*')[0]
#         germline_list=perCDR.query.filter(perCDR.cdr_germline.contains(queryname)).all()
#         tsvFile=f'downloads/text-files/{queryname}.tab'
#         write_text_file(germline_list,tsvFile,True)

#         return render_template('germline.html',queryname=queryname,germline_list=germline_list,tsvFile=tsvFile)


if __name__ == '__main__':
    app.run(debug=True)
