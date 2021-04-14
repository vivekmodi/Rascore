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
import pandas as pd
from natsort import natsorted
import pickle

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
    gene_name=db.Column(db.Text)
    protein_name=db.Column(db.Text)
    bound_protein=db.Column(db.Text)
    bound_peptide=db.Column(db.Text)
    mutation_status=db.Column(db.Text)
    nucleotide_class=db.Column(db.Text)
    drug_class=db.Column(db.Text)
    bound_protein_class=db.Column(db.Text)
    homodimer_status=db.Column(db.Text)
    switch1_cluster=db.Column(db.Text)
    switch2_cluster=db.Column(db.Text)
    conformational_state=db.Column(db.Text)
    residue_range=db.Column(db.Text)
    nucleotide=db.Column(db.Text)
    drug=db.Column(db.Text)
    other_ligands=db.Column(db.Text)
    experiment_type=db.Column(db.Text)
    resolution=db.Column(db.Text)
    rfactor=db.Column(db.Text)
    crystal_form=db.Column(db.Text)
    deposit_date=db.Column(db.Text)
    pmid=db.Column(db.Text)

    
    def __init__(self,pdb=pdb,gene_name=gene_name,protein_name=protein_name,bound_protein=bound_protein,bound_peptide=bound_peptide,mutation_status=mutation_status,nucleotide_class=nucleotide_class,drug_class=drug_class,\
                 bound_protein_class=bound_protein_class,homodimer_status=homodimer_status,switch1_cluster=switch1_cluster,\
                 switch2_cluster=switch2_cluster,conformational_state=conformational_state,residue_range=residue_range,nucleotide=nucleotide,\
                 drug=drug,other_ligands=other_ligands,experiment_type=experiment_type,resolution=resolution,rfactor=rfactor,\
                 crystal_form=crystal_form,deposit_date=deposit_date,pmid=pmid):

        self.pdb=pdb;self.gene_name=gene_name;self.protein_name=protein_name;self.bound_protein=bound_protein;self.bound_peptide=bound_peptide;
        self.mutation_status=mutation_status;self.nucleotide_class=nucleotide_class;
        self.drug_class=drug_class;self.bound_protein_class=bound_protein_class;self.homodimer_status=homodimer_status;
        self.switch1_cluster=switch1_cluster;self.switch2_cluster=switch2_cluster;self.conformational_state=conformational_state;
        self.residue_range=residue_range;self.nucleotide=nucleotide;self.drug=drug;self.other_ligands=other_ligands;self.experiment_type=experiment_type;
        self.resolution=resolution;self.rfactor=rfactor;self.crystal_form=crystal_form;self.deposit_date=deposit_date;self.pmid=pmid;

    def __repr__(self):
        return f'{self.pdb} {self.gene_name} {self.protein_name} {self.bound_protein} {self.bound_peptide} {self.mutation_status} {self.nucleotide_class} {self.drug_class} {self.bound_protein_class}\
            {self.homodimer_status} {self.switch1_cluster} {self.switch2_cluster} {self.conformational_state} {self.residue_range} {self.nucleotide}\
            {self.drug} {self.other_ligands} {self.experiment_type} {self.resolution} {self.rfactor} {self.crystal_form} {self.deposit_date} {self.pmid}'

def create_lists():
    pdb_set=set();conformation_set=set();bound_protein_set=set();bound_protein_class_set=set();nucleotide_class_set=set();drug_class_set=set()
    mutation_set=set();state1_form1_list=dict();confChainCount=dict()
            
    for items in Clusters.query.all():
        pdb_set.add(items.pdb[:4])
        conformation_set.add(items.conformational_state.strip())
        nucleotide_class_set.add(items.nucleotide_class.strip())
        drug_class_set.add(items.drug_class.strip())
        bound_protein_class_set.add(items.bound_protein_class.strip())   #Check if these variables can also have ',' in them
        
        if ',' in items.bound_protein:
            for protein in items.bound_protein.split(','):
                bound_protein_set.add(protein.strip())
        else:
            bound_protein_set.add(items.bound_protein.strip())
        
        if ',' in items.mutation_status:
            for mutation in items.mutation_status.split(','):
                mutation_set.add(mutation.strip())
        else:
            mutation_set.add(items.mutation_status.strip())
           
    count_entire_set=dict()
   
    for item in Clusters.query.all():
        if ',' in item.bound_protein:
            for proteins in item.bound_protein.split(','):
                state1_form1_list[f'{item.conformational_state};{proteins}']=Clusters.query.filter(Clusters.conformational_state.contains(item.conformational_state),Clusters.bound_protein.contains(proteins)).count()
        else:
            state1_form1_list[f'{item.conformational_state};{item.bound_protein}']=Clusters.query.filter(Clusters.conformational_state.contains(item.conformational_state),Clusters.bound_protein.contains(item.bound_protein)).count()
    
    for conf in conformation_set:
        confChainCount[conf]=Clusters.query.filter(Clusters.conformational_state.contains(conf)).count()

    return natsorted(pdb_set),natsorted(conformation_set),natsorted(drug_class_set),natsorted(bound_protein_class_set),natsorted(mutation_set),confChainCount,count_entire_set

def count_entries(Clusters):    #generalize this function
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

@app.route('/formSearchPDB', methods=['GET','POST'])
def formSearchPDB():
    (pdbListDb,conformation_set,drug_class_set,bound_protein_class_set,mutation_set,confChainCount,count_entire_set)=create_lists()

    if request.method=='POST':
        inputString=request.form['pdb_select']
        return redirect(url_for('uniqueQuery',queryname=inputString,settings='PDB'))
       
    
    return render_template('search.html',pdbListDb=pdbListDb,confChainCount=confChainCount,count_entire_set=count_entire_set,conformation_set=conformation_set,\
                           drug_class_set=drug_class_set,bound_protein_class_set=bound_protein_class_set,mutation_set=mutation_set)

@app.route('/formSearchConf', methods=['GET','POST'])
def formSearchConf():
    (pdbListDb,conformation_set,drug_class_set,bound_protein_class_set,mutation_set,confChainCount,count_entire_set)=create_lists()

    if request.method=='POST':
        inputString=request.form['conf_select']
        return redirect(url_for('uniqueQuery',queryname=inputString,settings='conformation'))
       
    
    return render_template('search.html',pdbListDb=pdbListDb,confChainCount=confChainCount,count_entire_set=count_entire_set,conformation_set=conformation_set,\
                           drug_class_set=drug_class_set,bound_protein_class_set=bound_protein_class_set,mutation_set=mutation_set)


@app.route('/formSearchMut', methods=['GET','POST'])
def formSearchMut():
    (pdbListDb,conformation_set,drug_class_set,bound_protein_class_set,mutation_set,confChainCount,count_entire_set)=create_lists()

    if request.method=='POST':
        inputString=request.form['mutation_select']
        return redirect(url_for('uniqueQuery',queryname=inputString,settings='mutation'))
       
    
    return render_template('search.html',pdbListDb=pdbListDb,confChainCount=confChainCount,count_entire_set=count_entire_set,conformation_set=conformation_set,\
                           drug_class_set=drug_class_set,bound_protein_class_set=bound_protein_class_set,mutation_set=mutation_set)

@app.route('/formSearchDrug', methods=['GET','POST'])
def formSearchDrug():
    (pdbListDb,conformation_set,drug_class_set,bound_protein_class_set,mutation_set,confChainCount,count_entire_set)=create_lists()

    if request.method=='POST':
        inputString=request.form['drug_class_select']
        return redirect(url_for('uniqueQuery',queryname=inputString,settings='drug_class'))
       
    
    return render_template('search.html',pdbListDb=pdbListDb,confChainCount=confChainCount,count_entire_set=count_entire_set,conformation_set=conformation_set,\
                           drug_class_set=drug_class_set,bound_protein_class_set=bound_protein_class_set,mutation_set=mutation_set)

@app.route('/formSearchBP', methods=['GET','POST'])
def formSearchBP():
    (pdbListDb,conformation_set,drug_class_set,bound_protein_class_set,mutation_set,confChainCount,count_entire_set)=create_lists()

    if request.method=='POST':
        inputString=request.form['bound_protein_class_select']
        return redirect(url_for('uniqueQuery',queryname=inputString,settings='bound_protein_class'))
       
    
    return render_template('search.html',pdbListDb=pdbListDb,confChainCount=confChainCount,count_entire_set=count_entire_set,conformation_set=conformation_set,\
                           drug_class_set=drug_class_set,bound_protein_class_set=bound_protein_class_set,mutation_set=mutation_set)


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
        protein_name=Clusters.query.filter(Clusters.pdb.contains(queryname)).first().protein_name
        experiment_type=Clusters.query.filter(Clusters.pdb.contains(queryname)).first().experiment_type
        resolution=Clusters.query.filter(Clusters.pdb.contains(queryname)).first().resolution
        rfactor=Clusters.query.filter(Clusters.pdb.contains(queryname)).first().rfactor
        crystal_form=Clusters.query.filter(Clusters.pdb.contains(queryname)).first().crystal_form
        homodimer_status=Clusters.query.filter(Clusters.pdb.contains(queryname)).first().homodimer_status
        residue_range=Clusters.query.filter(Clusters.pdb.contains(queryname)).first().residue_range
        pmid=Clusters.query.filter(Clusters.pdb.contains(queryname)).first().pmid
        
        return render_template('pdbs.html',queryname=queryname,pdb_list=pdb_list,protein_name=protein_name,\
                               experiment_type=experiment_type,resolution=resolution,rfactor=rfactor,crystal_form=crystal_form,\
                               homodimer_status=homodimer_status,residue_range=residue_range,pmid=pmid)

    if settings=='All':
        retrieve_str=dict();chain_count=dict();entry_count=dict()
        retrieve_str['All']=Clusters.query.all();chain_count['All']=Clusters.query.count()
        retrieve_str['HRAS']=Clusters.query.filter(Clusters.protein_name.contains('HRAS'));chain_count['HRAS']=Clusters.query.filter(Clusters.protein_name.contains('HRAS')).count()
        retrieve_str['KRAS']=Clusters.query.filter(Clusters.protein_name.contains('KRAS'));chain_count['KRAS']=Clusters.query.filter(Clusters.protein_name.contains('KRAS')).count()
        retrieve_str['NRAS']=Clusters.query.filter(Clusters.protein_name.contains('NRAS'));chain_count['NRAS']=Clusters.query.filter(Clusters.protein_name.contains('NRAS')).count()
        entry_count=count_entries(Clusters)
    
        return render_template('browse.html',retrieve_str=retrieve_str,chain_count=chain_count,entry_count=entry_count)

    if settings=='conformation':
        retrieve_str=dict();chain_count=dict();entry_count=dict()
        retrieve_str['All']=Clusters.query.filter(Clusters.conformational_state.contains(queryname));chain_count['All']=Clusters.query.filter(Clusters.conformational_state.contains(queryname)).count()
        retrieve_str['HRAS']=Clusters.query.filter(Clusters.protein_name.contains('HRAS'),Clusters.conformational_state.contains(queryname));chain_count['HRAS']=Clusters.query.filter(Clusters.conformational_state.contains(queryname),Clusters.protein_name.contains('HRAS')).count()
        retrieve_str['KRAS']=Clusters.query.filter(Clusters.protein_name.contains('KRAS'),Clusters.conformational_state.contains(queryname));chain_count['KRAS']=Clusters.query.filter(Clusters.conformational_state.contains(queryname),Clusters.protein_name.contains('KRAS')).count()
        retrieve_str['NRAS']=Clusters.query.filter(Clusters.protein_name.contains('NRAS'),Clusters.conformational_state.contains(queryname));chain_count['NRAS']=Clusters.query.filter(Clusters.conformational_state.contains(queryname),Clusters.protein_name.contains('NRAS')).count()
        entry_count=count_entries(Clusters)   #generalize this function
    
        return render_template('conformation.html',queryname=queryname,retrieve_str=retrieve_str,chain_count=chain_count,entry_count=entry_count)
    
    if settings=='bound_protein_class':
        retrieve_str=dict();chain_count=dict();entry_count=dict()
        retrieve_str['All']=Clusters.query.filter(Clusters.bound_protein_class.contains(queryname));chain_count['All']=Clusters.query.filter(Clusters.bound_protein_class.contains(queryname)).count()
        retrieve_str['HRAS']=Clusters.query.filter(Clusters.protein_name.contains('HRAS'),Clusters.bound_protein_class.contains(queryname));chain_count['HRAS']=Clusters.query.filter(Clusters.bound_protein_class.contains(queryname),Clusters.protein_name.contains('HRAS')).count()
        retrieve_str['KRAS']=Clusters.query.filter(Clusters.protein_name.contains('KRAS'),Clusters.bound_protein_class.contains(queryname));chain_count['KRAS']=Clusters.query.filter(Clusters.bound_protein_class.contains(queryname),Clusters.protein_name.contains('KRAS')).count()
        retrieve_str['NRAS']=Clusters.query.filter(Clusters.protein_name.contains('NRAS'),Clusters.bound_protein_class.contains(queryname));chain_count['NRAS']=Clusters.query.filter(Clusters.bound_protein_class.contains(queryname),Clusters.protein_name.contains('NRAS')).count()
        entry_count=count_entries(Clusters)   #generalize this function
    
        return render_template('bound_protein_class.html',queryname=queryname, retrieve_str=retrieve_str,chain_count=chain_count,entry_count=entry_count)

    if settings=='mutation':
        retrieve_str=dict();chain_count=dict();entry_count=dict()
        retrieve_str['All']=Clusters.query.filter(Clusters.mutation_status.contains(queryname));chain_count['All']=Clusters.query.filter(Clusters.mutation_status.contains(queryname)).count()
        retrieve_str['HRAS']=Clusters.query.filter(Clusters.protein_name.contains('HRAS'),Clusters.mutation_status.contains(queryname));chain_count['HRAS']=Clusters.query.filter(Clusters.mutation_status.contains(queryname),Clusters.protein_name.contains('HRAS')).count()
        retrieve_str['KRAS']=Clusters.query.filter(Clusters.protein_name.contains('KRAS'),Clusters.mutation_status.contains(queryname));chain_count['KRAS']=Clusters.query.filter(Clusters.mutation_status.contains(queryname),Clusters.protein_name.contains('KRAS')).count()
        retrieve_str['NRAS']=Clusters.query.filter(Clusters.protein_name.contains('NRAS'),Clusters.mutation_status.contains(queryname));chain_count['NRAS']=Clusters.query.filter(Clusters.mutation_status.contains(queryname),Clusters.protein_name.contains('NRAS')).count()
        entry_count=count_entries(Clusters)   #generalize this function
    
        return render_template('mutation.html',queryname=queryname,retrieve_str=retrieve_str,chain_count=chain_count,entry_count=entry_count)
    
    if settings=='drug_class':
        retrieve_str=dict();chain_count=dict();entry_count=dict()
        retrieve_str['All']=Clusters.query.filter(Clusters.drug_class.contains(queryname));chain_count['All']=Clusters.query.filter(Clusters.drug_class.contains(queryname)).count()
        retrieve_str['HRAS']=Clusters.query.filter(Clusters.protein_name.contains('HRAS'),Clusters.drug_class.contains(queryname));chain_count['HRAS']=Clusters.query.filter(Clusters.drug_class.contains(queryname),Clusters.protein_name.contains('HRAS')).count()
        retrieve_str['KRAS']=Clusters.query.filter(Clusters.protein_name.contains('KRAS'),Clusters.drug_class.contains(queryname));chain_count['KRAS']=Clusters.query.filter(Clusters.drug_class.contains(queryname),Clusters.protein_name.contains('KRAS')).count()
        retrieve_str['NRAS']=Clusters.query.filter(Clusters.protein_name.contains('NRAS'),Clusters.drug_class.contains(queryname));chain_count['NRAS']=Clusters.query.filter(Clusters.drug_class.contains(queryname),Clusters.protein_name.contains('NRAS')).count()
        entry_count=count_entries(Clusters)   #generalize this function
    
        return render_template('drug_class.html',queryname=queryname, retrieve_str=retrieve_str,chain_count=chain_count,entry_count=entry_count)



        

if __name__ == '__main__':
    app.run(debug=True)
