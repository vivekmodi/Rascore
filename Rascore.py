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
    mutation_set=set();nucleotide_class_count_dict=dict()
            
    for items in Clusters.query.all():
        pdb_set.add(items.pdb[:4])
        conformation_set.add(items.conformational_state.strip())
        nucleotide_class_set.add(items.nucleotide_class.strip('"'))
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
           
    conformation_set.add('All') ;mutation_set.add('All');drug_class_set.add('All');bound_protein_class_set.add('All')
  
   
    for nucleo in nucleotide_class_set:    #For count of chains per nucleotide class in Search menu
        for conf in conformation_set:
            if conf=='All':
                nucleotide_class_count_dict[nucleo,conf]=Clusters.query.filter(Clusters.nucleotide_class.contains(nucleo)).count()
            else:    
                nucleotide_class_count_dict[nucleo,conf]=Clusters.query.filter(Clusters.conformational_state.contains(conf),Clusters.nucleotide_class.contains(nucleo)).count()
        
        for mut in mutation_set:
            if mut=='All':
                nucleotide_class_count_dict[nucleo,conf]=Clusters.query.filter(Clusters.nucleotide_class.contains(nucleo)).count()
            else:
                nucleotide_class_count_dict[nucleo,mut]=Clusters.query.filter(Clusters.mutation_status.contains(mut),Clusters.nucleotide_class.contains(nucleo)).count()
        
        for drug in drug_class_set:
            if drug=='All':                                
                nucleotide_class_count_dict[nucleo,conf]=Clusters.query.filter(Clusters.nucleotide_class.contains(nucleo)).count()
            else:
                nucleotide_class_count_dict[nucleo,drug]=Clusters.query.filter(Clusters.drug_class.contains(drug),Clusters.nucleotide_class.contains(nucleo)).count()
        
        for bp in bound_protein_class_set:
            if bp=='All':
                nucleotide_class_count_dict[nucleo,conf]=Clusters.query.filter(Clusters.nucleotide_class.contains(nucleo)).count()
            else:
                nucleotide_class_count_dict[nucleo,bp]=Clusters.query.filter(Clusters.bound_protein_class.contains(bp),Clusters.nucleotide_class.contains(nucleo)).count()
    
    
    return natsorted(pdb_set),natsorted(conformation_set),natsorted(drug_class_set),natsorted(bound_protein_class_set),natsorted(mutation_set,key=lambda x: x.replace('All', 'A00')),nucleotide_class_set,nucleotide_class_count_dict

def count_entries(sublist):      #Count number of PDB entries
    count=set()
    for items in sublist:
        count.add(items.pdb[0:4])
    
    return len(count)
    
    
def write_text_file(sublist,tsvFile):
    fhandle_textFile=open(f'{pwd}/static/{tsvFile}','w')
    
    fhandle_textFile.write('PDB ID\tGene Name\tProtein Name\tBound Protein\tBound Peptide\tMutation Status\tNucleotide Class\tDrug Class\tBound Protein Class\tHomodimer Status\tSwitch 1 Cluster\tSwitch 2 Cluster\tConformational State\tResidue Range\tNucleotide\tDrug\tOther Ligands\tExperiment Type\tResolution\tR-Factor\tCrystal Form\tDeposit Date\tPMID\n')
    for item in sublist:
        fhandle_textFile.write(f'{item.pdb}\t{item.gene_name}\t{item.protein_name}\t{item.bound_protein}\t{item.bound_peptide}\t{item.mutation_status}\t{item.nucleotide_class}\t{item.drug_class}\t{item.bound_protein_class}\
            {item.homodimer_status}\t{item.switch1_cluster}\t{item.switch2_cluster}\t{item.conformational_state}\t{item.residue_range}\t{item.nucleotide}\
            {item.drug}\t{item.other_ligands}\t{item.experiment_type}\t{item.resolution}\t{item.rfactor}\t{item.crystal_form}\t{item.deposit_date}\t{item.pmid}\n')

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
    (pdbListDb,conformation_set,drug_class_set,bound_protein_class_set,mutation_set,nucleotide_class_set,nucleotide_class_count_dict)=create_lists()

    if request.method=='POST':
        inputString=request.form['pdb_select']
        return redirect(url_for('uniqueQuery',queryname=inputString,settings='PDB'))
       
    
    return render_template('search.html',pdbListDb=pdbListDb,conformation_set=conformation_set,\
                           drug_class_set=drug_class_set,bound_protein_class_set=bound_protein_class_set,mutation_set=mutation_set,nucleotide_class_set=nucleotide_class_set,nucleotide_class_count_dict=nucleotide_class_count_dict)

@app.route('/formSearchConf', methods=['GET','POST'])
def formSearchConf():
    (pdbListDb,conformation_set,drug_class_set,bound_protein_class_set,mutation_set,nucleotide_class_set,nucleotide_class_count_dict)=create_lists()

    if request.method=='POST':
        inputString=request.form['conf_select']
        if inputString=='All':
            return redirect(url_for('uniqueQuery',settings='All',queryname='All'))
        else:
            return redirect(url_for('uniqueQuery',queryname=inputString,settings='conformation'))
       
    
    return render_template('search.html',pdbListDb=pdbListDb,conformation_set=conformation_set,\
                           drug_class_set=drug_class_set,bound_protein_class_set=bound_protein_class_set,mutation_set=mutation_set,nucleotide_class_set=nucleotide_class_set,nucleotide_class_count_dict=nucleotide_class_count_dict)


@app.route('/formSearchMut', methods=['GET','POST'])
def formSearchMut():
    (pdbListDb,conformation_set,drug_class_set,bound_protein_class_set,mutation_set,nucleotide_class_set,nucleotide_class_count_dict)=create_lists()

    if request.method=='POST':
        inputString=request.form['mutation_select']
        if inputString=='All':
            return redirect(url_for('uniqueQuery',settings='All',queryname='All'))
        else:
            return redirect(url_for('uniqueQuery',queryname=inputString,settings='mutation'))
       
    
    return render_template('search.html',pdbListDb=pdbListDb,conformation_set=conformation_set,\
                           drug_class_set=drug_class_set,bound_protein_class_set=bound_protein_class_set,mutation_set=mutation_set,nucleotide_class_set=nucleotide_class_set,nucleotide_class_count_dict=nucleotide_class_count_dict)

@app.route('/formSearchDrug', methods=['GET','POST'])
def formSearchDrug():
    (pdbListDb,conformation_set,drug_class_set,bound_protein_class_set,mutation_set,nucleotide_class_set,nucleotide_class_count_dict)=create_lists()

    if request.method=='POST':
        inputString=request.form['drug_class_select']
        if inputString=='All':
            return redirect(url_for('uniqueQuery',settings='All',queryname='All'))
        else:
            return redirect(url_for('uniqueQuery',queryname=inputString,settings='drug_class'))
       
    
    return render_template('search.html',pdbListDb=pdbListDb,conformation_set=conformation_set,\
                           drug_class_set=drug_class_set,bound_protein_class_set=bound_protein_class_set,mutation_set=mutation_set,nucleotide_class_set=nucleotide_class_set,nucleotide_class_count_dict=nucleotide_class_count_dict)

@app.route('/formSearchBP', methods=['GET','POST'])
def formSearchBP():
    (pdbListDb,conformation_set,drug_class_set,bound_protein_class_set,mutation_set,confChainCount,count_entire_set,nucleotide_class_set,nucleotide_class_count_dict)=create_lists()

    if request.method=='POST':
        inputString=request.form['bound_protein_class_select']
        if inputString=='All':
            return redirect(url_for('uniqueQuery',settings='All',queryname='All'))
        else:
            return redirect(url_for('uniqueQuery',queryname=inputString,settings='bound_protein_class'))
       
    
    return render_template('search.html',pdbListDb=pdbListDb,conformation_set=conformation_set,\
                           drug_class_set=drug_class_set,bound_protein_class_set=bound_protein_class_set,mutation_set=mutation_set,nucleotide_class_set=nucleotide_class_set,nucleotide_class_count_dict=nucleotide_class_count_dict)


@app.route('/webserver', methods=['GET','POST'])
def webserver():
    return render_template('webserver.html')

@app.route('/download')
def download():
    return render_template('download.html')

@app.route('/contact')
def contact():
    return render_template('contact.html')

@app.route('/dunbrackLab')
def dunbrackLab():
    return redirect("http://dunbrack.fccc.edu")


@app.route('/<settings>/<queryname>')
def uniqueQuery(settings,queryname):
    tsvFile=dict();pymolSessionRe=dict();pymolScriptRe=dict();pymolSession=dict();pymolScript=dict()
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
        
        tsvFile=f'downloads/text-files/{queryname}.tsv'
        write_text_file(pdb_list,tsvFile)
        
        return render_template('pdbs.html',queryname=queryname,pdb_list=pdb_list,protein_name=protein_name,\
                               experiment_type=experiment_type,resolution=resolution,rfactor=rfactor,crystal_form=crystal_form,\
                               homodimer_status=homodimer_status,residue_range=residue_range,pmid=pmid,tsvFile=tsvFile,\
                               pymolSession=pymolSession,pymolScript=pymolScript)

    if settings=='All':
        retrieve_str=dict();chain_count=dict();entry_count=dict()
        retrieve_str['All']=Clusters.query.all();chain_count['All']=Clusters.query.count()
        entry_count['All']=count_entries(retrieve_str['All'])
        tsvFile['All']='downloads/text-files/All_RAS.tsv'
        write_text_file(retrieve_str['All'],tsvFile['All'])
        pymolSession['All']='downloads/pymolSessions/All_RAS.pse.zip';pymolScript['All']='downloads/pymolScripts/All_RAS.pml.zip'
        pymolSessionRe['All']='downloads/pymolSessionsRe/Repr_RAS.pse.zip';pymolScriptRe['All']='downloads/pymolScriptsRe/Repr_RAS.pml.zip'
        
        for gene in ('HRAS','KRAS','NRAS'):
            retrieve_str[gene]=Clusters.query.filter(Clusters.gene_name.contains(gene));chain_count[gene]=Clusters.query.filter(Clusters.gene_name.contains(gene)).count()
            entry_count[gene]=count_entries(retrieve_str[gene])
            tsvFile[gene]=f'downloads/text-files/All_{gene}.tsv'
            write_text_file(retrieve_str[gene],tsvFile[gene])
            pymolSession[gene]=f'downloads/pymolSessions/All_{gene}.pse.zip';pymolScript[gene]=f'downloads/pymolScripts/All_{gene}.pml.zip'
            pymolSessionRe[gene]=f'downloads/pymolSessionsRe/Repr_{gene}.pse.zip';pymolScriptRe[gene]=f'downloads/pymolScriptsRe/Repr_{gene}.pml.zip'
        
        return render_template('browse.html',retrieve_str=retrieve_str,chain_count=chain_count,entry_count=entry_count,tsvFile=tsvFile,\
                               pymolSessionRe=pymolSessionRe,pymolScriptRe=pymolScriptRe,pymolSession=pymolSession,pymolScript=pymolScript)

    if settings=='conformation':
        retrieve_str=dict();chain_count=dict();entry_count=dict()
        retrieve_str['All']=Clusters.query.filter(Clusters.conformational_state.contains(queryname));chain_count['All']=Clusters.query.filter(Clusters.conformational_state.contains(queryname)).count()
        entry_count['All']=count_entries(retrieve_str['All'])
        tsvFile['All']=f'downloads/text-files/All_{queryname}.tsv'
        write_text_file(retrieve_str['All'],tsvFile['All'])
        pymolSession['All']=f'downloads/pymolSessions/All_{queryname}.pse.zip';pymolScript['All']=f'downloads/pymolScripts/All_{queryname}.pml.zip'
        pymolSessionRe['All']=f'downloads/pymolSessionsRe/Repr_{queryname}.pse.zip';pymolScriptRe['All']=f'downloads/pymolScriptsRe/Repr_{queryname}.pml.zip'
        
        for gene in ('HRAS','KRAS','NRAS'):
            retrieve_str[gene]=Clusters.query.filter(Clusters.gene_name.contains(gene),Clusters.conformational_state.contains(queryname));chain_count[gene]=Clusters.query.filter(Clusters.conformational_state.contains(queryname),Clusters.gene_name.contains(gene)).count()
            retrieve_str[gene]=Clusters.query.filter(Clusters.gene_name.contains(gene),Clusters.conformational_state.contains(queryname));chain_count[gene]=Clusters.query.filter(Clusters.conformational_state.contains(queryname),Clusters.gene_name.contains(gene)).count()
            retrieve_str[gene]=Clusters.query.filter(Clusters.gene_name.contains(gene),Clusters.conformational_state.contains(queryname));chain_count[gene]=Clusters.query.filter(Clusters.conformational_state.contains(queryname),Clusters.gene_name.contains(gene)).count()
            entry_count[gene]=count_entries(retrieve_str[gene])
            tsvFile[gene]=f'downloads/text-files/{gene}_{queryname}.tsv'
            write_text_file(retrieve_str[gene],tsvFile[gene])
            pymolSession[gene]=f'downloads/pymolSessions/All_{gene}_{queryname}.pse.zip';pymolScript[gene]=f'downloads/pymolScripts/All_{gene}_{queryname}.pml.zip'
            pymolSessionRe[gene]=f'downloads/pymolSessionsRe/Repr_{gene}_{queryname}.pse.zip';pymolScriptRe[gene]=f'downloads/pymolScriptsRe/Repr_{gene}_{queryname}.pml.zip'
    
        return render_template('conformation.html',queryname=queryname,retrieve_str=retrieve_str,chain_count=chain_count,entry_count=entry_count,tsvFile=tsvFile,\
                               pymolSessionRe=pymolSessionRe,pymolScriptRe=pymolScriptRe,pymolSession=pymolSession,pymolScript=pymolScript)
    
    if settings=='bound_protein_class':
        retrieve_str=dict();chain_count=dict();entry_count=dict()
        retrieve_str['All']=Clusters.query.filter(Clusters.bound_protein_class.contains(queryname));chain_count['All']=Clusters.query.filter(Clusters.bound_protein_class.contains(queryname)).count()
        entry_count['All']=count_entries(retrieve_str['All'])
        tsvFile['All']=f'downloads/text-files/All_{queryname}.tsv'
        write_text_file(retrieve_str['All'],tsvFile['All'])
        pymolSession['All']=f'downloads/pymolSessions/All_{queryname}.pse.zip';pymolScript['All']=f'downloads/pymolScripts/All_{queryname}.pml.zip'
        pymolSessionRe['All']=f'downloads/pymolSessionsRe/Repr_{queryname}.pse.zip';pymolScriptRe['All']=f'downloads/pymolScriptsRe/Repr_{queryname}.pml.zip'
        
        for gene in ('HRAS','KRAS','NRAS'):
            retrieve_str[gene]=Clusters.query.filter(Clusters.gene_name.contains(gene),Clusters.bound_protein_class.contains(queryname));chain_count[gene]=Clusters.query.filter(Clusters.bound_protein_class.contains(queryname),Clusters.gene_name.contains(gene)).count()
            entry_count[gene]=count_entries(retrieve_str[gene])
            tsvFile[gene]=f'downloads/text-files/{gene}_{queryname}.tsv'
            write_text_file(retrieve_str[gene],tsvFile[gene])
            pymolSession[gene]=f'downloads/pymolSessions/All_{gene}_{queryname}.pse.zip';pymolScript[gene]=f'downloads/pymolScripts/All_{gene}_{queryname}.pml.zip'
            pymolSessionRe[gene]=f'downloads/pymolSessionsRe/Repr_{gene}_{queryname}.pse.zip';pymolScriptRe[gene]=f'downloads/pymolScriptsRe/Repr_{gene}_{queryname}.pml.zip'
            
    
        return render_template('bound_protein_class.html',queryname=queryname, retrieve_str=retrieve_str,chain_count=chain_count,entry_count=entry_count,tsvFile=tsvFile,\
                               pymolSessionRe=pymolSessionRe,pymolScriptRe=pymolScriptRe,pymolSession=pymolSession,pymolScript=pymolScript)

    if settings=='mutation':
        retrieve_str=dict();chain_count=dict();entry_count=dict()
        retrieve_str['All']=Clusters.query.filter(Clusters.mutation_status.contains(queryname));chain_count['All']=Clusters.query.filter(Clusters.mutation_status.contains(queryname)).count()
        entry_count['All']=count_entries(retrieve_str['All'])
        tsvFile['All']=f'downloads/text-files/All_{queryname}.tsv'
        write_text_file(retrieve_str['All'],tsvFile['All'])
        pymolSession['All']=f'downloads/pymolSessions/All_{queryname}.pse.zip';pymolScript['All']=f'downloads/pymolScripts/All_{queryname}.pml.zip'
        pymolSessionRe['All']=f'downloads/pymolSessionsRe/Repr_{queryname}.pse.zip';pymolScriptRe['All']=f'downloads/pymolScriptsRe/Repr_{queryname}.pml.zip'
        
        for gene in ('HRAS','KRAS','NRAS'):
            retrieve_str[gene]=Clusters.query.filter(Clusters.gene_name.contains(gene),Clusters.mutation_status.contains(queryname));chain_count[gene]=Clusters.query.filter(Clusters.mutation_status.contains(queryname),Clusters.gene_name.contains(gene)).count()
            entry_count[gene]=count_entries(retrieve_str[gene])
            tsvFile[gene]=f'downloads/text-files/{gene}_{queryname}.tsv'
            write_text_file(retrieve_str[gene],tsvFile[gene])
            pymolSession[gene]=f'downloads/pymolSessions/All_{gene}_{queryname}.pse.zip';pymolScript[gene]=f'downloads/pymolScripts/All_{gene}_{queryname}.pml.zip'
            pymolSessionRe[gene]=f'downloads/pymolSessionsRe/Repr_{gene}_{queryname}.pse.zip';pymolScriptRe[gene]=f'downloads/pymolScriptsRe/Repr_{gene}_{queryname}.pml.zip'
    
        return render_template('mutation.html',queryname=queryname,retrieve_str=retrieve_str,chain_count=chain_count,entry_count=entry_count,tsvFile=tsvFile,\
                               pymolSessionRe=pymolSessionRe,pymolScriptRe=pymolScriptRe,pymolSession=pymolSession,pymolScript=pymolScript)
    
    if settings=='drug_class':
        retrieve_str=dict();chain_count=dict();entry_count=dict()
        retrieve_str['All']=Clusters.query.filter(Clusters.drug_class.contains(queryname));chain_count['All']=Clusters.query.filter(Clusters.drug_class.contains(queryname)).count()
        entry_count['All']=count_entries(retrieve_str['All'])
        tsvFile['All']=f'downloads/text-files/All_{queryname}.tsv'
        write_text_file(retrieve_str['All'],tsvFile['All'])
        pymolSession['All']=f'downloads/pymolSessions/All_{queryname}.pse.zip';pymolScript['All']=f'downloads/pymolScripts/All_{queryname}.pml.zip'
        pymolSessionRe['All']=f'downloads/pymolSessionsRe/Repr_{queryname}.pse.zip';pymolScriptRe['All']=f'downloads/pymolScriptsRe/Repr_{queryname}.pml.zip'
        
        for gene in ('HRAS','KRAS','NRAS'):
            retrieve_str[gene]=Clusters.query.filter(Clusters.gene_name.contains(gene),Clusters.drug_class.contains(queryname));chain_count[gene]=Clusters.query.filter(Clusters.drug_class.contains(queryname),Clusters.gene_name.contains(gene)).count()
            entry_count[gene]=count_entries(retrieve_str[gene])
            tsvFile[gene]=f'downloads/text-files/{gene}_{queryname}.tsv'
            write_text_file(retrieve_str[gene],tsvFile[gene])
            pymolSession[gene]=f'downloads/pymolSessions/All_{gene}_{queryname}.pse.zip';pymolScript[gene]=f'downloads/pymolScripts/All_{gene}_{queryname}.pml.zip'
            pymolSessionRe[gene]=f'downloads/pymolSessionsRe/Repr_{gene}_{queryname}.pse.zip';pymolScriptRe[gene]=f'downloads/pymolScriptsRe/Repr_{gene}_{queryname}.pml.zip'
    
        return render_template('drug_class.html',queryname=queryname, retrieve_str=retrieve_str,chain_count=chain_count,entry_count=entry_count,tsvFile=tsvFile,\
                               pymolSessionRe=pymolSessionRe,pymolScriptRe=pymolScriptRe,pymolSession=pymolSession,pymolScript=pymolScript)



        

if __name__ == '__main__':
    app.run(debug=True)
