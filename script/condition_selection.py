# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 20:15:43 2023

@author: Yonglan Liu
"""
import pandas as pd

def gene(database, genes):
    """filter case by gene name (Hugo_Symbol)"""
    
    database_r = pd.DataFrame()
    for gene in genes:
        database_r = pd.concat([database_r, database[database["Hugo_Symbol"] == gene]], axis = 0)
    return database_r

def varclass(database, varclass):
    
    """filter case by Variant Classification:
        Missense_Mutation
        Silent
        Nonsense_Mutation
        Intron
        Nonstop_Mutation
        Splice_Region
        Intron
        Frame_Shift_Del
        Splice_Site
        Frame_Shift_Ins
        In_Frame_Del
        In_Frame_Ins
        5'UTR 
        Translation_Start_Site
        3'UTR
        RNA
        Targeted_Region
    """
    
    database_r = pd.DataFrame()
    
    for vc in varclass:
        
        database_r = pd.concat([database_r, database[database["Variant_Classification"] == vc]], axis = 0)
                                
    return database_r
    
def patient(database, patientID):
    """filter case by Patient ID"""
    database = database[database["Patient_ID"] == patientID]
    return database

def chromosome(database, chromosome):
    """
    filter case by chromosome
    """
    database = database[database["Chromosome"] == chromosome]
    return database

def vartype(database, vartype):
    """ 
    filter case by variation type:
    SNP: Single nucleotide polymorphism -- a substitution in one nucleotide
    DNP: Double nucleotide polymorphism -- a substitution in two consecutive nucleotides
    TNP: Triple nucleotide polymorphism -- a substitution in three consecutive nucleotides
    ONP: Oligo-nucleotide polymorphism -- a substitution in more than three consecutive nucleotides
    DEL: Deletion -- the removal of nucleotides
    INS: Insertion -- the addition of nucleotides
    
    SNP, DNP, TNP, and ONP are the missense mutations
    """
    database  = database[database["Variant_Type"] == vartype]
    return database

def sample(database, sampleID):
    """
    filter case by tumor barcode
    """
    database = database[database["Sample_ID"] == sampleID]
    return database

def vaf(database, vaf_cutoff=0.125):
    
    """
    filter case with Variant allele frequency (VAF) less than a cutoff
    VAF = t_alt_count/(t_ref_count + t_alt_count)
    
    VAF = (number of reads with variant / total number of reads)
    """
    database = database[database["VAF"] > vaf_cutoff]
    return database

def freq_filter(database, freq=3):
    """filter by mutation frequency"""
    db_freq = database.groupby(["HGVSp_Short"], as_index = False).size()
    db_freq = db_freq[db_freq["size"] > 3]
    protein_mut = db_freq["HGVSp_Short"].values
    database = database[database["HGVSp_Short"].isin(protein_mut)]
    return database

def tumor_filter(database, tumorIDs):
    database = database[database["Tumor_ID"].isin(tumorIDs)]
    return database