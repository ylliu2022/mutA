# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 20:15:43 2023

@author: Yonglan Liu
"""

import pandas as pd
import numpy as np


GENIE_COLUMN_NAMES = ["Hugo_Symbol", 
                  "Chromosome",
                  "Start_Position",
                  "End_Position",
                  "Variant_Classification",
                  "Variant_Type",
                  "Reference_Allele",	
                  "Tumor_Seq_Allele1",
                  "Tumor_Seq_Allele2",
                  "Tumor_Sample_Barcode",
                  "Mutation_Status",
                  "t_depth",
                  "t_ref_count",
                  "t_alt_count",
                  "n_depth",
                  "n_ref_count",	
                  "n_alt_count",
                  "HGVSc",
                  "HGVSp_Short",
                  "Transcript_ID",
                  "RefSeq",
                  "Protein_position",
                  "SIFT_Prediction",
                  "SIFT_Score",
                  "Polyphen_Prediction",
                  "Polyphen_Score"]

def genieM(gene):
    
    """load genie mutation database"""
    try:
        db = pd.read_pickle("../database/genie_refined.pkl")
    except:
        db = pd.read_csv("../database/genie/data_mutations_extended.txt", 
                     sep = "\t", 
                     usecols= GENIE_COLUMN_NAMES)
        db.to_pickle("../database/genie_refined.pkl")
    
    # select data for target gene
    db = db[db["Hugo_Symbol"] == gene]
    db=db[db["Variant_Classification"]=="Missense_Mutation"]
    
    # data clean
    db["HGVSp_Short"] = db["HGVSp_Short"].apply(lambda x: x[2:] if x is not np.nan else np.nan)
    db.rename({"Tumor_Sample_Barcode": "Sample_ID"}, axis=1, inplace=True)
    db["VAF"] = db["t_alt_count"]/db["t_depth"]
    return db


def genieClinicS():
    
    """load clinical sample data from Genie database"""
    try:
        db = pd.read_pickle("../database/genie_clinc.pkl")
    except:
        db = pd.read_csv("../database/genie/data_clinical_sample.txt", sep = "\t")
        db.to_pickle("../database/genie_clinc.pkl")
    
    db.columns = db.iloc[3]
    db = db.iloc[4:, ::]
    
    db = db.drop(["SEQ_ASSAY_ID", 
                  "CANCER_TYPE_DETAILED",
                  "SAMPLE_TYPE_DETAILED",
                  "SAMPLE_TYPE",
                  "CANCER_TYPE"], axis=1)
    
    db.rename({"PATIENT_ID":"Patient_ID",
               "SAMPLE_ID": "Sample_ID",
               "AGE_AT_SEQ_REPORT": "Age",
               "ONCOTREE_CODE": "Oncotree_Code"
               }, axis=1, inplace=True)
    return db


###################################################################
#                         Load TCGA database                      #
###################################################################

MAF_COLUMN_NAMES = ["Hugo_Symbol", 
                  "Chromosome",
                  "Start_Position",
                  "End_Position",
                  "Variant_Classification",
                  "Variant_Type",
                  "Reference_Allele",	
                  "Tumor_Seq_Allele1",
                  "Tumor_Seq_Allele2",
                  "Tumor_Sample_Barcode",
                  "Mutation_Status",
                  "t_depth",
                  "t_ref_count",
                  "t_alt_count",
                  "n_depth",
                  "n_ref_count",	
                  "n_alt_count",
                  "HGVSc",
                  "HGVSp_Short",
                  "Transcript_ID",
                  "RefSeq",
                  "Protein_position",
                  "SIFT",
                  "PolyPhen"]
    
def tcgaMC3MAF(gene):

    """load genie mutation database"""
    try:
        db = pd.read_pickle("../database/tcga_refined.pkl")
    except:
        db = pd.read_csv("../database/tcga/mc3.v0.2.8.PUBLIC.maf", 
                    sep = "\t", 
                    usecols= MAF_COLUMN_NAMES)
        db.to_pickle("../database/tcga_refined.pkl")
    # select data for target gene
    db = db[db["Hugo_Symbol"] == gene]
    db=db[db["Variant_Classification"]=="Missense_Mutation"]
    
    # data clean
    db["HGVSp_Short"] = db["HGVSp_Short"].apply(lambda x: x[2:] if x is not np.nan else np.nan)
    db.rename({"Tumor_Sample_Barcode": "Sample_ID"}, axis=1, inplace=True)
    db["VAF"] = db["t_alt_count"]/db["t_depth"]
    
    # add patient ID
    db["Patient_ID"] = db["Sample_ID"].str.split("-")
    db["Patient_ID"] = db["Patient_ID"].apply(lambda x: x[0] + "-" + x[1] + "-" + x[2])
    
    # add columns that is consistant with GENIE
    db["SIFT"] = db["SIFT"].str.replace("(", ",").str.replace(")", "").str.split(",")
    db["SIFT_Prediction"] = db["SIFT"].apply(lambda x: x[0] if len(x)==2 else np.nan)
    db["SIFT_Score"] = db["SIFT"].apply(lambda x: float(x[1]) if len(x) ==2 else np.nan)
    db["PolyPhen"] = db["PolyPhen"].str.replace("(", ",").str.replace(")", "").str.split(",")
    db["Polyphen_Prediction"] = db["PolyPhen"].apply(lambda x: x[0] if len(x)==2 else np.nan)
    db["Polyphen_Score"] = db["PolyPhen"].apply(lambda x: float(x[1]) if len(x) ==2 else np.nan)
    db = db.drop(["SIFT", "PolyPhen"], axis=1)
    
    return db


def tcgaClinicP():
    """load clinical patient data from Genie database"""
    try:
        db = pd.read_pickle("../database/tcga_clinic.pkl")
    except:
        db = pd.read_csv("../database/tcga/tcgaClinic.txt", sep = "\t")
        db.to_pickle("../database/tcga_clinic.pkl")
    return db


def tumor(db, tumorIDs):
    
    if db == "GENIE":
        """load genie mutation database"""
        
        db = pd.read_csv("../database/genie/data_mutations_extended.txt", 
                         sep = "\t", 
                         usecols= GENIE_COLUMN_NAMES)
    elif db == "TCGA":
        
        """load tcga mutation database"""
        db = pd.read_csv("../database/tcga/mc3.v0.2.8.PUBLIC.maf", 
                             sep = "\t", 
                             usecols= MAF_COLUMN_NAMES)
    
    # data clean
    db["HGVSp_Short"] = db["HGVSp_Short"].apply(lambda x: x[2:] if x is not np.nan else np.nan)
    db.rename({"Tumor_Sample_Barcode": "Sample_ID"}, axis=1, inplace=True)
    db["VAF"] = db["t_alt_count"]/db["t_depth"]
    db = db[db["Sample_ID"].isin(tumorIDs)]
        
    return db