# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 20:02:23 2023

@author: Yonglan Liu
"""
import pandas as pd

source_file = pd.read_excel("../database/oncotree/oncotree_tumor_tissue.xlsx")

all_oncotree_code = source_file["oncotree_code"].values

def code_to_disease(code):
    
    """map oncotree code to tumor name"""
    if (code in all_oncotree_code) == True:
        tumor = source_file[source_file["oncotree_code"] == code].values[0,3]
    else:
        tumor = "Unknown"
        
    return tumor
        
def code_to_tissue(code):
    
    """map oncotree code to tissue"""
    if (code in all_oncotree_code) == True:
        tissue = source_file[source_file["oncotree_code"] == code].values[0,4]
    else:
        tissue = "Unknown"
    return tissue


