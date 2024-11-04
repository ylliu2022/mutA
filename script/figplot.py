# -*- coding: utf-8 -*-
"""
Created on Sat Jun  3 10:32:16 2023

@author: fyong
"""

# access to Numpy and pandas functionality
import pandas as pd


def pp_lollipop(databases):
    
    """
    This is for draw a protein frequency figure using ProteinPaint: https://proteinpaint.stjude.org/
    The required columns: Hugo_Symbol, isoform, Chromosome, start, aachange
    Hugo_Symbol: gene name
    isoform: NCBI gene annotation (RefSeq), like NM_????. This can be obtained in bioportal on the top-right position
    Chromosome: can be keep 
    start: gene start_position
    aachange: Protein Change
    """
    # create a empty DataFrame
    data_df = pd.DataFrame()

    # required columns, except for the four required columns. Hugo_Symbol and isoform are needed
    #cols = ["Hugo_Symbol", "HGVSp_Short", "Variant_Classification", "Chromosome", "Start_Position", "RefSeq"] # columns to be selected
    cols = ["HGVSp_Short", "Protein_position", "Variant_Classification"]

    # load downloaded database and combine all databases
    for db in databases:
        	data_df = db[cols]

    # modify required data
    data_df["Variant_Classification"] = "M"
    return data_df
