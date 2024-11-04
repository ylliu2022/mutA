# -*- coding: utf-8 -*-
"""
Yonglan Liu

This is a temporary script file.
"""
import pandas as pd
import database
import condition_selection
import figplot
import mapping
import numpy as np

OUT_PATH = "../mtor_out"

# load clinic data
genieClinic = database.genieClinicS()
tcgaClinic = database.tcgaClinicP()

# select target gene data
gene = "MTOR"
genieM = database.genieM(gene)
tcgaM = database.tcgaMC3MAF(gene)

# select mutation type data, it is a list
mutType = ["Missense_Mutation"]
genie_missense =condition_selection.varclass(genieM, mutType)
tcga_missense = condition_selection.varclass(tcgaM, mutType)

# merge mutation data and clinic data
genie_merged = pd.merge(genie_missense, genieClinic, on = "Sample_ID", how = "left")
tcga_merged = pd.merge(tcga_missense, tcgaClinic, on = "Patient_ID", how = "left")

# merge tcga and genie
merged = pd.concat([tcga_merged, genie_merged], axis=0)

# filter out samples with mutilple substitutions
merged_updated = merged[~(merged["HGVSp_Short"].str.contains("delins"))]

# drop out number of reads with variant is 0 or NA
merged_updated = merged_updated[(merged_updated["t_ref_count"].notna()) &
                                           (merged_updated["t_ref_count"] != 0)]

# filter out mutations with VAF less than 0.125 to ensure >25% reads contain mutations
merged_updated = condition_selection.vaf(merged_updated, 0.125)

# filter out mutations with frequence less than or equal to 3
merged_updated = condition_selection.freq_filter(merged_updated, 3)

# save file for Lollipop plot
databases = [merged_updated] # databases to be analyzed 
proteinpaint_data = figplot.pp_lollipop(databases)

# save data
proteinpaint_data.to_csv(OUT_PATH + "/" + "ProteinPaint.txt", sep = "\t", index=False)

# mapping tissue and disease
merged_updated["Tissue"] = merged_updated["Oncotree_Code"].apply(mapping.code_to_tissue)
merged_updated["Disease"] = merged_updated["Oncotree_Code"].apply(mapping.code_to_disease)
merged_updated.to_excel("../mtor_out/mutation_data.xlsx", index = False)

# create a table for mapping tissue analysis
by_tissue = merged_updated.groupby(["Tissue"], as_index=False).size()
by_tissue.rename({"size": "Count"}, axis=1, inplace=True)
#by_tissue = by_tissue[by_tissue["Count"] > 5]
by_tissue.to_excel("../mtor_out/mutation_count_by_tissue.xlsx", index = False)

by_tissue_mut = merged_updated.groupby(["Tissue", "HGVSp_Short"], as_index=False).size()
by_tissue_mut.to_excel("../mtor_out/tissue_mut.xlsx", index = False)

unique_mutations = merged_updated["HGVSp_Short"].unique()
unique_mutation_df = np.zero(76)
