#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import codecs


# Gene list of reduced features

# In[2]:


overall_27k = pd.read_csv("data/genes_long_format/sure/27k_All.csv")
#print(ANOVA_RF_df.head)
overall_27k_genes = set(list(overall_27k["Gene_Symbol"]))

overall_450k = pd.read_csv("data/genes_long_format/sure/450k_All.csv")
#print(ANOVA_RF_df.head)
overall_450k_genes = set(list(overall_450k["Gene_Symbol"]))

ANOVA_RF_df_450s = pd.read_csv("data/genes_long_format/sure/450_sim_AnovaRF_markerData_wSMOTE.csv")
#print(ANOVA_RF_df.head)
ANOVA_RF_genes_450s = set(list(ANOVA_RF_df_450s["Gene_Symbol"]))

ANOVA_RF_df_450wos = pd.read_csv("data/genes_long_format/sure/450_sim_AnovaRF_markerData_woSMOTE.csv")
#print(ANOVA_RF_df.head)
ANOVA_RF_genes_450wos = set(list(ANOVA_RF_df_450wos["Gene_Symbol"]))
#print(len(ANOVA_RF_genes_450wos))
ANOVA_RF_df_27s = pd.read_csv("data/genes_long_format/sure/27_sim_AnovaRF_markerData_wSMOTE.csv")
#print(ANOVA_RF_df.head)
ANOVA_RF_genes_27s = set(list(ANOVA_RF_df_27s["Gene_Symbol"]))

ANOVA_RF_df_27wos = pd.read_csv("data/genes_long_format/sure/27_sim_AnovaRF_markerData_woSMOTE.csv")
#print(ANOVA_RF_df.head)
ANOVA_RF_genes_27wos = set(list(ANOVA_RF_df_27wos["Gene_Symbol"]))


# Common genes in 27k and 450k dataset

# In[3]:


common_27_450_smote_genes = list(set(ANOVA_RF_genes_450s).intersection(ANOVA_RF_genes_27s))
print(len(ANOVA_RF_genes_450s))
print(len(ANOVA_RF_genes_27s))
print("common genes of 27k smote and 450k smote:",len(common_27_450_smote_genes))

common_27_450_wosmote_genes = list(set(ANOVA_RF_genes_450wos).intersection(ANOVA_RF_genes_27wos))
print(len(ANOVA_RF_genes_450wos))
print(len(ANOVA_RF_genes_27wos))
print("common genes of 27k without smote and 450k without smote:", len(common_27_450_wosmote_genes))


# Gene list from COSMIC

# In[5]:


census = pd.read_csv("data/Census_All[22180]_Long.csv")
census_genes = list(census["Gene Symbol"])
census_genes_synonym = list(census["Synonym"])

print("Gene list from Gene_Symbol:",len(set(census_genes)))
print("Gene list from Synonym:",len(set(census_genes_synonym))) #used this column for analysis


# Gene list from TSGene database

# In[6]:


doc = codecs.open('data/tumor_suppressor.txt','rU')
tumor_suppressor_df = pd.read_csv(doc, sep='\t')
#print(tumor_suppressor)
ts_genes = list(tumor_suppressor_df["GeneSymbol"])
print("Total genes in TSGene database:", len(ts_genes))


# In[7]:


#Combine all genes from two databases
combined_genes_TS_COSMIC = set(census_genes_synonym + ts_genes)
print("Combined oncogenes between two databases:",len(combined_genes_TS_COSMIC))


# Common genes between genes of reduced features and combined oncogenes of two databases

# In[9]:


common_450k_all = list(set(overall_450k_genes).intersection(combined_genes_TS_COSMIC))
print("Total genes in overall 450k dataset:", len(overall_450k_genes))
print("Common genes between overall 450k dataset and TS+COSMIC oncogenes:", len(common_450k_all))
print("\n")
common_450k_smote = list(set(ANOVA_RF_genes_450s).intersection(common_450k_all))
print("Total genes in 450k SMOTE dataset:", len(ANOVA_RF_genes_450s))
print("Common genes between 450k Smote dataset and TS+COSMIC oncogenes:", len(common_450k_smote))
print("\n")
common_450k_wosmote = list(set(ANOVA_RF_genes_450wos).intersection(common_450k_all))
print("Total genes in 450k without SMOTE dataset:", len(ANOVA_RF_genes_450wos))
print("Common genes between 450k wothout Smote dataset and TS+COSMIC oncogenes:", len(common_450k_wosmote))
print("\n")
common_27k_all = list(set(overall_27k_genes).intersection(combined_genes_TS_COSMIC))
print("Total genes in overall 27k dataset:", len(overall_27k_genes))
print("Common genes between overall 27k dataset and TS+COSMIC oncogenes:", len(common_27k_all))
print("\n")
common_27k_smote = list(set(ANOVA_RF_genes_27s).intersection(common_27k_all))
print("Total genes in 27k SMOTE dataset:", len(ANOVA_RF_genes_27s))
print("Common genes between 27k Smote dataset and TS+COSMIC oncogenes:", len(common_27k_smote))
print("\n")
common_27k_wosmote = list(set(ANOVA_RF_genes_27wos).intersection(common_27k_all))
print("Total genes in 27k without SMOTE dataset:", len(ANOVA_RF_genes_27wos))
print("Common genes between 27k wothout Smote dataset and TS+COSMIC oncogenes:", len(common_27k_wosmote))


# In[10]:


smote_27_450_db = list(set(common_450k_smote).intersection(common_27k_smote))
print(len(common_450k_smote))
print(len(common_27k_smote))
print("Common between 27k smote, 450k smote genes, and COSMIC+TSgene:", (smote_27_450_db))
wosmote_27_450_db = list(set(common_450k_wosmote).intersection(common_27k_wosmote))
print("Common between 27k wosmote, 450k wosmote genes, and COSMIC+TSgene:", (wosmote_27_450_db))


# Collected 100 cancer related genes from several resources (https://www.breastcancer.org/risk/risk-factors/genetics, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4478970/)

# In[13]:


sample_gene_list = ["AKT1","CCND1","ERBB2","ESR1","ETV6","FOXA1","GATA3","IRS4","MAP2K4","MAP3K1","MAP3K13","NOTCH1","NTRK3","PIK3CA","PPM1D","SALL4","TBX3","TP53"]
print(len(sample_gene_list))
additional_genes = ["BRCA1", "BRCA2","PALB2","PTEN","ATM","CDH1","CHEK2","NBN","NF1","STK11","BARD1","BRIP1","MLH1", 
                    "MSH2", "MSH6", "PMS2", "EPCAM","RAD51C", "RAD51D","CASP8","CTLA4","NBN","CYP19A1","TERT","XRCC3","AKT1", "APC", "ATR", 
                    "BABAM1", "BAP1", "BARD1", "BMPR1A", "CDK4", "CDKN2A", "CHEK1", "CTNNA1", "EPCAM", "FAM175A", "MLH1", "MRE11A", 
                    "MSH2", "MSH6", "MUTYH", "NBN", "PIK3CA", "PMS2", "POLD1", "POLE", "PRSS1", "RAD50", "RAD51", "RAD51C", "RAD51D", 
                    "RET", "SDHB", "SDHC", "SDHD", "SMAD4", "TP53BP1", "VHL", "XRCC2","BMPR1A", "EPCAM", "GALNT12", "GREM1", 
                    "MLH1", "MSH2", "MSH6", "MUTYH", "PMS2", "SMAD4","BARD1", "MRE11A", "MUTYH", "NBN", "NF1", "RAD50", "RAD51C",
                    "RAD51D","EPCAM", "MRE11A", "BLM", "FANCC", "HOXB13", "MLH1", "MRE11A", "XRCC2","AXIN2", "BLM", "FANCC", 
                    "PALLD", "PMS2", "RAD50", "RAD51C", "RAD51D", "SMAD4", "VHL", "XRCC2", "PMS1", "CDC73", "CDKN1C", "CDKN2A", "FH", "FLCN", 
                    "GPC3", "MAX", "MEN1", "MET", "NF1", "NF2", "PHOX2B", "PRKAR1A", "PTCH1", "RAD51C", "RAD51D", "RET", "SDHAF2", "SDHB", 
                    "SDHC", "SDHD", "SUFU", "TMEM127", "VHL", "WT1","ERCC4", "RAD51C", "RAD51D", "XRCC2", "FAM175A", "FANCA", "FANCC", "FANCD2", 
                    "FANCE", "FANCF", "FANCG", "MEN1", "MLH1", "MRE11A", "MSH2", "MSH3", "MSH6", "PTCH1", "RAD50", "RAD51C", "RAD51D", "SLX4", 
                    "UIMC1", "XRCC2"]
sample_gene_list.extend(additional_genes)
total_sample_genes = set(sample_gene_list)
print(len(total_sample_genes))


# Compared 100 sample genes with our reduced gene sets

# In[15]:


db_common_450k_smote = list(set(overall_450k_genes).intersection(total_sample_genes))
print("Common genes between 450k all dataset and sample genes:", len(db_common_450k_smote))
db_common_450k_smote = list(set(ANOVA_RF_genes_450s).intersection(total_sample_genes))
print("Common genes between 450k Smote dataset and sample genes:", len(db_common_450k_smote))
db_common_450k_wosmote = list(set(ANOVA_RF_genes_450wos).intersection(total_sample_genes))
print("Common genes between 450k without Smote dataset and sample genes:", len(db_common_450k_wosmote))
db_common_450k_smote = list(set(overall_27k_genes).intersection(total_sample_genes))
print("Common genes between 27k all dataset and sample genes:", len(db_common_450k_smote))
db_common_27k_smote = list(set(ANOVA_RF_genes_27s).intersection(total_sample_genes))
print("Common genes between 27k Smote dataset and sample genes:", len(db_common_27k_smote))
db_common_27k_wosmote = list(set(ANOVA_RF_genes_27wos).intersection(total_sample_genes))
print("Common genes between 27k without Smote dataset and sample genes:", (db_common_27k_wosmote))
db_and_sample = list(set(combined_genes_TS_COSMIC).intersection(total_sample_genes))
print("\nCommon genes between 2 databases and sample genes:", len(db_and_sample))

