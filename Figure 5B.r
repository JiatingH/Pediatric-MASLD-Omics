################################################################################
# Project: MASLD Multi-omics Analysis
# Task: Metabolomics Preprocessing and Differential Abundance Analysis
# Aim: Identify significant metabolites across HC, OB, and MASLD groups
################################################################################

# --- 1. Load Libraries ---
library(readxl)
library(gdata)
library(reshape2)
library(Hmisc)
library(pheatmap)
library(tableone)
library(ggpubr)
library(stringr)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(matrixStats) 
library(vegan)
library(corrplot)
library(ade4)
library(readr)

# --- 2. Load Custom Functional Scripts ---
# Keeping original source filenames as requested
source('./scripts/functions/KW.r')
source('./scripts/functions/Wilcox.r')
source('./scripts/functions/FoldChange.r')
source('./scripts/functions/DataMining.r')

# --- 3. Data Import ---
# Using generalized relative paths for public deployment
data_dir <- "./data/metabolomics/"
output_dir <- "./results/metabolites/"

# Load Metabolite Abundance Data
metabolites = read.csv(paste0(data_dir, 'metabolites_merged_batches.csv'), sep = ',', header = 1, row.names = 1)
# Create unique rownames combining KEGG ID and original name
rownames(metabolites) = paste(metabolites$KEGG.ID, rownames(metabolites), sep = '_')

# Load Metabolite Sample Mapping Info
sample.info = read.table(paste0(data_dir, 'metabolites_sample_info.csv'), sep = ',', header = 1, row.names = 1)

# Load Clinical Metadata
clinical = read.table(paste0(data_dir, 'clinical_metadata_standardized.csv'), sep = ',', header = 1, row.names = 1)
# Simplify subgroup classification (Merging MASLD Fibrosis and MASH)
clinical$Class_Two_Subgroup = ifelse(clinical$Class_Three_Subgroup %in% c('MASLD_Fibrosis', 'MASH'), 'MASH_Fibrosis', clinical$Class_Three_Subgroup)

# Extract matching clinical and metabolomic subsets
metabolites.clinical = clinical[sample.info$Hospital.Num, ]
metabolites.data = metabolites[sample.info$Submit]

# --- 4. Data Preprocessing ---
# Step 1: Remove metabolites with >40% missing values (NAs)
metabolites.filter = Remove_HighNA(metabolites.data, axis = 'row', threshold = 0.4)

# Step 2: Transpose for sample-wise processing
metabolites.filter.data = t(metabolites.filter) %>% as.data.frame

# Step 3: Impute missing values using 1/2 of the minimum value detected for that metabolite
metabolites.fillna = Fill_MissingValues(metabolites.filter.data, method = 'half_min')

# Step 4: Finalize data frame and re-attach external identifiers (HMDB ID)
metabolites.fillna.data = t(metabolites.fillna) %>% as.data.frame
metabolites.fillna.data$HMDBID = metabolites[rownames(metabolites.fillna.data), 'HMDB.ID']

# --- 5. Statistical Analysis ---
# Perform Kruskal-Wallis test across the three groups (HC, OB, MASLD)
metabolites.kw = KW(data = metabolites.fillna, compare = colnames(metabolites.fillna), class = 'Group')

# Filter for metabolites with nominal significance (p < 0.05)
metabolites.sig = filter(metabolites.kw, p < 0.05)
metabolites.sig.name = rownames(metabolites.sig)

# --- 6. Visualization ---
# Source the fold change analysis and plotting script
source('./scripts/plotting/FoldChangeAnalysis.r')

# Generate fold change objects for the significant metabolites
metabolites.obj = FoldChangeAnalysis(data = metabolites.fillna, 
                                     com = metabolites.sig.name, 
                                     class = 'Group', 
                                     group_label = c('HC', 'OB', 'MASLD'))

# Generate and save the clustered heatmap of differential metabolites
heatmap.FoldChangeAnalysis(metabolites.obj, 
                           save_filename = paste0(output_dir, 'Differential_Metabolites_Clustering_Heatmap.pdf'))