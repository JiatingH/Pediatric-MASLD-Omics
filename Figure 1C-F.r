################################################################################
# Project: MASLD Metagenomics Analysis
# Task: Phylum-level Taxonomic Comparison
# Aim: Identify and visualize significantly different phyla across groups
################################################################################

# --- 1. Load Libraries ---
library(vegan)
library(Hmisc)
library(viridis)
library(dplyr)
library(readr)
library(readxl)
library(gdata)
library(reshape2)
library(pheatmap)
library(tableone)
library(ggpubr)
library(stringr)
library(RColorBrewer)

# --- 2. Data Import & Cleaning ---
# Replace absolute local paths with relative project paths for portability
data_dir <- "./data/taxonomy/"
output_dir <- "./results/plots/"
script_dir <- "./scripts/functions/"

# Load Phylum-level abundance data
phylum = read.table(paste0(data_dir, 'all_phylum_abundance.csv'), 
                    sep = ',', header = 1, row.names = 1, check.names = FALSE)

# Extract abundance matrix (excluding the Group column)
phylum.data = phylum[setdiff(colnames(phylum), 'Group')]

# Standardize group nomenclature (Update 'NAFLD' to 'MASLD')
phylum$Group = ifelse(phylum$Group == 'NAFLD', 'MASLD', phylum$Group)

# Replace spaces with underscores in Phylum names to avoid syntax issues in R
colnames(phylum) = gsub(' ', '_', colnames(phylum))
colnames(phylum.data) = gsub(' ', '_', colnames(phylum.data))

# --- 3. Load Custom Functional Scripts ---
# Keeping original source filenames as requested
source(paste0(script_dir, 'KW.r'))
source(paste0(script_dir, 'Wilcox.r'))
source(paste0(script_dir, 'FoldChange.r'))
source(paste0(script_dir, 'BoxPlot.r'))

# --- 4. Statistical Analysis ---
# Prepare the list of phyla names for comparison
phylum.name = colnames(phylum.data)

# Perform Kruskal-Wallis (KW) test to find differences across groups
# This custom function identifies features that vary significantly
phylum.kw = KW(data = phylum, compare = phylum.name, class = 'Group')

# Filter for significant phyla using a False Discovery Rate (FDR) threshold < 0.05
sig.phylum.name = rownames(filter(phylum.kw, fdr < 0.05))

# --- 5. Visualization ---
# Generate a multi-panel PDF containing boxplots for all significant phyla
# 


MultipleBoxPlot(phylum, 
                x = 'Group', 
                y = sig.phylum.name,
                group_elements = c('HC', 'OB', 'MASLD'),
                color_values = c('#CD6155', '#040676', '#5F93CE'),
                pdf_path = paste0(output_dir, 'Significant_Phyla_Boxplots.pdf'))