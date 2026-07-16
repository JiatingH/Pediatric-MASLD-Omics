################################################################################
# Project: MASLD Metagenomics - Functional Core Analysis
# Task: TCA Cycle (ko00020) Taxonomic Contribution (Figure 5C)
# Aim: Identify microbial species driving differential TCA cycle genes
################################################################################

# --- 1. Load Libraries ---
library(datasets); library(caret); library(pROC); library(readxl)
library(gdata); library(reshape2); library(Hmisc); library(pheatmap)
library(tableone); library(ggpubr); library(stringr); library(dplyr)
library(RColorBrewer); library(viridis); library(matrixStats); library(vegan)
library(tidyr); library(KEGGREST)

# --- 2. Data Import & Configuration ---
# Define standardized relative paths for public use
input_dir <- "./data/functional_analysis/"
result_dir <- "./results/pathway_contribution/"
script_dir <- "./scripts/microbiome_functions/"

# Load metagenomic KO abundance and sample metadata
all.ko.temp = read.csv(paste0(input_dir, 'all_ko_abundance_temp.csv'), sep = ',', header = 1, row.names = 1, check.names = FALSE)
group.info = read.csv(paste0(input_dir, 'metagenomics_sample_info.csv'), sep = ',', header = 1, row.names = 1, check.names = FALSE)

# Standardize group labeling (Merging NAFL/NASH into MASLD)
group.info$ClassNote = ifelse(group.info$ClassNote %in% c('NAFL','NASH'), 'MASLD', group.info$ClassNote)

# Load summary of differential KOs within KEGG pathways
diff.ko.summary.in.path = read.table(paste0(input_dir, 'diff_ko_pathway_summary.csv'), sep = ',', header = 1, row.names = 1)

# --- 3. Pathway-Specific Analysis: Citrate cycle (TCA cycle) ---
# Example: Figure 5C logic
# Step A: Extract all KO IDs associated with the TCA Cycle
ko00020 = strsplit(diff.ko.summary.in.path['Citrate cycle (TCA cycle)', 'KO'], ',')[[1]]

# Step B: Generate a stacked barplot showing overall species contribution to this pathway
# 
Pathway_Top5_Contributed_StackedBar(
    data = all.ko.temp, 
    ko_name = ko00020, 
    group.info = group.info,
    save_path = paste0(result_dir, 'TCA_Cycle_Pathway_Contribution_StackedBar.pdf')
)

# Step C: Extract only the significantly Differential KOs within the TCA Cycle
diff.ko00020 = strsplit(diff.ko.summary.in.path['Citrate cycle (TCA cycle)', 'Diff_KO'], ',')[[1]]

# Step D: Generate detailed Pie Plots (contribution) and Box Plots (abundance comparison)
# for the top 5 species contributing to these specific differential KOs
# 
output_subdir <- paste0(result_dir, 'TCA_Detailed_Analysis/')
if(!dir.exists(output_subdir)) dir.create(output_subdir, recursive = TRUE)

Path_KO_PiePlot_BoxPlot(
    data = all.ko.temp, 
    ko_name = diff.ko00020, 
    group.info = group.info, 
    n = 5,
    save_path1 = output_subdir, # Path for Pie Plots
    save_path2 = output_subdir  # Path for Box Plots
)