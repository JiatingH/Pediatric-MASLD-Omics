################################################################################
# Project: MASLD Multi-omics Research
# Task: Taxonomic Contribution to Functionally-Relevant KOs
# Aim: Identify top 5 microbial contributors to Carbohydrate/Glycan metabolic KOs
################################################################################

# --- 1. Load Libraries ---
library(randomForest); library(datasets); library(caret); library(pROC)
library(readxl); library(gdata); library(reshape2); library(Hmisc)
library(pheatmap); library(tableone); library(ggpubr); library(stringr)
library(RColorBrewer); library(viridis); library(dplyr); library(matrixStats)
library(vegan); library(tidyr); library(KEGGREST)

# --- 2. Data Import & Configuration ---
# Define public-ready paths
input_dir <- "./data/metagenomics/"
output_dir <- "./results/contribution_analysis/"
script_dir <- "./scripts/functions/"

# Load KO abundance and sample metadata
all.ko.temp = read.csv(paste0(input_dir, 'all_ko_abundance.csv'), sep = ',', header = 1, row.names = 1, check.names = FALSE)
group.info = read.csv(paste0(input_dir, 'sample_metadata.csv'), sep = ',', header = 1, row.names = 1, check.names = FALSE)

# Standardize group labeling (Combining NAFL/MASH into MASLD)
group.info$ClassNote = ifelse(group.info$ClassNote %in% c('NAFL','MASH'), 'MASLD', group.info$ClassNote)

# Source custom functions (Keeping original filenames)
source(paste0(script_dir, 'KW.r'))
source(paste0(script_dir, 'KO_Contribution.r'))

# Load correlation results from previous clinical association step
ko.carbohydrate.glycan.cor = read.table(paste0(input_dir, 'diff_ko_clinical_correlations.csv'), sep = ',', header = 1, row.names = 1)

# --- 3. Run Contribution Analysis ---
# Extract top 5 contributing species for the targeted KOs
ko.carbohydrate.glycan.heatmap.pre = KO_HeatMap(all.ko.temp, ko_name = ko.carbohydrate.glycan.cor$KO, group.info = group.info, n = 5)

# Extract contributor abundance and significance (P-values)
ko.carbohydrate.glycan.heatmap = ko.carbohydrate.glycan.heatmap.pre$Contribute
ko.carbohydrate.glycan.heatmap.p = ko.carbohydrate.glycan.heatmap.pre$P[colnames(ko.carbohydrate.glycan.heatmap)]

# --- 4. Subgroup Analysis: Positive vs. Negative Correlations ---
# Filter KOs by clinical correlation direction (ALT, Fibro_Pen, or kPa)
ko.carbohydrate.glycan.cor.negative = ko.carbohydrate.glycan.cor %>% filter(ALT < 0 | Fibro_Pen < 0 | kPa < 0)
ko.carbohydrate.glycan.cor.positive = ko.carbohydrate.glycan.cor %>% filter(ALT > 0 | Fibro_Pen > 0 | kPa > 0)

# Extract specific contribution data for both subsets
ko.carbohydrate.glycan.cor.negative.pre = KO_HeatMap(all.ko.temp, ko_name = rownames(ko.carbohydrate.glycan.cor.negative), group.info = group.info, n = 5)
ko.carbohydrate.glycan.cor.positive.pre = KO_HeatMap(all.ko.temp, ko_name = rownames(ko.carbohydrate.glycan.cor.positive), group.info = group.info, n = 5)

# Prepare matrices for Positive correlation contributors
ko.carbohydrate.glycan.cor.positive.heatmap = ko.carbohydrate.glycan.cor.positive.pre$Contribute
ko.carbohydrate.glycan.cor.positive.p = ko.carbohydrate.glycan.cor.positive.pre$P[colnames(ko.carbohydrate.glycan.cor.positive.heatmap)]

# Prepare matrices for Negative correlation contributors
ko.carbohydrate.glycan.cor.negative.heatmap = ko.carbohydrate.glycan.cor.negative.pre$Contribute
ko.carbohydrate.glycan.cor.negative.p = ko.carbohydrate.glycan.cor.negative.pre$P[colnames(ko.carbohydrate.glycan.cor.negative.heatmap)]

# --- 5. Visualization Preparation (Annotation & Colors) ---
# Color coding based on KEGG Pathway Class
ko.carbohydrate.glycan.annotation.col = ko.carbohydrate.glycan.cor[, 'Class', drop = FALSE]
ko.carbohydrate.glycan.col.color = AssignColor(unique(ko.carbohydrate.glycan.annotation.col$Class))

# Assign colors for the positive and negative subsets
ko.carbohydrate.glycan.annotation.positive = ko.carbohydrate.glycan.annotation.col[rownames(ko.carbohydrate.glycan.cor.positive), , drop = FALSE]
ko.carbohydrate.glycan.annotation.positive.color = AssignColor(unique(ko.carbohydrate.glycan.annotation.positive$Class))

ko.carbohydrate.glycan.annotation.negative = ko.carbohydrate.glycan.annotation.col[rownames(ko.carbohydrate.glycan.cor.negative), , drop = FALSE]
ko.carbohydrate.glycan.annotation.negative.color = AssignColor(unique(ko.carbohydrate.glycan.annotation.negative$Class))

# Replace zeros with NA for better heatmap visualization (white background for zeros)
ko.carbohydrate.glycan.heatmap[ko.carbohydrate.glycan.heatmap == 0] = NA
ko.carbohydrate.glycan.cor.positive.heatmap[ko.carbohydrate.glycan.cor.positive.heatmap == 0] = NA
ko.carbohydrate.glycan.cor.negative.heatmap[ko.carbohydrate.glycan.cor.negative.heatmap == 0] = NA

# --- 6. Plotting and Saving Heatmaps ---
library(viridis)

# Global Heatmap
pheatmap(ko.carbohydrate.glycan.heatmap, 
         cluster_row = FALSE, cluster_col = FALSE, 
         color = plasma(100), na_col = "white", 
         fontsize_row = 10, cellwidth = 10, cellheight = 10,
         annotation_col = ko.carbohydrate.glycan.annotation.col,
         annotation_colors = list(Class = ko.carbohydrate.glycan.col.color),
         # Add '*' where species abundance differs significantly between groups (p < 0.05)
         display_numbers = matrix(ifelse(ko.carbohydrate.glycan.heatmap.p < 0.05, '*', ''), 
                                  nrow = nrow(ko.carbohydrate.glycan.heatmap.p)), 
         filename = paste0(output_dir, "KO_Species_Contribution_Global.pdf"))

# Positive Correlation Heatmap (Potential Pathogenic Factors)
pheatmap(ko.carbohydrate.glycan.cor.positive.heatmap, 
         cluster_row = FALSE, cluster_col = FALSE, 
         color = plasma(100), na_col = "white", 
         fontsize_row = 10, cellwidth = 10, cellheight = 10,
         annotation_col = ko.carbohydrate.glycan.annotation.positive,
         annotation_colors = list(Class = ko.carbohydrate.glycan.annotation.positive.color),
         display_numbers = matrix(ifelse(ko.carbohydrate.glycan.cor.positive.p < 0.05, '*', ''), 
                                  nrow = nrow(ko.carbohydrate.glycan.cor.positive.p)), 
         filename = paste0(output_dir, "KO_Species_Contribution_Positive.pdf"))

# Negative Correlation Heatmap (Potential Protective Factors)
pheatmap(ko.carbohydrate.glycan.cor.negative.heatmap, 
         cluster_row = FALSE, cluster_col = FALSE, 
         color = plasma(100), na_col = "white", 
         fontsize_row = 10, cellwidth = 10, cellheight = 10,
         annotation_col = ko.carbohydrate.glycan.annotation.negative,
         annotation_colors = list(Class = ko.carbohydrate.glycan.annotation.negative.color),
         display_numbers = matrix(ifelse(ko.carbohydrate.glycan.cor.negative.p < 0.05, '*', ''), 
                                  nrow = nrow(ko.carbohydrate.glycan.cor.negative.p)), 
         filename = paste0(output_dir, "KO_Species_Contribution_Negative.pdf"))