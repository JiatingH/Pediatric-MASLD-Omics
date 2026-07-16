################################################################################
# Project: MASLD Multi-omics Integration
# Task: Adonis (PERMANOVA) Analysis for Clinical Factor Association (Figure S1)
# Aim: Determine the effect size (R2) of clinical traits on multi-omics variance
################################################################################

# --- 1. Load Libraries ---
library(readxl); library(gdata); library(reshape2); library(Hmisc)
library(pheatmap); library(tableone); library(ggpubr); library(stringr)
library(dplyr); library(RColorBrewer); library(viridis)
library(matrixStats); library(vegan); library(corrplot); library(ade4); library(readr)

# Set working directory to project root (Generalize path for portability)
# setwd('./Project_Root_Directory/')

# --- 2. Data Preparation ---

# Load Clinical Data
all.clinical = read.csv('./Data/Clinical/All.Clinical.Metadata.csv', sep = ',', header = 1, row.names = 1) 
rownames(all.clinical) = all.clinical$Hospital.Num

# Load Metagenomic Data (Taxonomy and Functional Tiers)
all.taxonomy = read.csv('./Data/Metagenomics/Taxonomy_Table.csv', header = 1, row.names = 1, check.names = FALSE)
all.species = read.csv('./Data/Metagenomics/Species_Level_Abundance.csv', header = 1, row.names = 1, check.names = FALSE)
all.kegg = read.csv('./Data/Metagenomics/KEGG_Pathways.csv', header = 1, row.names = 1, check.names = FALSE)
all.metacyc = read.csv('./Data/Metagenomics/MetaCyc_Pathways.csv', header = 1, row.names = 1, check.names = FALSE)
all.rxn = read.csv('./Data/Metagenomics/Enzymatic_Reactions.csv', header = 1, row.names = 1, check.names = FALSE)
all.ko = read.csv('./Data/Metagenomics/KO_Abundance.csv', header = 1, row.names = 1, check.names = FALSE)

# Load and Preprocess CAZyme Data
all.cazyme = read.csv('./Data/Metagenomics/CAZyme_Abundance.csv', header = 1, row.names = 1, check.names = FALSE)
all.cazyme[all.cazyme == 0] = NA
# Fill_MissingValues is a custom function from previous steps
all.cazyme = Fill_MissingValues(all.cazyme, method = 'half_min')

# Load Sample Mapping Info
all.feces.sample.info = read.csv('./Data/Metadata/Sample_Mapping_Info.csv', sep = ',', header = 1, row.names = 1, check.names = FALSE)

# Align Taxonomy data with clinical IDs
all.taxonomy.data = all.taxonomy[all.feces.sample.info$MicrobioID,]
rownames(all.taxonomy.data) = all.feces.sample.info$Hospital.Num

# Align Functional and Metabolomic data
all.kegg.data = t(all.kegg[all.feces.sample.info$Hospital.Num]) %>% as.data.frame
all.ko.data = t(all.ko[all.feces.sample.info$Hospital.Num]) %>% as.data.frame
all.metacyc.data = t(all.metacyc[all.feces.sample.info$Hospital.Num]) %>% as.data.frame
all.rxn.data = t(all.rxn[all.feces.sample.info$Hospital.Num]) %>% as.data.frame
all.cazyme.data = t(all.cazyme[intersect(all.feces.sample.info$Hospital.Num, colnames(all.cazyme))]) %>% as.data.frame

# Load Metabolomics Data
metabolites = read.csv('./Data/Metabolomics/Metabolite_Concentration.csv', sep = ',', header = 1, row.names = 1)
metabolites.info = read.csv('./Data/Metadata/Metabolite_Sample_Info.csv')

# Sync Clinical data for both Metagenomics and Metabolomics cohorts
metabolites.clinical = all.clinical[rownames(metabolites),]
taxonomy.clinical = all.clinical[rownames(all.taxonomy.data),]

# --- 3. Clinical Factor Categorization ---
# Load Adonis Analysis custom function
source('./Scripts/Functions/Adonis.r')

# Define clinical indices by category
basic.index = c('Age', 'Gender', 'Height', 'Weight', 'BMI', 'BMI_Zscore', 'Waist', 'Hip', 'Waist.Hip', 'SP', 'DP')
blood.index = c('WBC', 'Lymphocyte_Count', 'Neutrophil_Count', 'Monocyte_Count', 'Eosinophil_Count', 'Basophil_Count', 'HCT', 'RBC', 'Hemoglobin', 'Platelet_Count')
kidney.index = c('ALB', 'GLB')
liver.index = c('ALT', 'AST', 'kPa', 'Fibro_Pen', 'AKP', 'TBIL', 'DBIL', 'IBIL', 'GGT', 'ADA')
lipid.index = c('TG', 'TC', 'HDL', 'LDL', 'non.HDL')
glucose.index = c('Glu0', 'Glu30', 'Glu60', 'Glu120', 'Insulin0', 'Insulin30', 'Insulin60', 'Insulin120',
                  'Cpeptide0', 'Cpeptide30', 'Cpeptide60', 'Cpeptide120', 'HbA1c', 'HOMA_IR', 'GluAUC', 'CpAUC', 'InsAUC')

compare.index = c(basic.index, blood.index, kidney.index, liver.index, lipid.index, glucose.index)

# Exclude kPa (FibroScan) from general analysis due to missing values (analyzed separately)
adonis.clinical.index <- setdiff(compare.index, "kPa")
adonis.index = c('Group', 'Class_MASH', 'Class_Fibrosis_2', 'Class_Three_Subgroup', adonis.clinical.index)

# Create combined DataFrames for Adonis Analysis
all.taxonomy.anova = cbind(all.taxonomy.data, taxonomy.clinical[adonis.index])
all.kegg.anova = cbind(all.kegg.data, taxonomy.clinical[adonis.index])
all.ko.anova = cbind(all.ko.data, taxonomy.clinical[adonis.index])
all.metacyc.anova = cbind(all.metacyc.data, taxonomy.clinical[adonis.index])
all.rxn.anova = cbind(all.rxn.data, taxonomy.clinical[adonis.index])
all.metabolites.anova = cbind(metabolites, metabolites.clinical[adonis.index])
all.cazyme.anova = cbind(all.cazyme.data, taxonomy.clinical[rownames(all.cazyme.data), adonis.index])

# Separate subset analysis for kPa (FibroScan) - Filtering for non-missing values
kpa.taxonomy.anova = filter(cbind(all.taxonomy.data, taxonomy.clinical), kPa > 0)
kpa.taxonomy.data = all.taxonomy.data[rownames(kpa.taxonomy.anova),]
# (Repeat same filtering logic for other omics types if kPa analysis is required)

# --- 4. Perform Adonis Analysis ---
taxonomy.adonis.out = perform_adonis_analysis(all.taxonomy.data, all.taxonomy.anova, adonis.index)
kegg.adonis.out = perform_adonis_analysis(all.kegg.data, all.kegg.anova, adonis.index)
ko.adonis.out = perform_adonis_analysis(all.ko.data, all.ko.anova, adonis.index)
metacyc.adonis.out = perform_adonis_analysis(all.metacyc.data, all.metacyc.anova, adonis.index)
rxn.adonis.out = perform_adonis_analysis(all.rxn.data, all.rxn.anova, adonis.index)
metabolites.adonis.out = perform_adonis_analysis(metabolites, all.metabolites.anova, adonis.index)
cazyme.adonis.out = perform_adonis_analysis(all.cazyme.data, all.cazyme.anova, adonis.index)

# Special analysis for kPa (FibroScan stiffness)
kpa.taxonomy.adonis.out = perform_adonis_analysis(kpa.taxonomy.data, kpa.taxonomy.anova, c('kPa'))
# ... (perform for other omics)

# --- 5. Compile Results Table ---

# Extract R2 (Variance Explained)
adonis.r = data.frame(Taxonomy = taxonomy.adonis.out$R2, KEGG = kegg.adonis.out$R2, KO = ko.adonis.out$R2,
                      MetaCyc = metacyc.adonis.out$R2, RXN = rxn.adonis.out$R2, 
                      Metabolites = metabolites.adonis.out$R2, CAZyme = cazyme.adonis.out$R2)

# Extract P-values
adonis.p = data.frame(Taxonomy = taxonomy.adonis.out$`Pr(>F)`, KEGG = kegg.adonis.out$`Pr(>F)`, 
                      KO = ko.adonis.out$`Pr(>F)`, MetaCyc = metacyc.adonis.out$`Pr(>F)`, 
                      RXN = rxn.adonis.out$`Pr(>F)`, Metabolites = metabolites.adonis.out$`Pr(>F)`, 
                      CAZyme = cazyme.adonis.out$`Pr(>F)`)

rownames(adonis.r) = adonis.index
rownames(adonis.p) = adonis.index

# Add FibroScan results to the final table
# (Assuming kpa.r and kpa.p were extracted as per above logic)
adonis.R2 = rbind(adonis.r, kpa.adonis.r)
adonis.P = rbind(adonis.p, kpa.adonis.p)

# Format R2 as percentage with 2 decimal places
adonis.R2 = round(adonis.R2 * 100, 2)

# Define final visualization order for rows (Clinical traits)
order = c('Group', 'Class_MASH', 'Class_Fibrosis_2', 'Class_Three_Subgroup', 'Age', 'Gender', 
          'Height', 'Weight', 'BMI', 'BMI_Zscore', 'Waist', 'Hip', 'Waist.Hip', 'SP', 'DP', 
          'WBC', 'Lymphocyte_Count', 'Neutrophil_Count', 'Monocyte_Count', 'Eosinophil_Count', 
          'Basophil_Count', 'HCT', 'RBC', 'Hemoglobin', 'Platelet_Count', 'ALB', 'GLB', 'ALT', 
          'AST', 'FibroScan', 'Fibro_Pen', 'AKP', 'TBIL', 'DBIL', 'IBIL', 'GGT', 'ADA', 'TG', 'TC', 'HDL', 
          'LDL', 'non.HDL', 'Glu0', 'Glu30', 'Glu60', 'Glu120', 'Insulin0', 'Insulin30', 
          'Insulin60', 'Insulin120', 'Cpeptide0', 'Cpeptide30', 'Cpeptide60', 'Cpeptide120', 
          'HbA1c', 'HOMA_IR', 'GluAUC', 'CpAUC', 'InsAUC')

adonis.R2.table = adonis.R2[order,]
adonis.P.table = adonis.P[order,]

# Annotate significance with asterisks (*) for p < 0.05
adonis.R2.df <- as.data.frame(
  mapply(function(r2, p) ifelse(p < 0.05, paste0(r2, "*"), as.character(r2)),
         as.data.frame(lapply(adonis.R2.table, as.character)), adonis.P.table, SIMPLIFY = FALSE)
)

# --- 6. Final Visualization (Figure S1 Heatmap) ---

bk = c(seq(0, 20, by = 1))
colors <- colorRampPalette(c("white", "#E64B35CC"))(20)

pheatmap(adonis.R2.table,
         cluster_row = FALSE,
         cluster_col = FALSE,
         cellwidth = 45,
         cellheight = 20,
         color = colors,
         legend_breaks = seq(0, 20, 5),
         breaks = bk,
         display_numbers = format(adonis.R2.df),
         fontsize_number = 13,
         color_number = "black",
         filename = './Results/Heatmaps/FigureS1_Adonis_Omics_Associations.pdf'
)