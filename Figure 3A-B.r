################################################################################
# Project: Multi-omics Analysis of Metabolic Dysfunction (MASLD)
# Task: KEGG Pathway Enrichment & Differential KO Analysis (Figures 3A-B)
# Aim: Identify specific microbial metabolic shifts in disease progression
################################################################################

# --- 1. Load Libraries ---
library(randomForest); library(datasets); library(caret); library(pROC)
library(readxl); library(gdata); library(reshape2); library(Hmisc)
library(pheatmap); library(tableone); library(ggpubr); library(stringr)
library(RColorBrewer); library(viridis); library(dplyr); library(matrixStats)
library(vegan); library(tidyr); library(KEGGREST); library(xgboost); library(shapr)
library(ComplexHeatmap); library(circlize)

# --- 2. Data Import & Preprocessing ---
# Define public paths
input_dir <- "./Project_Data/KEGG_Abundance/"
output_dir <- "./Project_Results/Function/"

# Load KEGG and Clinical Metadata
kegg <- read.table(paste0(input_dir, 'All_KEGG_Abundance.csv'), sep =',', header = 1, row.names = 1, check.names = FALSE)
clinical <- read.table(paste0(input_dir, 'Clinical_Metadata.csv'), sep = ',', header = 1, row.names = 1, check.names = FALSE)

# Filter for Metabolic Pathways (KEGG Level 1 == 'Metabolism')
kegg.metabolism <- filter(kegg, KEGGLevel1 == 'Metabolism')

# Intersection of samples between clinical and KEGG data
common_samples <- intersect(colnames(kegg), rownames(clinical))
kegg.data <- kegg[, common_samples]
clinical.sub <- clinical[common_samples, ]

# Data Transformation: Samples as rows, Pathways as columns
kegg.df <- as.data.frame(t(kegg.data))
kegg.df <- kegg.df[, rownames(kegg.metabolism)]

# --- 3. Statistical Analysis (Kruskal-Wallis) ---
# Filter pathways with low prevalence (remove if < 50% across samples)
source('./Scripts/DataMining_Utils.r') # Custom cleaning function
kegg.df <- Remove_HighNA_Zero(kegg.df, threshold = 0.5, axis = 'col')
kegg.df$Group <- clinical.sub$Group

# Run multi-group comparison (KW Test)
source('./Scripts/Statistical_Tests.r')
kegg.kw <- KW_Test(kegg.df, compare = setdiff(colnames(kegg.df), 'Group'), class = 'Group')
kegg.p <- filter(kegg.kw, fdr < 0.05) # Filter by FDR < 0.05

# --- 4. Enrichment Analysis (Hypergeometric Test) ---
# Calculate if a pathway is enriched with significant KOs

# Load KO-level abundance
ko.abundance <- read.csv(paste0(input_dir, 'All_KO_Data.csv'), header = 1, row.names = 1, check.names = FALSE)
ko.table <- as.data.frame(t(ko.abundance[, common_samples]))
ko.table$Group <- clinical.sub$Group

# Identify significant KOs (KW Test)
ko.kw <- KW_Test(ko.table, compare = colnames(ko.table)[-ncol(ko.table)], class = 'Group')
ko.sig <- filter(ko.kw, fdr < 0.05)

# Hypergeometric test function: phyper(q, m, n, k)
# q: Significant KOs in pathway, m: Total significant KOs, 
# n: Total non-significant KOs, k: Total KOs in pathway
calculate_phyper <- function(df, total_ko, total_sig_ko) {
  df$p_value <- phyper(df$Diff_KO_Count - 1, total_sig_ko, total_ko - total_sig_ko, df$Total_Count, lower.tail = FALSE)
  return(df)
}

# --- 5. Visualization: Figure 3A (Pathway Enrichment Barplot) ---

pdf(file = paste0(output_dir, "Figure3A_KEGG_Enrichment.pdf"))
ggplot(diff.KEGG.KO.summary, aes(x = reorder(Class, enrichment_ratio), y = enrichment_ratio, fill = -log10(phyper_p_value))) + 
  geom_bar(stat = "identity") + 
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  geom_text(aes(label = Diff_KO_Count), vjust = -0.5) + 
  theme_classic() + 
  labs(x = "KEGG Pathway Class", y = "Enrichment Ratio", title = "Pathway Enrichment Analysis") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# --- 6. Visualization: Figure 3B (Circular Heatmap of KO Fold Changes) ---

# Focus on Carbohydrate and Glycan Metabolism
# Prepare Fold Change (FC) data for HC vs OB vs MASLD
pdf(file = paste0(output_dir, "Figure3B_Circular_KO_Heatmap.pdf"))
circos.par(gap.after = c(30)) # Create visual gap for legend

# Track 1: Disease Subgroup FC (Advanced vs Simple MASLD)
circos.heatmap(subgroup_fc_data, col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), 
               track.height = 0.05, rownames.side = "outside")

# Track 2: Global 3-Group FC (HC/OB/MASLD)
circos.heatmap(global_fc_data, col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), 
               track.height = 0.2, cell.border = "white")

# Track 3: Pathway Classification mapping
circos.heatmap(class_mapping, col = category_colors, track.height = 0.05)

# Add Legends
draw(Legend(title = "Log10(FC)", col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))), 
     x = unit(0.2, "npc"), y = unit(0.5, "npc"))
circos.clear()
dev.off()