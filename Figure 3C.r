################################################################################
# Project: Multi-omics Analysis of MASLD Progression
# Task: Functional KO & Clinical Phenotype Correlation (Figure 3C)
# Aim: Identify specific metabolic KOs significantly associated with liver fibrosis
################################################################################

# --- 1. Load Libraries ---
library(readxl); library(gdata); library(reshape2); library(Hmisc)
library(tableone); library(ggpubr); library(stringr); library(dplyr)
library(tidyr); library(RColorBrewer); library(viridis); library(vegan)
library(corrplot); library(ade4); library(readr); library(ggplot2)
library(gridExtra); library(ggalluvial)

# --- 2. Data Import ---
# Using public/generalized paths
input_dir <- "./Project_Data/Functional_Analysis/"
output_dir <- "./Project_Results/Correlation_Models/"

# Load Differential KO abundance and metadata
diff_ko_data <- read.table(paste0(input_dir, 'Differential_KO_Abundance.csv'), sep = ',', header = 1, row.names = 1, check.names = FALSE)
ko_description <- read.csv(paste0(input_dir, 'KO_KEGG_Descriptions.csv'), sep = ',', header = 1, row.names = 1, check.names = FALSE)
clinical_data <- read.table(paste0(input_dir, 'Clinical_Indices.csv'), sep = ',', header = 1, row.names = 1)

# Define Clinical Index Groupings
liver_indices <- c('ALT', 'AST', 'kPa', 'Fibro_Pen', 'GGT') # Liver damage & Fibrosis
lipid_indices <- c('TG', 'TC', 'HDL', 'LDL')               # Lipid profile
glucose_indices <- c('Glu0', 'HbA1c', 'HOMA_IR')           # Glucose metabolism

# Target indices for correlation
target_clinical <- c(liver_indices, lipid_indices, glucose_indices)

# Merge KO data with target clinical variables (Adjusted for Age, Gender, BMI)
combined_data <- cbind(diff_ko_data[, -ncol(diff_ko_data)], 
                       clinical_data[rownames(diff_ko_data), c(target_clinical, 'Age', 'Gender', 'BMI')])

# --- 3. Statistical Modeling (Spearman + GLM) ---
source('./Functions/CorrelationPlusGLM.r')

# Initialize Correlation Model
cor_model <- CorrelationModel$new(combined_data)

# Step A: Perform Spearman Correlation (Filter by p < 0.05)
# Focus on key liver fibrosis markers: ALT, Fibro_Pen (Fibrosis score), kPa (Stiffness)
spearman_results <- cor_model$SpearmanCor(p.cutoff = 0.05, 
                                          row = colnames(diff_ko_data), 
                                          col = c('ALT', 'Fibro_Pen', 'kPa'))

# Step B: GLM Adjustment (Accounting for Age & Gender)
# Ensure the associations remain significant after controlling for covariates
glm_results <- cor_model$GLM_Model(elements = c('ALT', 'Fibro_Pen', 'kPa'), 
                                   covariates = c('Age', 'Gender'))

# Step C: Filtering Spearman R values by GLM P-values
# Only retain correlations where the GLM p-value is < 0.05
glm_p_matrix <- glm_results %>% select(ends_with('.P')) %>% mutate(across(everything(), ~replace_na(., 1)))
cor_r_matrix <- spearman_results$R

# Mask non-significant GLM results (set R to 0)
filtered_r_matrix <- cor_r_matrix
filtered_r_matrix[glm_p_matrix >= 0.05] <- 0

# --- 4. Filtering for Carbohydrate & Glycan Metabolism ---
# Focus analysis on specific metabolic pathways of interest
target_pathways <- read.table(paste0(input_dir, 'Metabolism_Class_Mapping.csv'), sep = ',', header = 1, row.names = 1)

# Intersection of targeted pathways and significant correlations
final_cor_data <- filtered_r_matrix[intersect(rownames(target_pathways), rownames(filtered_r_matrix)), ]
final_cor_data <- final_cor_data[rowSums(final_cor_data) != 0, ] # Remove all-zero rows

# Prepare data for Alluvial Plot (Melt and Clean)
final_cor_data$Class <- target_pathways[rownames(final_cor_data), 'Class']
final_cor_data$KO <- rownames(final_cor_data)
plot_df <- subset(melt(final_cor_data, id.vars = c("KO", "Class")), value != 0)

# Set plotting attributes
plot_df$Abs_R <- abs(plot_df$value)
plot_df$Direction <- ifelse(plot_df$value > 0, 'Positive', 'Negative')

# --- 5. Visualization: Figure 3C (Alluvial Plot) ---



# Color Palette for various metabolic classes
metabolic_colors <- c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf')

pdf(file = paste0(output_dir, "Figure3C_KO_Clinical_Alluvial.pdf"), width = 10, height = 8)
ggplot(data = plot_df, aes(axis1 = KO, axis2 = variable, y = Abs_R)) +
  geom_alluvium(aes(fill = Direction), alpha = 0.7) +
  geom_stratum(aes(fill = Class), color = "grey") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2) +
  scale_x_discrete(limits = c("Microbial KO", "Clinical Index"), expand = c(0.15, 0.05)) +
  scale_fill_manual(values = c("Positive" = "#E64B35", "Negative" = "#4DBBD5", metabolic_colors)) +
  theme_minimal() +
  labs(title = "Associations between Carbohydrate/Glycan KOs and Liver Status",
       y = "Correlation Strength (Abs Spearman R)") +
  theme(panel.grid = element_blank())
dev.off()