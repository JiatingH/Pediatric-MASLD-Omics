################################################################################
# Project: NAFLD/MASLD Omics Analysis
# Task: Correlation Analysis (Spearman) + GLM Adjustment (Age, Gender, BMI)
# Visualization: Multi-panel Heatmap (ComplexHeatmap)
################################################################################

# --- 1. Load Required Libraries ---
library(readxl)      # Excel data import
library(gdata)       # Data manipulation
library(reshape2)    # Data reshaping
library(Hmisc)       # Statistical functions
library(tableone)    # Summary statistics
library(ggpubr)      # Publication ready plots
library(stringr)     # String manipulation
library(dplyr)       # Data wrangling
library(tidyr)       # Data tidying
library(RColorBrewer)# Color palettes
library(viridis)     # Scientific color maps
library(matrixStats) # Row-wise statistics
library(vegan)       # Ecological statistics
library(corrplot)    # Correlation visualization
library(ade4)        # Multivariate analysis
library(readr)       # Fast data reading
library(ComplexHeatmap) # Advanced heatmaps
library(circlize)    # Circular visualization and color mapping

# --- 2. Data Loading ---

# Load significantly differential species (previously identified)
# Format: Rows = Samples, Columns = Species
diff.species <- read.table('[Project_Path]/Temp_Data/Diff.Species.Table.csv', 
                           sep = ',', header = 1, row.names = 1, check.names = FALSE)
species.name <- colnames(diff.species)

# Load clinical metadata
clinical <- read.table('[Project_Path]/Data/All.Clinical.Original.Add.Class.20240806.csv', 
                       sep = ',', header = 1, row.names = 1)

# --- 3. Clinical Index Categorization ---

# Grouping clinical parameters for structured analysis and annotation
blood.index <- c('WBC', 'Lymphocyte_Count', 'Neutrophil_Count', 'Monocyte_Count', 
                 'Eosinophil_Count', 'Basophil_Count', 'HCT', 'RBC', 'Hemoglobin', 'Platelet_Count')
kiney.index <- c('ALB', 'GLB')
liver.index <- c('ALT', 'AST', 'kPa', 'Fibro_Pen', 'AKP', 'TBIL', 'DBIL', 'IBIL', 'GGT', 'ADA')
lipid.index <- c('TG', 'TC', 'HDL', 'LDL', 'non.HDL')
glucose.index <- c('Glu0','Glu30','Glu60','Glu120','Insulin0','Insulin30','Insulin60','Insulin120',
                  'Cpeptide0','Cpeptide30','Cpeptide60','Cpeptide120','HbA1c','HOMA_IR','GluAUC','CpAUC','InsAUC')

# Final selection of clinical indices for correlation heatmap
cor.index <- c(blood.index, kiney.index, 
               'Glu0','Glu120','HbA1c','HOMA_IR','GluAUC','CpAUC','InsAUC',
               lipid.index, 'AKP', 'TBIL', 'DBIL', 'IBIL', 'GGT', 'ADA', 'AST','ALT','Fibro_Pen','kPa')

# Align clinical data rows with species data rows and select target columns + covariates
taxonomy.adjusted.cor.clinical <- clinical[rownames(diff.species), 
                                           c(cor.index, 'Age', 'Gender', 'BMI_Zscore', 'BMI')]

# Combine species and clinical data into one dataframe for analysis
diff.species.clinical <- cbind(diff.species, taxonomy.adjusted.cor.clinical)

# --- 4. Correlation and GLM Modeling ---

# Load custom utility scripts
source('[Script_Path]/CorrelationPlusGLM.r') # Contains CorrelationModel class
source('[Script_Path]/MergeDataFrame.r')     # Contains merge_multiple_dataframes function

# Initialize the correlation model
species.clinical.model <- CorrelationModel$new(diff.species.clinical)

# Calculate standard Spearman correlations (Unadjusted)
# Note: Excluding 'kPa' from the column list here as it's often handled separately in GLM
species.clinical.cor <- species.clinical.model$SpearmanCor(p.cutoff = 0.05, 
                                                          row = species.name, 
                                                          col = setdiff(cor.index, c('kPa')))

# --- 5. GLM Adjustment (Confounder Control) ---

# Model A: Adjust for Age + Gender
# Calculating residuals/adjusted correlations for key liver markers
species.clinical.alt.age.gender.glm <- species.clinical.model$GLM_Model(elements = c('ALT'), covariates = c("Age", "Gender"))
species.clinical.fibropen.age.gender.glm <- species.clinical.model$GLM_Model(elements = c('Fibro_Pen'), covariates = c("Age", "Gender"))
species.clinical.kpa.age.gender.glm <- species.clinical.model$GLM_Model(elements = c('kPa'), covariates = c("Age", "Gender"))

# Model B: Adjust for Age + BMI_Zscore
species.clinical.alt.age.bmi.glm <- species.clinical.model$GLM_Model(elements = c('ALT'), covariates = c('Age','BMI_Zscore'))
species.clinical.fibropen.age.bmi.glm <- species.clinical.model$GLM_Model(elements = c('Fibro_Pen'), covariates = c('Age','BMI_Zscore'))
species.clinical.kpa.age.bmi.glm <- species.clinical.model$GLM_Model(elements = c('kPa'), covariates = c('Age','BMI_Zscore'))

# --- 6. Extraction and Merging ---

# Extraction logic: Extract correlation coefficient (.Cor) and P-value (.P)
# Renaming columns to reflect the adjustment model used (Age_Gender vs Age_BMIZscore)

# [Example of extraction for ALT - Logic repeated for FibroPen and kPa in the full merge]
species.glm.alt.age.gender.cor <- species.clinical.alt.age.gender.glm %>% dplyr::select(ends_with(".Cor")) %>% rename_with(~ gsub("Cor", "Age_Gender", .x))
species.glm.alt.age.gender.p <- species.clinical.alt.age.gender.glm %>% dplyr::select(ends_with('.P')) %>% rename_with(~ gsub("\\.P", ".Age_Gender", .x))

# Merging all Correlation coefficients into one master table
# species.clinical.R and species.clinical.P are assumed to be output objects from $SpearmanCor()
merge.species.cor <- merge_multiple_dataframes(list(
    species.clinical.R, 
    species.glm.alt.age.gender.cor, species.glm.alt.age.bmi.cor,
    species.glm.fibropen.age.gender.cor, species.glm.fibropen.age.bmi.cor,
    species.glm.kpa.age.gender.cor, species.glm.kpa.age.bmi.cor
)) %>% mutate(across(everything(), ~ifelse(is.na(.), 0, .x)))

# Merging all P-values into one master table
merge.species.p <- merge_multiple_dataframes(list(
    species.clinical.P, 
    species.glm.alt.age.gender.p, species.glm.alt.age.bmi.p,
    species.glm.fibropen.age.gender.p, species.glm.fibropen.age.bmi.p,
    species.glm.kpa.age.gender.p, species.glm.kpa.age.bmi.p
)) %>% mutate(across(everything(), ~ifelse(is.na(.), 1, .x)))

# --- 7. Heatmap Preparation ---

# Split merged correlation matrix for panel-wise plotting
merge.species.cor1 <- merge.species.cor[, cor.index] # Main correlation panel
merge.species.cor2.1 <- merge.species.cor %>% dplyr::select(starts_with('ALT.'))
merge.species.cor3.1 <- merge.species.cor %>% dplyr::select(starts_with('Fibro_Pen.'))
merge.species.cor4.1 <- merge.species.cor %>% dplyr::select(starts_with('kPa.'))

# Prepare Top Annotation (Clinical Categories)
Clinical.Class <- c(rep('CBC',10), rep('Protein',2), 
                    rep('Glucose Homeostasis',7), rep('Lipid Metabolism', 5), 
                    rep('Liver Function',10))
annotation1 <- data.frame(Class = Clinical.Class)

# Prepare Right Annotation (Microbial Taxonomy)
taxonomy.class <- read.table('[Project_Path]/Data/Taxonomy.Class.Table.csv', 
                             sep = ',', header = 1, row.names = 1)
diff.species.class <- filter(taxonomy.class, Species %in% rownames(merge.species.cor1))
annotation.class <- diff.species.class[, c('Phylum'), drop = FALSE]
rownames(annotation.class) <- rownames(merge.species.cor1)

# Define Color Mappings
col1 <- colorRamp2(c(-0.45, 0, 0.45), c("#4DBBD5CC", "white", "#E64B35CC")) # Spearman
col2 <- colorRamp2(c(-1, 0, 1), c('#B692B9', '#E0E0E0', '#4DAF4A'))         # GLM Panels

# --- 8. Heatmap Generation ---

# Panel 1: Main Spearman Correlations
ht1 = Heatmap(as.matrix(merge.species.cor1), 
              name = "Spearman Correlation", col = col1, 
              cluster_columns = TRUE, border = TRUE,
              rect_gp = gpar(col = "white", lwd = 1.2),
              top_annotation = HeatmapAnnotation(df = annotation1,
                  col = list(Class = c('CBC'='#1B9E77','Protein'='#D95F02','Glucose Homeostasis'='#7570B3',
                                     'Lipid Metabolism'='#E7298A','Liver Function'='#66A61E'))),
              right_annotation = rowAnnotation(df = annotation.class,
                  col = list(Phylum = c('Actinobacteria'='#66C2A5', 'Bacteroidetes'='#FC8D62', 'Firmicutes'='#8DA0CB')))
)

# Panel 2: ALT Adjusted Results
ht2 = Heatmap(as.matrix(merge.species.cor2.1), name = 'ALT GLM', col = col2, 
              cluster_columns = FALSE, width = unit(1, "cm"), border = TRUE,
              rect_gp = gpar(col = "white", lwd = 1.2))

# Panel 3: FibroPen Adjusted Results
ht3 = Heatmap(as.matrix(merge.species.cor3.1), name = 'Fibro_Pen GLM', col = col2, 
              cluster_columns = FALSE, width = unit(1, "cm"), border = TRUE,
              rect_gp = gpar(col = "white", lwd = 1.2))

# Panel 4: kPa Adjusted Results
ht4 = Heatmap(as.matrix(merge.species.cor4.1), name = 'kPa GLM', col = col2, 
              cluster_columns = FALSE, width = unit(1, "cm"), border = TRUE,
              rect_gp = gpar(col = "white", lwd = 1.2))

# --- 9. Export Final Figure ---

pdf(file = '[Project_Path]/Results/Species.Clinical.ComplexHeatmap.pdf', width = 17, height = 12)
draw(ht1 + ht2 + ht3 + ht4) # Combine panels side-by-side
dev.off()

circos.clear() # Reset circular plotting parameters