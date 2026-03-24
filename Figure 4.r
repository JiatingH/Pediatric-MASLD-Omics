################################################################################
# Project: MASLD Metagenomics - CAZyme Functional Analysis
# Task: Differential CAZyme Families and Taxonomic Contribution (Figure 4)
################################################################################

# --- 1. Load Libraries ---
library(randomForest); library(datasets); library(caret); library(pROC)
library(readxl); library(gdata); library(reshape2); library(Hmisc)
library(pheatmap); library(tableone); library(ggpubr); library(stringr)
library(RColorBrewer); library(viridis); library(dplyr); library(matrixStats)
library(vegan); library(tidyr); library(KEGGREST)
library(ComplexHeatmap); library(circlize)

# --- 2. Data Import & Initial Processing ---
# Standardize paths to public/relative directories
data_dir <- "./Project_Data/CAZyme_Analysis/"
result_dir <- "./Project_Results/Figures/"
script_dir <- "./Project_Scripts/Functions/"

# Load Sample Info and KO Abundance
sampleinfo = read.table(paste0(data_dir, 'metagenomics_sample_info.csv'), sep = ',', header = 1, check.names = FALSE)
ko.temp = read.table(paste0(data_dir, 'all_ko_abundance_temp.csv'), sep = ',', header = 1, row.names = 1, check.names = FALSE)
ko = read.table(paste0(data_dir, 'all_ko_abundance.csv'), sep = ',', header = 1, row.names = 1)

# Extract EC Numbers from KO descriptions
ko.class.data <- ko.class %>%
  mutate(EC_Numbers = sapply(str_extract_all(Description, "EC:[0-9.]+"), 
                             function(x) paste(x, collapse = ", ")))

ko.class.data$EC_Numbers = gsub(':', " ", ko.class.data$EC_Numbers)
ko.class.data$KO = rownames(ko.class.data)

# Load CAZyme annotation data
old.cazyme = read.csv(paste0(data_dir, 'CAZyme_V1.txt'), sep = '\t', header = 1, check.names = FALSE, row.names = 1)
new.cazyme = read.csv(paste0(data_dir, 'CAZyme_V2.txt'), sep = '\t', header = 1, check.names = FALSE, row.names = 1)

# Merge CAZyme datasets using custom function
source(paste0(script_dir, 'MergeDataFrame.r'))
all.cazyme = merge_dataframes(old.cazyme, new.cazyme)
all.cazyme.class = all.cazyme[c(1,2)]

# Synchronize CAZyme data with sample metadata
rownames(sampleinfo) = sampleinfo$MicrobioID
all.cazyme.sampleinfo = sampleinfo[intersect(colnames(all.cazyme), sampleinfo$MicrobioID),]
all.cazyme.data = all.cazyme[all.cazyme.sampleinfo$MicrobioID]
colnames(all.cazyme.data) = all.cazyme.sampleinfo$Hospital.Num

# Re-aligning IDs and handling zero values
intersect.id = intersect(sampleinfo$MicrobioID, colnames(all.cazyme))
all.cazyme.data = t(all.cazyme[intersect.id]) %>% as.data.frame
all.cazyme.data[all.cazyme.data == 0] = NA

# --- 3. Data Cleaning and Statistics ---
source(paste0(script_dir, 'DataMining.r'))
source(paste0(script_dir, 'KW.r'))

# Remove low-prevalence features and fill missing values (Half-minimum method)
all.cazyme.data.removeNA = Remove_HighNA(all.cazyme.data, axis = 'col', threshold = 0.2)
all.cazyme.data.fillNA = Fill_MissingValues(all.cazyme.data.removeNA, method = 'half_min')

# Calculate Sums for CAZyme Classes (AA, CBM, CE, GH, GT, PL)
cazyme.class.sum = data.frame(
  AA = rowSums(all.cazyme.data.fillNA[grep('AA', colnames(all.cazyme.data.fillNA))]),
  CBM = rowSums(all.cazyme.data.fillNA[grep('CBM', colnames(all.cazyme.data.fillNA))]),
  CE = rowSums(all.cazyme.data.fillNA[grep('CE', colnames(all.cazyme.data.fillNA))]),
  GH = rowSums(all.cazyme.data.fillNA[grep('GH', colnames(all.cazyme.data.fillNA))]),
  GT = rowSums(all.cazyme.data.fillNA[grep('GT', colnames(all.cazyme.data.fillNA))]),
  PL = rowSums(all.cazyme.data.fillNA[grep('PL', colnames(all.cazyme.data.fillNA))])
)

# Define Groups: HC, OB, MASLD
cazyme.class.sum$Group = cazyme.sample.info$ClassNote
cazyme.class.sum$Group = ifelse(cazyme.class.sum$Group %in% c('NASH','NAFL'), 'MASLD', cazyme.class.sum$Group)

# --- 4. Differential Abundance Analysis (Figure 4A-F) ---
source(paste0(script_dir, 'BoxPlot.r'))

MultipleBoxPlot(cazyme.class.sum, x = 'Group', y = setdiff(colnames(cazyme.class.sum), 'Group'),
                group_elements = c('HC','OB','MASLD'),
                color_values = c('#CD6155','#040676','#5F93CE'),
                pdf_path = paste0(result_dir, 'Figure4_CAZyme_Class_Abundance.pdf'))

# Family-level comparison with Bonferroni correction
all.cazyme.data.fillNA$Group = cazyme.sample.info$ClassNote
all.cazyme.data.fillNA$Group = ifelse(all.cazyme.data.fillNA$Group %in% c('NAFL','NASH'), 'MASLD', all.cazyme.data.fillNA$Group)

all.cazyme.pval.bonfer = KW(all.cazyme.data.fillNA, compare = setdiff(colnames(all.cazyme.data.fillNA),'Group'), class = 'Group')
all.cazyme.fdr.bonfer = filter(all.cazyme.pval.bonfer, fdr < 0.05)

# --- 5. Mapping CAZyme Families to KOs via EC Numbers ---
diff.cazyme.bonfer = all.cazyme[rownames(all.cazyme.fdr.bonfer), c('Family','Activities_in_Family')]

# Extract EC numbers from significant families
diff.cazyme.bonfer.data <- diff.cazyme.bonfer %>%
  mutate(EC_Numbers = sapply(str_extract_all(Activities_in_Family, "EC [0-9.]+"), 
                             function(x) paste(x, collapse = ", ")))
diff.cazyme.bonfer.data$ID = rownames(diff.cazyme.bonfer.data)

# Flatten and match with KO data
split_ecs <- diff.cazyme.bonfer.data %>%
  mutate(EC_Numbers = strsplit(EC_Numbers, ", ")) %>%
  unnest(EC_Numbers)

diff.cazyme.ec.ko = filter(split_ecs, EC_Numbers %in% intersect(ko.class.data$EC_Numbers, unique(unlist(str_extract_all(diff.cazyme.bonfer.data$EC_Numbers, "EC [0-9.]+")))))

# Loop to assign KO strings back to ECs
ec.ko.vector = c()
for(n in seq_along(diff.cazyme.ec.ko$EC_Numbers)){
    current_ec = diff.cazyme.ec.ko$EC_Numbers[n]
    ec.ko = filter(ko.class.data, EC_Numbers == current_ec)
    ec.ko.data = if (nrow(ec.ko) > 0) paste(unique(ec.ko$KO), collapse = ", ") else NA
    ec.ko.vector = c(ec.ko.vector, ec.ko.data)
}
diff.cazyme.ec.ko$KO = ec.ko.vector

# Aggregating KOs by Family
diff.cazyme.ko_expanded <- diff.cazyme.ec.ko %>%
  group_by(ID) %>%
  summarize(across(everything(), ~ paste(unique(.), collapse = ", ")), .groups = 'drop') %>%
  separate_rows(KO, sep = ", ")

# --- 6. Taxonomic Contribution (Figure 4G-H) ---
source(paste0(script_dir, 'KO_Contribution.r'))

# Generate stacked bar plots for the top 5 contributing species of each CAZyme class
# (Repeating for CE, GH, GT, PL using unique KOs found in step 5)
classes_to_plot <- c("CE", "GH", "GT", "PL")
for(cls in classes_to_plot) {
  target_kos <- diff.cazyme.ko_expanded %>% filter(grepl(cls, ID)) %>% pull(KO) %>% unique()
  Pathway_Top5_Contributed_StackedBar(data = ko.temp, ko_name = target_kos, group.info = group.info, n = 500, 
                                      save_path = paste0(result_dir, cls, "_Contribution_StackedBar.pdf"))
}

# --- 7. Heatmap and Subgroup Analysis (Figure 4K) ---
# Calculate Z-scores and Fold Changes between subgroups (Simple MASLD vs. Advanced)
source(paste0(script_dir, 'Wilcox.r'))

# Heatmap setup for Glycoside Hydrolases (GH)
# Comparison: HC vs MASLD and MASLD Subgroup Comparison
diff.cazyme.zscore.mean = data.frame(HC = colMeans(hc_mat), OB = colMeans(ob_mat), MASLD = colMeans(masld_mat))
# [Calculated Fold Changes (FC) and Wilcoxon P-values for MASLD vs Advanced Subgroups]

# Heatmap Visualization with Barplot Annotations
col1 = colorRamp2(c(-0.5, 0, 0.5), viridis(3))
col2 = colorRamp2(c(0, 30, 60), plasma(3))

# Final Combination Heatmap: Abundance (left) + Taxon Contribution (right)
# Trace specific species contribution (e.g., Bacteroides vulgatus) to GH families
gh1 = Heatmap(as.matrix(diff.cazyme.zscore.mean[rownames(cazyme.class4),]), 
              name = "Z-score", col = col1, cluster_columns = FALSE)

gh2 = Heatmap(as.matrix(GH_top10_species_contribution_matrix), 
              name = "Contribution %", col = col2, cluster_columns = FALSE)

pdf(file = paste0(result_dir, 'Figure4K_GH_Contribution_Heatmap.pdf'), width = 17, height = 12)
gh1 + gh2
dev.off()