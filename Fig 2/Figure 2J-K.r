################################################################################
# Project: NAFLD/MASLD Microbial Meta-Analysis
# Task: Cross-Cohort Network Integration (Venn & Prevalence Analysis)
# Aim: Identify robust microbial "signatures" conserved across different populations
################################################################################

# --- 1. Load Required Libraries ---
library(readxl); library(gdata); library(reshape2); library(Hmisc)
library(tableone); library(ggpubr); library(stringr); library(dplyr); library(tidyr)
library(RColorBrewer); library(viridis); library(matrixStats); library(vegan)
library(corrplot); library(ade4); library(readr); library(VennDiagram)

# --- 2. Data Import ---
# Loading SparCC results from different cohorts (formatted as Var1, Var2, cor, p)
# Set your public base path here
base_path <- "./Public_Project/SparCC_Results/"

# Cohort A: China Adult Study (Different stages)
adult_china.NAFL = read.table(paste0(base_path, 'china_adult_NAFL.csv'), sep = ',', header = 1)
adult_china.All  = read.table(paste0(base_path, 'china_adult_All.csv'), sep = ',', header = 1)

# Cohort B: USA Adult Study
adult_usa.All    = read.table(paste0(base_path, 'usa_adult_All.csv'), sep = ',', header = 1)

# Cohort C: China 
adult_china688.NAFLD = read.table(paste0(base_path, 'china_688_NAFLD.csv'), sep = ',', header = 1)

# Cohort D: UK Pediatric Study
child_uk.NAFLD   = read.table(paste0(base_path, 'uk_child_MASLD.csv'), sep = ',', header = 1)

# Internal Cohort
masld_internal   = read.table(paste0(base_path, 'internal_MASLD.csv'), sep = ',', header = 1)

# --- 3. Filtering and Preprocessing ---

# Retain only significant correlations (p < 0.05)
hc.sig    = filter(masld_internal, p < 0.05)
# ... [Repeat filtering for all cohorts] ...

# Remove redundant/mirrored rows for undirected network (A-B is the same as B-A)
# We sort Var1 and Var2 alphabetically to identify duplicates
remove_mirrored_edges <- function(df) {
  df[!duplicated(t(apply(df[, c("Var1", "Var2")], 1, sort))), ]
}

masld.unique = remove_mirrored_edges(filter(masld_internal, p < 0.05))
uk.unique    = remove_mirrored_edges(filter(child_uk.NAFLD, p < 0.05))
# ... [Apply to all datasets] ...

# Standardize species names and calculate absolute correlation (ABS)
standardize_data <- function(data){
  data %>% mutate(
    Var1 = gsub("\\.", " ", Var1),
    Var2 = gsub("\\.", " ", Var2),
    ABS = abs(cor)
  )
}

# --- 4. Cross-Cohort Intersection Function ---

calculate_intersections_and_plot <- function(dataframes, cols = c("Var1", "Var2")) {
  
  # Validate presence of required columns
  validate_dfs <- function(df) {
    if (is.null(df) || nrow(df) == 0 || !all(cols %in% colnames(df))) return(NULL)
    return(df)
  }
  
  # Create a unique ID for each interaction: "SpeciesA_SpeciesB" (alphabetical)
  standardize_pairs <- function(df) {
    apply(df[cols], 1, function(row) paste(sort(row), collapse = "_"))
  }
  
  valid_dfs <- Filter(Negate(is.null), lapply(dataframes, validate_dfs))
  std_sets  <- lapply(valid_dfs, standardize_pairs)
  names(std_sets) <- paste0("Cohort_", seq_along(std_sets))
  
  # Plot Venn Diagram
  colors <- brewer.pal(min(length(std_sets), 12), "Set3")
  venn_plot <- venn.diagram(
    x = std_sets, filename = NULL, main = "Conserved Interactions",
    fill = colors, alpha = 0.5, cex = 1.5, fontface = "bold"
  )
  
  # Statistics: How many cohorts share each interaction?
  counts <- as.data.frame(table(unlist(std_sets)))
  colnames(counts) <- c("Combination", "Count")
  counts <- counts %>% arrange(desc(Count))
  
  return(list(venn_plot = venn_plot, combination_counts = counts))
}

# Execute Analysis across 5 major datasets
result = calculate_intersections_and_plot(list(
  standardize_data(masld.unique), 
  standardize_data(uk.unique), 
  standardize_data(adult_usa.All.unique),
  standardize_data(adult_china.All.unique),
  standardize_data(adult_china688.NAFLD.unique)
))

# --- 5. Final Visualization ---

# Save Venn Diagram (Figure 2L)
pdf('./Results/Venn_Overlap_Network.pdf')
grid.draw(result$venn_plot)
dev.off()

# Plot Prevalence Barplot (Figure 2M)
# Focus on interactions appearing in >3 datasets
core_interactions = filter(result$combination_counts, Count > 3)

pdf('./Results/Core_Interaction_Prevalence.pdf')
ggplot(core_interactions, aes(x = reorder(Combination, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Microbial Interaction (Species Pair)", y = "Cohort Prevalence Count")
dev.off()

# --- 6. Export for Cytoscape ---
# Extract full edge data for interactions appearing in multiple datasets
# This CSV can be imported into Cytoscape for final network visualization
core_edge_list <- standardize_data(masld.unique) %>% 
  filter(apply(.[c("Var1", "Var2")], 1, function(x) paste(sort(x), collapse = "_")) %in% core_interactions$Combination)

write.csv(core_edge_list, './Results/Conserved_Network_Edges.csv', row.names = FALSE)