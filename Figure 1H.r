################################################################################
# Project: NAFLD/MASLD Omics Analysis
# Task: Differential Abundance & Sensitivity Analysis (Figure 1H)
# Logic: 1. KW-test -> 2. Binary Weighted Sorting -> 3. Random Subsampling (100x)
################################################################################

# Load statistical and machine learning libraries
library(randomForest)
library(datasets)
library(caret)
library(pROC)

# Load data manipulation and visualization libraries
library(readxl)      # Read Excel files
library(gdata)
library(reshape2)    # Reshape data
library(Hmisc)
library(pheatmap)    # Heatmap visualization
library(tableone)    # Baseline characteristics
library(ggpubr)      # Publication ready plots
library(stringr)

# Load color palettes and data wrangling tools
library(RColorBrewer)
library(viridis)
library(dplyr)       # Data manipulation

# Load specialized bioinformatics and statistical tools
library(matrixStats) # Row statistics
library(vegan)       # Ecology/Community analysis
library('tidyr')
library('KEGGREST')

# Load matching and imputation libraries
library(Matching)
library(tableone)
library(mice)
library(MatchIt)

# --- Path Configuration (Placeholders for Public Repository) ---
# Replace [Project_Data_Path] with your local data directory or project root
DATA_DIR <- "[Project_Data_Path]/Data"
SCRIPT_DIR <- "[Project_Script_Path]/Scripts"
OUTPUT_DIR <- "[Project_Output_Path]/Results/Sensitivity"

# --- Data Loading and Pre-processing ---

# Load Species abundance data (Example: Species level)
species = read.table(file.path(DATA_DIR, 'All.Species.csv'), sep =',', header = 1, row.names = 1, check.names = FALSE)

# Update group naming from NAFLD to MASLD and set factor levels
species$Group = ifelse(species$Group == 'NAFLD', 'MASLD', species$Group)
species$Group = factor(species$Group, levels = c('HC','OB','MASLD'))

# Load taxonomy classification table
taxonomy.class = read.table(file.path(DATA_DIR, 'Taxonomy.Class.Table.csv'), sep = ',', header = 1, row.names = 1, check.names = FALSE)

# Import custom statistical functions
source(file.path(SCRIPT_DIR, 'KW.r'))         # Kruskal-Wallis Function
source(file.path(SCRIPT_DIR, 'Wilcox.r'))     # Wilcoxon Function
source(file.path(SCRIPT_DIR, 'FoldChange.r')) # Fold Change Function

# Extract species elements for comparison
species.data = species[setdiff(colnames(species), 'Group')]
species.name = colnames(species.data)

# Perform Kruskal-Wallis test across the three groups
species.kw = KW(data = species, compare = species.name, class = 'Group')
# Filter significant species (FDR < 0.05)
sig.species.name = rownames(filter(species.kw, fdr < 0.05))

# --- Custom Sorting Logic ---

# Define function to calculate FoldChange and P-values, then sort by significance weights
FoldChange_Order <- function(data, sig, group) {
     sig.data <- data[, c(sig, group)]
    
     # Calculate Fold Change for pairwise comparisons
     hc.ob.fd <- FD(sig.data, com = sig, class = group, 'HC', 'OB')
     ob.masld.fd <- FD(sig.data, com = sig, class = group, 'OB', 'MASLD')
     hc.masld.fd <- FD(sig.data, com = sig, class = group, 'HC', 'MASLD')
    
     # Calculate P-values for pairwise comparisons
     hc.ob.p <- Wilcox(sig.data, com = sig, class = group, 'HC', 'OB')
     ob.masld.p <- Wilcox(sig.data, com = sig, class = group, 'OB', 'MASLD')
     hc.masld.p <- Wilcox(sig.data, com = sig, class = group, 'HC', 'MASLD')
    
     # Merge Fold Change results
     fd <- cbind(hc.ob.fd, ob.masld.fd, hc.masld.fd)
     colnames(fd) <- c('HC.OB', 'OB.MASLD', 'HC.MASLD')
    
     # Merge P-value results
     p <- data.frame(HC.OB = hc.ob.p$p, OB.MASLD = ob.masld.p$p, HC.MASLD = hc.masld.p$p)
     rownames(p) <- sig
    
     # Sort by significance binary values
     # Convert P-values to binary (-1 for non-significant, 1 for significant)
     species_p_value <- p %>%
          mutate_all(~ifelse(. < 0.05, 1, ifelse(. > 0.05, -1, 0)))
    
     # Convert Fold Change to binary (-1 for decrease, 1 for increase)
     species_fd_value <- fd %>%
          mutate_all(~ifelse(. < 0, 1, ifelse(. > 0, -1, 0)))
    
     # Apply custom weights to Fold Change values for sorting priority
     species_fd_value$HC.OB <- ifelse(species_fd_value$HC.OB == -1, -2, ifelse(species_fd_value$HC.OB == 1, 2, species_fd_value$HC.OB))
     species_fd_value$HC.MASLD <- ifelse(species_fd_value$HC.MASLD == -1, -2, ifelse(species_fd_value$HC.MASLD == 1, 2, species_fd_value$HC.MASLD))
     species_fd_value$OB.MASLD <- ifelse(species_fd_value$OB.MASLD == -1, -3, ifelse(species_fd_value$OB.MASLD == 1, 3, species_fd_value$OB.MASLD))
    
     # Apply custom weights to P-values for sorting priority
     species_p_value$HC.OB <- ifelse(species_p_value$HC.OB == 1, 2, species_p_value$HC.OB)
     species_p_value$HC.MASLD <- ifelse(species_p_value$HC.MASLD == 1, 2, species_p_value$HC.MASLD)
     species_p_value$OB.MASLD <- ifelse(species_p_value$OB.MASLD == 1, 3, species_p_value$OB.MASLD)
    
     # Calculate row sums for sorting rank
     species_fd_value$Sum <- rowSums(species_fd_value[, c(1:3)])
     species_p_value$Sum <- rowSums(species_p_value[, c(1:3)])
    
     # Sort data frame based on weighted sums
     species_fd_order <- species_fd_value[order(species_p_value$Sum, species_fd_value$Sum, abs(species_fd_value$HC.MASLD), decreasing = TRUE),]
    
     # Reorder raw FD and P-value data frames based on calculated order
     species_fd_plot <- fd[rownames(species_fd_order),]
     species_p_plot <- p[rownames(species_fd_order),]
    
     return(list(fd_plot = species_fd_plot,    
                                   p_plot = species_p_plot,
                                   fd_value = species_fd_value,    
                                   p_value = species_p_value))
}

# Execute sorting function
species.fold.change.order = FoldChange_Order(species, sig.species.name, 'Group')

# Extract plotting data
species.p = species.fold.change.order$p_plot
species.fd = species.fold.change.order$fd_plot

# --- Taxonomic Class Sorting ---

# Reorder results based on Phylum assignment ('Bacteroidetes', 'Firmicutes', 'Actinobacteria')
species.taxo.class = filter(taxonomy.class, Species %in% rownames(species.p)) %>% distinct(Species, .keep_all = TRUE)
rownames(species.taxo.class) = species.taxo.class$Species

species.p$Phylum = species.taxo.class[rownames(species.p), 'Phylum']
species.p$Phylum = factor(species.p$Phylum, levels = c('Bacteroidetes','Firmicutes','Actinobacteria'))
species.p.order = species.p[order(species.p$Phylum),]

species.fd$Phylum = species.taxo.class[rownames(species.fd), 'Phylum']
species.fd$Phylum = factor(species.fd$Phylum, levels = c('Bacteroidetes','Firmicutes','Actinobacteria'))
species.fd.order = species.fd[order(species.fd$Phylum),]

# --- Visualization: Figure 1H Left (Heatmap) ---

# Define color scale breaks
bk <- c(seq(-1, 1, by = 0.01))
BK <- bk[!duplicated(bk)]

# Generate significance matrix (stars for p < 0.05)
species.display_numbers <- matrix(ifelse(species.p.order < 0.05, '*', ''), nrow = nrow(species.p.order), ncol = ncol(species.p.order))
    
# Draw Heatmap using pheatmap
pheatmap(
          species.order[c(1,2,3)],
          cluster_row = FALSE,
          cluster_col = FALSE,
          cellwidth = 20,
          cellheight = 20,
          color = c(colorRampPalette(colors = c("#4DBBD5CC", "white"))(100),
                                   colorRampPalette(colors = c("white", "#E64B35CC"))(100)),
          legend_breaks = seq(-1, 1, 1),
          breaks = BK,
          display_numbers = species.display_numbers,
          fontsize_number = 10,
          fontsize = 10,
          annotation_row = species.order[, 4, drop = FALSE],
          filename = file.path(OUTPUT_DIR, 'My.Diff.Species.pdf')
)


# --- Sensitivity Analysis: Figure 1H Right ---

# Define function for Random Subsampling Kruskal-Wallis test (100 iterations)
Random_KW <- function(df) {
     if (!("Group" %in% colnames(df))) {
          stop("The input dataframe must have a 'Group' column.")
     }
    
     required_groups <- c("HC", "OB", "NAFLD")
     if (!all(required_groups %in% unique(df$Group))) {
          stop("The 'Group' column must contain 'HC', 'OB', and 'NAFLD'.")
     }
    
     result.list = list()
     set.seed(123) # Ensure reproducibility
    
     # Repeat subsampling 100 times
     for (i in 1:100) {
          # Randomly select 55 NAFLD samples
          NAFLD_samples <- df[df$Group == "NAFLD", ]
          selected_NAFLD <- NAFLD_samples[sample(nrow(NAFLD_samples), 55), ]
         
          # Combine with original HC and OB samples
          other_samples <- df[df$Group %in% c("HC", "OB"), ]
          combined_samples <- rbind(selected_NAFLD, other_samples)
         
          # Execute KW test for each feature
          result = KW(combined_samples, compare = setdiff(colnames(combined_samples),'Group'), class = 'Group', adjusted_method = 'BH')
         
          colnames(result) = c(paste0('P_', i), paste0('FDR_', i))
          result.list[[i]] = result    
     }
     
     # Merge results and count significance occurrences
     sig.df = merge_multiple_dataframes(result.list)
     
     p.sig.df = sig.df[, grep("^P_", names(sig.df))]
     fdr.sig.df = sig.df[, grep('^FDR_', names(sig.df))]
     
     p.count = apply(p.sig.df, 1, function(x) sum(x < 0.05))
     fdr.count = apply(fdr.sig.df, 1, function(x) sum(x < 0.05))
     
     sig.df$P.Count = p.count
     sig.df$FDR.Count = fdr.count
     
     return(sig.df)
}

# Define function for Random Subsampling Wilcoxon test (100 iterations)
Random_Wilcox <- function(df, g1, g2) {
     if (!("Group" %in% colnames(df))) {
          stop("The input dataframe must have a 'Group' column.")
     }
    
     required_groups <- c("HC", "OB", "NAFLD")
     if (!all(required_groups %in% unique(df$Group))) {
          stop("The 'Group' column must contain 'HC', 'OB', and 'NAFLD'.")
     }
    
     result.list <- list()
     set.seed(123) 
    
     for (i in 1:100) {
          NAFLD_samples <- df[df$Group == "NAFLD", ]
          if (nrow(NAFLD_samples) < 55) {
               stop("Not enough NAFLD samples to randomly select 55.")
          }
         
          selected_NAFLD <- NAFLD_samples[sample(nrow(NAFLD_samples), 55), ]
          other_samples <- df[df$Group %in% c("HC", "OB"), ]
          combined_samples <- rbind(selected_NAFLD, other_samples)
         
          # Pairwise comparison via Wilcoxon test
          result <- Wilcox(combined_samples,    
                                                      com = setdiff(colnames(combined_samples), 'Group'),    
                                                      class = 'Group',    
                                                      g1 = g1,    
                                                      g2 = g2,    
                                                      adjusted_method = 'BH')
         
          colnames(result) <- c(sprintf("P_%d", i), sprintf("FDR_%d", i))
          result.list[[i]] <- result
     }
    
     # Merge results and count significance
     sig.df <- merge_multiple_dataframes(result.list)
    
     p.sig.df <- sig.df[, grep("^P_", names(sig.df))]
     fdr.sig.df <- sig.df[, grep("^FDR_", names(sig.df))]
    
     p.count <- apply(p.sig.df, 1, function(x) sum(x < 0.05))
     fdr.count <- apply(fdr.sig.df, 1, function(x) sum(x < 0.05))
    
     sig.df$P.Count <- p.count
     sig.df$FDR.Count <- fdr.count
    
     return(sig.df)
}

# Execute Random Analysis
my.species.Random.KW = Random_KW(my.species.data)
my.species.Random.Wilcox.HC.MASLD = Random_Wilcox(my.species.data, g1 = 'HC', g2 = 'NAFLD')
my.species.Random.Wilcox.OB.MASLD = Random_Wilcox(my.species.data, g1 = 'OB', g2 = 'NAFLD')

# Extract significant counts from random iterations
my.species.Random.count = my.species.Random.KW[c('P.Count','FDR.Count')]
my.species.Random.Wilcox.HC.MASLD.count = my.species.Random.Wilcox.HC.MASLD[rownames(species.fd.order), c('P.Count','FDR.Count')]
my.species.Random.Wilcox.OB.MASLD.count = my.species.Random.Wilcox.OB.MASLD[rownames(species.fd.order), c('P.Count','FDR.Count')]

my.species.Random.KW.count = my.species.Random.count[rownames(species.fd.order),]

# Set specific column names for merging
colnames(my.species.Random.KW.count) = paste0(colnames(my.species.Random.KW.count), '_KW')
colnames(my.species.Random.Wilcox.HC.MASLD.count) = paste0(colnames(my.species.Random.Wilcox.HC.MASLD.count), '_HC.MASLD_Wilcox')
colnames(my.species.Random.Wilcox.OB.MASLD.count) = paste0(colnames(my.species.Random.Wilcox.OB.MASLD.count), '_OB.MASLD_Wilcox')

# Combine summary results
species.Random.results.summary = cbind(my.species.Random.KW.count, my.species.Random.Wilcox.HC.MASLD.count, my.species.Random.Wilcox.OB.MASLD.count)

# Normalize counts to 100%
species.Random.result.data = species.Random.results.summary[, c(2, 3, 5)]/100

# Visualization: Figure 1H Right (Stability Pie Plot)
pdf(file = file.path(OUTPUT_DIR, 'My.Diff.Species.Random.pdf'))
corrplot(corr = as.matrix(species.Random.result.data), method = 'pie',
         col = viridis(50), cl.lim = c(-1, 1),
         tl.col = "black",
         tl.cex = 0.8 )
dev.off()