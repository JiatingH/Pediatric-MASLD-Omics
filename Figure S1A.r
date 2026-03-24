################################################################################
# Project: NAFLD Omics Analysis
# Task: Taxonomic Composition Analysis (Top 10 Species Stacked Bar/Donut Chart)
# Note: This script uses relative paths via the 'here' package for portability.
################################################################################

# Load required libraries
library(readxl)      # Read Excel files
library(gdata)       # Data manipulation tools
library(reshape2)    # Data reshaping
library(Hmisc)       # Miscellaneous data analysis functions
library(pheatmap)    # Heatmap visualization
library(tableone)    # Table 1 generation
library(ggpubr)      # Publication-ready plots
library(stringr)     # String manipulation
library(ggplot2)     # Core plotting system
library(RColorBrewer)# Color palettes
library(viridis)     # Color scales
library(dplyr)       # Data manipulation and merging
library(here)        # IMPORTANT: Enables relative paths from project root

#-------------------------------------------------------------------------------
# 1. Data Loading (Using Relative Paths)
#-------------------------------------------------------------------------------

# Load Pediatric WGS Species Data
# Assumes 'All.Species.csv' is in your project's data subfolder
my.species.data <- read.table(here("Data", "All.Species.csv"), 
                             sep = ',', header = TRUE, row.names = 1, check.names = FALSE)

# Load Taxonomy data using relative path from project root
child.uk.wgs <- read.table(here("Public_Taxonomy", "Relative_Abundance", "Children_NAFLD_WGS", "PRJNA328258_Taxonomy.csv"), 
                          sep = ',', header = TRUE)

# Transpose and convert to data frame
child.uk.wgs.data <- t(child.uk.wgs[, c(8:40)]) %>% as.data.frame()
colnames(child.uk.wgs.data) <- child.uk.wgs$Species

# Load Adult WGS Species Data (Relative paths)
adult_china_wgs.species <- read.table(here("Public_Taxonomy", "Relative_Abundance", "Adult_NAFLD_WGS", "NMDC10018157_Species.csv"), 
                                     sep = ',', header = TRUE, row.names = 1, check.names = FALSE) %>% t() %>% as.data.frame()

adult_usa_wgs.species <- read.table(here("Public_Taxonomy", "Relative_Abundance", "Adult_NAFLD_WGS", "PRJNA373901_Species.csv"), 
                                   sep = ',', header = TRUE, row.names = 1, check.names = FALSE) %>% t() %>% as.data.frame()

adult_china688_wgs.species <- read.table(here("Public_Taxonomy", "Relative_Abundance", "Adult_NAFLD_WGS", "PRJNA686835_Species.csv"), 
                                        sep = ',', header = TRUE, row.names = 1, check.names = FALSE) %>% t() %>% as.data.frame()

#-------------------------------------------------------------------------------
# 2. Function Definition
#-------------------------------------------------------------------------------

# Function to plot Top 10 taxonomic composition as a stacked bar/pie chart
plot_top10_stacked_bar <- function(data, percent_col, genus_col, top_n = 10) {
  colnames(data)[colnames(data) == percent_col] <- "Percent"
  colnames(data)[colnames(data) == genus_col] <- "Genus"
  
  data_processed <- data %>%
    arrange(desc(Percent)) %>%
    mutate(Genus = ifelse(row_number() > top_n, "others", Genus)) %>%
    group_by(Genus) %>%
    summarise(Percent = sum(Percent)) %>%
    ungroup() %>%
    arrange(desc(Percent)) %>%
    mutate(Genus = factor(Genus, levels = c(Genus[Genus != "others"], "others")))
  
  plot <- ggplot(data_processed, aes(x = "", y = Percent, fill = Genus)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_brewer(palette = "Set3") + 
    theme_minimal() + 
    theme(
      panel.border = element_rect(fill = NA, colour = "black", size = 1),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    ) +
    coord_polar(theta = "y")

  return(plot)
}

#-------------------------------------------------------------------------------
# 3. Data Processing: Calculate Relative Abundance
#-------------------------------------------------------------------------------

my.species.100 <- data.frame(
  Percent = colSums(my.species.data[setdiff(colnames(my.species.data), 'Group')]) / sum(colSums(my.species.data[setdiff(colnames(my.species.data), 'Group')])) * 100, 
  Species = setdiff(colnames(my.species.data), 'Group')
)

adult.china.wgs.100 <- data.frame(
  Percent = colSums(adult_china_wgs.species[setdiff(colnames(adult_china_wgs.species), c('Group', 'Class'))]) / sum(colSums(adult_china_wgs.species[setdiff(colnames(adult_china_wgs.species), c('Group', 'Class'))])) * 100, 
  Species = setdiff(colnames(adult_china_wgs.species), c('Group', 'Class'))
)

adult.usa.wgs.100 <- data.frame(
  Percent = colSums(adult_usa_wgs.species[setdiff(colnames(adult_usa_wgs.species), c('Group', 'Class'))]) / sum(colSums(adult_usa_wgs.species[setdiff(colnames(adult_usa_wgs.species), c('Group', 'Class'))])) * 100, 
  Species = setdiff(colnames(adult_usa_wgs.species), c('Group', 'Class'))
)

adult.china688.wgs.100 <- data.frame(
  Percent = colSums(adult_china688_wgs.species[setdiff(colnames(adult_china688_wgs.species), c('Group', 'Class'))]) / sum(colSums(adult_china688_wgs.species[setdiff(colnames(adult_china688_wgs.species), c('Group', 'Class'))])) * 100, 
  Species = setdiff(colnames(adult_china688_wgs.species), c('Group', 'Class'))
)

child.uk.wgs.data[is.na(child.uk.wgs.data)] <- 0
child.uk.wgs.100 <- data.frame(
  Percent = colSums(child.uk.wgs.data[setdiff(colnames(child.uk.wgs.data), c('Group', 'Class'))]) / sum(colSums(child.uk.wgs.data[setdiff(colnames(child.uk.wgs.data), c('Group', 'Class'))])) * 100, 
  Species = setdiff(colnames(child.uk.wgs.data), c('Group', 'Class'))
)

#-------------------------------------------------------------------------------
# 4. Export Visualizations (Figure 1B)
#-------------------------------------------------------------------------------

# Create a 'plots' directory if it doesn't exist locally
if(!dir.exists(here("Top_Abun"))) dir.create(here("Top_Abun"))

# Save plots to PDF using the 'here' relative path
pdf(file = here("Top_Abun", "my.species.top10.pdf"))
plot_top10_stacked_bar(my.species.100, 'Percent', 'Species')
dev.off()

pdf(file = here("Top_Abun", "Adult.China.Species.Top10.pdf"))
plot_top10_stacked_bar(adult.china.wgs.100, 'Percent','Species')
dev.off()

pdf(file = here("Top_Abun", "Adult.USA.Species.Top10.pdf"))
plot_top10_stacked_bar(adult.usa.wgs.100, 'Percent', 'Species')
dev.off()

pdf(file = here("Top_Abun", "Adult.China688.Species.Top10.pdf"))
plot_top10_stacked_bar(adult_china688_wgs.species.100 , 'Percent', 'Species')
dev.off()

pdf(file = here("Top_Abun", "Child.UK.Species.Top10.pdf"))
plot_top10_stacked_bar(child.uk.wgs.100, 'Percent', 'Species')
dev.off()