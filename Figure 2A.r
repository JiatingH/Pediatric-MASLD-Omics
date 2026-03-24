################################################################################
# Project: NAFLD/MASLD Omics Analysis
# Task: Co-occurrence Network Analysis (SparCC)
# Logic: 1. Filter by Significance -> 2. Network Construction -> 3. Topological Ranking
################################################################################

# --- 1. Load Required Libraries ---
library(readxl)
library(gdata)
library(reshape2)
library(Hmisc)
library(tableone)
library(ggpubr)
library(stringr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(viridis)
library(matrixStats)
library(vegan)
library(corrplot)
library(ade4)
library(readr)
library(igraph) # Essential for network/graph theory

# --- 2. Data Import (SparCC Results) ---
# Note: SparCC is computationally intensive and usually run on a Linux server.
# Here we load the pre-calculated correlation (Cor) and significance (P) tables.

# Load Correlation Matrices
hc.sparcc.cor = read.table('[Path]/All.hc.sparcc.cor.table.csv', sep = ',', header = 1, row.names = 1, check.names = FALSE)
ob.sparcc.cor = read.table('[Path]/All.ob.sparcc.cor.table.csv', sep = ',', header = 1, row.names = 1, check.names = FALSE)
nafld.sparcc.cor = read.table('[Path]/All.nafld.sparcc.cor.table.csv', sep = ',', header = 1, row.names = 1, check.names = FALSE)

# Load P-value Matrices
hc.sparcc.p = read.table('[Path]/All.hc.sparcc.p.table.csv', sep = ',', header = 1, row.names = 1, check.names = FALSE)
ob.sparcc.p = read.table('[Path]/All.ob.sparcc.p.table.csv', sep = ',', header = 1, row.names = 1, check.names = FALSE)
nafld.sparcc.p = read.table('[Path]/All.nafld.sparcc.p.table.csv', sep = ',', header = 1, row.names = 1, check.names = FALSE)

# --- 3. Significance Filtering ---

# Handle missing values (NA) by setting p-value to 1 (non-significant)
hc.sparcc.p[is.na(hc.sparcc.p)] = 1
ob.sparcc.p[is.na(ob.sparcc.p)] = 1
nafld.sparcc.p[is.na(nafld.sparcc.p)] = 1

# Stringent Filtering: Set correlation to 0 if p >= 0.05
# This ensures only statistically robust edges remain in the network
hc.sparcc.cor[hc.sparcc.p >= 0.05] = 0
ob.sparcc.cor[ob.sparcc.p >= 0.05] = 0
nafld.sparcc.cor[nafld.sparcc.p >= 0.05] = 0

# --- 4. Network Object Construction ---

# Create adjacency graphs from matrices (undirected, weighted by absolute correlation)
hc.g = graph.adjacency(as.matrix(abs(hc.sparcc.cor)), weight = TRUE, mode = "undirected")
ob.g = graph.adjacency(as.matrix(abs(ob.sparcc.cor)), weight = TRUE, mode = "undirected")
nafld.g = graph.adjacency(as.matrix(abs(nafld.sparcc.cor)), weight = TRUE, mode = "undirected")

# Clean Networks: Remove isolated nodes (species with Degree = 0)
hc.g.deleted = delete_vertices(hc.g, names(degree(hc.g)[degree(hc.g) == 0]))
ob.g.deleted = delete_vertices(ob.g, names(degree(ob.g)[degree(ob.g) == 0]))
nafld.g.deleted = delete_vertices(nafld.g, names(degree(nafld.g)[degree(nafld.g) == 0]))

# Export to GraphML format for advanced visualization in Gephi
write_graph(hc.g.deleted, '[Path]/hc.network.graphml', format = 'graphml')
write_graph(ob.g.deleted, '[Path]/ob.network.graphml', format = 'graphml')
write_graph(nafld.g.deleted, '[Path]/nafld.network.graphml', format = 'graphml')

