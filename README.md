# Pediatric-MASLD-Omics
Multi-omic Integration (Metagenomics + Metabolomics) and Meta-analysis of Pediatric MASLD Patients and Public Datasets

Project Overview
This repository contains the integrated bioinformatic workflow and R scripts for the Pediatric-MASLD-Omics project. The study explores the complex interplay between the gut microbiome and host metabolism in children with Metabolic Dysfunction-Associated Steatotic Liver Disease (MASLD).

By combining primary clinical samples with high-quality publicly available datasets, this project identifies core microbial signatures and metabolic pathways associated with pediatric MASLD.

1) Metagenomics: Taxonomic profiling (Phylum to Species), functional annotation (KO, KEGG, CAZyme), and diversity analysis.
2) Metabolomics: Preprocessing of serum metabolic profiles and identification of differential metabolites.
3) Cross-Study Validation: Meta-analysis incorporating external public cohorts to ensure the robustness of the identified biomarkers.
4) Multi-omic Integration (Adonis/PERMANOVA): Quantifying the variance explained ($R^2$) by clinical traits (e.g., BMI, ALT, Lipid profiles) across all omics layers.
5) Taxonomic Contribution Analysis: Identifying "Driver Species" for specific metabolic functions (e.g., TCA cycle, GAG degradation) using contribution decomposition to link microbial members to functional outputs.
6) Clinical Correlation (Spearman): Evaluating the strength and direction of associations between differential microbes/metabolites and host clinical parameters (e.g., FibroScan kPa, HOMA-IR).

Citation & Contact
If you find this code or the associated datasets helpful for your research, please cite our work.
