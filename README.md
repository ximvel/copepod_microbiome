# Feeding strategy and dietary preference shape the microbiome of epipelagic copepods in a warm nutrient-impoverished ecosystem
Copepod microbiome analysis

## Bioinformatics
The obtained sequences were processed with R v.4.1.0 package "biostring” and DADA2 pipeline [dada2 tutorial](https://benjjneb.github.io/dada2/tutorial_1_8.html). Primer sequences were trimmed from the paired-end reads using ”biostring” and then, DADA2 was used for sequence assembly, quality filtration, chimera removal, preparation of the amplicon sequence variants (ASVs), and ASVs taxonomical assignment by Silva 16S database v. 138. See the following R script: `16s_DADA2_Biostring_pipeline_copepod_microbiome.R`  
The contaminant sequences were filtered using the “decontam” package, based on the negative controls. The presumed contaminant ASVs that were not detected by the prevalence method [decontam tutorial](https://benjjneb.github.io/decontam/vignettes/decontam_intro.html) and were present in the negative controls were also filtered. Eukaryotes, chloroplast, and mitochondria-derived reads were removed from the dataset using the “dplyr” package in R. See following R scripts: `Decontam_data_filtration.R`. The ASVs with ≤20 read counts throughout the entire dataset were removed and the remaining data was normalized by total sum scaling (TSS). this refers to the following R script: `phyloseq_data_analysis.R`. 


## Statistical analysis and data visualization
The following data was used for the anaylsis`taxa_copepod_microbiome_2020.csv`, `asvabundance_copepod_microbiome_2020.csv`, `metadata_copepod_microbiome_2020`. To see sequences of the ASVs pleasee check: `asvabundance_taxa_merged_copepod_microbiome_2020.csv`. 

Data analyses and visualization were conducted using the “phyloseq”, “vegan” and “ggplot2” packages in R, see following R script:`phyloseq_data_analysis.R`. 
Alpha diversity indices (Shannon diversity and Chao richness) were calculated for the copepod microbiomes and seawater for each season, and significant differences across samples and seasons were estimated by two-way ANOVA. Post-hoc Tukey analyses were performed to compare the differences between the host-associated and the seawater microbiota, and between seasons. The major associated microbial families were illustrated in a heatmap, and a dendrogram analysis was performed with hierarchical clustering based on Bray–Curtis dissimilarity using the “pheatmap” package to analyze the microbial composition and similarity. Based on the observed family-level clusters, the associated microbiota per season were represented by chord diagram plots using the “circlize” package.

A principal coordinate analysis (PCoA) was created based on the Bray-Curtis distances to assess host specificity. To test whether the microbial diversity is species-specific, permutational multivariate analysis of variance (PERMANOVA) with pairwise comparisons was applied. Environmental factors were checked for collinearity using the Spearman correlation. The factors that showed a correlation value of <0.75 or > −0.75 were considered non-collinear and used to estimate their effect on the copepod microbial composition and the seawater. The significance of each variable was tested using an ANOVA-like permutational test for constrained correspondence analysis (anova.cca) and the significant environmental factors (p-value <0.05) were represented in a Bray-Curtis distance-based Redundancy Analysis (db-RDA). The explanatory value (in %) of the factors was calculated with a variation partitioning analysis.

The associated microbial core taxa of each copepod species and the seawater were identified using the “microbiome” package. Taxa that were prevalent in >30% of the samples (in each copepod species, or the seawater) with relative abundances >0.005 were considered as core microbiome (see following data:core_ponticus.csv`, `core_stylifera.csv`, `core_nana.csv`). The ASVs that were detected in >70% of the samples (relative abundances >0.005) were identified per each copepod species and denominated as temporal core. To determine which specific core ASVs were significantly different between the seasons, the Kruskal-Wallis test was performed, and the generated p-values were corrected for multiple testing using the False Discovery Rate (FDR) Benjamini-Hochberg algorithm. 

The metabolic potential of the core microbiota in the seawater and copepod species was predicted from the 16S rRNA amplicon data using Tax4Fun2 [Tax4Fun tutorial](https://github.com/bwemheu/Tax4Fun2). We searched metabolic pathways that may be relevant to the metabolic interaction between the copepods and the associated microbiota (see the follwing data:`selected_ponticus_pathways.csv`, `selected_stylifera_pathways.csv`,`selected_nana_pathways.csv`, `selected_seawater_pathways.csv`) and performed a PCoA analysis (Bray-Curtis distances) to assess the specificity of the predicted functions to the specific core ASVs. The differentially abundant metabolic pathways between copepods (mean relative abundance across copepod species) and seawater were determined by the Kruskal-Wallis test and the p-values were adjusted for multiple testing. See following R script: `Tax4Fun2_metabolic_analysis`. 

## Scripts, data and methods are part of the following study:
Velasquez, X., Morov, A.R., Belkin, N., Kurt-Terbiyik, T., Rubin-Blum, M., Meron, D., Tchernov, D., & Guy-Haim, T. (2022). Feeding strategy and dietary preference shape the microbiome of epipelagic copepods in a warm nutrient-impoverished ecosystem. *Environmental DNA*, 00:1–18. https://doi.org/10.1002/edn3.357



