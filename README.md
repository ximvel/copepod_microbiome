# Feeding strategy and dietary preference shape the microbiome of epipelagic copepods in a warm nutrient-impoverished ecosystem
Copepod microbiome analysis

## Bioinformatics
The obtained sequences were processed with R v.4.1.0 package "biostring” and DADA2 pipeline [dada2 tutorial](https://benjjneb.github.io/dada2/tutorial_1_8.html). Primer sequences were trimmed from the paired-end reads using ”biostring” and then, DADA2 was used for sequence assembly, quality filtration, chimera removal, preparation of the amplicon sequence variants (ASVs), and ASVs taxonomical assignment by Silva 16S database v. 138. See the following R code `16s_DADA2_Biostring_pipeline_Copepod_Microbiome.R`  
The contaminant sequences were filtered using the “decontam” package, based on the negative controls. The presumed contaminant ASVs that were not detected by the prevalence method [decontam tutorial] (https://benjjneb.github.io/decontam/vignettes/decontam_intro.html) and were present in the negative controls were also filtered. Eukaryotes, chloroplast, and mitochondria-derived reads were removed from the dataset using the “dplyr” package in R. See following R code `Decontam_analysis.R`. The ASVs with ≤20 read counts throughout the entire dataset were removed and the remaining data was normalized by total sum scaling (TSS). this refers to the following R code: `phyloseq_data_analysis.R`. 
The following data was used for the stadistical anaylsis

## Statistical analysis and data visualization
The following data was used for the stadistical anaylsis`taxa_data_copeod_microbiome_2020.csv`, 

`taxa_data_new_data_2020.csv`

```
this is
very
long
code
```
I used the following 3 R scripts to generate the ASV tables:
* `code_a.R`
* `code_b.R`
* `code_c.R`
This is in **copepod** *copepod*.

## smaller
### more small
#### even smaller


## Citation
If you use the code in this repository, please cite the following paper:
