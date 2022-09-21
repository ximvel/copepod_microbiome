library(decontam)

#input data 
asv_df = read.csv(file="DADA2_ASV_abundance_copepod_microbiome_2020.csv", header=T, na.strings="NA",row.names = 1, stringsAsFactors=F, check.names=F, sep=",")
taxa_df = read.csv(file="DADA2_taxonomy_copepod_microbiome_2020.csv", header=T, na.strings="NA", row.names = 9, stringsAsFactors=F, check.names=F, sep=",")
metadata_df = read.csv(file="metadata_for_decontam.csv", header=T, na.strings="NA", row.names = 1,stringsAsFactors=F, check.names=F, sep=",")

#taxa_df =taxa_df[,-8]

dim(asv_df)
dim(metadata_df)
all(colnames(asv_df) == rownames(metadata_df))
all(rownames(taxa_df) == rownames(asv_df))

#create matrix
asv_mat = as.matrix(asv_df)
taxa_mat = as.matrix(taxa_df)

#Put all 3 tables into phyloseq object:
asv = otu_table(asv_mat, taxa_are_rows = TRUE)
tax = tax_table(taxa_mat)
metadata = sample_data(metadata_df)
copmicrobio = phyloseq(asv, tax, metadata)
copmicrobio


#Decontam: Identify Contaminants - Prevalence (based on the negative control):
sample_data(copmicrobio)$is.neg = sample_data(copmicrobio)$sample_negative == "negative"
contamdf.prev = isContaminant(copmicrobio, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
#More aggresive filtering:
contamdf.prev = isContaminant(copmicrobio, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev$contaminant)
#the locations (row number of the contaminants)
which(contamdf.prev$contaminant)
#This list now holds the row numbers with the contaminations
contam_row_numbers = which(contamdf.prev$contaminant)
contam_row_numbers_asv = paste("asv", as.character(contam_row_numbers), sep = "")

#Select all the non-contam ASVs:
asv_df_filt = asv_df %>% filter (!rownames(asv_df) %in% contam_row_numbers_asv)
taxa_df_filt = taxa_df %>% filter (!rownames(taxa_df) %in% contam_row_numbers_asv)

#Combine the taxonomy and the abundance tables:
combined_taxa_abundance = cbind(taxa_df_filt, asv_df_filt)

#Save the filtered data to the computer:
write.table(asv_df_filt, file="DADA2_ASV_abundance_copepod_microbiome_2020_Filter.csv", sep = ",", quote = F, row.names = T, col.names = T)
write.table(taxa_df_filt, file="DADA2_taxonomy_copepod_microbiome_2020_Filter.csv", sep = ",", quote = F, row.names = T, col.names = T)
write.table(combined_taxa_abundance, file="DADA2_taxa_abundance_copepod_microbiome_2020_Filter.csv", sep = ",", quote = F, row.names = T, col.names = T)


cont = read.csv(file="cont.csv", header=T, na.strings="NA", stringsAsFactors=F, check.names=F, sep=",")
asv_df_filt1 = read.csv(file="DADA2_ASV_abundance_copepod_microbiome_2020_Filter.csv", header=T, na.strings="NA", stringsAsFactors=F, check.names=F, sep=",")
taxa_df_filt1 = read.csv(file="DADA2_taxonomy_copepod_microbiome_2020_Filter.csv", header=T, na.strings="NA", stringsAsFactors=F, check.names=F, sep=",")

asv_df_filt2 = asv_df_filt1 %>% filter(!asv_df_filt1$asvID %in% cont$asvID)
taxa_df_filt2 = taxa_df_filt1 %>% filter (!taxa_df_filt1$asvID %in% cont$asvID)
asv_df_filt3= asv_df_filt2[,-1]

combined_taxa_abundance1 = cbind(taxa_df_filt2, asv_df_filt3)
#combined_taxa_abundance1 = combined_taxa_abundance1[,-9]

combined_taxa_abundance2 = combined_taxa_abundance1 %>% filter (Family != "Mitochondria" | is.na(combined_taxa_abundance1$Family)) # is.na funtion was added to prevent the automatic filtration of na values from all the taxonomical labes
combined_taxa_abundance3 = combined_taxa_abundance2 %>% filter (Order != "Chloroplast" | is.na(combined_taxa_abundance2$Order))
combined_taxa_abundance4 = combined_taxa_abundance3 %>% filter (Class != "Chloroplast" | is.na(combined_taxa_abundance3$Class))
combined_taxa_abundance5 = combined_taxa_abundance4 %>% filter (Kingdom != "Eukaryota")# | is.na(combined_taxa_abundance4$Kingdom))

combined_taxa_abundance6 = combined_taxa_abundance5 %>% select(-Neg1,-Neg2,-Neg3,-Neg4,-Neg5)
combined_taxa_abundance7 = combined_taxa_abundance6 %>% select(-DOO1,-DOO2,-DOO3,-DOO4,-DOO5,-DON1,-DON2,-DON3,-DON4,-DON5,-SBN1,-SBN2,-SBN3) 

taxa_data = combined_taxa_abundance7[, c(1:8)]
asv_data = combined_taxa_abundance7[,c(1,10:length(combined_taxa_abundance7))]

write.table(asv_data, file="asvabundance_copepod_microbiome_2020.csv", sep = ",", quote = F, row.names = F, col.names = T)
write.table(taxa_data, file="taxa_copepod_microbiome_2020.csv", sep = ",", quote = F, row.names = F, col.names = T)
write.table(combined_taxa_abundance7, file="asvabundance_taxa_merged_copepod_microbiome_2020.csv", sep = ",", quote = F, row.names = F, col.names = T)





