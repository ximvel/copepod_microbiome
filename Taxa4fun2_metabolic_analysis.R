library(Tax4Fun2)

###Build the default reference database
#buildReferenceData(path_to_working_directory = ".", use_force = FALSE, install_suggested_packages = TRUE)
###Install dependencies
#buildDependencies(path_to_reference_data = "./Tax4Fun2_ReferenceData_v2", install_suggested_packages = TRUE)


###Making functional predictions using the default reference data only
# 1. Run the reference blast
runRefBlast(path_to_otus = "./Taxa4Fun2_core_microbiome/core_ponticus_seq.fasta", path_to_reference_data = "./Tax4Fun2_ReferenceData_v2", path_to_temp_folder = "./Taxa4Fun2_core_microbiome/temp_core_ponticus", database_mode = "Ref100NR", use_force = T, num_threads = 4)
#run for all the other copepod species and seawater

# 2) Predicting functional profiles
makeFunctionalPrediction(path_to_otu_table = "./Taxa4Fun2_core_microbiome/core_ponticus.txt", path_to_reference_data = "./Tax4Fun2_ReferenceData_v2", path_to_temp_folder = "./Taxa4Fun2_core_microbiome/temp_core_water", database_mode = "Ref100NR", normalize_by_copy_number = TRUE, min_identity_to_reference = 0.97, normalize_pathways = FALSE)
#run for all the other copepod species and seawater

######KW
pathway_ponticus= read.csv(file="./Taxa4Fun2_core_microbiome/selected_ponticus_pathways.csv", header=T, na.strings="NA", stringsAsFactors=F,row.names = 1, check.names=F, sep=",")
pathway_nana= read.csv(file="./Taxa4Fun2_core_microbiome/selected_nana_pathways.csv", header=T, na.strings="NA", stringsAsFactors=F,row.names = 1, check.names=F, sep=",")
pathway_stylifera= read.csv(file="./Taxa4Fun2_core_microbiome/selected_stylifera_pathways.csv", header=T, na.strings="NA", stringsAsFactors=F,row.names = 1, check.names=F, sep=",")
pathway_water= read.csv(file="./Taxa4Fun2_core_microbiome/selected_water_pathways.csv", header=T, na.strings="NA", stringsAsFactors=F,row.names = 1, check.names=F, sep=",")

#combine all KO numbers with levels and add as rowname
rownames(pathway_ponticus) = paste(rownames(pathway_ponticus), pathway_ponticus$level1, pathway_ponticus$level2, sep="_")
rownames(pathway_nana) = paste(rownames(pathway_nana), pathway_nana$level1, pathway_nana$level2, sep="_")
rownames(pathway_stylifera) = paste(rownames(pathway_stylifera), pathway_stylifera$level1, pathway_stylifera$level2, sep="_")
rownames(pathway_water) = paste(rownames(pathway_water), pathway_water$level1, pathway_water$level2, sep="_")

###delete the KO levels from the columns
pathway_ponticus = pathway_ponticus[,-1:-3]
pathway_nana = pathway_nana[,-1:-3]
pathway_stylifera = pathway_stylifera[,-1:-3]
pathway_water = pathway_water[,-1:-3]

#replace all spaces with "_" in the rownames and add this to a new col "Ko"
pathway_ponticus$Ko = gsub(" ","_", rownames(pathway_ponticus))
pathway_nana$Ko = gsub(" ","_", rownames(pathway_nana))
pathway_stylifera$Ko = gsub(" ","_", rownames(pathway_stylifera))
pathway_water$Ko = gsub(" ","_", rownames(pathway_water))

#pathway_ponticus$Ko = rownames(pathway_ponticus)
#pathway_nana$Ko = rownames(pathway_nana)
#pathway_stylifera$Ko = rownames(pathway_stylifera)
#pathway_water$Ko = rownames(pathway_water)

#Merge the 4 tables by "Ko" column
merged1 = merge(pathway_ponticus, pathway_nana, by = "Ko", all = T)
merged2 = merge(merged1, pathway_stylifera, by = "Ko", all = T)
merged3 = merge(merged2, pathway_water, by = "Ko", all = T)
rownames(merged3) = merged3$Ko
merged3 = merged3[,-1] #removed the KO column
merged3 = merged3[,c(-23,-25,-26,-36,-39)]
dist_merged3 = vegdist(as.data.frame(t(merged3)), method="bray", binary = F, na.rm=T)

write.table(merged3, file="./Taxa4Fun2_core_microbiome/combined_patways_new.csv", sep = "\t", row.names=F, col.names=T, quote = F)

PCoA = cmdscale(dist_merged3, k=2, eig=T)
exp=round(100*(abs(PCoA$eig))/(sum(abs(PCoA$eig))),1)
df=as.data.frame(PCoA[[1]]) #eigvalues are in 1 index of the list.

metadatadf_v2 = metadatadf[c(-23,-25,-26,-36,-39),]


ggplot(df, aes(df[,1], df[,2], color=metadatadf_v2$Cop_sp)) + #, label=factor_map$Pa #, label=factor_map$"Sample_MAG-WGS_ID"
  geom_point(aes(color = metadatadf_v2$Cop_sp), size=5, stroke=1.3) +
  scale_color_manual(values=c("gray78","darkslategray3","gray31", "deepskyblue1")) + #For discrete colors for each point (categories) #For abx-project
  #scale_shape_manual(values=c(17, 16)) + #17/16 - filled triangle/circle; 2/1 - opened triangle/circle
  #geom_text(size=(3),vjust=1.2,hjust=1.2) +
  stat_ellipse(type = "t",linetype = 2) +
  labs(x=paste0(exp[1],"%"),y=paste0(exp[2],"%")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title = element_text(size=18, face="bold")) +
  theme(legend.title = element_text(size=15), legend.text = element_text(size=13)) +
  theme(legend.position= "right") +
  theme(legend.title = element_blank())+
  #theme(legend.background = element_rect(colour = 'black', fill = 'grey93', size = 0.5, linetype='solid')) +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(panel.background = element_rect(fill='white', colour='black')) #white background, black borderline
ggsave("./Taxa4Fun2_core_microbiome/PCoA_functions.pdf", width=12, height=9, dpi = 300)

merged3 = as.data.frame(t(merged3))
metadatadf_v2$Cop_vs_SW = c(rep("cop", 49), rep("sw", 12))
merged3 = merged3[,c(-9,-11,-36)]


test_res = c() #an empty vector to store the p-values from the KW test
#mean_ponticus = c() #an empty vector to store the mean abundances of each feature in specific group
#mean_nana = c()
#mean_stylifera = c()
mean_water = c()
mean_cop = c()

for (i in 1:length(merged3)){
  print (colnames(merged3)[i])
  test_res[i] = kruskal.test(merged3[,i], metadatadf_v2$Cop_vs_SW)$p.value
  #mean_ponticus[i] = mean(merged3[metadatadf_v2$Cop_sp=="C. ponticus", i])
  #mean_nana[i] = mean(merged3[metadatadf_v2$Cop_sp=="O. nana", i])
  #mean_stylifera[i] = mean(merged3[metadatadf_v2$Cop_sp=="T. stylifera", i])
  mean_water[i] = mean(merged3[metadatadf_v2$Cop_vs_SW=="sw", i])
  mean_cop[i] = mean(merged3[metadatadf_v2$Cop_vs_SW=="cop", i])
}

#Combine the result vector of p-values with features names and mean values to make a df: (add as many mean_group items as needed)
test_res_pathway = cbind.data.frame(colnames(merged3),mean_cop , mean_water,test_res, stringsAsFactors=FALSE)
colnames(test_res_pathway) = c("Ko","mean_cop", "mean_water", "p_value")

#Run p.adjust() on the raw p-valuess to correct for multiple-testing, using FDR (BH, Benjamini-Hochberg): 
test_res_pathway$p_fdr_BH = p.adjust(test_res_pathway$p_value, method = "BH")
#Keep all the significant features after p-value FDR<0.05
test_res_pathway = test_res_pathway[test_res_pathway$p_fdr_BH < 0.05, ]

#Save the results
write.table(test_res_pathway, file="./Taxa4Fun2_core_microbiome/KW_test_pathways_meancop_new.csv", sep = "\t", row.names=F, col.names=T, quote = F)

test_res_pathway$foldchange = test_res_pathway$mean_cop/test_res_pathway$mean_water
test_res_pathway$foldchange_log =log10(test_res_pathway$foldchange)
test_res_pathway = test_res_pathway[order(test_res_pathway[,7]), ] #sort by log-FC
test_res_pathway$group =  c(rep("sw", 30), rep("cop", 21))
test_res_pathway = test_res_pathway[,-2:-6]
write.table(test_res_pathway, file="./Taxa4Fun2_core_microbiome/KW_pathways_organize.csv", sep = "\t", row.names=F, col.names=T, quote = F)

#optional repalce spaces or other characters
test_res_pathway$Ko = gsub("_", " ", test_res_pathway$Ko)

pathways_cop_sw = read.csv(file="./Taxa4Fun2_core_microbiome/barplot_pathways_1.csv", header=T, na.strings="NA", stringsAsFactors=F, check.names=F, sep=",")
#pathways_cop_sw$pathway = factor(pathways_cop_sw$pathway, levels = c("Xenobiotics biodegradation and metabolism","Metabolism of other amino acids","Metabolism of cofactors and vitamins","Energy metabolism","Carbohydrate metabolism","Amino acid metabolism"))
pathways_cop_sw$group = factor(pathways_cop_sw$group)
pathways_cop_sw$Ko = factor(pathways_cop_sw$Ko, levels = pathways_cop_sw$Ko[order(pathways_cop_sw$pathway)])

ggplot(pathways_cop_sw, aes(x = foldchange_log, y =Ko, fill =group)) + 
  geom_bar(stat = "identity") + #coord_flip() + 
  scale_fill_manual(values=c("orange","deepskyblue1"))+#= colorRampPalette(rev(brewer.pal(, "Blues")))(22)) +
  theme(axis.title.y = element_text(size=14, face="bold")) +
  theme(axis.text.y  = element_text(size=12, face="bold")) +
  theme(axis.title.x  = element_text(size=14, face="bold"))+
  theme(axis.text.x  = element_text(size=12, face="bold"))+
  theme(legend.title = element_blank(),#(size = 18, face= "bold"),
        legend.text = element_text(size = 14, face= "bold"))+
  labs(x= "fold change [log10]", y="KO pathways")+
  theme(panel.grid.major.y = element_line(color = "black"))+
  theme(panel.background = element_rect(fill='white', colour='black'))
ggsave("./Taxa4Fun2_core_microbiome/patways_copvswater.pdf", width=20, height=10, dpi = 300)

