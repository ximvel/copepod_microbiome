library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(gplots)
library(pheatmap)
library(reshape2)
library(RVAideMemoire)
library(GGally)
library(Hmisc)
library(corrplot)

#Phyloseq objects need to have row.names(asvID)
asvdf = read.csv(file="asvabundance_copepod_microbiome_2020.csv", header=T, na.strings="NA", stringsAsFactors=F,row.names = 1, check.names=F, sep=",")
taxadf = read.csv(file="taxa_copepod_microbiome_2020.csv", header=T, na.strings="NA", stringsAsFactors=F, row.names = 1, check.names=F, sep=",")
metadatadf = read.csv(file="metadata_copepod_microbiome_2020.csv", header=T, na.strings="NA", stringsAsFactors=F, row.names = 1, check.names=F, sep=",")

asvdf = asvdf %>% select(-ONO4)
metadatadf = metadatadf %>% filter (rownames(metadatadf) != "ONO4")

#Change the order of specific metadata columns
metadatadf$Cop_sp = factor(metadatadf$Cop_sp, levels = c("C. ponticus","O. nana","T. stylifera","Sea water"))
metadatadf$Season = factor(metadatadf$Season, levels = c("Winter","Spring","Summer","Autumn"))
#metadatadf$Month = factor(metadatadf$Month, levels = c("February","April","July","October"))


dim(asvdf)
dim(metadatadf)
dim (taxadf)

all(colnames(asvdf) == rownames(metadatadf))
all(rownames(taxadf) == rownames(asvdf))

#create matrix
asv_mat = as.matrix(asvdf)
taxa_mat = as.matrix(taxadf)

#Put all 3 tables into phyloseq object:
asv = otu_table(asv_mat, taxa_are_rows = TRUE)
tax = tax_table(taxa_mat)
metadata = sample_data(metadatadf)
copmicro = phyloseq(asv, tax, metadata)
copmicro

# Remove singletons (taxa with less than 20-counts are removed)!
copmicro_filt = prune_taxa(taxa_sums(copmicro) > 20, copmicro) 
copmicro_filt

# Make a data frame with a column for the read counts of each sample
sample_sum_df = data.frame(sum = sample_sums(copmicro_filt))
sample_sum_df$asvnumb = colSums(asv_table_filt != 0)
write.table(sample_sum_df, file="read_counts_final.csv", sep = ",", quote = F, row.names = TRUE, col.names = T)

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

# mean, median, max and min of sample read counts
smin = min(sample_sums(copmicro_filt))
smean = mean(sample_sums(copmicro_filt))
smax = max(sample_sums(copmicro_filt))
smedian = median(sample_sums(copmicro_filt))

#check the data before analysis
asv_table_filt = data.frame(otu_table(copmicro_filt))
taxa_table_filt = data.frame(tax_table(copmicro_filt))
combined_filt = cbind(taxa_table_filt, asv_table_filt)
write.table(combined_filt, file="combined_copmicro_filt_avs_taxa_check_2020.csv", sep = ",", quote = F, row.names = T, col.names = T)


#Normalize each sample read counts into relative abundance units:
copmicro_ra  = transform_sample_counts(copmicro_filt, function(x) x / sum(x)) 
sample_sums(copmicro_ra) #check that each sample is sum to 1 (100%)
asv_table_ra = data.frame(otu_table(copmicro_ra))
taxa_table_ra = data.frame(tax_table(copmicro_ra))
colSums(asv_table_ra) #check that each sample is sum to 1 (100%)
copmicro_ra 


#prepare data for heatmap
copmicro_ra_agg_family = tax_glom(copmicro_ra, taxrank = "Family")
copmicro_ra_agg_family_filt = filter_taxa(copmicro_ra_agg_family, function(x) mean(x) > 0.005, prune=TRUE)
copmicro_ra_agg_family_filt

#Calculate distance and save as a matrix
asv_family_df = data.frame (otu_table(copmicro_ra_agg_family_filt))
taxa_family_df = data.frame (tax_table(copmicro_ra_agg_family_filt))
rownames(asv_family_df)= taxa_family_df$Family

colormap1 = c("mediumblue","orange","red","forestgreen")
names(colormap1) = levels(metadatadf$Season)
colormap2 = c("white","white","white","white") #"gray0","gray31","gray","gray94", khaki1","hotpink2","mediumorchid1","steelblue2
names(colormap2) = levels(metadatadf$Cop_sp)
mycolors = list(Season = colormap1, Cop_sp = colormap2)

#metadata_season = data.frame(metadatadf$Season)
#rownames(metadata_season) = rownames(metadatadf)
#colnames(metadata_season)[1] = "Season"
metadata_season_sp = data.frame(metadatadf$Season, metadatadf$Cop_sp)
rownames(metadata_season_sp) = rownames(metadatadf)
colnames(metadata_season_sp) = c("Season","Cop_sp")

distmat = vegdist(asv_family_df, method="bray", binary = F) #distance mat for rows (OTUs)
clust_row = hclust(distmat, method = "ward.D2") #clustering OTUs based on composition
distmat = vegdist(t(asv_family_df), method="bray", binary = F) #distance mat for columns (samples)
clust_col = hclust(distmat, method = "ward.D2") #clustering samples based on composition


#### heatmap plot ####
ph = pheatmap(log10(asv_family_df+0.0001), scale = "none", cluster_rows=clust_row, cluster_cols=clust_col, show_colnames=T, show_rownames = T, #clustering_method = "average",
              annotation_col=metadata_season_sp, annotation_row=NA, annotation_names_col=T, annotation_names_row=T, legend = T,
              annotation_colors = mycolors, fontsize = 12, cutree_rows = NA, annotation_legend = T, display_numbers = F,
              treeheight_col=100, treeheight_row=75, #clustering_distance_cols = "euclidean", #clustering_distance_rows = "euclidean",
              color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100), fontsize_row =16,fontsize_col=16, cellwidth = 20, cellheight = 20)
ggsave("Heatmap_microbiome_family_2020_final.pdf", width=25, height=15, plot=ph, dpi = 300)


####subset samples (remove the sea water samples for the PCoA and Simper)####
copmicro_ra_cop = subset_samples(copmicro_ra, Cop_sp != "Sea water")
sample_data(copmicro_ra_cop)$Cop_sp
rownames(sample_data(copmicro_ra_cop))
copmicro_ra_cop

#remove the asv with 0 and were only from the sea water 
asv_cop_df = data.frame(otu_table(copmicro_ra_cop))
asv_cop_df_nop = asv_cop_df[rowSums(asv_cop_df)!=0, ] # remove rows that contain only zeros 
meta_cop_df = data.frame(sample_data(copmicro_ra_cop))
write.table(asv_cop_df_nop, file="onlycop.csv", sep = ",", quote = F, row.names = T, col.names = T)

meta_cop_df$Season = factor(meta_cop_df$Season, levels = c("Winter","Spring","Summer","Autumn"))
meta_cop_df$Cop_sp = factor(meta_cop_df$Cop_sp, levels = c("C. ponticus","O. nana","T. stylifera","Sea water"))

#create matrix
asv_mat1 = as.matrix(asv_cop_df_nop)
taxa_mat1 = as.matrix(taxadf)


#Put all 3 tables into phyloseq object:
asv_cop = otu_table(asv_mat1, taxa_are_rows = TRUE)
tax_cop = tax_table(taxa_mat1)
metadata_cop = sample_data(meta_cop_df)
copmicro_ra_cop = phyloseq(asv_cop, tax_cop, metadata_cop)
copmicro_ra_cop


#### PCoA (unconstrained ordination analysis) #### 
#only copepods sp. 
copmicro_pcoa = ordinate(physeq = copmicro_ra_cop, method = "PCoA", distance = "bray")

# Plot 
#only Copepods sp.

plot_ordination(physeq =copmicro_ra_cop , ordination = copmicro_pcoa, color = "Season")+
  scale_color_manual(values = c("mediumblue","orange","red","forestgreen"))+
  geom_point(aes(color = Season), alpha = 0.5, size = 8)+
  #stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t",linetype = 2) +
  theme(axis.title.x = element_text(size=12, face="bold"))+
  theme(axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 14,face= "plain"))+
  theme(panel.background = element_rect(fill='white', colour='black')) +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())
ggsave("PCoA_betweenCop_season.pdf", width=8, height=5, dpi = 300)


plot_ordination(physeq =copmicro_ra_cop , ordination = copmicro_pcoa, color= "Cop_sp", )+ #shape = "Cop_sp"
  scale_color_manual(values = c("gray66","darkslategray3","gray31"))+
  #scale_shape_manual(values = c(15,17,18))+
  geom_point(aes(color = Cop_sp), size = 8, alpha = 0.7) + #alpha = 0.5 , shape = Cop_sp
  #stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t",linetype = 2) +
  theme(axis.title.x = element_text(size=12, face="bold"))+
  theme(axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 12,face= "italic"))+
  theme(panel.background = element_rect(fill='white', colour='black')) +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())
ggsave("PCoA_betweenCop_sp.pdf", width=8, height=5, dpi = 300)

plot_ordination(physeq =copmicro_ra_cop , ordination = copmicro_pcoa, color = "Season", shape = "Cop_sp")+
  scale_color_manual(values = c("mediumblue","orange","red","forestgreen"))+ #"dodgerblue3","gold","firebrick2", "chartreuse3"))
  geom_point(aes(color = Season, shape = Cop_sp), alpha = 0.5, size = 5)+
  #stat_ellipse(type = "norm", linetype = 2) +
  #stat_ellipse(type = "t") +
  theme(axis.title.x = element_text(size=12, face="bold"))+
  theme(axis.title.y = element_text(size=12, face="bold"))+
  scale_shape_manual(values = c(15,17,18))+
  theme(legend.title = element_text(size = 10, face= "bold"),
        legend.text = element_text(size = 10,face= "italic" ))+
  theme(panel.background = element_rect(fill='white', colour='black')) +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())
ggsave("PCoA_Season_Cop_sp_2020.png", width=8, height=5, dpi = 300)

####Permanova test_cop_sp#### 
# Calculate bray curtis distance matrix (only opepods sp.)
copmicro_bray = phyloseq::distance(copmicro_ra_cop, method = "bray")
# make a data frame from the sample_data (only copepod sp.)
metadata_filt = data.frame(sample_data(copmicro_ra_cop))

# Adonis test
#only Copepod sp.
adonis(copmicro_bray~ Season, data = metadata_filt)
adonis(copmicro_bray~ Cop_sp, data = metadata_filt)
adonis(copmicro_bray~ Season*Cop_sp, data = metadata_filt)
#adonis(copmicro_bray~ Season*Cop_sp, data = metadata_filt, permutations = 999)

#pairwise 
#res = pairwise.adonis(copmicro_bray, metadata_filt$Season, F = TRUE, R2=TRUE)
res = pairwise.perm.manova(copmicro_bray, metadata_filt$Season, F = TRUE, R2=TRUE)
res$R2
res$p.value
res$F.value
res

res1=pairwise.perm.manova(copmicro_bray, metadata_filt$Cop_sp, F = TRUE, R2=TRUE)
res1$R2
res1$p.value
res1$F.value


#### Alpha Diversity indexes ####
alpha_diversity = estimate_richness(copmicro_filt)
write.table(alpha_diversity, file="index_diversity_2020.csv", sep = ",", quote = F, row.names = T, col.names = T)

plot_richness(copmicro_filt, x="Cop_sp", measures=c("Shannon")) +
  facet_grid(.~ Season, scales="free")+
  geom_point(size = 0.1, stroke = 0.1)+
  geom_boxplot(aes(fill=Cop_sp), fatten = 1, lwd=0.5, alpha =0.5)+ #lwd=0.4, alpha=0.5, fatten = 1)+
  #scale_y_continuous(breaks=c(0,1,2,3,4,5.6), limits=c(0,6))+
  #geom_jitter(position=position_jitter(0.2), aes(color = Cop_sp), stroke=0.4)+
  scale_fill_manual(values=c("gray78","darkslategray3","gray31", "deepskyblue1")) + #this is for the boxplots "orange","red","green3","blue"
  theme(strip.background = element_rect(color="black", fill="black", size=1.5, linetype="solid"))+ # Change the apperance of the rectangle around facet label
  theme(strip.text.x = element_text(size=14, color="white", face="bold")) + # Change facet text font. Possible values for the font style: 'plain', 'italic', 'bold', 'bold.italic'.
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x  = element_blank()) + #element_text(size=10, face="bold.italic", color = "black"))+
  theme(axis.ticks.x = element_blank())+
  theme(axis.title.y = element_text(size=12, face="bold")) +
  theme(axis.text.y  = element_text(size=10, face="bold", color = "black")) +
  theme(legend.position='right')+
  theme(legend.text = element_text(size=12, face = "italic"),legend.title = element_blank())+
  labs(y = "Alpha diveristy (Shannon)") +
  theme(panel.background = element_rect(fill='white', colour='black'))  
ggsave("Shannon_diversity_2020.pdf", width=10, height=5, dpi = 500)

plot_richness(copmicro_filt, x="Cop_sp", measures=c("Chao1")) + 
  facet_grid(.~ Season, scales="free")+
  geom_point(size = 0.1, stroke = 0.1)+
  geom_boxplot(aes(fill=Cop_sp), fatten = 1, lwd=0.5, alpha =0.5)+
  #scale_y_continuous(breaks=c(0,1,2,3,4,5.6), limits=c(0,6))+
  #geom_jitter(position=position_jitter(0.2), aes(color = Cop_sp), stroke=0.4)+
  scale_fill_manual(values=c("gray78","darkslategray3","gray31","deepskyblue1")) + #this is for the boxplots "orange","red","green3","blue"
  theme(strip.background = element_rect(color="black", fill="black", size=1.5, linetype="solid"))+ # Change the apperance of the rectangle around facet label
  theme(strip.text.x = element_text(size=14, color="white", face="bold")) + # Change facet text font. Possible values for the font style: 'plain', 'italic', 'bold', 'bold.italic'.
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x  = element_blank()) + #element_text(size=10, face="bold.italic", color = "black"))+
  theme(axis.ticks.x = element_blank())+
  theme(axis.title.y = element_text(size=12, face="bold")) +
  theme(axis.text.y  = element_text(size=10, face="bold", color = "black")) +
  theme(legend.position="right")+
  theme(legend.text = element_text(size=12, face = "italic"),legend.title = element_blank())+
  labs(y = "Richness (Chao1)") +
  theme(panel.background = element_rect(fill='white', colour='black'))  
ggsave("Richness_2020.pdf", width=10, height=5, dpi = 500)

#2 ways-Anova 
index = read.csv(file="index_diversity_2020.csv", header=T, na.strings="NA", stringsAsFactors=F,row.names = 1, check.names=F, sep=",")

res.aov = aov(Shannon ~ Season + Cop_sp, data = index)
summary(res.aov)

res.aov1 = aov(Chao1 ~ Season + Cop_sp, data = index)
summary(res.aov1)

#TukeyHSD (res.aov)

TukeyHSD(res.aov, which = "Cop_sp")
TukeyHSD(res.aov, which = "Season")

TukeyHSD(res.aov1, which = "Cop_sp")
TukeyHSD(res.aov1, which = "Season")

#normality test 
plot(res.aov, 2)
plot(res.aov1,2)
leveneTest(Shannon ~ Season*Cop_sp, data = index)
leveneTest(Chao1 ~ Season*Cop_sp, data = index)

#Extract the residuals
aov_residuals = residuals(object = res.aov)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals)


####Check collinearity of environmental metadata#### 
meta_cop_env = data.frame(sample_data(copmicro_ra_cop))
meta_cop_env1 = meta_cop_df [,-1:-3] 


##correlation between enviromental factors##
rcorr_env_method2 = rcorr(as.matrix(meta_cop_env1), type="spearman")
rcorr_env_method2 = rcorr_env_method2$r

#correlation matrix
pdf("collinearity_env.pdf", width=8, height=5)
#png("collinearity_env.png", width=8, height=5, units = "in", res=500)
corrplot(rcorr_env_method2, type = "lower", tl.col = "black", tl.srt = 45, tl.pos='lt', diag = F,
         cl.pos = "b", tl.cex = 0.8, cl.cex = 0.8, addCoef.col = "black", number.digits = 2, number.cex = 0.6, method = "color", COL1("Blues", 10))
dev.off()


####Constrained Ordinations####
#dbRDA ordinate

dbrda_ord = ordinate(physeq = copmicro_ra_cop, method = "CAP", distance = copmicro_bray,
                     formula = ~ Temperature)

anova.cca(dbrda_ord, by = "terms")

dbrda_ord1= ordinate(physeq = copmicro_ra_cop, method = "CAP", distance = copmicro_bray,
                     formula = ~ PNEA)

anova.cca(dbrda_ord1, by = "terms")

dbrda_ord2 = ordinate(physeq = copmicro_ra_cop, method = "CAP", distance = copmicro_bray,
                      formula = ~ BP)

anova.cca(dbrda_ord2, by = "terms")

dbrda_ord3 = ordinate(physeq = copmicro_ra_cop, method = "CAP", distance = copmicro_bray,
                      formula = ~Temperature + PNEA + BP)

summary (dbrda_ord3)
dbrda_ord3


## dbRDA plot_copepods ##
rda_plot = plot_ordination(physeq=copmicro_ra_cop, ordination= dbrda_ord3, color="Season", shape="Cop_sp") + 
  scale_color_manual(values = c("mediumblue","orange","red","forestgreen"))+
  geom_point(aes(color = Season), alpha = 0.5, size = 5)+
  scale_shape_manual(values = c(15,17,18,19))+
  theme(panel.background = element_rect(fill='white', colour='black')) +
  theme(axis.title.x = element_text(size=12, face="bold"))+
  theme(axis.title.y = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 10,face= "italic" ))+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())

# Now add the environmental variables as arrows into a matrix
arrowmat = vegan::scores(dbrda_ord3, display = "bp")

#transform matrix into a dataframe, and add labels
arrowdf = data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map = aes(xend = CAP1, yend = CAP2, x = 0, y = 0, shape = NULL, color = NULL, label = labels)
label_map = aes(x = 1.3 * CAP1, y = 1.3 * CAP2, shape = NULL, color = NULL, label = labels)
arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
rda_plot + geom_segment(mapping = arrow_map, size = .8, data = arrowdf, color = "black", arrow = arrowhead) + 
  geom_text(mapping = label_map, size = 3.5, vjust=-0.5,hjust=1.2, data = arrowdf, show.legend = F)
ggsave("dbRDA_Temperature + PNEA + BP_2020.pdf", width=8, height=5, dpi = 300)


####Partition the Variation with explanatory environmental factors####
asv_cop = data.frame(otu_table(copmicro_ra_cop))

#transpose the ASV table:
asv_cop_df_t = as.data.frame(t(asv_cop))
meta_cop_env 

varpart_res = varpart(asv_cop_df_t, ~ Temperature, ~PNEA, ~BP, ~Cop_sp, data = meta_cop_env)
varpart_res

pdf("VP_for_dbRDA_enviromental_factors_cop_final_1.pdf", width=8, height=5)
#png("VP_for_dbRDA_enviromental_factors_1.png", width=8, height=5, units = "in", res=500)
plot (varpart_res, digits = 2, Xnames = c("Temp.", "PNEA", "BP", "Cop sp."), bg = c("tomato", "green","paleturquoise1","blue"), cex = 0.8, id.size = 0.8, cutoff = 1)
dev.off()

##dbRDA only sea water##
copmicro_ra_water = subset_samples(copmicro_ra, Cop_sp == "Sea water")
sample_data(copmicro_ra_water)$Cop_sp
rownames(sample_data(copmicro_ra_water))
copmicro_ra_water

#remove the asv with 0 and were only from the sea water 
copmicro_ra_water_df = data.frame(otu_table(copmicro_ra_water))
copmicro_ra_water_nop = copmicro_ra_water_df[rowSums(copmicro_ra_water_df)!=0, ] # remove rows that contain only zeros 
meta_copmicro_ra_water= data.frame(sample_data(copmicro_ra_water))
#write.table(copmicro_ra_ponticus_nop, file="onlyponticus.csv", sep = ",", quote = F, row.names = T, col.names = T)

#create matrix
asv_ra_water = as.matrix(copmicro_ra_water_nop)
taxa_ra_water = as.matrix(taxadf)


#Put all 3 tables into phyloseq object:
asv_ra_water_new = otu_table(asv_ra_water, taxa_are_rows = TRUE)
tax_ra_water = tax_table(taxa_ra_water)
meta_ra_water = sample_data(meta_copmicro_ra_water)
copmicro_water = phyloseq(asv_ra_water_new,tax_ra_water,meta_ra_water)
copmicro_water

# Calculate bray curtis distance matrix (only opepods sp.)
copmicro_water_bray = phyloseq::distance(copmicro_water, method = "bray")

#dbRDA ordinate 

ordi_water = ordinate(physeq = copmicro_water, method = "CAP", distance = copmicro_water_bray,
                      formula = ~ Temperature)

anova.cca(ordi_water, by = "terms")

ordi_water1= ordinate(physeq = copmicro_water, method = "CAP", distance = copmicro_water_bray,
                      formula = ~ PNEA)

anova.cca(ordi_water1, by = "terms")

ordi_water2 = ordinate(physeq = copmicro_water, method = "CAP", distance = copmicro_water_bray,
                       formula = ~ BP)

anova.cca(ordi_water2, by = "terms")


ordi_water3 = ordinate(physeq = copmicro_water , method = "CAP", distance = copmicro_water_bray, 
                       formula = ~Temperature + PNEA)

summary (ordi_water3)
ordi_water3
anova.cca(ordi_water3, by = "terms")

## dbRDA plot ##
rda_plot_water = plot_ordination(physeq = copmicro_water, ordination = ordi_water3, color="Season", shape="Cop_sp") + 
  scale_color_manual(values = c("mediumblue","orange","red","forestgreen"))+
  geom_point(aes(color = Season), alpha = 0.5, size = 5)+
  scale_shape_manual(values = c(19))+
  theme(panel.background = element_rect(fill='white', colour='black')) +
  theme(axis.title.x = element_text(size=12, face="bold"))+
  theme(axis.title.y = element_text(size=12, face="bold")) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 10,face= "italic" ))+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())

# Now add the environmental variables as arrows into a matrix
arrowmat = vegan::scores(ordi_water3, display = "bp")

#transform matrix into a dataframe, and add labels
arrowdf = data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map = aes(xend = CAP1, yend = CAP2, x = 0, y = 0, shape = NULL, color = NULL, label = labels)
label_map = aes(x = 1.3 * CAP1, y = 1.3 * CAP2, shape = NULL, color = NULL, label = labels)
arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
rda_plot_water + geom_segment(mapping = arrow_map, size = .8, data = arrowdf, color = "black", arrow = arrowhead) + 
  geom_text(mapping = label_map, size = 3.5, vjust=-1.2,hjust=0.5, data = arrowdf, show.legend = F)
ggsave("dbRDA_Temperature + PNEA_water.pdf", width=8, height=5, dpi = 300)

#Partition the Variation with explanatory environmental factors#
asv_water_pv = data.frame(otu_table(copmicro_water))

#transpose the ASV table:
asv_water_pv_t = as.data.frame(t(asv_water_pv))
meta_copmicro_ra_water

varpart_water = varpart(asv_water_pv_t, ~ Temperature, ~PNEA, data = meta_copmicro_ra_water)
varpart_water

pdf("VP_for_dbRDA_enviromental_factors_water_final_1.pdf", width=8, height=5)
#png("VP_for_dbRDA_clust_fam_Temperature + PNEA_water1_1.png", width=8, height=5, units = "in", res=500)
plot(varpart_water, digits = 2, Xnames = c("Temp.", "PNEA"), bg = c("tomato", "green"), cex = 0.8, id.size = 0.8,cutoff = 1)
dev.off()


#### ChordDiagram ####
library(RColorBrewer)
library(circlize)

display.brewer.all() #display all available color palletes

#this part need to be run again when different cluster is done
asv_family_df = data.frame (otu_table(copmicro_ra_agg_family_filt))
rownames(asv_family_df)= taxa_family_df$Family


#February
cp = rowMeans(asv_family_df[,1:5]) # average per each asv_family across the group
on = rowMeans(asv_family_df[,21:25])
sw = rowMeans(asv_family_df[,55:57])
#April
cp = rowMeans(asv_family_df[,6:10])
on = rowMeans(asv_family_df[,26:30])
ts = rowMeans(asv_family_df[,40:44])
sw = rowMeans(asv_family_df[,58:60])
#July
cp = rowMeans(asv_family_df[,11:15])
on = rowMeans(asv_family_df[,31:35])
ts = rowMeans(asv_family_df[,45:49])
sw = rowMeans(asv_family_df[,61:63])
#October
cp = rowMeans(asv_family_df[,16:20])
on = rowMeans(asv_family_df[,36:39])
ts = rowMeans(asv_family_df[,50:54])
sw = rowMeans(asv_family_df[,64:66])

family_list = c("Flavobacteriaceae","Saprospiraceae","Rhodobacteraceae","Cyanobiaceae","Halieaceae","Vibrionaceae","Saccharospirillaceae","Cyclobacteriaceae","Pseudoalteromonadaceae","Streptococcaceae","Pirellulaceae") # families of the cluster II 
family_list = c("Beijerinckiaceae","Hymenobacteraceae","Planococcaceae","Pseudomonadaceae","Halomicrobiaceae","Bacillaceae","Moraxellaceae","Micrococcaceae","Sphingobacteriaceae","Chitinophagaceae","Oxalobacteraceae","Comamonadaceae","Rubritaleaceae","Corynebacteriaceae","Sphingomonadaceae")# families cluster III

asv_family_df = asv_family_df %>% filter (rownames(asv_family_df) %in% family_list) # take only the families for specific cluster

# create variable with all the average x each group as a matrix
fdf = as.matrix(cbind(cp, on, sw)) 
fdf = as.matrix(cbind(cp, on, ts, sw))
fdf

#circos.clear()

# for defining the color of the ribons of the chord diagram plots
grid.col = setNames(colorRampPalette(brewer.pal(n=8, name = "Set3"))(nrow(fdf)), rownames(fdf)) ##automatically defined colors for the the ASVs_families cluster II
grid.col = setNames(colorRampPalette(brewer.pal(n=8, name = "Paired"))(nrow(fdf)), rownames(fdf))##automatically defined colors for the the ASVs_families cluster III

#manually adding colors for the sample groups (columns)
grid.col[["cp"]] = "gray78"
grid.col[["on"]] = "darkslategray3"
grid.col[["ts"]] = "gray31"
grid.col[["sw"]] = "deepskyblue1"  

chordDiagram(fdf, grid.col = grid.col)

pdf("ChordDiagram_October_ClusterIII.pdf", width=20, height=10)
#png("ChordDiagram_October_ClusterII.png", width=20, height=10, units = "in", res=350)
# now, the image with rotated labels
chordDiagram(fdf, annotationTrack = "grid", preAllocateTracks = 1, grid.col = grid.col, transparency = 0.3)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + 3, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.8)
  circos.axis(h = "top", labels.cex = 0.8, major.tick.length =mm_y(1), sector.index = sector.name, track.index = 2)
}, bg.border = NA)
#legend("left",pch=20,legend=rownames(fdf),
#col=grid.col[rownames(fdf)],bty="n",
#cex=2,pt.cex=5, xjust=-20) # Set legend
dev.off()

####core microbiome####
library(microbiome)
library(microbiomeutilities)
library(eulerr)
library ("ggpolypath")

cop_sp_names = unique(as.character(meta(copmicro_ra)$Cop_sp))
list_core = c()

for (n in cop_sp_names){ # for each variable n in cop_sp_names
  cop_sub = subset_samples(copmicro_ra, Cop_sp == n) # Choose sample from Cop_sp by n
  cop_core = core_members(cop_sub, detection = 0.005, prevalence = 30/100) #detection (minimal abundance detection limite), prevalence in which % of the samples appear 
  print(paste0("No. of core taxa in ", n, " : ", length(cop_core))) # print core taxa identified in each Cop_sp
  list_core[[n]] = cop_core # add to a list core taxa for each group.
}

# Specify colors and plot venn, supplying colors in the order they appear in cop_sp_names
mycols = c("C. ponticus"="gray78", "O. nana"="darkslategray3", "T. stylifera"="gray31", "Sea water"= "deepskyblue1") 
pdf("Venn.dagram_core.pdf", width=10, height=10)
#png("Venn.dagram_core_cop.png", width=10, height=10, units = "in", res=350)
eulerr_options(edges = list(col = "black",lwd = 2),
               quantities = list(fontsize = 22))
plot(venn(list_core), fills = mycols, labels = F)#,legend = list(side = "right"))
dev.off()

#Separate objects for core members in each species
cop_bysp = subset_samples(copmicro_ra_cop, Cop_sp == "T. stylifera")
cop_bysp_core_stylifera = core(cop_bysp, detection = 0.005, prevalence = 30/100) #detection (minimal abundance detection limite), prevalence in which % of the samples appear 
cop_bysp = subset_samples(copmicro_ra_cop, Cop_sp == "C. ponticus")
cop_bysp_core_ponticus = core(cop_bysp, detection = 0.005, prevalence = 30/100)
cop_bysp = subset_samples(copmicro_ra_cop, Cop_sp == "O. nana")
cop_bysp_core_nana = core(cop_bysp, detection = 0.005, prevalence = 30/100)
cop_bysp = subset_samples(copmicro_ra, Cop_sp == "Sea water")
cop_bysp_core_Seawater = core(cop_bysp, detection = 0.005, prevalence = 30/100)


#for unique_core_microbiome
copmicro_ra_cop

asv_for_core = data.frame(otu_table(copmicro_ra_cop))
taxa_for_core = data.frame(tax_table(copmicro_ra_cop))
combined_core = cbind(taxa_for_core, asv_for_core)
combined_core$asvID = rownames(combined_core) # add a new row for the asvID 

copmicro_ra_water_df = data.frame(otu_table(copmicro_ra_water))
copmicro_ratax_water_df = data.frame(tax_table(copmicro_ra_water))
combined_water = cbind(copmicro_ratax_water_df, copmicro_ra_water_df)
combined_water$asvIDwater = rownames(combined_water)

core = read.csv(file="core.csv", header=T, na.strings="NA", stringsAsFactors=F, check.names=F, sep=",")
combined_core1 = combined_core %>% filter(combined_core$asvID %in% core$asvID)
combined_core1 = combined_core1[,-62]

combined_water = combined_water %>% filter(combined_water$asvIDwater %in% core$asvIDwater)
combined_water = combined_water[,-20]

write.table(combined_core1, file="combined_core_cop_2020.csv", sep = ",", quote = F, row.names = T, col.names = T)
write.table(combined_water,file="combined_seawater_2020.csv", sep = ",", quote = F, row.names = T, col.names = T)

forseq = read.csv(file="asvabundance_taxa_merged_copepod_microbiome_2020.csv", header=T, na.strings="NA", stringsAsFactors=F, check.names=F, sep=",")
forseq1 = forseq %>% filter(forseq$asvID %in% core$asvID)
forseq1 = forseq[,c(1:9,10:75)]
forseq2 = forseq %>% filter(forseq$asvID %in% core$asvIDwater)
forseq2 = forseq[,c(1:9,65:76)]

write.table(forseq1, file="forseqcore_cop.csv", sep = ",", quote = F, row.names = T, col.names = T)
write.table(forseq2, file="forseqcore_seawater.csv", sep = ",", quote = F, row.names = T, col.names = T)

### Metadata (samples groups-metadata in columns. Samples in the rows in the same order as in 'features' table)
uniq_core_ponticus= read.csv(file="core_ponticus.csv", header=T, na.strings="NA", stringsAsFactors=F,row.names = 1, check.names=F, sep=",")
uniq_core_nana= read.csv(file="core_nana.csv", header=T, na.strings="NA", stringsAsFactors=F,row.names = 1, check.names=F, sep=",")
uniq_core_stylifera= read.csv(file="core_stylifera.csv", header=T, na.strings="NA", stringsAsFactors=F,row.names = 1, check.names=F, sep=",")

asv_core_df = uniq_core_ponticus[,-1:-7]
taxa_core_df = uniq_core_ponticus[,1:7]
meta_ponticus
rownames(asv_core_df) = paste0(rownames(taxa_core_df), "_", taxa_core_df$Genus)

asv_core_df = uniq_core_nana[,-1:-7]
taxa_core_df = uniq_core_nana[,1:7]
meta_nana
rownames(asv_core_df) = paste0(rownames(taxa_core_df), "_", taxa_core_df$Genus)

asv_core_df = uniq_core_stylifera[,-1:-7]
taxa_core_df = uniq_core_stylifera[,1:7]
meta_stylifera
rownames(asv_core_df) = paste0(rownames(taxa_core_df), "_", taxa_core_df$Genus)

core_micro = as.data.frame(t(asv_core_df))

all(rownames(core_micro) == rownames(meta_ponticus))

#Run KW test with a for loop for each feature in the table:
test_res = c() #an empty vector to store the p-values from the KW test
mean_February_Winter = c() #an empty vector to store the mean abundances of each feature in specific group
mean_April_Spring = c()
mean_July_Summer = c()
mean_October_Autumn = c()

for (i in 1:length(core_micro)){
  print (colnames(core_micro)[i])
  test_res[i] = kruskal.test(core_micro[,i], meta_ponticus$Season)$p.value
  mean_February_Winter[i] = mean(core_micro[meta_ponticus$Season=="Winter", i])
  mean_April_Spring[i] = mean(core_micro[meta_ponticus$Season=="Spring", i])
  mean_July_Summer[i] = mean(core_micro[meta_ponticus$Season=="Summer", i])
  mean_October_Autumn[i] = mean(core_micro[meta_ponticus$Season=="Autumn", i])
}

#Combine the result vector of p-values with features names and mean values to make a df: (add as many mean_group items as needed)
test_res_df = cbind.data.frame(colnames(core_micro),mean_February_Winter,mean_April_Spring,mean_July_Summer, mean_October_Autumn, test_res, stringsAsFactors=FALSE)
colnames(test_res_df) = c("ASV", "mean_February_Winter", "mean_April_Spring", "mean_July_Summer", "mean_October_Autumn", "p_value")

#Run p.adjust() on the raw p-valuess to correct for multiple-testing, using FDR (BH, Benjamini-Hochberg): 
test_res_df$p_fdr_BH = p.adjust(test_res_df$p_value, method = "BH")
#Keep all the significant features after p-value FDR<0.05
test_res_df = test_res_df[test_res_df$p_fdr_BH < 0.05, ]

#Save the results
write.table(test_res_df, file="KW_test_core_microbiome_ponticus_month_FDR_pval_res_1.csv", sep = ",", row.names=F, col.names=T, quote = F)

###----------------------------------------------------------------------------------------------------------###

#### Create a long-format dataframe for plotting ####
library(reshape2)

core_micro$Sample_ID = rownames(core_micro)
core_micro$Season = meta_ponticus$Season
core_micro_l = melt(core_micro, id.vars = c("Sample_ID","Season"),
                    value.name = "Relative_abundance",
                    variable.name = "ASV") #convert from wide into long format

#core_micro_l$Cop_sp = factor(core_micro_l$Cop_sp, levels = c("C. ponticus","O. nana","T. stylifera","Sea water"))

ggplot(core_micro_l, aes(x=Season, y=Relative_abundance, group=Season)) +
  geom_boxplot(aes(fill=Season), outlier.shape=NA, color="black", lwd=0.6, alpha=0.5) +
  geom_jitter(position=position_jitter(0.2), aes(color = Season), stroke=0.5) +
  scale_fill_manual(values=c("mediumblue","orange","red","forestgreen")) + #this is for the boxplots "mediumblue"
  scale_color_manual(values=c("mediumblue","orange","red","forestgreen")) + #this is for the jitter points
  facet_wrap(~ ASV, ncol=3, nrow=3, scales="free",  strip.position="bottom") + #By default, all the panels have the same scales (scales="fixed"). They can be made independent, by setting scales to free  
  theme(strip.background = element_rect(colour="black", fill="black", size=1.2, linetype="solid")) + # Change the apperance of the rectangle around facet label
  theme(strip.text.x = element_text(size=9, color="white", face="bold.italic")) + # Change facet text font. Possible values for the font style: 'plain', 'italic', 'bold', 'bold.italic'.
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x  = element_blank()) + #element_text(size=10, face="bold.italic", color = "black"))+
  theme(axis.title.y = element_text(size=12, face="bold")) +
  theme(axis.text.y  = element_text(size=10, face="bold", color = "black")) +
  theme(legend.position= "0")+
  #theme(legend.text = element_text(size=10, face = "italic" ))+
  labs(y ="Relative abundance") +
  theme(axis.ticks.x = element_blank()) +
  theme(panel.background = element_rect(fill='white', colour='black'))
ggsave("core_ponticus_KW.png", width=12, height=8, dpi = 350)

##only significant kw
KW_sig= read.csv(file="significant_KW.csv", header=T, na.strings="NA", stringsAsFactors=F,row.names = 1, check.names=F, sep=",")
asv_kw = KW_sig[,-1:-7]
taxa_kw = KW_sig[,1:7]
rownames(asv_kw) = paste0(rownames(taxa_kw), "_", taxa_kw$Genus)
metakw = meta_cop_env[c(1:20,40:54),]
metakw = metakw[,c(2:3)]

kw_core = as.data.frame(t(asv_kw))

kw_core$Sample_ID = rownames(kw_core)
kw_core$Season = metakw$Season
kw_core_l = melt(kw_core, id.vars = c("Sample_ID","Season"),
                 value.name = "Relative_abundance",
                 variable.name = "ASV") #convert from wide into long format

kw_core_l$Season = factor(kw_core_l$Season, levels = c("February","April","July","October"))


ggplot(kw_core_l,aes(x=Season, y=Relative_abundance, group=Season)) +
  geom_boxplot(aes(fill=Season), outlier.shape=NA, color="black", lwd=0.6, alpha=0.5) +
  geom_jitter(position=position_jitter(0.2), aes(color = Season), stroke=0.5) +
  scale_fill_manual(values=c("mediumblue","orange","red","forestgreen")) + #this is for the boxplots "mediumblue"
  scale_color_manual(values=c("mediumblue","orange","red","forestgreen")) + #this is for the jitter points
  facet_wrap(~ ASV, ncol=5, nrow=2, scales="free",  strip.position="bottom") + #By default, all the panels have the same scales (scales="fixed"). They can be made independent, by setting scales to free  
  theme(strip.background = element_rect(colour="black", fill="black", size=1.2, linetype="solid")) + # Change the apperance of the rectangle around facet label
  theme(strip.text.x = element_text(size=9, color="white", face="bold.italic")) + # Change facet text font. Possible values for the font style: 'plain', 'italic', 'bold', 'bold.italic'.
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x  = element_blank()) + #element_text(size=10, face="bold.italic", color = "black"))+
  theme(axis.title.y = element_text(size=12, face="bold")) +
  theme(axis.text.y  = element_text(size=10, face="bold", color = "black")) +
  theme(legend.position= "0")+
  theme(legend.text = element_text(size=12, face = "plain"),legend.title = element_blank())+ 
  labs(y ="Relative abundance") +
  theme(axis.ticks.x = element_blank()) +
  theme(panel.background = element_rect(fill='white', colour='black'))
ggsave("significant_KW.pdf", width=12, height=6, dpi = 300)


##temporal core
temp_co= read.csv(file="temporal_core.csv", header=T, na.strings="NA", stringsAsFactors=F,row.names = 1, check.names=F, sep=",")
asv_temp = temp_co[,-1:-7]
taxa_temp = temp_co[,1:7]
rownames(asv_temp) = paste0(rownames(taxa_temp), "_", taxa_temp$Genus)
meta_temp = meta_cop_env[,c(1:3)]

temp_core = as.data.frame(t(asv_temp))

temp_core$Sample_ID = rownames(temp_core)
temp_core$Season = meta_temp$Season
temp_core_l = melt(temp_core, id.vars = c("Sample_ID","Season"),
                   value.name = "Relative_abundance",
                   variable.name = "ASV") #convert from wide into long format


ggplot(temp_core_l,aes(x=Season, y=Relative_abundance, group=Season)) +
  geom_boxplot(aes(fill=Season), outlier.shape=NA, color="black", lwd=0.6, alpha=0.5) +
  geom_jitter(position=position_jitter(0.2), aes(color = Season), stroke=0.5) +
  scale_fill_manual(values=c("mediumblue","orange","red","forestgreen")) + #this is for the boxplots "mediumblue"
  scale_color_manual(values=c("mediumblue","orange","red","forestgreen")) + #this is for the jitter points
  facet_wrap(~ ASV, ncol=4, nrow=2, scales="free",  strip.position="bottom") + #By default, all the panels have the same scales (scales="fixed"). They can be made independent, by setting scales to free  
  theme(strip.background = element_rect(colour="black", fill="black", size=1.2, linetype="solid")) + # Change the apperance of the rectangle around facet label
  theme(strip.text.x = element_text(size=9, color="white", face="bold.italic")) + # Change facet text font. Possible values for the font style: 'plain', 'italic', 'bold', 'bold.italic'.
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x  = element_blank()) + #element_text(size=10, face="bold.italic", color = "black"))+
  theme(axis.title.y = element_text(size=12, face="bold")) +
  theme(axis.text.y  = element_text(size=10, face="bold", color = "black")) +
  theme(legend.position= "right")+
  theme(legend.text = element_text(size=12, face = "plain"),legend.title = element_blank())+  
  labs(y ="Relative abundance") +
  theme(axis.ticks.x = element_blank()) +
  theme(panel.background = element_rect(fill='white', colour='black'))
ggsave("temporalcore.pdf", width=12, height=6, dpi = 300)

##not significant
not_sig= read.csv(file="not_significant.csv", header=T, na.strings="NA", stringsAsFactors=F,row.names = 1, check.names=F, sep=",")
asv_not = not_sig[,-1:-7]
taxa_not = not_sig[,1:7]
rownames(asv_not) = paste0(rownames(taxa_not), "_", taxa_not$Genus)
meta_not = meta_cop_env[,c(1:3)]

not_sing_core = as.data.frame(t(asv_not))

not_sing_core$Sample_ID = rownames(not_sing_core)
not_sing_core$Season = meta_temp$Season
not_sing_core_l = melt(not_sing_core, id.vars = c("Sample_ID","Season"),
                       value.name = "Relative_abundance",
                       variable.name = "ASV") #convert from wide into long format


ggplot(not_sing_core_l,aes(x=Season, y=Relative_abundance, group=Season)) +
  geom_boxplot(aes(fill=Season), outlier.shape=NA, color="black", lwd=0.6, alpha=0.5) +
  geom_jitter(position=position_jitter(0.2), aes(color = Season), stroke=0.5) +
  scale_fill_manual(values=c("mediumblue","orange","red","forestgreen")) + #this is for the boxplots "mediumblue"
  scale_color_manual(values=c("mediumblue","orange","red","forestgreen")) + #this is for the jitter points
  facet_wrap(~ ASV, ncol=3, nrow=4, scales="free",  strip.position="bottom") + #By default, all the panels have the same scales (scales="fixed"). They can be made independent, by setting scales to free  
  theme(strip.background = element_rect(colour="black", fill="black", size=1.2, linetype="solid")) + # Change the apperance of the rectangle around facet label
  theme(strip.text.x = element_text(size=10, color="white", face="bold.italic")) + # Change facet text font. Possible values for the font style: 'plain', 'italic', 'bold', 'bold.italic'.
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x  = element_blank()) + #element_text(size=10, face="bold.italic", color = "black"))+
  theme(axis.title.y = element_text(size=12, face="bold")) +
  theme(axis.text.y  = element_text(size=10, face="bold", color = "black")) +
  theme(legend.position= "right")+
  theme(legend.text = element_text(size=12, face = "plain"),legend.title = element_blank())+  
  labs(y ="Relative abundance") +
  theme(axis.ticks.x = element_blank()) +
  theme(panel.background = element_rect(fill='white', colour='black'))
ggsave("not_significant_core.pdf", width=10, height=12, dpi = 300)


#### Average % of total core microbiome in each species ####
core_stylifera_asv = data.frame(otu_table(cop_bysp_core_stylifera))
mean(colSums(core_stylifera_asv))
core_ponticus_asv = data.frame(otu_table(cop_bysp_core_ponticus))
mean(colSums(core_ponticus_asv))
core_nana_asv = data.frame(otu_table(cop_bysp_core_nana))
mean(colSums(core_nana_asv))
core_water_asv = data.frame(otu_table(cop_bysp_core_Seawater))
mean(colSums(core_water_asv))

##take off the share ASVs with the water and other copepods
core_stylifera_asv = core_stylifera_asv[c(-1),]
mean(colSums(core_stylifera_asv))*100
core_ponticus_asv = core_ponticus_asv[c(-1,-7),]
mean(colSums(core_ponticus_asv))*100
core_nana_asv =  core_nana_asv[c(-1,-6),]
mean(colSums(core_nana_asv))*100


#### Percentage of unique core microbiome in each species/month ####
core_stylifera_asv_df = as.data.frame(colSums(core_stylifera_asv))
core_ponticus_asv_df = as.data.frame(colSums(core_ponticus_asv))
core_nana_asv_df = as.data.frame(colSums(core_nana_asv))
colnames(core_stylifera_asv_df) = "core_micro"
colnames(core_ponticus_asv_df) = "core_micro"
colnames(core_nana_asv_df) = "core_micro"

#Percentage of unique core microbiome in each month
mean(core_stylifera_asv_df$core_micro[meta_stylifera$Month == "April"])*100
mean(core_stylifera_asv_df$core_micro[meta_stylifera$Month == "July"])*100
mean(core_stylifera_asv_df$core_micro[meta_stylifera$Month == "October"])*100
mean(core_ponticus_asv_df $core_micro[meta_ponticus$Month == "February"])*100
mean(core_ponticus_asv_df $core_micro[meta_ponticus$Month == "April"])*100
mean(core_ponticus_asv_df $core_micro[meta_ponticus$Month == "July"])*100
mean(core_ponticus_asv_df $core_micro[meta_ponticus$Month == "October"])*100
mean(core_nana_asv_df $core_micro[meta_nana$Month == "February"])*100
mean(core_nana_asv_df $core_micro[meta_nana$Month == "April"])*100
mean(core_nana_asv_df $core_micro[meta_nana$Month == "July"])*100
mean(core_nana_asv_df $core_micro[meta_nana$Month == "October"])*100

all_core_asv_df = rbind(core_ponticus_asv_df, core_nana_asv_df, core_stylifera_asv_df)
all(rownames(all_core_asv_df) == rownames(meta_cop_env))
all_core_asv_df$Cop_sp = meta_cop_env$Cop_sp
all_core_asv_df$Season = meta_cop_env$Season
all_core_asv_df$core_micro1 = all_core_asv_df$core_micro*100


ggplot(all_core_asv_df, aes(x=Cop_sp, y=core_micro1, group=Cop_sp)) +
  facet_grid(.~ Season, scales="free")+
  geom_boxplot(aes(fill=Cop_sp), outlier.shape=NA, color="black", lwd=0.6, alpha=0.5) +
  geom_jitter(position=position_jitter(0.2), aes(color = Cop_sp), stroke=0.5) +
  scale_y_continuous(breaks=c(0,20,40,60,80,100), limits=c(0,100))+
  scale_fill_manual(values=c("gray78","darkslategray3","gray31")) + #this is for the boxplots
  scale_color_manual(values=c("gray78","darkslategray3","gray31")) + #this is for the jitter points
  theme(strip.background = element_rect(color="black", fill="black", size=1.5, linetype="solid"))+ # Change the apperance of the rectangle around facet label
  theme(strip.text.x = element_text(size=12, color="white", face="bold")) + # Change facet text font. Possible values for the font style: 'plain', 'italic', 'bold', 'bold.italic'.
  theme(axis.title.x = element_blank()) + 
  theme(axis.ticks.x = element_blank())+
  theme(axis.text.x  = element_blank()) + #element_text(size=10, face="bold.italic", color = "black"))+
  theme(axis.title.y = element_text(size=12, face="bold")) +
  theme(axis.text.y  = element_text(size=12, face="bold", color = "black")) +
  theme(legend.position="right")+
  labs(y = "% specific core taxa") +
  theme(legend.text = element_text(size=12, face = "italic"), legend.title = element_blank())+
  theme(panel.background = element_rect(fill='white', colour='black'))  
ggsave("% of core_cop.pdf", width=10, height=5, dpi = 300)



####bubble chart####

require(gridExtra)

cop_bysp_core_ponticus
cop_bysp_core_nana
cop_bysp_core_stylifera
meta_ponticus = meta_cop_env [1:20,1:3]

asv_ponticus= data.frame(otu_table(cop_bysp_core_ponticus))
taxa_ponticus = data.frame(tax_table(cop_bysp_core_ponticus))
asv_ponticus = asv_ponticus[c(-1,-7),]
taxa_ponticus =taxa_ponticus[c(-1,-7),]
rownames(asv_ponticus) = paste0(rownames(taxa_ponticus), "_", taxa_ponticus$Genus)

ponticus = as.data.frame(t(asv_ponticus))
ponticus$Cop_sp = meta_ponticus$Cop_sp
ponticus$Season = meta_ponticus$Season
ponticus$SampleID = rownames(meta_ponticus)


ponticus_l = melt(ponticus, id.vars = c("Season","Cop_sp","SampleID"),
                  value.name = "Relative_abundance",
                  variable.name = "ASV") #convert from wide into long format

#ponticus_l$rela_log = log10(ponticus_l$Relative_abundance+0.0001)

p = ggplot(ponticus_l, aes(x = SampleID, y = ASV)) + 
  geom_point(aes(size = Relative_abundance, fill = ASV), alpha = 0.75, shape = 21) +
  #scale_fill_manual(values = colorRampPalette(brewer.pal(n = 9, name = "Greys"))(8)) +
  scale_fill_manual(values = rep("gray78", 8)) +
  scale_size_continuous(limits = c(0.0002, 0.5), breaks = c(0.0002,0.001,0.005,0.01,0.05,0.1,0.5))+
  facet_grid(.~ Season, scales="free")+
  theme(strip.background = element_rect(color="black", fill="black", size=1.5, linetype="solid"))+ # Change the apperance of the rectangle around facet label
  theme(strip.text.x = element_text(size=12, color="white", face="bold"))+ # Change facet text font. Possible values for the font style: 'plain', 'italic', 'bold', 'bold.italic'.
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x  = element_text(size=8, face="bold", color = "black", angle = 90))+
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y  = element_text(size=10, face="bold", color = "black")) +
  theme(legend.position="0")+
  #theme(legend.text = element_text(size=6, face = "italic" ))+
  theme(panel.background = element_rect(fill='white', colour='black'))
p


asv_nana= data.frame(otu_table(cop_bysp_core_nana))
taxa_nana = data.frame(tax_table(cop_bysp_core_nana))
asv_nana = asv_nana[c(-1,-6),]
taxa_nana =taxa_nana[c(-1,-6),]
rownames(asv_nana) = paste0(rownames(taxa_nana), "_", taxa_nana$Genus)
meta_nana = meta_cop_env [21:39,1:3]

nana = as.data.frame(t(asv_nana))
nana$Cop_sp = meta_nana$Cop_sp
nana$Season = meta_nana$Season
nana$SampleID = rownames(meta_nana)

nana_l = melt(nana, id.vars = c("Season","Cop_sp","SampleID"),
              value.name = "Relative_abundance",
              variable.name = "ASV") #convert from wide into long format


n = ggplot(nana_l, aes(x = SampleID, y = ASV)) + 
  geom_point(aes(size = Relative_abundance, fill = ASV), alpha = 0.75, shape = 21) +
  scale_fill_manual(values = rep("darkslategray3", 7)) +
  scale_size_continuous(limits = c(0.0002, 0.5), breaks = c(0.0002,0.001,0.005,0.01,0.05,0.1,0.5))+
  facet_grid(.~ Season, scales="free")+
  theme(strip.background = element_rect(color="black", fill="black", size=1.5, linetype="solid"))+ # Change the apperance of the rectangle around facet label
  theme(strip.text.x = element_text(size=12, color="white", face="bold"))+ # Change facet text font. Possible values for the font style: 'plain', 'italic', 'bold', 'bold.italic'.
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x  = element_text(size=8, face="bold", color = "black",angle = 90))+
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y  = element_text(size=10, face="bold", color = "black")) +
  theme(legend.position="0")+
  #theme(legend.text = element_text(size=6, face = "italic" ))+
  theme(panel.background = element_rect(fill='white', colour='black'))
n


asv_stylifera= data.frame(otu_table(cop_bysp_core_stylifera))
taxa_stylifera = data.frame(tax_table(cop_bysp_core_stylifera))
asv_stylifera = asv_stylifera[c(-1),]
taxa_stylifera =taxa_stylifera[c(-1),]
rownames(asv_stylifera) = paste0(rownames(taxa_stylifera), "_", taxa_stylifera$Genus)
meta_stylifera = meta_cop_env [40:54,1:3]

stylifera = as.data.frame(t(asv_stylifera))
stylifera$Cop_sp = meta_stylifera$Cop_sp
stylifera$Season = meta_stylifera$Season
stylifera$SampleID = rownames(meta_stylifera)

stylifera_l = melt(stylifera, id.vars = c("Season","Cop_sp","SampleID"),
                   value.name = "Relative_abundance",
                   variable.name = "ASV") #convert from wide into long format

s = ggplot(stylifera_l, aes(x = SampleID, y = ASV)) + 
  geom_point(aes(size = Relative_abundance, fill = ASV), alpha = 0.75, shape = 21) +
  scale_fill_manual(values = rep("gray31", 19)) +
  scale_size_continuous(limits = c(0.0002, 0.5), breaks = c(0.0002,0.001,0.005,0.01,0.05,0.1,0.5))+
  facet_grid(.~ Season, scales="free")+
  theme(strip.background = element_rect(color="black", fill="black", size=1.5, linetype="solid"))+ # Change the apperance of the rectangle around facet label
  theme(strip.text.x = element_text(size=12, color="white", face="bold"))+ # Change facet text font. Possible values for the font style: 'plain', 'italic', 'bold', 'bold.italic'.
  theme(axis.title.x = element_blank()) + 
  theme(axis.text.x  = element_text(size=8, face="bold", color = "black",angle = 90))+
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y  = element_text(size=10, face="bold", color = "black")) +
  theme(legend.position= "0")+
  #theme(legend.text = element_text(size = 10)) +
  theme(panel.background = element_rect(fill='white', colour='black'))
s
#ggsave("for_legent.pdf", plot=s, width=10, height=12, dpi = 300)

gr= grid.arrange(p, n, s, ncol=1, nrow=3)
ggsave("core_cop_microbiome.pdf", plot=gr, width=10, height=12, dpi = 300)