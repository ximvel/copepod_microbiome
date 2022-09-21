library(Biostrings)
library(dada2)

#### Remove primers with Biostrings ####

path="C:/Rspace/Microbiome_project/Data_seq_2020_copepod_microbiome/" ##input row data
names_F=list.files(path=path, pattern="_R1.", recursive=T)
names_R=list.files(path=path, pattern="_R2.", recursive=T)

#Your primer sequence used in the amplification
fp=DNAString("GTGYCAGCMGCCGCGGTAA") #Forward primer for 16S
rp=DNAString("GGACTACNVGGGTWTCTAAT") #Reverse primer for 16S 
stats.R=stats.F=data.frame(Raw.num=numeric(),With.primer=numeric(),final.num=numeric()) #Table for reads statistics (numbers) after each trimming step:

#Create a modified paths to fastq files without the subdirectories included:
names_F_mod = sapply(strsplit(basename(names_F), "-"),"[", 2)
names_R_mod = sapply(strsplit(basename(names_R), "-"),"[", 2)

#Loop to trim (remove reads with >2 mismatches to the primer seq and cut the primer squence from the reads) all the R1-R2 files in your dir:
for (i in 1:length(names_F)){
  
  FR=readQualityScaledDNAStringSet(paste0(path,names_F[i]),quality.scoring = "illumina")
  RR=readQualityScaledDNAStringSet(paste0(path,names_R[i]),quality.scoring = "illumina")
  stats.R[i,1]=length(RR)
  stats.F[i,1]=length(FR)
  print(paste("read files", i))
  keepF=which(vcountPattern(fp,FR,max.mismatch=2,fixed=FALSE)>0) #check if the primer sequence is found in the forward reads
  keepR=which(vcountPattern(rp,RR,max.mismatch=2,fixed=FALSE)>0) #check if the primer sequence is found in the reverse reads
  keep=intersect(keepF,keepR) #keep ONLY seqs that have both forward (on FR) and reverse (on RR) primers
  FR=FR[keep]
  RR=RR[keep]
  
  stats.F[i,2]=length(keepF) #How many redas left after removing of the the reads with with >2 mismatches
  stats.R[i,2]=length(keepR) #How many redas left after removing of the the reads with with >2 mismatches
  
  ### Now remove the actual primers sequence from the reads:
  FR.s=subseq(FR,19,width(FR)) #strip off seq primers, the number is the length of your primer
  RR.s=subseq(RR,20,width(RR)) #strip off seq primers, the number is the length of your primer
  
  stats.F[i,3] = length(FR.s) # How many reads left after removing the primers 
  stats.R[i,3] = length(RR.s) # How many reads left after removing the primers 
  rownames(stats.F)[i]=rownames(stats.R)[i]=gsub("_R1.fastq","",names_F_mod[i])
  
  #Save in the computer the new trimmed FASTQ files:
  writeQualityScaledXStringSet(FR.s,paste0(path,"primers_trimmed/",names_F_mod[i]),compress=TRUE)
  writeQualityScaledXStringSet(RR.s,paste0(path,"primers_trimmed/",names_R_mod[i]),compress=TRUE)
}

#For reads statistics (# of reads after each step)
#Add more columns with reads fraction after each step
stats.F$per.with.primer = stats.F$With.primer / stats.F$Raw.num
stats.F$per.final = stats.F$final / stats.F$Raw.num
stats.R$per.with.primer = stats.R$With.primer / stats.R$Raw.num
stats.R$per.final = stats.R$final / stats.R$Raw.num

write.csv(stats.R,paste0(path,"stats.R.csv"))
write.csv(stats.F,paste0(path,"stats.F.csv"))


#### DADA2 Analysis ####

path="C:/Rspace/Microbiome_project/Data_seq_2020_microbiome/November_seawater_dalit/primers_trimmed/"
fnFs = sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs = sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names = sapply(strsplit(basename(fnFs), "_"), "[", 1)

#visualizing the quality profiles of the forward and reverse reads (before quality trimming):
plotQualityProfile(fnFs[1:length(fnFs)]) #forward
plotQualityProfile(fnRs[1:length(fnFs)]) #reverse

### Filter and trim

# Save the filtered files in filtered/ subdirectory
filtFs = file.path(path, "filtered_tight", paste0(sample.names, "_F_filt.fastq.gz")) #Assign the filenames for the filtered fastq.gz files
filtRs = file.path(path, "filtered_tight", paste0(sample.names, "_R_filt.fastq.gz")) #Assign the filenames for the filtered fastq.gz files

#Trimming itself:
out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,190),
                    maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, multithread=FALSE)

out=as.data.frame(out) #convert to df the reads stats table

#visulize reads quality after trimming:
plotQualityProfile(filtFs[1:4]) #forward
plotQualityProfile(filtRs[1:4]) #reverse

#Learn the Error Rates for your specific samples
errF = learnErrors(filtFs, multithread=FALSE)
errR = learnErrors(filtRs, multithread=FALSE)

#visualize the estimated error rates:
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)


### Dereplication
#Dereplication combines all identical sequencing reads into "unique sequences" with a corresponding "abundance" equal to the number of reads with that unique sequence. Dereplication substantially reduces computation time by eliminating redundant comparisons.
derepFs = derepFastq(filtFs, verbose=TRUE)
derepRs = derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) = sample.names
names(derepRs) = sample.names


### Sample Inference
#We are now ready to apply the core sample inference algorithm to the dereplicated data
#Now run the dada magic - the error models (errR and errF) are used to seperate true seqeuences from sequencing errors. The derpelicated seqeunecs (derepF,derepR) are the inputs + the error models (errF,errR)

dadaFs=dada(derepFs, err=errF, multithread = F)
dadaRs=dada(derepRs, err=errR, multithread = F)

#inspect the DADA object (that contains the actual amplicon sequence variants [ASV]), see how many sequence variants were inferred in each sample
dadaFs[3] #[n] - the number of samples


### Merge paired reads
#We now merge the forward and reverse reads together to obtain the full denoised sequences (forward + reverse complement="contig"):
mergers = mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = 8, maxMismatch = 0, verbose=TRUE) #this is default values
head(mergers[1]) #inspect the results - sequence variants and summary table


### Construct sequence table
#We can now construct an amplicon sequence variant table (ASV) table, a higher-resolution version of the OTU table produced by traditional methods

seqtab = makeSequenceTable(mergers) # Make data to a classic sample abundance table. Note that the colnames are the actual DNA seqeunces! 
dim(seqtab) #number of samples, number of ASVs

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#to remove sequences that are shorter or longer then expected because of unspecific primings:
#seqtab2 = seqtab[,nchar(colnames(seqtab)) %in% seq(250,256)


### Remove chimeras
seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=F, verbose=TRUE)
dim(seqtab.nochim)
#Percent of ASVs remained after chimera removeal
ncol(seqtab.nochim) / ncol(seqtab) * 100

saveRDS(seqtab.nochim,paste0(path,"seqtab.nochim.rds")) #this is saved as a dada2 object - so you can reload it later to R if you want, for example,to merge it with another dada table.


### Track reads through the pipeline
#As a final check of our progress, we'll look at the number of reads that made it through each step in the pipeline:
getN = function(x) sum(getUniques(x))
track = cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim") #denoised=dereplicated
rownames(track) = sample.names
write.csv(track, paste0(path,"DADA2_reads_tracking_all_steps.csv"))


### Assign taxonomy
taxa = assignTaxonomy(seqtab.nochim, "C:/Rspace/Microbiome_project/SILVA_db/silva_nr99_v138_train_set.fa.gz", multithread=F, tryRC=F)
taxa = addSpecies(taxa, "C:/Rspace/Microbiome_project/SILVA_db/silva_species_assignment_v138.fa.gz", allowMultiple = FALSE) #Note many seqs will have several possible species matches. If you set allowMultiple to TRUE, you get all possible matches. If you set it to FALSE, then sequences that have more than one possible match will be assigned only down to Genus level

taxa.print = taxa # Removing sequence rownames for display only
rownames(taxa.print) = NULL
head(taxa.print)

#Rearrange objects (taxonomy sequence and abundance of each ASV in each sample) for easy use:
identical(rownames(taxa),colnames(seqtab.nochim)) #verify the objects are aligned properly
taxa = as.data.frame(taxa) #covert to a dataframe
taxa$Sequence = rownames(taxa)
rownames(taxa) = NULL #More efficient - DNA info is retained in "Sequence" column
taxa$asvID = paste0("asv", seq(1,nrow(taxa))) # assign an "asv" (accurate sequence variant) ID for each seqeunce
colnames(seqtab.nochim) = taxa$asvID #replace seqs with IDs in the sample abundance table
seqtab.nochim.t = as.data.frame(t(seqtab.nochim))
identical(rownames(seqtab.nochim.t), taxa$asvID) #confirm that the number of rows (ASVs) is the same in the taxonomy and in the abunadance table
combined_taxa_abundance = cbind(taxa, seqtab.nochim.t)

write.csv(taxa, paste0(path,"DADA2_taxonomy_copepod_microbiome_2020.csv"),row.names = FALSE) #save the the sequence and taxa information
write.csv(seqtab.nochim.t, paste0(path,"DADA2_ASV_abundance_copepod_microbiome_2020.csv"))  #save the sample abundance information
write.csv(combined_taxa_abundance, paste0(path,"DADA2_ASV_taxa_abundance_merged_copepod_microbiome_2020.csv"))  #save combined taxa+abundance table


