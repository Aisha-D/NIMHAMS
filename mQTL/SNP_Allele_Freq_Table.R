## This script creates an allele frquency chart for different population based on selected number of SNPs

####      Read in files     ####
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringi)

setwd('/gpfs/ts0/scratch/and202/NIMHAMS/mQTL/Output/')
mydir = "/gpfs/ts0/scratch/and202/NIMHAMS/mQTL/Output/"
myfiles = list.files(pattern="*.txt")
myfiles = myfiles[-1]

#Read multiple files from each chr
for(i in 1:22) {
  chr <- myfiles[i]
  chr_ob <- strsplit(chr, ".txt") #removes the .txt taken from file name
  tmp <- read.table(paste("Ind_chr",i,".txt", sep=""), stringsAsFactors = F, header = T) #read in the file as a dataframe
  chr_ob <- chr_ob[[1]] #get the dataframe from the list
  assign(chr_ob, tmp) #assigns the file name tmp to the object chr_ob which is the dataframe
  chr_ob <- NULL
  tmp <- NULL}
rm(chr, chr_ob, i, mydir, myfiles, tmp)

#now we rbind all chr as a dataframe in a loop
all_chr <- Ind_chr1
for (i in 2:22){
  tmp <- paste("Ind_chr",i, sep = "")
  tmp <- get(tmp)
  all_chr  <- rbind(all_chr, tmp)
}
rm(list=setdiff(ls(), 'all_chr'))

##Add annotation files
epicManifest <- read.csv("/gpfs/ts0/scratch/and202/NIMHAMS/epicManifest_hg38.csv", 
                         stringsAsFactors = F, header = T, sep =" ")

genotypeManifest <- read.csv("/gpfs/ts0/scratch/gn261/Nimhams/Infinium_Global_Screening-24_v2.0_manifest.csv")

##Read in the european samples that had mQTLs run on
load("/gpfs/ts0/scratch/and202/UKDataset/AllAutosomeResults_mQTLs_CER_Anno_with2PCs.Rdata")






####      Find SNPs of interest     ####
#check which cpgs are in both datasets
length(intersect(res.cer.anno$gene, all_chr$gene)) #2286 cpgs
#Take the snps that do not intersect in both datasets

#SNPs only in India. Which snps not in european mqtl
unique_ind_snps = all_chr[-which(all_chr$SNP %in% res.cer.anno$SNP),] #8885
#remove the a/c/G/T on the SNP name then read out as text.file
unique_ind_snps$SNP <- gsub("_A", "", unique_ind_snps$SNP)
unique_ind_snpsSNP <- gsub("_T", "", unique_ind_snps$SNP)
unique_ind_snps$SNP <- gsub("_G", "", unique_ind_snps$SNP)
unique_ind_snps$SNP <- gsub("_C", "", unique_ind_snps$SNP)
unique_ind_snps <- unique_ind_snps[order(unique_ind_snps$FDR),]
ind_snps = unique_ind_snps$SNP

# write.table(ind_snps, 'SNPs_only_in_india_mqtl.txt', row.names = F, quote = F, col.names =F)

#save the snps that are present in ind_snps. These snps will be used to load the allele frequency in differe
#populations. The code to do this is: cat SNPs_only_in_india_mqtl.txt | while read R ; do wget -q -O - "https://www.ncbi.nlm.nih.gov/snp/${R}?download=frequency" | grep -E '^(#Study|1000Genomes)' | sed "s/^/${R}\t/" >> snp_freq_1000G.txt ; done

####      Visualise Allele Frequency In SA and Euro populations       #####

##read in allele frequency
ind_snp_freq = read.table("snp_allele_freq/snp_freq_1000G.txt",
                          sep = "\t", fill = TRUE, stringsAsFactors = F) #36756
ind_snp_freq <- ind_snp_freq[!duplicated(ind_snp_freq), ]

# make an index eg. every 7th to be removed as there are blank rows and remove from dataframe
ind <- seq(1, nrow(ind_snp_freq), by=7)
ind_snp_freq = ind_snp_freq[-ind, ] #31505

rownames(ind_snp_freq) <- c(1:nrow(ind_snp_freq))
colnames(ind_snp_freq) <- c("SNP", "Study", "Population", "Group", "SampleSize", "RefAllele",
                            "AltAllele", "BioProjectID", "BiosSampleID")
ind_snp_freq<- ind_snp_freq[which(ind_snp_freq$Population == c("Europe" , "South Asian")),]
#create two tables for each population then combine them later
ind_snp_freq_euro <- ind_snp_freq[which(ind_snp_freq$Population == "Europe"),]
ind_snp_freq_sa <- ind_snp_freq[which(ind_snp_freq$Population == "South Asian"),]


##make combinations of allele frequency in respective tables. This will be joined to the table later on.
#Euro
ind_snp_freq_euro$Freq <- rep(NA, nrow(ind_snp_freq_euro))
ind_snp_freq_euro$Freq = paste( gsub("[^0-9.]",'',ind_snp_freq_euro$RefAllele),
                                gsub("[^0-9.]",'',ind_snp_freq_euro$AltAllele),
                                sep = "/")
#SA
ind_snp_freq_sa$Freq <- rep(NA, nrow(ind_snp_freq_sa))
ind_snp_freq_sa$Freq = paste( gsub("[^0-9.]",'',ind_snp_freq_sa$RefAllele),
                              gsub("[^0-9.]",'',ind_snp_freq_sa$AltAllele),
                              sep = "/")


##Now merge two population to one
ind_snp_freq_tbl <- merge(ind_snp_freq_euro,ind_snp_freq_sa, by = "SNP", all= TRUE)

ind_snp_freq_tbl$AlleleComb = rep(NA, nrow(ind_snp_freq_tbl))
ind_snp_freq_tbl$AlleleComb = paste(stri_extract_first_regex(ind_snp_freq_tbl$RefAllele.x, ".{1}"),
                                    stri_extract_first_regex(ind_snp_freq_tbl$AltAllele.x, ".{1}"),
                                    sep = "/")
colnames(ind_snp_freq_tbl)[colnames(ind_snp_freq_tbl) == 'Freq.x'] <- 'Europe'
colnames(ind_snp_freq_tbl)[colnames(ind_snp_freq_tbl) == 'Freq.y'] <- 'South.Asia'
rem = c("Population.x", "Population.y", "SampleSize.x", "SampleSize.y",
        "BioProjectID.x", "BiosSampleID.x", "BioProjectID.y", "BiosSampleID.y",
        "Group.x", "Group.y", "Study.y", "Study.x", "RefAllele.x", "RefAllele.y", 
        "AltAllele.x", "AltAllele.y")
ind_snp_freq_tbl = ind_snp_freq_tbl[, -which(names(ind_snp_freq_tbl) %in% rem)]


ind_snp_freq_tbl2<- cbind(ind_snp_freq_tbl, genotypeManifest[match(ind_snp_freq_tbl$SNP, genotypeManifest$Name), 
                                                             c("Chr", "MapInfo")])
ind_snp_freq_tbl2$Chr <- as.numeric(as.character(ind_snp_freq_tbl2$Chr))
ind_snp_freq_tbl2$Location <- rep(NA, nrow(ind_snp_freq_tbl2))
ind_snp_freq_tbl2$Location <- paste(ind_snp_freq_tbl2$Chr, ind_snp_freq_tbl2$MapInfo, sep = ":")
ind_snp_freq_tbl2 = ind_snp_freq_tbl2[, -which(names(ind_snp_freq_tbl2) %in% c("Chr", "MapInfo"))]

ind_snp_freq_tbl2 = ind_snp_freq_tbl2[c("SNP", "Location", "AlleleComb", "Europe", "South.Asia")]
colnames(ind_snp_freq_tbl2)[colnames(ind_snp_freq_tbl2) == 'AlleleComb'] <- 'Allele'
colnames(ind_snp_freq_tbl2)[colnames(ind_snp_freq_tbl2) == 'South.Asia'] <- 'South Asia'
rownames(ind_snp_freq_tbl2) <- NULL
pdf('Ind SNP frequency in 1KG.pdf', w = 8, h = 8)
grid.table(ind_snp_freq_tbl2[1:20,], rows = NULL)
dev.off()
