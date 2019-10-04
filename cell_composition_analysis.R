####Cell Composition in Prefrontal Cortex India brians ###

library(minfi)
require(IlluminaHumanMethylationEPICmanifest)
library(wateRmelon)
library(dplyr)
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("FlowSorted.DLPFC.450k")

##loading data in and sorting data frames
setwd("/mnt/data1/NIMHANS/Down_stream_analysis/")
load("NIMHAM_FPC_Normalised_snpsremoved.rdat")

pheno <- read.csv("/mnt/data1/NIMHANS/EPIC/NIMHAM_SamplesPassedQC.csv")
pheno$Basename <-paste(pheno$Chip, pheno$Position, sep = "_") #adding basename column to pheno file
pheno <- pheno[grep("f", pheno$SampleType),] 
pheno <- select(pheno, -pFilterPass)


#### Making RGset with prefrontal only idats for cell composotion count ######

#Basename <- pheno$Basename #making vector of basenames
# 
# #all idats in vector
#idat_subset <- list.files("/mnt/data1/NIMHANS/epigenetic_unzip/idats/")
# 
# #selecting prefronal cortex idats
#idat_subset2 <- unique(grep(paste(Basename,collapse="|"), 
 #                            idat_subset, value=TRUE))
# 
# #copying selected idats into new directory called prefrontal cortex idats
#setwd("/mnt/data1/NIMHANS/epigenetic_unzip/idats/")
#file.copy(idat_subset2, "/mnt/data1/NIMHANS/Down_stream_analysis/Prefrontal_cortex_idats/")
# 
# #making RGset
#idatPath<-c("/mnt/data1/NIMHANS/Down_stream_analysis/Prefrontal_cortex_idats/")
#RGset <- read.metharray.exp(base = idatPath, targets = pheno, force=T)
#save(RGset, file="NIMHAMS_prefrontal_cortex_RGset.rdat")


## Load RGset for cell ccounts 
load("NIMHAMS_prefrontal_cortex_RGset.rdat")

#Turn everything in the data frame into characters so the estimate cell count fucntion works
colData(RGset)$Sample_Name <- as.character(colData(RGset)$Sample_Name)
colData(RGset)$Position <- as.character(colData(RGset)$Position)
colData(RGset)$Sex1 <- as.character(colData(RGset)$Sex1)
colData(RGset)$Sex <- as.character(colData(RGset)$Sex)
colData(RGset)$predictedSex <- as.character(colData(RGset)$predictedSex)
colData(RGset)$SampleType <- as.character(colData(RGset)$SampleType)
colData(RGset)$ID <- as.character(colData(RGset)$ID)
colData(RGset)$Intensity <- as.character(colData(RGset)$Intensity)
colData(RGset)$predictedduplicates <- as.character(colData(RGset)$predictedduplicates)

#Making cell counts file using referance data set, this takes a while so has already been done
#counts <- estimateCellCounts(RGset, compositeCellType ="DLPFC", cellTypes = c("NeuN_neg", "NeuN_pos"), processMethod = "auto")
#save(counts, file="EXTEND16_CellCounts.rdat")

load("(EXTEND16_CellCounts.rdat")
counts<-counts[match(pheno$Basename, rownames(counts)),]
png("Boxplot_of_cell_counts.png")
boxplot(counts, use.cols=T, main="Boxplot of Estimated Cell Proportions in Prefrontal Cortex Samples")
dev.off()
