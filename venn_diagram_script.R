#This script produces venn diagrams of the siginficant probes across brain regions

##########################
####   LOAD DATA    ######
##########################
#The output for the linear models will be loaded

#Cerebellum India
cer_ind <- read.csv("res_c_chip_included_manifest.csv",stringsAsFactors = F, header=T, row.names = 1)
#Frontal India
fron_ind <- read.csv("res_f__chip_included_manifest.csv", stringsAsFactors = F, header=T, row.names = 1)

#Cerebellum UK
cer_uk <- read.csv("LM_MLM_Age_Sex_CovariedChipNeuronalProportion_MLMethod.csv", header= T, row.names = 1)
#Frontal UK
fron_uk <- read.csv("LM_AgeSq_Age_Sex_Cov_Chip_NeuronalProportion.csv", header =T, row.names = 1)



#######   AGE    #########


#India frontal + cereb
fron_ind_sig <-  fron_ind[which(fron_ind$Age_P <= 0.05),]
cer_ind_sig <- cer_ind[which(cer_ind$Age_P <= 0.05),] 

#UK frontal _ cereb
fron_uk_sig <-  fron_uk[which(fron_uk$Age.P <= 0.05),]
cer_uk_sig <- cer_uk[which(cer_uk$Age.P <= 0.05),]

library(limma)
library(VennDiagram)
set1 <- rownames(fron_ind_sig)
set2 <- rownames(cer_ind_sig)
set3 <- rownames(fron_uk_sig)
set4 <- rownames(cer_uk_sig)

universe <- sort(unique(c(set1, set2, set3, set4)))

Counts <- matrix(0, nrow=length(universe), ncol=4)
# Populate the said matrix
for (i in 1:length(universe)) {
  Counts[i,1] <- universe[i] %in% set1
  Counts[i,2] <- universe[i] %in% set2
  Counts[i,3] <- universe[i] %in% set3
  Counts[i,4] <- universe[i] %in% set4
}

# Name the columns with the sample names
colnames(Counts) <- c("India Frontal","India Cerebellum", "UK Frontal", "UK Cerebellum")

# Specify the colors for the sets
cols<-c("Red", "Green", "Blue", "Yellow")
pdf('venn diagram Age.pdf')
vennDiagram(vennCounts(Counts), circle.col=cols, main = "Significant Probes across age")
dev.off()



#######   SEX    #########


#India frontal + cereb
fron_ind_sig <-  fron_ind[which(fron_ind$Sex_P <= 0.05),]
cer_ind_sig <- cer_ind[which(cer_ind$Sex_P <= 0.05),] 

#UK frontal _ cereb
fron_uk_sig <-  fron_uk[which(fron_uk$Sex.P <= 0.05),]
cer_uk_sig <- cer_uk[which(cer_uk$Sex.P <= 0.05),]

library(limma)
library(VennDiagram)
set1 <- rownames(fron_ind_sig)
set2 <- rownames(cer_ind_sig)
set3 <- rownames(fron_uk_sig)
set4 <- rownames(cer_uk_sig)

universe <- sort(unique(c(set1, set2, set3, set4)))

Counts <- matrix(0, nrow=length(universe), ncol=4)
# Populate the said matrix
for (i in 1:length(universe)) {
  Counts[i,1] <- universe[i] %in% set1
  Counts[i,2] <- universe[i] %in% set2
  Counts[i,3] <- universe[i] %in% set3
  Counts[i,4] <- universe[i] %in% set4
}

# Name the columns with the sample names
colnames(Counts) <- c("India Frontal","India Cerebellum", "UK Frontal", "UK Cerebellum")

# Specify the colors for the sets
cols<-c("Red", "Green", "Blue", "Yellow")
pdf('venn diagram Sex.pdf')
vennDiagram(vennCounts(Counts), circle.col=cols, main = "Significant Probes across sex")
dev.off()
