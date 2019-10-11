library(lme4)
setwd("/mnt/data1/NIMHANS/Down_stream_analysis/NIMHAMSg/EWAS/")
load("/mnt/data1/NIMHANS/Down_stream_analysis/NIMHAM_Normalised_snpsremoved.rdat")

pheno <- pheno[order(rownames(pheno)),]
betas <- betas[,order(colnames(betas))]
print(identical(rownames(pheno), colnames(betas)))

pheno <-pheno[-grep("203490640029_R02C01", rownames(pheno)),]
betas <- betas[,-grep("203490640029_R02C01", colnames(betas))]


pheno$ID <-as.character(pheno$ID)
pheno$Age <- as.character(as.character(pheno$Age))


##Mixed effect model for sample type ####

res<-matrix(data = NA, ncol = 3, nrow = nrow(betas))
rownames(res)<-rownames(betas)
colnames(res)<-c("SampleType_Beta", "SampleType_SE", "SampleType_P")
for(i in 1:nrow(betas)){

  modellmer<-lmer(betas[i,] ~ pheno$SampleType + pheno$Sex + pheno$Age + as.factor(pheno$Chip) + (1|pheno$ID), REML = FALSE, subset = which(is.na(pheno$SampleType) == FALSE))
  res[i,1]<-fixef(modellmer)["pheno$SampleTypef"]
  res[i,2]<-summary(modellmer)$coefficients["pheno$SampleTypef",2]
  model.null<-lmer(betas[i,] ~  pheno$Sex + pheno$Age + as.factor(pheno$Chip) + (1|pheno$ID), REML = FALSE, subset = which(is.na(pheno$SampleType) == FALSE))
  res[i,3]<-anova(modellmer,model.null)["modellmer","Pr(>Chisq)"]
}
write.csv(res, file ="EWAS_Brain_Regions.csv")

# As the beta values here represent proportion of DNA methylation (i.e. they lie between 0 and 1), the regression coefficients represent the change in proportion. 
# Typically we report our findings on the % scale therefore we will multiple the regression coefficients and SE by 100. 
# This will need to be done for all variables for which you saved the results.


res[,"SampleType_Beta"]<-res[,"SampleType_Beta"]*100
res[,"SampleType_SE"]<-res[,"SampleType_SE"]*100

### Annotate the output

# At this stage we will be interested in adding information regarding where each probe is located, and what genes or regulatory features it overlaps with. 
# There is a lot of annotation available, the code below only takes a subset of columns which I think may be most relevant. 


epicManifest<-read.csv("/mnt/data1/IndiaWorkshop/epicManifest.csv", stringsAsFactors = F, header=T)


epicManifest<-epicManifest[match(rownames(res), epicManifest$Name),]
res<-cbind(res, epicManifest[,c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DNase_Hypersensitivity_NAME", "OpenChromatin_NAME", "TFBS_Evidence_Count", "Methyl450_Loci", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")]) ## this can be edited to include additional columns or exclude as required


## QQ plot 
pdf("/mnt/data1/NIMHANS/Down_stream_analysis/NIMHAMSg/EWAS/Plots/QQplot_of_EWAS_SampleType.pdf")
qq(res$SampleType_P)
dev.off()


##Manhattan Plot
res<-res[which(res$CHR != "Y"),] ## to also exclude the X chromosome repeat and edit this line

res$CHR<-as.character(res$CHR)
res$CHR[which(res$CHR == "X")]<-23 ## to also recode Y chromosome repeat and edit this line
res$CHR<-as.numeric(res$CHR)
res<-res[which(res$CHR != ""),]

bonfP<-0.05/nrow(res)

pdf("/mnt/data1/NIMHANS/Down_stream_analysis/NIMHAMSg/EWAS/Plots/Manhattan_plot_Mixed_Effect_Model_SampleType.pdf")
manhattan(res, p = "Status_P", bp = "MAPINFO", chr = "CHR", genomewide = -log10(bonfP), suggestiveline = -log10(5e-5), logp=T, col=c("blue","yellow"), ylim=c(0,90))
dev.off()

### Mixed effect model of sex #####

res1<-matrix(data = NA, ncol = 3, nrow = nrow(betas))
rownames(res1)<-rownames(betas)
colnames(res1)<-c("Sex_Beta", "Sex_SE", "Sex_P")
  
for(i in 1:nrow(betas)){
  modellmer1<-lmer(betas[i,] ~ pheno$SampleType + pheno$Sex + pheno$Age + as.factor(pheno$Chip) + (1|pheno$ID), REML = FALSE, subset = which(is.na(pheno$Sex) == FALSE))
  res1[i,1]<-fixef(modellmer1)["pheno$SexM"]
  res1[i,2]<-summary(modellmer1)$coefficients["pheno$SexM",2]
  model.null1<-lmer(betas[i,] ~  pheno$SampleType + pheno$Age + as.factor(pheno$Chip) + (1|pheno$ID), REML = FALSE, subset = which(is.na(pheno$Sex) == FALSE))
  res1[i,3]<-anova(modellmer1,model.null1)["modellmer1","Pr(>Chisq)"]
}

write.csv(res1, file="/mnt/data1/NIMHANS/Down_stream_analysis/NIMHAMSg/EWAS/EWAS_Sex_analysis.csv")

# As the beta values here represent proportion of DNA methylation (i.e. they lie between 0 and 1), the regression coefficients represent the change in proportion. 
# Typically we report our findings on the % scale therefore we will multiple the regression coefficients and SE by 100. 
# This will need to be done for all variables for which you saved the results.


res[,"Sex_Beta"]<-res[,"Sex_Beta"]*100
res[,"Sex_SE"]<-res[,"Sex_SE"]*100

### Annotate the output

# At this stage we will be interested in adding information regarding where each probe is located, and what genes or regulatory features it overlaps with. 
# There is a lot of annotation available, the code below only takes a subset of columns which I think may be most relevant. 


#epicManifest<-read.csv("/mnt/data1/IndiaWorkshop/epicManifest.csv", stringsAsFactors = F, header=T)


epicManifest<-epicManifest[match(rownames(res1), epicManifest$Name),]
res1<-cbind(res1, epicManifest[,c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DNase_Hypersensitivity_NAME", "OpenChromatin_NAME", "TFBS_Evidence_Count", "Methyl450_Loci", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")]) ## this can be edited to include additional columns or exclude as required


## QQ plot 
pdf("/mnt/data1/NIMHANS/Down_stream_analysis/NIMHAMSg/EWAS/Plots/QQplot_of_EWAS_Sex.pdf")
qq(res1$Sex_P)
dev.off()


##Manhattan Plot
res1<-res1[which(res1$CHR != "Y"),] ## to also exclude the X chromosome repeat and edit this line
res1<-res1[which(res1$CHR != "X"),] ## to also exclude the X chromosome repeat and edit this line


# res1$CHR<-as.character(res1$CHR)
# res1$CHR[which(res1$CHR == "X")]<-23 ## to also recode Y chromosome repeat and edit this line
res1$CHR<-as.numeric(res1$CHR)
res1<-res1[which(res1$CHR != ""),]

bonfP<-0.05/nrow(res1)

pdf("/mnt/data1/NIMHANS/Down_stream_analysis/NIMHAMSg/EWAS/Plots/Manhattan_plot_of_Sex.pdf")
manhattan(res1, p = "Status_P", bp = "MAPINFO", chr = "CHR", genomewide = -log10(bonfP), suggestiveline = -log10(5e-5), logp=T, col=c("blue","yellow"), ylim=c(0,50))
dev.off()

### Mixed effect model of Age #####

res2<-matrix(data = NA, ncol = 3, nrow = nrow(betas))
rownames(res2)<-rownames(betas)
colnames(res2)<-c("Age_Beta", "Age_SE", "Age_P")

for(i in 1:nrow(betas)){
  modellmer2<-lmer(betas[i,] ~ pheno$SampleType + pheno$Sex + pheno$Age + as.factor(pheno$Chip) + (1|pheno$ID), REML = FALSE, subset = which(is.na(pheno$Sex) == FALSE))
  res2[i,1]<-fixef(modellmer2)["pheno$Age"]
  res2[i,2]<-summary(modellmer2)$coefficients["pheno$SAge",2]
  model.null2<-lmer(betas[i,] ~  pheno$SampleType + pheno$Sex + as.factor(pheno$Chip) + (1|pheno$ID), REML = FALSE, subset = which(is.na(pheno$Sex) == FALSE))
  res2[i,3]<-anova(modellmer2,model.null2)["modellmer2","Pr(>Chisq)"]
}

write.csv(res1, file="/mnt/data1/NIMHANS/Down_stream_analysis/NIMHAMSg/EWAS/EWAS_Age_analysis.csv")

# As the beta values here represent proportion of DNA methylation (i.e. they lie between 0 and 1), the regression coefficients represent the change in proportion. 
# Typically we report our findings on the % scale therefore we will multiple the regression coefficients and SE by 100. 
# This will need to be done for all variables for which you saved the results.


res[,"Age_Beta"]<-res[,"Age_Beta"]*100
res[,"Age_SE"]<-res[,"Age_SE"]*100

### Annotate the output

# At this stage we will be interested in adding information regarding where each probe is located, and what genes or regulatory features it overlaps with. 
# There is a lot of annotation available, the code below only takes a subset of columns which I think may be most relevant. 


#epicManifest<-read.csv("/mnt/data1/IndiaWorkshop/epicManifest.csv", stringsAsFactors = F, header=T)


epicManifest<-epicManifest[match(rownames(res2), epicManifest$Name),]
res2<-cbind(res2, epicManifest[,c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DNase_Hypersensitivity_NAME", "OpenChromatin_NAME", "TFBS_Evidence_Count", "Methyl450_Loci", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")]) ## this can be edited to include additional columns or exclude as required


## QQ plot 
pdf("/mnt/data1/NIMHANS/Down_stream_analysis/NIMHAMSg/EWAS/Plots/QQplot_of_EWAS_Age.pdf")
qq(res1$Age_P)
dev.off()


##Manhattan Plot
res2<-res2[which(res2$CHR != "Y"),] ## to also exclude the X chromosome repeat and edit this line

res2$CHR<-as.character(res2$CHR)
res2$CHR[which(res2$CHR == "X")]<-23 ## to also recode Y chromosome repeat and edit this line
res2$CHR<-as.numeric(res2$CHR)
res2<-res2[which(res2$CHR != ""),]

bonfP<-0.05/nrow(res2)

pdf("/mnt/data1/NIMHANS/Down_stream_analysis/NIMHAMSg/EWAS/Plots/Manhattan_Plots_of_Age.pdf")
manhattan(res2, p = "Status_P", bp = "MAPINFO", chr = "CHR", genomewide = -log10(bonfP), suggestiveline = -log10(5e-5), logp=T, col=c("blue","yellow"), ylim=c(0,30))
dev.off()

