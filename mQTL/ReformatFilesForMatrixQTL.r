setwd("/gpfs/ts0/scratch/and202/NIMHAMS/")

load("NIMHAM_CER_Normalised_snpsremoved.rdat")
betasc <- betas

pheno <- phenoc
betas <-  betasc

crosshyb<-read.table("/gpfs/ts0/scratch/and202/EPIC_reference/CrossHydridisingProbes_McCartney.txt", stringsAsFactors = FALSE)
tofilter<-read.csv("/gpfs/ts0/scratch/and202/EPIC_reference/EPICArrayProbesToFilter.csv", stringsAsFactors = FALSE)
snpProbes<-read.table("/gpfs/ts0/scratch/and202/EPIC_reference/SNPProbes_McCartney.txt", stringsAsFactors = FALSE, header = TRUE)
crosshyb2<-read.csv("/gpfs/ts0/scratch/and202/EPIC_reference/Pidsley_SM1.csv", stringsAsFactors = FALSE)
snpProbes2<-read.csv("/gpfs/ts0/scratch/and202/EPIC_reference/Pidsley_SM4.csv", stringsAsFactors = FALSE)
snpProbes3<-read.csv("/gpfs/ts0/scratch/and202/EPIC_reference/Pidsley_SM5.csv", stringsAsFactors = FALSE)
snpProbes4<-read.csv("/gpfs/ts0/scratch/and202/EPIC_reference/Pidsley_SM6.csv", stringsAsFactors = FALSE)

snpProbes<-snpProbes[which(snpProbes$DIST_FROM_MAPINFO < 10 & snpProbes$SAS_AF > 0.01),]
snpProbes2<-snpProbes2[which(snpProbes2$SAS_AF > 0.01),]
snpProbes3<-snpProbes3[which(snpProbes3$SAS_AF > 0.01),]
snpProbes4<-snpProbes4[which(snpProbes4$SAS_AF > 0.01),]

dist<-cbind(abs(snpProbes4$VARIANT_END - snpProbes4$MAPINFO), abs(snpProbes4$VARIANT_START - snpProbes4$MAPINFO))
dist<-apply(dist, 1, min)
snpProbes4<-snpProbes4[which(dist <=10),]

remove<-unique(c(tofilter$IlmnID, crosshyb[,1], snpProbes$IlmnID, snpProbes2$PROBE, snpProbes3$PROBE, snpProbes4$PROBE))

dat <- betas            
dat<-dat[-match(intersect(remove, rownames(dat)), rownames(dat)),]
#dat<-dat[-grep("rs", rownames(dat)),] #no rs probes

epicManifest<-read.csv("/gpfs/ts0/scratch/and202/EPIC_reference/MethylationEPIC_v-1-0_B2.csv", skip = 7)
epicManifest<-epicManifest[match(rownames(dat), epicManifest$IlmnID),]

write.table(cbind(paste(pheno$Sample_Name, pheno$Sample_Name, sep = "_"),
                  paste(pheno$Sample_Name, pheno$Sample_Name, sep = "_")), 
            "DNAmSamples.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

chr<-22
geno<-read.table(paste("mQTL/Ind_maf_geno_mind_hwe_chr", chr, ".raw", sep = ""), header = T, stringsAsFactors = FALSE)
geno.map<-read.table(paste("mQTL/Ind_maf_geno_mind_hwe_chr", chr, ".bim", sep = ""), stringsAsFactors = FALSE)
colnames(geno.map)<-c("Chr", "id", "cm", "bp", "A1", "A2")

### match methylation & genotypes files
##The genotype file contains 48 samples which means the failed qc samples have not been used
#The pheno file includes passed samples only - we need to remove the failsed samples from geno
pheno$SampleType <- as.character(pheno$SampleType)
pheno$Sample_Name <- as.character(pheno$Sample_Name)

geno2 <- pheno$Sample_Name[!duplicated(geno$IID)] #remove samples from geno that did not pass QC
geno_2 <- geno[which(geno$IID %in% geno2),]
geno3 <- geno_2$IID[!duplicated(geno2)] #remove the sample 7c which has no genotype data but has methylation
pheno_2 <- pheno[which(pheno$Sample_Name %in% geno3),]

geno <- geno_2
pheno <- pheno_2
identical(pheno$Sample_Name, geno$IID) #check that they have the same samples

#the beta/dat methyaltion data needs to match the pheno data
dat <- dat[,which(colnames(dat) %in% rownames(pheno))]

write.table(dat, "Ind_Methylation.txt", sep = "\t", quote = FALSE)
write.table(epicManifest[,c("IlmnID", "CHR", "MAPINFO")], "Ind_Methylation_Info.txt", sep = "\t", quote = FALSE)

#PCA??
pcs<-read.table("/gpfs/ts0/scratch/and202/NIMHAMS/nimhams2.pca.eigenvec")
pcs$V2 <- as.character(pcs$V2)
#there are some outliers so remove these - most likelt the samples that are not in pheno
#without remocing the samples it looks like 48 and 49 samples are outliers by ~1%
pcs<- pcs[which(pcs$V2 %in% geno3),]
pcs$V2 <- as.factor(pcs$V2)
pdf("sample pca2.pdf")
plot(pcs[,3], pcs[,4], xlab = "PC1", ylab = "PC2",
     col = rainbow(nlevels(pcs$V2))[pcs$V2]) #third and fourth values are PC1 and Pc2 in matrix
legend("topright", levels(pcs$V2), col = rainbow(nlevels(pcs$V2)), pch = 10, cex=0.5)
dev.off()

#Make sex numeric
pheno$Sex <- as.character(pheno$Sex)
pheno$nsex  <- rep(NA, nrow(pheno))
pheno$nsex[pheno$Sex == 'F'] <- 2
pheno$nsex[pheno$Sex == 'M'] <- 1

#make a list of covariates
cov<-t(cbind((as.numeric(pheno$nsex)-1),pheno$Age, pcs[,3:12]))
colnames(cov)<-pheno$Sample_Name
rownames(cov)<-c("Sex", "Age", paste("PC", 1:10, sep = ""))
write.table(cov, "mQTL/cov_with10PCs.txt", sep = "\t", quote = FALSE)

cov<-t(cbind((as.numeric(pheno$nsex)-1),pheno$Age))
colnames(cov)<-pheno$Sample_Name
rownames(cov)<-c("Sex", "Age")
write.table(cov, "mQTL/cov_agesex.txt", sep = "\t", quote = FALSE)

geno.ids<-geno$FID
geno<-geno[,-c(1:6)]
geno<-t(geno)
colnames(geno)<-pheno$Sample_Name

#remove the rows/snps that appear minimum less than 5 times
a<-NULL
for(i in 1:nrow(geno)){
	if(min(table(geno[i,])) < 5){
		a<-append(a, i)
		}
}
geno<-geno[-a,]
geno.map<-geno.map[-a,]
write.table(geno, "mQTL/Ind_Genotypes_Imputed_MinGenoCount5_Chr22.txt", sep = "\t", quote = FALSE)
write.table(cbind(geno.map$id, geno.map$Chr, geno.map$bp), "mQTL/Ind_Genotypes_Imputed_MapInfo_MinGenoCount5_Chr22.txt", sep = "\t", quote = FALSE)

for(chr in 1:21){

	geno<-read.table(paste("mQTL/Ind_maf_geno_mind_hwe_chr", chr, ".raw", sep = ""), header = T, stringsAsFactors = FALSE)
	geno.map<-read.table(paste("mQTL/Ind_maf_geno_mind_hwe_chr", chr, ".bim", sep = ""), stringsAsFactors = FALSE)
	colnames(geno.map)<-c("Chr", "id", "cm", "bp", "A1", "A2")
	
	#match geno and pheno data
	geno2 <- pheno$Sample_Name[!duplicated(geno$IID)] #remove samples from geno that did not pass QC
	geno_2 <- geno[which(geno$IID %in% geno2),]
	geno3 <- geno_2$IID[!duplicated(geno2)] #remove the sample 7c which has no genotype data but has methylation
	pheno_2 <- pheno[which(pheno$Sample_Name %in% geno3),]
	
	geno <- geno_2
	pheno <- pheno_2
	identical(pheno$Sample_Name, geno$IID) #check that they have the same samples
	
	geno<-geno[,-c(1:6)]
	geno<-t(geno)
	colnames(geno)<-pheno$Sample_Name


	a<-NULL
	for(i in 1:nrow(geno)){
		if(min(table(geno[i,])) < 5){
			a<-append(a, i)
			}
	}
	geno<-geno[-a,]
	geno.map<-geno.map[-a,]
	write.table(geno, paste("mQTL/Ind_Genotypes_Imputed_MinGenoCount5_Chr", chr, ".txt", sep = ""), sep = "\t", quote = FALSE)
	write.table(cbind(geno.map$id, geno.map$Chr, geno.map$bp), paste("mQTL/Ind_Genotypes_Imputed_MapInfo_MinGenoCount5_Chr", chr, ".txt", sep = ""), sep = "\t", quote = FALSE)
	
}