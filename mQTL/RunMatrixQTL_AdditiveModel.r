args <- commandArgs(trailingOnly = TRUE)
chr<-args

setwd("/gpfs/ts0/scratch/and202/NIMHAMS/")
library(MatrixEQTL)
base.dir = find.package("MatrixEQTL")
useModel = modelLINEAR; # modelANOVA or modelLINEAR or modelLINEAR_CROSS 

covariates_file_name = "/gpfs/ts0/scratch/and202/NIMHAMS/mQTL/cov_agesex.txt";

##load covariate data
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
cvrt$LoadFile(covariates_file_name);

SNP_file_name = paste("mQTL/Ind_Genotypes_Imputed_MinGenoCount5_Chr", chr, ".txt", sep = "")
#snps_location_file_name = 
	
## threshold of results to save
pvOutputThreshold = 1e-8; #save data less than this threshold
cisDist<-500000
## set error covarriance. NOTE rarely used instead set to numeric()
errorCovariance = numeric();

 ## load SNP data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA" ; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile(SNP_file_name);


expression_file_name<-"/gpfs/ts0/scratch/and202/NIMHAMS/Ind_Methylation.txt"

#gene_location_file_name = paste(base.dir, "/data/geneloc.txt", sep="");
## if no covariates set to  character()
output_file_name = paste("mQTL/Output/Ind_chr", chr, ".txt", sep = "")

## load methylation data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
gene$LoadFile(expression_file_name);

## Run the analysis
#snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
#snpspos[,1]<-rownames(snps)

### run eqtls
me = Matrix_eQTL_main(
	snps = snps,
	gene = gene,
	cvrt = cvrt,
	output_file_name = output_file_name,
	pvOutputThreshold = pvOutputThreshold,
	useModel = useModel, 
	errorCovariance = errorCovariance, 
	verbose = TRUE,
	pvalue.hist = "qqplot",
	min.pv.by.genesnp = FALSE,
	noFDRsaveMemory = FALSE)

