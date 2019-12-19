library(data.table)

tab<-matrix(data = NA, nrow = 22, ncol = 8)
colnames(tab)<-c("nmQTLs", "nmQTLs_PCs", "overlap", "mean effect mQTLs", 
                 "mean effect mQTLs with PCs", "mean effect overlap mQTLs", 
                 "mean effect mQTLs not with PCs", "mean effect mQTLs only with PCs")

for(chr in 1:22){
	output_file_name = paste("mQTL/Output/without_cov_PCs/Ind_chr", chr, ".txt", sep = "")
	tmp<-fread(output_file_name)
	output_file_name = paste("mQTL/Output/Ind_chr", chr, ".txt", sep = "")
	tmp.2<-fread(output_file_name)
	
	tab[chr,1]<-nrow(tmp)
	tab[chr,2]<-nrow(tmp.2)
	overlap<-table(paste(tmp$SNP, tmp$gene) %in% paste(tmp.2$SNP, tmp.2$gene))
	tab[chr,3]<-overlap[1]
	tab[chr,4]<-mean(abs(tmp$beta))
	tab[chr,5]<-mean(abs(tmp.2$beta))
	tab[chr,6]<-mean(abs(tmp$beta[paste(tmp$SNP, tmp$gene) %in% paste(tmp.2$SNP, tmp.2$gene)]))
	tab[chr,7]<-mean(abs(tmp$beta[!paste(tmp$SNP, tmp$gene) %in% paste(tmp.2$SNP, tmp.2$gene)]))
	tab[chr,8]<-mean(abs(tmp.2$beta[!paste(tmp.2$SNP, tmp.2$gene) %in% paste(tmp$SNP, tmp$gene)]))
}
	
write.table(tab, file = "/gpfs/ts0/scratch/and202/NIMHAMS/mQTL/Output/Cerebellum_mQTL.txt")