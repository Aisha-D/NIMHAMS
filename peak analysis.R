###Sensitivity of probes type1,type2 

#Load in mset and beta data
setwd('/mnt/data1/PRISM/')
load('/mnt/data1/PRISM/PRISMmsetEPIC.rdat')

betas <- betas(msetEPIC2)
onetwo<-fData(msetEPIC2)$DESIGN
mat<-betas(msetEPIC2)

#visualise the probes
plotmset_density<-function(mset, study=""){
  onetwo<-fData(mset)$DESIGN
  mat<-betas(mset)
  
  plot(density(mat[onetwo=="I",1], na.rm=T, bw=0.03), cex.main=0.8, main=paste(study, "Betas"), ylim=c(0, 5.2), xlab="")
  lines(density(mat[onetwo=="II",1], na.rm=T, bw=0.03), col="red")
  
  for(j in 2:ncol(mat)){
    lines(density(mat[onetwo=="I",j], na.rm=T, bw=0.03))
    lines(density(mat[onetwo=="II",j], na.rm=T, bw=0.03), col="red")
  }
  
  legend("topright", legend=c("Type I", "Type II"), lty=1, col=c("black", "red")) 
}
plotmset_density(msetEPIC2, study = 'PRISM')

pk <- matrix(data = NA, nrow = ncol(mat), ncol = 4)
for(j in 1:ncol(mat)){
  one <- density(mat[onetwo=="I",j], na.rm=T, bw=0.03)
  one <- one[1:2]
  one <- as.data.frame(one)
  oneunmeth <- one[which(one$x < 0.5),]
  onemeth <- one[which(one$x > 0.5),]
  one_unme <- oneunmeth$y
  one_meth <- onemeth$y
  pk[j,1] <- max(one_unme)
  pk[j,2] <- max(one_meth)
  two <- density(mat[onetwo=="II",j], na.rm=T, bw=0.03)
  two <- two[1:2]
  two <- as.data.frame(two)
  twounmeth <- two[which(two$x < 0.5),]
  twometh <- two[which(two$x > 0.5),]
  two_unme <- twounmeth$y
  two_meth <- twometh$y
  pk[j,3] <- max(two_unme)
  pk[j,4] <- max(two_meth)
}

#1=type 1 unmeth, 2=type 1 meth, 3=type 2 unmeth, 4=type 2 meth
boxplot(pk[,1], pk[,2], pk[,3], pk[,4] )
sdpk1 <- sd(pk[,1])
sdpk2 <- sd(pk[,2])
sdpk3 <- sd(pk[,3])
sdpk4 <- sd(pk[,4])



#normalised data
library(wateRmelon)
msetEPIC.dasen<-dasen(msetEPIC2)
onetwo<-fData(msetEPIC.dasen)$DESIGN
mat<-betas(msetEPIC.dasen)
pk2 <- matrix(data = NA, nrow = ncol(mat), ncol = 4)
for(j in 1:ncol(mat)){
  one <- density(mat[onetwo=="I",j], na.rm=T, bw=0.03)
  one <- one[1:2]
  one <- as.data.frame(one)
  oneunmeth <- one[which(one$x < 0.5),]
  onemeth <- one[which(one$x > 0.5),]
  one_unme <- oneunmeth$y
  one_meth <- onemeth$y
  pk2[j,1] <- max(one_unme)
  pk2[j,2] <- max(one_meth)
  two <- density(mat[onetwo=="II",j], na.rm=T, bw=0.03)
  two <- two[1:2]
  two <- as.data.frame(two)
  twounmeth <- two[which(two$x < 0.5),]
  twometh <- two[which(two$x > 0.5),]
  two_unme <- twounmeth$y
  two_meth <- twometh$y
  pk2[j,3] <- max(two_unme)
  pk2[j,4] <- max(two_meth)
}

boxplot(pk2[,1], pk2[,2], pk2[,3], pk2[,4] )
sdpk2_1 <- sd(pk2[,1])
sdpk2_2 <- sd(pk2[,2])
sdpk2_3 <- sd(pk2[,3])
sdpk2_4 <- sd(pk2[,4])


mnpk <- mean(pk[,1])



