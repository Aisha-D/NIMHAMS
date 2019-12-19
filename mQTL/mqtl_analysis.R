#### READ IN ALL CHR FILES ####
library(ggplot2)
library(dplyr)
library(tidyr)

###### Cerebellum Aisha  #####
setwd('/gpfs/ts0/scratch/and202/NIMHAMS/mQTL/Output/')

cereb_mqtl <- read.table("/gpfs/ts0/scratch/and202/NIMHAMS/mQTL/Output/Cerebellum_mQTL.txt")

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


##Add CPG annotation data
epicManifest <- read.csv("/gpfs/ts0/scratch/and202/NIMHAMS/epicManifest_hg38.csv", 
                         stringsAsFactors = F, header = T, sep =" ")

genotypeManifest <- read.csv("/gpfs/ts0/scratch/gn261/Nimhams/Infinium_Global_Screening-24_v2.0_manifest.csv")



meth.all <- cbind(all_chr, epicManifest[match(all_chr$gene, epicManifest$IlmnID),
                                        c("CHR", "MAPINFO")])

#Recode X and Y chromosomes and remove snps that do not have chromosome locations
meth.all <- meth.all[which(meth.all$CHR != "X"),]
meth.all <- meth.all[which(meth.all$CHR != "Y"),]
meth.all<-meth.all[which(!is.na(meth.all$CHR)),] #5416 to 5416
meth.all$CHR<-as.numeric(as.character(meth.all$CHR))

meth.all <- meth.all[, c("SNP", "CHR", "MAPINFO")]

##Add SNP annotation
all_chr$SNP <- gsub("_A", "", all_chr$SNP)
all_chrSNP <- gsub("_T", "", all_chr$SNP)
all_chr$SNP <- gsub("_G", "", all_chr$SNP)
all_chr$SNP <- gsub("_C", "", all_chr$SNP)

snps.all <- cbind(all_chr, genotypeManifest[match(all_chr$SNP, genotypeManifest$Name), 
                                            c("Chr", "MapInfo")])

snps.all <-snps.all[which(!is.na(snps.all$Chr)),] #8071 to 5632
snps.all$Chr<-as.numeric(as.character(snps.all$Chr))

snps.all <- snps.all[, c("SNP", "Chr", "MapInfo")]

chrEnd.meth<-vector(length = 22)
chrEnd.snps<-vector(length = 22)
for(i in 1:22){
  chrEnd.meth[i]<-max(meth.all[which(meth.all$CHR == i),3])
  chrEnd.snps[i]<-max(snps.all[which(snps.all[,2] == i),3])
}

chrStart.meth<-vector(length = 22)
chrStart.snps<-vector(length = 22)
for(i in 1:22){
  chrStart.meth[i]<-min(meth.all[which(meth.all$CHR == i),3])
  chrStart.snps[i]<-min(snps.all[which(snps.all[,2] == i),3])
}

chrSize.snps<-chrEnd.snps-chrStart.snps
chrSize.meth<-chrEnd.meth-chrStart.meth

snp.coord<-vector(length = nrow(snps.all))
meth.coord<-vector(length = nrow(meth.all))
for(i in 1:22){
  snp.coord[which(snps.all[,2] == i)]<-(snps.all[which(snps.all[,2] == i),3]-chrStart.snps[i])/chrSize.snps[i]+i-1
  meth.coord[which(meth.all[,2] == i)]<-(meth.all[which(meth.all[,2] == i),3]-chrStart.meth[i])/chrSize.meth[i]+i-1
}

snps.all<-cbind(snps.all, snp.coord)
meth.all<-cbind(meth.all, meth.coord)

points.y<-snps.all[match(gsub("_.", "", all_chr$SNP), snps.all[,1]),4]
points.x<-meth.all[match(all_chr$SNP, meth.all$SNP),4]

### take absolute values of effect size
logP<-abs(all_chr$beta)

### see source data for resulting data points.

logP_col_palette<-cbind(seq(from = 0, to = 0.61, by = 0.01) , colorRampPalette(c("tan1", "orange",  "red", "red3", "red4",  "darkred"))(62))

pdf("MQTL_cerebellum.pdf")
layout(matrix(c(1,2), ncol = 1), height = c(1,0.1))
op <- par(oma=c(5,7,1,1))
par(op)

points.col<-logP_col_palette[match(round(logP,2), logP_col_palette[,1]),2]
plot(points.x, points.y, pch = 15, cex =0.75, col = points.col, type = "n", ylab = "Position of SNP (Chr)", xlab = "Position of DNA methylation site (Chr)", main = "",  axes = FALSE, xlim = c(0,22), ylim = c(0,22), xaxs = "i", yaxs = "i",cex.lab = 1, cex.axis = 1)
axis(1, at = c(-1,seq(0.5,21.5, by = 1), 23), c("", 1:22, ""), cex.axis = 1)
axis(2, at = c(-1,seq(0.5,21.5, by = 1), 23), c("", 1:22, ""), cex.axis = 1, las = 2)
points(points.x[length(points.x):1], points.y[length(points.x):1], pch = 15, cex =0.75, col = points.col[length(points.x):1])

### plot legend
par(op)
plot(as.numeric(logP_col_palette[,1]), rep(1, nrow(logP_col_palette)), col = logP_col_palette[,2], pch = 15, cex = 1, axes = FALSE, main  = "", ylab = "", ylim = c(0.9999,1.001),xlab = "", cex.lab = 1, cex.axis = 1, cex.main = 1, adj = 0)
axis(1, at = c(0,0.2,0.4,0.6),  c(0,0.2,0.4,0.6)*100, cex.axis = 1)
dev.off()


############# ggplot ################
pdf("MQTL_ggplot.pdf")
ggplot(allchr3, aes(CHR, Chr)) +
  geom_jitter(aes(color = allchr3$beta)) +
  scale_color_gradient(low = "orange", high = "red")
dev.off()



########### Cis ###############
### output is table of all bonferonni significant mQTLs
### see associated data source for relevant columns

## calculate distance between SNP and methylation probe
all_chr <-  cbind(all_chr, epicManifest[match(all_chr$gene, epicManifest$IlmnID),
                                        c("CHR", "MAPINFO")])
all_chr <-  cbind(all_chr, genotypeManifest[match(all_chr$SNP, genotypeManifest$Name), 
                                            c("Chr", "MapInfo")])

dist<-abs(all_chr$MapInfo - all_chr$MAPINFO)
dist[which(as.character(all_chr$Chr) != as.character(all_chr$CHR))]<--9

signedDist<-(all_chr$MapInfo-all_chr$MAPINFO)/1000

par(op)
pdf("Cis_mQTL_Cerebellum.pdf")
plot(signedDist[which(dist != -9)], abs(all_chr$beta[which(dist != -9)])*100, main = "", 
     xlab = "Distance (Mb)", ylab = "Effect size (% difference in DNA methylation per minor allele)",
     pch = 16, cex = 0.8, ylim = c(-0,45), xlim = c(-500,500), cex.axis = 1, cex.lab = 0.8, axes = FALSE, xaxs = "i", 
     yaxs = "i")
axis(1, seq(-1, 1, 0.5), at = seq(-1000, 1000, 500), cex.axis = 1, cex.lab = 1)
axis(2, las = 2, cex.axis = 1, cex.lab = 1)
dev.off()


########### Trans ##############
library(RCircos)
#all_chr.og <- all_chr
all_chr2 <- na.omit(all_chr, x = Chr)
all_chr2[which(all_chr2$Chr == 0)] <- NA
all_chr2 <- na.omit(all_chr2, x = Chr)
all_chr2<-all_chr2[which(all_chr2$CHR != "X"),] 
all_chr2<-all_chr2[which(all_chr2$CHR != "Y"),] 
all_chr2<-all_chr2[which(all_chr2$Chr != "X"),] 
all_chr2<-all_chr2[which(all_chr2$Chr != "Y"),] 
all_chr2$Chr <- as.character(all_chr2$Chr)
all_chr <- all_chr2

### res is table of all mQTLs
### identify transcisTrans<-rep("cis", nrow(res))
cisTrans<-rep("cis", nrow(all_chr))
cisTrans[which(all_chr$Chr != all_chr$CHR | abs(all_chr$MapInfo - all_chr$MAPINFO) > 500000)]<-"trans"


### set up plot
data(UCSC.HG38.Human.CytoBandIdeogram)


chr.exclude <- c("chrX", "chrY")
cyto.info <- UCSC.HG38.Human.CytoBandIdeogram;
tracks.inside <- 10;
tracks.outside <- 0;
RCircos.Set.Core.Components(cyto.info, chr.exclude,tracks.inside, tracks.outside);

rcircos.params <- RCircos.Get.Plot.Parameters();
rcircos.cyto <- RCircos.Get.Plot.Ideogram();
rcircos.position <- RCircos.Get.Plot.Positions();

out.file <- "RCircosDemoHumanGenome.pdf";
pdf(file=out.file, height=8, width=8, compress=TRUE);
RCircos.Set.Plot.Area();

par(mai=c(0.25, 0.25, 0.25, 0.25));
plot.new();
plot.window(c(-2.5,2.5), c(-2.5, 2.5));
RCircos.Chromosome.Ideogram.Plot();

track.num <- 11;
data(RCircos.Link.Data)
RCircos.Link.Plot(RCircos.Link.Data, 0, FALSE); #Error here!!
data(RCircos.Ribbon.Data);
RCircos.Ribbon.Plot(ribbon.data=RCircos.Ribbon.Data, track.num=11, by.chromosome=FALSE, twist=FALSE);


tmp<-all_chr[which(cisTrans == "trans"), c(7,8,8, 9,10,10)]
tmp[,1]<-paste("chr", tmp[,1], sep = "")
tmp[,4]<-paste("chr", tmp[,4], sep = "")
tmp <- tmp[which(tmp$Chr != 'chr0'),]


params <- RCircos.Get.Plot.Parameters();
params$text.size <- 2;
RCircos.Reset.Plot.Parameters(params)


pdf(file=out.file, height=12, width=12);
RCircos.Set.Plot.Area();

#par(mai=c(0.25, 0.25, 0.25, 0.25));
#plot.new();
#plot.window(c(-2.5,2.5), c(-2.5, 2.5));
RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
RCircos.Pos <- RCircos.Get.Plot.Positions()
RCircos.Par <- RCircos.Get.Plot.Parameters()
right.side <- nrow(RCircos.Pos)/2
outer.location <- RCircos.Par$chr.ideog.pos + RCircos.Par$chrom.width
inner.location <- RCircos.Par$chr.ideog.pos
chroms <- unique(RCircos.Cyto$Chromosome)
for (a.chr in 1:length(chroms)) {
  the.chr <- RCircos.Cyto[RCircos.Cyto$Chromosome == chroms[a.chr],
                          ]
  start <- the.chr$StartPoint[1] 
  end <- the.chr$EndPoint[nrow(the.chr)]
  mid <- round((end - start)/2, digits = 0) + start
  chr.color <- the.chr$ChrColor[nrow(the.chr)]
  pos.x <- c(RCircos.Pos[start:end, 1] * outer.location,
             RCircos.Pos[end:start, 1] * inner.location)
  pos.y <- c(RCircos.Pos[start:end, 2] * outer.location,
             RCircos.Pos[end:start, 2] * inner.location)
  polygon(pos.x, pos.y)
  chr.name <- sub(pattern = "chr", replacement = "", chroms[a.chr])
  text(RCircos.Pos[mid, 1] * RCircos.Par$chr.name.pos,
       RCircos.Pos[mid, 2] * RCircos.Par$chr.name.pos, label = chr.name,
       srt = RCircos.Pos$degree[mid], cex = 2)
  lines(RCircos.Pos[start:end, ] * RCircos.Par$highlight.pos,
        col = chr.color, lwd = RCircos.Par$highlight.width)
}
for (a.band in 1:nrow(RCircos.Cyto)) {
  a.color <- RCircos.Cyto$BandColor[a.band]
  if (a.color == "white") {
    next
  }
  start <- RCircos.Cyto$ChromStart[a.band] 
  end <- RCircos.Cyto$ChromEnd[a.band]
  pos.x <- c(RCircos.Pos[start:end, 1] * outer.location,
             RCircos.Pos[end:start, 1] * inner.location)
  pos.y <- c(RCircos.Pos[start:end, 2] * outer.location,
             RCircos.Pos[end:start, 2] * inner.location)
  polygon(pos.x, pos.y, col = a.color, border = NA)
}


### get colours of chromosome
RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
ChrColMap<-unique(RCircos.Cyto[,c(1,7)])
LineCexMap<-cbind(seq(from = 0, to = 0.4, by = 0.01), seq(from = 0.5, to = 4, by = 0.0875))
dev.off()

#RCircos.Link.Plot(tmp, 1, FALSE);
RCircos.Pos <- RCircos.Get.Plot.Positions()
RCircos.Par <- RCircos.Get.Plot.Parameters()
tmp2 <- RCircos.Validate.Genomic.Data(tmp, plot.type = "link") ###ISSUE IS HERE!!!
one.track <- RCircos.Par$track.height + RCircos.Par$track.padding
start <- RCircos.Par$track.in.start - (1 - 1) * one.track
base.positions <- RCircos.Pos * start
data.points <- matrix(rep(0, nrow(tmp) * 2), ncol = 2)
for (a.link in 1:nrow(tmp)) {
  data.points[a.link, 1] <- RCircos.Data.Point(tmp[a.link,
                                                   1], tmp[a.link, 2])
  data.points[a.link, 2] <- RCircos.Data.Point(tmp[a.link,
                                                   4], tmp[a.link, 5])
  if (data.points[a.link, 1] == 0 || data.points[a.link,
                                                 2] == 0) {
    print("Error in chromosome locations ...")
    break
  }
}
#   link.colors <- RCircos.Get.Link.Colors(tmp, by.chromosome)
link.colors<-ChrColMap[match(tmp[,4], ChrColMap[,1]),2]
link.cex.line<-LineCexMap[match(round(abs(res$Fetal.beta[which(cisTrans == "trans")]), 2), LineCexMap[,1]),2]

for (a.link in 1:nrow(data.points)) {
  point.one <- data.points[a.link, 1]
  point.two <- data.points[a.link, 2]
  if (point.one > point.two) {
    point.one <- data.points[a.link, 2]
    point.two <- data.points[a.link, 1]
  }
  P0 <- as.numeric(base.positions[point.one, ])
  P2 <- as.numeric(base.positions[point.two, ])
  links <- RCircos.Link.Line(P0, P2)
  lines(links$pos.x, links$pos.y, type = "l", col = link.colors[a.link], lwd = link.cex.line[a.link])
}






