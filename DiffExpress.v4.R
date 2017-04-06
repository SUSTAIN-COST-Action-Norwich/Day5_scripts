# Upload previous results from OTUsummary.vX.R
setwd("/home/student1/Day5/02_ORFsummary/")
plantORFsF <- read.table("plantORFsF.txt", header=T, sep="\t")
pathogenORFsF <- read.table("pathogenORFsF.txt", header=T, sep="\t")
uniORFsummary <- read.table("uniORFsummary.txt", header=T, sep="\t")

############################################################
# Upload Cufflinks results calculated on the whole dataset #
############################################################
setwd("/home/student1/Day5/01_Data/")
h21RPKM <- read.table("Hpc21.cufflinks.tab", sep="\t", header=T)
summary(h21RPKM)
h31RPKM <- read.table("Hpc31.cufflinks.tab", sep="\t", header=T)
h42RPKM <- read.table("Hpc42.cufflinks.tab", sep="\t", header=T)
R13RPKM <- read.table("Reseda13.cufflinks.tab", sep="\t", header=T)
R21RPKM <- read.table("Reseda21.cufflinks.tab", sep="\t", header=T)
R31RPKM <- read.table("Reseda31.cufflinks.tab", sep="\t", header=T)


##########################################################
# Dotplot RPKM values of plant ORFs for each sample pair #
##########################################################
setwd("/home/student1/Day5/05_DiffExpress/")
pdf("RPKMcomparison.pdf", width=8, height=16)
old.par <- par(mfrow=c(5,3))
plot(((subset(h21RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, ((subset(h31RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, pch=".", log="xy", xlim=c(0.01,50000), ylim=c(0.01,50000), xlab="Hpc#2-1", ylab="Hpc#3-1", main="H31_H21")
plot(((subset(h21RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, ((subset(h42RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, pch=".", log="xy", xlim=c(0.01,50000), ylim=c(0.01,50000), xlab="Hpc#2-1", ylab="Hpc#4-2", main="H42_H21")
plot(((subset(h31RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, ((subset(h42RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, pch=".", log="xy", xlim=c(0.01,50000), ylim=c(0.01,50000), xlab="Hpc#3-1", ylab="Hpc#4-2", main="H42_H31")
plot(((subset(h21RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, ((subset(R13RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, pch=".", log="xy", xlim=c(0.01,50000), ylim=c(0.01,50000), xlab="Hpc#2-1", ylab="Reseda#1-3", main="R13_H21")
plot(((subset(h21RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, ((subset(R21RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, pch=".", log="xy", xlim=c(0.01,50000), ylim=c(0.01,50000), xlab="Hpc#2-1", ylab="Reseda#2-1", main="R21_H21")
plot(((subset(h21RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, ((subset(R31RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, pch=".", log="xy", xlim=c(0.01,50000), ylim=c(0.01,50000), xlab="Hpc#2-1", ylab="Reseda#3-1", main="R31_H21")
plot(((subset(h31RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, ((subset(R13RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, pch=".", log="xy", xlim=c(0.01,50000), ylim=c(0.01,50000), xlab="Hpc#3-1", ylab="Reseda#1-3", main="R13_H31")
plot(((subset(h31RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, ((subset(R21RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, pch=".", log="xy", xlim=c(0.01,50000), ylim=c(0.01,50000), xlab="Hpc#3-1", ylab="Reseda#2-1", main="R21_H31")
plot(((subset(h31RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, ((subset(R31RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, pch=".", log="xy", xlim=c(0.01,50000), ylim=c(0.01,50000), xlab="Hpc#3-1", ylab="Reseda#3-1", main="R31_H31")
plot(((subset(h42RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, ((subset(R13RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, pch=".", log="xy", xlim=c(0.01,50000), ylim=c(0.01,50000), xlab="Hpc#4-2", ylab="Reseda#1-3", main="R13_H42")
plot(((subset(h42RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, ((subset(R21RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, pch=".", log="xy", xlim=c(0.01,50000), ylim=c(0.01,50000), xlab="Hpc#4-2", ylab="Reseda#2-1", main="R21_H42")
plot(((subset(h42RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, ((subset(R31RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, pch=".", log="xy", xlim=c(0.01,50000), ylim=c(0.01,50000), xlab="Hpc#4-2", ylab="Reseda#3-1", main="R31_H42")
plot(((subset(R13RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, ((subset(R21RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, pch=".", log="xy", xlim=c(0.01,50000), ylim=c(0.01,50000), xlab="Reseda#1-3", ylab="Reseda#2-1", main="R21_R13")
plot(((subset(R13RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, ((subset(R31RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, pch=".", log="xy", xlim=c(0.01,50000), ylim=c(0.01,50000), xlab="Reseda#1-3", ylab="Reseda#3-1", main="R31_R13")
plot(((subset(R21RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, ((subset(R31RPKM, gene_id %in% plantORFsF$orfID))$FPKM)+0.01, pch=".", log="xy", xlim=c(0.01,50000), ylim=c(0.01,50000), xlab="Reseda#2-1", ylab="Reseda#3-1", main="R31_R21")
par(old.par)
dev.off()

##################################################################################
# Merge RPKM results to calculate mean RPKM of infected and non-infected samples #
##################################################################################
RPKMmerge <- as.data.frame(matrix(0,nrow=dim(h21RPKM)[1], ncol=7))
colnames(RPKMmerge) <- c("orfID", "h21FPKM", "h31FPKM", "h42FPKM", "R13FPKM", "R21FPKM", "R31FPKM")
RPKMmerge$orfID <- h21RPKM$gene_id
RPKMmerge$h21FPKM <- h21RPKM$FPKM
for(i in 1:dim(RPKMmerge)[1]){
  s0 <- (subset(h21RPKM, gene_id %in% RPKMmerge[i,"orfID"]))
  s1 <- (subset(h31RPKM, gene_id %in% RPKMmerge[i,"orfID"]))
  RPKMmerge[i,"h31FPKM"] <- s1[1,"FPKM"]
  s2 <- (subset(h42RPKM, gene_id %in% RPKMmerge[i,"orfID"]))
  RPKMmerge[i,"h42FPKM"] <- s2[1,"FPKM"]
  s3 <- (subset(R13RPKM, gene_id %in% RPKMmerge[i,"orfID"]))
  RPKMmerge[i,"R13FPKM"] <- s3[1,"FPKM"]
  s4 <- (subset(R21RPKM, gene_id %in% RPKMmerge[i,"orfID"]))
  RPKMmerge[i,"R21FPKM"] <- s4[1,"FPKM"]
  s5 <- (subset(R31RPKM, gene_id %in% RPKMmerge[i,"orfID"]))
  RPKMmerge[i,"R31FPKM"] <- s5[1,"FPKM"]
  h1 <- c(s0$FPKM, s1$FPKM,s2$FPKM)
  r1 <- c(s3$FPKM, s4$FPKM,s5$FPKM)
  RPKMmerge[i,"meanHpc.FPKM"] <- mean(h1)
  RPKMmerge[i,"stdvHpc.FPKM"] <- sd(h1)
  RPKMmerge[i,"meanReseda.FPKM"] <- mean(r1)
  RPKMmerge[i,"stdvReseda.FPKM"] <- sd(r1)
  rm(s0, s1, s2, s3, s4, s5, h1, r1)
}
#RPKMmerge <- subset(RPKMmerge, orfID %in% uniORFsummary$orfID)

##########################################################################
# Plot mean RPKM values of infected samples against non-infected samples #
##########################################################################
AtORFs <- as.character((subset(uniORFsummary, HpaHomolog == 0 & AtHomolog == 1))$orfID)  # unique ORFs with homologs only in A. thaliana
HpaORFs <- as.character((subset(uniORFsummary, HpaHomolog == 1 & AtHomolog == 0))$orfID) # unique ORFs with homologs only in H. arabidopsidis
AtHpaORFs <- as.character((subset(uniORFsummary, HpaHomolog == 1 & AtHomolog == 1))$orfID) # unique ORFs with homologs in both, A. thaliana and H. arabidopsidis
pdf("meanRPKMcomparison.Arabidopsis.Hpa.pdf", width=8, height=8)
plot(((subset(RPKMmerge, orfID %in% AtORFs))$meanHpc.FPKM)+0.01, ((subset(RPKMmerge, orfID %in% AtORFs))$meanReseda.FPKM)+0.01, pch=20, log="xy", xlim=c(0.01,50000), ylim=c(0.01,50000), xlab="infected", ylab="non-infected", main="Mean FPKMs plant, pathogen and shared ORFs", col="green3")
points(((subset(RPKMmerge, orfID %in% HpaORFs))$meanHpc.FPKM)+0.01, ((subset(RPKMmerge, orfID %in% HpaORFs))$meanReseda.FPKM)+0.01, pch=20, log="xy", xlim=c(0.01,50000), ylim=c(0.01,50000), col="gray48")
points(((subset(RPKMmerge, orfID %in% AtHpaORFs))$meanHpc.FPKM)+0.01, ((subset(RPKMmerge, orfID %in% AtHpaORFs))$meanReseda.FPKM)+0.01, pch=20, log="xy", xlim=c(0.01,50000), ylim=c(0.01,50000), col="orange")
dev.off()
#dev.print(pdf, "meanRPKMcomparison.Arabidopsis.Hpa.pdf", width=8, height=8)

#################################
# Add OTU info to summary table # !!!verbose script, needs improvement!!!
#################################
setwd("/home/student1/Day5/03_OTUpicking/")
Bclust <- read.table("uniORFsBlastClust.txt", sep=" ", header=F, fill=T)
Bclust$otuID <- c(paste0("OTU_", 1:(dim(Bclust)[1])))
#Bclust.list <- as.list(as.data.frame(t(read.table("uniORFsBlastClust.txt", sep=" ", header=F, fill=T, row.names=paste0("OTU_", 1:(dim(Bclust)[1]))))))
for(i in 1:dim(uniORFsummary)[1]){
  v1 <- subset(Bclust, V1 == as.character(uniORFsummary$orfID[i]))
  v2 <- subset(Bclust, V2 == as.character(uniORFsummary$orfID[i]))
  v3 <- subset(Bclust, V3 == as.character(uniORFsummary$orfID[i]))
  v4 <- subset(Bclust, V4 == as.character(uniORFsummary$orfID[i]))
  v5 <- subset(Bclust, V5 == as.character(uniORFsummary$orfID[i]))
  v6 <- subset(Bclust, V6 == as.character(uniORFsummary$orfID[i]))
  v7 <- subset(Bclust, V7 == as.character(uniORFsummary$orfID[i]))
  v8 <- subset(Bclust, V8 == as.character(uniORFsummary$orfID[i]))
  v9 <- subset(Bclust, V9 == as.character(uniORFsummary$orfID[i]))
  v10 <- subset(Bclust, V10 == as.character(uniORFsummary$orfID[i]))
  v11 <- subset(Bclust, V11 == as.character(uniORFsummary$orfID[i]))
  if(dim(v1)[1] == 1){
    uniORFsummary[i,"otuID"] <- v1[1,"otuID"]
  }
  if(dim(v2)[1] == 1){
    uniORFsummary[i,"otuID"] <- v2[1,"otuID"]
  }
  if(dim(v3)[1] == 1){
    uniORFsummary[i,"otuID"] <- v3[1,"otuID"]
  }
  if(dim(v4)[1] == 1){
    uniORFsummary[i,"otuID"] <- v4[1,"otuID"]
  }
  if(dim(v5)[1] == 1){
    uniORFsummary[i,"otuID"] <- v5[1,"otuID"]
  }
  if(dim(v6)[1] == 1){
    uniORFsummary[i,"otuID"] <- v6[1,"otuID"]
  }
  if(dim(v7)[1] == 1){
    uniORFsummary[i,"otuID"] <- v7[1,"otuID"]
  }
  if(dim(v8)[1] == 1){
    uniORFsummary[i,"otuID"] <- v8[1,"otuID"]
  }
  if(dim(v9)[1] == 1){
    uniORFsummary[i,"otuID"] <- v9[1,"otuID"]
  }
  if(dim(v10)[1] == 1){
    uniORFsummary[i,"otuID"] <- v10[1,"otuID"]
  }
  if(dim(v11)[1] == 1){
    uniORFsummary[i,"otuID"] <- v11[1,"otuID"]
  }
  rm(v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11)
}


#######################################################################################################
# Analysis of differential expression of plant OTUs comparing infected and non-infected plant samples #
#######################################################################################################
# Upload htseqcount data from samples
setwd("/home/student1/Day5/01_Data/")
hpc21Count <- read.table("Hpc21_htseq.tab", sep="\t", header=F, col.names=c("orfID", "count"))
hpc31Count <- read.table("Hpc31_htseq.tab", sep="\t", header=F, col.names=c("orfID", "count"))
hpc42Count <- read.table("Hpc42_htseq.tab", sep="\t", header=F, col.names=c("orfID", "count"))
R13Count <- read.table("Reseda13_htseq.tab", sep="\t", header=F, col.names=c("orfID", "count"))
R21Count <- read.table("Reseda21_htseq.tab", sep="\t", header=F, col.names=c("orfID", "count"))
R31Count <- read.table("Reseda31_htseq.tab", sep="\t", header=F, col.names=c("orfID", "count"))

# Add OTU IDs to count tables
#Hpc21
d1 <- hpc21Count
colnames(d1) <- c("orfID", "countH21")
# Add OTU ID to count table
for(i in 1:dim(d1)[1]){
  d2 <- subset(uniORFsummary, orfID %in% d1[i,"orfID"])
  d1[i,"otuID"] <- d2[1,"otuID"]
  rm(d2)
}
hpc21Count <- d1
rm(d1)

#Hpc31
d1 <- hpc31Count
colnames(d1) <- c("orfID", "countH31")
# Add OTU ID to count table
for(i in 1:dim(d1)[1]){
  d2 <- subset(uniORFsummary, orfID %in% d1[i,"orfID"])
  d1[i,"otuID"] <- d2[1,"otuID"]
  rm(d2)
}
hpc31Count <- d1
rm(d1)

#Hpc42
d1 <- hpc42Count
colnames(d1) <- c("orfID", "countH42")
# Add OTU ID to count table
for(i in 1:dim(d1)[1]){
  d2 <- subset(uniORFsummary, orfID %in% d1[i,"orfID"])
  d1[i,"otuID"] <- d2[1,"otuID"]
  rm(d2)
}
hpc42Count <- d1
rm(d1)

d2 <- merge(hpc21Count[,c("orfID", "countH21")], hpc31Count[,c("orfID", "countH31")], by.x="orfID", by.y="orfID", all=T)
d3 <- merge(d2, hpc42Count[,c("orfID", "countH42", "otuID")], by.x="orfID", by.y="orfID", all=T)
hpcCount <- d3
rm(d2,d3)

# Reseda13
d1 <- R13Count
colnames(d1) <- c("orfID", "countR13")
# Add OTU ID to count table
for(i in 1:dim(d1)[1]){
  d2 <- subset(uniORFsummary, orfID %in% d1[i,"orfID"])
  d1[i,"otuID"] <- d2[1,"otuID"]
  rm(d2)
}
R13Count <- d1

# Reseda21
d1 <- R21Count
colnames(d1) <- c("orfID", "countR21")
# Add OTU ID to count table
for(i in 1:dim(d1)[1]){
  d2 <- subset(uniORFsummary, orfID %in% d1[i,"orfID"])
  d1[i,"otuID"] <- d2[1,"otuID"]
  rm(d2)
}
R21Count <- d1

# Reseda31
d1 <- R31Count
colnames(d1) <- c("orfID", "countR31")
# Add OTU ID to count table
for(i in 1:dim(d1)[1]){
  d2 <- subset(uniORFsummary, orfID %in% d1[i,"orfID"])
  d1[i,"otuID"] <- d2[1,"otuID"]
  rm(d2)
}
R31Count <- d1

d2 <- merge(R13Count[,c("orfID", "countR13")], R21Count[,c("orfID", "countR21")], by.x="orfID", by.y="orfID", all.x=T)
d3 <- merge(d2, R31Count[,c("orfID", "countR31", "otuID")], by.x="orfID", by.y="orfID", all.x=T)
RCount <- d3
rm(d2, d3, d1)

# Create new count table with IDs and summed counts of plant OTUs
plantOTUs <- data.frame(otuID = (subset(uniORFsummary, HpaHomolog == 0 & AtHomolog == 1))$otuID) # unique OTUs with homologs only in A. thaliana but not Hpa
plantOTUs <- plantOTUs[!duplicated(plantOTUs$otuID),]
allOTUs <- (uniORFsummary[!duplicated(uniORFsummary$otuID),])$otuID
df1 <- as.data.frame(matrix(0,nrow=length(plantOTUs), ncol=7))
colnames(df1) <- c("otuID", "countSumH21", "countSumH31", "countSumH42", "countSumR13", "countSumR21", "countSumR31")
for(i in 1:length(plantOTUs)){
  d4 <- subset(hpcCount, otuID == as.character(plantOTUs[i]))
  d5 <- subset(RCount, otuID == as.character(plantOTUs[i]))
  if(dim(d4)[1] > 0){
    df1[i,"otuID"] <- as.character(plantOTUs[i])
    df1[i,"countSumH21"] <- sum(d4$countH21)
    df1[i,"countSumH31"] <- sum(d4$countH31)
    df1[i,"countSumH42"] <- sum(d4$countH42)
  }
  rm(d4)
  if(dim(d5)[1] > 0){
    df1[i,"otuID"] <- as.character(plantOTUs[i])
    df1[i,"countSumR13"] <- sum(d5$countR13)
    df1[i,"countSumR21"] <- sum(d5$countR21)
    df1[i,"countSumR31"] <- sum(d5$countR31)
  }
  rm(d5)   
}
SumCountPlant <- df1
rm(df1)

# Save OTU raw count tables to be used by DESeq2
setwd("/home/student1/Day5/05_DiffExpress")
write.table(SumCountPlant[,c("otuID", "countSumH21")], "Hpc21.plantOTU.csv", sep="\t", row.names=F, quote=F, col.names=F)
write.table(SumCountPlant[,c("otuID", "countSumH31")], "Hpc31.plantOTU.csv", sep="\t", row.names=F, quote=F, col.names=F)
write.table(SumCountPlant[,c("otuID", "countSumH42")], "Hpc42.plantOTU.csv", sep="\t", row.names=F, quote=F, col.names=F)
write.table(SumCountPlant[,c("otuID", "countSumR13")], "Reseda13.plantOTU.csv", sep="\t", row.names=F, quote=F, col.names=F)
write.table(SumCountPlant[,c("otuID", "countSumR21")], "Reseda21.plantOTU.csv", sep="\t", row.names=F, quote=F, col.names=F)
write.table(SumCountPlant[,c("otuID", "countSumR31")], "Reseda31.plantOTU.csv", sep="\t", row.names=F, quote=F, col.names=F)


###################
# DESeq2 analysis #
###################
source("http://www.bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
# FYI, you can always pull up the vignette appropriate for your version by loading it from your R session
# vignette("DESeq2")

##############################################################################################
# Upload htseqCount results as a RangedSummarizedExperiment object for analysis using DESeq2 #
##############################################################################################
# Htseq Sample Table for data upload in DESeq2 #
sampleID=c("Hpc21","Hpc31","Hpc42","Reseda13","Reseda21","Reseda31")
fileID=c("Hpc21.plantOTU.csv", "Hpc31.plantOTU.csv", "Hpc42.plantOTU.csv", "Reseda13.plantOTU.csv", "Reseda21.plantOTU.csv", "Reseda31.plantOTU.csv")
condition=c("infected","infected","infected","noninfected","noninfected","noninfected")
htseqTable=data.frame(sampleID,fileID,condition) 

# Create RangedSummarizedExperiment object
setwd("/home/student1/Day5/05_DiffExpress/")
dds <- DESeqDataSetFromHTSeqCount(sampleTable=htseqTable, directory="/home/student1/Day5/05_DiffExpress/", design = ~ condition)

# Pre-filtering of low count genes
dds <- dds[ rowSums(counts(dds)) > 1, ]

#########################################
# Heatmap of sample-to-sample distances #
#########################################
# Regularised log transformation
# The function rlog, stands for regularized log, transforming the original count data to the log2 scale 
# by fitting a model with a term for each sample and a prior distribution on the coefficients which is estimated from the data.
rld <- rlog(dds, blind=FALSE)
# Here, we apply the dist function to the transpose of the transformed count matrix to get sample-to-sample distances.
sampleDists <- dist(t(assay(rld)))
# Plot heatmap
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, colnames(rld), sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("SampleDistances.pdf", width=6, height=4)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

# Principal component plot of the samples
setwd("/home/student1/Day5/05_DiffExpress/")
pdf("PCA_OTUs.pdf", width=8, height=6)
plotPCA(rld, intgroup="condition")
dev.off()
# Alternatively you can customize the PCA plot using the ggplot function
# data <- plotPCA(rld, intgroup="condition", returnData=TRUE)
# percentVar <- round(100 * attr(data, "percentVar"))
# library("ggplot2")
# ggplot(data, aes(PC1, PC2, color=condition, shape=condition)) +
#   geom_point(size=3) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#   coord_fixed()


###############################################
# Differential Expression analysis including: #
# - estimating size factors                   #
# - estimating dispersions                    #
# - gene-wise dispersion estimates            #
# - mean-dispersion relationship              #
# - final dispersion estimates                #
# - fitting model and testing                 #
###############################################
dds <- DESeq(dds)
res <- results(dds)
summary(res)

# Order results by log2 fold change
res <- res[order(res$log2FoldChange),]

# Filter DESeq2 results by different adjusted p-value thresholds:
res0.05 <- subset(res, padj < 0.05)
summary(res0.05)
res0.01 <- subset(res, padj < 0.01)
summary(res0.01)
res0.001 <- subset(res, padj < 0.001)
summary(res0.001)
res0.001
sum(res$padj < 0.01, na.rm=TRUE)
Top20padj0.01 <- rbind(as.data.frame(head(res0.01, n=20)),as.data.frame(tail(res0.01, n=20)))
Top20padj0.01$otuID <- rownames(Top20padj0.01)
Top20padj0.01 <- merge(Top20padj0.01, uniORFsummary[,c("otuID", "AtQuery", "orfID")], by.x="otuID", by.y="otuID", all.x=T)
Top20padj0.01 <- Top20padj0.01[order(Top20padj0.01$log2FoldChange),]
setwd("/home/student1/Day5/05_DiffExpress/")
write.table(Top20padj0.01, "Top20_up_down_regulatedOTUs.padj0.01.csv", sep="\t", row.names=F)

# MA-Plot: showing the log2 fold changes attributable to a given variable over the mean of normalized counts.
# Points will be colored red if the adjusted p value is less than 0.1. 
# Points which fall out of the window are plotted as open triangles pointing either up or down.
setwd("/home/student1/Day5/05_DiffExpress/")
pdf("MAplot.pdf", width=8, height=6)
plotMA(res, main="Plant OTUs - non-infected vs. infected", ylim=c(-6,6), xlab="mean FPKM")
dev.off()

