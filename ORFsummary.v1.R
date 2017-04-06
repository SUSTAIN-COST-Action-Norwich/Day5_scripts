#####################
# ORF Summary table #
#####################
setwd("/home/student1/Day5/01_Data/")
# Upload DATA (tables with ORFs collected from blast output)
library("seqinr")
R13At <- read.table("Reseda13.AtORFs.tab", sep="\t", header=T)
R13Hpa <- read.table("Reseda13.HpaORF.tab", sep="\t", header=T)
R21At <- read.table("Reseda21.AtORFs.tab", sep="\t", header=T)
R21Hpa <- read.table("Reseda21.HpaORF.tab", sep="\t", header=T)
R31At <- read.table("Reseda31.AtORFs.tab", sep="\t", header=T)
R31Hpa <- read.table("Reseda31.HpaORF.tab", sep="\t", header=T)
h21At <- read.table("Hpc21.AtORFs.tab", sep="\t", header=T)
h21Hpa <- read.table("Hpc21.HpaORFs.tab", sep="\t", header=T)
h31At <- read.table("Hpc31.AtORFs.tab", sep="\t", header=T)
h31Hpa <- read.table("Hpc31.HpaORFs.tab", sep="\t", header=T)
h42At <- read.table("Hpc42.AtORFs.tab", sep="\t", header=T)
h42Hpa <- read.table("Hpc42.HpaORFs.tab", sep="\t", header=T)

# Add sample + blast identifier to dataframe 
R13At$SampleBlastIdent <- "R13At" 
R13Hpa$SampleBlastIdent <- "R13Hpa"
R21At$SampleBlastIdent <- "R21At" 
R21Hpa$SampleBlastIdent <- "R21Hpa"
R31At$SampleBlastIdent <- "R31At" 
R31Hpa$SampleBlastIdent <- "R31Hpa"
h21At$SampleBlastIdent <- "h21At" 
h21Hpa$SampleBlastIdent <- "h21Hpa"
h31At$SampleBlastIdent <- "h31At" 
h31Hpa$SampleBlastIdent <- "h31Hpa"
h42At$SampleBlastIdent <- "h42At" 
h42Hpa$SampleBlastIdent <- "h42Hpa"

# Concatenate dataframes 
ORFs <- rbind(R13At, R13Hpa, R21At, R21Hpa, R31At, R31Hpa, h21At, h21Hpa, h31At, h31Hpa, h42At, h42Hpa)

# Filter ORFs for length >= 150 bp and presence of stop codon
ORFsF <- subset(ORFs, STOP.codon. %in% "YES" & ORF.length >= 150)
#ORFsF <- ORFsF[!duplicated(ORFsF[,c("Hit.ID", "ORF.length", "STOP.codon.", "Hit.to.Query", "Sequence", "SampleBlastIdent")]),]

# Remove identical sequences
#dim(ORFsF[!duplicated(ORFsF$Sequence),])
ORFsU <- ORFsF[!duplicated(ORFsF$Sequence),]

# Add unique identifier for each unique ORF
for(i in 1:dim(ORFsU)[1]){
  ORFsU[i,"orfID"] <- paste0("ORF_", i)
}
summary(ORFsU)

# Add infos for which sample contained which ORF
uniORFsummary <- ORFsU
uniORFsummary[,c("H21_Hpa", "H31_Hpa", "H42_Hpa","H21_At", "H31_At", "H42_At","R13_Hpa", "R21_Hpa", "R31_Hpa", "R13_At", "R21_At", "R31_At")] <- 0    # representing Blast search hits of respective database (H21,...) against Hpa or At proteome; 1 if unique ORF is in this sample, 0 if not
uniORFsummary <- uniORFsummary[,c("orfID", "ORF.length", "STOP.codon.","H21_Hpa", "H31_Hpa", "H42_Hpa","H21_At", "H31_At", "H42_At","R13_Hpa", "R21_Hpa", "R31_Hpa", "R13_At", "R21_At", "R31_At", "Sequence")]
for (i in 1:dim(uniORFsummary)[1]){
  if(uniORFsummary[i,"Sequence"] %in% h21Hpa$Sequence){
    uniORFsummary[i,"H21_Hpa"] <- 1
  }
  if(uniORFsummary[i,"Sequence"] %in% h21At$Sequence){
    uniORFsummary[i,"H21_At"] <- 1
  }
  if(uniORFsummary[i,"Sequence"] %in% h31Hpa$Sequence){
    uniORFsummary[i,"H31_Hpa"] <- 1
  }
  if(uniORFsummary[i,"Sequence"] %in% h31At$Sequence){
    uniORFsummary[i,"H31_At"] <- 1
  }
  if(uniORFsummary[i,"Sequence"] %in% h42Hpa$Sequence){
    uniORFsummary[i,"H42_Hpa"] <- 1
  }
  if(uniORFsummary[i,"Sequence"] %in% h42At$Sequence){
    uniORFsummary[i,"H42_At"] <- 1
  }
  if(uniORFsummary[i,"Sequence"] %in% R13Hpa$Sequence){
    uniORFsummary[i,"R13_Hpa"] <- 1
  }
  if(uniORFsummary[i,"Sequence"] %in% R13At$Sequence){
    uniORFsummary[i,"R13_At"] <- 1
  }
  if(uniORFsummary[i,"Sequence"] %in% R21Hpa$Sequence){
    uniORFsummary[i,"R21_Hpa"] <- 1
  }
  if(uniORFsummary[i,"Sequence"] %in% R21At$Sequence){
    uniORFsummary[i,"R21_At"] <- 1
  }
  if(uniORFsummary[i,"Sequence"] %in% R31Hpa$Sequence){
    uniORFsummary[i,"R31_Hpa"] <- 1
  }
  if(uniORFsummary[i,"Sequence"] %in% R31At$Sequence){
    uniORFsummary[i,"R31_At"] <- 1
  }
}

#Pick full length ORFs with homologs in either A. thaliana or H. arabidopsidis 
setwd("/home/student1/Day5/02_ORFsummary/")
plantORFsF <- subset(ORFsF, SampleBlastIdent == "R13At" | SampleBlastIdent == "R21At" | SampleBlastIdent == "R31At" | SampleBlastIdent == "h21At" | SampleBlastIdent == "h31At" | SampleBlastIdent == "h42At")
plantORFsF <- plantORFsF[!duplicated(plantORFsF$Sequence),]
plantORFsF <- merge(plantORFsF, uniORFsummary[,c("orfID", "Sequence")], by.x="Sequence", by.y="Sequence", all.x=T)
write.table(plantORFsF, "plantORFsF.txt", sep="\t", row.names=F, quote=F)
pathogenORFsF <- subset(ORFsF, SampleBlastIdent == "R13Hpa" | SampleBlastIdent == "R21Hpa" | SampleBlastIdent == "R31Hpa" | SampleBlastIdent == "h21Hpa" | SampleBlastIdent == "h31Hpa" | SampleBlastIdent == "h42Hpa" )
pathogenORFsF <- pathogenORFsF[!duplicated(pathogenORFsF$Sequence),]
pathogenORFsF <- merge(pathogenORFsF, uniORFsummary[,c("orfID", "Sequence")], by.x="Sequence", by.y="Sequence", all.x=T)
write.table(pathogenORFsF, "pathogenORFsF.txt", sep="\t", row.names=F, quote=F)

#Add homology information to summary table
uniORFsummary[,c("HpaHomolog", "AtHomolog")] <- 0
for (i in 1:dim(uniORFsummary)[1]){
  if(as.character(uniORFsummary[i,"Sequence"]) %in% as.character(pathogenORFsF$Sequence)){
    uniORFsummary[i,"HpaHomolog"] <- 1
  }
  if(as.character(uniORFsummary[i,"Sequence"]) %in% as.character(plantORFsF$Sequence)){
    uniORFsummary[i,"AtHomolog"] <- 1
  }
}
t1 <- merge(uniORFsummary, pathogenORFsF[,c("orfID", "Hit.to.Query")], by.x="orfID", by.y="orfID", all.x=T)
colnames(t1) <- c("orfID","ORF.length","STOP.codon.","H21_Hpa","H31_Hpa","H42_Hpa","H21_At","H31_At","H42_At","R13_Hpa","R21_Hpa","R31_Hpa","R13_At","R21_At","R31_At","Sequence","HpaHomolog","AtHomolog","HpaQuery")
t1 <- merge(t1, plantORFsF[,c("orfID", "Hit.to.Query")], by.x="orfID", by.y="orfID", all.x=T)
colnames(t1) <- c("orfID","ORF.length","STOP.codon.","H21_Hpa","H31_Hpa","H42_Hpa","H21_At","H31_At","H42_At","R13_Hpa","R21_Hpa","R31_Hpa","R13_At","R21_At","R31_At","Sequence","HpaHomolog","AtHomolog","HpaQuery","AtQuery")
uniORFsummary <- t1
rm(t1)
setwd("/home/student1/Day5/02_ORFsummary/")
write.table(uniORFsummary, "uniORFsummary.txt", sep="\t", row.names=F)


#Plot sequence length histogram of unique ORFs
homAt <- data.frame(length = (subset(uniORFsummary, HpaHomolog == 0 & AtHomolog == 1))$ORF.length)  # unique ORFs with homologs only in A. thaliana
homHpa <- data.frame(length = (subset(uniORFsummary, HpaHomolog == 1 & AtHomolog == 0))$ORF.length)  # unique ORFs with homologs only in H. arabidopsidis
homBoth <- data.frame(length = (subset(uniORFsummary, HpaHomolog == 1 & AtHomolog == 1))$ORF.length)  # unique ORFs with homologs in both, A. thaliana and H. arabidopsidis
# combine dataframes into one for density plot
homAt$Homology <- "At"
homAt$Colour <- "green"
homHpa$Homology <- "Hpa"
homHpa$Colour <- "grey"
homBoth$Homology <- "Both"
homBoth$Colour <- "orange"
ORFplot <- rbind(rbind(homBoth, homHpa), homAt)
library("ggplot2", lib.loc="~/Library/R/3.2/library")
Plot1 <- ggplot(ORFplot, aes(length)) + theme_bw() +
  geom_histogram(data=subset(ORFplot,Homology=="At"), fill="green3", alpha=0.5, bins = 250) +
  geom_histogram(data=subset(ORFplot,Homology=="Hpa"), fill="gray48", alpha=0.75, bins = 250) +
  geom_histogram(data=subset(ORFplot,Homology=="Both"), fill="orange", alpha=1, bins = 250)
ggsave(filename ="uniORFsLength.hist.pdf", plot = Plot1, width=6, height=4, device="pdf")       # Prints histogram to current working directory


#Write fasta file with unique ORFs for ORF expression analysis
setwd("/home/student1/Day5/04_Mapping/")
library("seqinr")
write.fasta(sequences=as.list(uniORFsummary$Sequence), names=as.list(uniORFsummary$orfID), file.out="uniORFs.fasta", open="w", nbchar=100)

#Create gtf file for read mapping with cufflinks
orfGTF <- as.data.frame(matrix(0,nrow=dim(ORFsU)[1], ncol=9))
colnames(orfGTF) <- c("seqname", "source",  "feature",	"start",	"end",	"score",	"strand",	"frame",	"attribute")
orfGTF$seqname <- paste0("ref",ORFsU$orfID)
orfGTF$source <- "ORFscreen"
orfGTF$feature <- "gene"
orfGTF$start <- 1
orfGTF$end <- ORFsU$ORF.length
orfGTF$score <- 0
orfGTF$strand <- "+"
orfGTF$frame <- 0
for (i in 1:dim(orfGTF)[1]){
  orfGTF[i,"attribute"] <- paste0("gene_id \"", as.character(ORFsU[i,"orfID"]), "\";","transcript_id \"", as.character(ORFsU[i,"orfID"]), "\";")
}
write.table(orfGTF, "uniORFs.gtf", sep="\t", row.names=F, col.names=F, quote=F)