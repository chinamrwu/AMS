rm(list=ls())
library(data.table)
library(sqldf)

setwd('F:/projects/STAD')
probInf <- fread('F:/projects/allData/TCGA/probeInf.txt',sep="\t",header=T,stringsAsFactors=F)
probInf <- as.data.frame(probInf)
probInf <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]

DMR <- read.csv('output/STAD_DMR_20190612.csv',header=T,stringsAsFactors=F)

meta27k   <-  read.csv('F:/projects/allData/illumina_humanmethylation27_content.csv',header=T,stringsAsFactors=F)
meta450k  <- read.csv('F:/projects/allData/HumanMethylation450_15017482_v1-2.csv',skip=7,header=T,stringsAsFactors=F)[,c(2,11,12,13,15,16)]
meta450k  <- meta450k[!duplicated(meta450k$Name),]
rownames(meta450k) <- meta450k[,1]
##############################
matGSE30601 <- fread('data/matGSE30601.txt',sep="\t",header=T,stringsAsFactors=F)
matGSE30601 <- as.data.frame(matGSE30601)
matGSE25869 <- fread('data/matGSE25869.txt',sep="\t",header=T,stringsAsFactors=F)
matGSE25869 <- as.data.frame(matGSE25869)

########################################
for(i in 1:dim(DMR)[1]){
  


}