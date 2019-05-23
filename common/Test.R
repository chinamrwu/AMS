rm(list=ls())
library(pROC)
library(sqldf)
library(data.table)
setwd('F:/projects/TCGA')

probInf <- read.table('F:/projects/TCGA/data/TCGA_450K_sample.txt',sep="\t",header=T,stringsAsFactors=F)[,-2]
colnames(probInf)[1] <- "probeId"
probInf <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]

COAD <- fread('F:/projects/TCGA/data/COAD_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
READ <- fread('F:/projects/TCGA/data/READ_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
mat450 <- cbind(as.data.frame(COAD),as.data.frame(READ[,-1]));

rm(COAD)
rm(READ)
rownames(mat450) <- mat450$probeId
mat450 <- mat450[,-1]
mat450 <- data.frame(t(mat450),stringsAsFactors=F,check.names=F)
mat450 <- mat450[grepl('-01A-|-11A-',rownames(mat450)),]

mat450$label <- as.character(sapply(rownames(mat450),function(v){a <- strsplit(v,"-")[[1]][4];ifelse(a=='01A','cancer','normal')}))

SDC2  <- sqldf("SELECT * FROM probInf where Gene_Symbol like '%;SDC2;%'")[,-c(6:8,10)]
TFPI2 <- sqldf("SELECT * FROM probInf where Gene_Symbol like '%;TFPI2;%'")[,-c(6:8,10)]

matESCA <- fread('data/ESCA_450k.txt',sep="\t",header=T,stringsAsFactors=F,check.names=F)
matESCA <- as.data.frame(matESCA)
rownames(matESCA) <- matESCA$probeId
matESCA <- matESCA[,-1]
matESCA <- data.frame(t(matESCA),check.name=F,stringsAsFactors=F)
matESCA$label <- as.character(sapply(rownames(matESCA),function(v){a <- strsplit(v,"-")[[1]][4];ifelse(a=='01A','cancer','normal')}))

tmp <- getProbeBlocks('DIDO1',matESCA)
system.time(tmp <- getDMR(mat450))