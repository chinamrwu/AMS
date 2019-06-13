rm(list=ls())
library(data.table)
library(sqldf)
library(pROC)
library(ggplot2)

setwd('F:/projects/STAD')
source('F:/projects/common/islandDMR.R')
matFiles       <- c('matGSE99553.txt','matGSE103186.txt')
sampleFiles    <- c('GSE99553_sampleInf.txt','GSE103186_sampleInf.txt')

GSE99553   <- fread(paste0('data/',matFiles[1]),sep="\t",header=T,stringsAsFactors=F)
GSE99553   <- as.data.frame(GSE99553)
rownames(GSE99553) <- GSE99553[,1]
GSE99553 <- GSE99553[,-1]
GSE99553 <- data.frame(t(GSE99553))

GSE103186  <- fread(paste0('data/',matFiles[2]),sep="\t",header=T,stringsAsFactors=F)
GSE103186  <- as.data.frame(GSE103186)
rownames(GSE103186) <- GSE103186[,1]
GSE103186 <- GSE103186[,-1]
GSE103186 <- data.frame(t(GSE103186))

sampleInf99553  <- read.table(paste0('data/',sampleFiles[1]),header=T,sep="\t",stringsAsFactors=F)
rownames(sampleInf99553) <- sampleInf99553$acc

sampleInf103186 <- read.table(paste0('data/',sampleFiles[2]),header=T,sep="\t",stringsAsFactors=F)
rownames(sampleInf103186) <- sampleInf103186$acc


probeIds <- colnames(GSE99553)
GSE99553$label <- sampleInf99553[rownames(GSE99553),'status']
GSE99553 <- GSE99553[,c('label',probeIds)]
GSE99553$label['control'==GSE99553$label] <- 'normal'
GSE99553$label['case'==GSE99553$label] <- 'cancer'

probeIds <- colnames(GSE103186)
GSE103186$label <- sampleInf103186[rownames(GSE103186),'label']
GSE103186 <- GSE103186[,c('label',probeIds)]
GSE103186$label[GSE103186$label!='normal'] <- 'cancer'


TCGA <- fread('F:/projects/allData/TCGA/STAD_450k.txt',sep="\t",header=T,stringsAsFactors=F)
TCGA <- as.data.frame(TCGA)
rownames(TCGA) <- TCGA$probeId
TCGA <- TCGA[,-1]
TCGA <- data.frame(t(TCGA),check.names=F)
TCGA <- TCGA[grepl('-01A-|-11A-',rownames(TCGA)),]
probeIds <- colnames(TCGA)
TCGA$label <- as.character(sapply(rownames(TCGA),function(v){a <- strsplit(v,"-")[[1]][4];ifelse(a=='01A','cancer','normal')}))
TCGA <- TCGA[,c('label',probeIds)]
TCGA <- TCGA[TCGA$label=='cancer',]

####################################################################
if(F){
		probInf <- fread('F:/projects/allData/TCGA/probeInf.txt',sep="\t",header=T,stringsAsFactors=F)
		probInf <- as.data.frame(probInf)
		probInf <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]
		probeIds01 <- intersect(probInf$probeId,colnames(GSE99553))
		probInf <- probInf[probInf$probeId %in% probeIds01,]
		tmp01 <- GSE99553[,c('label',probeIds01)]
		DMR1 <- getDMR(GSE99553[,c('label',probeIds01)])
		#############################################################################
		probInf <- fread('F:/projects/allData/TCGA/probeInf.txt',sep="\t",header=T,stringsAsFactors=F)
		probInf <- as.data.frame(probInf)
		probInf <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]
		probeIds02 <- intersect(probInf$probeId,colnames(GSE103186))
		probInf <- probInf[probInf$probeId %in% probeIds02,]
		DMR2 <- getDMR(GSE103186[,c('label',probeIds02)])
}
source('F:/projects/common/common.R')
R0 <- apply(GSE103186,2,function(v){sum(is.na(v))})
library(umap)

L <- dim(GSE103186)[2]
k <- 10
color2 <- c('normal'='blue','cancer'='red')
N <- 20
plots1 <- list()
for(i in 1:N){
  index <- sample(2:L,k,replace=F)
  plots1[[length(plots1)+1]] <- drawUMAP(GSE103186[,c(1,index)],color2,rowNormalization=F,strTitle=sprintf('GSE103186:%s',paste0(colnames(GSE103186)[index],collapse=",")))
}

library(glmnet)

##################################
probes <- intersect(colnames(TCGA),colnames(GSE99553))

tmp1 <-   GSE99553[,probes]
tmp1$label <- 'normal'
tmp2 <-       TCGA[,probes]
tmp <- rbind(tmp1,tmp2)
L <- dim(tmp)[2]
k <- 10
color2 <- c('normal'='blue','cancer'='red')
N <- 20
plots2 <- list()
for(i in 1:N){
  index <- sample(2:L,k,replace=F)
  plots[[length(plots)+1]] <- drawUMAP(tmp[,c(1,index)],color2,rowNormalization=F,strTitle='TCGA+GSE99553')
}
#######################################
probes      <- intersect(colnames(GSE99553),colnames(GSE103186))
tmp1        <- GSE103186[GSE103186$label=='normal',probes]
tmp2        <- GSE99553[,probes]
tmp2$label  <- 'cancer'
tmp         <- rbind(tmp1,tmp2)

L <- dim(tmp)[2]
k <- 10
color2 <- c('normal'='blue','cancer'='red')
N <- 20
plots3 <- list()
for(i in 1:N){
  index <- sample(2:L,k,replace=F)
  plots3[[length(plots3)+1]] <- drawUMAP(tmp[,c(1,index)],color2,rowNormalization=F,strTitle='GSE103186+GSE99553')
}
