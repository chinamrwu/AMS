rm(list=ls())
library(data.table)
library(sqldf)
library(pROC)
rm(list=ls())
library(data.table)
library(randomForest)
library(sqldf)
library(pROC)
library(caret)

setwd('F:/projects/UCEC')
mat450 <- fread('F:/projects/allData/TCGA/UCEC_450k.txt',sep="\t",header=T,stringsAsFactors=F,check.names=F)
mat450 <- as.data.frame(mat450)
rownames(mat450) <- mat450$probeId
mat450 <- mat450[,-1]
mat450 <- mat450[,grepl('-01A-|-11A-',colnames(mat450))] ##474:428+46
Labels <- as.character(sapply(colnames(mat450),function(v){a <- strsplit(v,"-")[[1]][4];ifelse('01A'==a,'cancer','normal')}))
mat450 <- data.frame(t(mat450),check.names=F)
probeIds <- colnames(mat450)
mat450$label <- Labels
mat450 <- mat450[,c('label',probeIds)]

probeInf <- fread('F:/projects/allData/TCGA/probeInf.txt',sep="\t",header=T,stringsAsFactors=F,check.names=F)
probeInf <- as.data.frame(probeInf)

BHLHE22 <- sqldf("SELECT * from probeInf where Chromosome='chr8' and (Start >= 64578738 and End  <= 64581936)")[,1]
CDO1    <- sqldf("SELECT * from probeInf where Chromosome='chr5' and (Start >= 115814209 and End <= 115818208)")[,1]
CELF4   <- sqldf("SELECT * from probeInf where Chromosome='chr18' and (Start >= 37563538 and End <= 37567537)")[,1]
ZNF662  <- sqldf("SELECT * from probeInf where Chromosome='chr3' and (Start >= 42904666 and End <= 42908665)")[,1]


genes <- list()
genes[[length(genes)+1]] <- BHLHE22; 
genes[[length(genes)+1]] <- CDO1;
genes[[length(genes)+1]] <- CELF4;
genes[[length(genes)+1]] <- ZNF662;
names(genes) <- c('BHLHE22','CDO1','CELF4','ZNF662')


### ROC for single gene
getROC <- function(probes){
  v0 <- as.numeric(apply(mat450[,probes],1,mean,na.rm=T))
  roc1 <- roc('controls'=v0[mat450$label=='normal'],'cases'=v0[mat450$label=='cancer'])
  roc1
}

ROCs <- sapply(genes,getROC,simplify = F)

# BHLHE22和CDO1
ROCs[[length(ROCs)+1]] <- getROC(c(BHLHE22,CDO1))

#CDO1和CELF4
ROCs[[length(ROCs)+1]] <- getROC(c(CDO1,CELF4))

#BHLHE22和CELF4
ROCs[[length(ROCs)+1]] <- getROC(c(BHLHE22,CELF4))

#BHLHE22、CDO1和CELF4
ROCs[[length(ROCs)+1]] <- getROC(c(BHLHE22,CDO1,CELF4))

#CDO1和ZNF662
ROCs[[length(ROCs)+1]] <- getROC(c(CDO1,ZNF662))

names(ROCs) <- c('BHLHE22','CDO1','CELF4','ZNF662','CHLHE22+CDO1','CDO1+CELF4','BHLHE22+CELF4','BHLHE22+CDO1+CELF4','CDO1+ZNF662')

pdf('output/UCEC_ROC_5genes_20190605.pdf',width=10,height=10)
   par(mar=c(0.5,0.5,2,0.5))
	layout(matrix(1:4,ncol=2,byrow=T),heights=c(5,5,5,5))
   for(i in 1:4){
	 plot(ROCs[[i]],print.auc=TRUE,print.thres="best",col='blue',legacy.axes = TRUE,main=names(ROCs)[i],print.auc.cex=2.0,print.thres.cex=1.5)
	}
   par(mar=c(0.5,0.5,2,0.5))
	layout(matrix(1:4,ncol=2,byrow=T),heights=c(5,5,5,5))
   for(i in 5:8){
	 plot(ROCs[[i]],print.auc=TRUE,print.thres="best",col='blue',legacy.axes = TRUE,main=names(ROCs)[i],print.auc.cex=2.0,print.thres.cex=1.5)
	}
	par(mar=c(0.5,0.5,2,0.5))
	layout(matrix(1:4,ncol=2,byrow=T),heights=c(5,5,5,5))
	 plot(ROCs[[9]],print.auc=TRUE,print.thres="best",col='blue',legacy.axes = TRUE,main=names(ROCs)[i],print.auc.cex=2.0,print.thres.cex=1.5)
dev.off()
