rm(list=ls())
library(data.table)
library(readr)
setwd("F:/projects/CESC")
options('width'=500)
options('digits'=5)

mat1 <- fread("GEO/01/matrix01.txt",sep="\t",header = T,stringsAsFactors = F)
mat1 <- as.data.frame(mat1)
rownames(mat1) <- mat1[,1]

mat2 <- fread("GEO/02/matrix02.txt",sep="\t",header = T,stringsAsFactors = F)
mat2 <- as.data.frame(mat2)
rownames(mat2) <- mat2[,1]

mat3 <- fread("F:/projects/TCGA/data/CESC_450k.txt",sep="\t",header = T,stringsAsFactors = F)
mat3 <- as.data.frame(mat3)
rownames(mat3) <- mat3$probeId

probInf              <- fread(file='F:/projects/TCGA/data/TCGA_450k_sample.txt',sep="\t",header=T,check.names=F,stringsAsFactors=F)[,-2];
probInf              <- as.data.frame(probInf)
colnames(probInf)[1] <- "probeId"
rownames(probInf)    <- probInf$probeId
probInf              <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]
print('450K probe information loaded!')

int1 <- intersect(mat1$ID_REF,probInf$probeId)
int2 <- intersect(mat2$ID_REF,probInf$probeId)
int3 <- intersect(mat3$probeId,probInf$probeId)
ints <- intersect(intersect(int1,int2),int3)

mat1 <- mat1[ints,-1]
mat2 <- mat2[ints,-1]
mat3 <- mat3[ints,-1]


mat <- rbind(t(mat1),t(mat2))
mat <- data.frame(rbind(mat,t(mat3)))

R0  <- apply(mat,2,function(v){sum(is.na(v))})
mat <- mat[,R0==0]
clnames <- colnames(mat)

sampleInf1 <- read.table('GEO/01/sampleInf.txt',sep='\t',header=T,stringsAsFactors=F)
sampleInf1$label <- tolower(sampleInf1$label)

sampleInf2 <- read.table('GEO/02/sampleInf.txt',sep='\t',header=T,stringsAsFactors=F)
sampleInf2$label <- tolower(sampleInf2$label)

sampleInf <- rbind(sampleInf1[,c(1,2)],sampleInf2[,c(1,2)])
colnames(sampleInf) <- c('sampleId','label')
#############################################################################################################################
rownames(sampleInf) <- sampleInf$sampleId


mat$label <- c(rep('paper1',dim(mat1)[2]),rep('paper2',dim(mat2)[2]),rep('paper3',dim(mat3)[2]))
mat <- mat[,c('label',clnames)]

library(umap)
source('F:/src/common.R')
f1 <- function(probes=20){
   tmp <- apply(mat[,sample(2:dim(mat)[2],probes,replace=F)],1,function(v){(v-mean(v))/sd(v)})
	tmp <- data.frame(t(tmp))
	tmp$label <- mat$label
	tmp
}
f2 <- function(probes=20){
   mat[,c(1,sample(2:dim(mat)[2],probes,replace=F))]
} 

f3 <- function(M,featureNumber=20){
     indx <- which(colnames(M)=='label')
     index <- 1:dim(M)[1];
	  indx0 <- sample(index[-indx],featureNumber,replace=F);
     M[,c(indx,indx0)]
}


M0 <- f2(20)
M1 <- as.matrix(M0[,-1])
M1 <- log2(M1/(1-M1))
M1[is.infinite(M1)] <- min(M1[!is.infinite(M1)])
M1 <- data.frame(M1)
M1$label <- M0$label
M1 <- M1[,c('label',colnames(M0)[-1])]

color3 <- c('paper1'='red','paper2'='black','paper3'='blue')
drawUMAP(M0,color3,strTitle="UMAP:raw data",rowNormalization=F)
drawUMAP(M0,color3,strTitle="UMAP:row normalized",rowNormalization=T)
drawUMAP(M0,color3,strTitle="UMAP:row normalized",rowNormalization=T)
drawUMAP(M1,color3,strTitle="UMAP:raw M-value",rowNormalization=F)
drawUMAP(M1,color3,strTitle="UMAP: normalized M-value",rowNormalization=F)

drawUMAP(f2(20),color3)
