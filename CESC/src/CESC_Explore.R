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

mat3 <- fread("GEO/03/matrix03.txt",sep="\t",header = T,stringsAsFactors = F)
mat3 <- as.data.frame(mat3)
rownames(mat3) <- mat3[,1]

probInf              <- fread(file='F:/projects/TCGA/data/TCGA_450k_sample.txt',sep="\t",header=T,check.names=F,stringsAsFactors=F)[,-2];
probInf              <- as.data.frame(probInf)
colnames(probInf)[1] <- "probeId"
rownames(probInf)    <- probInf$probeId
probInf              <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]
print('450K probe information loaded!')

int1=intersect(mat1$ID_REF,probInf$probeId)
int2=intersect(mat2$ID_REF,probInf$probeId)
int3=intersect(mat3$ID_REF,probInf$probeId)
ints <- intersect(intersect(int1,int2),int3)

mat1 <- mat1[ints,-1]
mat2 <- mat2[ints,-1]
mat3 <- mat3[ints,-1]

mat <- rbind(t(mat1),t(mat2))
mat <- data.frame(rbind(mat,t(mat3)))
clnames <- colnames(mat)

sampleInf1 <- read.table('GEO/01/sampleInf.txt',sep='\t',header=T,stringsAsFactors=F)
sampleInf1$label <- tolower(sampleInf1$label)

sampleInf2 <- read.table('GEO/02/sampleInf.txt',sep='\t',header=T,stringsAsFactors=F)
sampleInf2$label <- tolower(sampleInf2$label)
sampleInf3 <- read.table('GEO/03/sampleInf.txt',sep='\t',header=T,stringsAsFactors=F)
sampleInf3$label <- tolower(sampleInf3$label)
colnames(sampleInf3)[1] <- 'sampleId';

sampleInf <- rbind(sampleInf1[,c(1,2)],sampleInf2[,c(1,2)])
sampleInf <- rbind(sampleInf,sampleInf3[,c(1,2)])
colnames(sampleInf) <- c('sampleId','label')
sampleInf$src <- c(rep('paper1',dim(sampleInf1)[1]),rep('paper2',dim(sampleInf2)[1]),rep('paper3',dim(sampleInf3)[1]))
#############################################################################################################################
rownames(sampleInf) <- sampleInf$sampleId
mat$label <- sampleInf[rownames(mat),'src']
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
color3 <- c('paper1'='red','paper2'='black','paper3'='blue')
drawUMAP(f1(20),color3)

