rm(list=ls())
library(randomForest)
library(glmnet)
library(data.table)
library(readr)

library(pROC)

#setwd('F:/projects/TCGA')
setwd("/home/zcwu/")
DMR <- read.table('output/CRCA_DMR_0424.txt',sep='\t',header=T,stringsAsFactors=F,check.names=F);
DMR$rgId <- as.character(sapply(1:dim(DMR)[1],function(i){paste0('rg',i)}))

COAD <- fread('data/COAD_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
READ <- fread('data/READ_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
mat450 <- cbind(as.data.frame(COAD),as.data.frame(READ[,-1]));

rm(COAD)
rm(READ)
rownames(mat450) <- mat450$probeId
mat450 <- mat450[,-1]
mat450 <- data.frame(t(mat450),stringsAsFactors=F,check.names=F)

probInf              <- fread(file='data/TCGA_450k_sample.txt',sep="\t",header=T,check.names=F,stringsAsFactors=F)[,-2];
probInf              <- as.data.frame(probInf)
colnames(probInf)[1] <- "probeId"
rownames(probInf)    <- probInf$probeId
probInf              <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]
print('450K probe information loaded!')

mat450 <- mat450[,probInf$probeId]
R0  <- apply(mat450,2,function(v){sum(is.na(v))})
mat450 <- mat450[,R0==0]
probeIds <- colnames(mat450)
probInf <- probInf[probeIds,]

sampleTypes <- as.character(sapply(rownames(mat450),function(x){x1 <- strsplit(x,"-")[[1]][4];x1}))
mat450 <- mat450[sampleTypes %in% c('01A','11A'),]


X <- sapply(DMR$probeIds,function(v){
     a <- strsplit(v,";")[[1]]
	  as.numeric(apply(mat450[,a],1,mean))
})

colnames(X) <- DMR$rgId
rownames(X) <- rownames(mat450)
X           <- data.frame(X,stringsAsFactors=F)
X$label     <-  as.character(sapply(rownames(X),function(v){a <- strsplit(v,"-")[[1]][4];ifelse('01A'==a,'Cancer','Normal')}))
X <- X[,c('label',DMR$rgId)]
X$label <- as.factor(X$label)




#mat450$label <- as.character(sapply(rownames(mat450),function(v){a <- strsplit(v,"-")[[1]][4];ifelse('01A'==a,'Cancer','Normal')}))
#mat450 <- mat450[,c('label',probeIds)]
#mat450$label <- as.factor(mat450$label)

if(F){ # example code to using LASSO

	cvfit  <- cv.glmnet(as.matrix(X[,-1]),X$label,family='binomial',alpha=0.1,type.measure='class');
	cf <- coef(cvfit,s='lambda.min')
	cf <- rownames(cf)[cf[,1]!=0][-1]
   
	X1 <- X[,c('label',cf)]
	objRF <- randomForest(label ~ . ,data=X1,importance=T,ntree=1000,nodesize=7,na.action=na.exclude)
	imps  <- data.frame(importance(objRF));
	impScore <- imps$MeanDecreaseAccuracy
	imps <- imps[order(impScore,decreasing=T),]
	importance <- rownames(imps)

   #############

   cvfit1  <- cv.glmnet(as.matrix(mat450),X$label,family='binomial',alpha=0.1,type.measure='class');
	cf1 <- coef(cvfit1,s='lambda.min')
	cf1 <- rownames(cf1)[cf1[,1]!=0][-1]

   
	X1 <- mat450[,cf1]
	X1$label <- X$label

	objRF <- randomForest(label ~ . ,data=X1,importance=T,ntree=1000,nodesize=7,na.action=na.exclude)
	imps  <- data.frame(importance(objRF));
	impScore <- imps$MeanDecreaseAccuracy
	imps <- imps[order(impScore,decreasing=T),]
	importance1 <- rownames(imps)



}
