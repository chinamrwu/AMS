rm(list=ls())
library(randomForest)
library(pROC)
library(data.table)
library(caret)
library(sqldf)
setwd('F:/projects/CRC')

COAD <- fread('F:/projects/TCGA/data/COAD_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
READ <- fread('F:/projects/TCGA/data/READ_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
mat450 <- cbind(as.data.frame(COAD),as.data.frame(READ[,-1]));

rm(COAD)
rm(READ)
rownames(mat450) <- mat450$probeId
mat450 <- mat450[,-1]
mat450 <- data.frame(t(mat450),stringsAsFactors=F,check.names=F)

probInf              <- fread(file='F:/projects/TCGA/data/TCGA_450k_sample.txt',sep="\t",header=T,check.names=F,stringsAsFactors=F)[,-2];
probInf              <- as.data.frame(probInf)
colnames(probInf)[1] <- "probeId"
rownames(probInf)    <- probInf$probeId
probInf              <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]

mat450 <- mat450[,probInf$probeId]
R0  <- apply(mat450,2,function(v){sum(is.na(v))})
mat450 <- mat450[,R0==0]
probeIds <- colnames(mat450)
probInf <- probInf[probeIds,]

sampleTypes <- as.character(sapply(rownames(mat450),function(x){x1 <- strsplit(x,"-")[[1]][4];x1}))
mat450 <- mat450[sampleTypes %in% c('01A','11A'),]

probs <- unique(read.table('CRC_DMR_fromXihui_20190506.txt',sep="\t",header=T,stringsAsFactors=F))
DMR0 <- read.table('F:/projects/CRC/output/CRCA_DMR_0424.txt',sep="\t",header=T,stringsAsFactors=F)
tmp <- sapply(DMR0$probeIds,function(v){strsplit(v,";")[[1]]})
index <- 0
for(obj in unique(probs$GeneSymbol)){
   v0 <- probs$probeId[probs$GeneSymbol==obj]
	for(k in 1:length(tmp)){
     if(length(intersect(v0,tmp[[k]])) >0 ) {index <- c(index,k);break}
	}
}
DMR <- DMR0[index,]
rownames(DMR) <- 1:dim(DMR)[1]
DMR$GeneSymbol[10] <- 'MARCH11'
DMR$GeneSymbol[11] <- 'FLI1+SENCR'
DMR$GeneSymbol[13] <- 'PTPRT'
DMR <- rbind(DMR,DMR0[DMR0$GeneSymbol=='SDC2',])

label <-  as.character(sapply(rownames(mat450),function(v){a <- strsplit(v,"-")[[1]][4];ifelse('01A'==a,'Cancer','Normal')}))

X <- sapply(DMR$probeIds,function(v){
     a <- strsplit(v,";")[[1]]
	  as.numeric(apply(mat450[,a],1,mean))
})

colnames(X) <- DMR$GeneSymbol
rownames(X) <- rownames(mat450)
X           <- data.frame(X,stringsAsFactors=F)
X$label     <-  label
X <- X[,c(dim(X)[2],1:(dim(X)[2]-1))]

###############################

library(caret)
for(gene in colnames(X)[- c(1,14)]){
	set.seed(1);
	X1 <- X[,c('label','SDC2',gene)]
	predictions <- matrix(0,nrow=dim(X1)[1],ncol=2)
	rownames(predictions) <- rownames(X1)
	for(i in 1:100){
		folds <- createFolds(X1$label,10)
		for(fold in folds){
			valids <- X1[fold,]
			trains <- X1[setdiff(1:dim(X1)[1],fold),]
			trains$label <- as.factor(trains$label)
			tmpRF <- randomForest(label ~ . ,data=trains,importance=T,ntree=1000,nodesize=5)
			predicted <- predict(tmpRF,valids,type='prob')
			predictions[rownames(predicted),] <- predictions[rownames(predicted),]+predicted
		}
	}
	colnames(predictions) <- colnames(predicted)

	predicts <- t(apply(predictions,1,function(v){v/sum(v)}))
	colnames(predicts) <- colnames(predicted)
	predicts <- data.frame(predicts,check.names=F)
	predicts$predicted <- apply(predicts,1,function(v){names(v)[max(v)==v]})
	predicts$observed <- X1$label

	ROC <- roc(ifelse(predicts$observed=="Cancer", "Cancer", "Normal"), as.numeric(predicts$Cancer))
	png(sprintf('output/CRCA_SDC2+%s_ROC_20190506.png',gene))
	plot.roc(ROC,print.auc=T,col = "blue3",ylim=c(0,1), print.thres="no",	
	  main=sprintf('SDC2+%s',gene),legacy.axes =TRUE,print.auc.cex=1.2)
   dev.off()
}
