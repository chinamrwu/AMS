#1. SDC2 Segment1：96493351 and 96493635 ："cg13096260" "cg18719750" "cg24732574" "cg08979737" "cg25070637"
#2. SDC2  segment2: 96493701 and 96493935: "cg08979737" "cg25070637" "cg14538332" "cg16935295"
#3  SDC2     Seg 3: 96494006 and 96494335: "cg14538332" "cg16935295"

#TFPI2 segment1: 93890630 and 93890991:"cg24531255" "cg17338208" "cg14775114" "cg16934178" "cg26739865" "cg22441533" "cg14377593"
#TFPI2 segment2: 93889918 and 93890319: "cg19784477" "cg12973591" "cg22799321"

rm(list=ls())
library(data.table)
library(pROC)
library(randomForest)
library(caret)
setwd("F:/projects/CRC")

SDCProbes <- c('cg13096260;cg18719750;cg24732574;cg08979737;cg25070637','cg08979737;cg25070637;cg14538332;cg16935295')
SDCProbes <- c(SDCProbes,'cg14538332;cg16935295')
TFPI2Probes <- c('cg24531255;cg17338208;cg26739865;cg22441533;cg14377593')
TFPI2Probes <- c(TFPI2Probes,'cg12973591;cg22799321')
probes <- unique(as.character(unlist(sapply(c(SDCProbes,TFPI2Probes),function(v){strsplit(v,";")[[1]]}))))




COAD <- fread('F:/projects/TCGA/data/COAD_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
READ <- fread('F:/projects/TCGA/data/READ_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
mat450 <- cbind(as.data.frame(COAD),as.data.frame(READ[,-1]));

rm(COAD)
rm(READ)
rownames(mat450) <- mat450$probeId
mat450 <- mat450[,-1]
mat450 <- data.frame(t(mat450),stringsAsFactors=F,check.names=F)
mat450 <- mat450[,probes]

mat450 <- mat450[grepl('-01A-|-11A-',rownames(mat450)),] #438 samples:393 cancer vs 45 normal
Label <- as.character(sapply(rownames(mat450),function(v){a <- strsplit(v,"-")[[1]][4];ifelse('01A'==a,'Cancer','Normal')}))
mat450$label <- Label
mat450 <- mat450[,c('label',probes)]


clinic <- read.table('F:/projects/TCGA/data/clinical.tsv',sep="\t",header=T,stringsAsFactors=F)
submitters <- as.character(sapply(rownames(mat450),function(v){paste0(strsplit(v,"-")[[1]][1:3],collapse="-")}))
clinic <- clinic[clinic$submitter_id %in% submitters,c('submitter_id','site_of_resection_or_biopsy','tissue_or_organ_of_origin')] 
rownames(clinic) <- clinic$submitter_id
sampleSite <- data.frame('sampleId'=rownames(mat450),stringsAsFactors=F)
sampleSite$site =as.character(sapply(sampleSite$sampleId,function(v){a <- paste0(strsplit(v,"-")[[1]][1:3],collapse="-");
clinic[a,'site_of_resection_or_biopsy']}))
rownames(sampleSite) <- sampleSite[,1]
sampleSite$label <- "NO"
sampleSite$label[grepl('-01A-',rownames(sampleSite))] <- 'cancer'
sampleSite$label[grepl('-11A-',rownames(sampleSite))] <- 'normal'

############################################################################
getAUC <- function(M0,iter=5){
	#M0 <- mat[,c(1:8)]
	set.seed(1)
	iter=5
	predictions <- matrix(0,nrow=dim(M0)[1],ncol=2)
	rownames(predictions) <- rownames(M0)
	for(i in 1:iter){
		folds <- createFolds(M0$label,5)
		for(fold in folds){
			valids <- M0[fold,]
			trains <- M0[setdiff(1:dim(M0)[1],fold),]
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
	predicts$observed <- M0$label
	ROC <- roc(ifelse(predicts$observed=="Cancer", "Cancer", "Normal"), as.numeric(predicts$Cancer),direction='>')
	v=abs(ROC$thresholds-0.5)
	indx <- which(v==min(v))

	pIndex <- which(M0$label=='Cancer')
	nIndex <- which(M0$label=='Normal')

	TP <- sum(predicts$predicted[pIndex]=='Cancer')
	TN <- sum(predicts$predicted[nIndex]=='Normal')

	sensi <- TP/length(pIndex)
	speci <- TN/length(nIndex)
	print(sprintf('%5.3f %5.3f %5.3f',sensi,speci,ROC$auc))
	return(list('sensitivity'=sensi,'specificity'=speci,'ROC'=ROC))
}
########################################
getPrediction <- function(M0,iter=5){
	#M0 <- mat[,c(1:8)]
	set.seed(1)
	iter=5
	predictions <- matrix(0,nrow=dim(M0)[1],ncol=2)
	rownames(predictions) <- rownames(M0)
	for(i in 1:iter){
		folds <- createFolds(M0$label,5)
		for(fold in folds){
			valids <- M0[fold,]
			trains <- M0[setdiff(1:dim(M0)[1],fold),]
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
	predicts$observed <- M0$label
	predicts$result <- predicts$observed==predicts$predicted
	predicts$site <- clinic[as.character(sapply(rownames(predicts),function(v){a <- strsplit(v,"-")[[1]][1:3];paste0(a,collapse="-")})),2]
	predicts <- predicts[order(predicts$result),]
	return(predicts)
}


###############################
getAUC1 <- function(M0){
   tmp <- data.frame('label'=M0$label,'value'=apply(M0[,colnames(M0)!='label'],1,mean))
   obj <- roc(label ~ value,data=tmp)
   #c(mean(obj$sensitivities),mean(obj$specificities),obj$auc)
	return(obj)
}

#####################################################################################

sitePrediction <- function(site,probeIds,iter=5){
	set.seed(1)
	iter=5
	M0 <- mat450[sampleSite$sampleId[sampleSite$site==site],c('label',probeIds)]
   M0$label <- tolower(M0$label)
   M0[is.na(M0)] <- 0
   
	predictions <- matrix(0,nrow=dim(M0)[1],ncol=2)
	rownames(predictions) <- rownames(M0)
	for(i in 1:iter){
		folds <- createFolds(M0$label,5)
		for(fold in folds){
			valids <- M0[fold,]
			trains <- M0[setdiff(1:dim(M0)[1],fold),]
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
	predicts$observed <- M0$label
	positive <- which(M0$label=='cancer')

	return(c(dim(M0)[1],sum(M0$label=='cancer'),sum(predicts$predicted[positive]=='cancer')))
}
sdc2Names <- c('DMR1','DMR2','DMR3')
indexs <- combn(1:3,2)

for(i in 1:3){ sdc2Names <- c(sdc2Names,sprintf('DMR%d+DMR%d',indexs[1,i],indexs[2,i]))}

sitePred <- c()

for(seg in SDCProbes){
   pb <- strsplit(seg,";")[[1]]
   for(site in sites){
     sitePred <- rbind(sitePred,sitePrediction(site,pb))
	}
}

for(i in 1:3){ #### SDC2 3 个 DMR之间的组合
  pb1 <- strsplit(SDCProbes[indexs[1,i]],";")[[1]]
  pb2 <- strsplit(SDCProbes[indexs[2,i]],";")[[1]]

  k1 <- apply(mat450[,pb1],1,mean,na.rm=T)
  k2 <- apply(mat450[,pb2],1,mean,na.rm=T)
  k0 <- data.frame('label'=mat450$label,'value1'=k1,value2=k2)
   sdc2Prediction[[length(sdc2Prediction)+1]] <- getPrediction(k0)
}

names(sdc2Prediction) <- sdc2Names
for(i in 1: length(sdc2Prediction)){
  obj <- sdc2Prediction[[i]]
  write.table(obj,file=sprintf("output/CRC_Prediction_SDC2_%s_20190522.txt",names(sdc2Prediction)[i]),sep="\t",col.names=T,row.names=F,quote=F)
}


###########################################################################

tfPrediction <- list()
for(seg in TFPI2Probes){
   pb <- strsplit(seg,";")[[1]]
   k0 <- apply(mat450[,pb],1,mean,na.rm=T)
	k1 <- data.frame('label'=mat450$label,'value'=k0)
   tfPrediction[[length(tfPrediction)+1]] <- getPrediction(k1)
}

pb1 <- strsplit(TFPI2Probes[1],";")[[1]]
pb2 <- strsplit(TFPI2Probes[2],";")[[1]]
k1 <- apply(mat450[,pb1],1,mean,na.rm=T)
k2 <- apply(mat450[,pb2],1,mean,na.rm=T)
k0 <- data.frame('label'=mat450$label,'value1'=k1,value2=k2)
tfPrediction[[length(tfPrediction)+1]] <- getPrediction(k0)

tfNames <- c('DMR1','DMR2','DMR1+DMR2')
names(sdc2Prediction) <- sdc2Names
for(i in 1: length(sdc2Prediction)){
  obj <- sdc2Prediction[[i]]
  write.table(obj,file=sprintf("output/CRC_Prediction_SDC2_%s_20190522.txt",names(sdc2Prediction)[i]),sep="\t",col.names=T,row.names=F,quote=F)
}
############################## SDC2 3个DMR 及其组合 在不同site的检出率

sdc2Names <- c('DMR1','DMR2','DMR3')
indexs <- combn(1:3,2)
for(i in 1:3){ sdc2Names <- c(sdc2Names,sprintf('DMR%d+DMR%d',indexs[1,i],indexs[2,i]))}
sdc2Prediction <- list()
for(seg in SDCProbes){
   
   pb <- strsplit(seg,";")[[1]]
   k0 <- apply(mat450[,pb],1,mean,na.rm=T)
	k1 <- data.frame('label'=mat450$label,'value'=k0,stringsAsFactors=F)
	#k1$result <- k1$label==k1$predicted
   sdc2Prediction[[length(sdc2Prediction)+1]] <- getPrediction(k1)
}

for(i in 1:3){ #### SDC2 3 个 DMR之间的组合
  pb1 <- strsplit(SDCProbes[indexs[1,i]],";")[[1]]
  pb2 <- strsplit(SDCProbes[indexs[2,i]],";")[[1]]

  k1 <- apply(mat450[,pb1],1,mean,na.rm=T)
  k2 <- apply(mat450[,pb2],1,mean,na.rm=T)
  k0 <- data.frame('label'=mat450$label,'value1'=k1,value2=k2)
   sdc2Prediction[[length(sdc2Prediction)+1]] <- getPrediction(k0)
}

names(sdc2Prediction) <- sdc2Names
for(i in 1: length(sdc2Prediction)){
  obj <- sdc2Prediction[[i]]
  write.table(obj,file=sprintf("output/CRC_Prediction_SDC2_%s_20190522.txt",names(sdc2Prediction)[i]),sep="\t",col.names=T,row.names=F,quote=F)
}

