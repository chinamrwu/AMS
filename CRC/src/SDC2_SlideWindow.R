
rm(list=ls())
library(sqldf)
library(data.table)
library(randomForest)
library(caret)
library(pROC)
setwd("F:/projects/TCGA")

probInf              <- fread(file='data/TCGA_450k_sample.txt',sep="\t",header=T,check.names=F,stringsAsFactors=F)[,-2];
probInf              <- as.data.frame(probInf)
colnames(probInf)[1] <- "probeId"
rownames(probInf)    <- probInf$probeId
probInf              <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.',]# & probInf$Feature_Type=='Island',]

SDC2 <- sqldf("SELECT * FROM probInf where Gene_Symbol like '%SDC2%' and Chromosome='chr8' and Feature_Type !='.' order by Start")[,-c(6:8)]
SDC2 <- SDC2[order(SDC2$Start),]
#w0 <- SDC2$End[2:dim(SDC2)[1]] - SDC2$Start[1:(dim(SDC2)[1]-1)]
#SDC2 <- SDC2[1:which(w0 > 500)[1],]   ## 


COAD <- fread('data/COAD_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
READ <- fread('data/READ_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
mat450 <- cbind(as.data.frame(COAD),as.data.frame(READ[,-1]));

mat <- mat450[mat450$probeId %in% SDC2$probeId,]
R0  <- apply(mat,2,function(v){sum(is.na(v))})
mat <- mat[,R0==0]
rownames(mat) <- mat$probeId
mat <- data.frame(t(mat[,-1]))
sampleType <- as.character(sapply(rownames(mat),function(v){strsplit(v,"-")[[1]][4]}))
mat <- mat[sampleType %in% c('01A','11A'),]
probes <- SDC2$probeId
mat$label <- as.character(sapply(rownames(mat),function(v){
 ifelse(strsplit(v,"-")[[1]][4]=='01A',"Cancer","Normal")
}))
mat <- mat[,c('label',probes)]
set.seed(1)

getAUC <- function(M0,iter=5){
      #M0 <- mat[,c(1:8)]
		#iter=5
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
		ROC <- roc(ifelse(predicts$observed=="Cancer", "Cancer", "Normal"), as.numeric(predicts$Cancer))
      v=abs(ROC$thresholds-0.5)
      indx <- which(v==min(v))

		pIndex <- which(M0$label=='Cancer')
		nIndex <- which(M0$label=='Normal')

		TP <- sum(predicts$predicted[pIndex]=='Cancer')
		TN <- sum(predicts$predicted[nIndex]=='Normal')

      sensi <- TP/length(pIndex)
		speci <- TN/length(nIndex)
       

		print(sprintf('%5.3f %5.3f %5.3f',sensi,speci,ROC$auc))
		return(c(sensi,speci,ROC$auc))
}

getAUC1 <- function(M0){
  tmp  <- data.frame('label'=M0$label,'value'=apply(M0[,colnames(M0)!='label'],2,mean))
  obj  <- roc(label ~ value,data=tmp)
  sensitivity <- mean(obj$sensitivities)
  specificity <- mean(obj$specificities)
  return(c(sensitivity,specificity,obj$auc)
}
fillGap <- function(M0){
   




}
#aucs <- sapply(probes,function(v){getAUC(mat[,c('label',v)])})
L <- length(probes)
probes <- probes[order(probInf[probes,'Start'])];
k0 <- unlist(sapply(2:8,function(w){
   k1 <- sapply(1:(L-w+1),function(i){
	  pbs <- probes[i:(i+w-1)]; 
     #ac <-  getAUC(mat[,c('label',pbs)])
	  ac  <- getAUC1(mat[,c('label',pbs)])
	  sprintf("%s|%5.3f %5.3f %5.3f",paste0(pbs,collapse="_"),ac[1],ac[2],ac[3])
   })
	k1
}))
k1 <- data.frame(t(sapply(k0,function(v){strsplit(v," |\\|")[[1]];})),stringsAsFactors=F)
rownames(k1) <- 1:dim(k1)[1]
k1$X2 <- as.numeric(k1$X2)
k1 <- k1[order(k1$X2,decreasing=T),]
colnames(k1) <- c('probeIds','AUC')
k1$Start <- sapply(k1$probeIds,function(v){ 
  a <- strsplit(v,"_")[[1]];
  v1 <- probInf[probInf$probeId %in% a,'Start']
  min(v1)
})
k1$End <- sapply(k1$probeIds,function(v){ 
  a <- strsplit(v,"_")[[1]];
  v1 <- probInf[probInf$probeId %in% a,'End']
  max(v1)
})
k1$Width <- abs(k1$End - k1$Start)
pb0 <- strsplit(k1$probeIds[1],"_")[[1]][1]
chr <- probInf[pb0,'Chromosome']
k1$Chr <- rep(chr,dim(k1)[1])
k1$probeNumber <- sapply(k1$probeIds,function(v){length(strsplit(v,"_")[[1]])})

k1 <- k1[,c('Chr','Start','End','probeIds','probeNumber','Width','AUC')]

