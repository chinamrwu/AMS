rm(list=ls())
library(data.table)

setwd('F:/projects/TCGA')

mat450 <- fread(file='data/ESCA_450k.txt',sep="\t",header=T,check.names=F,stringsAsFactors=F)
mat450 <- as.data.frame(mat450)
clinic <- read.table('ESAD/data/clinical.tsv',sep="\t",header=T,stringsAsFactors=F)
submitter <- clinic$submitter_id[grepl('squamous',tolower(clinic$primary_diagnosis))]
samples <- as.character(sapply(colnames(mat450)[-1],function(v){a <- strsplit(v,"-")[[1]][1:3];paste0(a,collapse="-")}))
index <- match(submitter,samples)+1
samples <- colnames(mat450)[index]
normalSamples <- colnames(mat450)[grepl('-11A-',colnames(mat450))]
samples <- c(samples,normalSamples)
persons <- as.character(sapply(samples,function(v){a <- strsplit(v,"-")[[1]][1:3];paste0(a,collapse="-")}))
samples <- samples[!duplicated(persons)]
Label <-  as.character(sapply(samples,function(v){strsplit(v,"-")[[1]][4]}))
samples <- samples[Label %in% c('01A','11A')]
Label <-  as.character(sapply(samples,function(v){strsplit(v,"-")[[1]][4]}))

mat <- mat450[,c('probeId',samples)]
#DMR <- findDMR(tmp,probInfPath='F:/projects/TCGA/data/TCGA_450K_sample.txt')
DMR <- read.table("ESAD/output/ESAD_DMR_0509.txt",sep="\t",header = T,stringsAsFactors = F)
library(pROC)
library(randomForest)
library(caret)

rownames(mat) <- mat$probeId
mat <- mat[,-1]
mat <- data.frame(t(mat),check.names=F)

getAUC <- function(M0,iter=5){
	#M0 <- mat[,c(1:8)]
	iter=5
	M0 <- tmp
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

getAUC1 <- function(M0){
   tmp <- data.frame('label'=M0$label,'value'=apply(M0[,colnames(M0)!='label'],1,mean))
   obj <- roc(label ~ value,data=tmp)
   #c(mean(obj$sensitivities),mean(obj$specificities),obj$auc)
	return(obj)
}

CDO1   <- strsplit(DMR$probeIds[grepl('CDO1',DMR$GeneSymbol)],";")[[1]]
NKX2   <- strsplit(DMR$probeIds[grepl('NKX2-6',DMR$GeneSymbol)],";")[[1]]
ZNF454 <- strsplit(DMR$probeIds[grepl('ZNF454',DMR$GeneSymbol)],";")[[1]]

ds <- list('CDO1'=CDO1,'NKX2-6'=NKX2,'ZNF454'=ZNF454)
Label <- as.character(sapply(rownames(mat),function(v){a <- strsplit(v,"-")[[1]][4];ifelse('01A'==a,'Cancer','Normal')}))
mats <- list()
for(obj in ds){
   tmp <- as.numeric(apply(mat[,obj],1,mean,na.rm=T))
	tmp <- data.frame('beta'=tmp,'label'=Label)
	mats[[length(mats)+1]] <- tmp
}

plots <- list()
names(mats) <- names(ds)
for(nm in names(mats)){
  plot(roc(label ~ beta,data=mats[[nm]]),print.auc=TRUE,col='blue',main=sprintf("ESAD: %s",nm))
}

genes <- names(mats)
cbn <- combn(1:3,2)
for(i in 1:3){
  indx0 <- cbn[1,i]
  indx1 <- cbn[2,i]
  tmp <- data.frame(cbind(mats[[indx0]][,c(2,1)],mats[[indx1]][,1]),stringsAsFactors=F)
  colnames(tmp) <- c('label','beta1','beta2')
  obj <- getAUC(tmp)
  plot(obj$ROC,print.auc=TRUE,col='blue',main=sprintf("ESAD: %s",paste0(genes[c(indx0,indx1)],collapse="+")))
}