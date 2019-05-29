#requirement:
# ZNF454: 
#      chr5:178367763-178368288  chr5:178940762-178941287: cg02165355;cg03234732;cg10575261;cg23037403;cg16536329;cg20778451;cg17840719
#      chr5:178367621-178368725  chr5:178940620-178941724
#		 chr5:178367763-178368725  chr5:178940762-178941724
#		 chr5:178367621-178368288  chr5:178940620-178941287
# CDO1:
#      chr5:115151911-115152561  chr5:115816214-115816864


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

probeInf <- fread('F:/projects/allData/TCGA/probeInf.txt',sep="\t",header=T,stringsAsFactors=F,check.names=F)[,-c(2,7:9)]
probeInf <- as.data.frame(probeInf)

ZNF454 <- probeInf[grepl('ZNF454;|;ZNF454',probeInf$Gene_Symbol) & probeInf$Feature_Type=='Island',]
ZNF454 <- ZNF454[order(ZNF454$Start),]
CDO1   <- probeInf[grepl('CDO1;|;CDO1',probeInf$Gene_Symbol) & probeInf$Feature_Type=='Island',]
CDO1   <- CDO1[order(CDO1$Start),]

dmr1 <-              c("cg24843380","cg02165355","cg03234732","cg10575261","cg23037403","cg16536329","cg20778451","cg17840719")
dmr2 <- c("cg17741986","cg24843380","cg02165355","cg03234732","cg10575261","cg23037403","cg16536329","cg20778451","cg17840719","cg03355526","cg10902717")
dmr3 <-              c("cg24843380","cg02165355","cg03234732","cg10575261","cg23037403","cg16536329","cg20778451","cg17840719","cg03355526","cg10902717")
dmr4 <- c("cg17741986","cg24843380","cg02165355","cg03234732","cg10575261","cg23037403","cg16536329","cg20778451","cg17840719")
dmr5 <- c("cg07405021","cg16265906","cg12880658","cg16707405","cg02792792","cg14470895","cg23180938","cg08516516","cg11036833")

usedProbes <- unique(c(dmr1,dmr2,dmr3,dmr4,dmr5)) # 19 probes

mat450 <- mat450[,c('label',usedProbes)]

getROC <- function(dmr){
        v0 <- apply(mat450[,dmr],1,mean)
		  obj <- roc('controls'=v0[mat450$label=='normal'],'cases'=v0[mat450$label!='normal'])
		  obj
}
getRFROC <- function(M0,iter=5){
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
	ROC <- roc(ifelse(predicts$observed=="cancer", "cancer", "normal"), as.numeric(predicts$cancer),direction='>')
   return(ROC)
}
getCombineROC <- function(pbs1,pbs2){
      v1 <- apply(mat450[,pbs1],1,mean)
		v2 <- apply(mat450[,pbs2],1,mean)
      v3 <- data.frame('label'=mat450$label,'v1'=v1,'v2'=v2)
		obj <- getRFROC(v3)
		obj
}
#############################################################
roc1  <- getROC(dmr1)
roc2  <- getROC(dmr2)
roc3  <- getROC(dmr3)
roc4  <- getROC(dmr4)
roc15 <- getCombineROC(dmr1,dmr5)
roc25 <- getCombineROC(dmr2,dmr5)
roc35 <- getCombineROC(dmr3,dmr5)
roc45 <- getCombineROC(dmr4,dmr5)

pdf('output/UCEC_ROC_0528_ChenMY.pdf',width=10,height=12)
   par(mar=c(0.5,0.5,2,0.5))
	layout(matrix(1:4,ncol=2,byrow=T),heights=c(5,5,5,5))
   plot(roc1,print.auc=TRUE,print.thres="best",col='blue',legacy.axes = TRUE,main="ZNF454:DMR1")
	plot(roc2,print.auc=TRUE,print.thres="best",col='blue',legacy.axes = TRUE,main="ZNF454:DMR2")
	plot(roc3,print.auc=TRUE,print.thres="best",col='blue',legacy.axes = TRUE,main="ZNF454:DMR3")
	plot(roc4,print.auc=TRUE,print.thres="best",col='blue',legacy.axes = TRUE,main="ZNF454:DMR4")

   par(mar=c(0.5,0.5,2,0.5))
	layout(matrix(1:4,ncol=2,byrow=T),heights=c(5,5,5,5))
   plot(roc15,print.auc=TRUE,print.thres="best",col='blue',legacy.axes = TRUE,main="ZNF454:DMR1+CDO1")
	plot(roc25,print.auc=TRUE,print.thres="best",col='blue',legacy.axes = TRUE,main="ZNF454:DMR2+CDO1")
	plot(roc35,print.auc=TRUE,print.thres="best",col='blue',legacy.axes = TRUE,main="ZNF454:DMR3+CDO1")
	plot(roc45,print.auc=TRUE,print.thres="best",col='blue',legacy.axes = TRUE,main="ZNF454:DMR4+CDO1")
dev.off()


