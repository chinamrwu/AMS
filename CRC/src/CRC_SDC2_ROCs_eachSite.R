#1. SDC2 Segment1£º96493351 and 96493635 £º"cg13096260" "cg18719750" "cg24732574" "cg08979737" "cg25070637"
#2. SDC2  segment2: 96493701 and 96493935: "cg08979737" "cg25070637" "cg14538332" "cg16935295"
#3  SDC2     Seg 3: 96494006 and 96494335: "cg14538332" "cg16935295"

#TFPI2 segment1: 93890630 and 93890991:"cg24531255" "cg17338208" "cg14775114" "cg16934178" "cg26739865" "cg22441533" "cg14377593"
#TFPI2 segment2: 93889918 and 93890319: "cg19784477" "cg12973591" "cg22799321"

rm(list=ls())
library(data.table)
library(pROC)
library(randomForest)
library(caret)


SDCProbes <- c('cg13096260;cg18719750;cg24732574;cg08979737;cg25070637','cg08979737;cg25070637;cg14538332;cg16935295')
SDCProbes <- c(SDCProbes,'cg14538332;cg16935295')
TFPI2Probes <- c('cg24531255;cg17338208;cg26739865;cg22441533;cg14377593')
TFPI2Probes <- c(TFPI2Probes,'cg12973591;cg22799321')
probes <- unique(as.character(unlist(sapply(c(SDCProbes,TFPI2Probes),function(v){strsplit(v,";")[[1]]}))))




COAD <- fread('F:/projects/allData/TCGA/COAD_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
READ <- fread('F:/projects/allData/TCGA/READ_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
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

clinic <- read.table('F:/projects/allData/TCGA/clinical.tsv',sep="\t",header=T,stringsAsFactors=F)
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
#########################################################################################
getROC <- function(M0,iter=5){
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
   return(ROC)
}

###############################################

normalSamples <- rownames(mat450)[grepl('-11A-',rownames(mat450))]

siteROC <- function(site,DMRs){
	sampleIds <- sampleSite$sampleId[sampleSite$site==site]
	sampleIds <- unique(c(sampleIds,normalSamples))
   
	k0 <- data.frame(sapply(DMRs,function(v){
       pbs <- strsplit(v,";")[[1]]
		 as.numeric(apply(mat450[sampleIds,pbs],1,function(v){mean(v,na.rm=T)}))
	}))
	k0$label <- as.factor(mat450[sampleIds,'label'])

   return(getROC(k0))

}


############################################################################
ROCs <- list()
sites <- T2$site
combs <- list(1,2,3,c(1,2),c(1,3),c(2,3))

pdf('E:/tmp/CRC_site_ROC.pdf',width=10,height=15)
for(site in sites ){
   par(mar=c(0.5,0.5,2,0.5))
	layout(matrix(c(1,1,2:7),ncol=2,byrow=T),heights=c(1,4,4,4))
	plot.new()
	text(0.5,0.5,site,cex=2,font=2)
	for(cbs in combs){
      obj <- siteROC(site,SDCProbes[cbs])
		#str1 <- sprintf("%s:%s",site,paste0('DMR',cbs,collapse="+"))
		str1 <- paste0('DMR',cbs,collapse="+");
      plot(obj,print.auc=TRUE,print.thres="best",col='blue',legacy.axes = TRUE,main=str1)
	}
	#mtext(site,,side=3,line=-1.5)
}
dev.off()

pdf('E:/tmp/CRC_site_ROC.pdf',width=10,height=15)
for(rc in ROCs){
  par(mfrow=c(3,2))
   for(obj in rc){
         plot(obj,print.auc=TRUE,print.thres="best",col='blue',legacy.axes = TRUE,print.auc.cex=par("cex"))
	}
}
dev.off()

pdf("F:/tmp/CRCA_ROC.pdf")
for(i in 1:length(tfAUC)) {
   obj <- tfAUC[[i]]
   plot(obj$ROC,print.auc=TRUE,print.thres="best",col='blue',legacy.axes = TRUE,main=sprintf("CRCA TFPI2: %s",tfNames[i]))
}

for(i in 1:length(sdc2AUC)) {
   obj <- sdc2AUC[[i]]
   plot(obj$ROC,print.auc=TRUE,print.thres="best",col='blue',legacy.axes = TRUE,main=sprintf("CRCA SDC2: %s",sdc2Names[i]))
}
dev.off()