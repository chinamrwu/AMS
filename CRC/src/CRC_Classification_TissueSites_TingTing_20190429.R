rm(list=ls())
library(randomForest)
library(glmnet)
library(data.table)
library(readr)

library(pROC)

setwd('F:/projects/TCGA')
#setwd("/home/zcwu/")
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
probInf              <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.',]# & probInf$Feature_Type=='Island',]
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



manifest <- 'F:/projects/TCGA/data/gdc_all450K_manifest.2019-04-03.txt' # manifest file that contains the id of files 
mf       <- read.table(manifest,sep="\t",header=T,stringsAsFactors=F)
rnames   <- rownames(X)
sampleInf <- data.frame('sampleId'=rnames,'label'=X$label,
                        'submiter'=as.character(sapply(rnames,function(v){paste0(strsplit(v,"-")[[1]][1:3],collapse="-")})),
								stringsAsFactors=F)
library(sqldf)
tmp <- sqldf("SELECT A.*,B.sampleId FROM clinic A,sampleInf B where A.submitter_id=B.submiter")
if(F){ # below codes to draw ROC plot using random forest

   library(caret)
	features <- list('SCD2'='rg149','TFPI2'='rg196','SCD2+TFPI2'=c('rg149','rg196'))
   clinic <- read.table('F:/projects/TCGA/data/clinical.tsv',sep="\t",header=T,stringsAsFactors=F)
   
	submitters <- as.character(sapply(rownames(X),function(v){paste0(strsplit(v,"-")[[1]][1:3],collapse="-")}))
	clinic <- clinic[clinic$submitter_id %in% submitters,c('submitter_id','site_of_resection_or_biopsy','tissue_or_organ_of_origin')] 
   # site_of_resection_or_biopsy  tissue_or_organ_of_origin
	preds <- list();
   
   for(nm in names(features)){
		X1 <- X[,c('label',features[[nm]])]
      set.seed(1)
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
      predicts$submitter <- submitters

		t0 <- sqldf('SELECT A.*,B.site_of_resection_or_biopsy,B.tissue_or_organ_of_origin FROM predicts A,clinic B WHERE A.submitter=B.submitter_id')
      rownames(t0) <- rownames(predicts)
		preds[[length(preds)+1]] <- t0
	 }
	 names(preds) <- names(features)
    
	 sites1 <- unique(clinic$site_of_resection_or_biopsy)
	 sites2 <- unique(clinic$tissue_or_organ_of_origin)
    
	 siteOf <- function(mPred,site){
      tmp <- mPred[mPred$site_of_resection_or_biopsy==site,]
	   pIndx <- which(tmp$observed=='Cancer')
      nIndx <- which(tmp$observed=='Normal')
		TP <- sum(tmp$predicted[pIndx]=='Cancer')
		TN <- sum(tmp$predicted[nIndx]=='Normal')
		c(dim(tmp)[1],length(pIndx),TP,TP/length(pIndx)*100,length(nIndx),TN,TN/length(nIndx)*100)
	 }
     
	 obj1 <- c()
	 for(site in sites1){
     obj1 <- cbind(obj1,sapply(preds,siteOf,site=site))
	 }
    obj1 <- t(obj1)
	 colnames(obj1) <- c("Total","Cancer","predictedCancer","Sensitivity","Normal","predictedNormal","Speficity")
	 obj1 <- data.frame(obj1)
	  L <- dim(obj1)[1]
	 index <- seq(1,L,by=3)
	 obj1 <- obj1[c(index,index+1,index+2),]
	 site <- as.character(sapply(sites1,function(v){paste0(strsplit(v," |,|;")[[1]],collapse="_")}));
	 obj1$site <- c(site,site,site)
	 obj1$GeneSymbol <- c(rep('SDC2',length(sites1)),rep('TFPI2',length(sites1)),rep('SDC2+TFPI2',length(sites1)))
    obj1 <- obj1[,c('GeneSymbol','site','Total','Cancer',"predictedCancer","Sensitivity","Normal","predictedNormal","Speficity")]
   
	  tissueOf <- function(mPred,site){
      tmp <- mPred[mPred$tissue_or_organ_of_origin==site,]
	   pIndx <- which(tmp$observed=='Cancer')
      nIndx <- which(tmp$observed=='Normal')
		TP <- sum(tmp$predicted[pIndx]=='Cancer')
		TN <- sum(tmp$predicted[nIndx]=='Normal')
		c(dim(tmp)[1],length(pIndx),TP,TP/length(pIndx)*100,length(nIndx),TN,TN/length(nIndx)*100)
	 }

    obj2 <- c()
	 for(site in sites2){
     obj2 <- cbind(obj2,sapply(preds,tissueOf,site=site))
	 }
    obj2 <- t(obj2)
	 colnames(obj2) <- c("Total","Cancer","predictedCancer","Sensitivity","Normal","predictedNormal","Speficity")
	 obj2 <- data.frame(obj2)
	 L <- dim(obj2)[1]
	 index <- seq(1,L,by=3)
	 obj2 <- obj2[c(index,index+1,index+2),]
    tissue <- as.character(sapply(sites2,function(v){paste0(strsplit(v," |,|;")[[1]],collapse="_")}))
	 obj2$tissue <- c(tissue,tissue,tissue)
	 obj2$GeneSymbol <- c(rep('SDC2',length(sites2)),rep('TFPI2',length(sites2)),rep('SDC2+TFPI2',length(sites2)))
    obj2 <- obj2[,c('GeneSymbol','tissue','Total','Cancer',"predictedCancer","Sensitivity","Normal","predictedNormal","Speficity")]
  
  write.table(obj1,file="site_of_resection_or_biopsy.txt",sep="\t",quote=F,col.names=T,row.names=F)
  write.table(obj2,file="tissue_or_organ_of_origin.txt",sep="\t",quote=F,col.names=T,row.names=F)




}
   #############
