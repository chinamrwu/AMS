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



clinic <- read.table('F:/projects/TCGA/data/clinical.tsv',sep="\t",header=T,stringsAsFactors=F)
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

		ROC <- roc(ifelse(predicts$observed=="Cancer", "Cancer", "Normal"), as.numeric(predicts$Cancer))
		#pdf(sprintf('output/CRCA_%s_ROC.pdf',nm))
		plot.roc(ROC,print.auc=T,col = "blue3",ylim=c(0,1), print.thres="no",	
		  main=sprintf('ROC for %s in CRCA',nm),legacy.axes =TRUE,print.auc.cex=1.2)
		#dev.off()
	}
   #############




}
