# 回应李婷婷 《结直肠癌投稿论文需要的生信分析结果》
# 1）所有SDC2低甲基化样本（不区分部位）与癌旁样本DMR结果；
# 2)直肠(Rectum) SDC2低甲基化样本与癌旁样本DMR结果；
# 3)结肠(Colon) SDC2低甲基化样本与癌旁样本DMR结果；
# 4)乙状结肠(sigmoid colon) SDC2低甲基化样本与癌旁样本DMR结果；
# 5)降结肠(Descending colon) SDC2低甲基化样本与癌旁样本DMR结果；
# 6)直乙交界(Rectosigmoid junction) SDC2低甲基化样本与癌旁样本DMR结果；
#
#
remove(list=ls())
library(data.table)
library(sqldf)
library(pROC)

setwd('F:/projects/TCGA');
source('F:/projects/TCGA/src/findDMR.R')
COAD <- fread('F:/projects/TCGA/data/COAD_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
COAD <- as.data.frame(COAD)
rownames(COAD) <- COAD$probeId
COAD <- COAD[,-1]

READ <- fread('F:/projects/TCGA/data/READ_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
READ <- as.data.frame(READ)
rownames(READ) <- READ$probeId
READ <- READ[,-1]

matALL <- cbind(COAD,READ);


clinic <- read.table('data/clinical.tsv',sep="\t",header = T,stringsAsFactors = F)
clinic <- clinic[clinic$project_id %in% c('TCGA-READ','TCGA-COAD'),]
probInf <- read.table('F:/projects/TCGA/data/TCGA_450K_sample.txt',sep="\t",header=T,stringsAsFactors=F)[,-2]
colnames(probInf)[1] <- "probeId"
probInf <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]
DMR <- read.table('F:/projects/CRC/output/CRCA_DMR_0424.txt',sep="\t",header=T,stringsAsFactors=F)
SDC2Probes <- strsplit(DMR$probeIds[grepl('SDC2',DMR$GeneSymbol)][1],";")[[1]]

normalSamples <- colnames(matALL)[grepl('-11A-',colnames(matALL))]
allPersons <- as.character(sapply(colnames(matALL),function(v){a <- strsplit(v,"-")[[1]][1:3];paste0(a,collapse="-")}))

SDC2 <- apply(matALL[SDC2Probes,],2,mean,na.rm=T)

############################################################################################

getSDC2LowSamples <- function(PersonIds =NULL){
    rtv <- NULL
	 if('all'==PersonIds){  # 1
       rtv <- names(SDC2)[SDC2 < 0.2]
	 }else if(PersonIds=='COAD'){ #2 
      m0 <- apply(COAD[SDC2Probes,],2,mean,na.rm=T)
		rtv  <- names(m0)[m0 < 0.2]
	 } else if(PersonIds=='READ'){ # 3
      m0 <- apply(READ[SDC2Probes,],2,mean,na.rm=T)
		rtv  <- names(m0)[m0 < 0.2] # 4
	 }else    if('Sigmoid colon'==PersonIds){
       submiters <- clinic$submitter_id[clinic$site_of_resection_or_biopsy =='Sigmoid colon']
		 index <- match(submiters,allPersons)
       m0   <- apply(matALL[SDC2Probes,index],2,mean,na.rm=T)
		 rtv  <- names(m0)[m0 < 0.2]
	 }else	 if('Rectosigmoid junction'==PersonIds){
       submiters <- clinic$submitter_id[clinic$site_of_resection_or_biopsy =='Rectosigmoid junction']
		 index <- match(submiters,allPersons)
       m0   <- apply(matALL[SDC2Probes,index],2,mean,na.rm=T)
		 rtv  <- names(m0)[m0 < 0.2]
	 }else if('Descending colon'==PersonIds){
       submiters <- clinic$submitter_id[clinic$site_of_resection_or_biopsy =='Descending colon']
		 index <- match(submiters,allPersons)
       m0   <- apply(matALL[SDC2Probes,index],2,mean,na.rm=T)
		 rtv  <- names(m0)[m0 < 0.25]
	 }
	 rtv
}
#################################################################

workNames <- c('all','COAD','READ','Sigmoid colon','Rectosigmoid junction')
probPath='F:/projects/TCGA/data/TCGA_450K_sample.txt'
for( nm in workNames){
   nm <- 'Descending colon'
   sampleIds <- unique(c(normalSamples,getSDC2LowSamples(nm)))
	print(sprintf("%s:%d samples",nm,length(sampleIds)))
	mat <- matALL[,sampleIds]
	DMR <- findDMR(mat,probPath)
   write.table(DMR,file=sprintf("CRC_%s_DMR_20190509.txt",nm),col.names=T,row.names=F,quote=F)
}
