rm(list=ls())
library(data.table)
library(randomForest)
library(caret)
library(rmarkdown)
library(knitr)
library(pROC)
library(sqldf)
source('F:/projects/common/islandDMR.R')

#########################
inputMat450    <- 'F:/projects/allData/TCGA/UCEC_450k.txt'     # The matrix file from TCGA 450k
inputProbeInf  <- 'F:/projects/allData/TCGA/probeInf.txt' # probe information file
inputClinic    <- 'F:/projects/allData/TCGA/clinical.tsv' # clinical information all the patients
projectId      <- 'TCGA-UCEC' # TCGA project code,like 'TCGA-BRCA for breast cance\COAD for colon adenom
#####################################################################################################
print('Loading methylation data......')
mat450 <- fread(inputMat450,sep="\t",header=T,stringsAsFactors=F,check.names=F)
mat450 <- as.data.frame(mat450)
rownames(mat450) <- mat450$probeId
mat450 <- mat450[,-1]
mat450 <- data.frame(t(mat450),check.name=F,stringsAsFactors=F)
mat450 <- mat450[grepl('-01A-|-11A-',rownames(mat450)),]
print('Reading probe information ......')
probInf <- fread(inputProbeInf,header=T,stringsAsFactors=F,check.names=F)[,-c(2,7:9)]
probInf <- as.data.frame(probInf)
probInf <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]
mat450  <- mat450[,probInf[,1]]
lbls <- as.character(sapply(rownames(mat450),function(v){a <- strsplit(v,'-')[[1]][4];ifelse('01A'==a,'cancer','normal')}))
probeIds <- colnames(mat450)
mat450$label <- lbls
mat450 <- mat450[,c('label',probeIds)]
#####################################################################################################
print('Reading clinic information ......')

clinic <- read.table(inputClinic,sep="\t",header=T,stringsAsFactors=F,check.names=F)
clinic <- clinic[clinic$project_id==projectId,]

sites <- unique(clinic$site_of_resection_or_biopsy)
siteSample <- c()
for(site in sites){
	siteSample <- rbind(siteSample,c(site,'cancer',0))
	siteSample <- rbind(siteSample,c(site,'normal',0))
}
siteSample <- data.frame(siteSample,stringsAsFactors=F)
colnames(siteSample) <- c('site','sample_type','sample_number')
siteSample$sample_number <- as.numeric(siteSample$sample_number)

tmp1 <- data.frame('sampleId'  = rownames(mat450),'label'=lbls,
                  'patientId' = as.character(sapply(rownames(mat450),function(v){a <- strsplit(v,'-')[[1]][1:3];paste0(a,collapse="-")})))
tmp2 <- clinic[,c('submitter_id','site_of_resection_or_biopsy')]
colnames(tmp2) <- c('patientId','site')

tmp3 <- sqldf('SELECT tmp2.patientId,site,label FROM tmp1,tmp2 where tmp1.patientId=tmp2.patientId')

T1   <- sqldf('SELECT site,label,count(*) sampleNumber FROM tmp3 group by site,label')
T2   <- sqldf('SELECT site,count(*) CT from tmp3 group by site order by count(*) desc')
T2$ord <- 1:dim(T2)[1]
T2   <- sqldf('SELECT A.*,T2.ord FROM siteSample A, T2 where A.site=T2.site order by T2.ord')

sql  <- 'SELECT T2.site,T2.sample_type,T1.sampleNumber FROM T2 left join T1 '
sql  <- paste0(sql,'ON T2.site=T1.site AND T2.sample_type=T1.label order by T2.ord')
T2   <- sqldf(sql)
T2[is.na(T2)] <- 0
colnames(T2) <- c('site','SampleType','SampleNumber')
T2[dim(T2)[1]+1,] <- c('TOTAL','cancer',sum(T2$SampleNumber[T2$SampleType=='cancer']))
T2 <- data.frame(T2,stringsAsFactors=F)
T2$SampleNumber <- as.numeric(T2$SampleNumber)
T2[dim(T2)[1]+1,] <- c('TOTAL','normal',sum(T2$SampleNumber[T2$SampleType=='normal']))
T2 <- data.frame(T2,stringsAsFactors=F)
T2$SampleNumber <- as.numeric(T2$SampleNumber)
Report.tables.siteSample <- T2  #####
rm(list=c('tmp1','tmp2','tmp3','T1','T2','siteSample','site'))

##################################################################################################
patientInf <- data.frame('sampleId'=rownames(mat450),'label'=lbls,stringsAsFactors=F,
'patientId'=as.character(sapply(rownames(mat450),function(v){a <- strsplit(v,"-")[[1]][1:3];paste0(a,collapse="-")})))
patientInf <- sqldf("SELECT A.*,B.site_of_resection_or_biopsy site FROM patientInf A,clinic B where A.submitter_id=B.patientId")


DMR.overall <- getDMR(mat450)
for(site in sites){
    tmp <- mat450[patientInf$sampleId[patient

}