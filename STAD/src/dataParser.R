rm(list=ls())
library(data.table)
library(sqldf)

setwd('F:/projects/STAD')
matInf <- read.table('data/GSE99553_series_matrix.txt',sep="\n",nrows=72,stringsAsFactors=F)
accs <- matInf[grepl('!Sample_geo_accession',matInf[,1]),]
accs <- strsplit(accs,"\t")[[1]][-1]

tmp <- matInf[grepl('!Sample_characteristics_ch1',matInf[,1]),]
status <- strsplit(tmp[2],"\t")[[1]][-1]
status <- tolower(trimws(as.character(sapply(status,function(v){strsplit(v,":")[[1]][2]}))))

infection <-  strsplit(tmp[3],"\t")[[1]][-1]
infection <-  tolower(trimws(as.character(sapply(infection,function(v){strsplit(v,":")[[1]][2]}))))
GSE99553_sampleInf <- data.frame('acc'=accs,'status'=status,'infection'=infection,stringsAsFactors=F)
write.table(GSE99553_sampleInf,file='data/GSE99553_sampleInf.txt',sep="\t",col.names=T,row.names=F,quote=F)
mat0 <- read.table('data/GSE99553_series_matrix.txt',sep="\n",header=T,stringsAsFactors=F)
indx1 <- which(grepl('!series_matrix_table_begin',mat0[,1]))+1
indx2 <- which(grepl('!series_matrix_table_end',mat0[,1]))-1
mat0 <- mat0[indx1:indx2,]
write.table(mat0,file='data/matGSE99553.txt',sep="\n",col.names=F,row.names=F,quote=F)
###################################
matInf <- read.table('data/GSE103186_series_matrix.txt',sep="\n",stringsAsFactors=F)
accs <- matInf[grepl('!Sample_geo_accession',matInf[,1]),]
accs <- strsplit(accs,"\t")[[1]][-1]

tmp <- matInf[grepl('!Sample_characteristics_ch1',matInf[,1]),]
tmp <- strsplit(tmp,"\t")[[1]][-1]
tmp <- as.character(sapply(tmp,function(v){strsplit(v,": ")[[1]][2]}))
tmp1 <- as.character(sapply(tmp,function(v){
  a<-   ifelse(v=='intestinal metaplasia biopsy from gastric antrum',      paste0('IM','_antrum'),
	ifelse(v=='intestinal metaplasia biopsy from gastric body',        paste0('IM','_body'),
	ifelse(v=='intestinal metaplasia biopsy from gastric cardia',      paste0('IM','_cardia'),
	ifelse(v=='mild intestinal metaplasia biopsy from gastric antrum', paste0('MIM','_antrum'),
	ifelse(v=='normal biopsy from gastric antrum',                     paste0('normal','_antrum'),
	ifelse(v=='normal biopsy from gastric body',                       paste0('normal','_body'),
	ifelse(v=='normal biopsy from gastric cardia',                     paste0('normal','_cardia'),NA)))))))
	a
}))
tmp1 <- sapply(tmp1,function(v){strsplit(v,"_")[[1]]})
tmp1 <- data.frame(t(tmp1),stringsAsFactors=F)
rownames(tmp1) <- 1:dim(tmp1)[1]
colnames(tmp1) <- c('label','site')
tmp1$acc <- accs
sampleInf <- tmp1[,c('acc','label','site')]
write.table(sampleInf,file='data/GSE103186_sampleInf.txt',sep="\t",col.names=T,row.names=F,quote=F)

indx1 <- which(grepl('!series_matrix_table_begin',matInf[,1]))+1
indx2 <- which(grepl('!series_matrix_table_end',matInf[,1]))-1
mat0  <- matInf[indx1:indx2,]
write.table(mat0,file='data/matGSE103186.txt',sep="\n",col.names=F,row.names=F,quote=F)

#########################################
tmp <- read.table('data/GSE25869_series_matrix.txt',sep="\n",header=F,stringsAsFactors=F)
acc <- strsplit(tmp[grepl('!Sample_geo_accession',tmp[,1]),1],"\t")[[1]][-1]
src <- strsplit(tmp[grepl('!Sample_source_name_ch1',tmp[,1]),1],"\t")[[1]][-1]
sampleCharacters <- tmp[grepl('!Sample_characteristics_ch1',tmp[,1]),1]
tissueType <- strsplit(sampleCharacters[2],"\t")[[1]][2:65]
tissueType <- as.character(sapply(tissueType,function(v){strsplit(v,": gastric ")[[1]][2]}))

status <- strsplit(sampleCharacters[3],"\t")[[1]][2:65]
status <- as.character(sapply(status,function(v){strsplit(v,": ")[[1]][2]}))
indx1 <- which(grepl('!series_matrix_table_begin',tmp[,1]))+1
indx2 <- which(grepl('!series_matrix_table_end',tmp[,1]))-1
mat0  <- tmp[indx1:indx2,]
write.table(mat0,file='data/matGSE25869.txt',sep="\n",col.names=F,row.names=F,quote=F)
mat0 <- read.table('data/matGSE25869.txt',sep="\t",header=T,stringsAsFactors=F)
rownames(mat0) <- mat0$ID_REF
mat0 <- mat0[,-1]
mat0 <- data.frame(t(mat0))
mat0 <- mat0[1:64,]
probes <- colnames(mat0)
mat0$label <- status
mat0 <- mat0[,c('label',probes)]
write.table(mat0,file='data/matGSE25869.txt',sep="\t",col.names=F,row.names=F,quote=F)

##############################################################################################

tmp <- read.table('data/GSE30601_series_matrix.txt',sep="\n",header=F,stringsAsFactors=F)
acc <- strsplit(tmp[grepl('!Sample_geo_accession',tmp[,1]),1],"\t")[[1]][-1]
sampleCharacters <- tmp[grepl('!Sample_characteristics_ch1',tmp[,1]),1]

status <- strsplit(sampleCharacters[2],"\t")[[1]][-1]
status <- as.character(sapply(status,function(v){a <- strsplit(v,": ")[[1]][2];ifelse('normal'==a,'normal','cancer')}))

indx1 <- which(grepl('!series_matrix_table_begin',tmp[,1]))+1
indx2 <- which(grepl('!series_matrix_table_end',tmp[,1]))-1
mat0  <- tmp[indx1:indx2,]
write.table(mat0,file='data/matGSE30601.txt',sep="\n",col.names=F,row.names=F,quote=F)
mat0 <- read.table('data/matGSE30601.txt',sep="\t",header=T,stringsAsFactors=F)
rownames(mat0) <- mat0$ID_REF
mat0 <- mat0[,-1]
mat0 <- data.frame(t(mat0))
probes <- colnames(mat0)
mat0$label <- status
mat0 <- mat0[,c('label',probes)]
write.table(mat0,file='data/matGSE30601.txt',sep="\t",col.names=F,row.names=F,quote=F)

