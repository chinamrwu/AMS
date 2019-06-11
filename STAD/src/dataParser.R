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
tmp <- sapply(tmp,function(v){


})