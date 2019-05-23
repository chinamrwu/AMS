rm(list=ls())
library(data.table)
library(readr)
setwd('F:/projects/CESC')


df0 <- read.table('GEO/01/GSE46306_series_matrix.txt',sep='\n',header=F);
geoAcc <- strsplit(as.character(df0[grepl("!Sample_geo_accession",df0[,1]),1]),'\t')[[1]][-1]

index <- which(grepl('!Sample_characteristics_ch1',df0[,1]))
label <- as.character(sapply(strsplit(as.character(df0[index[3],1]),"\t")[[1]][-1],function(v){
         a=strsplit(v,'status:')[[1]][2];
			trimws(a)
			}))
label[is.na(label)] <- 'cin3'

HPV <- as.character(sapply(strsplit(as.character(df0[index[4],1]),"\t")[[1]][-1],function(v){
         a=strsplit(v,'status:')[[1]][2];
			trimws(a)
			}))
sampleInf <- data.frame('sampleId'=geoAcc,'label'=label,'HPV'=HPV,stringsAsFactors=F)
write.table(sampleInf,file='GEO/01/sampleInf.txt',sep='\t',col.names=T,row.names=F,quote=F)



index01 <- which(grepl('!series_matrix_table_begin',df0[,1]))
index02 <- which(grepl('!series_matrix_table_end',  df0[,1]))
write.table(df0[(index01+1):(index02-1),],file='GEO/01/matrix01.txt',sep="\n",col.names=F,row.names=F,quote=F)

####################################### parse datasets for paper2
df0 <- read.table('GEO/02/GSE99511_series_matrix.txt',sep='\n',header=F);
Series    <- df0[grepl('!Series_',df0[,1]),]
sampleTitle <- strsplit(as.character(df0[grepl('!Sample_title',df0[,1]),1]),'\t')[[1]][-1]

index01 <- which(grepl('!series_matrix_table_begin',df0[,1]))
index02 <- which(grepl('!series_matrix_table_end',  df0[,1]))
# write.table(df0[(index01+1):(index02-1),],file='GEO/02/matrix02.txt',sep="\n",col.names=F,row.names=F,quote=F)
sampleNames <- strsplit(as.character(df0[index01+1,1]),'\t')[[1]][-1]
sampleInf <- data.frame('sampleId'=sampleNames,
                       'label'=as.character(sapply(sampleTitle,function(v){strsplit(v,"_")[[1]][1]})),stringsAsFactors=F)
#write.table(sampleInf,file='GEO/02/sampleInf.txt',sep='\t',col.names=T,row.names=F,quote=F)

#####################################################################################################################################
df0 <- fread('GEO/03/GSE68339-GPL13534_series_matrix.txt',sep='\n',header=F);
df0 <- as.data.frame(df0);
geoAcc <- strsplit(as.character(df0[grepl("!Sample_geo_accession",df0[,1]),1]),'\t')[[1]][-1]
geoAcc <- as.character(sapply(geoAcc,function(v){substr(v,2,nchar(v)-1)}))

sampleSrc <-  strsplit(as.character(df0[grepl("!Sample_source_name_ch1",df0[,1]),1]),'\t')[[1]][-1]
sampleSrc <- as.character(sapply(sampleSrc,function(v){a=substr(v,2,nchar(v)-1);strsplit(a," ")[[1]][2]}))

stage <-  strsplit(as.character(df0[grepl("!Sample_characteristics_ch1",df0[,1]),1]),'\t')[[1]][-1]
stage <- as.character(sapply(stage,function(v){a=substr(v,2,nchar(v)-1);b=strsplit(a,"stage:")[[1]][2];trimws(b)}))

sampleInf <- data.frame('geoAcc'=geoAcc,'label'=sampleSrc,'figo_stage'=stage,stringsAsFactors=F)
write.table(sampleInf,file='GEO/03/sampleInf.txt',sep='\t',col.names=T,row.names=F,quote=F)
index01 <- which(grepl('!series_matrix_table_begin',df0[,1]))
index02 <- which(grepl('!series_matrix_table_end',  df0[,1]))
write.table(df0[(index01+1):(index02-1),],file='GEO/03/matrix01.txt',sep="\n",col.names=F,row.names=F,quote=F)

