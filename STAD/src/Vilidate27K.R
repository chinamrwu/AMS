rm(list=ls())
library(data.table)
library(sqldf)

setwd('F:/projects/STAD')
probInf <- fread('F:/projects/allData/TCGA/probeInf.txt',sep="\t",header=T,stringsAsFactors=F)
probInf <- as.data.frame(probInf)
probInf <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]

DMR <- read.csv('output/STAD_DMR_20190612.csv',header=T,stringsAsFactors=F)

meta27k   <-  read.csv('F:/projects/allData/illumina_humanmethylation27_content.csv',header=T,stringsAsFactors=F)[,c(1,3,4,11,12,20)]
meta450k  <- read.csv('F:/projects/allData/HumanMethylation450_15017482_v1-2.csv',skip=7,header=T,stringsAsFactors=F)[,c(2,12,13,16)]
meta450k  <- meta450k[!duplicated(meta450k$Name),]
rownames(meta450k) <- meta450k[,1]
meta450k$Coordinate_36 <- as.integer(meta450k$Coordinate_36)
meta450k <- meta450k[!is.na(meta450k$Coordinate_36),]
colnames(meta450k) <- c('probeId','chr','start450','start27')
##############################
matGSE30601 <- fread('data/matGSE30601.txt',sep="\t",header=T,stringsAsFactors=F)
matGSE30601 <- as.data.frame(matGSE30601)
matGSE25869 <- fread('data/matGSE25869.txt',sep="\t",header=T,stringsAsFactors=F)
matGSE25869 <- as.data.frame(matGSE25869)
matGSE25869$label[matGSE25869$label!='normal'] <- 'cancer'
########################################
getAreaPlot <- function(mat1){
   features <- colnames(mat1)[colnames(mat1)!='label']
	mat2 <- mat1[c(which(mat1$label=='normal'),which(mat1$label!='normal')),]
   L <- dim(mat2)[1]
	
	dfX <- c()
	for( feature in features){
      df1 <- data.frame('X'=1:L,'Beta'=mat2[,feature],'probeId'=rep(feature,L),'sampleType'=mat2$label,stringsAsFactors=F)
		dfX <- rbind(dfX,df1)
	}
   
	p1 <- ggplot(dfX,aes(x=X,y=Beta,fill=sampleType))+geom_area()+facet_grid(probeId ~ .)+
			ggtitle(paste0(features,collapse=";"))+
			scale_x_discrete(expand = c(0.01,0))+
			theme(legend.position="bottom",
					panel.grid.major = element_blank(),
					panel.grid.minor = element_blank(),
					panel.border = element_blank(),
					axis.line.x = element_line(color="black", size = 0.25),
					axis.line.y = element_line(color="black", size = 0.25),
					plot.title   = element_text(size=16),
					panel.background = element_blank()
			)
	return(p1)
}
#####################################
getProbes27k <- function(probes450k){
   probes <- strsplit(probes450k,";")[[1]];
	chr   <- meta450k[probes[1],'chr']
	Start <- meta450k[probes[1],'start27']-100
	End   <- meta450k[probes[length(probes)],'start27']+100
	meta27k$Name[meta27k$Chr==chr & meta27k$MapInfo >= Start &meta27k$MapInfo <= End]
}

library(ggplot2)

pbs270 <- c()
for(i in 1:500){
   pbs <- getProbes27k(DMR$probeIds[i]);
   if(length(pbs)>0){ pbs270 <- unique(c(pbs270,paste0(pbs,collapse=";")))}
}



plots <- list()
for(pbs in pbs270){
    fs <- strsplit(pbs,";")[[1]]
    plots[[length(plots)+1]] <- getAreaPlot(matGSE25869[,c('label',fs)])
}

