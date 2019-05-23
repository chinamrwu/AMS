rm(list=ls())
library(data.table)
library(pROC)
library(randomForest)
library(caret)
setwd("F:/projects/CRC")

SDCProbes <- c('cg13096260;cg18719750;cg24732574;cg08979737;cg25070637','cg08979737;cg25070637;cg14538332;cg16935295')
SDCProbes <- c(SDCProbes,'cg14538332;cg16935295')
TFPI2Probes <- c('cg24531255;cg17338208;cg26739865;cg22441533;cg14377593')
TFPI2Probes <- c(TFPI2Probes,'cg12973591;cg22799321')
probes <- unique(as.character(unlist(sapply(c(SDCProbes,TFPI2Probes),function(v){strsplit(v,";")[[1]]}))))




COAD <- fread('F:/projects/TCGA/data/COAD_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
READ <- fread('F:/projects/TCGA/data/READ_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
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


clinic <- read.table('F:/projects/TCGA/data/clinical.tsv',sep="\t",header=T,stringsAsFactors=F)
submitters <- as.character(sapply(rownames(mat450),function(v){paste0(strsplit(v,"-")[[1]][1:3],collapse="-")}))
clinic <- clinic[clinic$submitter_id %in% submitters,c('submitter_id','site_of_resection_or_biopsy','tissue_or_organ_of_origin')] 
rownames(clinic) <- clinic$submitter_id
########################################################################
 
 probes <- unique(strsplit(paste0(SDCProbes,collapse = ";"),";")[[1]])
   X1 <- c()
	for(pb in probes){
     X1 <- rbind(X1,data.frame('Beta'=mat450[,pb],'GeneSymbol'=rep(pb,dim(mat450)[1]),
	  'sampleType'=tolower(Label),  stringsAsFactors=F))
	}

	color2 <- c('normal'='blue','cancer'='red')

	p1 <- ggplot(X1, aes(GeneSymbol, Beta, color = sampleType,fill=sampleType))+geom_boxplot()+
	ggtitle('结直肠癌SDC2基因关键探针')+
	stat_boxplot(geom ='errorbar')+
	scale_x_discrete(labels=probes)+
	theme(
	    axis.text.x = element_text(angle = 90, hjust = 1,size=18,color="black"),
	    axis.text.y = element_text(size=18,color="black"),
		 axis.line.x = element_line(color="black", size = 0.25),
	    axis.line.y = element_line(color="black", size = 0.25),
	    panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
		 #panel.border = element_blank(),
		 panel.background = element_blank(),
		 plot.title = element_text(hjust = 0.5,size=18))+
	    scale_colour_manual(values=c('normal'='black','cancer'='black'))


###########

probes <- unique(strsplit(paste0(TFPI2Probes,collapse = ";"),";")[[1]])
   X1 <- c()
	for(pb in probes){
     X1 <- rbind(X1,data.frame('Beta'=mat450[,pb],'GeneSymbol'=rep(pb,dim(mat450)[1]),
	  'sampleType'=tolower(Label),  stringsAsFactors=F))
	}

	color2 <- c('normal'='blue','cancer'='red')

	p2 <- ggplot(X1, aes(GeneSymbol, Beta, color = sampleType,fill=sampleType))+geom_boxplot()+
	ggtitle('结直肠癌TFPI2基因关键探针')+
	stat_boxplot(geom ='errorbar')+
	scale_x_discrete(labels=probes)+
	theme(
	    axis.text.x = element_text(angle = 90, hjust = 1,size=18,color="black"),
	    axis.text.y = element_text(size=18,color="black"),
		 axis.line.x = element_line(color="black", size = 0.25),
	    axis.line.y = element_line(color="black", size = 0.25),
	    panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
		 #panel.border = element_blank(),
		 panel.background = element_blank(),
		 plot.title = element_text(hjust = 0.5,size=18))+
	    scale_colour_manual(values=c('normal'='black','cancer'='black'))

################################################################## 5 同4





################################################################## violinPlot for DMR
sampleSite <- data.frame('sampleId'=rownames(mat450),stringsAsFactors=F)
sampleSite$site =as.character(sapply(sampleSite$sampleId,function(v){a <- paste0(strsplit(v,"-")[[1]][1:3],collapse="-");
clinic[a,'site_of_resection_or_biopsy']}))
rownames(sampleSite) <- sampleSite[,1]
sampleSite$label <- "NO"
sampleSite$label[grepl('-01A-',rownames(sampleSite))] <- 'cancer'
sampleSite$label[grepl('-11A-',rownames(sampleSite))] <- 'normal'


X1 <- c()
for(i in 1:length(SDCProbes)){
    pb <- strsplit(SDCProbes[i],";")[[1]]
     X1 <- rbind(X1,data.frame('Beta'=as.numeric(apply(mat450[,pb],1,mean,na.rm=F)),
	  'island'=rep(paste0("CpG_island_",i),dim(mat450)[1]),
	  'sampleType'=tolower(Label),
	  'site'=sampleSite[rownames(mat450),2],
	    stringsAsFactors=F
	  ))
}
X1$site[X1$site=="Connective, subcutaneous and other soft tissues of abdomen"] <- 'Other soft tissues'

vioPlots <- list()
for(i in 1:3){
tmp <- ggplot(X1[X1$island==paste0('CpG_island_',i),],aes(x=site, y=Beta, color = sampleType,fill=sampleType))+geom_boxplot()+
	ggtitle(sprintf('结直肠不同区域在SDC2基因CpG-island-%d的甲基化水平',i))+
	stat_boxplot(geom ='errorbar')+
	#scale_x_discrete(labels=probes)+
	theme(
	    axis.text.x = element_text(angle = 90, hjust = 1,size=18,color="black"),
	    axis.text.y = element_text(size=18,color="black"),
		 axis.line.x = element_line(color="black", size = 0.25),
	    axis.line.y = element_line(color="black", size = 0.25),
	    panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
		 #panel.border = element_blank(),
		 panel.background = element_blank(),
		 plot.title = element_text(hjust = 0.5,size=18))+
	    scale_colour_manual(values=c('normal'='black','cancer'='black'))
  vioPlots[[length(vioPlots)+1]] <- tmp
 }
pdf("output/CRC_SDC2_sites_boxplot.pdf",family='GB1',width=12)
 for(i in 1:length(vioPlots)){print(vioPlots[[i]])}
dev.off()

###########################
X1 <- c()
for(i in 1:length(TFPI2Probes)){
    pb <- strsplit(TFPI2Probes[i],";")[[1]]
     X1 <- rbind(X1,data.frame('Beta'=as.numeric(apply(mat450[,pb],1,mean,na.rm=F)),
	  'island'=rep(paste0("CpG_island_",i),dim(mat450)[1]),
	  'sampleType'=tolower(Label),
	  'site'=sampleSite[rownames(mat450),2],
	    stringsAsFactors=F
	  ))
}
X1$site[X1$site=="Connective, subcutaneous and other soft tissues of abdomen"] <- 'Other soft tissues'

vioPlots2 <- list()
for(i in 1:length(TFPI2Probes)){
tmp <- ggplot(X1[X1$island==paste0('CpG_island_',i),],aes(x=site, y=Beta, color = sampleType,fill=sampleType))+geom_boxplot()+
	ggtitle(sprintf('结直肠不同区域在TFPI2基因CpG-island-%d的甲基化水平',i))+
	stat_boxplot(geom ='errorbar')+
	#scale_x_discrete(labels=probes)+
	theme(
	    axis.text.x = element_text(angle = 90, hjust = 1,size=18,color="black"),
	    axis.text.y = element_text(size=18,color="black"),
		 axis.line.x = element_line(color="black", size = 0.25),
	    axis.line.y = element_line(color="black", size = 0.25),
	    panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
		 #panel.border = element_blank(),
		 panel.background = element_blank(),
		 plot.title = element_text(hjust = 0.5,size=18))+
	    scale_colour_manual(values=c('normal'='black','cancer'='black'))
  vioPlots2[[length(vioPlots2)+1]] <- tmp
}

pdf("output/CRC_TFPI2_sites_boxplot.pdf",family='GB1',width=12)
print(vioPlots2[[1]])
print(vioPlots2[[2]])
dev.off()
