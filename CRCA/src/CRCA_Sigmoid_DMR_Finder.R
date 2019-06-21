#
# searching DMR for Sigmoid colon
# author: wuzhicheng
# date: 20190506
#

rm(list=ls())
library(sqldf)
library(data.table)
library(ggplot2)
setwd('F:/projects/CRC')

COAD <- fread('F:/projects/TCGA/data/COAD_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
READ <- fread('F:/projects/TCGA/data/READ_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
mat450 <- cbind(as.data.frame(COAD),as.data.frame(READ[,-1]));

rm(COAD)
rm(READ)
rownames(mat450) <- mat450$probeId
mat450 <- mat450[,-1]
mat450 <- data.frame(t(mat450),stringsAsFactors=F,check.names=F)


clinic <- read.table('F:/projects/TCGA/data/clinical.tsv',sep="\t",header=T,stringsAsFactors=F,check.names=F)#Squamous sample information
print('loading data completed')
patientId  <- as.character(sapply(rownames(mat450),function(v){a=strsplit(v,"-")[[1]][1:3];paste0(a,collapse="-")}))
clinic     <- clinic[clinic$submitter_id %in% patientId,]
Sigmoid    <- clinic$submitter_id[clinic$site_of_resection_or_biopsy=='Rectosigmoid junction']#'Sigmoid colon']

index      <- match(Sigmoid,patientId)
index      <- unique(c(index,which(grepl('-11A-',rownames(mat450)))))## Sigmoid + Normal

df0 <- mat450[index,]

cancerIndex <- which(grepl('-01A-',rownames(df0))) ## only Squamous
normIndex   <- which(grepl('-11A-',rownames(df0)))
df1 <- df0[c(normIndex,cancerIndex),]
normIndex   <- 1:length(normIndex)
cancerIndex <- (length(normIndex)+1):dim(df1)[1]
df0 <- df1
rm(df1)

probInf <- read.table('F:/projects/TCGA/data/TCGA_450K_sample.txt',sep="\t",header=T,stringsAsFactors=F)[,-2]
colnames(probInf)[1] <- "probeId"
probInf <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]



df0 <- df0[,probInf$probeId]
R0  <- apply(df0,2,function(v){sum(is.na(v))})
df0 <- df0[,R0==0]
probInf <- probInf[probInf$probeId %in% colnames(df0),]
rownames(probInf) <- probInf$probeId
###去除重复样本###
rnames <- as.character(sapply(rownames(df0),function(v){a=strsplit(v,"-|\\.")[[1]][1:7];paste0(a,collapse="-")}))
if(length(rnames)!=length(unique(rnames))){  
    index <- which(!duplicated(rnames));
    df0 <- df0[index,] 
	 rownames(df0) <- rnames[index]
}

df0 <- rbind(df0[which(grepl('-11A-',rownames(df0))),],df0[which(grepl('-01A-',rownames(df0))),])
normIndex <- which(grepl('-11A-',rownames(df0)))
cancerIndex <- which(grepl('-01A-',rownames(df0)))
k0 <- apply(df0[normIndex,],2,mean)
k1 <- apply(df0[cancerIndex,],2,mean)

####################################################


geneSymbolByIsland <- function(v){
    v0  <- probInf$Gene_Symbol[probInf$CGI_Coordinate==v]
    v1 <- unique(strsplit(paste0(v0,collapse=";"),";")[[1]])
	 v1 <- v1[!grepl("\\.",v1)]
	 paste0(v1,collapse=";")
}

isLands <- unique(probInf$CGI_Coordinate)
crs <- sapply(isLands,function(v){
   geneSym <- geneSymbolByIsland(v)

   prbs <- probInf[probInf$CGI_Coordinate==v,];
	prbs <- prbs[order(prbs$Start),]
	prbs <- prbs$probeId
   d0 <- k0[prbs]
	d1 <- abs(d0[2:length(d0)] - d0[1:(length(d0)-1)])
   Regions <- c()
	index <- which(d1 >= 0.2)
	if(length(index)>0){
			index <- c(1,index+1)
			if(index[length(index)] != length(d0)){ index <- c(index,length(d0)) }
			for(i in 2:length(index)){
           w <- index[i]-index[i-1]
			  if(w>=3){
			     probs <- paste0(prbs[index[i-1]:(index[i]-1)],collapse=";")
			     Width <- as.character(probInf[prbs[index[i]-1],'End']-probInf[prbs[index[i-1]],'Start'])
              Regions <- c(Regions,paste0(c(probs,Width,geneSym),collapse="#"))
			  }
			}
	}else if(length(prbs)>=2) { ## all the probes belongs to an whole continous region
	   
	    probs <- paste0(prbs,collapse=";")
		 Width <- as.character(probInf[prbs[length(prbs)],'End']- probInf[prbs[1],'Start'])
       Regions <- c(Regions,paste0(c(probs,Width,geneSym),collapse="#"))
	}
   Regions
})

DMRScores <- c()
for(obj0 in crs){
  for(obj1 in obj0){
   a <- strsplit(obj1,'#')[[1]]
	prbs <- strsplit(a[1],';')[[1]]
	m0 <- apply(df0[normIndex,prbs],2,mean)
	m1 <- apply(df0[cancerIndex,prbs],2,mean)
	DMRScores <- rbind(DMRScores,c(mean(m0),mean(m1),obj1))
	}
}

DMRScores <- data.frame(DMRScores,stringsAsFactors=F)
colnames(DMRScores) <- c('betaN','betaC','probeInf')
DMRScores$betaN <- as.numeric(DMRScores$betaN)
DMRScores$betaC <- as.numeric(DMRScores$betaC)
kk <- t(sapply(DMRScores$probeInf,function(v){a <- strsplit(v,'#')[[1]];c(a[3],a[2],a[1])}))
DMR <- cbind(DMRScores[,c(1,2)],kk)
colnames(DMR) <- c('betaN','betaC','GeneSymbol','BaseNumber','probeIds')
DMR$dltBeta  <- DMR$betaC-DMR$betaN
DMR <- DMR[, c('betaN','betaC','dltBeta','GeneSymbol','BaseNumber','probeIds')]
DMR <- DMR[order(DMR$dltBeta,decreasing=T),]
rownames(DMR) <- 1:dim(DMR)[1]

selection <- DMR[DMR$betaN <0.2 & DMR$betaC>0.5,]

locations <- apply(selection,1,function(v){
  
})

write.table(selection,file=sprintf("DMR_RectosigmoidJunction_20190507.txt",dim(selection)[1]),row.names=F,quote=F)

############################## boxplot

if(F){
   topK <- 5
	probs <- c()
   for(i in 1:50){
     probs <- c(probs,strsplit(DMR$probeIds[i],';')[[1]])
	}
	k0 <- df0[,probs]
	Labels  <- as.character(sapply(rownames(k0),function(v){a=strsplit(v,"-")[[1]][4];ifelse('11A'==a,'Normal','Cancer')}))
	
	k1 <- data.frame('Beta'=k0[,1],'probeId'=rep(probs[1],length(Labels)),'subtype'=Labels)
	for(i in 2:length(probs)){ k1 <- rbind(k1,data.frame('Beta'=k0[,i],'probeId'=rep(probs[i],length(Labels)),'subtype'=Labels))}
   
	XLabels <- as.character(apply(probInf[probs,],1,function(v){
	genes <- unique(strsplit(as.character(v[5]),";")[[1]]);
	if(length(genes)>=1){
	   genes <- genes[!grepl("\\.|-",genes)];
	   if(sum(grepl('PCDHGA',genes))>0){genes <- genes[1]}
		genes <- genes[!grepl("-",genes)]

	 }
	 paste0(v[1]," ",paste0(genes,collapse="|"))
	 }))
	color2 <- c('Normal'='blue','Cancer'='red')
	p1 <- ggplot(k1, aes(probeId, Beta, color = subtype,fill=subtype))+geom_boxplot()+ggtitle('鳞状食管癌前17候选DMR探针(92癌vs19癌旁)')+
	#labs(title="前50候选探针")+
	scale_x_discrete(labels=XLabels)+
	theme(axis.text.x = element_text(angle = 90, hjust = 1),
	    panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
		 plot.title = element_text(hjust = 0.5),
	    #panel.border = element_blank(),
	    axis.line.x = element_line(color="black", size = 0.25),
	    axis.line.y = element_line(color="black", size = 0.25),
	    #plot.title   = element_text(size=16),
	    panel.background = element_blank())+
	scale_colour_manual(values=color2)
	
	#pdf("
	
	#+stat_smooth(method = "rlm")

######################################### Heatmap
library(pheatmap);
library(RColorBrewer)

X <- DMR[DMR$betaN <0.2 & DMR$betaC>0.5 & !is.na(DMR$GeneSymbol),]
X1 <- apply(X,1,function(v){
   kk <- apply(df0[,strsplit(as.character(v[6]),";")[[1]]],1,mean);
})
colnames(X1) <- X$GeneSymbol
colnames(X1)[grepl('PCDHG',colnames(X1))] <- 'PCDHGA|BX'

spLabels <- data.frame('sampleType'=as.character(sapply(rownames(X1),function(v){a <- strsplit(v,"-")[[1]][4];ifelse('11A'==a,"Normal","Cancer")})),stringsAsFactors=F)
rownames(spLabels) <- rownames(X1)

pheatmap(X1,
   #     color= c('royalblue4','royalblue3','royalblue2','ivory2','ivory1','ivory','red4','red3','red2','red1'),
	#color= c('slateblue','salmon','#D70131'),
   #color=c('blue','slateblue1','#D70131'),
	#breaks=seq(min(X1),max(X1),by=(max(X1)-min(X1))/3),
	color =rev(brewer.pal(11,'RdYlBu')),
	#color=brewer.pal(
	show_colnames= T,
	show_rownames=F,
	annotation_row = spLabels,
	border_color=NA,
	filename =sprintf('ESAD_Top%d_DMR_10.pdf',dim(X1)[2]), 
	width = 9, 
	height = 12,
	main=sprintf("Top %d DMR for ESAD_1",dim(X1)[2])
)

}
