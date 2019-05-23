rm(list=ls())
setwd('G:/projects/COAD')
probInf <- read.table('G:/projects/TCGA_450K_sample.txt',sep="\t",header=T,stringsAsFactors=F)[,-2]
probInf <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]
colnames(probInf)[1] <- "probId"
df0 <- read.table('data/COAD_methy_mat_20190408.txt',sep='\t',header=T,stringsAsFactors=F,check.names=F)
df0 <- df0[probInf$probId,]
df0 <- t(df0)
R0  <- apply(df0,2,function(v){sum(is.na(v))})
df0 <- df0[,R0==0]
spt <- sapply(rownames(df0),function(v){strsplit(v,"-")[[1]][4]})
indx1 <- which(spt=='11A')
indx2 <- which(spt=='01A')
L     <- length(indx1)
df0 <- df0[c(indx1,indx2),]

probInf <- probInf[probInf$probId %in% colnames(df0),]
isLands <- unique(probInf$CGI_Coordinate)

X1  <- sapply(isLands,function(v){
      pbs <- probInf$probId[probInf$CGI_Coordinate==v]
		tmp <- c()
		if(length(pbs) >1){ tmp <- as.numeric(apply(df0[,pbs],1,mean))}
		else{ tmp <- df0[,pbs]}
		tmp
})
rownames(X1) <- rownames(df0)

scores <- apply(X1,2,function(v){
      v0 <- v[1:L]
		v1 <- v[(L+1):dim(X1)[1]]
		#mean(v1)/(mean(v0)+ sd(v1)/mean(v1))
		(mean(v1)-mean(v0))/(mean(v0)^2*10)
})
scores <- sort(scores,decreasing = T)
cancerBeta <- apply( X1[(L+1):dim(X1)[1],],2,mean)
normBeta <- apply( X1[1:L,],2,mean)

scores <- rbind(scores,cancerBeta[names(scores)])
scores <- rbind(scores,normBeta[colnames(scores)])
scores <- data.frame(t(scores))
colnames(scores) <- c('score','cancerBeta','normBetal')

geneSymbol <- function(v){
  v1 <- paste0(probInf$Gene_Symbol[probInf$CGI_Coordinate==v],collapse=";")
  v2 <- unique(strsplit(v1,";")[[1]])
  paste0(v2,collapse="|")
}
scores$Gene_Symbol <- sapply(rownames(scores),geneSymbol)
write.table(scores,"COAD_DMR_20190412.txt",sep="\t",col.names=T,row.names=T,quote=F)


############ HeatMap plot #####
library(pheatmap);

topK <- 50
X2 <- X1[,rownames(scores)[1:topK]]
spLabels <- data.frame('sampleType'=c(rep('Normal',L),rep('Cancer',dim(X1)[1]-L)),stringsAsFactors=F)
rownames(spLabels) <- rownames(X2)

pheatmap(X2,
        color= c('royalblue4','royalblue3','royalblue2',
	        'ivory2','ivory1','ivory',
		'red4','red3','red2','red1'),
	breaks=seq(min(X2),max(X2),by=0.1),
	show_colnames= F,
	show_rownames=F,
	annotation_row = spLabels,
	filename =sprintf('COAD_Top%d_DMR.pdf',topK), 
	width = 9, 
	height = 12,
	main=sprintf("Top %d DMR for colon cancer",topK)
)
################### ROC for each candidate islands
library(randomForest)
library(pROC)