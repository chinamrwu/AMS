rm(list=ls())
setwd('G:/working')
library(parallel)
library(sqldf)

projectName <- 'Test'
matrixFile  <- 'READSDC2NEG_matrix.txt'
probInf <- read.table('F:/projects/TCGA_450K_sample.txt',sep="\t",header=T,stringsAsFactors=F)[,-2]
colnames(probInf)[1] <- "probeId"

df0 <- NULL
if(file.exists(matrixFile)) {
  print(sprintf("Now loading the existing matrix file %s_matrix.txt,this will take some time......",projectName))
  system.time(df0 <- read.table(matrixFile,sep="\t",header=T,stringsAsFactors=F,row.names=1))
}else{
	f1 <- function(v){read.table(paste0("data/",v),sep="\t",header=T,stringsAsFactors=F,check.names=F)[,2]}
	fnames <- list.files("data")
	fnames <- fnames[grepl('hg38.txt',fnames)]
	print(sprintf("There are %d files to load,it will take some minutes,please waiting.....",length(fnames)))
	system.time(df0 <- data.frame(mclapply(fnames,f1)))
	colnames(df0) <- sapply(fnames,function(v){a=strsplit(v,"\\.")[[1]][6]})
	rownames(df0) <- probInf$probeId
	print(sprintf("Now writing the matrix to local file with name %s_matrix.txt ......",projectName))
	write.table(df0,file=sprintf('data/%s_matrix.txt',projectName),col.names=T,row.names=T,quote=F,sep="\t")
}
probInf <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]

df0 <- df0[probInf$probeId,]
R0  <- apply(df0,1,function(v){sum(is.na(v))})
df0 <- df0[R0==0,]

df0 <- t(df0)
rownames(df0) <- as.character(sapply(rownames(df0),function(v){a=strsplit(v,"\\.|-")[[1]];paste0(a,collapse="-")}))

normIndex    <- which(grepl('-11A-',rownames(df0)))
cancerIndex  <- which(grepl('-01A-',rownames(df0)))

stopifnot(
    {
           if(length(normIndex)*length(cancerIndex) ==0){
			    print("Sample sizes of both cancer and normal must great than one!")
			  }
           length(normIndex)*length(cancerIndex) > 0

   }
	)
df0 <- df0[c(normIndex,cancerIndex),]
normIndex   <- 1:length(normIndex)
cancerIndex <- (length(normIndex)+1):dim(df0)[1]

probScores <- apply(df0,2,function(v){
      m0 <- mean(v[normIndex])
		m1 <- mean(v[cancerIndex])
      c(m0,m1,m1-m0)
})
probScores <- t(probScores)
colnames(probScores) <- c('betaN','betaC','dltBeta')
probScores <- data.frame(probScores)

probScores$probeId <- rownames(probScores)
probScores <- probScores[,c('probeId','betaN','betaC','dltBeta')]
probScores <- sqldf('SELECT A.*,B.Chromosome,B.Start,B.End,B.CGI_Coordinate CGILocation FROM probScores A,probInf B WHERE A.probeId=B.probeId')


geneSymbol <- function(v){
  v1 <- paste0(probInf$Gene_Symbol[probInf$CGI_Coordinate==v],collapse=";")
  v2 <- unique(strsplit(v1,";")[[1]])
  paste0(v2,collapse="|")
}

probScores$GeneSymbol <- sapply(probScores$CGILocation,geneSymbol)
probScores <- probScores[order(probScores$dltBeta,decreasing=T),]

probInf <- probInf[probInf$probeId %in% colnames(df0),]
isLands <- unique(probInf$CGI_Coordinate)

X1  <- sapply(isLands,function(v){
      pbs <- probInf$probeId[probInf$CGI_Coordinate==v]
		tmp <- c()
		if(length(pbs) >1){ tmp <- as.numeric(apply(df0[,pbs],1,mean))}
		else{ tmp <- df0[,pbs]}
		tmp
})
rownames(X1) <- rownames(df0)

islandScores <- apply(X1,2,function(v){
      m0 <- mean(v[normIndex])
		m1 <- mean(v[cancerIndex])
      c(m0,m1,m1-m0)
})

islandScores <- t(islandScores)
colnames(islandScores) <- c('betaN','betaC','dltBeta')
islandScores <- data.frame(islandScores)
islandScores <- islandScores[order(islandScores$dltBeta,decreasing=T),]
islandScores$Gene_Symbol <- sapply(rownames(islandScores),geneSymbol)

#########################

islandScores$CGI <- rownames(islandScores)
islandScores$importance <- 1:dim(islandScores)[1]

tmp <- sqldf('SELECT A.*,B.Gene_Symbol,B.importance FROM probScores A,islandScores B WHERE A.CGILocation=B.CGI order by B.importance asc,A.Start')
strDate <- paste0(strsplit(as.character(Sys.Date()),"-")[[1]],collapse = "")
write.table(tmp[1:500,],file=sprintf('%s_DMisLands_%s.txt',projectName,strDate),sep="\t",col.names=T,row.names=F,quote=F)
write.table(probScores[1:5000,],file=sprintf('%s_DMProbes_%s.txt',projectName,strDate),sep="\t",col.names=T,row.names=F,quote=F)

############ HeatMap plot #####
if(F){
	library(pheatmap);

	topK <- 50
	X2 <- X1[,rownames(islandScores)[1:topK]]
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
}