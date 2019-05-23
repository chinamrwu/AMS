#ESAD
rm(list=ls())
library(sqldf)
library(data.table)
library(ggplot2)
setwd('F:/projects/ESAD')
df0 <- fread(file='data/ESAD_probe_matrix_20190417.txt',sep="\t",header=T,check.names=F,stringsAsFactors=F,nThread=2)
df0 <- data.frame(df0,check.names=F,stringsAsFactors=F)
clinic <- read.table('data/clinical.tsv',sep="\t",header=T,stringsAsFactors=F,check.names=F)#Squamous sample information
print('loading data completed')
patientId <- clinic$submitter_id
probes <- df0$probeId
samples <- as.character(sapply(colnames(df0)[-1],function(v){a=strsplit(v,"-")[[1]][1:3];paste0(a,collapse="-")}))
index <- match(patientId,samples)+1
cancers <- t(df0[,index]) ## only Squamous
nrms    <- t(df0[,which(grepl('-11A-',colnames(df0)))])
df1 <- rbind(nrms,cancers)
normIndex   <- 1:dim(nrms)[1]
cancerIndex <- (dim(nrms)[1]+1):dim(df1)[1]
df0 <- df1
rm(df1)
rm(cancers)
rm(nrms)
colnames(df0) <- probes
df0 <- data.frame(t(df0),check.names=F)

projectName <- 'ESAD'
probInf <- read.table('F:/projects/TCGA_450K_sample.txt',sep="\t",header=T,stringsAsFactors=F)[,-2]
colnames(probInf)[1] <- "probeId"
probInf <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]



df0 <- df0[probInf$probeId,]
R0  <- apply(df0,1,function(v){sum(is.na(v))})
df0 <- df0[R0==0,]
df0 <- data.frame(t(df0))
probInf <- probInf[probInf$probeId %in% colnames(df0),]

print("Now assess each probe ......")
probScores <- apply(df0,2,function(v){
      m0 <- mean(v[normIndex])
		m1 <- mean(v[cancerIndex])
		sd0 <- sd(v[normIndex])
		sd1 <- sd(v[cancerIndex])

		m2 <- (m1-m0)/(m0 * (0.5+sd0)*(0.5+sd1))
      c(m0,m1,m1-m0,m2)
})
probScores <- t(probScores)
colnames(probScores) <- c('betaN','betaC','dltBeta','score')
probScores <- data.frame(probScores)
probScores$probeId <- rownames(probScores)
probScores <- probScores[,c('probeId','betaN','betaC','dltBeta','score')]
probScores <- sqldf('SELECT A.*,B.Chromosome,B.Start,B.End,B.CGI_Coordinate CGILocation FROM probScores A,probInf B WHERE A.probeId=B.probeId')
geneSymbol <- function(v){
  v1 <- paste0(probInf$Gene_Symbol[probInf$CGI_Coordinate==v],collapse=";")
  v2 <- unique(strsplit(v1,";")[[1]])
  paste0(v2,collapse="|")
}

probScores$GeneSymbol <- sapply(probScores$CGILocation,geneSymbol)
probScores <- probScores[order(probScores$dltBeta,decreasing=T),]


######################################
probInf <- probInf[probInf$probeId %in% colnames(df0),]
isLands <- unique(probInf$CGI_Coordinate)
print("assessing each island......")
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

if(F){
	k0 <- df0[,match(tmp$probeId[1:50],colnames(df0))]
	Labels  <- as.character(sapply(rownames(k0),function(v){a=strsplit(v,"\\.")[[1]][4];ifelse('11A'==a,'Normal','Cancer')}))
	probs <- colnames(k0)
	
	k1 <- data.frame('Beta'=k0[,1],'probeId'=rep(probs[1],length(Labels)),'subtype'=Labels)
	for(i in 2:length(probs)){ k1 <- rbind(k1,data.frame('Beta'=k0[,i],'probeId'=rep(probs[i],length(Labels)),'subtype'=Labels))}
    
	XLabels <- as.character(apply(tmp[1:50,],1,function(v){paste0(v[c(10,1)],collapse=" ")}))
	color2 <- c('Normal'='blue','Cancer'='red')
	p1 <- ggplot(k1, aes(probeId, Beta, color = subtype,fill=subtype))+geom_boxplot()+ggtitle('ÁÛ×´Ê³¹Ü°©DMRÇ°50ºòÑ¡Ì½Õë(92°©vs19°©ÅÔ)')+
	#labs(title="Ç°50ºòÑ¡Ì½Õë")+
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
	
	
	#+stat_smooth(method = "rlm")
}