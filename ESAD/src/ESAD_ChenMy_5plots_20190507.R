rm(list=ls())
library(sqldf)
library(data.table)
library(ggplot2)
setwd('F:/projects/TCGA')
df0 <- fread(file='data/ESCA_450k.txt',sep="\t",header=T,check.names=F,stringsAsFactors=F,nThread=2)
df0 <- data.frame(df0,check.names=F,stringsAsFactors=F)
clinic <- read.table('F:/projects/ESAD/data/clinical.tsv',sep="\t",header=T,stringsAsFactors=F,check.names=F)#Squamous sample information
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

probInf <- read.table('data/TCGA_450K_sample.txt',sep="\t",header=T,stringsAsFactors=F)[,-2]
colnames(probInf)[1] <- "probeId"
probInf <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]



df0 <- df0[probInf$probeId,]
R0  <- apply(df0,1,function(v){sum(is.na(v))})
df0 <- df0[R0==0,unique(colnames(df0))]
df0 <- data.frame(t(df0),row.names=colnames(df0),check.names=F)
probInf <- probInf[probInf$probeId %in% colnames(df0),]
rownames(probInf) <- probInf$probeId
##########################
DMR <- read.csv('F:/projects/ESAD/ESAD_DMR_20190425.csv',header=T,stringsAsFactors=F)
DMR <- DMR[!is.na(DMR$GeneSymbol),]
DMR$GeneSymbol <- as.character(sapply(DMR$GeneSymbol,function(v){strsplit(v,";")[[1]][1]}))
######################################################## 1  heatmap for probes and samples
library(pheatmap)

probes <- as.character(unlist(sapply(DMR$probeIds,function(v){strsplit(v,";")[[1]]})));
tmp <- t(df0[,probes])
tmp <- data.frame(tmp,check.names=F)

Label <- as.character(sapply(colnames(tmp),function(v){a <- strsplit(v,"-")[[1]][4];ifelse('01A'==a,'cancer','normal')}));
colAnnotation <- data.frame("SampleType"=Label)
rownames(colAnnotation) <- colnames(tmp)

 pheatmap(as.matrix(tmp),
          show_colnames= F,
          color =rev(brewer.pal(11,'RdYlBu')),
          show_rownames=T,
          annotation_col = colAnnotation,
			 border_color = NA,
			filename = "../ESAD/output/ESAD_Heatmap_probes_0507.png",width =6,height =8) 
###################################################################################### 3 heatmap for DMR

normIndex   <- which(grepl('-11A-',rownames(df0)))
cancerIndex <- which(grepl('-01A-',rownames(df0)))

X <- sapply(DMR$probeIds,function(v){
      probes <- strsplit(v,";")[[1]]
		v0 <- apply(df0[,probes],1,mean)
      v0
})
colnames(X) <- DMR$GeneSymbol

Label <- as.character(sapply(rownames(X),function(v){a <- strsplit(v,"-")[[1]][4];ifelse('01A'==a,'cancer','normal')}));
sampleAnnotation <- data.frame("SampleType"=Label)
rownames(sampleAnnotation) <- rownames(X)

 pheatmap(as.matrix(X),
          color =rev(brewer.pal(11,'RdYlBu')),
          show_rownames=F,
			 show_colnames= T,
          annotation_row = sampleAnnotation,
			 border_color = NA,
			 main="DMR for ESAD",
			filename = "../ESAD/output/ESAD_Heatmap_DMR_0507.png",width =6,height =8) 
##################################################################################### 4.boxplot for DMR genes
 genes <- colnames(X)
   X1 <- c()
	for(gene in genes){
     X1 <- rbind(X1,data.frame('Beta'=X[,gene],'GeneSymbol'=rep(gene,dim(X)[1]),
	  'sampleType'=Label,
	  stringsAsFactors=F))
	}

	color2 <- c('normal'='blue','cancer'='red')

	p1 <- ggplot(X1, aes(GeneSymbol, Beta, color = sampleType,fill=sampleType))+geom_boxplot()+
	ggtitle('ÁÛ×´Ê³¹Ü°©DMR(91°©vs17°©ÅÔ)')+
	stat_boxplot(geom ='errorbar')+
	scale_x_discrete(labels=genes)+
	theme(
	    axis.text.x = element_text(angle = 90, hjust = 1,size=18,color="black"),
	    axis.text.y = element_text(size=18,color="black"),
		 axis.line.x = element_line(color="black", size = 0.25),
	    axis.line.y = element_line(color="black", size = 0.25),
	    panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
		 #panel.border = element_blank(),
		 panel.background = element_blank())+
		 plot.title = element_text(hjust = 0.5,size=18),
	    scale_colour_manual(values=c('normal'='black','cancer'='black'))
################################################################## 5 Í¬4





################################################################## violinPlot for DMR
X2            <-  X1
X2$GeneSymbol <-  as.character(apply(X2,1,function(v){paste0(v[2:3],collapse="/")}))
pdf(file="F:/projects/ESAD/output/ESAD_Vioplot_16DMR.pdf",width=25,height=8)
p2 <- ggplot(X2, aes(x=GeneSymbol, y=Beta,fill=sampleType)) + 
  geom_violin()+
stat_summary(fun.y=mean, geom="point", shape=23, size=2,color='red')+
theme(axis.text.x = element_text(angle = 90, hjust = 1),
	    panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
		 plot.title = element_text(hjust = 0.5),
	    #panel.border = element_blank(),
	    axis.line.x = element_line(color="black", size = 0.25),
	    axis.line.y = element_line(color="black", size = 0.25),
	    #plot.title   = element_text(size=16),
	    panel.background = element_blank())
print(p2)
dev.off()