# 回答卢老师关于宫颈癌的问题
#
rm(list=ls())
library(data.table)
library(readr)
library(pROC)
library(sqldf)
setwd("F:/projects/CESC")
options('width'=500)
options('digits'=5)
source('F:/projects/TCGA/src/islandDMR.R')
rm(list=ls())
library(data.table)
library(readr)
library(sqldf)
setwd('F:/projects/CESC')
options('width'=500)
options('digits'=5)

mat1 <- fread('GEO/01/matrix01.txt',sep="\t",header = T,stringsAsFactors = F)
mat1 <- as.data.frame(mat1)
rownames(mat1) <- mat1[,1]

mat2 <- fread("GEO/02/matrix02.txt",sep="\t",header = T,stringsAsFactors = F)
mat2 <- as.data.frame(mat2)
rownames(mat2) <- mat2[,1]

#mat3 <- fread("F:/projects/TCGA/data/CESC_450k.txt",sep="\t",header = T,stringsAsFactors = F)
#mat3 <- as.data.frame(mat3)
#rownames(mat3) <- mat3$probeId

probInf              <- fread(file='F:/projects/TCGA/data/TCGA_450k_sample.txt',sep="\t",header=T,check.names=F,stringsAsFactors=F)[,-2];
probInf              <- as.data.frame(probInf)
colnames(probInf)[1] <- "probeId"
rownames(probInf)    <- probInf$probeId
probInf              <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]
print('450K probe information loaded!')

PAX1 <- probInf[grepl('PAX1',probInf$Gene_Symbol),-c(5:8)]
PAX1 <- PAX1[order(PAX1$Start),]
SOX1 <- sqldf("SELECT * FROM probInf where Gene_Symbol ='SOX1' order by Start")[,-c(5:8)]
probes <- c(PAX1$probeId,SOX1$probeId)

int1 <- intersect(mat1$ID_REF,probes)
int2 <- intersect(mat2$ID_REF,probes)
ints <- intersect(int1,int2)

PAX1 <- PAX1[PAX1$probeId %in% ints,]
SOX1 <- SOX1[SOX1$probeId %in% ints,]

mat1 <- mat1[ints,-1]
mat2 <- mat2[ints,-1]
mat <- data.frame(rbind(t(mat1),t(mat2)))



R0  <- apply(mat,2,function(v){sum(is.na(v))})
mat <- mat[,R0==0]
clnames <- colnames(mat)

sampleInf1 <- read.table('GEO/01/sampleInf.txt',sep='\t',header=T,stringsAsFactors=F)
sampleInf1$label <- tolower(sampleInf1$label)

sampleInf2 <- read.table('GEO/02/sampleInf.txt',sep='\t',header=T,stringsAsFactors=F)
sampleInf2$label <- tolower(sampleInf2$label)

sampleInf <- rbind(sampleInf1[,c(1,2)],sampleInf2[,c(1,2)])
colnames(sampleInf) <- c('sampleId','label')
rownames(sampleInf) <- sampleInf$sampleId
sampleInf$label[sampleInf$label=='scc'] <- 'cancer'
mat$label <- sampleInf[rownames(mat),'label']


########## find DMR 
isLand1 <- sapply(unique(PAX1$CGI_Coordinate),function(v){PAX1$probeId[PAX1$CGI_Coordinate==v]})
isLand2 <- sapply(unique(SOX1$CGI_Coordinate),function(v){SOX1$probeId[SOX1$CGI_Coordinate==v]})


geneDMR <- function(islandList,labelPairs){
      isLands <- ifelse(is.list(islandList),islandList,ifelse(is.matrix(islandList),list(islandList[,1]),list()))
		result <- c()
		for(obj in isLands){
			tmp <- mat[mat$label %in% labelPairs,c('label',obj)]
			result <- rbind(result,bestIslandDMR(tmp))
		}
		result <- result[order(result$auc,decreasing=T),]
		#result <- strsplit(result$probeIds[1],";")[[1]]
		return(result[1,])
}

DMR_PAX1_normalCin3   <- geneDMR(isLand1,c('normal','cin3'))
DMR_PAX1_normalCancer <- geneDMR(isLand1,c('normal','cancer'))
DMR_PAX1_cancerCin3   <- geneDMR(isLand1,c('cancer','cin3'))

DMR_SOX1_normalCin3   <- geneDMR(isLand2,c('normal','cin3'))
DMR_SOX1_normalCancer <- geneDMR(isLand2,c('normal','cancer'))
DMR_SOX1_cancerCin3   <- geneDMR(isLand2,c('cancer','cin3'))


################# 2)不同等级宫颈癌与非癌甲基化水平比较（给出PAX1和SOX1平均水平数据和差异显著性统计数据)

cancerSamples <- sampleInf$sampleId[sampleInf$label=='cancer']
cin3Samples   <- sampleInf$sampleId[sampleInf$label=='cin3']
normalSample  <- sampleInf$sampleId[sampleInf$label=='normal']

cin3Mean1   <- apply(mat[cin3Samples,  strsplit(DMR_PAX1_normalCin3$probeIds,";")[[1]]],1,mean,na.rm=T);
normalMean1 <- apply(mat[normalSample, strsplit(DMR_PAX1_normalCin3$probeIds,";")[[1]]],1,mean,na.rm=T);

t.test(x=normalMean1,y=cin3Mean1)

cancerMean1   <- apply(mat[cancerSamples,  strsplit(DMR_PAX1_normalCancer$probeIds,";")[[1]]],1,mean,na.rm=T);
normalMean1 <- apply(mat[normalSample, strsplit(DMR_PAX1_normalCancer$probeIds,";")[[1]]],1,mean,na.rm=T);
t.test(x=normalMean1,y=cancerMean1)

cin3Mean2   <- apply(mat[cin3Samples,  strsplit(DMR_SOX1_normalCin3$probeIds,";")[[1]]],1,mean,na.rm=T);
normalMean2 <- apply(mat[normalSample, strsplit(DMR_SOX1_normalCin3$probeIds,";")[[1]]],1,mean,na.rm=T);
t.test(x=normalMean2,y=cin3Mean2)

cancerMean2   <- apply(mat[cancerSamples,strsplit(DMR_SOX1_normalCancer$probeIds,";")[[1]]],1,mean,na.rm=T);
normalMean2  <- apply(mat[normalSample, strsplit(DMR_SOX1_normalCancer$probeIds,";")[[1]]],1,mean,na.rm=T);
t.test(x=normalMean2,y=cancerMean2)

tmp <- rbind(tmp,data.frame('methylation'=c(normalMean2,cin3Mean2,cancerMean2),
'sampleType'=c(rep('normal',length(normalMean2)),
					rep('cin3',length(cin3Mean2)),
					rep('cancer',length(cancerMean2))),
'GeneSymbol'=rep('SOX1',length(normalMean2)+length(cin3Mean2)+length(cancerMean2))				
)

p1= ggplot(tmp, aes(sampleType, methylation, color = GeneSymbol,fill=sampleType))+geom_boxplot()+
stat_boxplot(geom ='errorbar')+scale_x_discrete(labels=c('cancer','cin3','normal'))+
	theme(
	    axis.text.x = element_text(size=18,color="black"),
	    axis.text.y = element_text(size=18,color="black"),
		 axis.line.x = element_line(color="black", size = 0.25),
	    axis.line.y = element_line(color="black", size = 0.25),
	    panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
		 panel.background = element_blank()
	)
## 3.	不同等级宫颈癌与非癌转录水平比较（给出平均水平数据和差异显著性统计数据）
tcgaCancer <- read.table('F:/tmp/GTEx/figshare/cesc-rsem-fpkm-tcga-t.txt',sep="\t",header=T,stringsAsFactors=F,check.names=F)
s1 <- as.character(sapply(colnames(tcgaCancer),function(v){strsplit(v,"-")[[1]][4]}))
index <- which(grepl('-01A-',colnames(tcgaCancer)))
tcgaCancer <- tcgaCancer[,c(1,2,index)]
tcgaNormal <- read.table('F:/tmp/GTEx/figshare/cesc-rsem-fpkm-tcga.txt',sep="\t",header=T,stringsAsFactors=F,check.names=F)
gtex   <- read.table('F:/tmp/GTEx/figshare/cervix-rsem-fpkm-gtex.txt',sep="\t",header=T,stringsAsFactors=F,check.names=F)

normal <- gt

tumor <- tcgaCancer[tcgaCancer[,1] %in% c('PAX1','SOX1'),-2]
rownames(tumor) <- tumor[,1]
tumor <- tumor[,-1]
tumor <- data.frame(t(tumor))

normal <- cbind(tcgaCancer[tcgaCancer[,1] %in% c('PAX1','SOX1'),-2],gtex[gtex[,1] %in% c('PAX1','SOX1'),-c(1,2)])
rownames(normal) <- normal[,1]
normal <- normal[,-1]
normal <- data.frame(t(normal))
##4. PAX1与SOX1在高、低甲基化样本的生存分析
# 见CESC_Survival.R

# 5. PAX1和SOX1基因代谢通路分析，输出代谢通路图。略

# 6. PAX1和SOX1基因以及双位点甲基化检测的对宫颈癌ROC曲线图分析，输出ROC曲线图及其相关数据。
  
cancerSamples <- sampleInf$sampleId[sampleInf$label=='cancer']
cin3Samples   <- sampleInf$sampleId[sampleInf$label=='cin3']
normalSample  <- sampleInf$sampleId[sampleInf$label=='normal']

cin3Mean1   <- apply(mat[cin3Samples,  strsplit(DMR_PAX1_normalCin3$probeIds,";")[[1]]],1,mean,na.rm=T);
normalMean1 <- apply(mat[normalSample, strsplit(DMR_PAX1_normalCin3$probeIds,";")[[1]]],1,mean,na.rm=T);
ROC <- roc('controls'=normalMean1,'cases'=cin3Mean1)
p1=plot(ROC,main="PAX1 in CESC :normal vs cin3",col='blue',print.thres='best',print.auc=T)

cancerMean1 <- apply(mat[cancerSamples,  strsplit(DMR_PAX1_normalCancer$probeIds,";")[[1]]],1,mean,na.rm=T);
normalMean1 <- apply(mat[normalSample, strsplit(DMR_PAX1_normalCancer$probeIds,";")[[1]]],1,mean,na.rm=T);
ROC <- roc('controls'=normalMean1,'cases'=cancerMean1)
p2=plot(ROC,main="PAX1 in CESC :normal vs cancer",col='blue',print.thres='best',print.auc=T)

cancerMean1 <- apply(mat[cancerSamples,  strsplit(DMR_PAX1_cancerCin3$probeIds,";")[[1]]],1,mean,na.rm=T);
cin3Mean1   <- apply(mat[cin3Samples,    strsplit(DMR_PAX1_cancerCin3$probeIds,";")[[1]]],1,mean,na.rm=T);
ROC <- roc('controls'=cin3Mean1,'cases'=cancerMean1)
p3=plot(ROC,main="PAX1 in CESC :cin3 vs cancer",col='blue',print.thres='best',print.auc=T)

############################
cin3Mean1   <- apply(mat[cin3Samples,  strsplit(DMR_SOX1_normalCin3$probeIds,";")[[1]]],1,mean,na.rm=T);
normalMean1 <- apply(mat[normalSample, strsplit(DMR_SOX1_normalCin3$probeIds,";")[[1]]],1,mean,na.rm=T);
ROC <- roc('controls'=normalMean1,'cases'=cin3Mean1)
p4 <- plot(ROC,main="SOX1 in CESC :normal vs cin3",col='blue',print.thres='best',print.auc=T)

cancerMean1 <- apply(mat[cancerSamples,  strsplit(DMR_SOX1_normalCancer$probeIds,";")[[1]]],1,mean,na.rm=T);
normalMean1 <- apply(mat[normalSample, strsplit(DMR_SOX1_normalCancer$probeIds,";")[[1]]],1,mean,na.rm=T);
ROC <- roc('controls'=normalMean1,'cases'=cancerMean1)
p5 <- plot(ROC,main="SOX1 in CESC :normal vs cancer",col='blue',print.thres='best',print.auc=T)

cancerMean1 <- apply(mat[cancerSamples,  strsplit(DMR_SOX1_cancerCin3$probeIds,";")[[1]]],1,mean,na.rm=T);
cin3Mean1   <- apply(mat[cin3Samples,    strsplit(DMR_SOX1_cancerCin3$probeIds,";")[[1]]],1,mean,na.rm=T);
ROC <- roc('controls'=cin3Mean1,'cases'=cancerMean1)
p6 <- plot(ROC,main="SOX1 in CESC :cin3 vs cancer",col='blue',print.thres='best',print.auc=T)