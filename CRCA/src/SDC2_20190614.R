rm(list=ls())
library(data.table)
library(pROC)
library(sqldf)

setwd("F:/projects/CRC")
COAD <- fread('F:/projects/allData/TCGA/COAD_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
READ <- fread('F:/projects/allData/TCGA/READ_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
mat450 <- cbind(as.data.frame(COAD),as.data.frame(READ[,-1]));

rm(COAD)
rm(READ)
rownames(mat450) <- mat450$probeId
mat450 <- mat450[,-1]
mat450 <- data.frame(t(mat450),stringsAsFactors=F,check.names=F)

mat450 <- mat450[grepl('-01A-|-11A-',rownames(mat450)),] #438 samples:393 cancer vs 45 normal
Label <- as.character(sapply(rownames(mat450),function(v){a <- strsplit(v,"-")[[1]][4];ifelse('01A'==a,'cancer','normal')}))
mat450$label <- Label

probInf <- fread('F:/projects/allData/TCGA/TCGA_450K_sample.txt',sep="\t",header=T,stringsAsFactors=F)[,-2]
probInf <- as.data.frame(probInf)
colnames(probInf)[1] <- "probeId"
probInf <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]

# chr8: 96493520-96495379
SDC2 <- sqldf("SELECT probeId FROM probInf where Chromosome='chr8' and Start between 96493519 and 96495380 order by Start")[,1]

#######################################
imputeMatrix <- function(M0){
      clss <- tolower(levels(as.factor(M0$label)))
		controlName  <- ifelse('normal' %in% clss,'normal',ifelse('cancer' %in% clss,clss[clss!='cancer'],NULL))
		normalIndex <- which(M0$label==controlName)
		cancerIndex <- which(M0$label!=controlName)
		
		R0 <- apply(M0,2,function(v){sum(is.na(v))})
	   M0[,R0!=0] <- sapply(names(R0)[R0!=0],function(v){
		      k0 <- M0[,v]
				v0 <- k0[normalIndex]
				v1 <- k0[cancerIndex]
				v0[is.na(v0)] <- mean(v0,na.rm=T)
				v1[is.na(v1)] <- mean(v1,na.rm=T)
				k0[normalIndex] <- v0
				k0[cancerIndex] <- v1
				k0
		})
		M0
}

getDMROfBlock <- function(block,M450){
   chr   <- probInf$Chromosome[probInf$probeId==block[1]][1]
	geneSymbol <- paste0(probInf$Gene_Symbol[probInf$probeId %in% block],collapse=";")
	geneSymbol <- paste0(unique(strsplit(geneSymbol,";")[[1]]),collapse=";")

	

   clss <- tolower(levels(as.factor(M450$label)))
	controlName  <- ifelse('normal' %in% clss,'normal',ifelse('cancer' %in% clss,clss[clss!='cancer'],NULL))
	normalIndex <- which(M450$label==controlName)
	cancerIndex <- which(M450$label!=controlName)
   Mat <- imputeMatrix(M450[,c('label',block)]) 
   result <- data.frame(matrix(nrow=0,ncol=12));
	colnames(result) <- c('geneSymbol','chr','Start','End','Width','betaN','betaC','dltBeta','senesitivity','specificity','AUC','probeIds')
	
	L=length(block)
	M0 <- Mat[,c('label',block)]
	for(W in 2:L){
	    for(indx in 1:(L-W+1)){
				 obj <- data.frame('label'=tolower(M0$label),'value'=apply(M0[,block[indx:(indx+W-1)]],1,mean,na.rm=T))
				 objROC <- roc(controls=obj$value[obj$label==controlName],cases=obj$value[obj$label!=controlName])
				 strFeatures <- paste0(block[indx:(indx+W-1)],collapse=";")
				 
				 #sensitivity   <- mean(objROC$sensitivities) ### Need modification here
				 #specificity <- mean(objROC$specificities)
            
				 x <- objROC$sensitivities+objROC$specificities
				 mxIndex <- which(x==max(x))
				 
				 #x <- abs(objROC$thresholds - 0.5)
             #mxIndex <- which(x==min(x))
				 
				 sensitivity <- objROC$sensitivities[mxIndex]
				 specificity <- objROC$specificities[mxIndex]
             #####################
				 Start <- probInf$Start[probInf$probeId==block[indx]]
				 End   <- probInf$End[probInf$probeId==block[indx+W-1]]
				 betaN <- mean(obj$value[normalIndex],na.rm=T)
				 betaC <- mean(obj$value[cancerIndex],na.rm=T)
				 dltBeta <- betaC - betaN
             
				 tdf <- data.frame(matrix(nrow=0,ncol=12));
				 colnames(tdf) <- c('geneSymbol','chr','Start','End','Width','betaN','betaC','dltBeta','sensitivity','specificity','AUC','probeIds')
             result[nrow(result)+1,] <- list(geneSymbol,chr,Start,End,abs(End-Start),betaN,betaC,dltBeta,
				 sensitivity,specificity,objROC$auc,strFeatures)
		 }
	 }
	
	result <- result[result$Width >=200,]
	result <- result[result$dltBeta > 0.2,]
	result <- result[result$betaN  < 0.2,]
	result <- result[order(result$dltBeta,decreasing=T),]
	return(result)

}
####################################
mat <- imputeMatrix(mat450[,c('label',SDC2)])
DMR <- getDMROfBlock(SDC2,mat)
DMR <- DMR[order(DMR$dltBeta,decreasing=T),]
DMR <- DMR[!is.na(DMR$probeIds),]

pdf('output/ROC_0614.pdf',width=12,height=12)
for(i in 1:52){
  probes <- strsplit(DMR$probeIds[i],";")[[1]]
  bta <- as.numeric(apply(mat[,probes],1,mean))
  objROC <- roc('controls'=bta[mat$label=='normal'],'cases'=bta[mat$label=='cancer'])
  if( (i-1) %% 4 ==0) {
       par(mar=c(0.5,0.5,2,0.5))
	    layout(matrix(1:4,ncol=2,byrow=T),heights=c(5,5,5,5))
  }
  	plot(objROC,print.auc=TRUE,print.thres="best",col='blue',legacy.axes = TRUE,main=sprintf('SDC2_DMR_%d',i),
	print.auc.cex=2.0,print.thres.cex=1.5)
}
