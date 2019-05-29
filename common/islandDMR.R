library(pROC)
library(data.table)

### find DMR in a CpG island
##  M0: data.frame, whose colnames are probes and contains a column named 'label' indicates a sample belong to cancer or not
##  Note: the probes in the columns of M0 must be orderred by their coordinates in chromosoal


getBlocksOfGene <- function(geneSymbol,Mat){
	siteInf <- probInf[ which(sapply(probInf$Gene_Symbol,function(v){ geneSymbol %in% strsplit(v,";")[[1]]})),]
	siteInf <- siteInf[order(siteInf$Start),]

	isLands <- unique(siteInf$CGI_Coordinate)
	M0 <- Mat[,siteInf$probeId]
	R0 <- apply(M0,2,function(v){sum(is.na(v))})
   blocks <- list() 
   L <- dim(Mat)[1]

	for(obj in isLands){
      probes <- siteInf$probeId[siteInf$CGI_Coordinate==obj]
		R1 <- R0[probes]
		index <- R1 < 0.05*L
		block <- c()
		for(i in 1:length(index)){
         if(index[i]){ block <- c(block,i)}
			else{
             if(length(block) >=2) { blocks[[length(blocks)+1]] <- names(R1)[block];}  
              block <- c()
			}
		}
		if(length(block) >=2){blocks[[length(blocks)+1]] <- names(R1)[block]}
	}
	blocks
}

getBlocksOfisLand <- function(island,M450){
   siteInf <- probInf[probInf$CGI_Coordinate==island,]
	siteInf <- siteInf[order(siteInf$Start),]
	if(dim(siteInf)[1] < 2 ){ return(NULL)}
	M0 <- M450[,siteInf$probeId]
	R0 <- apply(M0,2,function(v){sum(is.na(v))})
	blocks <- list();
   L <- dim(M0)[1]
   
	probes <- siteInf$probeId
	R1 <- R0[probes]
	index <- R1 < 0.05*L
	block <- c()
	for(i in 1:length(index)){
		if(index[i]){ block <- c(block,i)}
		else{
			 if(length(block) >=2) { blocks[[length(blocks)+1]] <- names(R1)[block];}  
			  block <- c()
		}
	}
	if(length(block) >=2){blocks[[length(blocks)+1]] <- names(R1)[block]}
	blocks
}

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

blockDltBeta <- function(block,M0){
	 tmp <- M0[,block]
	 normalIndex <- which(M0$label =='normal')
	 cancerIndex <- which(M0$label =='cancer')
	 k0 <- apply(tmp,1,mean,na.rm=T)
	 dlt <- mean(k0[cancerIndex])-mean(k0[normalIndex])
	 dlt
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
				 
				 sensitivity   <- mean(objROC$sensitivities) ### Need modification here
				 specificity <- mean(objROC$specificities)
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
	result <- result[result$dltBeta >0.1,]
	result <- result[order(result$AUC,decreasing=T),]
	return(result)

}
#############################
getDMROfisLand <- function(island,Mat450){
	chr   <- probInf$Chromosome[probInf$CGI_Coordinate==island][1]
	geneSymbol <- paste0(probInf$Gene_Symbol[probInf$CGI_Coordinate==island],collapse=";")
	geneSymbol <- paste0(unique(strsplit(geneSymbol,";")[[1]]),collapse=";")

	blocks <- getBlocksOfisLand(island,Mat450)
	if(length(blocks)==0){ return(NULL)}

   clss <- tolower(levels(as.factor(Mat450$label)))
	controlName  <- ifelse('normal' %in% clss,'normal',ifelse('cancer' %in% clss,clss[clss!='cancer'],NULL))
	normalIndex <- which(Mat450$label==controlName)
	cancerIndex <- which(Mat450$label!=controlName)
   Mat <- imputeMatrix(Mat450[,c('label',unlist(blocks))]) 

   result <- c()
	for(features in blocks){
	   L=length(features)
		M0 <- Mat[,c('label',features)]
	   for(W in 2:L){
	      for(indx in 1:(L-W+1)){
				 obj <- data.frame('label'=tolower(M0$label),'value'=apply(M0[,features[indx:(indx+W-1)]],1,mean,na.rm=T))
				 objROC <- roc(controls=obj$value[obj$label==controlName],cases=obj$value[obj$label!=controlName])
				 strFeatures <- paste0(features[indx:(indx+W-1)],collapse=";")
				 avgSensitivity   <- mean(objROC$sensitivities)
				 avgSpecificities <- mean(objROC$specificities)
				 Start <- probInf$Start[probInf$probeId==features[indx]]
				 End   <- probInf$End[probInf$probeId==features[indx+W-1]]
				 betaN <- mean(obj$value[normalIndex],na.rm=T)
				 betaC <- mean(obj$value[cancerIndex],na.rm=T)
				 dltBeta <- betaC - betaN
				 result <- rbind(result, sprintf('%s#%4.3f#%4.3f#%4.3f#%s#%d#%d#%d#%4.3f#%4.3f#%4.3f#%s',geneSymbol,betaN,betaC,dltBeta,chr,Start,End,abs(Start-End),avgSensitivity,avgSpecificities,objROC$auc,strFeatures))
		   }
	   }
	   result
	}
		k1 <- data.frame(t(apply(result,1,function(v){ a <- strsplit(v,'#')[[1]]})),stringsAsFactors=F)
		colnames(k1) <- c('geneSymbol','betaN','betaC','dltBeta','chr','Start','End','Width','avgSensitivities','avgSpecificities','auc','probeIds')
		k1$betaN <- as.numeric(k1$betaN)
		k1$betaC <- as.numeric(k1$betaC)
		k1$dltBata <- as.numeric(k1$dltBeta)
		k1$Width <- as.numeric(k1$Width)
		k1$Start <- as.numeric(k1$Start)
		k1$End   <- as.numeric(k1$End)
		k1$avgSensitivities <- as.numeric(k1$avgSensitivities)
		k1$avgSpecificities <- as.numeric(k1$avgSpecificities)
		k1$auc <- as.numeric(k1$auc)
		k1 <- k1[k1$Width >=200,]
		k1 <- k1[k1$dltBeta >0.1,]
		k1 <- k1[order(k1$auc,decreasing=T),]
		return(k1)
}

getDMROfGene <- function(geneSymbol,Mat450){
  genes <- as.logical(sapply(probInf$Gene_Symbol,function(v){ a <- strsplit(v,";")[[1]];geneSymbol %in% a}))
  isLands <- unique(probInf$CGI_Coordinate[genes])
  result <- c()
  for(obj in isLands){result <- rbind(result,getDMROfisLand(obj,Mat450))}
  
  if(! is.null(result)){
     result <- result[order(result$auc,decreasing=T),]
  }
  result
}

getDMR <- function(M450){
	 normalNumber <- sum(M450$label=='normal')
	 cancerNumber <- sum(M450$label=='cancer')
	 if(normalNumber < 5){ print(sprintf("%d normal samples are too fewer ,the statistical power will weak",normalNumber)) }
	 print(sprintf(" %d normal samples and %d cancer samples will be analysed",normalNumber,cancerNumber))
	 isLands <- unique(probInf$CGI_Coordinate)
    tmp <- sapply(isLands,getBlocksOfisLand,M450=M450)
	 blocks <- unlist(tmp,recursive=FALSE)
    dltBeta <- as.numeric(sapply(blocks,blockDltBeta,M0=M450))
	 blocks <- blocks[which(dltBeta>0.20)]
	 print(sprintf("There are %d CpG islands will be analyzed",length(blocks)))
	 k0 <- c()
	 for(obj in blocks){
       k0 <- rbind(k0,getDMROfBlock(obj,M450))
	 }
	 k0 <- k0[order(k0$AUC,decreasing=T),]
	 k0
}
findComplements <- function(matDMR){
}