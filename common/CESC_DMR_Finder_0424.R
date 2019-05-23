rm(list=ls())
library(data.table)
library(readr)
library(sqldf)
options('width'=400)
#搜寻连续CpG sites 构成的DMR区域
#mat450Path
findDMR <- function(mat450,probInfPath,topK=100){
  #mat450 <- as.data.frame(COAD)
  #probInfPath='F:/projects/TCGA_450K_sample.txt'
  if(is.null(mat450)){
     print("450K methylation matrix for all samples are required!")
	  return(NULL)
	}
    rownames(mat450) <- mat450$probeId
    mat450           <- mat450[,-1]

  	cnames <- as.character(sapply(colnames(mat450),function(v){a=strsplit(v,"-|\\.")[[1]][1:7];paste0(a,collapse="-")})) ## remove ducplicated samples
   if(length(cnames)!=length(unique(cnames))){  
		 index <- which(!duplicated(cnames));
		 mat450 <- mat450[,index] 
		 colnames(mat450) <- cnames[index]
		 print(sprintf('%d duplication been removed',length(cnames)-length(unique(cnames)))) 
    }
   normIndex    <- which(grepl('-11A-',colnames(mat450)))
   cancerIndex  <- which(grepl('-01A-',colnames(mat450)))
  
  if(length(normIndex)*length(cancerIndex) ==0) {
	  print("Sample sizes of both cancer and normal must great than one!")
	  return(NULL)
   }
  	print(sprintf('Analysis will be based on %d normal samples and %d cancer samples!',length(normIndex),length(cancerIndex)))
   
	probInf              <- fread(file=probInfPath,sep="\t",header=T,check.names=F,stringsAsFactors=F)[,-2];
	probInf              <- as.data.frame(probInf)
   colnames(probInf)[1] <- "probeId"
   probInf              <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]
   print('450K probe information loaded!')
  
   mat450 <- mat450[probInf$probeId,]
   R0  <- apply(mat450,1,function(v){sum(is.na(v))})
   mat450 <- mat450[R0==0,]
   mat450 <- data.frame(t(mat450))
   probInf <- probInf[probInf$probeId %in% colnames(mat450),]
	rownames(probInf) <- probInf$probeId
   
	geneSymbolByIsland <- function(v){
		 v0  <- probInf$Gene_Symbol[probInf$CGI_Coordinate==v]
		 v1 <- unique(strsplit(paste0(v0,collapse=";"),";")[[1]])
		 v2 <- v1[!grepl("\\.",v1)]
		 v3 <- unique(as.character(sapply(v1[grepl('\\.',v1)],function(s0){strsplit(s0,'\\.')[[1]]})))

		 paste0(c(v2,v3),collapse=";")
	}

	 ####################################
	 k0 <- apply(mat450[normIndex,],  2,mean)
    #k1 <- apply(mat450[cancerIndex,],2,mean)
   
	 ################################################################

	 isLands <- unique(probInf$CGI_Coordinate)
	 print(sprintf('Now searching for continuous CpG sites on each of %d islands',length(isLands)))
	 crs <- sapply(isLands,function(v){
				geneSym <- geneSymbolByIsland(v)

				prbs <- probInf[probInf$CGI_Coordinate==v,];
				prbs <- prbs[order(prbs$Start),]
				prbs <- prbs$probeId
				d0 <- k0[prbs]
				d1 <- abs(d0[2:length(d0)] - d0[1:(length(d0)-1)])
				Regions <- c()
				index <- which(d1 >= 0.1)  ##################################  key cutoff to break continuous CpG sites into regions
				if(length(index)>0){
					index <- c(1,index+1)
					if(index[length(index)] != length(d0)){ index <- c(index,length(d0)) }
					for(i in 2:length(index)){
					  w <- index[i]-index[i-1]
					  if(w>=3){
						  probs <- paste0(prbs[index[i-1]:(index[i]-1)],collapse=";")
						  Width <- as.character(probInf[prbs[index[i]-1],'End']-probInf[prbs[index[i-1]],'Start'])
						  Regions <- c(Regions,paste0(c(probs,Width,geneSym),collapse='#'))
					  }
					}
				}else if(length(prbs)>=2) { ## all the probes belongs to an whole continous region
					 probs <- paste0(prbs,collapse=";")
					 Width <- as.character(probInf[prbs[length(prbs)],'End']- probInf[prbs[1],'Start'])
					 Regions <- c(Regions,paste0(c(probs,Width,geneSym),collapse='#'))
				}
				Regions
      })

print('Assessing each continuous region......')
  DMRScores <- c()
      for(obj0 in crs){
        for(obj1 in obj0){
           a <- strsplit(obj1,'#')[[1]]
	        prbs <- strsplit(a[1],';')[[1]]
	        m0 <- apply(mat450[normIndex,prbs],2,mean)
	        m1 <- apply(mat450[cancerIndex,prbs],2,mean)
	       DMRScores <- rbind(DMRScores,c(mean(m0),mean(m1),obj1))
	   }
    }

	DMRScores <- data.frame(DMRScores,stringsAsFactors=F)
	colnames(DMRScores) <- c('betaN','betaC','probeInf')
	kk <- t(sapply(DMRScores$probeInf,function(v){a <- strsplit(v,'#')[[1]];c(a[3],a[2],a[1])}))
	DMR <- cbind(DMRScores[,c(1,2)],kk,stringsAsFactors=F)
	colnames(DMR) <- c('betaN','betaC','GeneSymbol','BaseNumber','probeIds')
	rownames(DMR) <- 1:dim(DMR)[1]
	DMR$betaN <- as.numeric(DMR$betaN)
	DMR$betaC <- as.numeric(DMR$betaC)

	DMR$dltBeta  <- DMR$betaC-DMR$betaN
	DMR <- DMR[, c('betaN','betaC','dltBeta','GeneSymbol','BaseNumber','probeIds')]
	DMR[,1] <- as.numeric(DMR[,1])
	DMR[,2] <- as.numeric(DMR[,2])
	DMR[,3] <- as.numeric(DMR[,3])
	DMR[,5] <- as.numeric(DMR[,5])

	DMR <- DMR[order(DMR$dltBeta,decreasing=T),]
	return(DMR)
}

################################################ Test findDMR

X <- fread('F:/projects/TCGA/data/CESC_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
X <- as.data.frame(X,stringsAsFactors=F)
probInfPath='F:/projects/TCGA/data/TCGA_450K_sample.txt'
DMR  <- findDMR(X,probInfPath)

probInf              <- fread(file=probInfPath,sep="\t",header=T,check.names=F,stringsAsFactors=F)[,-2];
probInf              <- as.data.frame(probInf)
colnames(probInf)[1] <- "probeId"
probInf              <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]
print('450K probe information loaded!')

X <- X[probInf$probeId,]
R0  <- apply(X,1,function(v){sum(is.na(v))})
X <- X[R0==0,]
X <- data.frame(t(X))
probInf <- probInf[probInf$probeId %in% colnames(X),]
rownames(probInf) <- probInf$probeId
DMR  <- findDMR(X,probInfPath)

topK <- 200
