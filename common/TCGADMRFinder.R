# Author: wuzhichen
# Date  : 2019-05-07

# Searching for DMR in TCGA 450K dataset
# mat450: matrix or data frame,whose rows are samples and column are probes of 450k,row names are the sample names from TCGA
# probInf: data.frame for probe information, each level3 450K sample files from TCGA is a probe information,each row is a 
# probe and the row names of probInf are the probeId
# breakCutOff: each DMR contains multiple probes of 450K array, when too many probes are closed, this cutoff is used to break long
#  DMR to ones with suitable length
# output: a data frame£¬which contains DMRs 
#

library(data.table)
TCGADMRFinder <- function(mat450,probInf,breakCutOff=0.2) {
    breakCutOff <- 0.2
    if(is.null(mat450) | is.null(probInf)){
         print("Missing one of input datasets")
			return(NULL)
	 }
    rn <- rownames(probInf)[1]
	 index <- 1
	 if(!startsWith(rn,'cg')){
	    v0 <- as.character(probInf[1,])
		 while(!startsWith(v0[index],'cg')){index <- index+1}
		 rownames(probInf) <- probInf[,index]
	 }
	 rm(rn)
    probInf <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]

		 

    mat450  <- mat450[,rownames(probInf)]
    mat450  <- mat450[grepl('-01A-|-11A-',rownames(mat450)),]
	 mat450  <- mat450[,rownames(probInf)]

	 label <- as.character(sapply(rownames(mat450),function(v){a <- strsplit(v,"-")[[1]][4];ifelse(a=='01A','cancer','normal')}))
	 normalIndex <- which(label == 'normal')
	 cancerIndex <- which(label == 'cancer')
    R0 <- apply(mat450,2,function(v){sum(is.na(v))})
	 mat450[,R0!=0] <- -breakCutOff
	
	 mn <- apply(mat450[normalIndex,],2,mean)
	 mc <- apply(mat450[cancerIndex,],2,mean)

	 ############################################
	 geneSymbolByIsland <- function(v){
			 v0  <- probInf$Gene_Symbol[probInf$CGI_Coordinate==v]
			 v1 <- unique(strsplit(paste0(v0,collapse=";"),";")[[1]])
			 paste0(v1,collapse=";")
    }
	
	###########################################
   isLands <- unique(probInf$CGI_Coordinate)

	dmr <- sapply(isLands,function(obj) {
	#for( obj in isLands){
		 geneSyb <-geneSymbolByIsland(obj);
		 chr <- probInf$probInf$Chromosome[probInf$CGI_Coordinate==obj][1]

       probes <- probInf[probInf$CGI_Coordinate==obj,c('Start','probeId')];
		 probes <- probes[order(probes$Start),]
		 probes <- probes$probeId
		 score  <- apply(rbind(mn[probes],mc[probes]),2,function(v){
            dlt <- v[2]-v[1];
				v[2]*dlt/(v[1]+0.0001)
		 })

		tf <- score >= breakCutOff
		i <- 1
		L <- length(tf)
		DMRs <- list()
		while(i < L){
			rg <- c()
			while(tf[i]==T & i <=L){rg <- c(rg,i); i <- i+1 }
			if(length(rg)>=2) { DMRs[[length(DMRs)+1]] <- rg}
			i <- i+1
		}
	   
	  s1 <- as.character(sapply(DMRs,function(v){
           iStart <- probInf[probes[v[1]],'Start']
			  iEnd   <- probInf[probes[v[length(v)]],'End']
			  iWidth <- abs(iEnd-iStart)+1
			  prbs   <- paste0(probes[v],collapse=";")
			  t0 <- paste0(c(geneSyb,chr,iStart,iEnd,iWidth,prbs),collapse="-")
			  t0
		}))
		s1
	})
	Len <- sapply(dmr,length)
	index <- which(Len>0)
	dmr[index]

}

############################################## Test it 
COAD <- fread('F:/projects/TCGA/data/COAD_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
READ <- fread('F:/projects/TCGA/data/READ_450k.txt',header=T,sep="\t",stringsAsFactors=F,check.names=F)
mat450 <- cbind(as.data.frame(COAD),as.data.frame(READ[,-1]));

rm(COAD)
rm(READ)
rownames(mat450) <- mat450$probeId
mat450 <- mat450[,-1]
mat450 <- data.frame(t(mat450),stringsAsFactors=F,check.names=F)

probInf <- read.table('F:/projects/TCGA/data/TCGA_450K_sample.txt',sep="\t",header=T,stringsAsFactors=F)[,-2]
colnames(probInf)[1] <- "probeId"
probInf <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]

system.time( result <- TCGADMRFinder(mat450,probInf))



