# This script generate a methylation matrix from samples file

rm(list=ls())
setwd("G:/projects/COAD")
srcDir <- "data/tcga/" # where the tcga methylation files exist
saveTo <- "data/COAD450K_filtered_20190410.txt";

fnames <- list.files(srcDir)

probInfo <- read.table(paste0(srcDir,fnames[1]),sep="\t",header=T,stringsAsFactors=F)[,-2]
colnames(probInfo)[1] <- c("probeId")
write.table(probInfo,file="probInfo.txt",sep="\t",



X      <- c()
index  <- 0
for(fn in fnames){
     X     <- cbind(X,read.table(paste0(srcDir,fn),sep="\t",header=T)[,2])
     index <- index+1
     print(sprintf("End processing %d/%d file",index,length(fnames)))
}
rownames(X) <- probInfo[,1]
colnames(X) <- as.character(sapply(fnames,function(v){strsplit(v,"\\.")[[1]][6]}))

probInfo <- probInfo[!probInfo$Chromosome %in% c('*','chrX','chrY') & probInfo$Gene_Symbol !='.' & probInfo$Feature_Type=='Island',]
rnames <- probInfo[,1]
X <- X[rnames,]

missingR <- as.integer(apply(X,1,function(v){sum(is.na(v))}))
X <- X[missingR==0,]

X <- data.frame(X,stringsAsFactors=F)

if(F){ 
	print("Now write data frame to local file......")
	write.table(X,file=saveTo,sep="\t",col.names=T,row.names=T,quote=F)
}


