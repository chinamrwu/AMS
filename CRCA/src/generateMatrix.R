rm(list=ls())
setwd("G:/projects/CRC")

fnames <- list.files("data/STAD")
X      <- c()
index  <- 0
for(fn in fnames){
     X     <- cbind(X,read.table(paste0("data/STAD/",fn),sep="\t",header=T,stringsAsFactors=F)[,2])
     index <- index+1
     print(sprintf("End processing %d/%d file",index,length(fnames)))
     #print(dim(X))
}
X <- data.frame(t(X))
probs <- read.table(paste0("data/STAD/",fnames[1]),sep="\t",header=T,stringsAsFactors=F)[,1]
colnames(X) <- probs
rownames(X) <- as.character(sapply(fnames,function(v){strsplit(v,"\\.")[[1]][6]}))
if(F){ 
	print("Now write data frame to local file......")
	write.table(X,file="data/STAD_cancerMatrix_0402.txt",sep="\t",col.names=T,row.names=T,quote=F)
}


