dataDir <- 'F:/projects/CRCA/output/Report/'
fnames  <- list.files(dataDir)
fnames  <- fnames[grepl('DMR',fnames)]

probes <- c()
for(nme in fnames){
  tmp <- read.csv(paste0(dataDir,nme),header=T,stringsAsFactors=F)
  for(dmr in tmp$probeIds){
     probes <- unique(c(probes,strsplit(dmr,";")[[1]]))
  }
}

library(data.table)
df1 <- fread('F:/projects/allData/TCGA/CRCA_450k.txt',header=T,sep="\t")
df1 <- as.data.frame(df1)
rownames(df1) <- df1$probeId
df1 <- df1[,-1]
df1 <- df1[,grepl('-01A-|-11A-',colnames(df1))]
df1 <- data.frame(t(df1))
df1 <- df1[,probes]
df1$label <- as.character(sapply(rownames(df1),function(v){a <- strsplit(v,"-")[[1]][4];ifelse('01A'==a,'cancer','normal')}))
df1$sampleId <- rownames(df1)
df1 <- df1[,c('sampleId','label',probes)]
rownames(df1) <- 1:dim(df1)[1]
write.table(format(df1,digits=4),file="output/Report/CRCA_reportMat.txt",sep="\t",col.names = T,row.names =F,quote = F)