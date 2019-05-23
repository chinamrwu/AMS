rm(list=ls())
library(glmnet)
library(randomForest)
library(pROC)
library(data.table)

setwd('F:/projects/TCGA')
disease <- 'THCA'

mat <- fread(sprintf('data/%s_450k.txt',disease),sep="\t",header=T,stringsAsFactors=F,check.names=F)
mat <- data.frame(mat,stringsAsFactors=F,check.names=F)
rownames(mat) <- mat$probeId
mat <- mat[,-1]
R0  <- apply(mat,1,function(v){sum(is.na(v))})
mat <- mat[R0==0,grepl('-01A-|-11A-',colnames(mat))]


probInf <- fread('data/probeInf.txt',sep="\t",header=T,stringsAsFactors=F)
probInf <- data.frame(probInf,stringsAsFactors=F,check.names=F)
probInf <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]

sharedProbs <- intersect(probInf$probeId,rownames(mat))
mat <- data.frame(t(mat[sharedProbs,]))
mat$label <- as.character(sapply(rownames(mat),function(v){a=strsplit(v,"-")[[1]][4];ifelse('01A'==a,'cancer','normal')}))
mat <- mat[,c('label',sharedProbs)]
probInf <- probInf[probInf$probeId %in% sharedProbs,]


index1 <- which(mat$label=='normal')
index2 <- which(mat$label=='cancer')

DMPrb <- apply(mat[,-1],2,function(v){
  m1 <- mean(v[index1])
  m2 <- mean(v[index2])
  c(m1,m2,m2-m1)
})

DMPrb <- data.frame(t(DMPrb))
colnames(DMPrb) <- c('betaN','betaC','dltB')

dmr <- DMPrb[DMPrb$dltB > 0.35,]
dmr <- dmr[order(dmr$dltB,decreasing=T),]

genes <- sapply(rownames(dmr),function(v){
   v0 <- probInf$Gene_Symbol[probInf$probeId==v]
	v0 <- strsplit(v0,";")[[1]];
	paste0(unique(v0),collapse=";")
})

#cvfit  <- cv.glmnet(as.matrix(mat[,-1]),as.factor(mat$label),family='binomial',alpha=0.5,type.measure='class');
#cf <- coef(cvfit,s='lambda.min')
#cf <- rownames(cf)[cf[,1]!=0][-1]

#########################  volcano
if(F){
drawVolcano <- function(df1,outFile=NULL){ 
   label <- df1$label
   lbls = unique(df1$label)
   table(df1$label) # A143,C75
   fc <- apply(2^df1[,colnames(df1)!='label'],2,function(x){
       log2( mean(na.omit(x[label==lbls[1]])) /mean(na.omit(x[label==lbls[2]])))
    })
   
   df1[is.na(df1)] <- 0
   pValue <- apply(df1[,colnames(df1)!='label'], 2, function(v) {
       p1 <- t.test(v[label == lbls[1]], v[label == lbls[2]], paired = F, var.equal = F)$p.value
       p1 <-  p.adjust(p1,method="BH")
       p1
    }) 

  pdf(pdfPath)
     plot(fc, -log10(pValue), col = '#00000033', pch = 19,xlab = 'log2(FC)', ylab = '-log10(p-value)', main = strTitle)
     abline(h = 1.3, v = c(-log2(1.5),log2(1.5)), lty = 2, lwd = 1)
  
     up  <- fc >= log2(1.5) & pValue <= 0.05
     points(fc[up], -log10(pValue[up]), col = 1,bg = brewer.pal(9,"YlOrRd")[6], pch = 21, cex = 2)
 
     down <- fc <= -log2(1.5) & pValue <= 0.05
     points(fc[down], -log10(pValue[down]), col = 1,bg = brewer.pal(11,"RdBu")[9], pch = 21, cex = 2)
  dev.off()
  
  name = list(up = data.frame(prot = colnames(df1[,colnames(df1)!='label'])[up],fc = fc[up], p_value = pValue[up],
                            type = rep("Upregulation",sum(up)),stringsAsFactors=F),
            down = data.frame(prot = colnames(df1[,colnames(df1)!='label'])[down],fc = fc[down], p_value = pValue[down],
                              type = rep("Downregulation",sum(down))),stringsAsFactors=F)
   name1 = rbind(name[[1]],name[[2]])
   rownames(name1) <- 1:dim(name1)[1]
   
   if( !is.null(outFile) ){
     write.table(name1,file=outFile,sep="\t",col.names=T,row.names=F,quote=F)
   }
   #write.xlsx(name1, "Table_1_TPD_volcano_prot_AC10_fc1.5_190131.xlsx",showNA = T,row.names = F)
   sum(up)
   sum(down)
   return(name1)
}
}