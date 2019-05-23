
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

