rm(list=ls())
library(data.table)
library(sqldf)
setwd('F:/projects/CESC')

clinic <- read.table('F:/projects/TCGA/data/clinical.tsv',sep="\t",header=T,stringsAsFactors=F)
clinic <- clinic[clinic$project_id=='TCGA-CESC',]
clinic <- clinic[,c(2,4,5,13,14,21,25)]
clinic$age_at_diagnosis <- as.numeric(clinic$age_at_diagnosis)
clinic <- clinic[clinic$submitter_id %in% patients,]
clinic$days_to_last_follow_up <- as.numeric(clinic$days_to_last_follow_up)
index <- which(!is.na(clinic$days_to_last_follow_up))
clinic <- clinic[index,]
clinic$vital_status <- sapply(clinic$vital_status,function(v){ifelse('alive'==v,0,1)})

probInf              <- fread(file='F:/projects/TCGA/data/TCGA_450k_sample.txt',sep="\t",header=T,check.names=F,stringsAsFactors=F)[,-2];
probInf              <- as.data.frame(probInf)
colnames(probInf)[1] <- "probeId"
rownames(probInf)    <- probInf$probeId
probInf              <- probInf[!probInf$Chromosome %in% c('*','chrX','chrY') & probInf$Gene_Symbol !='.' & probInf$Feature_Type=='Island',]

probePAX1 <- probInf$probeId[grepl(';PAX1;',probInf$Gene_Symbol)]
probeSOX1 <- sqldf("SELECT * FROM probInf where Gene_Symbol ='SOX1'")[,1]
probes <- unique(c(probePAX1,probeSOX1))


mat <- fread('F:/projects/TCGA/data/CESC_450k.txt',sep='\t',header=T,stringsAsFactors=F)
mat <- as.data.frame(mat)
rownames(mat) <- mat$probeId
mat <- mat[,-1]
mat <- data.frame(t(mat[probes,]))

submitters <- as.character(sapply(rownames(mat),function(v){a <- strsplit(v,"-")[[1]][1:3];paste0(a,collapse="-")}));
tmp <- submitters[which(!submitters %in% submitters[duplicated(submitters)])]
mat <- mat[match(tmp,submitters),]
rownames(mat) <- tmp

patients <- intersect(tmp,clinic$submitter_id)
mat <- mat[patients,]
mat$submitter <- rownames(mat)
R0 <- apply(mat,2,function(v){sum(is.na(v))})
mat <- mat[,R0==0]

probePAX1 <- intersect(probePAX1,colnames(mat))
probeSOX1 <- intersect(probeSOX1,colnames(mat))
######################################
m01    <-   mean(apply(mat[,probePAX1],2,mean),na.rm=T)
mv1    <-   apply(mat[,probePAX1],1,mean,na.rm=T)

lowIndex <- which(mv1 < 0.3) ##  16 samples
highIndex <- which(mv1 > 0.6) ## 44 samples

tmp <- data.frame('submitter'=names(mv1)[c(lowIndex,highIndex)],'methylation'=mv1[c(lowIndex,highIndex)],stringsAsFactors=F)
tmp$methylationLevel <- c(rep('hypomethylation',length(lowIndex)),rep('hypermethylation',length(highIndex)))
matPAX1 <- merge(x=clinic,y=tmp,by.x='submitter_id',by.y='submitter')


m02    <-  mean(apply(mat[,probeSOX1],2,mean,na.rm=T),na.rm=T)
mv2    <-  apply(mat[,probeSOX1],1,mean,na.rm=T)

lowIndex <-  which(mv2 < 0.4) ##  22 samples
highIndex <- which(mv2 > 0.7) ##  21 samples

tmp <- data.frame('submitter'=names(mv2)[c(lowIndex,highIndex)],'methylation'=mv2[c(lowIndex,highIndex)],stringsAsFactors=F)
tmp$methylationLevel <- c(rep('hypomethylation',length(lowIndex)),rep('hypermethylation',length(highIndex)))
matSOX1 <- merge(x=clinic,y=tmp,by.x='submitter_id',by.y='submitter')

#########################################################
library(survminer)
library(survival)

fit1 <- survfit(Surv(days_to_last_follow_up, vital_status) ~ methylationLevel,data=matPAX1)
p1 <- ggsurvplot(fit1, data = matPAX1, pval = T, risk.table = TRUE,title="CESC:PAX1",legend=c(0.80,0.2),
 legend.title = "", legend.labs = c("Hypomethlation", "Hypermethylation"))
p1$plot <- p1$plot + 
theme(legend.text = element_text(size = 16, color = "black"))
p1

fit2 <- survfit(Surv(days_to_last_follow_up, vital_status) ~ methylationLevel,data=matSOX1)
p2 <- ggsurvplot(fit2, data = matSOX1, pval = T, risk.table = TRUE,title="CESC:SOX1",legend=c(0.80,0.2),
 legend.title = "", legend.labs = c("Hypomethlation", "Hypermethylation"))
p2$plot <- p2$plot + 
theme(legend.text = element_text(size = 16, color = "black"))
p2


##########################
tmp <- data.frame('avgBeta'=c(mv1,mv2),'geneSymbol'=c(rep('PAX1',length(mv1)),rep('SOX1',length(mv2))))

p1= ggplot(tmp, aes(geneSymbol, avgBeta, fill=geneSymbol))+geom_boxplot()+
stat_boxplot(geom ='errorbar')+
theme(
	    axis.text.x = element_text(size=18,color="black"),
	    axis.text.y = element_text(size=18,color="black"),
		 axis.line.x = element_line(color="black", size = 0.25),
	    axis.line.y = element_line(color="black", size = 0.25),
	    panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
		 panel.background = element_blank()
	)