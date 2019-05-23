#Extract 450K files for a specified disease and contsruct matrix 
library(data.table)
library(readr)
extractMatrix <- function(disease){
		manifest <- "F:/projects/gdc_all450K_manifest.2019-04-03.txt" # manifest file that contains the id of files 
		mf       <- read.table(manifest,sep="\t",header=T,stringsAsFactors=F)
		target   <- "D:/allTCGAMethylation/" ## where to save the downloaded files

		dsList <- sapply(mf$filename,function(v){a=strsplit(v,"\\.")[[1]][2];b=strsplit(a,"_")[[1]][2];b})
		dmf <- mf[disease==dsList,]
		if('all'==tolower(disease)){dmf <- mf}
      
		fname <- dmf$filename[1]
      id <- dmf$id[dmf$filename==fname];
		tmp <- read.table(sprintf('%s%s/%s',target,id,fname),sep="\t",header=T,stringsAsFactors=F);
		clsses <- as.character(sapply(tmp,class))
		k0 <- tmp[,2]
		probeId <- tmp[,1]
      probeNumber=dim(tmp)[1]
		L=dim(dmf)[1]
		df0 <- sapply(dmf$filename[-1],function(fname){
		    id <- dmf$id[dmf$filename==fname];
		    as.data.frame(fread(sprintf('%s%s/%s',target,id,fname),sep="\t",header=T,stringsAsFactors=F,colClasses=clsses,nrow=probeNumber))[,2]
		})
		df0 <- cbind(k0,df0)
		colnames(df0)[1] <- fname;
		colnames(df0) <- as.character(sapply(colnames(df0),function(v){strsplit(v,"\\.")[[1]][6]}))
		clnames <- colnames(df0)
		df0 <- data.frame(df0,check.names=F)
		df0$probeId <- probeId
		df0 <- df0[,c('probeId',clnames)]
		return(df0)
}


####################################################
if(F){
		manifest <- "F:/projects/gdc_all450K_manifest.2019-04-03.txt" # manifest file that contains the id of files 
		mf       <- read.table(manifest,sep="\t",header=T,stringsAsFactors=F)
		target   <- "D:/allTCGAMethylation/" ## where to save the downloaded files

		ESCA <- extractMatrix("ESCA") # 202 食管癌
		write_tsv(ESCA,path='G:/ESCA_450k.txt');

		STAD <- extractMatrix("STAD") # 397 胃癌
		write_tsv(STAD,path='G:/STAD_450k.txt');

		PRAD <- extractMatrix("PRAD") # 553 前列腺
		write_tsv(PRAD,path='G:/PRAD_450k.txt');

		THCA <- extractMatrix("THCA") #
		write_tsv(THCA,path='G:/THCA_450k.txt');

		THYM <- extractMatrix("THYM")  # 胸腺
		write_tsv(THYM,path='G:/THYM_450k.txt');


		UCEC <- extractMatrix("UCEC")  # 485 samples 子宫内膜癌
		write_tsv(UCEC,path='G:/UCEC_450k.txt');

		HNSC <- extractMatrix("HNSC")
		write_tsv(HNSC,path='G:/HNSC_450k.txt');

		BLCA <- extractMatrix("BLCA")
		write_tsv(BLCA,path='G:/BLCA_450k.txt');
      
		LUAD <- extractMatrix("LUAD")
		write_tsv(BLCA,path='F:/projects/TCGA/data/LUAD_450k.txt');
		LUSC <- extractMatrix("LUSC")
		write_tsv(LUSC,path='F:/projects/TCGA/data/LUSC_450k.txt');
		LAML <- extractMatrix("LAML")
		write_tsv(LAML,path='F:/projects/TCGA/data/LAML_450k.txt');

		KICH <- extractMatrix("KICH")
		write_tsv(KICH,path='F:/projects/TCGA/data/KICH_450k.txt');
		KIRC <- extractMatrix("KIRC")
		write_tsv(KIRC,path='F:/projects/TCGA/data/KIRC_450k.txt');
		KIRP <- extractMatrix("KIRP")
		write_tsv(KIRP,path='F:/projects/TCGA/data/KIRP_450k.txt');

		##################################
		fnames <- list.files('F:/projects/TCGA/data/')
		fnames <- fnames[grepl('450k.txt',fnames)]
		Dones <- as.character(sapply(fnames,function(v){strsplit(v,"_")[[1]][1]})) 

		diseases <- unique(dsList)
		toDo <- setdiff(diseases,Dones)
      
		for(job in toDo){
        obj <- extractMatrix(job);
		  write_tsv(obj,path=sprintf('F:/projects/TCGA/data/%s_450k.txt',job));
		}



}




