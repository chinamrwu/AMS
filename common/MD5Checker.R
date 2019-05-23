#MD5Checker
rm(list=ls())
library(data.table)
disease <- 'all'
manifest <- 'F:/projects/TCGA/data/gdc_all450K_manifest.2019-04-03.txt' # manifest file that contains the id of files 
mf       <- read.table(manifest,sep="\t",header=T,stringsAsFactors=F)
target   <- "D:/allTCGAMethylation/" ## where to save the downloaded files

dsList <- sapply(mf$filename,function(v){a=strsplit(v,"\\.")[[1]][2];b=strsplit(a,"_")[[1]][2];b})
dmf <- mf[disease==dsList,]
if('all'==tolower(disease)){
 dmf <- mf
}

errIDs=c();
fnames <- list.files(target)
index=0;

for(id in dmf$id){
   fname=dmf$filename[dmf$id==id]
   k0=system(sprintf("fciv -md5 %s/%s/%s",target,id,fname),intern = T)
	if(length(k0) ==4){	
	  md5=strsplit(k0[4]," ")[[1]][1]
	  index <- index+1
	  
	  if(md5!= mf$md5[mf$filename==fname]){errIDs <- c(errIDs,mf$id[mf$filename==fname]);print(sprintf("ERROR:%s",id))}
	  else {
      print(sprintf('%4d %s',index,md5))
	  }
	}
	else{ errIDs <- c(errIDs,mf$id[mf$filename==fname])}
	
}

diseases <- as.character(sapply(mf$filename,function(v){a=strsplit(v,"\\.")[[1]][2];b=strsplit(a,"_")[[1]][2];b}))