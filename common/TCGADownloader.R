# This script download files from TCGA using gdc-client 
# input: A manifest file 
# output: downloaded files from TCGA
options("width"=400)
manifest <- "F:/projects/gdc_all450K_manifest.2019-04-03.txt" # manifest file that contains the id of files 
mf       <- read.table(manifest,sep="\t",header=T,stringsAsFactors=F)
target   <- "D:/allTCGAMethylation/" ## where to save the downloaded files
UUIDs    <- setdiff(mf$id,list.files(target))

L        <- length(UUIDs)
errIds <- c()
print(sprintf("%d files to download",L))
for(i in 1:L){
  uuid  <- UUIDs[i]
  fname <- mf$filename[mf$id==uuid]

  tryCatch({
            system(command=sprintf("gdc download -d %s %s",target,uuid),intern = T) 
            print(sprintf("%4d--%s  %s ",i,uuid,fname))
				#print(k0)
			  }, warning = function(w) {
              print(sprintf("A warning is issued when downloading %s",uuid))
				  print(w)
           }, error = function(e) {
             print(e)
				 errIds <- c(errIds,uuid)
           }, finally = {
             
           }
  )
}

if(length(errIds) >0){ ###
 for(uuid in errIds){
  tryCatch({
            k0 <- system(command=sprintf("gdc download -d %s %s",target,uuid),intern = T) 
            print(sprintf("remedy download --%s  ",uuid))
				print(k0)
			  }, warning = function(w) {
              print(sprintf("A warning is issued when downloading %s",uuid))
           }, error = function(e) {
             print(e)
				 errIds <- c(errIds,uuid)
           }, finally = {
             
           }
  )
 }
}
########################## check MD5 for data integrity
errIDs=c();
fnames <- list.files(target)

for(fn in fnames){
   k0=system(sprintf("fciv -md5 data/tcga/%s",fn),intern = T)
	if(length(k0) ==4){	
	  md5=strsplit(k0[4]," ")[[1]][1]
	  if(md5!= mf$md5[mf$filename==fn]){errIDs <- c(errIDs,mf$id[mf$filename==fn])}
	}
	else{ errIDs <- c(errIDs,mf$id[mf$filename==fn])}
}
################################
if(length(errIds) >0){ ###
 for(uuid in errIds){
  tryCatch({
            k0 <- system(command=sprintf("gdc download -d %s %s",target,uuid),intern = T) 
            print(sprintf("remedy download --%s  ",uuid))
				print(k0)
			  }, warning = function(w) {
              print(sprintf("A warning is issued when downloading %s",uuid))
           }, error = function(e) {
             print(e)
				 errIds <- c(errIds,uuid)
           }, finally = {
             
           }
  )
 }
}