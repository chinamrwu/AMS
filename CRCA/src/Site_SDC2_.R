
setwd('F:/projects/CRC')
fnames <- list.files("output")
fnames <- fnames[grepl('CRC_Prediction_SDC2',fnames)]

pred <- read.table(sprintf('output/%s',fnames[grepl('_DMR1_',fnames)]),sep="\t",header=T,stringsAsFactors=F)
v0 <- table(pred[pred$observed=='Cancer','site'])
v1 <- rep(0,length(v0))
names(v1) <- names(v0)
tmp <- table(pred[(pred$result==FALSE & pred$observed=='Cancer'),'site'])
v1[names(tmp)] <- tmp

total <- cbind(v0,v1)


pred <- read.table(sprintf('output/%s',fnames[grepl('_DMR2_',fnames)]),sep="\t",header=T,stringsAsFactors=F)
v1 <- rep(0,length(v0))
names(v1) <- names(v0)
tmp <- table(predDMR2[(pred$result==FALSE & pred$observed=='Cancer'),'site'])
v1[names(tmp)] <- tmp
total <- cbind(total,v1)

pred <- read.table(sprintf('output/%s',fnames[grepl('_DMR3_',fnames)]),sep="\t",header=T,stringsAsFactors=F)
v1 <- rep(0,length(v0))
names(v1) <- names(v0)
tmp <- table(predDMR2[(pred$result==FALSE & pred$observed=='Cancer'),'site'])
v1[names(tmp)] <- tmp
total <- cbind(total,v1)


pred <- read.table(sprintf('output/%s',fnames[grepl('DMR1\\+DMR2',fnames)]),sep="\t",header=T,stringsAsFactors=F)
v1 <- rep(0,length(v0))
names(v1) <- names(v0)
tmp <- table(predDMR2[(pred$result==FALSE & pred$observed=='Cancer'),'site'])
v1[names(tmp)] <- tmp
total <- cbind(total,v1)

pred <- read.table(sprintf('output/%s',fnames[grepl('DMR1\\+DMR3',fnames)]),sep="\t",header=T,stringsAsFactors=F)
v1 <- rep(0,length(v0))
names(v1) <- names(v0)
tmp <- table(predDMR2[(pred$result==FALSE & pred$observed=='Cancer'),'site'])
v1[names(tmp)] <- tmp
total <- cbind(total,v1)

pred <- read.table(sprintf('output/%s',fnames[grepl('DMR2\\+DMR3',fnames)]),sep="\t",header=T,stringsAsFactors=F)
v1 <- rep(0,length(v0))
names(v1) <- names(v0)
tmp <- table(predDMR2[(pred$result==FALSE & pred$observed=='Cancer'),'site'])
v1[names(tmp)] <- tmp
total <- cbind(total,v1)

total <- data.frame(total)
colnames(total) <- c("cancerNumber","DMR1","DMR2","DMR3","DMR1+DMR2","DMR1+DMR3","DMR2+DMR3")
total$site <- rownames(total)
total <- total[,c("site","cancerNumber","DMR1","DMR2","DMR3","DMR1+DMR2","DMR1+DMR3","DMR2+DMR3")]
total <- total[order(total$cancerNumber,decreasing=T),]
rownames(total) <- 1:dim(total)[1]
