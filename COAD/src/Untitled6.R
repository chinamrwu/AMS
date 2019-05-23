rm(list=ls())
library(glmnet)
library(randomForest)
library(data.table)

setwd("G:/projects/COAD")
tmp  <- read.table('data/COAD_methy_mat_20190403.txt',sep="\t",header=T,stringsAsFactors=F,row.names=1,nrows=1)
clss <- apply(tmp,2,class)