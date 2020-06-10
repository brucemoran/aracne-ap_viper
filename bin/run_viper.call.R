#! /usr/bin/env R

argsIn <- commandArgs(trailingOnly=TRUE)
source(argsIn[1])
NETWORK <- argsIn[2]
EXPRMAT <- argsIn[3]
METADATA <- argsIn[4]
TAG <- argsIn[5]

runMsViper(NETWORK, EXPRMAT, METADATA, TAG)
RDATAS <- grep("ens2ext",dir(pattern="RData"),invert=TRUE,value=T)
lapply(RDATAS, filterSigRes)
