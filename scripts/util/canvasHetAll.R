setwd("//ukch-prd-isi01/scratch_tmp/users/ccolombo/evaluation/EvaluateCoV")
rm(list=ls(all=TRUE))
library(plyr)

ttestDistPerc <- function(thrP, simNames) {
  print(thrP)
  res <- data.frame(thrP=thrP, id=simNames)
  res <- ddply(res, .(id), calcPvalue)
  return(res)
}

calcPvalue <- function(res) {
#   if (!(file.exists(paste0(res$id, "/", res$id, "_perc.bed"))))
#     print(paste0(res$id, "/", res$id, "_perc.bed"))
#   if(!(file.exists(paste0(res$id, "/", res$id, "_truth_perc.bed"))))
#     print(paste0(res$id, "/", res$id, "_truth_perc.bed"))
  dataPerc <- read.table(paste0(res$id, "/", res$id, "_perc.bed"), sep="\t",
                         col.names=c("chrV", "startV", "endV", "cnA", "cnB", "perc", "chrS", "startS", "endS", "cnCall", "obsCov", "expCov"))
  dataPerc$dist <- abs(dataPerc$obsCov - dataPerc$expCov)
  dataPerc$distRel <- dataPerc$dist/dataPerc$obsCov
  dataPerc$het <- dataPerc$perc <= res$thrP
  lev <- length(unique(dataPerc$het)) > 1
  res$pDist <-  ifelse(lev, t.test(dist ~ het, data=dataPerc, na.rm=T)$p.value, NA)
  res$pDistRel <- ifelse(lev, t.test(distRel ~ het, data=subset(dataPerc, distRel!=Inf), na.rm=T)$p.value, NA)
  
  dataPercTruth <- read.table(paste0(res$id, "/", res$id, "_truth_perc.bed"), sep="\t",
                         col.names=c("chrV", "startV", "endV", "cnA", "cnB", "perc", "chrS", "startS", "endS", "cnCall", "obsCov", "expCov"))
  dataPercTruth$dist <- abs(dataPercTruth$obsCov - dataPercTruth$expCov)
  dataPercTruth$distRel <- dataPercTruth$dist/dataPercTruth$obsCov
  dataPercTruth$het <- dataPercTruth$perc <= res$thrP
  lev <- length(unique(dataPercTruth$het)) > 1
  res$pTruthDist <-  ifelse(lev, t.test(dist ~ het, data=dataPercTruth, na.rm=T)$p.value, NA)
  res$pTruthDistRel <- ifelse(lev, t.test(distRel ~ het, data=subset(dataPercTruth, distRel!=Inf), na.rm=T)$p.value, NA)

  return(res)
}

covSimNames <- list.dirs(".", recursive=F, full.names=F)
covSimNames <- covSimNames[sapply(covSimNames, function(x) substr(x, 1, 7)=="cov_sim")]
covSimNames <- covSimNames[sapply(covSimNames, function(x) substr(x, 1, 11)!="cov_simNorm")]
resTot <- adply(c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8), 1, ttestDistPerc, covSimNames)

summary(subset(resTot, thrP==0.4))
by(resTot$pTruthDistRel, resTot$thrP, median, na.rm=T)
by(resTot$pTruthDist, resTot$thrP, median, na.rm=T)
subset(resTot, thrP==0.5 & pTruthDist > 0.5)
by(resTot$pTruthDist, resTot$thrP, function(x) sum(x<0.05, na.rm=T)/length(x))
by(resTot$pTruthDistRel, resTot$thrP, function(x) sum(x<0.05, na.rm=T)/length(x))
by(resTot$pDist, resTot$thrP, function(x) sum(x<0.05, na.rm=T)/length(x))
by(resTot$pDistRel, resTot$thrP, function(x) sum(x<0.05, na.rm=T)/length(x))

library(ggplot2)

