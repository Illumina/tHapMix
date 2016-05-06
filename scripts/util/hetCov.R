#setwd("//ukch-prd-isi01/scratch_tmp/users/ccolombo/evaluation/EvaluateCov")
setwd("/illumina/scratch/tmp/users/ccolombo/evaluation/EvaluateCov")
rm(list=ls(all=TRUE))
library(plyr)
library(ggplot2)
library(grid)

# id <- "cov_simGeL004_n2_m80_v-1"
# id <- "cov_simGeL004_n2_m80_v-1_truth"
# id <- "cov_simGeL004_n2_m60_v-1_truth"
# id <- "cov_simGeL004_n1_m-1_v-1"
# id <- "cov_simGeL004_n1_m-1_v-1_truth"
# id <- "cov_simHCC2218_n2_m70_v-1"
#id <- "cov_simHCC2218_n2_m70_v-1_truth"
#id = "cov_simHCC1187_n4_m50_v60_sl_truth"
# id = "cov_simNorm_truth"
# id = "cov_simNorm"
# id = "cov_simGeL004_n1_purity100_corr_truth"
# id = "cov_simHCC2218_n1_purity100_corr_truth"

#id <- "cov_simHCC2218_n1_m-1_v-1"
# id <- "cov_simGeL004_n1_purity100_truth"

id <- "cov_simGeL004_n3_m50_v-1"

args <- commandArgs(trailingOnly = T)
if (length(args) > 0)
  id <- args[1]

fileBaseName <- paste0(ifelse(substr(id, nchar(id)-5, nchar(id))=="_truth", substr(id, 1, nchar(id)-6), id), "/", id)

plotDir <- ifelse(substr(id, nchar(id)-5, nchar(id))=="_truth", paste0("./", substr(id, 1, nchar(id)-6), "/plots/truth"), paste0("./", id, "/plots/canvas"))

if (!(plotDir %in% list.dirs(paste0("./", ifelse(substr(id, nchar(id)-5, nchar(id))=="_truth", substr(id, 1, nchar(id)-6), id)), recursive=T)))
  dir.create(plotDir, recursive=T)
plotDir <- "/illumina/scratch/tmp/users/ccolombo/evaluation/EvaluateCov/new/"
#plotDir = "//ukch-prd-isi01/scratch_tmp/users/ccolombo/evaluation/EvaluateCov/new/"
tit <- substr(id, 5, nchar(id))
print(plotDir)

dataPerc <- read.table(paste0(fileBaseName, "_bin_perc_all.bed"), sep="\t", col.names=c("chrV", "startV", "endV", "cnA", "cnB", "perc", "chrB", "startB", "endB", "cnCall", "obsCov", "expCov"))
dataPerc = subset(dataPerc, !(cnA==1 & cnB==1) & chrV != "chrX" & chrV != "chrY")
dataPerc$dist <- abs(dataPerc$obsCov - dataPerc$expCov)
dataPerc$distRel <- dataPerc$dist/dataPerc$obsCov


# if (sum(dataPerc$perc <= 0.25 | dataPerc$perc >= 0.75)>0 & sum(dataPerc$perc < 0.75 & dataPerc$perc > 0.25)>0) {
#   dataPerc$bin <- as.factor(ifelse(dataPerc$perc <= 0.25 | dataPerc$perc >= 0.75, "non-heterogeneous", "heterogeneous"))
#   d_th <- theme(plot.title=element_text(size=12, vjust=1))
#   p <- ggplot(dataPerc, aes(x=bin, y=dist)) + geom_boxplot() +
#     ylim(0, 100) + ylab("|obsCov - expCov|") + d_th +
#     ggtitle(tit)
#   ggsave(paste0(plotDir, id, "_het_perc_dist_bp.png"), p, width=7, height=4, dpi=200)
#   p <- ggplot(dataPerc, aes(x=bin, y=distRel)) + geom_boxplot() +
#     ylim(0, 1) + ylab("|obsCov - expCov|/obsCov") + d_th +
#     ggtitle(tit)
#   ggsave(paste0(plotDir, id, "_het_perc_distRel_bp.png"), p, width=7, height=4, dpi=200)
#   
#   p <- ggplot(dataPerc, aes(color=bin, x=dist)) + geom_density() +
#     xlim(0, 100) + xlab("|obsCov - expCov|") + scale_color_discrete(name="variant %") + d_th +
#     ggtitle(tit)
#   ggsave(paste0(plotDir, id, "_het_perc_dist_dens.png"), p, width=7, height=5, dpi=200)
#   p <- ggplot(dataPerc, aes(color=bin, x=distRel)) + geom_density() +
#     xlim(0, 1) + xlab("|obsCov - expCov|/obsCov") + scale_color_discrete(name="variant %") + d_th +
#     ggtitle(tit)
#   ggsave(paste0(plotDir, id, "_het_perc_distRel_dens.png"), p, width=7, height=5, dpi=200)
# }
# 
# 
# plotDir <- "/illumina/scratch/tmp/users/ccolombo/evaluation/EvaluateCov/new1/"
# dataPerc$bin <- as.factor(ifelse(dataPerc$perc >=0.98, "non-heterogeneous", "heterogeneous"))
plotDir <- "/illumina/scratch/tmp/users/ccolombo/evaluation/EvaluateCov/new2/"
dataPerc$bin <- as.factor(ifelse(dataPerc$perc <= 0.1 | dataPerc$perc >= 0.9, "non-heterogeneous", "heterogeneous"))
  d_th <- theme(plot.title=element_text(size=12, vjust=1))
  dataPerc$dist <- abs(dataPerc$obsCov - dataPerc$expCov)
  dataPerc$distRel <- dataPerc$dist/dataPerc$obsCov
  p <- ggplot(dataPerc, aes(x=bin, y=dist)) + geom_boxplot() +
    ylim(0, 100) + ylab("|obsCov - expCov|") + d_th +
    ggtitle(tit)
  ggsave(paste0(plotDir, id, "_het_perc_dist_bp.png"), p, width=7, height=4, dpi=200)
  p <- ggplot(dataPerc, aes(x=bin, y=distRel)) + geom_boxplot() +
    ylim(0, 1) + ylab("|obsCov - expCov|/obsCov") + d_th +
    ggtitle(tit)
  ggsave(paste0(plotDir, id, "_het_perc_distRel_bp.png"), p, width=7, height=4, dpi=200)
  
  p <- ggplot(dataPerc, aes(color=bin, x=dist)) + geom_density() +
    xlim(0, 100) + xlab("|obsCov - expCov|") + scale_color_discrete(name="variant %") + d_th +
    ggtitle(tit)
  ggsave(paste0(plotDir, id, "_het_perc_dist_dens.png"), p, width=7, height=5, dpi=200)
  p <- ggplot(dataPerc, aes(color=bin, x=distRel)) + geom_density() +
    xlim(0, 1) + xlab("|obsCov - expCov|/obsCov") + scale_color_discrete(name="variant %") + d_th +
    ggtitle(tit)
  ggsave(paste0(plotDir, id, "_het_perc_distRel_dens.png"), p, width=7, height=5, dpi=200)

