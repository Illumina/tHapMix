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

args <- commandArgs(trailingOnly = T)
if (length(args) > 0)
  id <- args[1]

fileBaseName <- paste0(ifelse(substr(id, nchar(id)-5, nchar(id))=="_truth", substr(id, 1, nchar(id)-6), id), "/", id)

plotDir <- ifelse(substr(id, nchar(id)-5, nchar(id))=="_truth", paste0("./", substr(id, 1, nchar(id)-6), "/plots/truth"), paste0("./", id, "/plots/canvas"))

if (!(plotDir %in% list.dirs(paste0("./", ifelse(substr(id, nchar(id)-5, nchar(id))=="_truth", substr(id, 1, nchar(id)-6), id)), recursive=T)))
  dir.create(plotDir, recursive=T)
plotDir <- paste0(plotDir, "/")
tit <- substr(id, 5, nchar(id))
print(plotDir)

##########################################################################################
# Load data

dataBin <- read.table(paste0(fileBaseName, "_bin.bed"), sep="\t", col.names=c("chr", "start", "end", "CN", "obsCov", "expCov"))
dataBin$chr <- factor(dataBin$chr, levels=c(sapply(1:22, function(x) paste0("chr", x)), "chrX", "chrY"))
print(summary(dataBin$obsCov))
dataBin$dist <- abs(dataBin$obsCov - dataBin$expCov)
dataBin$distRel <- dataBin$dist/dataBin$obsCov

dataSegm <- read.table(paste0(fileBaseName, ".bed"), sep="\t", col.names=c("chr", "start", "end", "CN", "obsCov", "expCov"))
dataSegm$dist <- abs(dataSegm$obsCov - dataSegm$expCov)
dataSegm$distRel <- dataSegm$dist/dataSegm$obsCov

print("Finished reading data")


#########################################################################################
# Observed and expected coverage for each bin

# dataBinFlt <- dataBin[!dataBin$chr %in% c("chrX", "chrY", "chrM"), ]
# dataBinFlt$chr <- factor(dataBinFlt$chr, levels=sapply(1:22, function(x) paste0("chr", x)))
# obsByCN = by(dataBinFlt$obsCov, dataBinFlt$CN, median)
# expByCN = by(dataBinFlt$expCov, dataBinFlt$CN, median)
# # print(obsByCN-expByCN)
# # print(diff(obsByCN))
# p <- ggplot(dataBinFlt, aes(x=obsCov, color=as.factor(CN))) + geom_density() + xlim(min(0, min(dataBinFlt$obsCov)), min(450, max(dataBinFlt$obsCov))) + ylim(0, 0.2) +
#   geom_vline(xintercept=as.numeric(obsByCN), color="blue") + geom_vline(xintercept=as.numeric(expByCN)) + ggtitle(id)
# ggsave(paste0(plotDir, id, "_covDistrByCN.png"), p, width=10, height=6, dpi=350)
# rm(dataBinFlt)

dataBin$posBin <- dataBin$start + (dataBin$end - dataBin$start)/2
maxCN <- max(dataBin$CN)

# if (substr(id, nchar(id)-5, nchar(id))=="_truth") {
#   CNs <- unique(dataBin$CN)
#   cnColors <- c("0"="blue", "1"="cyan", "2"="green", "3"="yellow")
#   for (myCN in 0:2) {
#     myCNs = CNs[CNs>myCN & CNs<(myCN+1)]
#     if (length(myCNs)>0) {
#       gainColsPal <- colorRampPalette(c(cnColors[as.character(myCN)], cnColors[as.character(myCN+1)]))
#       gainCols <- gainColsPal(length(myCNs)+2)[2:(length(myCNs)+1)]
#       names(gainCols) <- sapply(myCNs, as.character)
#       cnColors <- c(cnColors, gainCols)
#     }
#   }
#   if (maxCN>3) {
#     gainColsPal <- colorRampPalette(c("orange", "purple4"))
#     gainCols <- gainColsPal(sum(CNs>3))
#     names(gainCols) <- sapply(CNs[CNs>3], as.character)
#     cnColors <- c(cnColors, gainCols)
#   }
# } else {
#   cnColors <- c("0"="blue", "1"="cyan", "2"="green", "3"="yellow", "4"="orange")
#   if (maxCN>4) {
#     gainColsPal <- colorRampPalette(c("orangered", "purple4"))
#     gainCols <- gainColsPal(maxCN-4)
#     names(gainCols) <- sapply(5:maxCN, as.character)
#     cnColors <- c(cnColors, gainCols)
#   }
#   cnColors <- cnColors[sapply(0:maxCN, as.character)]
# }
dataBin$CN <- factor(dataBin$CN)
# 
# ### Plot whole genome
# # print("Starting plotting CN calls for whole genome")
# # p <- ggplot(data=dataBin, aes(x=posBin, y=obsCov)) + geom_point(size=0.5) +
# #   geom_point(aes(x=posBin, y=expCov, colour=CN), size=0.8) + facet_wrap(~chr, scales="free_x") +
# #   ylab("coverage") + ggtitle(paste(tit, ": observed coverage vs CN calls")) + ylim(0, 1500) +
# #   scale_colour_manual(values=cnColors) + guides(colour=guide_legend(override.aes=list(size=5))) +
# #   theme(text=element_text(size=17))
# # ggsave(paste0(plotDir, id, ".png"), p, width=20, height=12, dpi=300)
# # print(paste0("Finished plotting CN calls for whole genome: ", id, ".png"))
# # 
# # #### Plot all genome, zoom in
# # print("Starting plotting CN calls for whole genome (zoom)")
# # p <- ggplot(data=dataBin, aes(x=posBin, y=obsCov)) + geom_point(size=0.5) +
# #   geom_point(aes(x=posBin, y=expCov, colour=CN), size=1) + facet_wrap(~chr, scales="free_x") +
# #   ylab("coverage") + ggtitle(paste(tit, ": observed coverage vs CN calls")) + ylim(0, 500) +
# #   scale_colour_manual(values=cnColors) + guides(colour=guide_legend(override.aes=list(size=5))) +
# #   theme(text=element_text(size=16), axis.text.x=element_blank())
# # ggsave(paste0(plotDir, id, "_zoom.png"), p, width=15, height=10, dpi=300)
# # print(paste0("Finished plotting CN calls for whole genome (zoom): ", id, "_zoom.png"))
# 
# ### Zoom in for single chromosomes
# # Create directories
# plotDirChr <- paste0(plotDir, "CNchr")
# if (!(plotDirChr %in% list.dirs(".")))
#   dir.create(plotDirChr)
# plotDirChr <- paste0(plotDirChr, "/")
# print("Plotting CN calls for single chromosomes (zoom):")
# # Plot each chr
# for (chrID in unique(dataBin$chr)) {
#   p <- ggplot(data=subset(dataBin, chr==chrID), aes(x=posBin, y=obsCov)) + geom_point(size=0.5) +
#     geom_point(aes(x=posBin, y=expCov, colour=CN), size=0.8) + ylim(0, 500) + ylab("coverage") +
#     scale_colour_manual(values=cnColors) + guides(colour=guide_legend(override.aes=list(size=3))) +
#     ggtitle(paste(tit, ",", chrID, ": observed coverage vs CN calls"))
#   ggsave(paste0(plotDirChr, id, "_", chrID, ".png"), p, width=8, height=4, dpi=300)
#   print(paste0("    Plot for ", chrID, ": ", id, "_", chrID, ".png"))
# }
# 
# #### Zoom in for single chromosomes, add percentages colors
# # Read data with perc per variants per bin
# dataBinPerc <- read.table(paste0(fileBaseName, "_bin_perc_all.bed"), sep="\t", col.names=c("chrV", "startV", "endV", "cnA", "cnB", "perc", "chrB", "startB", "endB", "cnCall", "obsCov", "expCov"))
# dataBinPerc$posBin <- dataBinPerc$startB + (dataBinPerc$endB - dataBinPerc$startB)/2
# dataBinPerc$CN <- dataBinPerc$cnA + dataBinPerc$cnB
# dataBinPerc$perc <- factor(dataBinPerc$perc, levels=sort(unique(dataBinPerc$perc)))
# # Create directories
# plotDirChrP <- paste0(plotDir, "CNchrPerc")
# if (!(plotDirChrP %in% list.dirs(".")))
#   dir.create(plotDirChrP)
# plotDirChrP <- paste0(plotDirChrP, "/")
# print(plotDirChrP)
# # Plot each chr
# percColsPal <- colorRampPalette(c("grey60", "black"))
# percColors <- percColsPal(length(levels(dataBinPerc$perc)))
# names(percColors) <- levels(dataBinPerc$perc)
# print("Plotting CN calls with clonal percentages for single chromosomes (zoom):")
# for (chrID in unique(dataBinPerc$chrV)) {
#   p <- ggplot(data=subset(dataBinPerc, chrV==chrID), aes(x=posBin, y=obsCov)) + geom_point(size=0.5, aes(color=perc)) +
#     scale_colour_manual(values=percColors, name="var %") + guides(colour=guide_legend(override.aes=list(size=3))) +
#     geom_point(aes(x=posBin, y=expCov), size=1, color="red") + ylim(0, 500) +
#     ylab("coverage") + ggtitle(paste(tit, ",", chrID, ": observed coverage vs CN call"))
#   ggsave(paste0(plotDirChrP, id, "_clperc_", chrID, ".png"), p, width=8, height=4, dpi=300)
#   print(paste0("    Plot for ", chrID, ": ", id, "_clperc_", chrID, ".png"))
# }
# #        
# # 
# 
#########################################################################################
# Observed and expected coverage for each segment/bin

d_th <- theme(plot.title=element_text(size=12, vjust=1))

### Absolute distance
print("RMSE per segment:")
dataSegm$len <- as.numeric(dataSegm$end - dataSegm$start)
print(rmse <- sqrt(sum(dataSegm$len *(dataSegm$obsCov-dataSegm$expCov)^2)/sum(dataSegm$len)))
print("Plotting stats for absolute distances")
# Histogram for segments
p <- ggplot(dataSegm, aes(x=dist)) + geom_histogram(binwidth=5) +
  xlab("|obsCov - expCov|") + ggtitle(paste(tit, ": absolute distance per segment"))+
  geom_text(aes(max(dist)-170, 70, label=paste0("RMSE =", round(rmse, 2))), size=4)
ggsave(paste0(plotDir, id, "_distSegm_hist.png"), p, width=6, height=4, dpi=200)
# Histogram for segments, fixed xlim, ylim
p <- ggplot(dataSegm, aes(x=dist)) + geom_histogram(binwidth=4)  +
  xlim(0, 1000) + ylim(0, 200) +
  xlab("|obsCov - expCov|") + ggtitle(paste(tit, ": absolute distance per segment"))
ggsave(paste0(plotDir, id, "_distSegm_hist_fixed.png"), p, width=6, height=4, dpi=200)
# Histogram for bins
dataBin$len <- as.numeric(dataBin$end - dataBin$start)
print(rmse <- sqrt(sum(dataBin$len *(dataBin$obsCov-dataBin$expCov)^2)/sum(dataBin$len)))
p <- ggplot(dataBin, aes(x=dist)) + geom_histogram(binwidth=5) +
  xlab("|obsCov - expCov|") + ggtitle(paste(tit, ": absolute distance per bin"))+
  geom_text(aes(max(dist)-170, 70, label=paste0("RMSE =", round(rmse, 2))), size=4)
ggsave(paste0(plotDir, id, "_distBin_hist.png"), p, width=6, height=4, dpi=200)
# Boxplots for bins, per chr
p <- ggplot(dataBin, aes(x=chr, y=dist)) + geom_boxplot() +
  ylab("|obsCov - expCov|") + ggtitle(paste(tit, ": absolute distance per bin"))
ggsave(paste0(plotDir, id, "_distBin_bp.png"), p, width=10, height=4, dpi=200)
p <- ggplot(dataBin, aes(x=chr, y=dist)) + geom_boxplot() + ylim(0, 200) +
  ylab("|obsCov - expCov|") + ggtitle(paste(tit, ": absolute distance per bin, zoom"))
ggsave(paste0(plotDir, id, "_distBin_bp_zoom.png"), p, width=10, height=4, dpi=200)
p <- ggplot(dataBin, aes(x=chr, y=dist)) + geom_boxplot() + ylim(0, 50) +
  ylab("|obsCov - expCov|") + ggtitle(paste(tit, ": absolute distance per bin, zoom"))
ggsave(paste0(plotDir, id, "_distBin_bp_zoom2.png"), p, width=10, height=4, dpi=200)
print("Finished abs dist plots")

### Relative distance
# Histogram for segments
print("Plotting stats for relative distances")
p <- ggplot(dataSegm, aes(x=distRel)) + geom_histogram(binwidth=0.05) +
  xlab("|obsCov - expCov|/obsCov") + ggtitle(paste(tit, ": relative distance per segment"))
ggsave(paste0(plotDir, id, "_distSegmRel_hist.png"), p, width=6, height=4, dpi=200)
p <- ggplot(dataSegm, aes(x=distRel)) + geom_histogram(binwidth=0.05) +
  xlim(0, 3) + ylim(0, 200) +
  xlab("|obsCov - expCov|/obsCov") + ggtitle(paste(tit, ": relative distance per segment"))
ggsave(paste0(plotDir, id, "_distSegmRel_hist_fixed.png"), p, width=6, height=4, dpi=200)
# Histogram for bins
p <- ggplot(dataBin, aes(x=distRel)) + geom_histogram(binwidth=0.05) +
  xlab("|obsCov - expCov|/obsCov") + ggtitle(paste(tit, ": relative distance per bin"))
ggsave(paste0(plotDir, id, "_distBinRel_hist.png"), p, width=6, height=4, dpi=200)
# Boxplots for bins, per chr
p <- ggplot(dataBin, aes(x=chr, y=distRel)) + geom_boxplot() +
  ylab("|obsCov - expCov|/obsCov") + ggtitle(paste(tit, ": relative distance per bin"))
ggsave(paste0(plotDir, id, "_distBinRel_bp.png"), p, width=10, height=4, dpi=200)
p <- ggplot(dataBin, aes(x=chr, y=distRel)) + geom_boxplot() + ylim(0, 5) + 
  ylab("|obsCov - expCov|/obsCov") + ggtitle(paste(tit, ": relative distance per bin, zoom"))
ggsave(paste0(plotDir, id, "_distBinRel_bp_zoom.png"), p, width=10, height=4, dpi=200)
p <- ggplot(dataBin, aes(x=chr, y=distRel)) + geom_boxplot() + ylim(0, 1) + 
  ylab("|obsCov - expCov|/obsCov") + ggtitle(paste(tit, ": relative distance per bin, zoom"))
ggsave(paste0(plotDir, id, "_distBinRel_bp_zoom2.png"), p, width=10, height=4, dpi=200)
print("Finished rel dist plots")
# 
# ### Relative distance for shifting CN
# if (length(list.files(id, paste0(id, "_CNshift_.+[.]bed")))>1) {
#   print("Plotting relative distance for shifting CN")
#   dataCNshift <- adply(c(-2, -1, 1, 2), 1, function(x) {
#     data <- read.table(paste0(fileBaseName, "_CNshift_", as.character(x), ".bed"), sep="\t",
#                        col.names=c("chr", "start", "end", "CN", "obsCov", "expCov"))
#     data$distRel <- abs(data$expCov - data$obsCov)/data$obsCov
#     return(cbind(data, data.frame(CNshift=x)))
#   })[-1]
#   dataCNshift0 <- cbind(dataSegm, data.frame(CNshift=0))
#   dataCNshift <- rbind(dataCNshift, cbind(dataCNshift0[, colnames(dataCNshift)]))
#   dataCNshift$CNshift <- as.factor(dataCNshift$CNshift)
#   p <- ggplot(dataCNshift, aes(x=CNshift, y=distRel)) + geom_boxplot() +
#     ylim(0, 1) + ggtitle(paste(tit, ": relative distance\nfor shifting CN Canvas calls, zoom"))
#   ggsave(paste0(plotDir, id, "_CNshift_distSegmRel.png"), p, width=7, height=4, dpi=200)
# }
# 
# 
# 
##########################################################################
# Compare heterogeneous variants: percentage vs distance

if (paste0(fileBaseName, "_perc.bed") %in% list.files(".", recursive=T)) {
  
  dataPerc <- read.table(paste0(fileBaseName, "_perc.bed"), sep="\t",
                         col.names=c("chrV", "startV", "endV", "cnA", "cnB", "perc", "chrS", "startS", "endS", "cnCall", "obsCov", "expCov"))
  dataPerc$perc <- as.factor(dataPerc$perc)
  dataPerc$dist <- abs(dataPerc$obsCov - dataPerc$expCov)
  dataPerc$distRel <- dataPerc$dist/dataPerc$obsCov
  print("Plotting variants percentage vs distances")
  
  p <- ggplot(dataPerc, aes(x=perc, y=dist)) + geom_boxplot() +
    ylim(0, 100) + ylab("|obsCov - expCov|") + xlab("variant %") + d_th +
    ggtitle(paste(tit, ": absolute distance vs % for each variant, zoom"))
  ggsave(paste0(plotDir, id, "_perc_dist_bp.png"), p, width=7, height=4, dpi=200)
  p <- ggplot(dataPerc, aes(x=perc, y=distRel)) + geom_boxplot() +
    ylim(0, 1) + ylab("|obsCov - expCov|/obsCov") +  xlab("variant %") + d_th +
    ggtitle(paste(tit, ": relative distance vs % for each variant, zoom"))
  ggsave(paste0(plotDir, id, "_perc_distRel_bp.png"), p, width=7, height=4, dpi=200)
  
  p <- ggplot(dataPerc, aes(color=perc, x=dist)) + geom_density() +
    xlim(0, 100) + xlab("|obsCov - expCov|") + scale_color_discrete(name="variant %") + d_th +
    ggtitle(paste(tit, ": absolute distance vs % for each variant, zoom"))
  ggsave(paste0(plotDir, id, "_perc_dist_dens.png"), p, width=7, height=5, dpi=200)
  p <- ggplot(dataPerc, aes(color=perc, x=distRel)) + geom_density() +
    xlim(0, 1) + xlab("|obsCov - expCov|/obsCov") + scale_color_discrete(name="variant %") + d_th +
    ggtitle(paste(tit, ": relative distance vs % for each variant, zoom"))
  ggsave(paste0(plotDir, id, "_perc_distRel_dens.png"), p, width=7, height=5, dpi=200)
  
  print("Finished plotting\n\n")
  
#   dataPerc$het <- as.numeric(as.character(dataPerc$perc)) <= 0.4
#   print(t.test(distRel ~ het, data=subset(dataPerc, distRel!=Inf), na.rm=T))
#   print(t.test(dist ~ het, data=dataPerc, na.rm=T))

  # tapply(dataPerc$distRel,dataPerc$perc, function(x) summary(x[x!=Inf]))
  # with(dataPerc[dataPerc$perc==1.0, ], print(summary(cnCall-cnA-cnB)))
  # with(dataPerc[dataPerc$perc==0.6, ], print(summary(cnCall-cnA-cnB)))
  # with(dataPerc[dataPerc$perc==0.4, ], print(summary(cnCall-cnA-cnB)))
  
}



#############################################################################
# dataPart <- read.table("cov_simGeL004_n1_m-1_v-1_truth_partitioned.bed", sep="\t", col.names=c("chr", "start", "end", "cov", "segm"))
# 
# 
# 
# dataPart <- read.table("//ukch-prd-isi01/scratch_tmp/users/ccolombo/Canvas/simNorm/Analysis/TempCNV_G1_P1/G1_P1", sep="\t", col.names=c("chr", "start", "end", "cov", "segm"))
# tit <- paste0("coverage for normal sample\n",
#               paste(names(summary(dataPart$cov)), collapse=" "), "\n",
#               paste(summary(dataPart$cov), collapse="    "))
# p <- ggplot(dataPart, aes(x=cov)) + geom_histogram(binwidth=5) + ggtitle(tit)
# ggsave("./Plots/norm_covHist.png", p, width=6, height=5, dpi=200)
# p <- ggplot(dataPart, aes(x=cov)) + geom_histogram(binwidth=5) + xlim(0, 200) + ggtitle(tit)
# ggsave("./Plots/norm_covHist.png", p, width=6, height=5, dpi=200)
# print(median(dataPart$cov))
# print(mean(dataPart$cov))
# print(summary(dataPart$cov))


# Cfr expected coverage vs observed coverage distributions
# FIXME: for normal simulation!!!!!!!!!
# cn=0:10
# expCov=57*cn + 24
# cndata = data.frame(cn=cn+1, expCov=expCov)
# p <- ggplot(dataBin, aes(x=CN, y=obsCov)) + geom_point(position = position_jitter(width = 0.2)) +
#   geom_point(data=cndata, aes(x=cn, y=expCov), color="red", size=4)
# dataSegm$CN <- factor(as.integer(dataSegm$CN), levels=0:10)
# p <- ggplot(dataSegm, aes(x=CN, y=obsCov)) + geom_boxplot() +
#   geom_point(data=cndata, aes(x=cn, y=expCov), color="red", size=4)
# p <- ggplot(dataBin, aes(x=as.factor(CN), y=obsCov)) + geom_boxplot() + geom_point(cndata, aes(x=cn, y=expCov), color="red", size=4)
# p <- ggplot(dataBin[dataBin$obsCov<1000 & dataBin$CN!=0,], aes(x=obsCov, color=as.factor(CN))) + geom_density()