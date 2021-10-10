setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
library(data.table)
load("../../main/requiredData/gsNFLD.RData")

# dt <- data.frame(readxl::read_excel("../denovo SNVs from 2018 An et al/aat6576_Table-S2.xlsx"), stringsAsFactors = F)
# denovo.hg38 <- fread("../SNVs/SSC_denovo/f_all_20200702.formatted.ann.final.ssc-dn.aid.txt", data.table = F)
dt <- read.delim("../SNVs/SSC_denovo/hg19.bed", stringsAsFactors = F, header = F)
ssc_rare <- read.delim("../SNVs/SSC_rare.gs.matrix.tsv", stringsAsFactors = F)
names(dt) <- c("chr", "start", "end", "sample", "varid")

# id.map <- read.delim("../../../MSSNGExpansion/data/SampleInfo/nygc_sfari_id_sample_map.data", sep=",", stringsAsFactors = F)
# id.map <- id.map[id.map$Sample.ID %in% ssc_rare$sample, ]
# dt <- dt[dt$SampleID %in% id.map$SFARI.ID, c(1:14)]
# dt$End <- as.numeric(dt$Pos) + 1
# dt$ID <- paste(dt$Chr, dt$Pos, dt$End, sep="#")
# write.table(dt[, c("Chr", "Pos", "End", "ID")], "../denovo SNVs from 2018 An et al/tobelifted_tohg19.txt", 
#             sep="\t", row.names=F, quote=F, col.names=F)
# hg19 <- read.delim("../denovo SNVs from 2018 An et al/hglft_genome_228b_5eea00.bed", stringsAsFactors = F, header = F)
# names(hg19) <- c("chr_hg19", "start_hg19", "end_hg19", "ID")
# 
# dt <- merge(dt, hg19, by = "ID", all.y = T)
# dt <- merge(dt, id.map, by.x = "SampleID", by.y = "SFARI.ID", all = F)
dt$varid <- paste(dt$chr, dt$start+1, dt$sample, sep="#")
# allvar <- fread("../SNVs/all.variants.tsv", data.table = F)
lof <- read.delim("../SNVs/SSC_rare_1pct_lof.variants.tsv", stringsAsFactors = F)
ms1 <- read.delim("../SNVs/SSC_rare_1pct_ms1.variants.tsv", stringsAsFactors = F)
ms2 <- read.delim("../SNVs/SSC_rare_1pct_ms2.variants.tsv", stringsAsFactors = F)

lof$ID <- paste(lof$chr, lof$start, lof$sample, sep = "#")
ms1$ID <- paste(ms1$chr, ms1$start, ms1$sample, sep = "#")
ms2$ID <- paste(ms2$chr, ms2$start, ms2$sample, sep = "#")

lof$denovo <- lof$ID %in% dt$varid
ms1$denovo <- ms1$ID %in% dt$varid
ms2$denovo <- ms2$ID %in% dt$varid

table(ms2$denovo)

####
iselement <- function(gs, enzid){
  return(enzid %in% gs)
}

getSummary <- function(dt, gs){
  totalVar <- aggregate(entrez_id ~ sample, dt, length)
  names(totalVar) <- c("sample", "Total")
  otherVar <- aggregate(.~sample, dt[, c("sample", names(gs))], sum, na.rm = T) 
  return(merge(totalVar, otherVar, by = "sample", all = T))
}

load("../../main/requiredData/gsNFLD.RData")

lof.feat <- getSummary(lof[lof$denovo, ], gsNFLD)
mis.t1.feat <- getSummary(ms1[ms1$denovo, ], gsNFLD)
mis.t2.feat <- getSummary(ms2[ms2$denovo, ], gsNFLD)

names(lof.feat)[-1] <- paste0(names(lof.feat)[-1], "_lof")
names(mis.t1.feat)[-1] <- paste0(names(mis.t1.feat)[-1], "_tier1_ms")
names(mis.t2.feat)[-1] <- paste0(names(mis.t2.feat)[-1], "_tier2_ms")

dt.out <- merge(ssc_rare[, 1:2], lof.feat, by = "sample", all = T)
dt.out <- merge(dt.out, mis.t1.feat, by = "sample", all = T)
dt.out <- merge(dt.out, mis.t2.feat, by = "sample", all = T)
dt.out[is.na(dt.out)] <- 0

write.table(dt.out, "../SNVs/SSC_denovo.gs.matrix.tsv", sep="\t", row.names=F, quote=F, col.names=T)
nfld <- read.delim("../../main/data/SNVs/denovo.gs.matrix.tsv", stringsAsFactors = F)

summary(dt.out$Total_lof)
summary(nfld$Total_lof)