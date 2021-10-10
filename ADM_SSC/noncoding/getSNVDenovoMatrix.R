setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

rm(list = ls())
library(data.table)

source("burdenRequireFunctions.R")

total.count <- read.delim("../SNVs/SSC_rare.gs.matrix.tsv", stringsAsFactors = F)
total.count <- total.count[, c("sample", "TotalRare")]

feat <- c("An_promoter", "DDD_promoter", "lnc_promoter", "tss_intolerance", "fetal_brain_intolerance", 
          "An_promoter_tier2", "DDD_promoter_tier2", "lnc_promoter_tier2", "tss_intolerance_tier2", 
          "fetal_brain_intolerance_tier2", "splicing_tier1", 
          "splicing_tier2", "utr3_tier1", "utr3_tier2", "utr5", "deepsea_tier1", "deepsea_tier2", "brain_h3k27ac", "brain_h3k27ac_2gerp", 
          "tad_tier1", "tad_tier2", "ctcf_tier1", "ctcf_tier2", "brain_enh_tier1", "brain_enh_tier2")       
dt <- data.frame(fread("../dataNC/snv.rare.annotated.tsv"), stringsAsFactors = F)
denovo <- read.delim("../SNVs/SSC_denovo/hg19.bed", stringsAsFactors = F, header = F)
names(denovo) <- c("chr", "start", "end", "sample", "varid")
denovo$varid <- paste(denovo$chr, denovo$start+1, denovo$sample, sep="#")
dt$varid <- paste(dt$chr, dt$start, dt$sample, sep="#")
dt <- dt[dt$varid %in% denovo$varid, ]

dt <- aggregate(.~sample, dt[, c("sample", feat)], sum, na.rm = T, na.action = NULL)

summary(dt$fetal_brain_intolerance)
dt <- merge(dt, total.count, by = "sample", all = T)
dt[is.na(dt)] <- 0

write.table(dt, "../dataNC/denovo.snv.matrix.tsv", sep="\t", row.names=F, quote=F, col.names=T)

