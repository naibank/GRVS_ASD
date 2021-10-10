setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

rm(list = ls())
library(data.table)
source("burdenRequireFunctions.R")

file.in <- read.delim("../CNVs/CNVs.SSC.freq_1percent.HQR.hg19.tsv", stringsAsFactors = F)
file.in <- file.in[file.in$SIZE <= 10000, ]
dt.out <- data.frame()
for(cnv in c("DEL", "DUP")){
  
  total.count <- data.frame(table(file.in$sample[file.in$SVTYPE == cnv]))
  names(total.count) <- c("Sample", "TotalRare")
  
  feat <- c("utr_5", "utr_3", "An_promoter", "DDD_promoter",
            "lnc_promoter", "brain_h3k27ac", "tad", "ctcf", "brain_enh")        
  dt <- data.frame(fread("../dataNC/cnv.annotated.tsv"), stringsAsFactors = F)
  dt <- dt[dt$CHROM %in% paste0("chr", 1:22), ]
  dt <- dt[dt$ts_overlap == 0, ]
  dt <- dt[dt$SIZE <= 10000 & dt$SVTYPE == cnv, ]
  feat <- feat[feat %in% names(dt)]
  dt <- aggregate(.~sample, dt[, c("sample", feat)], sum, na.rm = T, na.action = NULL)
  dt <- merge(dt, total.count, by.x = "sample", by.y = "Sample", all = T)
  dt[is.na(dt)] <- 0
  names(dt)[-1] <- paste(names(dt)[-1], cnv, sep="_")
  if(cnv == "DEL"){
    dt.out <- dt
  }else{
    dt.out <- merge(dt.out, dt, by = "sample", all = T)
    dt.out[is.na(dt.out)] <- 0
  }
}

write.table(dt.out, "../dataNC/small.cnv.matrix.tsv", sep="\t", row.names=F, quote=F, col.names=T)
