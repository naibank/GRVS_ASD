setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
library(data.table)
type <- c("rare")


### get gene-set matrix ###

for(th.type in type){
  # d <- data.frame(fread(sprintf("calls/%s/all.%s.variants.OCT2019.tsv", th.type, th.type)), stringsAsFactors = F)
  d <- fread("../SNVs/SSC_rare1pct_SNVs_INDELs_QCed.tsv", data.table = F)
  d <- d[d$chr %in% paste0("chr", 1:22), ]
  
  #lof rare = 5077 + 26 = 5103
  #mis2 rare = 20382 + 80= 20462
  
  total.rare.variants <- aggregate(entrez_id ~ sample, d, length)
  names(total.rare.variants) <- c("sample", "TotalRare")
  lof <- d[which(d$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | d$typeseq_priority == "splicing"), ]
  mis.t1 <- d[which(d$effect_priority == "nonsynonymous SNV"), ]
  mis.t2 <- mis.t1[rowSums(data.frame(
    mis.t1$sift_score < 0.05, 
    mis.t1$polyphen_score >= 0.90,
    mis.t1$ma_score >= 1.90, 
    mis.t1$phylopMam_avg >= 2.30,
    mis.t1$phylopVert100_avg >= 4.0, 
    mis.t1$CADD_phred >= 15,
    mis.t1$mt_score >= 0.5), na.rm = T) >= 5, ]
  
  mis.t1 <- mis.t1[rowSums(data.frame(
    mis.t1$sift_score < 0.05, 
    mis.t1$polyphen_score >= 0.90,
    mis.t1$ma_score >= 1.90, 
    mis.t1$phylopMam_avg >= 2.30,
    mis.t1$phylopVert100_avg >= 4.0, 
    mis.t1$CADD_phred >= 15,
    mis.t1$mt_score >= 0.5), na.rm = T) < 5, ]

  lof[, names(gsNFLD)] <- sapply(gsNFLD, iselement, lof$entrez_id)
  mis.t1[, names(gsNFLD)] <- sapply(gsNFLD, iselement, mis.t1$entrez_id)
  mis.t2[, names(gsNFLD)] <- sapply(gsNFLD, iselement, mis.t2$entrez_id)

  # > nrow(lof)/length(unique(d$sample))
  # [1] 41.53992
  # > nrow(mis.t1)/length(unique(d$sample))
  # [1] 343.6912
  # > nrow(mis.t2)/length(unique(d$sample))
  # [1] 74.71849
  table(lof$effect_priority)/length(unique(d$sample))
  table(lof$typeseq_priority)/length(unique(d$sample))
  
  write.table(lof, sprintf("../SNVs/SSC_%s_1pct_lof.variants.tsv", th.type), sep="\t", row.names=F, quote=F, col.names=T)
  write.table(mis.t1, sprintf("../SNVs/SSC_%s_1pct_ms1.variants.tsv", th.type), sep="\t", row.names=F, quote=F, col.names=T)
  write.table(mis.t2, sprintf("../SNVs/SSC_%s_1pct_ms2.variants.tsv", th.type), sep="\t", row.names=F, quote=F, col.names=T)
  
  lof.feat <- getSummary(lof, gsNFLD)
  mis.t1.feat <- getSummary(mis.t1, gsNFLD)
  mis.t2.feat <- getSummary(mis.t2, gsNFLD)
  
  names(lof.feat)[-1] <- paste0(names(lof.feat)[-1], "_lof")
  names(mis.t1.feat)[-1] <- paste0(names(mis.t1.feat)[-1], "_tier1_ms")
  names(mis.t2.feat)[-1] <- paste0(names(mis.t2.feat)[-1], "_tier2_ms")
  
  dt.out <- merge(total.rare.variants, lof.feat, by = "sample", all = T)
  dt.out <- merge(dt.out, mis.t1.feat, by = "sample", all = T)
  dt.out <- merge(dt.out, mis.t2.feat, by = "sample", all = T)
  dt.out[is.na(dt.out)] <- 0
  
  write.table(dt.out, sprintf("../SNVs/SSC_%s.gs.matrix.tsv", th.type), sep="\t", row.names=F, quote=F, col.names=T)
}
