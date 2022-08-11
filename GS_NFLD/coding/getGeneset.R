setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
library(data.table)
type <- c("rare", "denovo")
load("../../new_set/gsNFLD_2022.RData")

iselement <- function(gs, enzid){
  return(enzid %in% gs)
}

getSummary <- function(dt, gs){
  totalVar <- aggregate(entrez_id ~ Sample, dt, length)
  names(totalVar) <- c("Sample", "Total")
  otherVar <- aggregate(.~Sample, dt[, c("Sample", names(gs))], sum, na.rm = T) 
  return(merge(totalVar, otherVar, by = "Sample", all = T))
}

### get gene-set matrix ###
mpc <- read.delim("../../../mpc/NFLD_main_snv.rare.missenses.tsv.annovar.out_rev27.6_hg19.tsv", stringsAsFactors = F)
mpc <- mpc[which(mpc$MPC_score > 2), ]
mpc <- unique(paste(mpc$chr, mpc$start, mpc$ref_allele, mpc$alt_allele, sep="#"))

for(th.type in type){
  d <- data.frame(fread(sprintf("../../../../main/data/SNVs/calls/%s/all.%s.variants.Jun2022.tsv", th.type, th.type)), stringsAsFactors = F)
  d <- d[d$chr %in% paste0("chr", 1:22), ]
  d$varid <- paste(d$chr, d$start, d$ref_allele, d$alt_allele, sep="#")
  
  
  #lof rare = 5077 + 26 = 5103
  #mis2 rare = 20382 + 80= 20462
  
  total.rare.variants <- aggregate(entrez_id ~ Sample, d, length)
  names(total.rare.variants) <- c("Sample", "TotalRare")
  lof <- d[which(d$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain") | d$typeseq_priority == "splicing"), ]
  mis.t1 <- d[which(d$effect_priority == "nonsynonymous SNV"), ]
  mis.t2 <- mis.t1[rowSums(data.frame(
    mis.t1$sift_score < 0.05, 
    mis.t1$polyphen_score >= 0.90,
    mis.t1$ma_score >= 1.90, 
    mis.t1$phylopPMam_avg >= 2.30,
    mis.t1$phylopVert100_avg >= 4.0, 
    mis.t1$CADD_phred >= 15,
    mis.t1$mt_score >= 0.5), na.rm = T) >= 5, ]
  
  mis.t1 <- mis.t1[rowSums(data.frame(
    mis.t1$sift_score < 0.05, 
    mis.t1$polyphen_score >= 0.90,
    mis.t1$ma_score >= 1.90, 
    mis.t1$phylopPMam_avg >= 2.30,
    mis.t1$phylopVert100_avg >= 4.0, 
    mis.t1$CADD_phred >= 15,
    mis.t1$mt_score >= 0.5), na.rm = T) < 5, ]
  
  mis.mpc2 <- d[which(d$effect_priority == "nonsynonymous SNV" & d$varid %in% mpc), ]
  mis.mpc2 <- data.frame(table(mis.mpc2$Sample))
  names(mis.mpc2) <- c("Sample", "MisMPC_morethan2")
  
  lof[, names(gsNFLD)] <- sapply(gsNFLD, iselement, lof$entrez_id)
  mis.t1[, names(gsNFLD)] <- sapply(gsNFLD, iselement, mis.t1$entrez_id)
  mis.t2[, names(gsNFLD)] <- sapply(gsNFLD, iselement, mis.t2$entrez_id)

  # >   nrow(lof)/length(unique(d$Sample))
  # [1] 16.92333
  # >   nrow(mis.t1)/length(unique(d$Sample))
  # [1] 236.9067
  # >   nrow(mis.t2)/length(unique(d$Sample))
  # [1] 67.94
  # frameshift deletion frameshift insertion             stopgain 
  # 5.493333             2.560000             5.280000 
  # exonic exonic;splicing        splicing 
  # 13.31666667      0.01333333      3.59333333 
  
  write.table(lof, sprintf("../../../../main/data/SNVs/calls/%s/lof.variants.Jun2022.tsv", th.type), sep="\t", row.names=F, quote=F, col.names=T)
  write.table(mis.t1, sprintf("../../../../main/data/SNVs/calls/%s/ms1.variants.Jun2022.tsv", th.type), sep="\t", row.names=F, quote=F, col.names=T)
  write.table(mis.t2, sprintf("../../../../main/data/SNVs/calls/%s/ms2.variants.Jun2022.tsv", th.type), sep="\t", row.names=F, quote=F, col.names=T)
  
  lof.feat <- getSummary(lof, gsNFLD)
  mis.t1.feat <- getSummary(mis.t1, gsNFLD)
  mis.t2.feat <- getSummary(mis.t2, gsNFLD)
  
  names(lof.feat)[-1] <- paste0(names(lof.feat)[-1], "_lof")
  names(mis.t1.feat)[-1] <- paste0(names(mis.t1.feat)[-1], "_tier1_ms")
  names(mis.t2.feat)[-1] <- paste0(names(mis.t2.feat)[-1], "_tier2_ms")
  
  dt.out <- merge(total.rare.variants, lof.feat, by = "Sample", all = T)
  dt.out <- merge(dt.out, mis.t1.feat, by = "Sample", all = T)
  dt.out <- merge(dt.out, mis.t2.feat, by = "Sample", all = T)
  dt.out <- merge(dt.out, mis.mpc2, by = "Sample", all = T)
  dt.out[is.na(dt.out)] <- 0
  
  write.table(dt.out, sprintf("%s.gs.matrix.tsv", th.type), sep="\t", row.names=F, quote=F, col.names=T)
}
