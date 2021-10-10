setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
library(ggplot2)
p <- 0.1

check <- c("SSC06086", "SSC09081", "SSC09084", "SSC09086", "SSC09478")
check.fam <- c("12175", "13069", "13227", "14484")
n442 <- read.delim("../SNVs/n442.SSCaff_unaffsib_samplelist.tsv", stringsAsFactors = F)

# for(p in c(1, 0.5, 0.1, 0.05, 0.01)){
# for(var in c("duplication", "deletion", "lof", "tier1_ms", "tier2_ms")){
cnvs <- read.delim("../CNVs/SSC.gs.matrix.tsv", stringsAsFactors = F)
nccnvs <- read.delim("../dataNC/cnv.matrix.tsv", stringsAsFactors = F)

scnvs <- read.delim("../CNVs/smaller.SSC.gs.matrix.tsv", stringsAsFactors = F)
sncnvs <- read.delim("../dataNC/small.cnv.matrix.tsv", stringsAsFactors = F)

nccnvs <- nccnvs[nccnvs$sample %in% c(cnvs$sample, scnvs$sample), ]
sncnvs <- sncnvs[sncnvs$sample %in% c(cnvs$sample, scnvs$sample), ]

rSNV <- read.delim("../SNVs/SSC_rare.gs.matrix.tsv", stringsAsFactors = F)
rSNV.outliers <- rSNV$sample[rSNV$TotalRare > 3*IQR(rSNV$TotalRare) + quantile(rSNV$TotalRare, 0.75)]
eth <- read.delim("../SSC_data/Target.admixture.tagged.tsv", stringsAsFactors = F)
rSNV.outliers <- rSNV.outliers[rSNV.outliers %in% eth$IND[eth$ETH_TAG == "EUR"]]

dSNV <- read.delim("../SNVs/SSC_denovo.gs.matrix.tsv", stringsAsFactors = F)

nrsnv <- read.delim("../dataNC/rare.snv.matrix.tsv", stringsAsFactors = F)
ndsnv <- read.delim("../dataNC/denovo.snv.matrix.tsv", stringsAsFactors = F)

cnv.outliers <- readLines("CNVOutlier.txt")

outliers <- unique(c(rSNV.outliers, cnv.outliers))

# writeLines(outliers, "exceed_calls_outliers.txt")
scnvs <- scnvs[!scnvs$sample %in% outliers, ]
cnvs <- cnvs[!cnvs$sample %in% outliers, ]
rSNV <- rSNV[!rSNV$sample %in% outliers, ]
dSNV <- dSNV[!dSNV$sample %in% outliers, ]
nrsnv <- nrsnv[!nrsnv$sample %in% outliers, ]
ndsnv <- ndsnv[!ndsnv$sample %in% outliers, ]
nccnvs <- nccnvs[!nccnvs$sample %in% outliers, ]
sncnvs <- sncnvs[!sncnvs$sample %in% outliers, ]

coeff <- read.delim("../../main/scripts/ADM/coeff.ADM.from.all.samples.tsv", stringsAsFactors = F)
# coeff <- coeff[coeff$variant %in% c("coding_denovo_SNVs", "coding_rare_CNVs", "small_coding_rare_CNVs", "coding_rare_SNVs"), ]
coeff <- coeff[coeff$pvalue < p, ]
# coeff <- coeff[-grep("tier1_ms", coeff$set), ]
dt.out <- data.frame("sample" = cnvs$sample, "dummy1" = 1, stringsAsFactors = F)

for(i in 1:nrow(coeff)){
  if(coeff$variant[i] == "coding_rare_CNVs"){
    score <- cnvs
  }else if(coeff$variant[i] == "coding_denovo_SNVs"){
    score <- dSNV
  }else if(coeff$variant[i] == "coding_rare_SNVs"){
    score <- rSNV
  }else if(coeff$variant[i] == "noncoding_denovo_SNVs"){
    score <- ndsnv
  }else if(coeff$variant[i] == "noncoding_rare_CNVs"){
    score <- nccnvs
  }else if(coeff$variant[i] == "noncoding_rare_SNVs"){
    score <- nrsnv
  }else if(coeff$variant[i] == "small_coding_rare_CNVs"){
    score <- scnvs
  }else if(coeff$variant[i] == "small_noncoding_rare_CNVs"){
    score <- sncnvs
  }
  
  score <- score[, c("sample", coeff$set[i])]
  score[, sprintf("%s_%s", coeff$variant[i], i)] <- score[, 2]*coeff$coeff[i]
  dt.out <- merge(dt.out, score[, c(1,3)], by = "sample", all.x = T)
}


cnvs.col <- grep("coding_rare_CNVs", names(dt.out))
scnvs.col <- grep("small_coding_rare_CNVs", names(dt.out))
dsnvs.col <- grep("coding_denovo_SNVs", names(dt.out))
rsnvs.col <- grep("coding_rare_SNVs", names(dt.out))
nccnvs.col <- grep("noncoding_rare_CNVs", names(dt.out))
sncnvs.col <- grep("small_noncoding_rare_CNVs", names(dt.out))
ndsnvs.col <- grep("noncoding_denovo_SNVs", names(dt.out))
nrsnvs.col <- grep("noncoding_rare_SNVs", names(dt.out))

if(length(cnvs.col) > 1){
  dt.out$CNVs_score <- rowSums(dt.out[, cnvs.col])
}else{
  dt.out$CNVs_score <- dt.out[, cnvs.col]
}

if(length(scnvs.col) > 1){
  dt.out$smallCNVs_score <- rowSums(dt.out[, scnvs.col])
}else{
  dt.out$smallCNVs_score <- dt.out[, scnvs.col]
}

if(length(dsnvs.col) > 1){
  dt.out$dSNVs_score <- rowSums(dt.out[, dsnvs.col])
}else{
  dt.out$dSNVs_score <- dt.out[, dsnvs.col]
}

if(length(rsnvs.col) > 1){
  dt.out$rSNVs_score <- rowSums(dt.out[, rsnvs.col])
}else{
  dt.out$rSNVs_score <- dt.out[, rsnvs.col]
}

#### nc
if(length(nccnvs.col) > 1){
  dt.out$noncodeCNVs_score <- rowSums(dt.out[, nccnvs.col])
}else{
  dt.out$noncodeCNVs_score <- dt.out[, nccnvs.col]
}

if(length(sncnvs.col) > 1){
  dt.out$smallnoncodeCNVs_score <- rowSums(dt.out[, sncnvs.col])
}else{
  dt.out$smallnoncodeCNVs_score <- dt.out[, sncnvs.col]
}

if(length(ndsnvs.col) > 1){
  dt.out$noncodedSNVs_score <- rowSums(dt.out[, ndsnvs.col])
}else{
  dt.out$noncodedSNVs_score <- dt.out[, ndsnvs.col]
}

if(length(nrsnvs.col) > 1){
  dt.out$noncoderSNVs_score <- rowSums(dt.out[, nrsnvs.col])
}else{
  dt.out$noncoderSNVs_score <- dt.out[, nrsnvs.col]
}

dt.out <- dt.out[, c("sample", "CNVs_score", "smallCNVs_score", "dSNVs_score", "rSNVs_score",
                     "noncodeCNVs_score", "smallnoncodeCNVs_score", "noncodedSNVs_score", "noncoderSNVs_score")]
dt.out$totalScore <- rowSums(dt.out[, -1])

meta <- read.delim("../SSC_data/SSC_meta_data.tsv", stringsAsFactors = F)
dt.out <- merge(dt.out, meta[, -1], by.x = "sample", by.y = "Sample.ID", all.x = T)
dt.out$adm.dysmorphic <- as.character(dt.out$adm.dysmorphic)
dt.out$adm.dysmorphic[is.na(dt.out$adm.dysmorphic)] <- "unaff"
dt.out$adm.dysmorphic <- factor(dt.out$adm.dysmorphic, levels = c("unaff", "nondys", "dys"))
dt.out <- dt.out[dt.out$sample %in% n442$Sample.ID, ]
write.table(dt.out, sprintf("ssc.grs.p%s.tsv", p), sep="\t", row.names=F, quote=F, col.names=T)
# formatC(wilcox.test(score.out$total_risk[score.out$ADM == "1"],
#                     score.out$total_risk[score.out$ADM == "0"], "greater")$p.value, 
#         format = "e", digits = 1)
if(length(dt.out$CNVs_score) > 0)
  p1.p <- formatC(wilcox.test(dt.out$CNVs_score[dt.out$adm.dysmorphic == "dys"], 
            dt.out$CNVs_score[dt.out$adm.dysmorphic == "nondys"], alternative = "greater")$p.value, format = "e",  digits = 1)
if(length(dt.out$smallCNVs_score) > 0)
  p0.p <- formatC(wilcox.test(dt.out$smallCNVs_score[dt.out$adm.dysmorphic == "dys"], 
                            dt.out$smallCNVs_score[dt.out$adm.dysmorphic == "nondys"], alternative = "greater")$p.value, format = "e", digits = 1)
if(length(dt.out$rSNVs_score) > 0)
  p2.p <- formatC(wilcox.test(dt.out$rSNVs_score[dt.out$adm.dysmorphic == "dys"], 
                          dt.out$rSNVs_score[dt.out$adm.dysmorphic == "nondys"], alternative = "greater")$p.value, format = "e", digits = 1)
if(length(dt.out$dSNVs_score) > 0)
  p3.p <- formatC(wilcox.test(dt.out$dSNVs_score[dt.out$adm.dysmorphic == "dys"], 
                          dt.out$dSNVs_score[dt.out$adm.dysmorphic == "nondys"], alternative = "greater")$p.value, format = "e", digits = 1)

p4.p <- formatC(wilcox.test(dt.out$totalScore[dt.out$adm.dysmorphic == "dys"], 
                          dt.out$totalScore[dt.out$adm.dysmorphic == "nondys"], alternative = "greater")$p.value, format = "e", digits = 1)
p4.unaff.p <- formatC(wilcox.test(dt.out$totalScore[dt.out$adm.dysmorphic == "dys"], 
                            dt.out$totalScore[dt.out$adm.dysmorphic == "unaff"], alternative = "greater")$p.value, format = "e", digits = 1)
###nc
if(length(dt.out$noncodeCNVs_score) > 0)
  p5.p <- formatC(wilcox.test(dt.out$noncodeCNVs_score[dt.out$adm.dysmorphic == "dys"], 
                          dt.out$noncodeCNVs_score[dt.out$adm.dysmorphic == "nondys"], alternative = "greater")$p.value, format = "e", digits = 1)
if(length(dt.out$smallnoncodeCNVs_score) > 0)
  p6.p <- formatC(wilcox.test(dt.out$smallnoncodeCNVs_score[dt.out$adm.dysmorphic == "dys"], 
                            dt.out$smallnoncodeCNVs_score[dt.out$adm.dysmorphic == "nondys"], alternative = "greater")$p.value, format = "e", digits = 1)
if(length(dt.out$noncoderSNVs_score) > 0)
  p7.p <- formatC(wilcox.test(dt.out$noncoderSNVs_score[dt.out$adm.dysmorphic == "dys"], 
                            dt.out$noncoderSNVs_score[dt.out$adm.dysmorphic == "nondys"], alternative = "greater")$p.value, format = "e", digits = 1)
if(length(dt.out$noncodedSNVs_score) > 0)
  p8.p <- formatC(wilcox.test(dt.out$noncodedSNVs_score[dt.out$adm.dysmorphic == "dys"], 
                            dt.out$noncodedSNVs_score[dt.out$adm.dysmorphic == "nondys"], alternative = "greater")$p.value, format = "e", digits = 1)

perc.rank <- function(x) trunc(rank(x))/length(x)

# convert to percentile
if(length(dt.out$CNVs_score) > 0)
  dt.out$CNVs_score <- perc.rank(dt.out$CNVs_score)
if(length(dt.out$smallCNVs_score) > 0)
  dt.out$smallCNVs_score <- perc.rank(dt.out$smallCNVs_score)
if(length(dt.out$rSNVs_score) > 0)
  dt.out$rSNVs_score <- perc.rank(dt.out$rSNVs_score)
if(length(dt.out$dSNVs_score) > 0)
  dt.out$dSNVs_score <- perc.rank(dt.out$dSNVs_score)
###nc
if(length(dt.out$noncodeCNVs_score) > 0)
dt.out$noncodeCNVs_score <- perc.rank(dt.out$noncodeCNVs_score)
if(length(dt.out$smallnoncodeCNVs_score) > 0)
  dt.out$smallnoncodeCNVs_score <- perc.rank(dt.out$smallnoncodeCNVs_score)
if(length(dt.out$noncoderSNVs_score) > 0)
  dt.out$noncoderSNVs_score <- perc.rank(dt.out$noncoderSNVs_score)
if(length(dt.out$noncodedSNVs_score) > 0)
  dt.out$noncodedSNVs_score <- perc.rank(dt.out$noncodedSNVs_score)
dt.out$totalScore <- perc.rank(dt.out$totalScore)

p0 <- ggplot()
p1 <- ggplot()
p2 <- ggplot()
p3 <- ggplot()
p4 <- ggplot()
p5 <- ggplot()
p6 <- ggplot()
p7 <- ggplot()
p8 <- ggplot()

# dt.out$adm.dysmorphic <- factor(dt.out$adm.dysmorphic, levels = c("nondys", "dys"))
# dt.out$totalScore <- scale(dt.out$totalScore)
if(length(dt.out$CNVs_score) > 0)
    p1 <- ggplot(dt.out, aes(x = adm.dysmorphic, y = CNVs_score, fill = adm.dysmorphic)) + geom_violin() +
  annotate("text", label = sprintf("p=%s", p1.p), x = c(1.5), y = c(1.12), cex = 5) +
  annotate("line", x = c(1, 2), y = c(1.05, 1.05), group = c("A", "A")) +
  geom_boxplot(fill = "white", width = .2, alpha=.6) + xlab("") + theme_classic() +     
  theme(legend.position = "none", axis.title.x = element_blank()) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + 
  ylab("Coding CNVs")

if(length(dt.out$smallCNVs_score) > 0)
  p0 <- ggplot(dt.out, aes(x = adm.dysmorphic, y = smallCNVs_score, fill = adm.dysmorphic)) + geom_violin() +
  annotate("text", label = sprintf("p=%s", p0.p), x = c(1.5), y = c(1.12), cex = 5) +
  annotate("line", x = c(1, 2), y = c(1.05, 1.05), group = c("A", "A")) +
  geom_boxplot(fill = "white", width = .2, alpha=.6) + theme_classic() +     
  theme(legend.position = "none", axis.title.x = element_blank()) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + 
  ylab("Small coding CNVs")

if(length(dt.out$rSNVs_score) > 0)
  p2 <- ggplot(dt.out, aes(x = adm.dysmorphic, y = rSNVs_score, fill = adm.dysmorphic)) + geom_violin() +
  annotate("text", label = sprintf("p=%s", p2.p), x = c(1.5), y = c(1.12), cex = 5) +
  annotate("line", x = c(1, 2), y = c(1.05, 1.05), group = c("A", "A")) +
  geom_boxplot(fill = "white", width = .2, alpha=.6) + theme_classic() +     
  theme(legend.position = "none", axis.title.x = element_blank()) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + 
  ylab("Coding rare SNVs")
if(length(dt.out$dSNVs_score) > 0)
  p3 <- ggplot(dt.out, aes(x = adm.dysmorphic, y = dSNVs_score, fill = adm.dysmorphic)) + geom_violin() +
  annotate("text", label = sprintf("p=%s", p3.p), x = c(1.5), y = c(1.12), cex = 5) +
  annotate("line", x = c(1, 2), y = c(1.05, 1.05), group = c("A", "A")) +
  geom_boxplot(fill = "white", width = .2, alpha=.6) + theme_classic() +     
  theme(legend.position = "none", axis.title.x = element_blank()) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + 
  ylab("Coding de novo SNVs")
p4 <- ggplot(dt.out, aes(x = adm.dysmorphic, y = totalScore, fill = adm.dysmorphic)) + geom_violin() +
  annotate("text", label = sprintf("p=%s", p4.p), x = c(1.5), y = c(1.12), cex = 5) +
  annotate("line", x = c(1, 2), y = c(1.05, 1.05), group = c("A", "A")) +
  geom_boxplot(fill = "white", width = .2, alpha=.6) + theme_classic() +     
  theme(legend.position = "none", axis.title.x = element_blank()) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + 
  ylab("Total risk score")

##nc
if(length(dt.out$noncodeCNVs_score) > 0)
p5 <- ggplot(dt.out, aes(x = adm.dysmorphic, y = noncodeCNVs_score, fill = adm.dysmorphic)) + geom_violin() +
  annotate("text", label = sprintf("p=%s", p5.p), x = c(1.5), y = c(1.12), cex = 5) +
  annotate("line", x = c(1, 2), y = c(1.05, 1.05), group = c("A", "A")) +
  geom_boxplot(fill = "white", width = .2, alpha=.6) + theme_classic() +     
  theme(legend.position = "none", axis.title.x = element_blank()) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + 
  ylab("Noncoding CNVs")
if(length(dt.out$smallnoncodeCNVs_score) > 0)
  p6 <- ggplot(dt.out, aes(x = adm.dysmorphic, y = smallnoncodeCNVs_score, fill = adm.dysmorphic)) + geom_violin() +
  annotate("text", label = sprintf("p=%s", p6.p), x = c(1.5), y = c(1.12), cex = 5) +
  annotate("line", x = c(1, 2), y = c(1.05, 1.05), group = c("A", "A")) +
  geom_boxplot(fill = "white", width = .2, alpha=.6) +theme_classic() +     
  theme(legend.position = "none", axis.title.x = element_blank()) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + 
  ylab("Small noncoding CNVs")
if(length(dt.out$noncoderSNVs_score) > 0)
  p7 <- ggplot(dt.out, aes(x = adm.dysmorphic, y = noncoderSNVs_score, fill = adm.dysmorphic)) + geom_violin() +
  annotate("text", label = sprintf("p=%s", p7.p), x = c(1.5), y = c(1.12), cex = 5) +
  annotate("line", x = c(1, 2), y = c(1.05, 1.05), group = c("A", "A")) +
  geom_boxplot(fill = "white", width = .2, alpha=.6) + theme_classic() +     
  theme(legend.position = "none", axis.title.x = element_blank()) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + 
  ylab("Noncoding rare SNVs")
if(length(dt.out$noncodedSNVs_score) > 0)
  p8 <- ggplot(dt.out, aes(x = adm.dysmorphic, y = noncodedSNVs_score, fill = adm.dysmorphic)) + geom_violin() +
  annotate("text", label = sprintf("p=%s", p8.p), x = c(1.5), y = c(1.12), cex = 5) +
  annotate("line", x = c(1, 2), y = c(1.05, 1.05), group = c("A", "A")) +
  geom_boxplot(fill = "white", width = .2, alpha=.6) + theme_classic() +     
  theme(legend.position = "none", axis.title.x = element_blank()) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + 
  ylab("Noncoding de novo SNVs")

# if(length(dt.out$rSNVs_score) > 0){
#   cowplot::plot_grid(p1, p0, p2, p3, p4, nrow = 2)
# }else if(length(dt.out$dSNVs_score) > 0){
#   cowplot::plot_grid(p1, p0, p4, nrow = 2)
# }else{
#   cowplot::plot_grid(p1, p0, p3, p4, nrow = 2)
# }
cowplot::plot_grid(p1, p0, p5, p6, p2, p3, p7, p8, p4, nrow = 3)
# 
# message(mean(dt.out$CNVs_score[dt.out$adm.dysmorphic == "dys"]))
# message(mean(dt.out$CNVs_score[dt.out$adm.dysmorphic == "nondys"]))
# message(mean(dt.out$totalScore[dt.out$adm.dysmorphic == "dys"]))
# message(mean(dt.out$totalScore[dt.out$adm.dysmorphic == "nondys"]))
message("+++++++++++++++++++++++")
  
ggsave(sprintf("ssc.%s.score.pdf",p), width = 12, height = 8)
#}
dt.out$adm.dysmorphic <- as.character(dt.out$adm.dysmorphic)
dt.out$adm.dysmorphic <- gsub("unaff", "Unaffected sibling", dt.out$adm.dysmorphic)
dt.out$adm.dysmorphic <- gsub("dys", "ADM\nDysmorphic", dt.out$adm.dysmorphic)
dt.out$adm.dysmorphic <- gsub("nonADM\nDysmorphic", "ADM\nNondysmorphic", dt.out$adm.dysmorphic)
dt.out$adm.dysmorphic <- factor(dt.out$adm.dysmorphic, levels = c("Unaffected sibling", "ADM\nNondysmorphic", "ADM\nDysmorphic"))

ggplot(dt.out, aes(x = adm.dysmorphic, y = totalScore, fill = adm.dysmorphic)) + geom_violin(color = NA) +
  annotate("text", label = sprintf("p=%s", p4.p), x = c(2.5), y = c(1.12), cex = 5) +
  annotate("text", label = sprintf("p=%s", p4.unaff.p), x = c(2.5), y = c(1.27), cex = 5) +
  annotate("line", x = c(2, 3), y = c(1.05, 1.05), group = c("A", "A")) +
  annotate("line", x = c(1, 3), y = c(1.2, 1.2), group = c("A", "A")) +
  geom_boxplot(fill = "white", width = .2, alpha=.6) + theme_classic() +     
  theme(legend.position = "none", axis.title.x = element_blank()) + ylab("GRS percentile") + 
  scale_fill_manual(values = c("grey", "#009eed", "#F28602"))
ggsave("ssc.grs.pdf", width = 5, height = 4)
  #scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + 
  
# 
# sample CNVs_score dSNVs_score rSNVs_score totalScore adm.dysmorphic adm.autism_type
# 148 SSC06764   50.04089           0    8.439151   58.48004         nondys         complex
# 343 SSC10210   71.17565           0   10.139232   81.31488         nondys       essential
# 383 SSC10984   53.10212           0    8.463880   61.56600            dys         complex
