setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
library(data.table)
library(GSBurden)

load("../../main/requiredData/gsNFLD.RData")

appris <- read.delim("../../../ReferenceData/APPRIS/appris_data.principal.hg19.txt",
                     stringsAsFactors = F, header = F, col.names = c("gsymbol", "enzid", "isoid", "x1", "tier"))
appris$isoid <- sapply(sapply(appris$isoid, strsplit, "\\."), "[", 1)
appris <- appris[grep("PRINCIPAL", appris$tier), ]

trscrpt <- read.delim("../../../ReferenceData/RefGene_fromJohn/hg19_refGene_28_04_17.transcript.txt", 
                      stringsAsFactors = F, header = F,
                      col.names = c("chr", "start", "end", "isoid", "gsymbol", "enzid"))
appris <- trscrpt[which(trscrpt$isoid %in% appris$isoid), ]

dt <- read.delim("../CNVs/CNVs.SSC.freq_1percent.HQR.hg19.tsv", stringsAsFactors = F)
dt <- dt[which(dt$SIZE < 3000000 & dt$SIZE > 10000 & dt$Sample_QC == "ok" & dt$chrAnn %in% paste0("chr", 1:22)), ]
# dt <- dt[which(dt$sscCnvnPercFreq_50percRecipOverlap < 1 &
#                  dt$sscErdsPercFreq_50percRecipOverlap < 1), ]
ssc.counts <- data.frame(table(dt$sample, dt$SVTYPE), stringsAsFactors = F)
q3.del <- quantile(ssc.counts$Freq[ssc.counts$Var2 == "DEL"], 0.75)
iqr.del <- IQR(ssc.counts$Freq[ssc.counts$Var2 == "DEL"])
cutoff <- q3.del + 3*iqr.del

outliers <- unique(as.character(ssc.counts$Var1[ssc.counts$Freq > cutoff & ssc.counts$Var2 == "DEL"]))

eth <- read.delim("../SSC_data/Target.admixture.tagged.tsv", stringsAsFactors = F)
outliers <- outliers[outliers %in% eth$IND[eth$ETH_TAG == "EUR"]]
writeLines(outliers, "CNVOutlier.txt")

ssc_rare <- read.delim("../SNVs/SSC_rare.gs.matrix.tsv", stringsAsFactors = F)
dt <- dt[dt$sample %in% ssc_rare$sample & (!dt$sample %in% outliers), ]

cnvs <- CNVSet(dt$sample, dt$CHROM, dt$START, dt$END, dt$SVTYPE)
cnvs$size <- cnvs$end - cnvs$start + 1

gene.in <- na.omit(read.delim("../../main/requiredData/hg19_refGene_28_04_17.cds.txt", header = F, stringsAsFactors = F))
names(gene.in) <- c("chr", "start", "end", "isoid", "genesymbol", "entrezid")
genes <- GeneAnnotation(gene.in$entrezid, gene.in$chr, gene.in$start, gene.in$end, gene.in$genesymbol)

cnv.matrix <- getCNVGSMatrix(cnvs, genes, gsNFLD)

# 
# cnv.matrix <- getCNVGSMatrix(cnvs[cnvs$type == "deletion", ], genes, gsNFLD)
# cnv.dup.matrix <- getCNVDupGSMatrix(cnvs[cnvs$type == "duplication", ], genes, appris, gsNFLD)
# names(cnv.matrix)[-1] <- paste(names(cnv.matrix)[-1], "DEL", sep="_")
# 
# cnv.matrix <- merge(cnv.matrix, cnv.dup.matrix, by = "sample", all = T)
# cnv.matrix[is.na(cnv.matrix)] <- 0


write.table(cnv.matrix, "../CNVs/SSC.gs.matrix.tsv", sep="\t", row.names=F, quote=F, col.names=T)
nfld <- read.delim("../../main/scripts/ADM/cnv.coding.matrix.tsv", stringsAsFactors = F)

summary(cnv.matrix$cnv_count_DEL)
summary(nfld$cnv_count_DEL)

