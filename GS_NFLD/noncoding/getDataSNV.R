setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(data.table)
library(GenomicRanges)

rm(list = ls())
transcript <- read.delim("~/Documents/doc/working/ReferenceData/refflat/refFlat.hg19.txt",
                         stringsAsFactors = F, header = F)
names(transcript) <- c("gsymbol", "isoform", "chr", "strand", "start", "end", "cds.start", "cds.end", "exons",
                       "exon.starts", "exon.ends")
transcript$tss.2kb.start <- ifelse(transcript$strand == "+",
                              transcript$start - 2000, 
                              transcript$end)
transcript$tss.2kb.end <- ifelse(transcript$strand == "+",
                                 transcript$start,
                                 transcript$end + 2000)
  
transcript$tss.10kb.start <- ifelse(transcript$strand == "+",
                                   transcript$start - 10000, 
                                   transcript$end)
transcript$tss.10kb.end <- ifelse(transcript$strand == "+",
                                 transcript$start,
                                 transcript$end + 10000)

snv.rare <- fread("../data/SNVs/calls/rare/all.rare.variants.OCT2019.tsv")
snv.rare <- data.frame(snv.rare, stringsAsFactors = F)
snv.rare <- snv.rare[which(snv.rare$chr %in% paste0("chr", 1:22)), ]

snv.rare.g <- GRanges(snv.rare$chr, IRanges(snv.rare$start, snv.rare$end), "*")

### conserved promoter
cpromo <- read.delim("../dataNC/an.conserved.promoters.hg19.bed", stringsAsFactors = F, header = F)
cpromo.g <- GRanges(cpromo$V1, IRanges(cpromo$V2, cpromo$V3), "*")

olap <- findOverlaps(snv.rare.g, cpromo.g)
snv.rare$An_promoter <- 0
snv.rare$An_promoter[unique(olap@from)] <- 1
snv.rare$An_promoter_tier2 <- snv.rare$An_promoter
snv.rare$An_promoter_tier2[snv.rare$phastCons_placental < 0 |
                             snv.rare$phylopPMam_avg < 1.5] <- 0

### DDD genes
ddd <- readLines("../dataNC/DDD.genes.27092019.tsv")
ddd <- transcript[transcript$gsymbol %in% ddd, ]
ddd.g <- GRanges(ddd$chr, IRanges(ddd$tss.2kb.start, ddd$tss.2kb.end), "*")
olap <- findOverlaps(snv.rare.g, ddd.g)
snv.rare$DDD_promoter <- 0
snv.rare$DDD_promoter[unique(olap@from)] <- 1
snv.rare$DDD_promoter_tier2 <- snv.rare$DDD_promoter
snv.rare$DDD_promoter_tier2[snv.rare$phastCons_placental < 0 |
                              snv.rare$phylopPMam_avg < 1.5] <- 0

### lncRNA
lnc <- read.delim("../dataNC/GENCODE.lncRNA.tsv", stringsAsFactors = F)
lnc$tss.2kb.start <- ifelse(lnc$V7 == "+",
                            lnc$V4 - 2000, 
                            lnc$V5)
lnc$tss.2kb.end <- ifelse(lnc$V7 == "+",
                          lnc$V4,
                          lnc$V5 + 2000)
lnc.g <- GRanges(lnc$V1, IRanges(lnc$tss.2kb.start, lnc$tss.2kb.end), "*")
olap <- findOverlaps(snv.rare.g, lnc.g)
snv.rare$lnc_promoter <- 0
snv.rare$lnc_promoter[unique(olap@from)] <- 1
snv.rare$lnc_promoter_tier2 <- snv.rare$lnc_promoter
snv.rare$lnc_promoter_tier2[snv.rare$phastCons_placental < 0 |
                              snv.rare$phylopPMam_avg < 1.5] <- 0

### gnomad
intolerance <- read.delim("../../../ReferenceData/constraint/gnomAD_oe_lof0.35_v2.1_Feb14_b38_coordinate.txt",
                          stringsAsFactors = F)
snv.rare$tss_intolerance <- as.numeric(snv.rare$typeseq_priority == "upstream" & snv.rare$gene_symbol %in% intolerance$gene)
snv.rare$tss_intolerance_tier2 <- snv.rare$tss_intolerance
snv.rare$tss_intolerance_tier2[snv.rare$phastCons_placental < 0 |
                                 snv.rare$phylopPMam_avg < 1.5] <- 0

### fetal brain
fetal.brain <- read.delim("../dataNC/fetal.brain.epiRoadmap.27092019.tsv", stringsAsFactors = F)
fetal.brain.g <- GRanges(fetal.brain$seqnames, IRanges(fetal.brain$start, fetal.brain$end), "*")
olap <- findOverlaps(snv.rare.g, fetal.brain.g)

snv.rare$fetal_brain_intolerance <- 0
snv.rare$fetal_brain_intolerance[intersect(unique(olap@from),
                                       which(snv.rare$gene_symbol %in% intolerance$gene))] <- 1
snv.rare$fetal_brain_intolerance_tier2 <- snv.rare$fetal_brain_intolerance
snv.rare$fetal_brain_intolerance_tier2[snv.rare$phastCons_placental < 0 |
                                         snv.rare$phylopPMam_avg < 1.5] <- 0

### splicing
snv.rare$splicing_tier1 <- as.numeric(snv.rare$spx_dpsi < -3.5)
snv.rare$splicing_tier2 <- as.numeric(snv.rare$spx_dpsi < -5)

### 3UTR
snv.rare$utr3_tier1 <- as.numeric(snv.rare$typeseq_priority == "UTR3" &
                                    snv.rare$phastCons_placental > 0)

snv.rare$utr3_tier2 <- as.numeric(snv.rare$typeseq_priority == "UTR3" &
                                    snv.rare$phastCons_placental > 0 &
                                    snv.rare$phylopPMam_avg >= 1.5)

### 5UTR
snv.rare$utr5 <- as.numeric(snv.rare$typeseq_priority == "UTR5")

### DeepSEA
deepsea1 <- read.delim("../dataNC/DeepSEA.tier1.tsv", stringsAsFactors = F)[, 1:6]
deepsea2 <- read.delim("../dataNC/DeepSEA.tier2.tsv", stringsAsFactors = F)[, 1:6]

snv.rare$deepsea_tier1 <- as.numeric(paste(snv.rare$chr, snv.rare$start, snv.rare$ref_allele, snv.rare$alt_allele, sep="@") %in%
                                       paste(deepsea1$chr, deepsea1$pos, deepsea1$ref, deepsea1$alt, sep="@"))
snv.rare$deepsea_tier2 <- as.numeric(paste(snv.rare$chr, snv.rare$start, snv.rare$ref_allele, snv.rare$alt_allele, sep="@") %in%
                                       paste(deepsea2$chr, deepsea2$pos, deepsea2$ref, deepsea2$alt, sep="@"))

### histone + gerp
histone <- read.delim("../dataNC/H3K27ac.marks.bed", stringsAsFactors = F, header = F)
transcript.g <- GRanges(transcript$chr, IRanges(transcript$tss.10kb.start, transcript$tss.10kb.end), "*")
histone.g <- GRanges(paste0("chr", histone$V1), IRanges(histone$V2, histone$V3), "*")

olap.trans <- findOverlaps(snv.rare.g, transcript.g)
olap.histone <- findOverlaps(snv.rare.g, histone.g)
olap <- intersect(olap.trans@from, olap.histone@from)

snv.rare$brain_h3k27ac <- 0
snv.rare$brain_h3k27ac[olap] <- 1
snv.rare$brain_h3k27ac_2gerp_tier2 <- snv.rare$brain_h3k27ac_2gerp
snv.rare$brain_h3k27ac_2gerp_tier2[snv.rare$phastCons_placental < 0 |
                                     snv.rare$phylopPMam_avg < 1.5] <- 0
# write.table(snv.rare[olap, ], "../dataNC/variants.overlap.H3K27ac.marks.bed", sep="\t",
#             row.names=F, quote=F, col.names=T)
gerp.annovar <- read.delim("../dataNC/rare_with_gerp.csv", sep=",", stringsAsFactors = F)
gerp.annovar$GERP.._RS <- as.numeric(gerp.annovar$GERP.._RS)
gerp.annovar <- unique(gerp.annovar[which(gerp.annovar$GERP.._RS > 2), 1:5])

snv.rare$brain_h3k27ac_2gerp <- 0
snv.rare$brain_h3k27ac_2gerp[paste(snv.rare$chr, snv.rare$start, snv.rare$end, snv.rare$ref_allele, snv.rare$alt_allele, sep=";") %in%
                               paste(gerp.annovar$Chr, gerp.annovar$Start, gerp.annovar$End, gerp.annovar$Ref, gerp.annovar$Alt, sep=";")] <- 1

### TAD
tad.boundary <- read.delim("../dataNC/tad.boundary.score.txt", stringsAsFactors = F)
tad.g <- GRanges(tad.boundary$chr, IRanges(tad.boundary$start, tad.boundary$end), "*")
olap <- findOverlaps(snv.rare.g, tad.g)

snv.rare$tad_tier1 <- 0
snv.rare$tad_tier1[intersect(unique(olap@from), which(snv.rare$phastCons_placental > 0))] <- 1
snv.rare$tad_tier2 <- 0
snv.rare$tad_tier2[intersect(unique(olap@from), which(snv.rare$phastCons_placental > 0 & snv.rare$phylopPMam_avg >= 1.5))] <- 1

### CTCF
ctcf <- read.delim("../dataNC/ctcf.hg19.bed", stringsAsFactors = F, header = F)
ctcf.g <- GRanges(ctcf$V1, IRanges(ctcf$V2, ctcf$V3), "*")
olap <- findOverlaps(snv.rare.g, ctcf.g)

snv.rare$ctcf_tier1 <- 0
snv.rare$ctcf_tier1[intersect(unique(olap@from), which(snv.rare$phastCons_placental > 0))] <- 1
snv.rare$ctcf_tier2 <- 0
snv.rare$ctcf_tier2[intersect(unique(olap@from), which(snv.rare$phastCons_placental > 0 & snv.rare$phylopPMam_avg >= 1.5))] <- 1

### brain enhancer
brain.enh <- read.delim("../dataNC/brainEnh_epiRoadmap_20190927.tsv", stringsAsFactors = F)
brain.g <- GRanges(brain.enh$seqnames, IRanges(brain.enh$start, brain.enh$end), "*")
olap <- findOverlaps(snv.rare.g, brain.g)

snv.rare$brain_enh_tier1 <- 0
snv.rare$brain_enh_tier1[intersect(unique(olap@from), which(snv.rare$phastCons_placental > 0))] <- 1
snv.rare$brain_enh_tier2 <- 0
snv.rare$brain_enh_tier2[intersect(unique(olap@from), which(snv.rare$phastCons_placental > 0 & snv.rare$phylopPMam_avg >= 1.5))] <- 1
snv.rare <- snv.rare[rowSums(snv.rare[, 21:40], na.rm = T) > 0, ]

write.table(snv.rare, "../dataNC/snv.rare.annotated.tsv", sep="\t", row.names=F, quote=F, col.names=T)
