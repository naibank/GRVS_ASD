setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(data.table)
library(GenomicRanges)

rm(list = ls())
transcript <- read.delim("../../../ReferenceData/refflat/refFlat.hg19.txt",
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

transcript$utr5.start <- transcript$start
transcript$utr5.end <- transcript$cds.start
transcript$utr3.start <- transcript$cds.end
transcript$utr3.end <- transcript$end

utr5 <- transcript[transcript$strand == "+", ]
utr3 <- transcript[transcript$strand == "-", ]

utr5.g <- union(GRanges(utr5$chr, IRanges(utr5$utr5.start, utr5$utr5.end), "*"),
               GRanges(utr3$chr, IRanges(utr3$utr3.start, utr3$utr3.end), "*"))
utr3.g <- union(GRanges(utr3$chr, IRanges(utr3$utr5.start, utr3$utr5.end), "*"),
                GRanges(utr5$chr, IRanges(utr5$utr3.start, utr5$utr3.end), "*"))

ts.g <- GRanges(transcript$chr, IRanges(transcript$cds.start, transcript$cds.end), "*")

snv.rare <- fread("../CNVs/CNVs.SSC.freq_1percent.HQR.hg19.tsv", data.table = F)
# snv.rare <- data.frame(snv.rare, stringsAsFactors = F)
snv.rare <- snv.rare[!is.na(snv.rare$chr), ]

snv.rare.g <- GRanges(snv.rare$chrAnn, IRanges(snv.rare$START, snv.rare$END), "*")

olap <- findOverlaps(snv.rare.g, utr5.g)
snv.rare$utr_5 <- 1:nrow(snv.rare) %in% olap@from
olap <- findOverlaps(snv.rare.g, utr3.g)
snv.rare$utr_3 <- 1:nrow(snv.rare) %in% olap@from

olap <- findOverlaps(snv.rare.g, ts.g)  
snv.rare$ts_overlap <- 1:nrow(snv.rare) %in% olap@from

### conserved promoter
cpromo <- read.delim("../../main/dataNC/an.conserved.promoters.hg19.bed", stringsAsFactors = F, header = F)
cpromo.g <- GRanges(cpromo$V1, IRanges(cpromo$V2, cpromo$V3), "*")

olap <- findOverlaps(snv.rare.g, cpromo.g)
snv.rare$An_promoter <- 0
snv.rare$An_promoter[unique(olap@from)] <- 1

### DDD genes
ddd <- readLines("../../main/dataNC/DDD.genes.2709main.tsv")
ddd <- transcript[transcript$gsymbol %in% ddd, ]
ddd.g <- GRanges(ddd$chr, IRanges(ddd$tss.2kb.start, ddd$tss.2kb.end), "*")
olap <- findOverlaps(snv.rare.g, ddd.g)
snv.rare$DDD_promoter <- 0
snv.rare$DDD_promoter[unique(olap@from)] <- 1

### lncRNA
lnc <- read.delim("../../main/dataNC/GENCODE.lncRNA.tsv", stringsAsFactors = F)
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

### histone + gerp
histone <- read.delim("../../main/dataNC/H3K27ac.marks.bed", stringsAsFactors = F, header = F)
transcript.g <- GRanges(transcript$chr, IRanges(transcript$tss.10kb.start, transcript$tss.10kb.end), "*")
histone.g <- GRanges(paste0("chr", histone$V1), IRanges(histone$V2, histone$V3), "*")

olap.trans <- findOverlaps(snv.rare.g, transcript.g)
olap.histone <- findOverlaps(snv.rare.g, histone.g)
olap <- intersect(olap.trans@from, olap.histone@from)

snv.rare$brain_h3k27ac <- 0
snv.rare$brain_h3k27ac[olap] <- 1
snv.rare$brain_h3k27ac[snv.rare$utr_5 | snv.rare$utr_3] <- 0
# write.table(snv.rare[olap, ], "../dataNC/variants.overlap.H3K27ac.marks.bed", sep="\t",
#             row.names=F, quote=F, col.names=T)
snv.rare$utr_5[snv.rare$ts_overlap] <- F
snv.rare$utr_3[snv.rare$ts_overlap] <- F
snv.rare$utr_5 <- as.numeric(snv.rare$utr_5)
snv.rare$utr_3 <- as.numeric(snv.rare$utr_3)

### TAD
tad.boundary <- read.delim("../../main/dataNC/tad.boundary.score.txt", stringsAsFactors = F)
tad.g <- GRanges(tad.boundary$chr, IRanges(tad.boundary$start, tad.boundary$end), "*")
olap <- findOverlaps(snv.rare.g, tad.g)

snv.rare$tad <- 0
snv.rare$tad[intersect(unique(olap@from), which(snv.rare$phastCons_placental > 0))] <- 1

### CTCF
ctcf <- read.delim("../../main/dataNC/ctcf.hg19.bed", stringsAsFactors = F, header = F)
ctcf.g <- GRanges(ctcf$V1, IRanges(ctcf$V2, ctcf$V3), "*")
olap <- findOverlaps(snv.rare.g, ctcf.g)

snv.rare$ctcf <- 0
snv.rare$ctcf[intersect(unique(olap@from), which(snv.rare$phastCons_placental > 0))] <- 1

### brain enhancer
brain.enh <- read.delim("../../main/dataNC/brainEnh_epiRoadmap_main0927.tsv", stringsAsFactors = F)
brain.g <- GRanges(brain.enh$seqnames, IRanges(brain.enh$start, brain.enh$end), "*")
olap <- findOverlaps(snv.rare.g, brain.g)

snv.rare$brain_enh <- 0
snv.rare$brain_enh[intersect(unique(olap@from), which(snv.rare$phastCons_placental > 0))] <- 1

write.table(snv.rare, "../dataNC/cnv.annotated.tsv", sep="\t", row.names=F, quote=F, col.names=T)
