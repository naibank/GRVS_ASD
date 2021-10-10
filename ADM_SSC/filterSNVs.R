setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
gc();gc();
library(R.utils)
library(GenomicRanges)
files <- list.files("../SNVs/SSC/")
sel.columns <- c("chr", "start", "end", "REF", "ALT", "vartype", 
                 "typeseq_priority", "effect_priority", "gene_symbol", "entrez_id", 
                 "sift_score", "polyphen_score", "ma_score", "phylopMam_avg",
                 "phylopVert100_avg", "phastCons_placental", "CADD_phred", "spx_dpsi", "mt_score")
gc();
rowMax <- function(this.dt){
  sapplyRowMax <- function(row, this.dt){
    return(max(as.numeric(as.character(this.dt[row, ])), na.rm = T))
  }
  
  return(unlist(sapply(1:nrow(this.dt), sapplyRowMax, this.dt)))
}

# dt.out <- data.frame()
for(f in files){
  message(sprintf("processing file no: %s of %s at %s", which(files == f), length(files), Sys.time()))
  dt <- data.table::fread(sprintf("../SNVs/SSC/%s", f), data.table = F)
  dt <- dt[is.na(dt$SegDup), ]
  dt <- dt[dt$`#HG19_CHROM` %in% paste0("chr", c(1:22, "X", "Y")), ]
  
  # dt <- dt[which(dt$gene_type == "protein-coding"), ]
  dt <- dt[which(dt$FILTER == "PASS"), ]
  # dt <- dt[which(dt$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain", "synonymous SNV", "nonsynonymous SNV") |
  #                (dt$typeseq_priority == "splicing")), ]

  
  dt$chr <- dt$`#HG19_CHROM`
  dt$start <- dt$HG19_POS
  dt$end <- dt$start + nchar(dt$ALT)
  dt$vartype <- ifelse(nchar(dt$ALT) == nchar(dt$REF), "snp", "indel")
  
  dt <- dt[which(dt$freq_max < 0.01 & dt$cgi_max_int_freq < 0.01 & dt$ilm_max_int_freq < 0.01 & dt$A1000g_freq_max < 0.01),]
  
  dt$zygosity <- sapply(sapply(dt[, 12], strsplit, ":"), "[", 1)
  ad <- sapply(sapply(sapply(dt[, 11], strsplit, ":"), "==", "AD"), which)
  # dt$alt_frac <- sapply(sapply(dt[, 12], strsplit, ":"), "[", ad)
  
  dt.2 <- dt[ad == 2, ]
  dt.3 <- dt[ad == 3, ]
  
  dt.2$alt_frac <- sapply(sapply(dt.2[, 12], strsplit, ":"), "[", 2)
  dt.3$alt_frac <- sapply(sapply(dt.3[, 12], strsplit, ":"), "[", 3)
  
  dt <- rbind(dt.2, dt.3)
  rm(list=c("dt.2", "dt.3"))
  gc(); gc();
  
  dt$alt_frac <- 1-as.numeric(sapply(sapply(dt$alt_frac, strsplit, ","), "[", 1))/
    colSums(sapply(sapply(dt$alt_frac, strsplit, ","), as.numeric), na.rm = T)
  
  gq <- sapply(sapply(sapply(dt[, 11], strsplit, ":"), "==", "GQ"), which)
  
  dt.4 <- dt[gq == 4, ]
  dt.5 <- dt[gq == 5, ]
  dt.6 <- dt[gq == 6, ]
  dt.4$GQ <- sapply(sapply(dt.4[, 12], strsplit, ":"), "[", 4)
  dt.5$GQ <- sapply(sapply(dt.5[, 12], strsplit, ":"), "[", 5)
  dt.6$GQ <- sapply(sapply(dt.6[, 12], strsplit, ":"), "[", 6)
  
  dt <- rbind(dt.4, dt.5, dt.6)
  rm(list=c("dt.4", "dt.5", "dt.6"))
  gc(); gc();
  dt$GQ <- as.numeric(dt$GQ)
  dt <- dt[which(
                   ((dt$vartype == "indel" & dt$GQ >= 90) |
             (dt$vartype == "snp" & dt$zygosity %in% c("0/1", "1/0") & 
                 dt$alt_frac >= 0.3 & dt$alt_frac <= 0.7 & dt$GQ >= 99) |
                (dt$vartype == "snp" & dt$zygosity %in% c("1/1") & dt$alt_frac > 0.7 &
                   dt$GQ >= 25))), ]
  dt <- dt[, sel.columns]
  dt$sample <- gsub(".tsv.gz", "", f)
  
  write.table(dt, gsub(".gz", "", sprintf("../SNVs/filtered/%s", f)), sep="\t", row.names=F, quote=F, col.names=T)
  # dt.out <- rbind(dt.out, dt)
  #[1] "chr"            "ref_allele"     "alt_allele"     "varType"        "phylopPMam_avg"
}

# write.table(dt.out, "../SNVs/SSC_rare1pct_SNVs_INDELs_QCed.tsv", sep="\t", row.names=F, quote=F, col.names=T)
# write.table(dt.out, "../SNVs/SSC_rare1pct_SNVs_INDELs_QCed.tsv", sep="\t", row.names=F, quote=F, col.names=T)

rm(list=ls())
library(data.table)
library(GenomicRanges)
setwd("../SNVs/")
system("head -1 filtered/SS0013000.tsv > all.variants.tsv")
system("tail -qn +2 filtered/*.tsv >> all.variants.tsv")
# write.table(dt, "all.variants.tsv", sep="\t", row.names=F, quote=F, col.names=T)

###########################################################
# transcript <- read.delim("../../../ReferenceData/refflat/refFlat.hg19.txt",
#                          stringsAsFactors = F, header = F)
# names(transcript) <- c("gsymbol", "isoform", "chr", "strand", "start", "end", "cds.start", "cds.end", "exons",
#                        "exon.starts", "exon.ends")
# 
# transcript$tss.10kb.start <- ifelse(transcript$strand == "+",
#                                     transcript$start - 10000, 
#                                     transcript$end)
# transcript$tss.10kb.end <- ifelse(transcript$strand == "+",
#                                   transcript$start,
#                                   transcript$end + 10000)
# 
# histone <- read.delim("../../main/dataNC/H3K27ac.marks.bed", stringsAsFactors = F, header = F)
# transcript.g <- GRanges(transcript$chr, IRanges(transcript$tss.10kb.start, transcript$tss.10kb.end), "*")
# histone.g <- GRanges(paste0("chr", histone$V1), IRanges(histone$V2, histone$V3), "*")
# 
# files <- list.files("filtered/", full.names = T)
# 
# dt <- data.frame()
# dt.genic <- data.frame()
# for(f in files){
#   message(which(f == files))
#   tmp <- fread(f, data.table = F)
#   
#   dt.genic <- rbind(dt.genic, tmp[which(tmp$effect_priority %in% c("frameshift deletion", "frameshift insertion", "stopgain", "synonymous SNV", "nonsynonymous SNV") |
#                                           (tmp$typeseq_priority == "splicing")), ])
#   tmp <- tmp[, 1:5]
#   names(tmp) <- c("Chr", "Start", "End", "Ref", "Alt")
#   tmp.g <- GRanges(tmp$Chr, IRanges(tmp$Start, tmp$End), "*")
#   olap.trans <- findOverlaps(tmp.g, transcript.g)
#   olap.histone <- findOverlaps(tmp.g, histone.g)
#   olap <- intersect(olap.trans@from, olap.histone@from)
#   
#   tmp <- tmp[olap, ]
#   
#   
#   dt <- unique(rbind(dt, tmp))
# }
# 
# 
# write.table(dt.genic,
#             "../SNVs/SSC_rare1pct_SNVs_INDELs_QCed.tsv", sep="\t", row.names=F, quote=F, col.names=T)
# write.table(dt, "forWannoVar.tsv", sep="\t", row.names=F, quote=F, col.names=T)
# 
# ##########
# ########## convert tsv to vcf
# dt <- data.table::fread("../SNVs/forWannoVar.tsv", data.table = F)
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
# names(dt) <- c("#CHROM", "POS", "End", "REF", "ALT")
# 
# dt$ID <- paste(dt$Chr, dt$Start, dt$End, dt$Ref, dt$Alt, sep="#")
# dt$QUAL <- 99
# dt$FILTER <- "PASS"
# dt$INFO <- "."
# dt$FORMAT <- "GT"
# dt$SAMPLE <- "1/0"
# dt <- dt[, c("#CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	"FILTER",	"INFO",	"FORMAT",	"SAMPLE")]
# write.table(dt, "../SNVs/forWannoVar.vcf", sep="\t", row.names=F, quote=F, col.names=T)
