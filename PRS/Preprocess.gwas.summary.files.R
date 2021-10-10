#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(biomaRt) # biomaRt_2.30.0
library(snplist)

ipsych <- read.delim("iPSYCH/iPSYCH-PGC_ASD_Nov2017_P0.1_INFO0.9_SNPs.tsv", stringsAsFactors = F)
bmi <- read.delim("BMI/2020.02.12-BMI_SNPs_p0.2.txt", stringsAsFactors = F)
nfld <- read.delim("bed/NFLD.bim", stringsAsFactors = F, header = F)
table(ipsych$A1, ipsych$A2)
table(bmi$A1, bmi$A2)

#remove ambiguous SNPs i.e. A-T and C-G SNPs
ipsych <- ipsych[!(ipsych$A1 == "A" & ipsych$A2 == "T") &
                   !(ipsych$A1 == "C" & ipsych$A2 == "G") &
                   !(ipsych$A1 == "T" & ipsych$A2 == "A") &
                   !(ipsych$A1 == "G" & ipsych$A2 == "C"), ]
bmi <- bmi[!(bmi$A1 == "A" & bmi$A2 == "T") &
             !(bmi$A1 == "T" & bmi$A2 == "A") &
             !(bmi$A1 == "C" & bmi$A2 == "G") &
             !(bmi$A1 == "G" & bmi$A2 == "C"), ]
nfld <- nfld[nchar(nfld$V5) == 1 & nchar(nfld$V6) == 1, ]

nfld <- unique(nfld)
nfld.id <- paste(nfld$V1, nfld$V4, nfld$V5, nfld$V6, sep="#")

### clean iPSYCH data
ipsych.id <- paste(ipsych$CHR, ipsych$BP, ipsych$A1, ipsych$A2, sep="#")
ipsych.id.swap <- paste(ipsych$CHR, ipsych$BP, ipsych$A2, ipsych$A1, sep="#")

ipsych.match <- which(ipsych.id %in% nfld.id)
ipsych.swap <- which(ipsych.id.swap %in% nfld.id)

ipsych.swap <- setdiff(ipsych.swap, ipsych.match)

ipsych.swap.a1 <- ipsych$A2[ipsych.swap]
ipsych.swap.a2 <- ipsych$A1[ipsych.swap]
ipsych.or <- 1/ipsych$OR[ipsych.swap]

ipsych[ipsych.swap, c("A1", "A2", "OR")] <- data.frame(ipsych.swap.a1, ipsych.swap.a2, ipsych.or, stringsAsFactors = F)
ipsych <- ipsych[c(ipsych.match, ipsych.swap), ]

write.table(ipsych, "iPSYCH/iPSYCH-PGC_ASD_Nov2017_P0.1_INFO0.9_SNPs_clean.tsv",
            sep="\t", row.names=F, quote=F, col.names=T)

### clean bmi data
snp.map <- read.delim("../iPSYCH-PGC_ASD_BMI.bed", stringsAsFactors = F, header = F)
snp.map <- snp.map[, c("V1", "V3", "V4")]
names(snp.map) <- c("CHR", "BP", "SNP")

bmi <- merge(bmi, snp.map, by = "SNP", all.x = T)
bmi$b <- exp(bmi$b)
names(bmi) <- c("SNP", "A1", "A2", "FREQ", "OR", "SE", "P", "N", "CHR", "BP")
bmi <- bmi[, c("CHR", "SNP", "BP", "A1", "A2", "FREQ", "OR", "SE", "P", "N")]

bmi.id <- paste(bmi$CHR, bmi$BP, bmi$A1, bmi$A2, sep="#")
bmi.id.swap <- paste(bmi$CHR, bmi$BP, bmi$A2, bmi$A1, sep="#")

bmi.match <- which(bmi.id %in% nfld.id)
bmi.swap <- which(bmi.id.swap %in% nfld.id)

bmi.swap <- setdiff(bmi.swap, bmi.match)

bmi.swap.a1 <- bmi$A2[bmi.swap]
bmi.swap.a2 <- bmi$A1[bmi.swap]
bmi.or <- 1/bmi$OR[bmi.swap]

bmi[bmi.swap, c("A1", "A2", "OR")] <- data.frame(bmi.swap.a1, bmi.swap.a2, bmi.or, stringsAsFactors = F)
bmi <- bmi[c(bmi.match, bmi.swap), ]

write.table(bmi, "../BMI/2020.02.12-BMI_SNPs_p0.2_clean.tsv",
            sep="\t", row.names=F, quote=F, col.names=T)

