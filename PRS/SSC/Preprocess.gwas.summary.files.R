# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# library(biomaRt) # biomaRt_2.30.0
# library(snplist)
library(data.table)
ipsych <- fread("/hpf/largeprojects/tcagstor/projects/MSSNG_SSC_PRS/SSC_PRS/iPSYCH_Final_For_SSC_4685642_SNPs_P_01.txt", data.table = F)
nfld <- read.delim("/hpf/largeprojects/tcagstor/projects/MSSNG_SSC_PRS/SSC_PRS/SSC_chr1_22_Filter_Pass_iPSYCH_SNPs_7261_Subjects.bim", stringsAsFactors = F, header = F)
names(ipsych) <- c("SNP", "A1", "A2", "CHR", "BP", "OR", "UNK", "P")
ipsych$OR <- exp(ipsych$OR)
#remove ambiguous SNPs i.e. A-T and C-G SNPs
message(nrow(ipsych))
ipsych <- ipsych[!(ipsych$A1 == "A" & ipsych$A2 == "T") &
                   !(ipsych$A1 == "C" & ipsych$A2 == "G") &
                   !(ipsych$A1 == "T" & ipsych$A2 == "A") &
                   !(ipsych$A1 == "G" & ipsych$A2 == "C"), ]
message(nrow(ipsych))

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
nfld$id <- nfld.id
ipsych$id <- paste(ipsych$CHR, ipsych$BP, ipsych$A1, ipsych$A2, sep="#")
ipsych <- merge(ipsych, nfld[, c("V2", "id")], by="id", all.x = T)
names(ipsych)[c(2, 10)] <- c("uniqID", "SNP")
# ipsych$SNP <- paste(ipsych$CHR, ipsych$BP, ipsych$A1, ipsych$A2, sep="#")
write.table(ipsych, "iPSYCH-PGC_ASD_Nov2017_P0.1_INFO0.9_SNPs_clean_hg38_for_SSC.tsv",
            sep="\t", row.names=F, quote=F, col.names=T)

####
#process fam
tmp <- readLines("plinkfiles/SSC_chr1_22_Filter_Pass_iPSYCH_SNPs_7261_Subjects.fam")
writeLines(tmp, "plinkfiles/old_SSC_chr1_22_Filter_Pass_iPSYCH_SNPs_7261_Subjects.fam")
tmp[1:3000] <- gsub("-9", "1", tmp[1:3000])
tmp[3001:7260] <- gsub("-9", "2", tmp[3001:7260])
writeLines(tmp, "plinkfiles/SSC_chr1_22_Filter_Pass_iPSYCH_SNPs_7261_Subjects.fam")
