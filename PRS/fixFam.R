fam <- read.delim("bed/NFLD.fam", stringsAsFactors = F, header = F, sep = " ")
fam$V6 <- rep(1:2, nrow(fam))[1:nrow(fam)]
write.table(fam, "bed/NFLD.fam", sep=" ", row.names=F, quote=F, col.names=F)

