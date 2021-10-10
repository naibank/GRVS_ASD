setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
library(ggplot2)
scores <- read.delim("ssc.grs.p0.1.tsv", stringsAsFactors = F)
scores$totalScore <- scale(scores$totalScore)

fam <- data.frame(readxl::read_excel("../SSC_data/n422.SSCaff_unaffsib_samplelist.xlsx", sheet = 1), stringsAsFactors = F)

unaff.dis <- fam$Family_ID[which(fam$adm.dysmorphic == "dys")]
unaff.nondis <- fam$Family_ID[which(fam$adm.dysmorphic == "nondys")]

unaff.dis <- fam$Sample.ID[fam$Family_ID %in% unaff.dis & fam$adm.dysmorphic == "NA"]
unaff.nondis <- fam$Sample.ID[fam$Family_ID %in% unaff.nondis & fam$adm.dysmorphic == "NA"]

formatC(wilcox.test(scores$totalScore[scores$sample %in% unaff.dis], 
                                  scores$totalScore[scores$sample %in% unaff.nondis], alternative = "greater")$p.value, format = "e", digits = 1)
# p=2.7e-01
boxplot(scores$totalScore[scores$sample %in% unaff.dis],
        scores$totalScore[scores$sample %in% unaff.nondis], names = c("Unaff (Dys fam)", "Unaff (Nondys fam)"))
