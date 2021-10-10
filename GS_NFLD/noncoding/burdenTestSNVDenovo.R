setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

rm(list = ls())
library(data.table)
source("burdenRequireFunctions.R")

total.count <- read.delim("../data/SNVs/denovo.gs.matrix.tsv", stringsAsFactors = F)
total.count <- total.count[, c("Sample", "TotalRare")]

feat <- c("An_promoter", "DDD_promoter", "lnc_promoter", "tss_intolerance", "fetal_brain_intolerance", 
          "An_promoter_tier2", "DDD_promoter_tier2", "lnc_promoter_tier2", "tss_intolerance_tier2", 
          "fetal_brain_intolerance_tier2", "splicing_tier1", 
"splicing_tier2", "utr3_tier1", "utr3_tier2", "utr5", "deepsea_tier1", "deepsea_tier2", "brain_h3k27ac", "brain_h3k27ac_2gerp", 
"tad_tier1", "tad_tier2", "ctcf_tier1", "ctcf_tier2", "brain_enh_tier1", "brain_enh_tier2")        
dt <- data.frame(fread("../dataNC/snv.denovo.annotated.tsv"), stringsAsFactors = F)
feat <- feat[feat %in% names(dt)]
dt <- aggregate(.~Sample, dt[, c("Sample", feat)], sum, na.rm = T, na.action = NULL)

dt <- merge(dt, total.count, by = "Sample", all = T)
dt[is.na(dt)] <- 0

### get only columns that I use
info <- read.delim("../data/2021.02.01-updated NFLD phenotype table IQ ADM and 3 status.txt", stringsAsFactors = F)
info <- info[, c("WGS_ManuID", "DNA.source", "Relation", "Sex..CRV", "incl.aff.in.my.study", "WGS.status", "incld.cohort",  "Paternal.age",
                 "remove.family.for.rare.burden", "Platform", "Dysmorphology.classification", "Grouped_Dysmorphology", "Final.Dysmorphology.scores")]
info <- info[info$incl.aff.in.my.study == 1 & info$WGS.status == "WGS-ed, complete fam", ]

### read estimated PCA
pca <- read.delim("../estimated.pca.12March2019.tsv", stringsAsFactors = F)
pca <- pca[!pca$sample %in% c("3-0456-000", "3-0458-000"), ]
pca$sample <- gsub("A", "", pca$sample)

### merge PCA to phenotype table
info <- merge(info, pca[, c("sample", "pc1", "pc2", "pc3")], by.x = "WGS_ManuID", by.y = "sample", all.x =T)
info$Sample <- gsub("-", "_", info$WGS_ManuID)

sample.not.in <- info$Sample[!info$Sample %in% dt$Sample]
sample.not.in <- data.frame("Sample" = sample.not.in)
sample.not.in[names(dt)[-1]] <- 0

dt <- rbind(dt, sample.not.in)
merge <- merge(dt, info, by = "Sample", all.x = T)

library(GSBurden)

### load geneset
snv.matrix <- merge
### encode label variables used to build different regression models
snv.matrix$Grouped_Dysmorphology[snv.matrix$Grouped_Dysmorphology == "Equivocal"] <- "Non-Essential"
snv.matrix$MultiClass <- ifelse(snv.matrix$Dysmorphology.classification == "essential", 0, 1)
snv.matrix$MultiClass <- ifelse(snv.matrix$Dysmorphology.classification == "complex", 2, snv.matrix$MultiClass)

snv.matrix$ComplexVsEssential <- ifelse(snv.matrix$MultiClass == 1, NA, 0)
snv.matrix$ComplexVsEssential <- ifelse(snv.matrix$MultiClass == 2, 1, snv.matrix$ComplexVsEssential)

snv.matrix$MultiClass <- factor(snv.matrix$MultiClass, levels = 0:2)

snv.matrix$GroupedDysmorphology <- ifelse(snv.matrix$Grouped_Dysmorphology == "Non-Essential", 1, 0)

snv.matrix$DysmorphologyScore <- factor(snv.matrix$Final.Dysmorphology.scores, levels = sort(unique(as.numeric(snv.matrix$Final.Dysmorphology.scores))))
write.table(snv.matrix, "nc.snv.denovo.matrix.tsv", sep="\t", row.names = F, quote=F, col.names=T)

### define covariates
snv.matrix$Paternal.age <- as.numeric(snv.matrix$Paternal.age)

covariates <- c("Sex..CRV", "Platform", "pc1", "pc2", "pc3", "TotalRare")
perm.values <- data.frame()
for(test in c("GroupedDysmorphology", "ComplexVsEssential", "MultiClass", "DysmorphologyScore")){
  test.out <- CNVBurdenTest(snv.matrix[!is.na(snv.matrix[, test]),], feat, test, covariates, nperm = 1000)
  
  burden.out <- test.out$Test
  perm.temp <- test.out$Permutation.Test
  perm.values <- rbind(perm.values, perm.temp)
  write.table(burden.out, sprintf("../resultNC/denovo.%s.burden.tsv", test), sep="\t", row.names=F, quote=F, col.names=T)
}

for(test in c("GroupedDysmorphology", "ComplexVsEssential", "MultiClass", "DysmorphologyScore")){
  dt <- read.delim(sprintf("../resultNC/denovo.%s.burden.tsv", test), stringsAsFactors = F)
  dt <- dt[!is.na(dt$coefficient), ]
  dt$permFDR <- sapply(dt$pvalue, getPermPvalue, dt$pvalue, perm.values$pvalue)
  dt$permFDR <- sapply(1:nrow(dt), getPermFDR, dt)

  write.table(dt, sprintf("../resultNC/denovo.%s.burden.withPerm.tsv", test), sep="\t", row.names=F, quote=F, col.names=T)
}
