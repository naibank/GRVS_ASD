setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

rm(list = ls())
library(data.table)
source("burdenRequireFunctions.R")

file.in <- read.delim("../data/all.cnvs.may6.txt", stringsAsFactors = F)
file.in <- file.in[file.in$size > 10000 & file.in$size < 3000000, ]
dt.out <- data.frame()
perm.values <- data.frame()

for(cnv in c("deletion", "duplication")){
  
  total.count <- data.frame(table(file.in$sample[file.in$type == cnv]))
  names(total.count) <- c("Sample", "TotalRare")
  total.count$Sample <- gsub("-", "_", total.count$Sample)
  
  feat <- c("utr_5", "utr_3", "An_promoter", "DDD_promoter",
            "lnc_promoter", "brain_h3k27ac", "tad", "ctcf", "brain_enh")        
  dt <- fread("../dataNC/cnv.annotated.tsv", data.table = F)
  dt <- dt[dt$chr %in% paste0("chr", 1:22), ]
  dt <- dt[dt$size > 10000 & dt$size < 3000000, ]
  
  dt <- dt[dt$ts_overlap == 0, ]
  dt <- dt[dt$type == cnv, ]
  
  dt$Sample <- gsub("-", "_", dt$sample)
  feat <- feat[feat %in% names(dt)]
  dt <- aggregate(.~Sample, dt[, c("Sample", feat)], sum, na.rm = T, na.action = NULL)
  dt <- merge(dt, total.count, by = "Sample", all = T)
  dt[is.na(dt)] <- 0
  dt.tmp <- dt
  names(dt.tmp)[-1] <- paste(names(dt.tmp)[-1], cnv, sep="_")
  
  if(cnv == "deletion"){
    dt.out <- dt.tmp
  }else{
    dt.out <- merge(dt.out, dt.tmp, by = "Sample", all = T)
    dt.out[is.na(dt.out)] <- 0
  }
  write.table(dt.out, "../dataNC/nc.cnv.matrix.tsv", sep="\t", row.names=F, quote=F, col.names=T)
    
  ### get only columns that I use
  info <- read.delim("../data/2021.02.01-updated NFLD phenotype table IQ ADM and 3 status.txt", stringsAsFactors = F)
  info <- info[, c("WGS_ManuID", "DNA.source", "Relation", "Sex..CRV", "incl.aff.in.my.study", "incld.cohort", "Paternal.age",
                   "remove.family.for.rare.burden", "Platform", "Dysmorphology.classification", "Grouped_Dysmorphology", "Final.Dysmorphology.scores")]
  info <- info[info$incl.aff.in.my.study == 1 & info$remove.family.for.rare.burden != "remove", ]
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

  ### define covariates
  covariates <- c("Sex..CRV", "Platform", "pc1", "pc2", "pc3", "TotalRare")
  for(test in c("GroupedDysmorphology", "ComplexVsEssential", "MultiClass", "DysmorphologyScore")){
    test.out <- CNVBurdenTest(snv.matrix[!is.na(snv.matrix[, test]),], feat, test, covariates, nperm = 1000)

    burden.out <- test.out$Test
    perm.temp <- test.out$Permutation.Test
    perm.values <- rbind(perm.values, perm.temp)
    write.table(burden.out, sprintf("../resultNC/cnv.%s.%s.burden.tsv", cnv, test), sep="\t", row.names=F, quote=F, col.names=T)
  }
}

write.table(perm.values, "../resultNC/cnv.perm.pvalues.tsv", sep="\t", row.names=F, quote=F, col.names=T)

for(cnv in c("deletion", "duplication")){
  for(test in c("GroupedDysmorphology", "ComplexVsEssential", "MultiClass", "DysmorphologyScore")){
    dt <- read.delim(sprintf("../resultNC/cnv.%s.%s.burden.tsv", cnv, test), stringsAsFactors = F)
    dt <- dt[dt$geneset != "ts_overlap", ]
    dt <- dt[!is.na(dt$coefficient), ]

    dt$permFDR <- sapply(dt$pvalue, getPermPvalue, dt$pvalue, perm.values$pvalue)
    dt$permFDR <- sapply(1:nrow(dt), getPermFDR, dt)
    
    write.table(dt, sprintf("../resultNC/cnv.%s.%s.burden.withPerm.tsv", cnv, test), sep="\t", row.names=F, quote=F, col.names=T)
  }

}

