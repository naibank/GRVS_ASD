setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
library(ordinal)
library(ggplot2)
library(ggforce)
library(cowplot)
library(GenomicRanges)
library(fmsb)

### get score from combination of variant types

cor.jaccard <- function(p1, p2){
  return(sum(p1==p2 & p1 == TRUE)/(sum(p1==p2 & p1 == TRUE) + sum(p1!=p2 & (p1 == TRUE | p2 == TRUE))))
}

getCNVEvents <- function(all.cnv, set.sig, gsNFLD){
  set.tmp <- set.sig
  
  dt.out <- data.frame()
  for(j in 1:2){
    if(j == 1){
      set.sig <- set.tmp[set.tmp$coeff > 0, ]
    }else{
      set.sig <- set.tmp[set.tmp$coeff < 0, ]
    }
    
    del <- grep("_deletion", set.sig$set)
    dup <- grep("_duplication", set.sig$set)
    
    set <- gsub("_deletion|_duplication", "", set.sig$set)
    if(length(set) > 0){
      
        del <- all.cnv[all.cnv$type == "deletion", c("sample", set[del])]
        dup <- all.cnv[all.cnv$type == "duplication", c("sample", set[dup])]
        
        sample <- data.frame()
        if(!is.vector(del)){
          if(ncol(del) == 2){
            del <- del[del[, 2] > 0, ]
            del$allcount <- del[, 2]
            del$varcount <- as.numeric(del[, 2] > 0)
          }else{
            del <- del[rowSums(del[, -1]) > 0, ]
            del$allcount <- rowSums(del[, -1])
            del$varcount <- rowSums(del[, -1]) > 0
          }
          sample <- rbind(sample, del[, c("sample", "allcount", "varcount")])
        }
        if(!is.vector(dup)){
          if(ncol(dup) == 2){
            dup <- dup[dup[, 2] > 0, ]
            dup$allcount <- dup[, 2]
            dup$varcount <- as.numeric(dup[, 2] > 0)
          }else{
            dup <- dup[rowSums(dup[, -1]) > 0, ]
            dup$allcount <- rowSums(dup[, -1])
            dup$varcount <- rowSums(dup[, -1]) > 0
            
          }
          sample <- rbind(sample, dup[, c("sample", "allcount", "varcount")])
        }
        
        # out <- data.frame(table(sample))
        names(sample) <- c("sample", sprintf("%s_events", ifelse(j == 1, "increasing_severity", "decreasing_severity")),
                           sprintf("%s_variants", ifelse(j == 1, "increasing_severity", "decreasing_severity")))
        
        if(nrow(dt.out) == 0){
          dt.out <- sample
        }else{
          dt.out <- merge(dt.out, sample, by = "sample", all = T)
        }
    }
  }
  
  dt.out[is.na(dt.out)] <- 0
  return(dt.out)
}

getSNVEvents <- function(all.rare, set.sig){
  set.tmp <- set.sig
  
  dt.out <- data.frame()
  for(j in 1:2){
    if(j == 1){
      set.sig <- set.tmp[set.tmp$coeff > 0, ]
    }else{
      set.sig <- set.tmp[set.tmp$coeff < 0, ]
    }
    lof <- grep("_lof", set.sig$set)
    ms1 <- grep("_tier1_ms", set.sig$set)
    ms2 <- grep("_tier2_ms", set.sig$set)
    
    set <- gsub("_lof|_tier1_ms|_tier2_ms", "", set.sig$set)
    
    if(length(set) > 0){
        lof <- all.rare[all.rare$effect.tier == "lof", c("Sample", set[lof])]
        ms1 <- all.rare[all.rare$effect.tier == "tier1_ms", c("Sample", set[ms1])]
        ms2 <- all.rare[all.rare$effect.tier == "tier2_ms", c("Sample", set[ms2])]
        
        sample <- c()
        if(!is.vector(lof)){
          if(ncol(lof) == 2){
            lof <- lof[lof[, 2], ]
          }else{
            lof <- lof[rowSums(lof[, -1]) > 0, ]
          }
          sample <- c(sample, lof$Sample)
        }
        if(!is.vector(ms1)){
          if(ncol(ms1) == 2){
            ms1 <- ms1[ms1[, 2], ]
          }else{
            ms1 <- ms1[rowSums(ms1[, -1]) > 0, ]
          }
          sample <- c(sample, ms1$Sample)
        }
        if(!is.vector(ms2)){
          if(ncol(ms2) == 2){
            ms2 <- ms2[ms2[, 2], ]
          }else{
            ms2 <- ms2[rowSums(ms2[, -1]) > 0, ]
          }
          sample <- c(sample, ms2$Sample)
        }
        
        out <- data.frame(table(sample))
        out[, "tmp2"] <- out[, 2]
        names(out) <- c("sample", sprintf("%s_events", ifelse(j == 1, "increasing_severity", "decreasing_severity")),
                        sprintf("%s_variants", ifelse(j == 1, "increasing_severity", "decreasing_severity")))
        
        if(nrow(dt.out) == 0){
          dt.out <- out
        }else{
          dt.out <- merge(dt.out, out, by = "sample", all = T)
        }
    }
  }
  
  dt.out[is.na(dt.out)] <- 0
  return(dt.out)
}

getSigFeatures <- function(dt, global, covariates, test, cutoff = 1, class = "ADM"){
  coeff <- c()
  p <- c()
  r2 <- c()
  correlation <- data.frame()
  for(feat in test){
    for(feat2 in test){
      cor <- cor.jaccard(dt[, feat] > 0, dt[, feat2] > 0)
      if(feat != feat2)
        correlation <- rbind(correlation, data.frame(feat, feat2, cor))
    }
    if(sd(dt[, feat]) != 0)
      dt[, feat] <- scale(dt[, feat])
    ref <- sprintf("factor(%s) ~ %s + %s", class, paste(covariates, collapse=" + "), global)
    add <- sprintf("%s + %s", ref, feat)
    
    ref <- glm(ref, dt, family = binomial(link = "logit"))
    add <- glm(add, dt, family = binomial(link = "logit"))
    pvalue <- anova(ref, add, test = "Chisq")$`Pr(>Chi)`[2]
    
    # print(sprintf("%s's pvalue = %s", feat, pvalue))
    if(!is.na(pvalue) & pvalue <= cutoff){
      coeff <- c(coeff, add$coefficients[feat])
      p <- c(p, pvalue)
      r2 <- c(r2, NagelkerkeR2(add)$R2)
      
    }
  }
  
  if(!is.null(coeff)){
    coeff <- data.frame("set" = names(coeff), "coeff" = coeff, "r2" = r2, "pvalue" = p, stringsAsFactors = F)
    coeff <- coeff[order(coeff$r2, decreasing = T), ]
    
    correlation <- correlation[which(correlation$cor > 0.75), ]
    if(nrow(correlation) > 0)
    for(i in 1:nrow(correlation)){
      feat1 <- correlation$feat[i]
      feat2 <- correlation$feat2[i]
      
      feat1.r2 <- coeff$r2[coeff$set == as.character(feat1)]
      feat2.r2 <- coeff$r2[coeff$set == as.character(feat2)]
      
      if(length(feat1.r2) > 0 & length(feat2.r2) > 0){
        if(feat1.r2 > feat2.r2){
          coeff <- coeff[coeff$set != as.character(feat2), ]
        }else{
          coeff <- coeff[coeff$set != as.character(feat1), ]
        }
      }
    }
    
    result <- coeff$set
    dt.out <- list(result, coeff)
    names(dt.out) <- c("set", "coeff")
    return(dt.out)
  }
  return(NULL)
}

getScoreZ <- function(dt, coeff, score.name){
  if(nrow(coeff) == 0){
    dt <- data.frame("Sample" = dt$Sample, "tmp" = NA, stringsAsFactors = F)
    names(dt)[2] <- score.name
  }else{
    z
    for(i in 1:nrow(coeff)){
      dt[, coeff$set[i]] <- dt[, coeff$set[i]] * coeff$coeff[i]
    }
    
    if(ncol(dt) == 2)
      dt[, score.name] <- dt[, 2]
    else
      dt[, score.name] <- rowSums(dt[, 2:ncol(dt)])
  }
  
  return(dt[, c("Sample", score.name)])
}

getScore <- function(dt, target, baseline, test, score.name, class = "ADM"){
  if(length(test) > 0){
#     ref <- gsub("\n", "", sprintf("%s ~ offset(I(%s * (Sex..CRV == 'M') + 
#             %s * (Platform == 'Macrogen-Illumina HiSeqX') + %s * (Platform == 'Complete Genomics') + 
# %s * (Platform == 'TCAG-Illumina HiSeqX') +
# %s * (pc1) + 
# %s * (pc2) + 
#                                   %s * (pc3)))", class,
#                                   baseline$coefficients["Sex..CRVM"], baseline$coefficients["PlatformMacrogen-Illumina HiSeqX"], baseline$coefficients["PlatformComplete Genomics"], 
#                                   baseline$coefficients["PlatformTCAG-Illumina HiSeqX"], baseline$coefficients["pc1"], baseline$coefficients["pc2"], baseline$coefficients["pc3"]))
    ref <- gsub("\n", "", sprintf("%s ~ offset(I(%s * (Sex..CRV == 'M') + 
            %s * (Platform == 'Macrogen-Illumina HiSeqX') + 
%s * (Platform == 'TCAG-Illumina HiSeqX') +
%s * (pc1) + 
%s * (pc2) + 
                                  %s * (pc3)))", class,
                                  baseline$coefficients["Sex..CRVM"], baseline$coefficients["PlatformMacrogen-Illumina HiSeqX"], 
                                  baseline$coefficients["PlatformTCAG-Illumina HiSeqX"], baseline$coefficients["pc1"], baseline$coefficients["pc2"], baseline$coefficients["pc3"]))
    
    add <- sprintf("%s + %s", ref, paste(test, collapse = " + "))
    lm <- glm(add, data = dt)
    sm <- summary(lm)
    sig <- rownames(sm$coefficients)[(sm$coefficients[, 'Pr(>|t|)'] < 1)]
    
    target <- target[, c("Sample", test[test %in% sig])]
    
    if(is.vector(target)){
      target <- data.frame("Sample" = target, "tmp" = NA, stringsAsFactors = F)
      names(target)[2] <- score.name
    }else if(ncol(target) < 3){
      target[, 2] <- target[, 2] * lm$coefficients[names(target)[2]]
      target[, score.name] <- target[, 2]
    }else{
      for(i in 2:ncol(target)){
        target[, i] <- target[, i] * lm$coefficients[names(target)[i]]
      }
      
      target[, score.name] <- rowSums(target[, 2:ncol(target)])
    }
  }else{
    target[, score.name] <- NA
  }
  
  return(target[, c("Sample", score.name)])
}

load("../burden_analysis/new_set/gsNFLD_2022.RData")

cds <- read.delim("../../main/requiredData/hg19_refGene_28_04_17.cds.txt", stringsAsFactors = F, header = F)
cds.g <- GRanges(cds$V1, IRanges(cds$V2, cds$V3), "*")

all.cnv <- rbind(read.delim("../burden_analysis/coding/CNVs/all.cnvs.June2022.txt", stringsAsFactors = F),
                 read.delim("../burden_analysis/coding/CNVs/all.small.cnvs.June2022.txt", stringsAsFactors = F))
all.cnv.g <- GRanges(all.cnv$chr, IRanges(all.cnv$start, all.cnv$end), "*")

olap <- data.frame(findOverlaps(all.cnv.g, cds.g))
olap$enzid <- cds$V6[olap$subjectHits]
for(i in 1:length(gsNFLD)){
  olap[, names(gsNFLD)[i]] <- F
  olap[olap$enzid %in% gsNFLD[[i]], names(gsNFLD)[i]] <- T
}
olap <- unique(olap[, -c(2:3)])
olap <- aggregate(.~queryHits, olap, sum)
# olap[, -1] <- olap[, -1] > 0

all.cnv[, names(olap)[-1]] <- 0
all.cnv[olap$queryHits, names(olap)[-1]] <- olap[, -1]

rare.lof <- read.delim("../../main/data/SNVs/calls/rare/lof.variants.Jun2022.tsv", stringsAsFactors = F)
rare.ms1 <- read.delim("../../main/data/SNVs/calls/rare/ms1.variants.Jun2022.tsv", stringsAsFactors = F)
rare.ms2 <- read.delim("../../main/data/SNVs/calls/rare/ms2.variants.Jun2022.tsv", stringsAsFactors = F)
rare.lof$effect.tier <- "lof"
rare.ms1$effect.tier <- "tier1_ms"
rare.ms2$effect.tier <- "tier2_ms"
all.rare <- rbind(rare.lof, rbind(rare.ms1, rare.ms2))

denovo.lof <- read.delim("../../main/data/SNVs/calls/denovo/lof.variants.Jun2022.tsv", stringsAsFactors = F)
denovo.ms1 <- read.delim("../../main/data/SNVs/calls/denovo/ms1.variants.Jun2022.tsv", stringsAsFactors = F)
denovo.ms2 <- read.delim("../../main/data/SNVs/calls/denovo/ms2.variants.Jun2022.tsv", stringsAsFactors = F)
denovo.lof$effect.tier <- "lof"
denovo.ms1$effect.tier <- "tier1_ms"
denovo.ms2$effect.tier <- "tier2_ms"
all.denovo <- rbind(denovo.lof, rbind(denovo.ms1, denovo.ms2))

### read all 5 files
cnv.coding <- read.delim("../burden_analysis/coding/CNVs/ADM.cnv.coding.matrix.tsv", stringsAsFactors = F)
cnv.coding$sample <- gsub("-", "_", cnv.coding$sample)
cnv.coding$sample <- gsub("A|_A", "", cnv.coding$sample)
names(cnv.coding)[1] <- "Sample"

scnv.coding <- read.delim("../burden_analysis/coding/CNVs/ADM.smaller.cnv.coding.matrix.tsv", stringsAsFactors = F)
scnv.coding$sample <- gsub("-", "_", scnv.coding$sample)
scnv.coding$sample <- gsub("A|_A", "", scnv.coding$sample)
names(scnv.coding)[1] <- "Sample"

### noncoding cnv
cnv.ncoding <- read.delim("../../main/dataNC/nc.cnv.matrix.tsv", stringsAsFactors = F)
scnv.ncoding <- read.delim("../../main/dataNC/nc.smaller.cnv.matrix.tsv", stringsAsFactors = F)

snv.coding.rare <- read.delim("../burden_analysis/coding/SNVs/rare.coding.snv.tsv", stringsAsFactors = F)

meta <- read.delim("../../main/data/2021.02.01-updated NFLD phenotype table IQ ADM and 3 status.txt", stringsAsFactors = F)
meta$WGS_ManuID <- gsub("-", "_", meta$WGS_ManuID)
meta$WGS_ManuID <- gsub("A|_A", "", meta$WGS_ManuID)
meta <- meta[meta$incl.aff.in.my.study == 1, ]
meta$ADM[grep("0456|0458", meta$WGS_ManuID)] <- "Nondysmorphic"
meta <- meta[which(meta$ADM %in% c("Dysmorphic", "Nondysmorphic")), ]
meta$ADM <- ifelse(meta$ADM == "Dysmorphic", 1, 0)

#fix PCs for four samples in cnv.coding
cnv.coding <- merge(cnv.coding, meta[, c("WGS_ManuID", "ADM")], by.x = "Sample", by.y = "WGS_ManuID", all = F)
scnv.coding <- merge(scnv.coding, meta[, c("WGS_ManuID", "ADM")], by.x = "Sample", by.y = "WGS_ManuID", all = F)

cnv.ncoding <- merge(cnv.ncoding, meta[, c("WGS_ManuID", "ADM")], by.x = "Sample", by.y = "WGS_ManuID", all = F)
scnv.ncoding <- merge(scnv.ncoding, meta[, c("WGS_ManuID", "ADM")], by.x = "Sample", by.y = "WGS_ManuID", all = F)

snv.coding.rare <- merge(snv.coding.rare, meta[, c("WGS_ManuID", "ADM")], by.x = "Sample", by.y = "WGS_ManuID", all = F)

cnv.coding <- merge(cnv.coding[, -c(90:94)], snv.coding.rare[, c(1, 139:146)])
scnv.coding <- merge(scnv.coding[, -c(90:94)], snv.coding.rare[, c(1, 139:146)])

cnv.ncoding <- merge(cnv.ncoding[, -c(22)], snv.coding.rare[, c(1, 130, 135, 139:146)], all = F)
scnv.ncoding <- merge(scnv.ncoding[, -c(22)], snv.coding.rare[, c(1, 130, 135, 139:146)], all = F)

snv.coding.denovo <- read.delim("../burden_analysis/coding/SNVs/denovo.coding.snv.tsv", stringsAsFactors = F)
snv.coding.denovo <- merge(snv.coding.denovo, meta[, c("WGS_ManuID", "ADM")], by.x = "Sample", by.y = "WGS_ManuID", all = F)

snv.nc.rare <- read.delim("../../main/scriptsNC/nc.snv.rare.matrix.tsv", stringsAsFactors = F)
snv.nc.rare <- merge(snv.nc.rare, meta[, c("WGS_ManuID", "ADM")], by.x = "Sample", by.y = "WGS_ManuID", all = F)

snv.nc.denovo <- read.delim("../../main/scriptsNC/nc.snv.denovo.matrix.tsv", stringsAsFactors = F)
snv.nc.denovo <- merge(snv.nc.denovo, meta[, c("WGS_ManuID", "ADM")], by.x = "Sample", by.y = "WGS_ManuID", all = F)


baseline.feat <- c("Sample", "Sex..CRV",  "Platform", "ADM", "pc1", "pc2", "pc3")
samples <- na.omit(unique(rbind(rbind(cnv.coding[, baseline.feat], scnv.coding[, baseline.feat]), rbind(snv.coding.rare[, baseline.feat],
                                                    rbind(snv.coding.denovo[, baseline.feat],
                                                          rbind(snv.nc.rare[, baseline.feat],
                                                                snv.nc.denovo[, baseline.feat]))))))
samples <- na.omit(unique(rbind(cnv.ncoding[, baseline.feat], rbind(samples, scnv.ncoding[, baseline.feat]))))
# samples <- samples[samples$Sample %in% snv.coding.denovo$Sample, ]
samples <- samples[samples$Platform != "", ]
samples <- samples[which(samples$Platform != "Complete Genomics"), ]
# samples <- na.omit(unique(rbind(cnv.coding[, baseline.feat], snv.coding.rare[, baseline.feat])))
set.seed(1500)
score.out <- data.frame()
sig.set <- c()
coeff.all <- data.frame()
covariates <- c("Sex..CRV", "Platform", "pc1", "pc2", "pc3")

for(k in 1:30){
nfold <- 10
fold.dys <- split(sample(samples$Sample[samples$ADM == 1]), 1:nfold)
fold.nondys <- split(sample(samples$Sample[samples$ADM == 0]), 1:nfold)

for(i in 1:length(fold.dys)){
  disc.set <- c(unlist(fold.dys[-i]), unlist(fold.nondys[-i]))
  targ.set <- c(unlist(fold.dys[i]), unlist(fold.nondys[i]))
  
  baseline <- glm("ADM ~ Sex..CRV + Platform + pc1 + pc2 + pc3", data = samples[samples$Sample %in% disc.set, ])
  
  # getSigFeatures <- function(dt, global, covariates, test, cutoff = 1, class = "ADM"){
  ## get significant geneset
  cnv.del <- names(cnv.coding)[grep("DEL", names(cnv.coding))]
  cnv.dup <- names(cnv.coding)[grep("DUP", names(cnv.coding))]
  
  cnv1 <- getSigFeatures(cnv.coding[cnv.coding$Sample %in% disc.set, ], "gene_count_DUP", covariates, cnv.dup[-c(1:3)])
  cnv2 <- getSigFeatures(cnv.coding[cnv.coding$Sample %in% disc.set, ], "gene_count_DEL", covariates, cnv.del[-c(1:3)])
  # cnv3 <- getSigFeatures(cnv.coding[cnv.coding$Sample %in% disc.set, ], "gene_count_PartialDup", covariates, names(cnv.coding)[85:121])
  
  cnv.c.sig <- c(cnv1$set, cnv2$set)
  cnv.set <- rbind(cnv1$coeff, cnv2$coeff)
  
  if(length(cnv1$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(cnv1$coeff, "variant" = "coding_rare_CNVs", stringsAsFactors = F))
  if(length(cnv2$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(cnv2$coeff, "variant" = "coding_rare_CNVs", stringsAsFactors = F))
  
  scnv.del <- names(scnv.coding)[grep("DEL", names(scnv.coding))]
  scnv.dup <- names(scnv.coding)[grep("DUP", names(scnv.coding))]
  ### smaller cnvs
  scnv1 <- getSigFeatures(scnv.coding[scnv.coding$Sample %in% disc.set, ], "gene_count_DUP", covariates, scnv.dup[-c(1:3)])
  scnv2 <- getSigFeatures(scnv.coding[scnv.coding$Sample %in% disc.set, ], "gene_count_DEL", covariates, scnv.del[-c(1:3)])

  scnv.c.sig <- c(scnv1$set, scnv2$set)
  scnv.set <- rbind(scnv1$coeff, scnv2$coeff)
  
  if(length(scnv1$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(scnv1$coeff, "variant" = "small_coding_rare_CNVs", stringsAsFactors = F))
  if(length(scnv2$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(scnv2$coeff, "variant" = "small_coding_rare_CNVs", stringsAsFactors = F))
  #################
  ## get significant cnv nc
  names(cnv.ncoding) <- gsub("deletion", "DEL", names(cnv.ncoding))
  names(cnv.ncoding) <- gsub("duplication", "DUP", names(cnv.ncoding))
  ncnv.del <- names(cnv.ncoding)[grep("DEL", names(cnv.ncoding))]
  ncnv.dup <- names(cnv.ncoding)[grep("DUP", names(cnv.ncoding))]
  ncnv1 <- getSigFeatures(cnv.ncoding, "TotalRare_DUP", covariates, ncnv.dup)
  ncnv2 <- getSigFeatures(cnv.ncoding, "TotalRare_DEL", covariates, ncnv.del)
  # cnv3 <- getSigFeatures(cnv.coding, "gene_count_PartialDup", covariates, names(cnv.coding)[85:121])
  
  cnv.nc.sig <- c(ncnv1$set, ncnv2$set)
  ncnv.set <- rbind(ncnv1$coeff, ncnv2$coeff)
  
  if(length(ncnv.set$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(ncnv1$coeff, "variant" = "noncoding_rare_CNVs", stringsAsFactors = F))
  if(length(ncnv.set$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(ncnv2$coeff, "variant" = "noncoding_rare_CNVs", stringsAsFactors = F))
  
  ## get significant smaller cnv nc
  names(scnv.ncoding) <- gsub("deletion", "DEL", names(scnv.ncoding))
  names(scnv.ncoding) <- gsub("duplication", "DUP", names(scnv.ncoding))
  sncnv.del <- names(scnv.ncoding)[grep("DEL", names(scnv.ncoding))]
  sncnv.dup <- names(scnv.ncoding)[grep("DUP", names(scnv.ncoding))]
  sncnv1 <- getSigFeatures(scnv.ncoding, "TotalRare_DUP", covariates, sncnv.dup)
  sncnv2 <- getSigFeatures(scnv.ncoding, "TotalRare_DEL", covariates, sncnv.del)
  # cnv3 <- getSigFeatures(cnv.coding, "gene_count_PartialDup", covariates, names(cnv.coding)[85:121])
  
  scnv.nc.sig <- c(sncnv1$set, sncnv2$set)
  sncnv.set <- rbind(sncnv1$coeff, sncnv2$coeff)
  
  if(length(sncnv.set$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(sncnv1$coeff, "variant" = "small_noncoding_rare_CNVs", stringsAsFactors = F))
  if(length(sncnv.set$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(sncnv2$coeff, "variant" = "small_noncoding_rare_CNVs", stringsAsFactors = F))
  
  ##########
  
  snv1 <- getSigFeatures(snv.coding.rare[snv.coding.rare$Sample %in% disc.set, ], "Total_lof", covariates, names(snv.coding.rare)[4:43])
  snv2 <- getSigFeatures(snv.coding.rare[snv.coding.rare$Sample %in% disc.set, ], "Total_tier1_ms", covariates, names(snv.coding.rare)[c(45:84, 126)])
  snv3 <- getSigFeatures(snv.coding.rare[snv.coding.rare$Sample %in% disc.set, ], "Total_tier2_ms", covariates, names(snv.coding.rare)[86:125])
  snv.c.rare.sig <- c(snv1$set, snv2$set, snv3$set)
  rare.set <- rbind(snv1$coeff, rbind(snv2$coeff, snv3$coeff))
  
  if(length(snv1$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(snv1$coeff, "variant" = "coding_rare_SNVs", stringsAsFactors = F))
  if(length(snv2$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(snv2$coeff, "variant" = "coding_rare_SNVs", stringsAsFactors = F))
  if(length(snv3$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(snv3$coeff, "variant" = "coding_rare_SNVs", stringsAsFactors = F))
  
  snv1 <- getSigFeatures(snv.coding.denovo[snv.coding.denovo$Sample %in% disc.set, ], "Total_lof", covariates, names(snv.coding.denovo)[4:43])
  snv2 <- getSigFeatures(snv.coding.denovo[snv.coding.denovo$Sample %in% disc.set, ], "Total_tier1_ms", covariates, names(snv.coding.denovo)[c(45:84, 126)])
  snv3 <- getSigFeatures(snv.coding.denovo[snv.coding.denovo$Sample %in% disc.set, ], "Total_tier2_ms", covariates, names(snv.coding.denovo)[86:125])
  snv.c.denovo.sig <- c(snv1$set, snv2$set, snv3$set)
  denovo.set <- rbind(snv1$coeff, rbind(snv2$coeff, snv3$coeff))
  if(length(snv1$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(snv1$coeff, "variant" = "coding_denovo_SNVs", stringsAsFactors = F))
  if(length(snv2$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(snv2$coeff, "variant" = "coding_denovo_SNVs", stringsAsFactors = F))
  if(length(snv3$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(snv3$coeff, "variant" = "coding_denovo_SNVs", stringsAsFactors = F))
  
  snv.nc.rare.sig <- getSigFeatures(snv.nc.rare[snv.nc.rare$Sample %in% disc.set, ], "TotalRare", covariates, names(snv.nc.rare)[2:26])
  if(length(snv.nc.rare.sig$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(snv.nc.rare.sig$coeff, "variant" = "noncoding_rare_SNVs", stringsAsFactors = F))
  nc.rare.set <- snv.nc.rare.sig$coeff
  snv.nc.rare.sig <- snv.nc.rare.sig$set
  
  snv.nc.denovo.sig <- getSigFeatures(snv.nc.denovo[snv.nc.denovo$Sample %in% disc.set, ], "TotalRare", covariates, names(snv.nc.denovo)[2:26])
  if(length(snv.nc.denovo.sig$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(snv.nc.denovo.sig$coeff, "variant" = "noncoding_denovo_SNVs", stringsAsFactors = F))
  nc.denovo.set <- snv.nc.denovo.sig$coeff
  snv.nc.denovo.sig <- snv.nc.denovo.sig$set
  
  
  # 
  # all.event <- data.frame()
  # if(length(cnv.c.sig) > 0){
  #   cnv.event <- getCNVEvents(all.cnv, cnv.set, gsNFLD)
  #   names(cnv.event)[-1] <- paste(names(cnv.event)[-1], "rare_CNVs", sep="_")
  #   if(nrow(all.event) == 0)
  #     all.event <- cnv.event
  # }
  # if(length(snv.c.rare.sig) > 0){
  #   rare.event <- getSNVEvents(all.rare, rare.set)
  #   names(rare.event)[-1] <- paste(names(rare.event)[-1], "rare_SNVs", sep="_")
  #   if(nrow(all.event) == 0)
  #     all.event <- cnv.event
  #   else
  #     all.event <- merge(all.event, rare.event, by = "sample", all = T)
  # }
  # if(length(snv.c.denovo.sig) > 0){
  #   denovo.event <- getSNVEvents(all.denovo, denovo.set)
  #   names(denovo.event)[-1] <- paste(names(denovo.event)[-1], "denovo_SNVs", sep="_")
  #   if(nrow(all.event) == 0)
  #     all.event <- cnv.event
  #   else
  #     all.event <- merge(all.event, denovo.event, by = "sample", all = T)
  # }
  # 
  # all.event[is.na(all.event)] <- 0
  # if(length(grep("increasing_severity_events", names(all.event))) > 1)
  #   all.event$all_increasing_severity_events <- rowSums(all.event[, grep("increasing_severity_events", names(all.event))], na.rm = T)
  # if(length(grep("increasing_severity_events", names(all.event))) == 1)
  #   all.event$all_increasing_severity_events <- all.event[, grep("increasing_severity_events", names(all.event))]
  # if(length(grep("increasing_severity_events", names(all.event))) == 0)
  #   all.event$all_increasing_severity_events <- 0
  # 
  # if(length(grep("decreasing_severity_events", names(all.event))) > 1)
  #   all.event$all_decreasing_severity_events <- rowSums(all.event[, grep("decreasing_severity_events", names(all.event))], na.rm = T)
  # if(length(grep("decreasing_severity_events", names(all.event))) == 1)
  #   all.event$all_decreasing_severity_events <- all.event[, grep("decreasing_severity_events", names(all.event))]
  # if(length(grep("decreasing_severity_events", names(all.event))) == 0)
  #   all.event$all_decreasing_severity_events <- 0
  # 
  # 
  # if(length(grep("increasing_severity_variants", names(all.event))) > 1)
  #   all.event$all_increasing_severity_variants <- rowSums(all.event[, grep("increasing_severity_variants", names(all.event))], na.rm = T)
  # if(length(grep("increasing_severity_variants", names(all.event))) == 1)
  #   all.event$all_increasing_severity_variants <- all.event[, grep("increasing_severity_variants", names(all.event))]
  # if(length(grep("increasing_severity_variants", names(all.event))) == 0)
  #   all.event$all_increasing_severity_variants <- 0
  # 
  # if(length(grep("decreasing_severity_variants", names(all.event))) > 1)
  #   all.event$all_decreasing_severity_variants <- rowSums(all.event[, grep("decreasing_severity_variants", names(all.event))], na.rm = T)
  # if(length(grep("decreasing_severity_variants", names(all.event))) == 1)
  #   all.event$all_decreasing_severity_variants <- all.event[, grep("decreasing_severity_variants", names(all.event))]
  # if(length(grep("decreasing_severity_variants", names(all.event))) == 0)
  #   all.event$all_decreasing_severity_variants <- 0
  # 
  
  sig.set <- c(sig.set, 
               paste("CNVs",cnv.c.sig, sep=":"), 
               paste("small CNVs",scnv.c.sig, sep=":"), 
               paste("CNVs",cnv.nc.sig, sep=":"), 
               paste("small CNVs",scnv.nc.sig, sep=":"), 
               paste("rare SNVs", snv.c.rare.sig, sep=":"),
               paste("de novo SNVs", snv.c.denovo.sig, sep=":"), 
               paste("rare SNVs", snv.nc.rare.sig, sep=":"), 
               paste("de novo SNVs", snv.nc.denovo.sig, sep=":"))
  
  for(p in c(1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001)){
    
    # cnv.c.score <- getScoreZ(cnv.coding[cnv.coding$Sample %in% targ.set, ], cnv.set[cnv.set$pvalue < p, ], "coding_CNVs")
    # scnv.c.score <- getScoreZ(scnv.coding[scnv.coding$Sample %in% targ.set, ], scnv.set[scnv.set$pvalue < p,], "small_coding_CNVs")
    # cnv.nc.score <- getScoreZ(cnv.ncoding[cnv.ncoding$Sample %in% targ.set, ], ncnv.set[ncnv.set$pvalue < p,], "noncoding_CNVs")
    # scnv.nc.score <- getScoreZ(scnv.ncoding[scnv.ncoding$Sample %in% targ.set, ], sncnv.set[sncnv.set$pvalue < p,], "small_noncoding_CNVs")
    # snv.c.rare.score <- getScoreZ(snv.coding.rare[snv.coding.rare$Sample %in% targ.set, ], rare.set[rare.set$pvalue < p,], "coding_rare_SNVs")
    # snv.c.denovo.score <- getScoreZ(snv.coding.denovo[snv.coding.denovo$Sample %in% targ.set, ], denovo.set[denovo.set$pvalue < p,], "coding_denovo_SNVs")
    # snv.nc.rare.score <- getScoreZ(snv.nc.rare[snv.nc.rare$Sample %in% targ.set, ], nc.rare.set[nc.rare.set$pvalue < p,], "noncoding_rare_SNVs")
    # snv.nc.denovo.score <- getScoreZ(snv.nc.denovo[snv.nc.denovo$Sample %in% targ.set, ], nc.denovo.set[nc.denovo.set$pvalue < p,], "noncoding_denovo_SNVs")
    
    cnv.c.score <- getScore(cnv.coding[cnv.coding$Sample %in% disc.set, ],
                            cnv.coding[cnv.coding$Sample %in% targ.set, ], baseline, cnv.c.sig[cnv.set$pvalue < p], "coding_CNVs")

    scnv.c.score <- getScore(scnv.coding[scnv.coding$Sample %in% disc.set, ],
                            scnv.coding[scnv.coding$Sample %in% targ.set, ], baseline, scnv.c.sig[scnv.set$pvalue < p], "small_coding_CNVs")

    cnv.nc.score <- getScore(cnv.ncoding[cnv.ncoding$Sample %in% disc.set, ],
                            cnv.ncoding[cnv.ncoding$Sample %in% targ.set, ], baseline, cnv.nc.sig[ncnv.set$pvalue < p], "noncoding_CNVs")

    scnv.nc.score <- getScore(scnv.ncoding[scnv.ncoding$Sample %in% disc.set, ],
                             scnv.ncoding[scnv.ncoding$Sample %in% targ.set, ], baseline, scnv.nc.sig[sncnv.set$pvalue < p], "small_noncoding_CNVs")

    snv.c.rare.score <- getScore(snv.coding.rare[snv.coding.rare$Sample %in% disc.set, ],
                                 snv.coding.rare[snv.coding.rare$Sample %in% targ.set, ], baseline, snv.c.rare.sig[rare.set$pvalue < p], "coding_rare_SNVs")

    snv.c.denovo.score <- getScore(snv.coding.denovo[snv.coding.denovo$Sample %in% disc.set, ],
                            snv.coding.denovo[snv.coding.denovo$Sample %in% targ.set, ], baseline, snv.c.denovo.sig[denovo.set$pvalue < p], "coding_denovo_SNVs")

    snv.nc.rare.score <- getScore(snv.nc.rare[snv.nc.rare$Sample %in% disc.set, ],
                            snv.nc.rare[snv.nc.rare$Sample %in% targ.set, ], baseline, snv.nc.rare.sig[nc.rare.set$pvalue < p], "noncoding_rare_SNVs")

    snv.nc.denovo.score <- getScore(snv.nc.denovo[snv.nc.denovo$Sample %in% disc.set, ],
                            snv.nc.denovo[snv.nc.denovo$Sample %in% targ.set, ], baseline, snv.nc.denovo.sig[nc.denovo.set$pvalue < p], "noncoding_denovo_SNVs")
    
    final.score <- merge(merge(merge(merge(merge(merge(merge(cnv.c.score, scnv.c.score, by = "Sample", all = T), 
                                           snv.c.rare.score, by = "Sample", all = T), snv.c.denovo.score, by = "Sample", all = T), 
                               snv.nc.rare.score, by = "Sample", all = T), 
                         snv.nc.denovo.score, by = "Sample", all = T), cnv.nc.score, all = T), scnv.nc.score, all = T)
    
    final.score$total_risk <- rowSums(final.score[, -1], na.rm = T)
    # all.event$sample <- gsub("-", "_", all.event$sample)
    # final.score <- merge(final.score, all.event[, c("sample", "all_increasing_severity_events", "all_decreasing_severity_events",
    #                                                 "all_increasing_severity_variants", "all_decreasing_severity_variants")], by.x = "Sample", by.y = "sample", all.x = T)
    final.score$pvalue <- p
    # dt.plot <- merge(final.score, samples, by = "Sample", all.x = T)
    # dt.plot$DysmorphologyClassificaiton <- ifelse(dt.plot$MultiClass == 0, "Essential", "Complex")
    # dt.plot$DysmorphologyClassificaiton[dt.plot$MultiClass == 1] <- "Equivocal"
    # 
    # ggplot(dt.plot, aes(x = DysmorphologyClassificaiton, y = total_risk, fill = DysmorphologyClassificaiton)) +
    #   geom_boxplot() + theme_classic() + theme(legend.position = "none", axis.title.x = element_blank()) + coord_cartesian(ylim = c(-0.5, 2.5))
    #   
    score.out <- rbind(score.out, final.score)
  }
}
}

write.table(coeff.all, "ADM/coeff.from.30x.10fold.tsv", sep="\t", row.names=F, quote=F, col.names=T)
write.table(score.out, "ADM/score.from.30x.10fold.tsv", sep="\t", row.names=F, quote = F, col.names=T)
write.table(sig.set, "ADM/significant.set.from.30x.10fold.tsv", sep="\t", row.names=F, quote=F, col.names=T)
write.table(samples, "ADM/sample.with.score.tsv", sep="\t", row.names=F, quote=F, col.names=T)


##############
##############
samples <- read.delim("ADM/sample.with.score.tsv", stringsAsFactors = F)
denovo <- readLines("samples.with.denovo.txt")
samples <- samples[samples$Sample %in% denovo, ]

out.perf <- data.frame()
for(p in c(1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001)){
  score.out <- read.delim("ADM/score.from.30x.10fold.tsv", stringsAsFactors = F)
  score.out <- score.out[score.out$pvalue == p, ]
  # sig.set <- read.delim("ADM/significant.set.from.30x.10fold.tsv", stringsAsFactors = F)
  coeff.df <- read.delim("ADM/coeff.from.30x.10fold.tsv", stringsAsFactors = F)
  coeff.df <- coeff.df[coeff.df$pvalue < p, ]
  # coeff.df <- coeff.df[grep("duplication", coeff.df$set), ]
  agg.coeff.df <- merge(aggregate(.~ set + variant, coeff.df[-3], sum), as.data.frame(table(coeff.df$set)), by.x = "set", by.y = "Var1", all = T)
  agg.coeff.df$coeff <- agg.coeff.df$coeff/300
  agg.coeff.df
  tmp.out <- data.frame()
  for(s in unique(score.out$Sample)){
    tmp <- score.out[score.out$Sample == s, ]
    tmp <- colSums(tmp[, -1], na.rm = T)/nrow(tmp)
    
    tmp.out <- rbind(tmp.out, data.frame("Sample" = s, t(tmp)))
  }
  
  score.out <- tmp.out
  score.out <- merge(score.out, samples, by = "Sample", all = F)
  score.out$ADM <- factor(score.out$ADM)
  
  lm <- glm(ADM ~ total_risk, score.out, family = binomial(link = "logit"))
  
  # score.out$DysmorphologyClassificaiton <- ifelse(score.out$MultiClass == 0, "Essential", "Complex")
  # score.out$DysmorphologyClassificaiton[score.out$MultiClass == 1] <- "Equivocal"
  suffix <- ""
  
  # cli.snv <- readxl::read_excel("../data/clinicallysignificant_var_forBE.xlsx", sheet = 1)
  # cli.cnv <- readxl::read_excel("../data/clinicallysignificant_var_forBE.xlsx", sheet = 2)
  # score.out <- score.out[!score.out$Sample %in% gsub("-", "_", c(as.character(cli.snv$Case), as.character(cli.cnv$Sample))), ]
  # suffix <- ".noCRV"
  # score.out$total_risk <- rowSums(score.out[, 2:4], na.rm = T)
  # write.table(score.out, "score.sample.withno.CRV.tsv", sep="\t", row.names=F, quote=F, col.names=T)
  #write.table(score.out[order(score.out$total_risk, decreasing = T), ], "cross.validate.score.tsv", sep="\t", row.names=F, quote=F, col.names=T)
  perc.rank <- function(x) trunc(rank(x))/length(x)
  
  ### calculate percentile
  score.out$coding_CNVs_perc <- perc.rank(score.out$coding_CNVs)
  score.out$small_coding_CNVs_perc <- perc.rank(score.out$small_coding_CNVs)
  score.out$noncoding_CNVs_perc <- perc.rank(score.out$noncoding_CNVs)
  score.out$small_noncoding_CNVs_perc <- perc.rank(score.out$small_noncoding_CNVs)
  score.out$coding_rare_SNVs_perc <- perc.rank(score.out$coding_rare_SNVs)
  score.out$coding_denovo_SNVs_perc <- perc.rank(score.out$coding_denovo_SNVs)
  score.out$noncoding_rare_SNVs_perc <- perc.rank(score.out$noncoding_rare_SNVs)
  score.out$noncoding_denovo_SNVs_perc <- perc.rank(score.out$noncoding_denovo_SNVs)
  
  score.out$total_risk_perc <- perc.rank(score.out$total_risk)
  
  # ####################
  # ### trend
  # p1 <- ggplot(score.out, aes(x = all_increasing_severity_events, y = total_risk)) + 
  #   geom_smooth(data = score.out[score.out$ADM == "1",], method = "lm", fullrange=T, color = "red", fill = "red", alpha = .25, lty = 2, lwd = .4) + 
  #   geom_smooth(data = score.out[score.out$ADM == "0",], method = "lm", fullrange=T, color = "green", fill = "green", alpha = .25, lty = 2, lwd = .4) + 
  #   geom_point(shape = 1, aes(color = ADM, fill = ADM), show.legend = F, size = 1.5) + theme_bw() + 
  #   xlab("#genes impacted in positively correlated sets") + ylab("total score")
  # p2 <- ggplot(score.out, aes(x = all_decreasing_severity_events, y = total_risk)) + 
  #   geom_smooth(data = score.out[score.out$ADM == "1",], method = "lm", fullrange=T, color = "red", fill = "red", alpha = .25, lty = 2, lwd = .4) + 
  #   geom_smooth(data = score.out[score.out$ADM == "0",], method = "lm", fullrange=T, color = "green", fill = "green", alpha = .25, lty = 2, lwd = .4) + 
  #   geom_point(shape = 1, aes(color = ADM, fill = ADM), size = 1.5) + theme_bw() +
  #   xlab("#genes impacted in negatively correlated sets") + theme(legend.position = "bottom") + ylab("total score")
  # plot_grid(p1, p2, nrow =2, rel_heights = c(0.9, 1))
  # ggsave(sprintf("ADM/score.genes%s.pdf", suffix), width = 7, height = 10)
  # 
  # ### trend
  # p1 <- ggplot(score.out, aes(x = all_increasing_severity_variants, y = total_risk)) + 
  #   geom_smooth(data = score.out[score.out$ADM == "1",], method = "lm", fullrange=T, color = "red", fill = "red", alpha = .25, lty = 2, lwd = .4) + 
  #   geom_smooth(data = score.out[score.out$ADM == "0",], method = "lm", fullrange=T, color = "green", fill = "green", alpha = .25, lty = 2, lwd = .4) + 
  #   geom_point(shape = 1, aes(color = ADM, fill = ADM), show.legend = F, size = 1.5) + theme_bw() + 
  #   xlab("#variants in positively correlated sets") + ylab("total score")
  # p2 <- ggplot(score.out, aes(x = all_decreasing_severity_variants, y = total_risk)) + 
  #   geom_smooth(data = score.out[score.out$ADM == "1",], method = "lm", fullrange=T, color = "red", fill = "red", alpha = .25, lty = 2, lwd = .4) + 
  #   geom_smooth(data = score.out[score.out$ADM == "0",], method = "lm", fullrange=T, color = "green", fill = "green", alpha = .25, lty = 2, lwd = .4) + 
  #   geom_point(shape = 1, aes(color = ADM, fill = ADM), size = 1.5) + theme_bw() +
  #   xlab("#variants in negatively correlated sets") + theme(legend.position = "bottom") + ylab("total score")
  # plot_grid(p1, p2, nrow =2, rel_heights = c(0.9, 1))
  # ggsave(sprintf("ADM/score.variants%s.pdf", suffix), width = 7, height = 10)
  
  
  ########### plot ###########
  p1.pvalue <- paste0(sprintf("p=%s", formatC(wilcox.test(score.out$coding_CNVs[score.out$ADM == "1"],
                                   score.out$coding_CNVs[score.out$ADM == "0"], "greater")$p.value, 
                       format = "e", digits = 1)))
  
  p0.pvalue <- paste0(sprintf("p=%s", formatC(wilcox.test(score.out$small_coding_CNVs[score.out$ADM == "1"],
                                                          score.out$small_coding_CNVs[score.out$ADM == "0"], "greater")$p.value, 
                                              format = "e", digits = 1)))
  
  p2.pvalue <-  paste0(sprintf("p=%s", formatC(wilcox.test(score.out$coding_rare_SNVs[score.out$ADM == "1"],
                                                                   score.out$coding_rare_SNVs[score.out$ADM == "0"], "greater")$p.value, 
                                                       format = "e", digits = 1)))
  
  p3.pvalue <-  paste0(sprintf("p=%s", formatC(wilcox.test(score.out$coding_denovo_SNVs[score.out$ADM == "1"],
                                                                   score.out$coding_denovo_SNVs[score.out$ADM == "0"], "greater")$p.value, 
                                                       format = "e", digits = 1)))
  
  p4.pvalue <-  paste0(sprintf("p=%s", formatC(wilcox.test(score.out$noncoding_rare_SNVs[score.out$ADM == "1"],
                                                                   score.out$noncoding_rare_SNVs[score.out$ADM == "0"], "greater")$p.value, 
                                                       format = "e", digits = 1)))
  
  p5.pvalue <- paste0(sprintf("p=%s", formatC(wilcox.test(score.out$noncoding_denovo_SNVs[score.out$ADM == "1"],
                                                                  score.out$noncoding_denovo_SNVs[score.out$ADM == "0"], "greater")$p.value, 
                                                      format = "e", digits = 1)))
  
  p6.pvalue <-  paste0(sprintf("p=%s", formatC(wilcox.test(score.out$total_risk[score.out$ADM == "1"],
                                                                   score.out$total_risk[score.out$ADM == "0"], "greater")$p.value, 
                                                       format = "e", digits = 1)))
  
  p7.pvalue <-  paste0(sprintf("p=%s", formatC(wilcox.test(score.out$noncoding_CNVs[score.out$ADM == "1"],
                                                           score.out$noncoding_CNVs[score.out$ADM == "0"], "greater")$p.value, 
                                               format = "e", digits = 1)))
  
  p8.pvalue <-  paste0(sprintf("p=%s", formatC(wilcox.test(score.out$small_noncoding_CNVs[score.out$ADM == "1"],
                                                           score.out$small_noncoding_CNVs[score.out$ADM == "0"], "greater")$p.value, 
                                               format = "e", digits = 1)))
  
  out.perf <- rbind(out.perf, data.frame("cutoff" = p,
                                         "NagelkerkeR2" = signif(NagelkerkeR2(lm)$R2, digits = 3),
                                         "p" = p6.pvalue))
  
  message(sprintf("%s's R2 = %s, %s", p, signif(NagelkerkeR2(lm)$R2, digits = 3), p6.pvalue))
  
  score.out$ADM <- ifelse(score.out$ADM == 1, "dys", "nondys")
  score.out$ADM <- factor(score.out$ADM, levels = c("nondys", "dys"))
  
  p0 <- ggplot(score.out, aes(x = ADM, y = small_coding_CNVs_perc, fill = ADM)) +
    geom_violin() + geom_boxplot(fill = "white", alpha = .75, width = .2) + theme_classic() + 
    theme(legend.position = "none", axis.title.x = element_blank()) + 
    annotate("text", label = p0.pvalue, x = c(1.5), y = c(1.12), cex = 5) +
    # annotate("text", label = c("**"), x = c(1.5), y = c(1.12), cex = 5) +
    annotate("line", x = c(1, 2), y = c(1.05, 1.05), group = c("A", "A")) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + ylab("Small coding CNVs")
  
  p1 <- ggplot(score.out, aes(x = ADM, y = coding_CNVs_perc, fill = ADM)) +
    geom_violin() + geom_boxplot(fill = "white", alpha = .75, width = .2) + theme_classic() + 
    theme(legend.position = "none", axis.title.x = element_blank()) + 
    annotate("text", label = p1.pvalue, x = c(1.5), y = c(1.12), cex = 5) +
    # annotate("text", label = c("**"), x = c(1.5), y = c(1.12), cex = 5) +
    annotate("line", x = c(1, 2), y = c(1.05, 1.05), group = c("A", "A")) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + ylab("Coding CNVs")
  
  p01 <- ggplot(score.out, aes(x = ADM, y = small_noncoding_CNVs_perc, fill = ADM)) +
    geom_violin() + geom_boxplot(fill = "white", alpha = .75, width = .2) + theme_classic() + 
    theme(legend.position = "none", axis.title.x = element_blank()) + 
    annotate("text", label = p8.pvalue, x = c(1.5), y = c(1.12), cex = 5) +
    # annotate("text", label = c("**"), x = c(1.5), y = c(1.12), cex = 5) +
    annotate("line", x = c(1, 2), y = c(1.05, 1.05), group = c("A", "A")) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + ylab("Small noncoding CNVs")
  
  p12 <- ggplot(score.out, aes(x = ADM, y = noncoding_CNVs_perc, fill = ADM)) +
    geom_violin() + geom_boxplot(fill = "white", alpha = .75, width = .2) + theme_classic() + 
    theme(legend.position = "none", axis.title.x = element_blank()) + 
    annotate("text", label = p7.pvalue, x = c(1.5), y = c(1.12), cex = 5) +
    # annotate("text", label = c("**"), x = c(1.5), y = c(1.12), cex = 5) +
    annotate("line", x = c(1, 2), y = c(1.05, 1.05), group = c("A", "A")) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + ylab("Noncoding CNVs")
  
  p2 <- ggplot(score.out, aes(x = ADM, y = coding_rare_SNVs_perc, fill = ADM)) +
    geom_violin() + geom_boxplot(fill = "white", alpha = .75, width = .2) + theme_classic() + theme(legend.position = "none", axis.title.x = element_blank())+ 
    annotate("text", label = p2.pvalue, x = c(1.5), y = c(1.12), cex = 5) +
    annotate("line", x = c(1, 2), y = c(1.05, 1.05), group = c("A", "A")) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + ylab("Coding rare SNVs")
  
  p3 <- ggplot(score.out, aes(x = ADM, y = coding_denovo_SNVs_perc, fill = ADM)) +
    geom_violin() + geom_boxplot(fill = "white", alpha = .75, width = .2) + theme_classic() + 
    theme(legend.position = "none", axis.title.x = element_blank()) + 
    annotate("text", label = p3.pvalue, x = c(1.5), y = c(1.12), cex = 5) +
    annotate("line", x = c(1, 1), y = c(1.05, 1.05), group = c("A", "A")) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + ylab("Coding de novo SNVs")
  
  p4 <- ggplot(score.out, aes(x = ADM, y = noncoding_rare_SNVs_perc, fill = ADM)) +
    geom_violin() + geom_boxplot(fill = "white", alpha = .75, width = .2) + theme_classic() + 
    theme(legend.position = "none", axis.title.x = element_blank()) + 
    annotate("text", label = p4.pvalue, x = c(1.5), y = c(1.12), cex = 5) +
    annotate("line", x = c(1, 1), y = c(1.05, 1.05), group = c("A", "A")) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + ylab("Noncoding rare SNVs")
  
  p5 <- ggplot(score.out, aes(x = ADM, y = noncoding_denovo_SNVs_perc, fill = ADM)) +
    geom_violin() + geom_boxplot(fill = "white", alpha = .75, width = .2) + theme_classic() + 
    theme(legend.position = "none", axis.title.x = element_blank()) + 
    annotate("text", label = p5.pvalue, x = c(1.5), y = c(1.12), cex = 5) +
    annotate("line", x = c(1, 1), y = c(1.05, 1.05), group = c("A", "A")) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + ylab("Noncoding denovo SNVs")
  
  p6 <- ggplot(score.out, aes(x = ADM, y = total_risk_perc, fill = ADM)) +
    geom_violin() + geom_boxplot(fill = "white", alpha = .75, width = .2) + theme_classic() + 
    theme(legend.position = "none", axis.title.x = element_blank()) + 
    annotate("text", label = p6.pvalue, x = c(1.5), y = c(1.12), cex = 5) +
    annotate("line", x = c(1, 1), y = c(1.05, 1.05), group = c("A", "A")) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + ylab("Total risk score")
  
  if(suffix == ""){
    plot_grid(p1, p0, p12, p01, p2, p3, p4, p5, p6, nrow = 3)
    width = 12
    height = 8
  }else{
    p6
    width = 5
    height = 4
  }
  
  ggsave(sprintf("ADM/NFLD.%s_30x10fold.coding%sall.pdf", p, suffix), width = width, height = height)
}

## updated May 2021
# 1's R2 = 0.00124, p=3.4e-01
# 0.5's R2 = 0.0554, p=1.8e-03
# 0.1's R2 = 0.161, p=3.6e-06
# 0.05's R2 = 0.127, p=3.8e-05
# 0.01's R2 = 0.0321, p=3.7e-02
# 0.005's R2 = 0.00577, p=3.0e-02
# 0.001's R2 = 0.0018, p=8.4e-01


# 1's R2 = 0.000727, p=1.6e-01
# 0.5's R2 = 0.029, p=3.2e-02
# 0.1's R2 = 0.0521, p=9.5e-03
# 0.05's R2 = 0.0581, p=1.4e-02
# 0.01's R2 = 0.027, p=6.1e-02
# 0.005's R2 = 0.00257, p=1.9e-01
# 0.001's R2 = 0.00448, p=9.3e-01
write.table(out.perf, "performance.all.cutoffs.tsv", sep="\t", row.names=F)
out.perf <- read.delim("performance.all.cutoffs.tsv", stringsAsFactors = F)
out.perf$cutoff <- factor(out.perf$cutoff, levels = rev(out.perf$cutoff))
ggplot(out.perf, aes(x = paste0("p<",cutoff), y = NagelkerkeR2, label = p)) + 
  geom_bar(stat = "identity", color = "black", fill = "red") + xlab("p cut-off") + ylim(0, 0.23) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_text(angle = 90, nudge_y = 0.025)
ggsave("adm.performance.pdf", width = 4, height = 5)
# tmp.pvalue <- paste(paste0("coding CNVs\t", gsub("\n", "\t", p1.pvalue)),
#                     paste0("coding rare SNVs\t", gsub("\n", "\t", p2.pvalue)),
#                     paste0("coding denovo SNVs\t", gsub("\n", "\t", p3.pvalue)),
#                     paste0("noncoding rare SNVs\t", gsub("\n", "\t", p4.pvalue)),
#                     paste0("risk score\t", gsub("\n", "\t", p6.pvalue)), sep="\n")
# writeLines(tmp.pvalue, sprintf("pvalue.all%s.tsv", suffix))
# 
# sig.set.count <- data.frame(table(sig.set), stringsAsFactors = F)
# sig.set.count$set <- sapply(sapply(as.character(sig.set.count$sig.set), strsplit, ":"), "[", 2)
# sig.set.count$type <- sapply(sapply(as.character(sig.set.count$sig.set), strsplit, ":"), "[", 1)
# sig.set.count <- na.omit(sig.set.count)
# sig.set.count <- sig.set.count[order(sig.set.count$Freq, decreasing = T), ]
# 
# sig.set.count$set <- factor(sig.set.count$set, levels = unique(c(sig.set.count$set[sig.set.count$type == "CNVs"], 
#                                                           sig.set.count$set[sig.set.count$type == "rare SNVs"],
#                                                           sig.set.count$set[sig.set.count$type == "de novo SNVs"])))
# sig.set.count$region <- ifelse(sig.set.count$set %in% names(snv.nc.rare), "noncoding", "coding")
# 
# ggplot(sig.set.count, aes(x = set, y = Freq, fill = type)) + geom_bar(stat = "identity") + xlab("") + theme_classic() + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) + facet_grid(region~.)
# ggsave(sprintf("ADM/30x10fold_selected_genesets.coding%s.pdf", suffix), width = 7, height =5)
# 
# ### boxplot
# inc <- ggplot(score.out, aes(x = ADM, y = all_increasing_severity_variants)) + 
#   geom_boxplot(outlier.alpha = 0) + theme(axis.title.x = element_blank()) + ylab("#variants in positively correlated sets") +
#   geom_jitter(aes(fill = total_risk), shape = 21, size = 4, alpha = .8) + scale_fill_gradientn(values = c(0, 0.45, 0.9, 1.1), 
#                                                                                                colours = c("white", "yellow", "red", "red"))
# dec <- ggplot(score.out, aes(x = ADM, y = all_decreasing_severity_variants)) + 
#   geom_boxplot(outlier.alpha = 0) + theme(axis.title.x = element_blank()) + ylab("#variants in negatively correlated sets") +
#   geom_jitter(aes(fill = total_risk), shape = 21, size = 4, alpha = .8) + scale_fill_gradientn(values = c(0, 0.45, 0.9, 1.1), 
#                                                                                                colours = c("white", "yellow", "red", "red"))
# 
# plot_grid(inc, dec, nrow = 1)
# ggsave(sprintf("ADM/boxplot.numvariant%s.pdf",suffix), width = 14, height =7)
# 
# fig1 <- ggplot(score.out, aes(x = ADM, y = total_risk_perc, fill = ADM)) +
#   geom_violin(color = NA) + geom_boxplot(fill = "white", alpha = .75, width = .2) + theme_classic() + 
#   theme(legend.position = "none", axis.title.x = element_blank()) + 
#   annotate("text", label = c("**", "*"), x = c(2, 2.5), y = c(1.22, 1.12), cex = 5) +
#   annotate("line", x = c(1, 3, 2, 3), y = c(1.2, 1.2, 1.1, 1.1), group = c("A", "B", "A", "B")) +
#   scale_fill_manual(values = c("#D55E00", "#F0E442", "#0072B2")) +
#   scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + ylab("Total risk score (percentile)")
# 
# agg.coeff.df$variant <- gsub("_", " ", agg.coeff.df$variant)
# coeff <- ggplot(agg.coeff.df, aes(x= variant, y=log(coeff))) + geom_boxplot() + theme_classic() + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = c("#D55E00", "#F0E442", "#0072B2")) +
#   scale_x_discrete(limits = c("coding rare CNVs", "noncoding denovo SNVs",  "coding denovo SNVs", "coding rare SNVs"))
# 
# fig2 <- ggplot(score.out, aes(x = all_increasing_severity_variants, y = total_risk)) + 
#   geom_smooth(data = score.out[score.out$DysmorphologyClassificaiton == "Complex",], method = "lm", fullrange=T, color = "#D55E00", fill = "#D55E00", alpha = .35, lty = 2, lwd = .4) + 
#   geom_smooth(data = score.out[score.out$DysmorphologyClassificaiton == "Equivocal",], method = "lm", fullrange=T, color = "#F0E442", fill = "#F0E442", alpha = .35, lty = 2, lwd = .4) + 
#   geom_smooth(data = score.out[score.out$DysmorphologyClassificaiton == "Essential",], method = "lm", fullrange=T, color = "#0072B2", fill = "#0072B2", alpha = .35, lty = 2, lwd = .4) + 
#   geom_point(shape = 21, color = "black", aes(fill = DysmorphologyClassificaiton), size = 1.5, alpha = .7) + theme_bw() +
#   scale_fill_manual(values = c("#D55E00", "#F0E442", "#0072B2")) +
#   xlab("#variants in positively correlated sets") + theme(legend.position = "bottom", legend.title = element_blank()) + ylab("total score")
# 
# # fig3 <- ggplot(score.out, aes(x = all_decreasing_severity_variants, y = total_risk)) + 
# #   geom_smooth(data = score.out[score.out$DysmorphologyClassificaiton == "Complex",], method = "lm", fullrange=T, color = "red", fill = "red", alpha = .25, lty = 2, lwd = .4) + 
# #   geom_smooth(data = score.out[score.out$DysmorphologyClassificaiton == "Equivocal",], method = "lm", fullrange=T, color = "green", fill = "green", alpha = .25, lty = 2, lwd = .4) + 
# #   geom_smooth(data = score.out[score.out$DysmorphologyClassificaiton == "Essential",], method = "lm", fullrange=T, color = "blue", fill = "blue", alpha = .25, lty = 2, lwd = .4) + 
# #   geom_point(shape = 1, aes(color = DysmorphologyClassificaiton, fill = DysmorphologyClassificaiton), size = 1.5) + theme_bw() +
# #   xlab("#variants in negatively correlated sets") + theme(legend.position = "bottom", legend.title = element_blank()) + ylab("total score")
# plot_grid(fig1, fig2, coeff, nrow = 1)
# ggsave("ADM/score.with.crv.pdf", height = 4, width = 12)
# 
# 
# ### other variants
# 
# dt <- read.delim("../snvAnalysis/rare.snv.signal.tsv", stringsAsFactors = F)
# dt <- dt[dt$Sample %in% c("3_0269_000", "3_0406_000", "3_0391_000", "3_0392_000", "3_0368_000", "3_0207_000"), c("Sample", "chr", "start", "end", "var_type", "gene_symbol", "type", "geneset")]
# write.table(dt, "ADM/additional_variant_selected_samples.tsv", sep="\t", row.names=F, quote=F, col.names=T)

#### plot p 0.1
p <- 0.1
score.out <- read.delim("ADM/score.from.30x.10fold.tsv", stringsAsFactors = F)
score.out <- score.out[score.out$pvalue == p, ]
samples <- read.delim("ADM/sample.with.score.tsv", stringsAsFactors = F)
denovo <- readLines("samples.with.denovo.txt")
samples <- samples[samples$Sample %in% denovo, ]
tmp.out <- data.frame()
for(s in unique(score.out$Sample)){
  tmp <- score.out[score.out$Sample == s, ]
  tmp <- colSums(tmp[, -1], na.rm = T)/nrow(tmp)
  
  tmp.out <- rbind(tmp.out, data.frame("Sample" = s, t(tmp)))
}

score.out <- tmp.out
score.out <- merge(score.out, samples, by = "Sample", all = F)
score.out$ADM <- factor(score.out$ADM)

lm <- glm(ADM ~ total_risk, score.out, family = binomial(link = "logit"))

perc.rank <- function(x) trunc(rank(x))/length(x)
score.out$total_risk_perc <- perc.rank(score.out$total_risk)
p6.pvalue <-  paste0(sprintf("p=%s", formatC(wilcox.test(score.out$total_risk[score.out$ADM == "1"],
                                                         score.out$total_risk[score.out$ADM == "0"], "greater")$p.value, 
                                             format = "e", digits = 1)))
write.table(score.out, "nfld.adm.grs.p0.1.tsv", sep="\t",row.names=F, quote=F, col.names=T)

score.out$adm.dysmorphic <- as.character(score.out$ADM)
score.out$adm.dysmorphic <- ifelse(score.out$adm.dysmorphic == 1, "ADM\nDysmorphic", "ADM\nNondysmorphic")
score.out$adm.dysmorphic <- factor(score.out$adm.dysmorphic, levels = c("Unaffected sibling", "ADM\nNondysmorphic", "ADM\nDysmorphic"))

ggplot(score.out, aes(x = adm.dysmorphic, y = total_risk_perc, fill = adm.dysmorphic)) + geom_violin(color = NA) +
  annotate("text", label = p6.pvalue, x = c(1.5), y = c(1.12), cex = 5) +
  annotate("line", x = c(1, 2), y = c(1.05, 1.05), group = c("A", "A")) +
  geom_boxplot(fill = "white", width = .2, alpha=.6) + theme_classic() +     
  theme(legend.position = "none", axis.title.x = element_blank()) + ylab("RSV percentile") + 
  scale_fill_manual(values = c("#009eed", "#F28602"))
ggsave("nfld.adm.grs.pdf", width = 4, height = 4)

