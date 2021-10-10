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

getNCEvents <- function(snv, set.sig){
  set.tmp <- set.sig
  snv[is.na(snv)] <- 0
  dt.out <- data.frame()
  for(j in 1:2){
    if(j == 1){
      set.sig <- set.tmp[set.tmp$coeff > 0, ]
    }else{
      set.sig <- set.tmp[set.tmp$coeff < 0, ]
    }
    
    
    set <- set.sig$set
    if(length(set) > 0){
      
      snv.tmp <- snv[, c("Sample", set)]
      snv.tmp[is.na(snv.tmp)] <- 0
      sample <- data.frame()
      if(!is.vector(snv.tmp)){
        if(ncol(snv.tmp) == 2){
          snv.tmp <- snv.tmp[snv.tmp[, 2] > 0, ]
          snv.tmp$allcount <- snv.tmp[, 2]
          snv.tmp$varcount <- as.numeric(snv.tmp[, 2] > 0)
        }else{
          snv.tmp <- snv.tmp[rowSums(snv.tmp[, -1]) > 0, ]
          snv.tmp$allcount <- rowSums(snv.tmp[, -1])
          snv.tmp$varcount <- rowSums(snv.tmp[, -1]) > 0
        }
        sample <- snv.tmp[, c("Sample", "allcount", "varcount")]
      }
      
      # out <- data.frame(table(sample))
      names(sample) <- c("sample", sprintf("%s_events", ifelse(j == 1, "increasing_severity", "decreasing_severity")),
                         sprintf("%s_variants", ifelse(j == 1, "increasing_severity", "decreasing_severity")))
      if(nrow(sample) > 0)
        sample <- aggregate(.~sample, sample, sum)
      
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

getCNVEvents <- function(this.cnv, set.sig, gsNFLD){
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
      
        del <- this.cnv[this.cnv$type == "deletion", c("sample", set[del])]
        dup <- this.cnv[this.cnv$type == "duplication", c("sample", set[dup])]
        
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
        if(nrow(sample) > 0)
          sample <- aggregate(.~sample, sample, sum)
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

getSNVEvents <- function(snv, set.sig){
  set.tmp <- set.sig
  # "3_0046_000"
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
        lof <- snv[snv$effect.tier == "lof", c("Sample", set[lof])]
        ms1 <- snv[snv$effect.tier == "tier1_ms", c("Sample", set[ms1])]
        ms2 <- snv[snv$effect.tier == "tier2_ms", c("Sample", set[ms2])]
        
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
        
        if(nrow(out) > 0){
          out[, "tmp2"] <- out[, 2]
          names(out) <- c("sample", sprintf("%s_events", ifelse(j == 1, "increasing_severity", "decreasing_severity")),
                          sprintf("%s_variants", ifelse(j == 1, "increasing_severity", "decreasing_severity")))
          out <- aggregate(.~sample, out, sum)
        }else{
          out <- data.frame("x1" = 1, "x2" = 2, "x3" =  3)[-1, ]
          names(out) <- c("sample", sprintf("%s_events", ifelse(j == 1, "increasing_severity", "decreasing_severity")),
                          sprintf("%s_variants", ifelse(j == 1, "increasing_severity", "decreasing_severity")))
        }
    
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

getSigFeatures <- function(dt, global, covariates, test, cutoff = 1, class = "GroupedDysmorphology"){
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
    ref <- sprintf("%s ~ %s + %s", class, paste(covariates, collapse=" + "), global)
    add <- sprintf("%s + %s", ref, feat)
    
    ref <- glm(ref, data=dt)
    add <- glm(add, data=dt)
    pvalue <- anova(ref, add, test = "Chisq")$`Pr(>Chi)`[2]
    
    # print(sprintf("%s's pvalue = %s", feat, pvalue))
    if(!is.na(pvalue) & pvalue <= cutoff){
      coeff <- c(coeff, add$coefficients[feat])
      p <- c(p, pvalue)
      r2 <- c(r2, NagelkerkeR2(add)$R2)
      
      
    }
  }
  
  if(!is.null(coeff)){
    coeff <- data.frame("set" = names(coeff), "coeff" = coeff, "r2" = r2,  "pvalue" = p, stringsAsFactors = F)
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

getScore <- function(dt, target, baseline, test, score.name, class = "GroupedDysmorphology"){
  if(length(test) > 0){
    ref <- gsub("\n", "", sprintf("%s ~ offset(I(%s * (Sex..CRV == 'M') + 
            %s * (Platform == 'Macrogen-Illumina HiSeqX') + %s * (Platform == 'Complete Genomics') + 
%s * (Platform == 'TCAG-Illumina HiSeqX') +
%s * (pc1) + 
%s * (pc2) + 
                                  %s * (pc3)))", class,
                                  baseline$coefficients["Sex..CRVM"], baseline$coefficients["PlatformMacrogen-Illumina HiSeqX"], baseline$coefficients["PlatformComplete Genomics"], 
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

load("../requiredData/gsNFLD.RData")

cds <- read.delim("../requiredData/hg19_refGene_28_04_17.cds.txt", stringsAsFactors = F, header = F)
cds.g <- GRanges(cds$V1, IRanges(cds$V2, cds$V3), "*")

all.nc.cnv <- read.delim("../dataNC/cnv.annotated.tsv", stringsAsFactors = F)
all.nc.cnv <- all.nc.cnv[!all.nc.cnv$ts_overlap, ]
all.nc.cnv$sample <- gsub("-", "_", all.nc.cnv$sample)
all.nc.cnv$sample <- gsub("A|_A", "", all.nc.cnv$sample)

all.cnv <- rbind(read.delim("../data/all.cnvs.may6.txt", stringsAsFactors = F),
                 read.delim("../data/all.small.cnvs.may6.txt", stringsAsFactors = F))

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
all.cnv$sample <- gsub("-", "_", all.cnv$sample)
all.cnv$sample <- gsub("A|_A", "", all.cnv$sample)

rare.lof <- read.delim("../data/SNVs/calls/rare/lof.variants.DEC062019.tsv", stringsAsFactors = F)
rare.ms1 <- read.delim("../data/SNVs/calls/rare/ms1.variants.DEC062019.tsv", stringsAsFactors = F)
rare.ms2 <- read.delim("../data/SNVs/calls/rare/ms2.variants.DEC062019.tsv", stringsAsFactors = F)
rare.lof$effect.tier <- "lof"
rare.ms1$effect.tier <- "tier1_ms"
rare.ms2$effect.tier <- "tier2_ms"
all.rare <- rbind(rare.lof, rbind(rare.ms1, rare.ms2))
nc.rare <- read.delim("../dataNC/snv.rare.annotated.tsv", stringsAsFactors = F)
nc.denovo <- read.delim("../dataNC/snv.denovo.annotated.tsv", stringsAsFactors = F)

eth <- read.delim("../NFLD_ethnicity_admixture_and_eigensoft/nfld.eth.tag.tsv", stringsAsFactors = F)

denovo.lof <- read.delim("../data/SNVs/calls/denovo/lof.variants.DEC062019.tsv", stringsAsFactors = F)
denovo.ms1 <- read.delim("../data/SNVs/calls/denovo/ms1.variants.DEC062019.tsv", stringsAsFactors = F)
denovo.ms2 <- read.delim("../data/SNVs/calls/denovo/ms2.variants.DEC062019.tsv", stringsAsFactors = F)
denovo.lof$effect.tier <- "lof"
denovo.ms1$effect.tier <- "tier1_ms"
denovo.ms2$effect.tier <- "tier2_ms"
all.denovo <- rbind(denovo.lof, rbind(denovo.ms1, denovo.ms2))

### read all 5 files
cnv.coding <- read.delim("../cnv.coding.matrix.tsv", stringsAsFactors = F)
cnv.coding$sample <- gsub("-", "_", cnv.coding$sample)
cnv.coding$sample <- gsub("A|_A", "", cnv.coding$sample)
names(cnv.coding)[1] <- "Sample"

scnv.coding <- read.delim("../cnv.smaller.coding.matrix.tsv", stringsAsFactors = F)
scnv.coding$sample <- gsub("-", "_", scnv.coding$sample)
scnv.coding$sample <- gsub("A|_A", "", scnv.coding$sample)
names(scnv.coding)[1] <- "Sample"

cnv.ncoding <- read.delim("../dataNC/nc.cnv.matrix.tsv", stringsAsFactors = F)
scnv.ncoding <- read.delim("../dataNC/nc.smaller.cnv.matrix.tsv", stringsAsFactors = F)

snv.coding.rare <- read.delim("rare.coding.snv.tsv", stringsAsFactors = F)
snv.coding.rare <- snv.coding.rare[!is.na(snv.coding.rare$GroupedDysmorphology), ]
#fix PCs for four samples in cnv.coding
cnv.coding <- merge(cnv.coding[, -c(87:89)], snv.coding.rare[, c(1, 129:131)], by = "Sample", all.x = T)
scnv.coding <- merge(scnv.coding[, -c(87:89)], snv.coding.rare[, c(1, 129:131)], by = "Sample", all.x = T)
cnv.ncoding <- merge(cnv.ncoding, snv.coding.rare[, c(1, 120, 125, 129:135)], by = "Sample", all.x = T)
scnv.ncoding <- merge(scnv.ncoding, snv.coding.rare[, c(1, 120, 125, 129:135)], by = "Sample", all.x = T)

snv.coding.denovo <- read.delim("denovo.coding.snv.tsv", stringsAsFactors = F)
snv.coding.denovo <- snv.coding.denovo[!is.na(snv.coding.denovo$GroupedDysmorphology), ]

snv.nc.rare <- read.delim("../scriptsNC/nc.snv.rare.matrix.tsv", stringsAsFactors = F)
snv.nc.rare <- snv.nc.rare[!is.na(snv.coding.rare$GroupedDysmorphology), ]

snv.nc.denovo <- read.delim("../scriptsNC/nc.snv.denovo.matrix.tsv", stringsAsFactors = F)
snv.nc.denovo <- snv.nc.denovo[!is.na(snv.nc.denovo$GroupedDysmorphology), ]

baseline.feat <- c("Sample", "Sex..CRV",  "Platform", "GroupedDysmorphology", "pc1", "pc2", "pc3")
samples <- na.omit(unique(rbind(rbind(rbind(cnv.ncoding[, baseline.feat], scnv.ncoding[,baseline.feat]), 
                                      scnv.coding[, baseline.feat]), rbind(cnv.coding[, baseline.feat], 
                                                                           rbind(snv.coding.rare[, baseline.feat],
                                                    rbind(snv.coding.denovo[, baseline.feat],
                                                          rbind(snv.nc.rare[, baseline.feat],
                                                                snv.nc.denovo[, baseline.feat])))))))
samples <- samples[samples$Platform != "", ]
writeLines(samples$Sample[samples$Sample %in% snv.nc.denovo$Sample], "samples.with.denovo.txt")
###read clisig variants
cli.cnv <- read.delim("CRV.sample.all.cnvs.tsv", stringsAsFactors = F)
cli.rare <- read.delim("CRV.sample.all.rare.snvs.tsv", stringsAsFactors = F)
cli.denovo <- read.delim("CRV.sample.all.denovo.snvs.tsv", stringsAsFactors = F)

# samples <- na.omit(unique(rbind(cnv.coding[, baseline.feat], snv.coding.rare[, baseline.feat])))
set.seed(1500)
score.out <- data.frame()
sig.set <- c()
coeff.all <- data.frame()
covariates <- c("Sex..CRV", "Platform", "pc1", "pc2", "pc3")
cli.out <- data.frame()

for(k in 1:30){
nfold <- 10
fold.ess <- split(sample(samples$Sample[samples$GroupedDysmorphology == 0]), 1:nfold)
fold.equ <- split(sample(samples$Sample[samples$GroupedDysmorphology == 1]), 1:nfold)
# fold.com <- split(sample(samples$Sample[samples$MultiClass == 2]), 1:nfold)

for(i in 1:length(fold.ess)){
  disc.set <- c(unlist(fold.ess[-i]), unlist(fold.equ[-i])) #, unlist(fold.com[-i]))
  targ.set <- c(unlist(fold.ess[i]), unlist(fold.equ[i])) #, unlist(fold.com[i]))
  # samples$MultiClass <- factor(samples$MultiClass)
  baseline <- glm("GroupedDysmorphology ~ Sex..CRV + Platform + pc1 + pc2 + pc3", data = samples[samples$Sample %in% disc.set, ])
  
  ## get significant geneset
  cnv1 <- getSigFeatures(cnv.coding[cnv.coding$Sample %in% disc.set, ], "gene_count_deletion", covariates, names(cnv.coding)[c(45:81)])
  cnv2 <- getSigFeatures(cnv.coding[cnv.coding$Sample %in% disc.set, ], "gene_count_duplication", covariates, names(cnv.coding)[c(5:41)])
  cnv.c.sig <- c(cnv1$set, cnv2$set)
  cnv.set <- rbind(cnv1$coeff, cnv2$coeff)
  if(length(cnv1$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(cnv1$coeff, "variant" = "coding_rare_CNVs", stringsAsFactors = F))
  if(length(cnv2$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(cnv2$coeff, "variant" = "coding_rare_CNVs", stringsAsFactors = F))
  
  scnv1 <- getSigFeatures(scnv.coding[scnv.coding$Sample %in% disc.set, ], "gene_count_deletion", covariates, names(scnv.coding)[5:41])
  scnv2 <- getSigFeatures(scnv.coding[scnv.coding$Sample %in% disc.set, ], "gene_count_duplication", covariates, names(scnv.coding)[45:81])
  scnv.c.sig <- c(scnv1$set, scnv2$set)
  scnv.set <- rbind(scnv1$coeff, scnv2$coeff)
  if(length(scnv1$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(scnv1$coeff, "variant" = "small_coding_rare_CNVs", stringsAsFactors = F))
  if(length(scnv2$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(scnv2$coeff, "variant" = "small_coding_rare_CNVs", stringsAsFactors = F))
  
  ###nc cnv
  cnv1.nc <- getSigFeatures(cnv.ncoding[cnv.ncoding$Sample %in% disc.set, ], "TotalRare_deletion", covariates, names(cnv.ncoding)[2:10])
  cnv2.nc <- getSigFeatures(cnv.ncoding[cnv.ncoding$Sample %in% disc.set, ], "TotalRare_duplication", covariates, names(cnv.ncoding)[12:20])
  cnv.nc.sig <- c(cnv1.nc$set, cnv2.nc$set)
  cnv.nc.set <- rbind(cnv1.nc$coeff, cnv2.nc$coeff)
  if(length(cnv1.nc$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(cnv1.nc$coeff, "variant" = "noncoding_rare_CNVs", stringsAsFactors = F))
  if(length(cnv2.nc$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(cnv2.nc$coeff, "variant" = "noncoding_rare_CNVs", stringsAsFactors = F))
  
  scnv1.nc <- getSigFeatures(scnv.ncoding[scnv.ncoding$Sample %in% disc.set, ], "TotalRare_deletion", covariates, names(scnv.ncoding)[2:10])
  scnv2.nc <- getSigFeatures(scnv.ncoding[scnv.ncoding$Sample %in% disc.set, ], "TotalRare_duplication", covariates, names(scnv.ncoding)[12:20])
  scnv.nc.sig <- c(scnv1.nc$set, scnv2.nc$set)
  scnv.nc.set <- rbind(scnv1.nc$coeff, scnv2.nc$coeff)
  if(length(scnv1.nc$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(scnv1.nc$coeff, "variant" = "small_noncoding_rare_CNVs", stringsAsFactors = F))
  if(length(scnv2.nc$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(scnv2.nc$coeff, "variant" = "small_noncoding_rare_CNVs", stringsAsFactors = F))
  
  ###snvs
  snv1 <- getSigFeatures(snv.coding.rare[snv.coding.rare$Sample %in% disc.set, ], "Total_lof", covariates, names(snv.coding.rare)[4:40])
  snv2 <- getSigFeatures(snv.coding.rare[snv.coding.rare$Sample %in% disc.set, ], "Total_tier1_ms", covariates, names(snv.coding.rare)[42:78])
  snv3 <- getSigFeatures(snv.coding.rare[snv.coding.rare$Sample %in% disc.set, ], "Total_tier2_ms", covariates, names(snv.coding.rare)[80:116])
  snv.c.rare.sig <- c(snv1$set, snv2$set, snv3$set)
  rare.set <- rbind(snv1$coeff, rbind(snv2$coeff, snv3$coeff))
  
  if(length(snv1$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(snv1$coeff, "variant" = "coding_rare_SNVs", stringsAsFactors = F))
  if(length(snv2$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(snv2$coeff, "variant" = "coding_rare_SNVs", stringsAsFactors = F))
  if(length(snv3$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(snv3$coeff, "variant" = "coding_rare_SNVs", stringsAsFactors = F))
  
  snv1 <- getSigFeatures(snv.coding.denovo[snv.coding.denovo$Sample %in% disc.set, ], "Total_lof", covariates, names(snv.coding.denovo)[4:40])
  snv2 <- getSigFeatures(snv.coding.denovo[snv.coding.denovo$Sample %in% disc.set, ], "Total_tier1_ms", covariates, names(snv.coding.denovo)[42:78])
  snv3 <- getSigFeatures(snv.coding.denovo[snv.coding.denovo$Sample %in% disc.set, ], "Total_tier2_ms", covariates, names(snv.coding.denovo)[80:116])
  snv.c.denovo.sig <- c(snv1$set, snv2$set, snv3$set)
  denovo.set <- rbind(snv1$coeff, rbind(snv2$coeff, snv3$coeff))
  if(length(snv1$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(snv1$coeff, "variant" = "coding_denovo_SNVs", stringsAsFactors = F))
  if(length(snv2$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(snv2$coeff, "variant" = "coding_denovo_SNVs", stringsAsFactors = F))
  if(length(snv3$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(snv3$coeff, "variant" = "coding_denovo_SNVs", stringsAsFactors = F))
  
  snv.nc.rare.sig <- getSigFeatures(snv.nc.rare[snv.nc.rare$Sample %in% disc.set, ], "TotalRare", covariates, names(snv.nc.rare)[2:26])
  nc.rare.set <- snv.nc.rare.sig$coeff
  if(length(snv.nc.rare.sig$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(snv.nc.rare.sig$coeff, "variant" = "noncoding_rare_SNVs", stringsAsFactors = F))
  snv.nc.rare.sig <- snv.nc.rare.sig$set
  
  snv.nc.denovo.sig <- getSigFeatures(snv.nc.denovo[snv.nc.denovo$Sample %in% disc.set, ], "TotalRare", covariates, names(snv.nc.denovo)[2:26])
  nc.denovo.set <- snv.nc.denovo.sig$coeff
  if(length(snv.nc.denovo.sig$set) > 0)
    coeff.all <- rbind(coeff.all, data.frame(snv.nc.denovo.sig$coeff, "variant" = "noncoding_denovo_SNVs", stringsAsFactors = F))
  snv.nc.denovo.sig <- snv.nc.denovo.sig$set
  
  
  for(p in c(1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001)){
  
    all.event <- data.frame()
    if(sum(cnv.set$pvalue < p) > 0){
      cnv.event <- getCNVEvents(all.cnv[all.cnv$size > 10000 & all.cnv$sample %in% targ.set, ], cnv.set[cnv.set$pvalue < p, ], gsNFLD)
      names(cnv.event)[-1] <- paste(names(cnv.event)[-1], "rare_CNVs", sep="_")
      all.event <- cnv.event
    }
    if(sum(scnv.set$pvalue < p) > 0){
      scnv.event <- getCNVEvents(all.cnv[all.cnv$size <= 10000 & all.cnv$sample %in% targ.set, ], scnv.set[scnv.set$pvalue < p, ], gsNFLD)
      names(scnv.event)[-1] <- paste(names(scnv.event)[-1], "small_rare_CNVs", sep="_")
      if(nrow(all.event) == 0)
        all.event <- scnv.event
      else
        all.event <- merge(all.event, scnv.event, by = "sample", all = T)
    }
    if(sum(cnv.nc.set$pvalue < p) > 0){
      cnv.nc.event <- getCNVEvents(all.nc.cnv[all.nc.cnv$size > 10000 & all.nc.cnv$sample %in% targ.set, ], cnv.nc.set[cnv.nc.set$pvalue < p, ], gsNFLD)
      names(cnv.nc.event)[-1] <- paste(names(cnv.nc.event)[-1], "rare_noncoding_CNVs", sep="_")
      
      if(nrow(all.event) == 0)
        all.event <- cnv.nc.event
      else
        all.event <- merge(all.event, cnv.nc.event, by = "sample", all = T)
    }
    if(sum(scnv.nc.set$pvalue < p) > 0){
      scnv.nc.event <- getCNVEvents(all.nc.cnv[all.nc.cnv$size <= 10000 & all.nc.cnv$sample %in% targ.set, ], scnv.nc.set[scnv.nc.set$pvalue < p, ], gsNFLD)
      names(scnv.nc.event)[-1] <- paste(names(scnv.nc.event)[-1], "small_rare_noncoding_CNVs", sep="_")
      
      if(nrow(all.event) == 0)
        all.event <- scnv.nc.event
      else
        all.event <- merge(all.event, scnv.nc.event, by = "sample", all = T)
    }
    
    if(sum(rare.set$pvalue < p) > 0){
      rare.event <- getSNVEvents(all.rare[all.rare$Sample %in% targ.set, ], rare.set[rare.set$pvalue < p, ])
      names(rare.event)[-1] <- paste(names(rare.event)[-1], "rare_SNVs", sep="_")
      
      if(nrow(all.event) == 0)
        all.event <- rare.event
      else
        all.event <- merge(all.event, rare.event, by = "sample", all = T)
    }
    if(sum(denovo.set$pvalue < p) > 0){
      denovo.event <- getSNVEvents(all.denovo[all.denovo$Sample %in% targ.set, ], denovo.set[denovo.set$pvalue < p, ])
      names(denovo.event)[-1] <- paste(names(denovo.event)[-1], "denovo_SNVs", sep="_")
      
      if(nrow(all.event) == 0)
        all.event <- denovo.event
      else
        all.event <- merge(all.event, denovo.event, by = "sample", all = T)
    }
    if(sum(nc.rare.set$pvalue < p) > 0){
      nc.rare.event <- getNCEvents(nc.rare[nc.rare$Sample %in% targ.set, ], nc.rare.set[nc.rare.set$pvalue < p, ])
      names(nc.rare.event)[-1] <- paste(names(nc.rare.event)[-1], "rare_noncoding_SNVs", sep="_")
      
      if(nrow(all.event) == 0)
        all.event <- nc.rare.event
      else
        all.event <- merge(all.event, nc.rare.event, by = "sample", all = T)
    }
    if(sum(nc.denovo.set$pvalue < p) > 0){
      nc.denovo.event <- getNCEvents(nc.denovo[nc.denovo$Sample %in% targ.set, ], nc.denovo.set[nc.denovo.set$pvalue < p, ])
      names(nc.denovo.event)[-1] <- paste(names(nc.denovo.event)[-1], "denovo_noncoding_SNVs", sep="_")
      
      if(nrow(all.event) == 0)
        all.event <- nc.denovo.event
      else
        all.event <- merge(all.event, nc.denovo.event, by = "sample", all = T)
    }
    
    
    if(nrow(all.event) != 0){
      all.event[is.na(all.event)] <- 0
      
      if(length(grep("increasing_severity_events", names(all.event))) > 1)
        all.event$all_increasing_severity_events <- rowSums(all.event[, grep("increasing_severity_events", names(all.event))], na.rm = T)
      if(length(grep("increasing_severity_events", names(all.event))) == 1)
        all.event$all_increasing_severity_events <- all.event[, grep("increasing_severity_events", names(all.event))]
      if(length(grep("increasing_severity_events", names(all.event))) == 0)
        all.event$all_increasing_severity_events <- 0
      
      if(length(grep("decreasing_severity_events", names(all.event))) > 1)
        all.event$all_decreasing_severity_events <- rowSums(all.event[, grep("decreasing_severity_events", names(all.event))], na.rm = T)
      if(length(grep("decreasing_severity_events", names(all.event))) == 1)
        all.event$all_decreasing_severity_events <- all.event[, grep("decreasing_severity_events", names(all.event))]
      if(length(grep("decreasing_severity_events", names(all.event))) == 0)
        all.event$all_decreasing_severity_events <- 0
      
      
      if(length(grep("increasing_severity_variants", names(all.event))) > 1)
        all.event$all_increasing_severity_variants <- rowSums(all.event[, grep("increasing_severity_variants", names(all.event))], na.rm = T)
      if(length(grep("increasing_severity_variants", names(all.event))) == 1)
        all.event$all_increasing_severity_variants <- all.event[, grep("increasing_severity_variants", names(all.event))]
      if(length(grep("increasing_severity_variants", names(all.event))) == 0)
        all.event$all_increasing_severity_variants <- 0
      
      if(length(grep("decreasing_severity_variants", names(all.event))) > 1)
        all.event$all_decreasing_severity_variants <- rowSums(all.event[, grep("decreasing_severity_variants", names(all.event))], na.rm = T)
      if(length(grep("decreasing_severity_variants", names(all.event))) == 1)
        all.event$all_decreasing_severity_variants <- all.event[, grep("decreasing_severity_variants", names(all.event))]
      if(length(grep("decreasing_severity_variants", names(all.event))) == 0)
        all.event$all_decreasing_severity_variants <- 0
      
      
      sig.set <- c(sig.set, 
                   paste("CNVs",cnv.c.sig, sep=":"),
                   paste("small CNVs", scnv.c.sig, sep=":"), 
                   paste("noncoding CNVs",cnv.nc.sig, sep=":"),
                   paste("noncoding small CNVs", scnv.nc.sig, sep=":"), 
                   paste("rare SNVs", snv.c.rare.sig, sep=":"),
                   paste("de novo SNVs", snv.c.denovo.sig, sep=":"), 
                   paste("rare noncoding SNVs", snv.nc.rare.sig, sep=":"), 
                   paste("de novo noncoding SNVs", snv.nc.denovo.sig, sep=":"))
      if(sum(cnv.set$pvalue < p) > 0){
        cnv.c.cli <- getScore(cnv.coding[cnv.coding$Sample %in% disc.set, ],
                            cli.cnv[cli.cnv$Sample %in% targ.set, ], baseline, cnv.c.sig[cnv.set$pvalue < p], "coding_rare_CNVs")
      }else{
        cnv.c.cli <- data.frame("Sample" = NA, "crv_score" = NA)[-1, ]
      }
      if(sum(rare.set$pvalue < p & snv.c.rare.sig %in% names(cli.rare)) > 0){
        snv.rare.cli <- getScore(snv.coding.rare[snv.coding.rare$Sample %in% disc.set, ],
                                 cli.rare[cli.rare$Sample %in% targ.set, ], baseline, 
                                 snv.c.rare.sig[rare.set$pvalue < p & snv.c.rare.sig %in% names(cli.rare)], "coding_rare_SNVs")
      }else{
        snv.rare.cli <- data.frame("Sample" = NA, "crv_score" = NA)[-1, ]
      }
      if(sum(denovo.set$pvalue < p & snv.c.denovo.sig %in% names(cli.denovo)) > 0){
        snv.denovo.cli <- getScore(snv.coding.denovo[snv.coding.denovo$Sample %in% disc.set, ],
                                       cli.denovo[cli.denovo$Sample %in% targ.set, ], baseline, 
                                       snv.c.denovo.sig[denovo.set$pvalue < p & snv.c.denovo.sig %in% names(cli.denovo)], "coding_denovo_SNVs")
      }else{
        snv.denovo.cli <- data.frame("Sample" = NA, "crv_score" = NA)[-1, ]
      }
      names(cnv.c.cli)[2] <- names(snv.rare.cli)[2] <- names(snv.denovo.cli)[2] <- "crv_score"
      cli.tmp <- rbind(cnv.c.cli, rbind(snv.rare.cli, snv.denovo.cli))
      if(nrow(cli.tmp) > 0){
        cli.tmp$p <- p
        cli.out <- rbind(cli.out, cli.tmp)
      }
      
      cnv.c.score <- getScore(cnv.coding[cnv.coding$Sample %in% disc.set, ],
                              cnv.coding[cnv.coding$Sample %in% targ.set, ], baseline, cnv.c.sig[cnv.set$pvalue < p], "coding_rare_CNVs")
      
      scnv.c.score <- getScore(scnv.coding[scnv.coding$Sample %in% disc.set, ],
                              scnv.coding[scnv.coding$Sample %in% targ.set, ], baseline, scnv.c.sig[scnv.set$pvalue < p], "small_coding_rare_CNVs")
      
      cnv.nc.score <- getScore(cnv.ncoding[cnv.ncoding$Sample %in% disc.set, ],
                              cnv.ncoding[cnv.ncoding$Sample %in% targ.set, ], baseline, cnv.nc.sig[cnv.nc.set$pvalue < p], "noncoding_rare_CNVs")
      
      scnv.nc.score <- getScore(scnv.ncoding[scnv.ncoding$Sample %in% disc.set, ],
                               scnv.ncoding[scnv.ncoding$Sample %in% targ.set, ], baseline, scnv.nc.sig[scnv.nc.set$pvalue < p], "small_noncoding_rare_CNVs")
      
      snv.c.rare.score <- getScore(snv.coding.rare[snv.coding.rare$Sample %in% disc.set, ],
                                   snv.coding.rare[snv.coding.rare$Sample %in% targ.set, ], baseline, snv.c.rare.sig[rare.set$pvalue < p], "coding_rare_SNVs")
      
      snv.c.denovo.score <- getScore(snv.coding.denovo[snv.coding.denovo$Sample %in% disc.set, ],
                              snv.coding.denovo[snv.coding.denovo$Sample %in% targ.set, ], baseline, snv.c.denovo.sig[denovo.set$pvalue < p], "coding_denovo_SNVs")
      
      snv.nc.rare.score <- getScore(snv.nc.rare[snv.nc.rare$Sample %in% disc.set, ],
                              snv.nc.rare[snv.nc.rare$Sample %in% targ.set, ], baseline, snv.nc.rare.sig[nc.rare.set$pvalue < p], "noncoding_rare_SNVs")
    
      snv.nc.denovo.score <- getScore(snv.nc.denovo[snv.nc.denovo$Sample %in% disc.set, ],
                              snv.nc.denovo[snv.nc.denovo$Sample %in% targ.set, ], baseline, snv.nc.denovo.sig[nc.denovo.set$pvalue < p], "noncoding_denovo_SNVs")
      
      final.score <- merge(merge(merge(merge(merge(merge(merge(cnv.nc.score, scnv.nc.score, by = "Sample", all = T),
                                                   cnv.c.score, by = "Sample", all = T), 
                                                   scnv.c.score, by = "Sample", all = T), 
                                             snv.c.rare.score, by = "Sample", all = T), snv.c.denovo.score, by = "Sample", all = T), 
                                 snv.nc.rare.score, by = "Sample", all = T), 
                           snv.nc.denovo.score, by = "Sample", all = T)
      
      final.score$total_risk <- rowSums(final.score[, -1], na.rm = T)
      all.event$sample <- gsub("-", "_", all.event$sample)
      final.score <- merge(final.score, all.event[, c("sample", "all_increasing_severity_events", "all_decreasing_severity_events",
                                                      "all_increasing_severity_variants", "all_decreasing_severity_variants")], by.x = "Sample", by.y = "sample", all.x = T)
      final.score$pvalue <- p
      
      # dt.plot <- merge(final.score, samples, by = "Sample", all.x = T)
      # dt.plot$DysmorphologyClassificaiton <- ifelse(dt.plot$MultiClass == 0, "Essential", "Complex")
      # dt.plot$DysmorphologyClassificaiton[dt.plot$MultiClass == 1] <- "Equivocal"
      # 
      # ggplot(dt.plot, aes(x = DysmorphologyClassificaiton, y = total_risk, fill = DysmorphologyClassificaiton)) +
      #   geom_boxplot() + theme_classic() + theme(legend.position = "none", axis.title.x = element_blank()) + coord_cartesian(ylim = c(-0.5, 2.5))
      #   
      score.out <- rbind(score.out, final.score)
    }### all.event != 0
  }### p-value cutoff loop
}### 10 folds
}### 30 reps

write.table(cli.out, "cli.score.from.30x.10fold.2groups.new.tsv", sep="\t", row.names=F, quote=F, col.names=T)
write.table(coeff.all, "coeff.from.30x.10fold.2groups.new.tsv", sep="\t", row.names=F, quote=F, col.names=T)
write.table(score.out, "score.from.30x.10fold.2groups.new.tsv", sep="\t", row.names=F, quote = F, col.names=T)
# write.table(sig.set, "significant.set.from.30x.10fold.tsv", sep="\t", row.names=F, quote=F, col.names=T)
write.table(samples, "sample.with.score.2groups.new.tsv", sep="\t", row.names=F, quote=F, col.names=T)



########################################
######### test different cutoffs #######

samples <- read.delim("sample.with.score.2groups.new.tsv", stringsAsFactors = F)
denovo <- readLines("samples.with.denovo.txt")

samples <- samples[samples$Sample %in% denovo, ]

out.perf <- data.frame()
for(p in c(1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001)){
  score.out <- data.table::fread("score.from.30x.10fold.2groups.new.tsv", data.table = F)
  score.out <- score.out[score.out$pvalue == p, ]
  coeff.df <- read.delim("coeff.from.30x.10fold.2groups.new.tsv", stringsAsFactors = F)
  coeff.df <- coeff.df[coeff.df$pvalue < p, ]
  # coeff.df <- coeff.df[grep("duplication", coeff.df$set), ]
  agg.coeff.df <- merge(aggregate(.~ set + variant, coeff.df[-3], sum), as.data.frame(table(coeff.df$set)), by.x = "set", by.y = "Var1", all = T)
  agg.coeff.df$coeff <- agg.coeff.df$coeff/20
  agg.coeff.df
  tmp.out <- data.frame()
  for(s in unique(score.out$Sample)){
    tmp <- score.out[score.out$Sample == s, ]
    # message(nrow(tmp))
    tmp <- colSums(tmp[, -1], na.rm = T)/nrow(tmp)
    tmp.out <- rbind(tmp.out, data.frame("Sample" = s, t(tmp)))
  }
  
  score.out <- tmp.out
  score.out <- merge(score.out, samples, by = "Sample", all = F)
  # score.out$MultiClass <- factor(score.out$MultiClass)
  # score.out$total_risk <- score.out$coding_rare_CNVs + score.out$small_coding_rare_CNVs + 
  #   score.out$coding_rare_SNVs + score.out$coding_denovo_SNVs
  lm <- glm(GroupedDysmorphology ~ total_risk, data = score.out)

  p.value <- paste0(sprintf("p=%s", formatC(signif(summary(lm)$coefficients[2, "Pr(>|t|)"]), 2), 
                                            format = "e", digits = 1))
  out.perf <- rbind(out.perf, data.frame("cutoff" = p,
                                         "NagelkerkeR2" = signif(NagelkerkeR2(lm)$R2, digits = 3),
                                         "p" = p.value))
}

out.perf$cutoff <- factor(out.perf$cutoff, levels = rev(out.perf$cutoff))
ggplot(out.perf, aes(x = paste0("p<",cutoff), y = NagelkerkeR2, label = p)) + 
  geom_bar(stat = "identity", color = "black", fill = "red") + xlab("p cut-off") + ylim(0, 0.05) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_text(angle = 90, nudge_y = 0.008)
ggsave("performance.2subtypes.new.pdf", width = 4, height = 5)


########################################
### get cli samples
# new.cli.snv <- data.frame(readxl::read_excel("../data/Draft 6-Supplementary Tables.xlsx", sheet = 10, skip = 1), stringsAsFactors = F)
# new.cli.snv <- new.cli.snv[!new.cli.snv$Variant.Classification %in% c("ASD candidate-Category 2", "ASD candidate-Category 1"), ]
# 
# new.cli.cnv <- data.frame(readxl::read_excel("../data/Draft 6-Supplementary Tables.xlsx", sheet = 11, skip = 1), stringsAsFactors = F)
# new.cli.cnv <- new.cli.cnv[!new.cli.cnv$Variant.classification %in% c("ASD candidate-Category 2", "ASD candidate-Category 1"), ]
# # new.cli.cnv <- new.cli.cnv[1:40, ]
# new.cli.cnv <- new.cli.cnv[new.cli.cnv$Type %in% c("Deletion", "Duplication", "Deletion/unbalanced translocation",
#                                                    "Duplication (tandem)"), ]
# cli.samples <- unique(c(gsub("-", "_", new.cli.snv$Case),gsub("-", "_", new.cli.cnv$Sample)))
# writeLines(cli.samples, "sample.with.crv.txt")
########### original plots #############
########### p < 0.1 ###########
samples <- read.delim("sample.with.score.2groups.new.tsv", stringsAsFactors = F)
denovo <- readLines("samples.with.denovo.txt")
samples <- samples[samples$Sample %in% denovo, ]
score.out <- read.delim("score.from.30x.10fold.2groups.new.tsv", stringsAsFactors = F)

score.out <- score.out[score.out$pvalue == 0.1, ]
coeff.df <- read.delim("coeff.from.30x.10fold.2groups.new.tsv", stringsAsFactors = F)
coeff.df <- coeff.df[coeff.df$pvalue < 0.1, ]

cli.out <- read.delim("cli.score.from.30x.10fold.2groups.new.tsv", stringsAsFactors = F)
cli.out <- cli.out[cli.out$p == 0.1, ]
cli.out <- aggregate(crv_score ~ Sample, cli.out, sum)
cli.out$crv_score <- cli.out$crv_score/30
cli.out <- cli.out[cli.out$Sample %in% readLines("sample.with.crv.txt"), ]

agg.coeff.df <- merge(aggregate(.~ set + variant, coeff.df[-4], sum), as.data.frame(table(coeff.df$set)), by.x = "set", by.y = "Var1", all = T)
agg.coeff.df$coeff <- agg.coeff.df$coeff/300
mean.coeff <- aggregate(coeff ~ variant, agg.coeff.df, mean)
sd.coeff <- aggregate(coeff ~ variant, agg.coeff.df, sd)
names(mean.coeff)[2] <- "mean coeff"
names(sd.coeff)[2] <- "SD coeff"
agg.coeff.df <- merge(mean.coeff, sd.coeff, by = "variant")
agg.coeff.df$variant <- gsub("_", " ", agg.coeff.df$variant)
agg.coeff.df$variant <- factor(agg.coeff.df$variant, levels = agg.coeff.df$variant[order(agg.coeff.df$`mean coeff`, decreasing = T)])
coeff.plot <- ggplot(agg.coeff.df, aes(x = variant, y = `mean coeff`)) + geom_bar(stat = "identity", fill = "skyblue", color = "black", width = .5) +
  geom_point() +
  theme_classic() + geom_errorbar(aes(ymin = `mean coeff`-`SD coeff`, ymax = `mean coeff`+`SD coeff`), width = .2) +
  geom_hline(yintercept = 0, lty = 1) + xlab("") + ylab("Beta coefficient") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

tmp.out <- data.frame()
for(s in unique(score.out$Sample)){
  tmp <- score.out[score.out$Sample == s, ]
  tmp <- colSums(tmp[, -1], na.rm = T)/nrow(tmp)
  
  tmp.out <- rbind(tmp.out, data.frame("Sample" = s, t(tmp)))
}

score.out <- tmp.out
score.out <- merge(score.out, samples, by = "Sample", all = F)
score.out$group <- ifelse(score.out$GroupedDysmorphology == 0, "Nondysmorphic", "Dysmorphic")
score.out <- merge(score.out, cli.out, by = "Sample", all.x = T)
write.table(score.out[!is.na(score.out$crv_score), ], "crv.score.sample.tsv", sep="\t", row.names=F, quote=F, col.names=T)




# score.out$DysmorphologyClassificaiton[score.out$MultiClass == 1] <- "Equivocal"
# score.out$DysmorphologyClassificaiton <- factor(score.out$DysmorphologyClassificaiton, levels = c("Essential", "Equivocal", "Complex"))
suffix <- ""

# cli.snv <- data.frame(readxl::read_excel("../data/clinicallysignificant_var_forBE.xlsx", sheet = 1), stringsAsFactors = F)
# cli.cnv <- data.frame(readxl::read_excel("../data/clinicallysignificant_var_forBE.xlsx", sheet = 2), stringsAsFactors = F)
# score.out <- score.out[!score.out$Sample %in% gsub("-", "_", c(as.character(cli.snv$Case), as.character(cli.cnv$Sample))), ]
# suffix <- ".noCRV"
# score.out$total_risk <- rowSums(score.out[, 2:4], na.rm = T)
#write.table(score.out[order(score.out$total_risk, decreasing = T), ], "cross.validate.score.tsv", sep="\t", row.names=F, quote=F, col.names=T)

### calculate percentile
perc.rank <- function(x) trunc(rank(x))/length(x)

### calculate percentile
score.out$coding_CNVs_perc <- perc.rank(score.out$coding_rare_CNVs)
score.out$small_coding_CNVs_perc <- perc.rank(score.out$small_coding_rare_CNVs)
score.out$noncoding_CNVs_perc <- perc.rank(score.out$noncoding_rare_CNVs)
score.out$small_noncoding_CNVs_perc <- perc.rank(score.out$small_noncoding_rare_CNVs)
score.out$coding_rare_SNVs_perc <- perc.rank(score.out$coding_rare_SNVs)
score.out$coding_denovo_SNVs_perc <- perc.rank(score.out$coding_denovo_SNVs)
score.out$noncoding_rare_SNVs_perc <- perc.rank(score.out$noncoding_rare_SNVs)
score.out$noncoding_denovo_SNVs_perc <- perc.rank(score.out$noncoding_denovo_SNVs)

score.out$total_risk_perc <- perc.rank(score.out$total_risk)
####################
### trend
score.out$subtype <- ifelse(score.out$GroupedDysmorphology == 1, "Dysmorphic", "Nondysmorphic")
write.table(score.out, "nfld.original.grs.tsv", sep="\t", row.names=F, quote=F,col.names=T)

# p1 <- ggplot(score.out, aes(x = all_increasing_severity_events, y = total_risk)) + 
#   geom_smooth(data = score.out[score.out$subtype == "Dysmorphic",], method = "lm", fullrange=T, color = "#F28602", fill = "#F28602", alpha = .25, lty = 2, lwd = .4) + 
#   geom_smooth(data = score.out[score.out$subtype == "Nondysmorphic",], method = "lm", fullrange=T, color = "#009eed", fill = "#009eed", alpha = .25, lty = 2, lwd = .4) + 
#   geom_point(shape = 1, aes(color = subtype, fill = subtype), show.legend = F, size = 1.5) + theme_bw() + 
#   xlab("#genes impacted in positively correlated sets") + ylab("total score")
# p2 <- ggplot(score.out, aes(x = all_decreasing_severity_events, y = total_risk)) + 
#   geom_smooth(data = score.out[score.out$subtype == "Dysmorphic",], method = "lm", fullrange=T, color = "#F28602", fill = "#F28602", alpha = .25, lty = 2, lwd = .4) + 
#   geom_smooth(data = score.out[score.out$subtype == "Nondysmorphic",], method = "lm", fullrange=T, color = "#009eed", fill = "#009eed", alpha = .25, lty = 2, lwd = .4) + 
#   geom_point(shape = 1, aes(color = subtype, fill = subtype), show.legend = F, size = 1.5) + theme_bw() + 
#   xlab("#genes impacted in negatively correlated sets") + ylab("total score")
# plot_grid(p1, p2, nrow =2, rel_heights = c(0.9, 1))
# ggsave(sprintf("score.genes%s.pdf", suffix), width = 7, height = 10)

### trend
# p1 <- ggplot(score.out, aes(x = all_increasing_severity_variants, y = total_risk)) + 
#   geom_smooth(data = score.out[score.out$subtype == "Dysmorphic",], method = "lm", fullrange=T, color = "#F28602", fill = "#F28602", alpha = .25, lty = 2, lwd = .4) + 
#   geom_smooth(data = score.out[score.out$subtype == "Nondysmorphic",], method = "lm", fullrange=T, color = "#009eed", fill = "#009eed", alpha = .25, lty = 2, lwd = .4) + 
#   geom_point(aes(fill = subtype), show.legend = F, size = 1.5, shape = 21, color = "black") + theme_bw() + 
#   xlab("#variants impacted in positively correlated sets") + ylab("total score")
# # p2 <- ggplot(score.out, aes(x = all_decreasing_severity_variants, y = total_risk)) + 
# #   geom_smooth(data = score.out[score.out$subtype == "Dysmorphic",], method = "lm", fullrange=T, color = "#F28602", fill = "#F28602", alpha = .25, lty = 2, lwd = .4) + 
# #   geom_smooth(data = score.out[score.out$subtype == "Nondysmorphic",], method = "lm", fullrange=T, color = "#009eed", fill = "#009eed", alpha = .25, lty = 2, lwd = .4) + 
# #   geom_point(shape = 1, aes(color = subtype, fill = subtype), show.legend = F, size = 1.5) + theme_bw() + 
# #   xlab("#variants impacted in negatively correlated sets") + ylab("total score")
# plot_grid(p1, coeff.plot, nrow = 1)
# ggsave("num.variants_coeff.mean.pdf", width = 10, height = 5)


########### plot ###########
# p1.pvalue <- paste0(sprintf("Co-Es's p=%s", formatC(wilcox.test(score.out$coding_rare_CNVs[score.out$DysmorphologyClassificaiton == "Complex"],
#                                  score.out$coding_CNVs[score.out$DysmorphologyClassificaiton == "Essential"], "greater")$p.value, 
#                      format = "e", digits = 1)), "\n",
#                    sprintf("Co-Eq's p=%s", formatC(wilcox.test(score.out$coding_CNVs[score.out$DysmorphologyClassificaiton == "Complex"],
#                                                                score.out$coding_CNVs[score.out$DysmorphologyClassificaiton == "Equivocal"], "greater")$p.value, 
#                                                    format = "e", digits = 1)), "\n",
#                    sprintf("Eq-Es's p=%s", formatC(wilcox.test(score.out$coding_CNVs[score.out$DysmorphologyClassificaiton == "Equivocal"],
#                                                                score.out$coding_CNVs[score.out$DysmorphologyClassificaiton == "Essential"], "greater")$p.value, 
#                                                    format = "e", digits = 1)))
# 
# 
# p2.pvalue <-  paste0(sprintf("Co-Es's p=%s", formatC(wilcox.test(score.out$coding_rare_SNVs[score.out$DysmorphologyClassificaiton == "Complex"],
#                                                                  score.out$coding_rare_SNVs[score.out$DysmorphologyClassificaiton == "Essential"], "greater")$p.value, 
#                                                      format = "e", digits = 1)), "\n",
#                      sprintf("Co-Eq's p=%s", formatC(wilcox.test(score.out$coding_rare_SNVs[score.out$DysmorphologyClassificaiton == "Complex"],
#                                                                  score.out$coding_rare_SNVs[score.out$DysmorphologyClassificaiton == "Equivocal"], "greater")$p.value, 
#                                                      format = "e", digits = 1)), "\n",
#                      sprintf("Eq-Es's p=%s", formatC(wilcox.test(score.out$coding_rare_SNVs[score.out$DysmorphologyClassificaiton == "Equivocal"],
#                                                                  score.out$coding_rare_SNVs[score.out$DysmorphologyClassificaiton == "Essential"], "greater")$p.value, 
#                                                      format = "e", digits = 1)))
# 
# p3.pvalue <-  paste0(sprintf("Co-Es's p=%s", formatC(wilcox.test(score.out$coding_denovo_SNVs[score.out$DysmorphologyClassificaiton == "Complex"],
#                                                                  score.out$coding_denovo_SNVs[score.out$DysmorphologyClassificaiton == "Essential"], "greater")$p.value, 
#                                                      format = "e", digits = 1)), "\n",
#                      sprintf("Co-Eq's p=%s", formatC(wilcox.test(score.out$coding_denovo_SNVs[score.out$DysmorphologyClassificaiton == "Complex"],
#                                                                  score.out$coding_denovo_SNVs[score.out$DysmorphologyClassificaiton == "Equivocal"], "greater")$p.value, 
#                                                      format = "e", digits = 1)), "\n",
#                      sprintf("Eq-Es's p=%s", formatC(wilcox.test(score.out$coding_denovo_SNVs[score.out$DysmorphologyClassificaiton == "Equivocal"],
#                                                                  score.out$coding_denovo_SNVs[score.out$DysmorphologyClassificaiton == "Essential"], "greater")$p.value, 
#                                                      format = "e", digits = 1)))
# 
# p4.pvalue <-  paste0(sprintf("Co-Es's p=%s", formatC(wilcox.test(score.out$noncoding_rare_SNVs[score.out$DysmorphologyClassificaiton == "Complex"],
#                                                                  score.out$noncoding_rare_SNVs[score.out$DysmorphologyClassificaiton == "Essential"], "greater")$p.value, 
#                                                      format = "e", digits = 1)), "\n",
#                      sprintf("Co-Eq's p=%s", formatC(wilcox.test(score.out$noncoding_rare_SNVs[score.out$DysmorphologyClassificaiton == "Complex"],
#                                                                  score.out$noncoding_rare_SNVs[score.out$DysmorphologyClassificaiton == "Equivocal"], "greater")$p.value, 
#                                                      format = "e", digits = 1)), "\n",
#                      sprintf("Eq-Es's p=%s", formatC(wilcox.test(score.out$noncoding_rare_SNVs[score.out$DysmorphologyClassificaiton == "Equivocal"],
#                                                                  score.out$noncoding_rare_SNVs[score.out$DysmorphologyClassificaiton == "Essential"], "greater")$p.value, 
#                                                      format = "e", digits = 1)))
# 
# p5.pvalue <- paste0(sprintf("Co-Es's p=%s", formatC(wilcox.test(score.out$noncoding_denovo_SNVs[score.out$DysmorphologyClassificaiton == "Complex"],
#                                                                 score.out$noncoding_denovo_SNVs[score.out$DysmorphologyClassificaiton == "Essential"], "greater")$p.value, 
#                                                     format = "e", digits = 1)), "\n",
#                     sprintf("Co-Eq's p=%s", formatC(wilcox.test(score.out$noncoding_denovo_SNVs[score.out$DysmorphologyClassificaiton == "Complex"],
#                                                                 score.out$noncoding_denovo_SNVs[score.out$DysmorphologyClassificaiton == "Equivocal"], "greater")$p.value, 
#                                                     format = "e", digits = 1)), "\n",
#                     sprintf("Eq-Es's p=%s", formatC(wilcox.test(score.out$noncoding_denovo_SNVs[score.out$DysmorphologyClassificaiton == "Equivocal"],
#                                                                 score.out$noncoding_denovo_SNVs[score.out$DysmorphologyClassificaiton == "Essential"], "greater")$p.value, 
#                                                     format = "e", digits = 1)))
# p6.pvalue <-  paste0(sprintf("Co-Es's p=%s", formatC(wilcox.test(score.out$total_risk[score.out$DysmorphologyClassificaiton == "Complex"],
#                                                                  score.out$total_risk[score.out$DysmorphologyClassificaiton == "Essential"], "greater")$p.value, 
#                                                      format = "e", digits = 1)), "\n",
#                      sprintf("Co-Eq's p=%s", formatC(wilcox.test(score.out$total_risk[score.out$DysmorphologyClassificaiton == "Complex"],
#                                                                  score.out$total_risk[score.out$DysmorphologyClassificaiton == "Equivocal"], "greater")$p.value, 
#                                                      format = "e", digits = 1)), "\n",
#                      sprintf("Eq-Es's p=%s", formatC(wilcox.test(score.out$total_risk[score.out$DysmorphologyClassificaiton == "Equivocal"],
#                                                                  score.out$total_risk[score.out$DysmorphologyClassificaiton == "Essential"], "greater")$p.value, 
#                                                      format = "e", digits = 1)))
# 
# p1 <- ggplot(score.out, aes(x = DysmorphologyClassificaiton, y = coding_CNVs_perc, fill = DysmorphologyClassificaiton)) +
#   geom_violin(color = NA) + geom_boxplot(fill = "white", alpha = .75, width = .2) + theme_classic() + 
#   theme(legend.position = "none", axis.title.x = element_blank()) + 
#   annotate("text", label = c("**", ""), x = c(2, 2.5), y = c(1.22, 1.12), cex = 5) +
#   annotate("line", x = c(1, 3, 2, 3), y = c(1.2, 1.2, 1.1, 1.1), group = c("A", "B", "A", "B")) +
#   scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + ylab("Coding CNVs (percentile)")
# p2 <- ggplot(score.out, aes(x = DysmorphologyClassificaiton, y = coding_rare_SNVs_perc, fill = DysmorphologyClassificaiton)) +
#   geom_violin(color = NA) + geom_boxplot(fill = "white", alpha = .75, width = .2) + theme_classic() + theme(legend.position = "none", axis.title.x = element_blank())+ 
#   annotate("text", label = c("*"), x = c(2.5), y = c(1.22), cex = 5) +
#   annotate("line", x = c(1, 3), y = c(1.2, 1.2), group = c("A", "A")) +
#   scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + ylab("Coding rare SNVs (percentile)")
# p3 <- ggplot(score.out, aes(x = DysmorphologyClassificaiton, y = coding_denovo_SNVs_perc, fill = DysmorphologyClassificaiton)) +
#   geom_violin(color = NA) + geom_boxplot(fill = "white", alpha = .75, width = .2) + theme_classic() + 
#   theme(legend.position = "none", axis.title.x = element_blank()) + 
#   annotate("text", label = c("*"), x = c(2.5), y = c(1.22), cex = 5) +
#   annotate("line", x = c(2, 3), y = c(1.2, 1.2), group = c("A", "A")) +
#   scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + ylab("Coding de novo SNVs (percentile)")
# p4 <- ggplot(score.out, aes(x = DysmorphologyClassificaiton, y = noncoding_rare_SNVs_perc, fill = DysmorphologyClassificaiton)) +
#   geom_violin(color = NA) + geom_boxplot(fill = "white", alpha = .75, width = .2) + theme_classic() + 
#   theme(legend.position = "none", axis.title.x = element_blank()) +
#   scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + ylab("Noncoding rare SNVs (percentile)")
# p5 <- ggplot(score.out, aes(x = DysmorphologyClassificaiton, y = noncoding_denovo_SNVs_perc, fill = DysmorphologyClassificaiton)) +
#   geom_violin(color = NA) + geom_boxplot(fill = "white", alpha = .75, width = .2) + theme_classic() + 
#   theme(legend.position = "none", axis.title.x = element_blank()) +
#   scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + ylab("Noncoding denovo SNVs (percentile)")
# p6 <- ggplot(score.out, aes(x = DysmorphologyClassificaiton, y = total_risk_perc, fill = DysmorphologyClassificaiton)) +
#   geom_violin(color = NA) + geom_boxplot(fill = "white", alpha = .75, width = .2) + theme_classic() + 
#   theme(legend.position = "none", axis.title.x = element_blank()) + ylab("GRS percentile") +
#   annotate("text", label = c("*"), x = c(1.5), y = c(1.12), cex = 5) +
#   annotate("line", x = c(1, 2), y = c(1.1, 1.1), group = c("A", "A")) +
#   scale_fill_manual(values = c("#009eed", "#fff700", "#fc0000"))
# p6  
# ggsave("grs.grouped.dys.p0.1.3subtypes.in.plot.pdf", width = 4, height = 4)
# score.out <- score.out[!score.out$Sample %in% cli.out$Sample, ]
# test <- wilcox.test(score.out$total_risk[score.out$subtype == "Dysmorphic"],
#                     score.out$total_risk[score.out$subtype == "Nondysmorphic"], "greater")$p.value
  # scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0th", "25th", "50th", "75th", "100th")) + ylab("Total risk score (percentile)")
# score.out$DysmorphologyClassificaiton <- as.character(score.out$DysmorphologyClassificaiton)
# 
# score.out$subtype <- factor(score.out$subtype,
#                                                 levels = c("Nondysmorphic", "Dysmorphic"))
# ggplot(score.out, aes(x = subtype, y = total_risk_perc, fill = subtype)) +
#   geom_violin(color = NA) + geom_boxplot(fill = "white", alpha = .75, width = .2) + theme_classic() + 
#   theme(legend.position = "none", axis.title.x = element_blank()) + ylab("GRS percentile") +
#   annotate("text", label = sprintf("p=%s", signif(test, digits=2)), x = c(1.5), y = c(1.12), cex = 5) +
#   annotate("line", x = c(1, 2), y = c(1.05, 1.05), group = c("A", "A")) +
#   scale_fill_manual(values = c("#009eed", "#F28602"))
# ggsave("grs.grouped.dys.p0.1.pdf", width = 4, height = 4)

# 
# if(suffix == ""){
#   plot_grid(p1, p2, p3, p5, p6, nrow = 2)
#   width = 8
#   height = 7
# }else{
#   p6
#   width = 5
#   height = 4
# }
# ggsave(sprintf("replication.0.01_30x10fold.coding%sall.pdf", suffix), width = width, height = height)
# 
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
# sig.set.count$set <- factor(sig.set.count$set, levels = c(sig.set.count$set[sig.set.count$type == "CNVs"], 
#                                                           sig.set.count$set[sig.set.count$type == "rare SNVs"],
#                                                           sig.set.count$set[sig.set.count$type == "de novo SNVs"]))
# sig.set.count$region <- ifelse(sig.set.count$set %in% names(snv.nc.rare), "noncoding", "coding")
# 
# ggplot(sig.set.count, aes(x = set, y = Freq, fill = type)) + geom_bar(stat = "identity") + xlab("") + theme_classic() + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) + facet_grid(region~.)
# ggsave(sprintf("30x10fold_selected_genesets.coding%s.pdf", suffix), width = 7, height =5)
# 
# ### boxplot
# inc <- ggplot(score.out, aes(x = DysmorphologyClassificaiton, y = all_increasing_severity_variants)) + 
#   geom_boxplot(outlier.alpha = 0) + theme(axis.title.x = element_blank()) + ylab("#variants in positively correlated sets") +
#   geom_jitter(aes(fill = total_risk), shape = 21, size = 4, alpha = .8) + scale_fill_gradientn(values = c(0, 0.45, 0.9, 1.1), 
#                                                                                                colours = c("white", "yellow", "red", "red"))
# dec <- ggplot(score.out, aes(x = DysmorphologyClassificaiton, y = all_decreasing_severity_variants)) + 
#   geom_boxplot(outlier.alpha = 0) + theme(axis.title.x = element_blank()) + ylab("#variants in negatively correlated sets") +
#   geom_jitter(aes(fill = total_risk), shape = 21, size = 4, alpha = .8) + scale_fill_gradientn(values = c(0, 0.45, 0.9, 1.1), 
#                                                                                                colours = c("white", "yellow", "red", "red"))
# 
# plot_grid(inc, dec, nrow = 1)
# ggsave(sprintf("boxplot.numvariant%s.pdf",suffix), width = 14, height =7)
# 
# fig1 <- ggplot(score.out, aes(x = DysmorphologyClassificaiton, y = total_risk_perc, fill = DysmorphologyClassificaiton)) +
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
# ggsave("score.with.crv.pdf", height = 4, width = 12)


### other variants

# dt <- read.delim("../snvAnalysis/rare.snv.signal.tsv", stringsAsFactors = F)
# dt <- dt[dt$Sample %in% c("3_0269_000", "3_0406_000", "3_0391_000", "3_0392_000", "3_0368_000", "3_0207_000"), c("Sample", "chr", "start", "end", "var_type", "gene_symbol", "type", "geneset")]
# write.table(dt, "additional_variant_selected_samples.tsv", sep="\t", row.names=F, quote=F, col.names=T)
# 
# 
# 
# ##### testing p-value cutoff for different var
# samples <- read.delim("sample.with.score.2groups.new.tsv", stringsAsFactors = F)
# scores <- c("noncoding_rare_CNVs", "small_noncoding_rare_CNVs", "coding_rare_CNVs", 
#             "small_coding_rare_CNVs", "coding_rare_SNVs", "coding_denovo_SNVs", "noncoding_rare_SNVs", "noncoding_denovo_SNVs")
# out.perf <- data.frame()
# for(p in c(1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001)){
#   for(score in scores){
#     score.out <- data.table::fread("score.from.30x.10fold.2groups.new.tsv", data.table = F)
#     score.out <- score.out[score.out$pvalue == p, ]
#     
#     tmp.out <- data.frame()
#     for(s in unique(score.out$Sample)){
#       tmp <- score.out[score.out$Sample == s, ]
#       tmp <- colSums(tmp[, -1], na.rm = T)/nrow(tmp)
#       
#       tmp.out <- rbind(tmp.out, data.frame("Sample" = s, t(tmp)))
#     }
#     
#     score.out <- tmp.out
#     score.out <- merge(score.out, samples, by = "Sample", all = T)
#     score.out$test_score <- score.out[, score]
#     # score.out$MultiClass <- factor(score.out$MultiClass)
#     # score.out$total_risk <- score.out$coding_rare_CNVs + score.out$small_coding_rare_CNVs + 
#     #   score.out$coding_rare_SNVs + score.out$coding_denovo_SNVs
#     lm <- glm(GroupedDysmorphology ~ test_score, data = score.out)
#     
#     if(nrow(summary(lm)$coefficients) > 1){
#       p.value <- paste0(sprintf("p=%s", formatC(signif(summary(lm)$coefficients[2, "Pr(>|t|)"]), 2), 
#                                 format = "e", digits = 1))
#       out.perf <- rbind(out.perf, data.frame("cutoff" = p, "score" = score,
#                                              "NagelkerkeR2" = signif(NagelkerkeR2(lm)$R2, digits = 3),
#                                              "p" = p.value))
#     }
#   }
# }
# 
# out.perf$cutoff <- factor(out.perf$cutoff, levels = unique(rev(out.perf$cutoff)))
# ggplot(out.perf, aes(x = paste0("p<",cutoff), y = NagelkerkeR2, label = p)) + facet_wrap(score~., nrow = 4, scales = "free_y") +
#   geom_bar(stat = "identity", color = "black", fill = "red") + xlab("p cut-off") + 
#   theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_text(nudge_y = 0.005)
