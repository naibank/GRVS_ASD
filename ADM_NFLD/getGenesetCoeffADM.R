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
    ref <- sprintf("factor(%s) ~ %s + %s", class, paste(covariates, collapse=" + "), global)
    add <- sprintf("%s + %s", ref, feat)
    
    ref <- glm(ref, dt, family = binomial(link = "logit"))
    add <- glm(add, dt, family = binomial(link = "logit"))
    pvalue <- anova(ref, add, test = "Chisq")$`Pr(>Chi)`[2]
    
    # print(sprintf("%s's pvalue = %s", feat, pvalue))
    if(!is.na(pvalue) & pvalue <= cutoff){
      r2 <- c(r2, NagelkerkeR2(add)$R2)
      
      add <- sprintf("factor(%s) ~ %s + %s + %s", class, paste(covariates, collapse=" + "), global, feat)
      add <- glm(add, dt, family = binomial(link = "logit"))
      
      coeff <- c(coeff, add$coefficients[feat])
      p <- c(p, pvalue)
      
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

load("../requiredData/gsNFLD.RData")

cds <- read.delim("../requiredData/hg19_refGene_28_04_17.cds.txt", stringsAsFactors = F, header = F)
cds.g <- GRanges(cds$V1, IRanges(cds$V2, cds$V3), "*")


all.cnv <- rbind(read.delim("../data/all.cnvs.may6.txt", stringsAsFactors = F),
                 read.delim("../data/all.small.cnvs.may6.txt", stringsAsFactors = F))
# all.cnv <- read.delim("../data/all.cnvs.nov27_10kb_3mb.txt", stringsAsFactors = F)

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

rare.lof <- read.delim("../data/SNVs/calls/rare/lof.variants.DEC062019.tsv", stringsAsFactors = F)
rare.ms1 <- read.delim("../data/SNVs/calls/rare/ms1.variants.DEC062019.tsv", stringsAsFactors = F)
rare.ms2 <- read.delim("../data/SNVs/calls/rare/ms2.variants.DEC062019.tsv", stringsAsFactors = F)
rare.lof$effect.tier <- "lof"
rare.ms1$effect.tier <- "tier1_ms"
rare.ms2$effect.tier <- "tier2_ms"
all.rare <- rbind(rare.lof, rbind(rare.ms1, rare.ms2))

denovo.lof <- read.delim("../data/SNVs/calls/denovo/lof.variants.DEC062019.tsv", stringsAsFactors = F)
denovo.ms1 <- read.delim("../data/SNVs/calls/denovo/ms1.variants.DEC062019.tsv", stringsAsFactors = F)
denovo.ms2 <- read.delim("../data/SNVs/calls/denovo/ms2.variants.DEC062019.tsv", stringsAsFactors = F)
denovo.lof$effect.tier <- "lof"
denovo.ms1$effect.tier <- "tier1_ms"
denovo.ms2$effect.tier <- "tier2_ms"
all.denovo <- rbind(denovo.lof, rbind(denovo.ms1, denovo.ms2))

### read all 5 files
cnv.coding <- read.delim("ADM/cnv.coding.matrix.tsv", stringsAsFactors = F)
cnv.coding$sample <- gsub("-", "_", cnv.coding$sample)
cnv.coding$sample <- gsub("A|_A", "", cnv.coding$sample)
names(cnv.coding)[1] <- "Sample"

scnv.coding <- read.delim("ADM/smaller.cnv.coding.matrix.tsv", stringsAsFactors = F)
scnv.coding$sample <- gsub("-", "_", scnv.coding$sample)
scnv.coding$sample <- gsub("A|_A", "", scnv.coding$sample)
names(scnv.coding)[1] <- "Sample"

### noncoding cnv
cnv.ncoding <- read.delim("../dataNC/nc.cnv.matrix.tsv", stringsAsFactors = F)
scnv.ncoding <- read.delim("../dataNC/nc.smaller.cnv.matrix.tsv", stringsAsFactors = F)

snv.coding.rare <- read.delim("rare.coding.snv.tsv", stringsAsFactors = F)

meta <- read.delim("../data/2021.02.01-updated NFLD phenotype table IQ ADM and 3 status.txt", stringsAsFactors = F)
meta$WGS_ManuID <- gsub("-", "_", meta$WGS_ManuID)
meta$WGS_ManuID <- gsub("A|_A", "", meta$WGS_ManuID)
meta <- meta[meta$incl.aff.in.my.study == 1, ]
meta$ADM[grep("0456|0458", meta$WGS_ManuID)] <- "Nondysmorphic"
meta <- meta[which(meta$ADM %in% c("Dysmorphic", "Nondysmorphic")), ]
meta$ADM <- ifelse(meta$ADM == "Dysmorphic", 1, 0)
meta <- meta[meta$Platform %in% c("BGI - Illumina HiSeq", "Macrogen-Illumina HiSeqX", "TCAG-Illumina HiSeqX"), ]
#fix PCs for four samples in cnv.coding
cnv.coding <- merge(cnv.coding, meta[, c("WGS_ManuID", "ADM")], by.x = "Sample", by.y = "WGS_ManuID", all = F)
snv.coding.rare <- merge(snv.coding.rare, meta[, c("WGS_ManuID", "ADM")], by.x = "Sample", by.y = "WGS_ManuID", all = F)

# cnv.coding <- merge(cnv.coding[, -c(124:128)], snv.coding.rare[, c(1, 129:136)])
cnv.coding <- merge(cnv.coding[, -c(84:88)], snv.coding.rare[, c(1, 129:136)])

scnv.coding <- merge(scnv.coding, meta[, c("WGS_ManuID", "ADM")], by.x = "Sample", by.y = "WGS_ManuID", all = F)
scnv.coding <- merge(scnv.coding[, -c(84:88)], snv.coding.rare[, c(1, 129:136)])

cnv.ncoding <- merge(cnv.ncoding, meta[, c("WGS_ManuID", "ADM")], by.x = "Sample", by.y = "WGS_ManuID", all = F)
scnv.ncoding <- merge(scnv.ncoding, meta[, c("WGS_ManuID", "ADM")], by.x = "Sample", by.y = "WGS_ManuID", all = F)

cnv.ncoding <- merge(cnv.ncoding[, -c(22)], snv.coding.rare[, c(1, 120, 125, 129:136)])
scnv.ncoding <- merge(scnv.ncoding[, -c(22)], snv.coding.rare[, c(1, 120, 125, 129:136)])

snv.coding.denovo <- read.delim("denovo.coding.snv.tsv", stringsAsFactors = F)
snv.coding.denovo <- merge(snv.coding.denovo, meta[, c("WGS_ManuID", "ADM")], by.x = "Sample", by.y = "WGS_ManuID", all = F)

snv.nc.rare <- read.delim("../scriptsNC/nc.snv.rare.matrix.tsv", stringsAsFactors = F)
snv.nc.rare <- merge(snv.nc.rare, meta[, c("WGS_ManuID", "ADM")], by.x = "Sample", by.y = "WGS_ManuID", all = F)

snv.nc.denovo <- read.delim("../scriptsNC/nc.snv.denovo.matrix.tsv", stringsAsFactors = F)
snv.nc.denovo <- merge(snv.nc.denovo, meta[, c("WGS_ManuID", "ADM")], by.x = "Sample", by.y = "WGS_ManuID", all = F)

# for(test in c("MultiClass", "ComplexVsEssential", "GroupedDysmorphology", "DysmorphologyScore")){
test <- "ADM"

baseline.feat <- c("Sample", "Sex..CRV",  "Platform", test, "pc1", "pc2", "pc3")
samples <- na.omit(unique(rbind(cnv.coding[, baseline.feat], rbind(snv.coding.rare[, baseline.feat],
                                                    rbind(snv.coding.denovo[, baseline.feat],
                                                          rbind(snv.nc.rare[, baseline.feat],
                                                                snv.nc.denovo[, baseline.feat]))))))

samples <- na.omit(unique(rbind(cnv.ncoding[, baseline.feat], rbind(samples, scnv.ncoding[, baseline.feat]))))
samples <- samples[samples$Sample %in% snv.coding.denovo$Sample, ]

# samples <- na.omit(unique(rbind(cnv.coding[, baseline.feat], snv.coding.rare[, baseline.feat])))
set.seed(1500)
score.out <- data.frame()
sig.set <- c()
coeff.all <- data.frame()
covariates <- c("Sex..CRV", "Platform", "pc1", "pc2", "pc3")

baseline <- glm(sprintf("%s ~ Sex..CRV + Platform + pc1 + pc2 + pc3", test), data = samples)

## get significant geneset
cnv.del <- names(cnv.coding)[grep("DEL", names(cnv.coding))]
cnv.dup <- names(cnv.coding)[grep("DUP", names(cnv.coding))]
cnv1 <- getSigFeatures(cnv.coding, "gene_count_DUP", covariates, cnv.dup[-c(1:3)])
cnv2 <- getSigFeatures(cnv.coding, "gene_count_DEL", covariates, cnv.del[-c(1:3)])
# cnv3 <- getSigFeatures(cnv.coding, "gene_count_PartialDup", covariates, names(cnv.coding)[85:121])

cnv.c.sig <- c(cnv1$set, cnv2$set)
cnv.set <- rbind(cnv1$coeff, cnv2$coeff)

if(length(cnv1$set) > 0)
  coeff.all <- rbind(coeff.all, data.frame(cnv1$coeff, "variant" = "coding_rare_CNVs", stringsAsFactors = F))
if(length(cnv2$set) > 0)
  coeff.all <- rbind(coeff.all, data.frame(cnv2$coeff, "variant" = "coding_rare_CNVs", stringsAsFactors = F))

scnv.del <- names(scnv.coding)[grep("DEL", names(scnv.coding))]
scnv.dup <- names(scnv.coding)[grep("DUP", names(scnv.coding))]
scnv1 <- getSigFeatures(scnv.coding, "gene_count_DUP", covariates, scnv.dup[-c(1:3)])
scnv2 <- getSigFeatures(scnv.coding, "gene_count_DEL", covariates, scnv.del[-c(1:3)])
# cnv3 <- getSigFeatures(cnv.coding, "gene_count_PartialDup", covariates, names(cnv.coding)[85:121])

scnv.c.sig <- c(cnv1$set, cnv2$set)
scnv.set <- rbind(cnv1$coeff, cnv2$coeff)

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
snv1 <- getSigFeatures(snv.coding.rare, "Total_lof", covariates, names(snv.coding.rare)[4:40], class = test)
snv2 <- getSigFeatures(snv.coding.rare, "Total_tier1_ms", covariates, names(snv.coding.rare)[42:78], class = test)
snv3 <- getSigFeatures(snv.coding.rare, "Total_tier2_ms", covariates, names(snv.coding.rare)[80:116], class = test)
snv.c.rare.sig <- c(snv1$set, snv2$set, snv3$set)
rare.set <- rbind(snv1$coeff, rbind(snv2$coeff, snv3$coeff))

if(length(snv1$set) > 0)
  coeff.all <- rbind(coeff.all, data.frame(snv1$coeff, "variant" = "coding_rare_SNVs", stringsAsFactors = F))
if(length(snv2$set) > 0)
  coeff.all <- rbind(coeff.all, data.frame(snv2$coeff, "variant" = "coding_rare_SNVs", stringsAsFactors = F))
if(length(snv3$set) > 0)
  coeff.all <- rbind(coeff.all, data.frame(snv3$coeff, "variant" = "coding_rare_SNVs", stringsAsFactors = F))

snv1 <- getSigFeatures(snv.coding.denovo, "Total_lof", covariates, names(snv.coding.denovo)[4:40], class = test)
snv2 <- getSigFeatures(snv.coding.denovo, "Total_tier1_ms", covariates, names(snv.coding.denovo)[42:78], class = test)
snv3 <- getSigFeatures(snv.coding.denovo, "Total_tier2_ms", covariates, names(snv.coding.denovo)[80:116], class = test)
snv.c.denovo.sig <- c(snv1$set, snv2$set, snv3$set)
denovo.set <- rbind(snv1$coeff, rbind(snv2$coeff, snv3$coeff))
if(length(snv1$set) > 0)
  coeff.all <- rbind(coeff.all, data.frame(snv1$coeff, "variant" = "coding_denovo_SNVs", stringsAsFactors = F))
if(length(snv2$set) > 0)
  coeff.all <- rbind(coeff.all, data.frame(snv2$coeff, "variant" = "coding_denovo_SNVs", stringsAsFactors = F))
if(length(snv3$set) > 0)
  coeff.all <- rbind(coeff.all, data.frame(snv3$coeff, "variant" = "coding_denovo_SNVs", stringsAsFactors = F))

snv.nc.rare.sig <- getSigFeatures(snv.nc.rare, "TotalRare", covariates, names(snv.nc.rare)[2:26], class = test)
if(length(snv.nc.rare.sig$set) > 0)
  coeff.all <- rbind(coeff.all, data.frame(snv.nc.rare.sig$coeff, "variant" = "noncoding_rare_SNVs", stringsAsFactors = F))
snv.nc.rare.sig <- snv.nc.rare.sig$set

snv.nc.denovo.sig <- getSigFeatures(snv.nc.denovo, "TotalRare", covariates, names(snv.nc.denovo)[2:26], class = test)
if(length(snv.nc.denovo.sig$set) > 0)
  coeff.all <- rbind(coeff.all, data.frame(snv.nc.denovo.sig$coeff, "variant" = "noncoding_denovo_SNVs", stringsAsFactors = F))
snv.nc.denovo.sig <- snv.nc.denovo.sig$set
  
write.table(coeff.all, sprintf("ADM/coeff.%s.from.all.samples.tsv", test), sep="\t", row.names=F, quote=F, col.names=T)


# }