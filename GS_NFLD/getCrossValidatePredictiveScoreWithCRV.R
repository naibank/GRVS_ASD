setwd("~/Documents/doc/working/NFLD/2019/scripts/")
rm(list=ls())
library(ordinal)
library(ggplot2)
library(ggforce)
library(cowplot)
library(GenomicRanges)
library(reshape)
library(magrittr)

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

getSigFeatures <- function(dt, global, covariates, test, cutoff = 0.01){
  result <- c() 
  coeff <- c()
  for(feat in test){
    ref <- sprintf("factor(MultiClass) ~ %s + %s", paste(covariates, collapse=" + "), global)
    add <- sprintf("%s + %s", ref, feat)
    
    ref <- clm(ref, data = dt)
    add <- clm(add, data = dt)
    pvalue <- anova(ref, add)$`Pr(>Chisq)`[2]
    
    # print(sprintf("%s's pvalue = %s", feat, pvalue))
    if(!is.na(pvalue) & pvalue <= cutoff){
      result <- c(result, feat)
      coeff <- c(coeff, add$coefficients[feat])
    }
  }
  
  coeff <- data.frame("set" = names(coeff), "coeff" = coeff)
  dt.out <- list(result, coeff)
  names(dt.out) <- c("set", "coeff")
  return(dt.out)
}

getScore <- function(dt, baseline, test){
  if(length(test) > 0){
    ref <- gsub("\n", "", sprintf("MultiClass ~ offset(I(%s * (Sex..CRV == 'M') + 
                                  %s * (Platform == 'Macrogen-Illumina HiSeqX') + %s * (Platform == 'Complete Genomics') + 
                                  %s * (Platform == 'TCAG-Illumina HiSeqX') +
                                  %s * (pc1) + 
                                  %s * (pc2) + 
                                  %s * (pc3)))", 
                                  baseline$coefficients["Sex..CRVM"], baseline$coefficients["PlatformMacrogen-Illumina HiSeqX"], baseline$coefficients["PlatformComplete Genomics"], 
                                  baseline$coefficients["PlatformTCAG-Illumina HiSeqX"], baseline$coefficients["pc1"], baseline$coefficients["pc2"], baseline$coefficients["pc3"]))
    
    add <- sprintf("%s + %s", ref, paste(test, collapse = " + "))
    lm <- glm(add, data = dt)
    sm <- summary(lm)
    sig <- rownames(sm$coefficients)[(sm$coefficients[, 'Pr(>|t|)'] < 1)]
    
    return(lm$coefficients[sig])
  }
}

load("../../data/GeneSet/gsNFLD.RData")

cds <- read.delim("~/Documents/CentOS_Documents/CNVs/data/annotateCNV/FromJohn/hg19_refGene_28_04_17.cds.txt", stringsAsFactors = F, header = F)
cds.g <- GRanges(cds$V1, IRanges(cds$V2, cds$V3), "*")

all.cnv <- read.delim("../data/all.cnvs.may6.txt", stringsAsFactors = F)
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

rare.lof <- read.delim("../data/SNVs/calls/rare/lof.variants.MAY062019.tsv", stringsAsFactors = F)
rare.ms1 <- read.delim("../data/SNVs/calls/rare/ms1.variants.MAY062019.tsv", stringsAsFactors = F)
rare.ms2 <- read.delim("../data/SNVs/calls/rare/ms2.variants.MAY062019.tsv", stringsAsFactors = F)
rare.lof$effect.tier <- "lof"
rare.ms1$effect.tier <- "tier1_ms"
rare.ms2$effect.tier <- "tier2_ms"
all.rare <- rbind(rare.lof, rbind(rare.ms1, rare.ms2))

denovo.lof <- read.delim("../data/SNVs/calls/denovo/lof.variants.MAY062019.tsv", stringsAsFactors = F)
denovo.ms1 <- read.delim("../data/SNVs/calls/denovo/ms1.variants.MAY062019.tsv", stringsAsFactors = F)
denovo.ms2 <- read.delim("../data/SNVs/calls/denovo/ms2.variants.MAY062019.tsv", stringsAsFactors = F)
denovo.lof$effect.tier <- "lof"
denovo.ms1$effect.tier <- "tier1_ms"
denovo.ms2$effect.tier <- "tier2_ms"
all.denovo <- rbind(denovo.lof, rbind(denovo.ms1, denovo.ms2))

### read all 5 files
cnv.coding <- read.delim("../cnv.coding.matrix.tsv", stringsAsFactors = F)
cnv.coding$sample <- gsub("-", "_", cnv.coding$sample)
cnv.coding$sample <- gsub("A|_A", "", cnv.coding$sample)
names(cnv.coding)[1] <- "Sample"

snv.coding.rare <- read.delim("rare.coding.snv.tsv", stringsAsFactors = F)
#fix PCs for four samples in cnv.coding
cnv.coding <- merge(cnv.coding[, -c(87:93)], snv.coding.rare[, c(1, 129:135)])

snv.coding.denovo <- read.delim("denovo.coding.snv.tsv", stringsAsFactors = F)
snv.nc.rare <- read.delim("../scriptsNC/nc.snv.rare.matrix.tsv", stringsAsFactors = F)
snv.nc.denovo <- read.delim("../scriptsNC/nc.snv.denovo.matrix.tsv", stringsAsFactors = F)

baseline.feat <- c("Sample", "Sex..CRV",  "Platform", "MultiClass", "pc1", "pc2", "pc3")
samples <- na.omit(unique(rbind(cnv.coding[, baseline.feat], rbind(snv.coding.rare[, baseline.feat],
                                                                   rbind(snv.coding.denovo[, baseline.feat],
                                                                         rbind(snv.nc.rare[, baseline.feat],
                                                                               snv.nc.denovo[, baseline.feat]))))))

crv.cnv <- read.delim("CRV.cnvs.tsv", stringsAsFactors = F)
crv.cnv$sample <- gsub("-", "_", crv.cnv$sample)

crv.rare.snv <- read.delim("CRV.rare.snvs.tsv", stringsAsFactors = F)
crv.denovo.snv <- read.delim("CRV.denovo.snvs.tsv", stringsAsFactors = F)

# samples <- na.omit(unique(rbind(cnv.coding[, baseline.feat], snv.coding.rare[, baseline.feat])))
set.seed(1400)
score.out <- data.frame()
sig.set <- c()
coeff.all <- data.frame()
covariates <- c("Sex..CRV", "Platform", "pc1", "pc2", "pc3")

score.tmp <- data.frame()
variant <- data.frame()
for(k in 1:30){
  nfold <- 10
  fold.ess <- split(sample(samples$Sample[samples$MultiClass == 0]), 1:nfold)
  fold.equ <- split(sample(samples$Sample[samples$MultiClass == 1]), 1:nfold)
  fold.com <- split(sample(samples$Sample[samples$MultiClass == 2]), 1:nfold)
  
  for(i in 1:length(fold.ess)){
    disc.set <- c(unlist(fold.ess[-i]), unlist(fold.equ[-i]), unlist(fold.com[-i]))
    targ.set <- c(unlist(fold.ess[i]), unlist(fold.equ[i]), unlist(fold.com[i]))
    
    baseline <- glm("MultiClass ~ Sex..CRV + Platform + pc1 + pc2 + pc3", data = samples[samples$Sample %in% disc.set, ])
    
    ## get significant geneset
    cnv1 <- getSigFeatures(cnv.coding[cnv.coding$Sample %in% disc.set, ], "gene_count_deletion", covariates, names(cnv.coding)[5:41])
    cnv2 <- getSigFeatures(cnv.coding[cnv.coding$Sample %in% disc.set, ], "gene_count_duplication", covariates, names(cnv.coding)[45:81])
    cnv.c.sig <- c(cnv1$set, cnv2$set)
    cnv.set <- rbind(cnv1$coeff, cnv2$coeff)
      
    
    cnv.tmp <- data.frame()
    if(nrow(cnv.set) > 0){
      cnv.coeff <- getScore(cnv.coding[cnv.coding$Sample %in% disc.set, ], baseline, cnv.set$set)
      cnv.coeff <- cnv.coeff[-1]
      for(j in 1:length(cnv.coeff)){
        tmp.set <- cnv.coeff[j]
        type <- strsplit(names(tmp.set), "_")[[1]]
        set <- paste(type[1:(length(type)-1)], collapse = "_")
        type <- type[length(type)]
        tmp.crv.cnv <- crv.cnv[crv.cnv$sample %in% targ.set & crv.cnv$type == type, c("sample", set)]
        tmp.crv.cnv$score <- tmp.crv.cnv[, set] * tmp.set
        
        tmp.var <- crv.cnv[crv.cnv$sample %in% targ.set & crv.cnv$type == type, c("sample", "chr", "start", "end", "type", set)]
        tmp.var$score <- tmp.var[, set] * tmp.set
        if(nrow(tmp.var) > 0){
        tmp.var$set <- set
        tmp.var <- tmp.var[, -6]
        variant <- rbind(variant, tmp.var)
        }
        cnv.tmp <- rbind(cnv.tmp, tmp.crv.cnv[, -2])
      }
    }
    
    if(nrow(cnv.tmp) > 0){
      cnv.tmp$variant_type <- "CNVs"
      names(cnv.tmp)[1] <- "Sample"
      score.tmp <- rbind(score.tmp, cnv.tmp)
    }
    
    snv1 <- getSigFeatures(snv.coding.rare[snv.coding.rare$Sample %in% disc.set, ], "Total_lof", covariates, names(snv.coding.rare)[4:40])
    snv2 <- getSigFeatures(snv.coding.rare[snv.coding.rare$Sample %in% disc.set, ], "Total_tier1_ms", covariates, names(snv.coding.rare)[42:78])
    snv3 <- getSigFeatures(snv.coding.rare[snv.coding.rare$Sample %in% disc.set, ], "Total_tier2_ms", covariates, names(snv.coding.rare)[80:116])
    snv.c.rare.sig <- c(snv1$set, snv2$set, snv3$set)
    rare.set <- rbind(snv1$coeff, rbind(snv2$coeff, snv3$coeff))
    
    rare.tmp <- data.frame()
    if(nrow(rare.set) > 0){
      rare.coeff <- getScore(snv.coding.rare[snv.coding.rare$Sample %in% disc.set, ], baseline, rare.set$set)
      rare.coeff <- rare.coeff[-1]
      for(j in 1:length(rare.coeff)){
        tmp.set <- rare.coeff[j]
        type <- strsplit(names(tmp.set), "_")[[1]]
        if(type[length(type) - 1] %in% c("tier1", "tier2")){
          set <- paste(type[1:(length(type)-2)], collapse = "_")
          type <- paste(type[(length(type)-1):(length(type))], collapse = "_")
        }else{
          set <- paste(type[1:(length(type)-1)], collapse = "_")
          type <- type[length(type)]
        }
        
        tmp.rare.snv <- crv.rare.snv[crv.rare.snv$Sample %in% targ.set & crv.rare.snv$effect.tier == type, c("Sample", set)]
        tmp.rare.snv$score <- tmp.rare.snv[, set] * tmp.set
        
        tmp.var <- crv.rare.snv[crv.rare.snv$Sample %in% targ.set & crv.rare.snv$effect.tier == type, c("Sample", "chr", "start", "end", "var_type",  set)]
        names(tmp.var) <- c("sample", "chr", "start", "end", "type", set)
        tmp.var$score <- tmp.var[, set] * tmp.set
        if(nrow(tmp.var) > 0){
        tmp.var$set <- set
        tmp.var <- tmp.var[, -6]
        variant <- rbind(variant, tmp.var)
        }
        rare.tmp <- rbind(rare.tmp, tmp.rare.snv[, -2])
      }
    }
    
    if(nrow(rare.tmp) > 0){
      rare.tmp$variant_type <- "rare SNVs"
      score.tmp <- rbind(score.tmp, rare.tmp)
    }
    
    
    snv1 <- getSigFeatures(snv.coding.denovo[snv.coding.denovo$Sample %in% disc.set, ], "Total_lof", covariates, names(snv.coding.denovo)[4:40])
    snv2 <- getSigFeatures(snv.coding.denovo[snv.coding.denovo$Sample %in% disc.set, ], "Total_tier1_ms", covariates, names(snv.coding.denovo)[42:78])
    snv3 <- getSigFeatures(snv.coding.denovo[snv.coding.denovo$Sample %in% disc.set, ], "Total_tier2_ms", covariates, names(snv.coding.denovo)[80:116])
    snv.c.denovo.sig <- c(snv1$set, snv2$set, snv3$set)
    denovo.set <- rbind(snv1$coeff, rbind(snv2$coeff, snv3$coeff))
    
    denovo.tmp <- data.frame()
    if(nrow(denovo.set) > 0){
      denovo.coeff <- getScore(snv.coding.denovo[snv.coding.denovo$Sample %in% disc.set, ], baseline, denovo.set$set)
      denovo.coeff <- denovo.coeff[-1]
      for(j in 1:length(denovo.coeff)){
        tmp.set <- denovo.coeff[j]
        type <- strsplit(names(tmp.set), "_")[[1]]
        if(type[length(type) - 1] %in% c("tier1", "tier2")){
          set <- paste(type[1:(length(type)-2)], collapse = "_")
          type <- paste(type[(length(type)-1):(length(type))], collapse = "_")
        }else{
          set <- paste(type[1:(length(type)-1)], collapse = "_")
          type <- type[length(type)]
        }
          
        tmp.denovo.snv <- crv.denovo.snv[crv.denovo.snv$Sample %in% targ.set & crv.denovo.snv$effect.tier == type, c("Sample", set)]
        tmp.denovo.snv$score <- tmp.denovo.snv[, set] * tmp.set
        
        tmp.var <- crv.denovo.snv[crv.denovo.snv$Sample %in% targ.set & crv.denovo.snv$effect.tier == type, c("Sample", "chr", "start", "end", "var_type",  set)]
        names(tmp.var) <- c("sample", "chr", "start", "end", "type", set)
        tmp.var$score <- tmp.var[, set] * tmp.set
        if(nrow(tmp.var) > 0){
        tmp.var$set <- set
        tmp.var <- tmp.var[, -6]
        variant <- rbind(variant, tmp.var)
        }
        denovo.tmp <- rbind(denovo.tmp, tmp.denovo.snv[, -2])
      }
    }
    
    if(nrow(denovo.tmp) > 0){
      denovo.tmp$variant_type <- "denovo SNVs"
      score.tmp <- rbind(score.tmp, denovo.tmp)
    }
  }
}

write.table(variant, "variant.with.score.tsv", sep="\t", row.names=F, quote=F, col.names=T)

score.tmp$id <- paste(score.tmp$Sample, score.tmp$variant_type, sep="#")
write.table(score.tmp, "CRV.score.tmp.tsv", sep="\t", row.names=F, quote=F, col.names=T)

score.tmp <- read.delim("CRV.score.tmp.tsv", stringsAsFactors = F)
tmp <- aggregate(score ~ id, score.tmp, sum)
tmp$score <- tmp$score/30

tmp$sample <- sapply(tmp$id, strsplit, "#") %>% sapply("[", 1)
tmp <- aggregate(score ~ sample, tmp, sum)

all.score <- read.delim("score.sample.with.CRV.tsv", stringsAsFactors = F)

tmp <- merge(tmp, all.score[, c("Sample", "total_risk", "DysmorphologyClassificaiton")], by.x = "sample", by.y = "Sample", all.x = T)
names(tmp) <- c("sample", "CRV_score", "total_score", "classification")

tmp <- tmp[tmp$CRV_score != 0, ]

tmp.mt <- melt(tmp, id.vars = c("sample", "classification"))
tmp.mt$classification <- factor(tmp.mt$classification, levels = c("Essential", "Equivocal", "Complex"))
tmp.mt$sample <- factor(tmp.mt$sample, levels = unique(tmp.mt$sample[order(tmp.mt$classification, tmp.mt$value)]))
ggplot(tmp.mt, aes(x = sample, y = value, fill = variable)) + geom_bar(stat = "identity", position = "identity", alpha = .5) +
  scale_fill_manual(values = c("red", "darkgrey")) + theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = 0)) +
  geom_vline(xintercept = c(6.5, 8.5), lty = 2) + ylab("Score")
ggsave("crv.contribution.pdf", height = 8, width = 8)
