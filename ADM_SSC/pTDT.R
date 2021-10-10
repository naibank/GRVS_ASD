setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# rm(list=ls())
library(ggplot2)
library(data.table)
##### SSC
check <- c("SSC06730", "SSC05954")
ssc.tmp <- fread("../PRS/SSC.PRS.all.score", data.table = F)
names(ssc.tmp) <- ssc.tmp[1, ]
ssc.tmp <- ssc.tmp[-1, ]

score <- read.delim("../../../MSSNG/pTDT/SSC_metadata+PRS_total.csv", stringsAsFactors = F, sep=",")
score <- merge(score, ssc.tmp[, c("IID", "0.1")], by.x = "SampleID", by.y="IID", all.x = T)
names(score)[23:24] <- c("pold", "SCORESUM")
plot(score$pold, score$SCORESUM)

score$ParentMean <- NA
score$Relation <- gsub("^other", "unaffected", score$Relation)
score <- score[score$Relation %in% c("father", "mother", "proband", "unaffected sibling"), ]
score$SCORESUM <- scale(score$SCORESUM)

for(i in 1:nrow(score)){
  tmp <- score[i, ]
  if(!tmp$Relation %in% c("father", "mother")){
    famid <- tmp$FamilyID
    parent <- score[score$FamilyID == famid & score$Relation %in% c("father", "mother"), ]
    
    if(nrow(parent) == 2){
      average <- mean(parent$SCORESUM)
      score$ParentMean[i] <- average
    }
    
  }
}

score$overtransmission <- score$SCORESUM-score$ParentMean
score <- score[!is.na(score$ParentMean), ]
n442 <- read.delim("../SNVs/n442.SSCaff_unaffsib_samplelist.tsv", stringsAsFactors = F)
score <- score[score$SampleID %in% n442$Sample.ID, ]
grs <- read.delim("ssc.grs.p0.1.tsv", stringsAsFactors = F)

mainmeta <- read.delim("../SSC_data/Target.admixture.tagged.tsv", stringsAsFactors = F)
mainmeta <- read.delim("../../../MSSNG/SampleInfo/SSC.all.metadata.tsv", stringsAsFactors = F)
mainmeta <- mainmeta[mainmeta$Dataset == "SSC", ]


meta <- read.delim("../SSC_data/SSC_meta_data.tsv", stringsAsFactors = F)
score <- merge(score, meta[, -1], by.x = "SampleID", by.y = "Sample.ID", all.x = T)
score$adm.dysmorphic <- as.character(score$adm.dysmorphic)
fam.ids <- score$FamilyID[!(is.na(score$adm.dysmorphic))]

score$adm.dysmorphic[is.na(score$adm.dysmorphic)] <- "unaff"
score$adm.dysmorphic <- factor(score$adm.dysmorphic, levels = c("unaff", "nondys", "dys"))
score <- score[score$FamilyID %in% fam.ids, ]
score.tmp <- score
score.tmp$tmp_dis <- as.character(score.tmp$adm.dysmorphic)
meta.out <- merge(score.tmp[, -c(27:28)], grs, by.x = "SampleID", by.y = "sample", all = T)
meta.out$adm.dysmorphic[is.na(meta.out$adm.dysmorphic)] <- meta.out$tmp_dis[is.na(meta.out$adm.dysmorphic)]
meta.out <- meta.out[, c(1, 24:26, 28:37)]

meta.out <- merge(meta.out, mainmeta[, c("Sample.ID", "Family.ID", "Relation", "Sex", 
                                         "Affection", "Platform", "Library.type", paste0("K", 1:5), "Predicted.ancestry")], by.x = "SampleID", by.y = "Sample.ID", all.x = T)
write.table(meta.out, "ssc.manifest.tsv", sep="\t", row.names=F, quote=F, col.names=T)
# meta <- read.delim("../SSC_data/SSC_meta_data.tsv", stringsAsFactors = F)
# clean.samples <- readLines("n442.final.SSC.samples.txt")
# meta <- meta[meta$Sample.ID %in% clean.samples, ]
# score <- merge(score, meta[, c("Sample.ID", "adm.dysmorphic", "adm.autism_type")], 
               # by.x = "SampleID", by.y = "Sample.ID")
# score <- merge(score, grs[, c("sample", "adm.dysmorphic")], by.x = "SampleID", by.y = "sample")

ggplot(score, 
       aes(x = adm.dysmorphic, y = overtransmission, color = adm.dysmorphic)) + 
  geom_boxplot(outlier.alpha = 0, show.legend = F, width = .4, fill = "white", alpha = .6) + 
  geom_jitter(shape = 21, width = .1, show.legend = F, alpha = .3) + theme_bw() +     
  ylab("PRS difference from the parents average") + geom_hline(yintercept = 0, lty = 2) + xlab("") +
  ggtitle("SSC")


table.out <- data.frame()
groups <- list("dys" = "dys",
               "nondys" = "nondys")
for(i in 1:length(groups)){
  message(names(groups)[i])
  group <- names(groups)[i]
  mean <- mean(score$SCORESUM[score$adm.dysmorphic %in% groups[[i]]], na.rm = T)
  sd <- sd(score$SCORESUM[score$adm.dysmorphic %in% groups[[i]]], na.rm = T)
  
  test.out <- data.frame(group, mean, sd)
  for(j in 1:length(groups)){
    vs_group <- names(groups)[j]
    test.out[, paste0("pvalue_vs_", vs_group)] <- t.test(score$SCORESUM[score$adm.dysmorphic %in% groups[[i]]],
                                                         score$SCORESUM[score$adm.dysmorphic %in% groups[[j]]])$p.value
  }
  test.out$test <- "direct_comparison"
  test.out$cohort <- "SSC"
  table.out <- rbind(table.out, test.out)
  
}


for(i in 1:length(groups)){
  message(names(groups)[i])
  group <- names(groups)[i]
  mean <- mean(score$overtransmission[score$adm.dysmorphic %in% groups[[i]]], na.rm = T)
  sd <- sd(score$overtransmission[score$adm.dysmorphic %in% groups[[i]]], na.rm = T)
  
  test.out <- data.frame(group, mean, sd)
  for(j in 1:length(groups)){
    vs_group <- names(groups)[j]
    test.out[, paste0("pvalue_vs_", vs_group)] <- t.test(score$overtransmission[score$adm.dysmorphic %in% groups[[i]]],
                                                         score$overtransmission[score$adm.dysmorphic %in% groups[[j]]])$p.value
  }
  test.out$test <- "parentalDeviationComparison"
  test.out$cohort <- "SSC"
  table.out <- rbind(table.out, test.out)
}

# t.test(score$SCORESUM[score$adm.dysmorphic == "dys"],
#        score$SCORESUM[score$adm.dysmorphic == "nondys"]) #p = 0.8
# t.test(score$overtransmission[score$adm.dysmorphic == "dys"],
#        score$overtransmission[score$adm.dysmorphic == "nondys"]) #p = 0.051
# 
t.dys <- t.test(score$overtransmission[score$adm.dysmorphic == "dys"], alternative = "greater", mu = 0)
t.nondys <- t.test(score$overtransmission[score$adm.dysmorphic == "nondys"], alternative = "greater", mu = 0)
t.unaff <- t.test(score$overtransmission[score$adm.dysmorphic == "unaff"], alternative = "greater", mu = 0)

score$Relation[score$adm.dysmorphic == "dys"] <- sprintf("ADM\nDysmorphic \n(n=%s)",
                                                         sum(score$adm.dysmorphic == "dys"))
score$Relation[score$adm.dysmorphic == "nondys"] <- sprintf("ADM\nNondysmorphic \n(n=%s)", 
                                                            sum(score$adm.dysmorphic == "nondys"))
score$Relation[score$adm.dysmorphic == "unaff"] <- sprintf("Unaffected siblings \n(n=%s)", 
                                                            sum(score$adm.dysmorphic == "unaff"))
write.table(score, "ssc.score.tsv", sep="\t", row.names=F, quote=T, col.names=T)
score$Relation <- factor(score$Relation, levels = unique(score$Relation)[c(2,1,3)])

ggplot(score, 
       aes(x = Relation, y = overtransmission, fill = Relation)) + geom_violin(color = NA) +
  geom_boxplot(fill = "white", width = .2, alpha=.6) + theme_classic() +     
  ylab("PRS difference from the parents average") + geom_hline(yintercept = 0, lty = 2) +
  annotate("text", x = c(1,2,3), y = c(3),
           label = paste("p=",c(signif(t.unaff$p.value, digits = 2),
                                signif(t.nondys$p.value, digits = 2),
                                signif(t.dys$p.value, digits = 2)), sep="")) + xlab("") +
  scale_fill_manual(values = c("grey", "#009eed", "#F28602")) + theme(legend.position = "none") +
  #scale_color_manual(values = c("grey", "#009eed", "#F28602")) + 
  ggtitle("SSC")
ggsave("ssc.adm.prs.pdf", width = 6, height = 4)



p2 <- ggplot(score, 
       aes(x = Relation, y = overtransmission, color = Relation)) + 
  geom_boxplot(outlier.alpha = 0, show.legend = F, width = .4, fill = "white", alpha = .6) + 
  geom_jitter(color = "black", 
              shape = 21, width = .1, show.legend = F, alpha = .7) + theme_bw() +     
  ylab("PRS difference from the parents average") + geom_hline(yintercept = 0, lty = 2) +
  annotate("text", x = c(1,2), y = c(3),
           label = paste("p=",c(signif(t.dys$p.value, digits = 2),
                     signif(t.nondys$p.value, digits = 2)), sep="")) +
  scale_color_manual(values = c("grey", "#009eed", "#F28602")) + xlab("") +
  scale_fill_manual(values = c("complex" = "red", "essential" = "grey", "equivocal" = "yellow")) +
  ggtitle("SSC") + coord_cartesian(ylim = c(-3, 3))
          
p2
ssc.aff <- score

#################
################ NFLD
score <- data.table::fread("../../main/common/ASD.PRS.2021.all.score", data.table = F)
# score <- data.table::fread("../../common/AdaPRS/ASD.PRS.all.score", data.table = F)

score$IID <- gsub("^0_", "", score$IID)
score$IID <- gsub("_", "-", score$IID)
score$IID <- gsub("--", "-", score$IID)
score$IID <- sapply(score$IID, substr, 1, 10)
names(score) <- c("FID", "WGS_ManuID", "PRS_0_1")
score$PRS_0_1 <- scale(score$PRS_0_1)
score <- score[!duplicated(score$WGS_ManuID), ]

check <- c("3-0066-001", "3-0116-001")
# score <- read.delim("NFLD.PRS.Jan2020.tsv", stringsAsFactors = F)
# score$PRS_0_1 <- scale(score$PRS_0_1)
meta <- read.delim("../../main/data/2021.02.01-updated NFLD phenotype table IQ ADM and 3 status.txt", stringsAsFactors = F)
meta <- meta[meta$incld.cohort == 1, ]
meta$ADM[grep("0456|0458", meta$WGS_ManuID)] <- "Nondysmorphic"
# adm <- read.delim("../scripts/data")
#names(score)[c(26:32)] <- c("PRS_0_1", "PRS_0_2", "PRS_0_3", "PRS_0_4", "PRS_0_5", "PRS_1", "var1") 
# write.table(score, "NFLD.PRS.Jan2020.tsv", sep="\t", row.names=F, quote=F, col.names=T)
score <- merge(score[, -1], meta, by = "WGS_ManuID", all.x = T)

# score <- read.delim("../../main/common/NFLD.PRS.Jan2020.tsv", stringsAsFactors = F)
score$ParentMean <- NA
score$Relation <- tolower(score$Relation)
score$Relation <- gsub("affectedsibling", "proband", score$Relation)
score <- score[score$Relation %in% c("father", "mother", "proband", "unaffected sibling"), ]
score <- score[!is.na(score$PRS_0_1),]
score$PRS_0_1 <- scale(score$PRS_0_1)

for(i in 1:nrow(score)){
  tmp <- score[i, ]
  if(!tmp$Relation %in% c("father", "mother")){
    famid <- tmp$WGS_FamilyID
    parent <- score[score$WGS_FamilyID == famid & score$Relation %in% c("father", "mother"), ]
    
    if(nrow(parent) == 2){
      average <- mean(parent$PRS_0_1)
      score$ParentMean[i] <- average
      
    }
  }
}

score$overtransmission <- score$PRS_0_1-score$ParentMean
score <- score[!is.na(score$ParentMean), ]
score <- score[score$Platform != "Complete Genomics", ]

# meta <- read.delim("../../main/data/2021.02.01-updated NFLD phenotype table IQ ADM and 3 status.txt", stringsAsFactors = F)
# meta <- meta[!is.na(meta$ADM), ]
# score <- merge(score, meta[, c("WGS_ManuID", "ADM")], 
#                by.x = "WGS_ManuID", by.y = "WGS_ManuID")

eth <- read.delim("../../main/NFLD_ethnicity_admixture_and_eigensoft/nfld.eth.tag.tsv", stringsAsFactors = F)
score <- score[score$WGS_ManuID %in% eth$sample[eth$clusterPop == "EUR"], ]

groups <- list("dys" = "Dysmorphic",
               "nondys" = "Nondysmorphic")
for(i in 1:length(groups)){
  message(names(groups)[i])
  group <- names(groups)[i]
  mean <- mean(score$PRS_0_1[score$ADM %in% groups[[i]]], na.rm = T)
  sd <- sd(score$PRS_0_1[score$ADM %in% groups[[i]]], na.rm = T)
  
  test.out <- data.frame(group, mean, sd)
  for(j in 1:length(groups)){
    vs_group <- names(groups)[j]
    test.out[, paste0("pvalue_vs_", vs_group)] <- t.test(score$PRS_0_1[score$ADM %in% groups[[i]]],
                                                         score$PRS_0_1[score$ADM %in% groups[[j]]])$p.value
  }
  test.out$test <- "direct_comparison"
  test.out$cohort <- "NFLD"
  table.out <- rbind(table.out, test.out)
  
}


for(i in 1:length(groups)){
  message(names(groups)[i])
  group <- names(groups)[i]
  mean <- mean(score$overtransmission[score$ADM %in% groups[[i]]], na.rm = T)
  sd <- sd(score$overtransmission[score$ADM %in% groups[[i]]], na.rm = T)
  
  test.out <- data.frame(group, mean, sd)
  for(j in 1:length(groups)){
    vs_group <- names(groups)[j]
    test.out[, paste0("pvalue_vs_", vs_group)] <- t.test(score$overtransmission[score$ADM %in% groups[[i]]],
                                                         score$overtransmission[score$ADM %in% groups[[j]]])$p.value
  }
  test.out$test <- "parentalDeviationComparison"
  test.out$cohort <- "NFLD"
  table.out <- rbind(table.out, test.out)
}

write.table(table.out, "NFLD_SSC_ADM_PRS_comparison.tsv", sep="\t",row.names=F, quote=F, col.names=T)

# t.test(score$PRS_0_1[score$ADM == "Dysmorphic"],
#        score$PRS_0_1[score$ADM == "Nondysmorphic"]) #p = 1
# t.test(score$overtransmission[score$ADM == "Dysmorphic"],
#        score$overtransmission[score$ADM == "Nondysmorphic"]) #p = 0.97
# 
score <- score[!is.na(score$ADM), ]
t.dys <- t.test(score$overtransmission[score$ADM == "Dysmorphic"], alternative = "greater", mu = 0)
t.nondys <- t.test(score$overtransmission[score$ADM == "Nondysmorphic"], alternative = "greater", mu = 0)
t.fn <- t.test(score$overtransmission[score$ADM == "Nondysmorphic-False negative"], alternative = "greater", mu = 0)

score$Relation[score$ADM == "Dysmorphic"] <- sprintf("ADM\nDysmorphic\n(n=%s)",
                                                         sum(score$ADM == "Dysmorphic"))
score$Relation[score$ADM == "Nondysmorphic"] <- sprintf("ADM\nNondysmorphic\n(n=%s)", 
                                                            sum(score$ADM == "Nondysmorphic"))
score$Relation[score$ADM == "Nondysmorphic-False negative"] <- sprintf("Nondysmorphic-False negative\n(n=%s)", 
                                                        sum(score$ADM == "Nondysmorphic-False negative"))
score <- score[!score$ADM %in% c("proband", "Nondysmorphic-False negative", ""), ]
score$Relation <- factor(score$Relation, levels = unique(score$Relation)[c(1,2)])

# score$Relation[score$Relation == "proband"] <- "removed FN\n(n=16)"
write.table(score, "nfld.score.tsv", sep="\t", row.names=F, quote=T, col.names=T)

ggplot(score, 
       aes(x = Relation, y = overtransmission, fill = Relation)) + geom_violin(color = NA) +
  geom_boxplot(fill = "white", width = .2, alpha=.6) + theme_classic() +     
  ylab("PRS difference from the parents average") + geom_hline(yintercept = 0, lty = 2) +
  annotate("text", x = c(1,2), y = c(4),
           label = paste("p=",c(signif(t.nondys$p.value, digits = 2),
                                signif(t.dys$p.value, digits = 2)), sep="")) + xlab("") +
  scale_fill_manual(values = c("#009eed", "#F28602")) + theme(legend.position = "none") +
  #scale_color_manual(values = c("grey", "#009eed", "#F28602")) + 
  ggtitle("NFLD")
ggsave("nfld.adm.prs.pdf", width = 4, height = 4)

p1 <- ggplot(score,
       aes(x = Relation, y = overtransmission, color = Relation)) +
  geom_boxplot(outlier.alpha = 0, width = .4, show.legend = F, fill = "white", alpha = .6) +
  geom_point(aes(fill = Relation),
             position = position_jitterdodge(jitter.width = .1, dodge.width = .4),
             color = "black",
              shape = 21, width = .1, show.legend = F, alpha = .7) + theme_bw() +
  ylab("PRS difference from the parents average") + geom_hline(yintercept = 0, lty = 2) +
  annotate("text", x = c(1,2), y = c(3),
           label = paste("p=",c(signif(t.dys$p.value, digits = 2),
                                signif(t.nondys$p.value, digits = 2)), sep="")) + xlab("") +
  scale_color_manual(values = c("red", "grey", "black"), guide = guide_legend()) +
  scale_fill_manual(values = c("complex" = "red", "essential" = "grey", "equivocal" = "orange")) +
  ggtitle("NFLD") + coord_cartesian(ylim = c(-3, 3))


# cowplot::plot_grid(p1, p2, nrow = 1)
# ggsave("PRS.ADM.color.pdf", width = 7.5, height = 5)

# ggplot(score[score$Relation != "proband", ], 
#        aes(x = Dysmorphology.classification, y = overtransmission, color = Dysmorphology.classification)) + 
#   geom_boxplot(outlier.alpha = 0, show.legend = F, width = .4, fill = "white", alpha = .6) + 
#   geom_jitter(shape = 21, width = .1, show.legend = F, alpha = .3) + theme_bw() +     
#   ylab("PRS difference from the parents average") + geom_hline(yintercept = 0, lty = 2) + xlab("") +
#   ggtitle("NFLD")

#### check sex bias
nfld.aff <- score[score$incl.aff.in.my.study == 1, ]
nfld.aff$Sex..CRV <- ifelse(nfld.aff$Sex..CRV == "M", "male", "female")

table(ssc.aff$Sex)
table(nfld.aff$Sex)

lm.ref <- lm(SCORESUM ~ adm.dysmorphic, ssc.aff)
lm.add <- lm(SCORESUM ~ adm.dysmorphic + Sex, ssc.aff)
anova(lm.ref, lm.add, test="Chisq")

lm.ref <- lm(PRS_0_1 ~ ADM, nfld.aff)
lm.add <- lm(PRS_0_1 ~ ADM + Sex..CRV, nfld.aff)
anova(lm.ref, lm.add, test="Chisq")

se1 <- ggplot(ssc.aff, aes(x = Sex, y = SCORESUM)) + geom_boxplot() + geom_jitter() + theme_bw() + xlab("Sex") + ylab("PRS") +
  ggtitle("SSC affected, p=0.97")
se2 <- ggplot(nfld.aff, aes(x = Sex..CRV, y = nfld.aff$PRS_0_1)) + geom_boxplot() + geom_jitter() + theme_bw() + 
  ggtitle("NFLD affected, p=0.09") +
  xlab("Sex") + ylab("PRS")
  
cowplot::plot_grid(se1, se2)
ggsave("PRS_Sex.png", height = 6, width = 6)
