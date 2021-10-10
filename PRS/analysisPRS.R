setwd("~/Documents/doc/working/NFLD/2019/common/")
rm(list=ls())
score <- read.delim("NFLD.PRS.Jan2020.tsv", stringsAsFactors = F)

#names(score)[c(26:32)] <- c("PRS_0_1", "PRS_0_2", "PRS_0_3", "PRS_0_4", "PRS_0_5", "PRS_1", "var1") 
write.table(score, "NFLD.PRS.Jan2020.tsv", sep="\t", row.names=F, quote=F, col.names=T)

library(ggplot2)
score.proband <- score[score$Platform != "Complete Genomics" & score$Dysmorphology.classification != "na-parent", ]
ggplot(score.proband, 
       aes(x = Dysmorphology.classification, y = PRS_0_1, color = Dysmorphology.classification)) + 
  geom_boxplot(outlier.alpha = 0, show.legend = F) + geom_jitter(shape = 21, show.legend = F) + theme_bw() +
  ylab("iPSYCH-derived PRS") + 
ggsave("PRS.dysmorphology.pdf", width = 8, height = 7)

t.test(score.proband$PRS_0_1[score.proband$Dysmorphology.classification == "complex"],
       score.proband$PRS_0_1[score.proband$Dysmorphology.classification == "essential"]) #p = 0.2031
t.test(score.proband$PRS_0_1[score.proband$Dysmorphology.classification == "equivocal"],
       score.proband$PRS_0_1[score.proband$Dysmorphology.classification == "essential"]) #p = 0.6285
t.test(score.proband$PRS_0_1[score.proband$Dysmorphology.classification == "complex"],
       score.proband$PRS_0_1[score.proband$Dysmorphology.classification == "equivocal"]) #p = 0.1538

ggplot(score, 
       aes(x = Platform, y = PRS_0_1, color = Platform)) + 
  geom_boxplot(outlier.alpha = 0, show.legend = F) + geom_jitter(shape = 21, show.legend = F) + theme_bw() +
  ylab("iPSYCH-derived PRS") +
ggsave("PRS.different.platforms.pdf", width = 8, height = 7)

summary(aov(PRS_0_1 ~ Dysmorphology.classification, 
            data = score[score$Platform != "Complete Genomics" & score$Dysmorphology.classification != "na-parent", ]))
summary(aov(PRS_0_1 ~ Platform, data = score))
summary(aov(PRS_0_1 ~ Platform, data = score[score$Platform != "Complete Genomics", ]))

#### calculate mid parent
score$ParentMean <- NA
for(i in 1:nrow(score)){
  tmp <- score[i, ]
  if(!tmp$Relation %in% c("Father", "Mother")){
    famid <- tmp$WGS_FamilyID
    parent <- score[score$WGS_FamilyID == famid & score$Relation %in% c("Father", "Mother"), ]
    
    if(nrow(parent) == 2){
      average <- mean(parent$PRS_0_1)
    }
    
    score$ParentMean[i] <- average
  }
}

score$overtransmission <- score$PRS_0_1-score$ParentMean
score.proband <- score[score$Platform != "Complete Genomics" & score$Dysmorphology.classification != "na-parent", ]
  
ggplot(score.proband, 
       aes(x = Dysmorphology.classification, y = overtransmission, color = Dysmorphology.classification)) + 
  geom_boxplot(outlier.alpha = 0, show.legend = F) + geom_jitter(shape = 21, show.legend = F) + theme_bw() +
  ylab("PRS difference from the parents average") + geom_hline(yintercept = 0, lty = 2) +
  annotate("text", x = c(1,2,3), y = c(0.00025),
           label = c("p=0.78", "p=0.07", "p=5e-4"))
ggsave("PRS.dysmorphology.transmission.eq.test.pdf", width = 8, height = 7)

summary(aov(PRS_0_1 ~ overtransmission, 
            data = score[score$Platform != "Complete Genomics" & score$Dysmorphology.classification != "na-parent", ]))


t.test(score.proband$overtransmission[score.proband$Dysmorphology.classification == "complex"],
       mu = 0) #p = 0.7775
t.test(score.proband$overtransmission[score.proband$Dysmorphology.classification == "equivocal"],
       mu = 0) #p = 0.07
t.test(score.proband$overtransmission[score.proband$Dysmorphology.classification == "essential"],
       mu = 0) #p = 0.0005719


########## PCA ##########
score.tmp <- readLines("cluster_NFLD.mds")
for(i in 1:10)
  score.tmp <- gsub("  ", " ", score.tmp)
writeLines(score.tmp, "NFLD.clusters.tsv")

pca <- read.delim("NFLD.clusters.tsv", stringsAsFactors = F, sep=" ")[, -c(1, 8)]
info <- read.delim("../data/2019.03.08-updated NFLD phenotype table full.txt", stringsAsFactors = F)

pca <- merge(pca, info, by.x = "IID", by.y = "WGS_ManuID", all.x = T)

library(ggplot2)
ggplot(pca[!is.na(pca$Platform), ], aes(x = C1, y = C2, color = Platform)) + geom_point()
ggsave("pca.common.variant.by.platform.pdf", width = 6, height = 6)
