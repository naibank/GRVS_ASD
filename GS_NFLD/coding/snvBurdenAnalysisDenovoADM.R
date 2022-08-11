### set working directory 237 subjects
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### install my burden analysis package
library(devtools)
#install_github("naibank/GSBurden")

### clear existing variables
rm(list=ls())

### load my package
library(GSBurden)

### read snv gene set
dt <- read.delim("denovo.gs.matrix.tsv", stringsAsFactors = F)

### read phenotype table
info <- read.delim("../../../../main/data/2021.02.01-updated NFLD phenotype table IQ ADM and 3 status.txt", stringsAsFactors = F)
# info <- read.delim("../data/2020.11.23-updated NFLD phenotype table IQ and ADM_bak_use_2021_file_instead.txt", stringsAsFactors = F)
# info$Dysmorphology.classification <- tolower(info$Dysmorphology.classification)
info <- info[!is.na(info$ADM), ]

### get only columns that I use
info <- info[, c("WGS_ManuID", "DNA.source", "Relation", "Sex..CRV", "incl.aff.in.my.study", "incld.cohort", "WGS.status", "Paternal.age",
                 "remove.family.for.rare.burden", "Platform", "Dysmorphology.classification", "Grouped_Dysmorphology", "Final.Dysmorphology.scores", "ADM")]
info <- info[info$incl.aff.in.my.study == 1 & info$WGS.status == "WGS-ed, complete fam", ]
info$Sample <- gsub("-", "_", info$WGS_ManuID)

sample.not.in <- info$Sample[!info$Sample %in% dt$Sample]
sample.not.in <- data.frame("Sample" = sample.not.in)
sample.not.in[names(dt)[-1]] <- 0


dt <- rbind(dt, sample.not.in)
### read estimated PCA
pca <- read.delim("../../../../NFLD2019/estimated.pca.12March2019.tsv", stringsAsFactors = F)
pca <- pca[!pca$sample %in% c("3-0456-000", "3-0458-000"), ]
pca$sample <- gsub("A", "", pca$sample)

### merge PCA to phenotype table
info <- merge(info, pca[, c("sample", "pc1", "pc2", "pc3")], by.x = "WGS_ManuID", by.y = "sample", all.x =T)
info$Sample <- gsub("-", "_", info$WGS_ManuID)

### limit to proband with following criterions and merge CNVs table and phenotype table
merge <- merge(dt, info, by = "Sample", all.x = T)


### load geneset
load("../../new_set/gsNFLD_2022.RData")

snv.matrix <- merge
snv.matrix$ADM <- factor(snv.matrix$ADM, levels = c("Nondysmorphic", "Dysmorphic"))

names(snv.matrix) <- gsub("Total_", "gene_count_", names(snv.matrix))
snv.matrix$Paternal.age <- as.numeric(snv.matrix$Paternal.age)

### separate gene-sets into 3 collections for multiple testing correction purpose
brainExp.gs <- grep("Expr|blueModule|Ilmn_BM|Bspan|FMR1|Synapse|Neurof|LoF|etal", names(gsNFLD))
phenotype.gs <- grep("Ph", names(gsNFLD))

### define covariates
covariates <- c("Sex..CRV", "Platform", "pc1", "pc2", "pc3")

### run burden test for each of 3 gene-set collections
perm.values <- data.frame()
nperm <- 1000
test <- "ADM"

snv.matrix[, test] <- factor(snv.matrix[, test])
brainExp.test.out <- CNVBurdenTest(snv.matrix[!is.na(snv.matrix[, test]),], gsNFLD[brainExp.gs], test, covariates, nperm = nperm)
phenotype.test.out <- CNVBurdenTest(snv.matrix[!is.na(snv.matrix[, test]),], gsNFLD[phenotype.gs], test, covariates, nperm = nperm)
# neurofunctions.test.out <- CNVBurdenTest(snv.matrix[!is.na(snv.matrix[, test]),], gsNFLD[neurofunctions.gs], test, covariates, nperm = 1000)
coeff.out <- CNVBurdenTest(snv.matrix[!is.na(snv.matrix[, test]),], gsNFLD, test, covariates,
                           permutation = F, correctGlobalBurden = F)$Test

###mpc test
ref.term <- sprintf("%s ~ %s + gene_count_tier1_ms", test, paste(covariates, collapse = " + "))
add.term <- sprintf("%s + %s", ref.term, "MisMPC_morethan2")
coeff.term <- sprintf("%s ~ %s + MisMPC_morethan2", test, paste(covariates, collapse = " + "))
  
  
ref.model <- glm(ref.term, snv.matrix, family = binomial(link = "logit"))
add.model <- glm(add.term, snv.matrix, family = binomial(link = "logit"))
coeff.model <- glm(coeff.term, snv.matrix, family = binomial(link = "logit"))
ano <- anova(ref.model, add.model, test = "Chisq")
  
names(ano)[length(names(ano))] <- "pvalue"
pvalue <- ano$pvalue[2]
coefficient <- add.model$coefficients["MisMPC_morethan2"]
sm <- summary(add.model)
  
if("MisMPC_morethan2" %in% rownames(sm$coefficients)){
  waldp <- sm$coefficients["MisMPC_morethan2", 4]
  stderr <- sm$coefficients["MisMPC_morethan2", 2]
  testval <- sm$coefficients["MisMPC_morethan2", 3]
}else{
  waldp <- stderr <- testval <- NA
}

conf <- confint(add.model)
conf.coeff <- confint(coeff.model)

perm.test.pvalues <- data.frame()
for(i in 1:nperm){
  snv.matrix$permlab <- sample(snv.matrix[, test])
  ref.term <- sprintf("permlab ~ %s + gene_count_tier1_ms", paste(covariates, collapse = " + "))
  add.term <- sprintf("%s + %s", ref.term, "MisMPC_morethan2")
  
  ref.model <- glm(ref.term, snv.matrix, family = binomial(link = "logit"))
  add.model <- glm(add.term, snv.matrix, family = binomial(link = "logit"))
  ano.perm <- anova(ref.model, add.model, test = "Chisq")
  
  names(ano.perm)[length(names(ano.perm))] <- "pvalue"
  perm.coeff <- add.model$coefficients["MisMPC_morethan2"]
  perm.test.pvalues <- rbind(perm.test.pvalues, data.frame("cnvtype" = "tier1_ms", "pvalue" = ano.perm$pvalue[2], "coeff" = perm.coeff))
}

coeff.out <- rbind(coeff.out, data.frame("geneset" = "MisMPC_morethan2", "type" = "tier1_ms", 
                                         "coefficient" = coeff.model$coefficients["MisMPC_morethan2"], 
                                         "coeff.upper" = conf.coeff["MisMPC_morethan2", 2],
                                         "coeff.lower" = conf.coeff["MisMPC_morethan2", 1], 
                                         "pvalue" = pvalue, stderr, testval, waldp, "permFDR" = NA))

#########################################

burden.test.out <- rbind(brainExp.test.out$Test, phenotype.test.out$Test)
burden.test.out <- rbind(burden.test.out, data.frame("geneset" = "MisMPC_morethan2", "type" = "tier1_ms", 
                                                     "coefficient" = coefficient, "coeff.upper" = conf["MisMPC_morethan2", 2],
                                                     "coeff.lower" = conf["MisMPC_morethan2", 2], 
                                                     "pvalue" = pvalue, stderr, testval, waldp, "permFDR" = NA))

burden.test.out <- merge(burden.test.out[, -c(3:5)], coeff.out[, 3:5], by = "row.names", all.x = T)
burden.test.out <- burden.test.out[, c(2:3, 9:11, 4:8)]

perm.temp <- rbind(brainExp.test.out$Permutation.Test, phenotype.test.out$Permutation.Test, perm.test.pvalues)
perm.temp$geneset <- gsub("[0-9]", "", rownames(perm.temp))
perm.values <- rbind(perm.values, perm.temp)

getPermPvalue <- function(target.pvalue, actual.pvalues, perm.pvalues){
  actual.times <- sum(actual.pvalues <= target.pvalue, na.rm = T)
  perm.times <- sum(perm.pvalues <= target.pvalue, na.rm = T)
  permFDR <- (perm.times/length(perm.pvalues))/(actual.times/length(actual.pvalues))
  permFDR <- ifelse(permFDR > 1, 1, permFDR)
  return(permFDR)
}

getPermFDR <- function(i, dt){
  pvalue <- dt$pvalue[i]
  return(min(dt$permFDR[dt$pvalue >= pvalue], na.rm = T))
}


perm.pvalues <- perm.values
dt <- burden.test.out

brainExp.test <- dt[dt$geneset %in% c("MisMPC_morethan2", names(gsNFLD)[brainExp.gs]), ]
phenotype.test <- dt[dt$geneset %in% names(gsNFLD)[phenotype.gs], ]

perm.pvalues$geneset <- gsub("_lof|_tier1_missense|_tier2_missense|_tier_ms|_all_tier_ms", "", perm.pvalues$geneset)
brain.perm <- perm.pvalues[perm.pvalues$geneset %in% c("MisMPC_morethan2", "MisMPC_morethan", names(gsNFLD)[brainExp.gs]), ] 
phenotype.perm <- perm.pvalues[perm.pvalues$geneset %in% names(gsNFLD)[phenotype.gs], ] 

brainExp.test$permFDR <- sapply(brainExp.test$pvalue, getPermPvalue, brainExp.test$pvalue, brain.perm$pvalue)
phenotype.test$permFDR <- sapply(phenotype.test$pvalue, getPermPvalue, phenotype.test$pvalue, phenotype.perm$pvalue)

brainExp.test$permFDR <- sapply(1:nrow(brainExp.test), getPermFDR, brainExp.test)
phenotype.test$permFDR <- sapply(1:nrow(phenotype.test), getPermFDR, phenotype.test)


dt.out <- rbind(brainExp.test, phenotype.test)
dt.out$coeff.lower[is.na(dt.out$coeff.lower) & !is.na(dt.out$coeff.upper)] <- -Inf
dt.out$coeff.upper[!is.na(dt.out$coeff.lower) & is.na(dt.out$coeff.upper)] <- Inf
  
write.table(dt.out, "burdenResult/snv.adm.denovo.gsburden.withPerm.tsv", sep="\t", row.names=F, quote=F, col.names=T)


### prepare data for plot
lof <- data.frame()
ms1 <- data.frame()
ms2 <- data.frame()
burden.test.out <- read.delim("burdenResult/snv.adm.denovo.gsburden.withPerm.tsv", stringsAsFactors = F)
burden.test.out$Test <- test
burden.test.out$FDRRange <- "NotSignificant"
burden.test.out$FDRRange[burden.test.out$permFDR < 0.20] <- "< 20%"
burden.test.out$FDRRange[burden.test.out$permFDR < 0.10] <- "< 10%"
burden.test.out$FDRRange[burden.test.out$permFDR < 0.05] <- "< 5%"
burden.test.out$FDRRange <- factor(burden.test.out$FDRRange, levels = c("NotSignificant", "< 20%", "< 10%", "< 5%"))

lof <- rbind(lof, burden.test.out[burden.test.out$type == "lof", ])
ms1 <- rbind(ms1, burden.test.out[burden.test.out$type == "tier1_ms", ])
ms2 <- rbind(ms2, burden.test.out[burden.test.out$type == "tier2_ms", ])


lof.significant <- unique(lof$geneset[lof$permFDR < 0.1])
ms1.significant <- unique(ms1$geneset[ms1$permFDR < 0.1])
ms2.significant <- unique(ms2$geneset[ms2$permFDR < 0.1])
### plot deletion result
lof$Test <- ifelse(lof$Test == "MultiClass", "Three subtypes", "Two subtypes")
lof$Test <- factor(lof$Test, levels = c("Two subtypes", "Three subtypes"))
lof$coeff.upper[is.na(lof$coeff.upper)] <- Inf
lof$coeff.lower[is.na(lof$coeff.lower)] <- -Inf
library(ggplot2)
ggplot(lof, aes(x = geneset, y = coefficient, fill = FDRRange)) +
  
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = coeff.lower, ymax = coeff.upper, group = Test), 
                position = position_dodge(width = 0.7), width = 0) +
  geom_point(aes(shape = Test), size = 3, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("NotSignificant" = "white",
                               "< 20%" = "yellow",
                               "< 10%" = "orange",
                               "< 5%" = "red")) +
  scale_shape_manual(values = c(21,23)) +
  theme_bw() + coord_cartesian(ylim=c(-0.5, 3)) +
  theme(axis.text.x = element_text(angle = 90, hjust =1)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) + ggtitle("LoF variants")

ggsave("burdenResult/gs.adm.denovo.lof.burden.pdf", width = 12, height = 5)

### plot ms1 result
ms1$Test <- ifelse(ms1$Test == "MultiClass", "Three subtypes", "Two subtypes")
ms1$Test <- factor(ms1$Test, levels = c("Two subtypes", "Three subtypes"))

ms1$coeff.upper[is.na(ms1$coeff.upper)] <- Inf
ms1$coeff.lower[is.na(ms1$coeff.lower)] <- -Inf

ggplot(ms1, aes(x = geneset, y = coefficient, fill = FDRRange)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = coeff.lower, ymax = coeff.upper, group = Test), 
                position = position_dodge(width = 0.7), width = 0) +
  geom_point(aes(shape = Test), size = 3, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("NotSignificant" = "white",
                               "< 20%" = "yellow",
                               "< 10%" = "orange",
                               "< 5%" = "red")) +
  scale_shape_manual(values = c(21,23)) +
  theme_bw() + coord_cartesian(ylim=c(-0.5, 3)) +
  theme(axis.text.x = element_text(angle = 90, hjust =1)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) + ggtitle("Tier 1 missenses")


ggsave("burdenResult/gs.adm.denovo.ms1.burden.pdf", width = 12, height = 5)


### plot ms2 result
ms2$Test <- ifelse(ms2$Test == "MultiClass", "Three subtypes", "Two subtypes")
ms2$Test <- factor(ms2$Test, levels = c("Two subtypes", "Three subtypes"))

ms2$coeff.upper[is.na(ms2$coeff.upper)] <- Inf
ms2$coeff.lower[is.na(ms2$coeff.lower)] <- -Inf

ggplot(ms2, aes(x = geneset, y = coefficient, fill = FDRRange)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = coeff.lower, ymax = coeff.upper, group = Test), 
                position = position_dodge(width = 0.7), width = 0) +
  geom_point(aes(shape = Test), size = 3, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("NotSignificant" = "white",
                               "< 20%" = "yellow",
                               "< 10%" = "orange",
                               "< 5%" = "red")) +
  scale_shape_manual(values = c(21,23)) +
  theme_bw() +coord_cartesian(ylim=c(-3, 3)) +
  theme(axis.text.x = element_text(angle = 90, hjust =1)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) + ggtitle("Tier 2 missenses")


ggsave("burdenResult/gs.adm.denovo.ms2.burden.pdf", width = 12, height = 5)

write.table(lof, "burdenResult/gs.adm.lof.denovo.result.table.tsv", sep="\t", row.names=F, quote=F, col.names=T)
write.table(ms1, "burdenResult/gs.adm.ms1.denovo.result.table.tsv", sep="\t", row.names=F, quote=F, col.names=T)
write.table(ms2, "burdenResult/gs.adm.ms2.denovo.result.table.tsv", sep="\t", row.names=F, quote=F, col.names=T)

####################################
######## global burden test ########
# snv.matrix$gene_count_lof <- snv.matrix$TotalRare
# snv.matrix$gene_count_tier1_ms <- snv.matrix$TotalRare
# snv.matrix$gene_count_tier2_ms <- snv.matrix$TotalRare
snv.matrix$cnv_count_lof <- snv.matrix$TotalRare
snv.matrix$cnv_count_tier1_ms <- snv.matrix$TotalRare
snv.matrix$cnv_count_tier2_ms <- snv.matrix$TotalRare
snv.matrix$cnv_size_lof <- snv.matrix$TotalRare
snv.matrix$cnv_size_tier1_ms <- snv.matrix$TotalRare
snv.matrix$cnv_size_tier2_ms <- snv.matrix$TotalRare

global.test <- data.frame()
covariates <- c("Sex..CRV", "Platform", "Paternal.age", "pc1", "pc2", "pc3")

tmp.test <- CNVGlobalTest(snv.matrix, test, covariates)
tmp.test$Test <- test
global.test <- rbind(global.test, tmp.test)

global.test <- global.test[global.test$global == "gene_count", ]
global.test$global <- "variants"
global.test$PvalueRange <- "NotSignificant"
global.test$PvalueRange[global.test$pvalue < 0.1] <- "< 0.1"
global.test$PvalueRange[global.test$pvalue < 0.05]  <- "< 0.05"
global.test$PvalueRange[global.test$pvalue < 0.01]  <- "< 0.01"
global.test$PvalueRange <- factor(global.test$PvalueRange, levels = c("NotSignificant", "< 0.1", "< 0.05", "< 0.01"))
global.test <- global.test[global.test$Test %in% c("GroupedDysmorphology", "MultiClass"), ]
global.test$Test <- ifelse(global.test$Test == "GroupedDysmorphology", "Two subtypes", "Three subtypes")
global.test$Test <- factor(global.test$Test, levels = c("Two subtypes", "Three subtypes"))
ggplot(global.test, aes(x = type, y = coefficient, fill = PvalueRange)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = lowerbound, ymax = upperbound, group = Test), 
                position = position_dodge(width = 0.7), width = 0) +
  geom_point(aes(shape = Test), size = 3, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("NotSignificant" = "white",
                               "< 0.1" = "yellow",
                               "< 0.05" = "orange",
                               "< 0.01" = "red")) +
  scale_shape_manual(values = c(21,23)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust =1)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) + ggtitle("Global Burden")

write.table(global.test, "burdenResult/global.adm.denovo.burden.tsv", sep="\t", row.names=F, quote=F, col.names=T)
ggsave("burdenResult/global.adm.denovo.burden.pdf", width = 6, height = 4)
