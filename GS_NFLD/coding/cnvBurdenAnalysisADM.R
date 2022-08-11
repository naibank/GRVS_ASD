### set working directory 301 subjects
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### install my burden analysis package
# library(devtools)
# install_github("naibank/GSBurden")

### clear existing variables
rm(list=ls())

library(ggplot2)
### read phenotype table
info <- read.delim("../../../../main/data/2021.02.01-updated NFLD phenotype table IQ ADM and 3 status.txt", stringsAsFactors = F)
info <- info[!is.na(info$ADM), ]

### get only columns that I use
info <- info[, c("WGS_ManuID", "DNA.source", "Relation", "Sex..CRV", "incl.aff.in.my.study", "incld.cohort", 
                 "Platform", "ADM")]

# unique(merge$sample[!merge$sample %in% info$WGS_ManuID])
# [1] "3-0741-000A" "3-0179-000A" "3-0450-000A" "3-0213-000A" "3-0730-000A"

### rename ID of re-sequenced samples, if I understand it correctly
info$WGS_ManuID <- gsub("3-0741-000", "3-0741-000A", info$WGS_ManuID)
info$WGS_ManuID <- gsub("3-0179-000", "3-0179-000A", info$WGS_ManuID)
info$WGS_ManuID <- gsub("3-0450-000", "3-0450-000A", info$WGS_ManuID)
info$WGS_ManuID <- gsub("3-0213-000", "3-0213-000A", info$WGS_ManuID)
info$WGS_ManuID <- gsub("3-0730-000", "3-0730-000A", info$WGS_ManuID)

### read estimated PCA
pca <- read.delim("../../../../NFLD2019/estimated.pca.12March2019.tsv", stringsAsFactors = F)
pca <- pca[!pca$sample %in% c("3-0456-000", "3-0458-000"), ]
pca$sample <- gsub("A", "", pca$sample)

### merge PCA to phenotype table
info <- merge(info, pca[, c("sample", "pc1", "pc2", "pc3")], by.x = "WGS_ManuID", by.y = "sample", all.x =T)
#0466 -> no population stratification done for this sample.

### limit to proband with following criterions and merge CNVs table and phenotype table
info <- info[info$incl.aff.in.my.study == 1, ]
# merge <- merge[merge$sample %in% info$WGS_ManuID, ]

### load my package
library(GSBurden)

### load geneset
load("../../new_set/gsNFLD_2022.RData")

### limit to CNVs smaller than 3MB
cnvs <- read.delim("all.cnvs.June2022.txt", stringsAsFactors = F)
# cnvs <- cnvs[which(cnvs$size > 10000), ]
### read coding exons table
gene.in <- na.omit(read.delim("../../../../main/requiredData/hg19_refGene_28_04_17.cds.txt", header = F, stringsAsFactors = F))
names(gene.in) <- c("chr", "start", "end", "isoid", "genesymbol", "entrezid")
genes <- GeneAnnotation(gene.in$entrezid, gene.in$chr, gene.in$start, gene.in$end, gene.in$genesymbol)

# appris <- read.delim("../../../ReferenceData/APPRIS/appris_data.principal.hg19.txt",
#                      stringsAsFactors = F, header = F, col.names = c("gsymbol", "enzid", "isoid", "x1", "tier"))
# appris$isoid <- sapply(sapply(appris$isoid, strsplit, "\\."), "[", 1)
# appris <- appris[grep("PRINCIPAL", appris$tier), ]
# 
# trscrpt <- read.delim("../../../ReferenceData/RefGene_fromJohn/hg19_refGene_28_04_17.transcript.txt", 
#                       stringsAsFactors = F, header = F,
#                       col.names = c("chr", "start", "end", "isoid", "gsymbol", "enzid"))
# appris <- trscrpt[which(trscrpt$isoid %in% appris$isoid), ]

### get gene-set matrix
cnvs$type <- ifelse(cnvs$type == "deletion", "DEL", "DUP")
cnv.matrix <- getCNVGSMatrix(cnvs, genes, gsNFLD)

# cnv.matrix <- getCNVGSMatrix(cnvs[cnvs$type == "deletion", ], genes, gsNFLD)
# cnv.dup.matrix <- getCNVDupGSMatrix(cnvs[cnvs$type == "duplication", ], genes, appris, gsNFLD)
# names(cnv.matrix)[-1] <- paste(names(cnv.matrix)[-1], "DEL", sep="_")
# 
# cnv.matrix <- merge(cnv.matrix, cnv.dup.matrix, by = "sample", all = T)
# cnv.matrix[is.na(cnv.matrix)] <- 0

### merge gene-set matrix with phenotype table
cnv.matrix <- merge(cnv.matrix, info[, c("WGS_ManuID", "Sex..CRV", "Platform", "ADM", "pc1", "pc2", "pc3")],
                    by.x = "sample", by.y = "WGS_ManuID", all.x = T)

write.table(cnv.matrix, "ADM.cnv.coding.matrix.tsv", sep="\t", row.names=F, quote=F, col.names=T)

### separate gene-sets into 3 collections for multiple testing correction purpose
brainExp.gs <- grep("Expr|blueModule|Ilmn_BM|Bspan|FMR1|Synapse|Neurof|LoF|etal", names(gsNFLD))
phenotype.gs <- grep("Ph", names(gsNFLD))

### define covariates
covariates <- c("Sex..CRV", "Platform", "pc1", "pc2", "pc3")

### run burden test for each of 3 gene-set collections
perm.values <- data.frame()

test <- "ADM"
cnv.matrix$ADM <- factor(cnv.matrix$ADM, levels = c("Nondysmorphic", "Dysmorphic"))
cnv.matrix <- cnv.matrix[which(cnv.matrix$Platform != "Complete Genomics"), ]
brainExp.test.out <- CNVBurdenTest(cnv.matrix[!is.na(cnv.matrix[, test]),], gsNFLD[brainExp.gs], test, covariates, nperm = 100)
phenotype.test.out <- CNVBurdenTest(cnv.matrix[!is.na(cnv.matrix[, test]),], gsNFLD[phenotype.gs], test, covariates, nperm = 100)

burden.test.out <- rbind(brainExp.test.out$Test, phenotype.test.out$Test)
perm.temp <- rbind(brainExp.test.out$Permutation.Test, phenotype.test.out$Permutation.Test)
perm.temp$geneset <- gsub("[0-9]", "", rownames(perm.temp))
perm.values <- rbind(perm.values, perm.temp)
perm.values$geneset <- perm.temp$geneset
### cal perm

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



dt <- burden.test.out
brainExp.test <- dt[dt$geneset %in% names(gsNFLD)[brainExp.gs], ]
phenotype.test <- dt[dt$geneset %in% names(gsNFLD)[phenotype.gs], ]

perm.pvalues <- perm.values
perm.pvalues$geneset <- gsub("_DEL|_FullDup|_PartialDup|_DUP", "", perm.pvalues$geneset)
brain.perm <- perm.pvalues[perm.pvalues$geneset %in% names(gsNFLD)[brainExp.gs], ] 
phenotype.perm <- perm.pvalues[perm.pvalues$geneset %in% names(gsNFLD)[phenotype.gs], ] 

brainExp.test$permFDR <- sapply(brainExp.test$pvalue, getPermPvalue, brainExp.test$pvalue, brain.perm$pvalue)
phenotype.test$permFDR <- sapply(phenotype.test$pvalue, getPermPvalue, phenotype.test$pvalue, phenotype.perm$pvalue)

brainExp.test$permFDR <- sapply(1:nrow(brainExp.test), getPermFDR, brainExp.test)
phenotype.test$permFDR <- sapply(1:nrow(phenotype.test), getPermFDR, phenotype.test)


dt.out <- rbind(brainExp.test, phenotype.test)

### prepare data for plot
del <- data.frame()
dup <- data.frame()

burden.test.out <- dt.out
burden.test.out$Test <- test
burden.test.out$FDRRange <- "NotSignificant"
burden.test.out$FDRRange[burden.test.out$permFDR < 0.20] <- "< 20%"
burden.test.out$FDRRange[burden.test.out$permFDR < 0.10] <- "< 10%"
burden.test.out$FDRRange[burden.test.out$permFDR < 0.05] <- "< 5%"
burden.test.out$FDRRange <- factor(burden.test.out$FDRRange, levels = c("NotSignificant", "< 20%", "< 10%", "< 5%"))

burden.test.out$type <- factor(burden.test.out$type, levels = c("DEL", "DUP"))
burden.test.out$geneset <- factor(burden.test.out$geneset, levels = names(gsNFLD)[names(gsNFLD) %in% burden.test.out$geneset])

####plot
ggplot(burden.test.out, aes(x = geneset, y = coefficient, fill = FDRRange)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = coeff.lower, ymax = coeff.upper, group = Test),
                position = position_dodge(width = 0.5), width = 0) +
  geom_point(shape = 21, size = 3, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("NotSignificant" = "white",
                               "< 20%" = "yellow",
                               "< 10%" = "orange",
                               "< 5%" = "red")) +
  theme_bw() + facet_wrap(.~type, nrow = 3) + coord_cartesian(ylim = c(-2, 2)) +
  theme(axis.text.x = element_text(angle = 90, hjust =1)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) + ggtitle("Deletions")

ggsave("burdenResult/ADM.NFLD.cnv.gs.burden.png", width=12, height = 7)
write.table(burden.test.out, "burdenResult/ADM.NFLD.cnv.gs.burden.tsv", sep="\t", row.names=F, quote=F, col.names=T)

# del <- rbind(del, burden.test.out[burden.test.out$type == "deletion", ])
# dup <- rbind(dup, burden.test.out[burden.test.out$type == "duplication", ])
# 
# del.significant <- unique(del$geneset[del$permFDR < 0.2])
# dup.significant <- unique(dup$geneset[dup$permFDR < 0.2])
# ### plot deletion result
# del$Test <- factor(del$Test)
# library(ggplot2)
# ggplot(del, aes(x = geneset, y = coefficient, fill = FDRRange)) +
#   
#   geom_hline(yintercept = 0) +
#   geom_errorbar(aes(ymin = coeff.lower, ymax = coeff.upper, group = Test), 
#                 position = position_dodge(width = 0.5), width = 0) +
#   geom_point(shape = 21, size = 3, position = position_dodge(width = 0.5)) +
#   scale_fill_manual(values = c("NotSignificant" = "white",
#                                "< 20%" = "yellow",
#                                "< 10%" = "orange",
#                                "< 5%" = "red")) +
#   scale_shape_manual(values = c("GroupedDysmorphology" = 21,
#                                "ComplexVsEssential" = 22,
#                                "MultiClass" = 23,
#                                "DysmorphologyScore" = 24)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust =1)) + 
#   guides(fill = guide_legend(override.aes=list(shape=21))) + ggtitle("Deletions")
# 
# ggsave("ADM/NFLD.gs.deletion.burden.png", width = 12, height = 5)
# 
# ### plot duplication result
# ggplot(dup, aes(x = geneset, y = coefficient, fill = FDRRange)) +
#   geom_hline(yintercept = 0) +
#   geom_errorbar(aes(ymin = coeff.lower, ymax = coeff.upper, group = Test), 
#                 position = position_dodge(width = 0.5), width = 0) +
#   geom_point(shape = 21, size = 3, position = position_dodge(width = 0.5)) +
#   scale_fill_manual(values = c("NotSignificant" = "white",
#                                "< 20%" = "yellow",
#                                "< 10%" = "orange",
#                                "< 5%" = "red")) +
#   scale_shape_manual(values = c("GroupedDysmorphology" = 21,
#                                 "ComplexVsEssential" = 22,
#                                 "MultiClass" = 23,
#                                 "DysmorphologyScore" = 24)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust =1)) +
#   guides(fill = guide_legend(override.aes=list(shape=21))) + ggtitle("Duplications")
# 
# 
# ggsave("ADM/NFLD.gs.duplication.burden.png", width = 12, height = 5)
