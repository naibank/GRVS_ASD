### set working directory 301 subjects
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

### install my burden analysis package
# library(devtools)
# install_github("naibank/GSBurden")

### clear existing variables
rm(list=ls())

### read data from R object
load("../data/2019.02.25-filteredCG_CNVs_cdsFixed_McgFixed.RData")
load("../data/2019.02.25-filteredSingles_CNVs_cdsFixed_McgFixed.RData")
load("../data/2019.02.25-filteredTrios_CNVs_cdsFixed_McgFixed.RData")

### get only column sample, chr, start, end, variant_type and size of CNVs
### match column and merge data from different tables
CG_filter <- CG_filter[, c(1:5, 7, 59)]
singles_filter <- singles_filter[, c(1:5, 7, 46)]
trios_filter <- trios_filter[, c(1:5, 7, 71)]
names(trios_filter)[4] <- "sample" 

merge <- rbind(singles_filter, trios_filter)
merge <- rbind(merge, CG_filter[, c(2:4, 1, 5:7)])
merge <- merge[merge$chr %in% paste0("chr", c(1:22)), ]
### read phenotype table
info <- read.delim("../data/2021.02.01-updated NFLD phenotype table IQ ADM and 3 status.txt", stringsAsFactors = F)

### get only columns that I use
info <- info[, c("WGS_ManuID", "DNA.source", "Relation", "Sex..CRV", "incl.aff.in.my.study", "incld.cohort", "remove.family.for.rare.burden",
                 "Platform", "Dysmorphology.classification", "Grouped_Dysmorphology", "Final.Dysmorphology.scores")]

# unique(merge$sample[!merge$sample %in% info$WGS_ManuID])
# [1] "3-0741-000A" "3-0179-000A" "3-0450-000A" "3-0213-000A" "3-0730-000A"

### rename ID of re-sequenced samples, if I understand it correctly
info$WGS_ManuID <- gsub("3-0741-000", "3-0741-000A", info$WGS_ManuID)
info$WGS_ManuID <- gsub("3-0179-000", "3-0179-000A", info$WGS_ManuID)
info$WGS_ManuID <- gsub("3-0450-000", "3-0450-000A", info$WGS_ManuID)
info$WGS_ManuID <- gsub("3-0213-000", "3-0213-000A", info$WGS_ManuID)
info$WGS_ManuID <- gsub("3-0730-000", "3-0730-000A", info$WGS_ManuID)

### read estimated PCA
pca <- read.delim("../estimated.pca.12March2019.tsv", stringsAsFactors = F)
pca <- pca[!pca$sample %in% c("3-0456-000", "3-0458-000"), ]
pca$sample <- gsub("A", "", pca$sample)

### merge PCA to phenotype table
info <- merge(info, pca[, c("sample", "pc1", "pc2", "pc3")], by.x = "WGS_ManuID", by.y = "sample", all.x =T)
#0466 -> no population stratification done for this sample.

### limit to proband with following criterions and merge CNVs table and phenotype table
info <- info[info$incl.aff.in.my.study == 1 & info$remove.family.for.rare.burden != "remove", ]
merge <- merge[merge$sample %in% info$WGS_ManuID, ]

### load my package
library(GSBurden)

### load geneset
load("../requiredData/gsNFLD.RData")

### limit to CNVs smaller than 10KB
merge <- merge[merge$CNV_size <= 10000, ]
cnvs <- CNVSet(merge$sample, merge$chr, merge$start, merge$end, merge$variant_type)
cnvs$size <- cnvs$end - cnvs$start + 1

write.table(cnvs, "../data/all.small.cnvs.may6.txt", sep="\t", row.names=F, quote=F, col.names=T)

### read coding exons table
gene.in <- na.omit(read.delim("../../../ReferenceData/geneInfo2019/hg19_refFlat_ext.txt", header = F, stringsAsFactors = F,
                              col.names = c("genesymbol", "isoid", "chr", "strand", "start", "end", "cds.start", "cds.end",
                                            "exons", "exon.starts", "exon.ends", "blank", "name2", "iso2", "var0", "var1", "var2", "entrezid", "blah")))
gene.in <- gene.in[, c("chr", "start", "end", "isoid", "genesymbol", "entrezid")]
genes <- GeneAnnotation(gene.in$entrezid, gene.in$chr, gene.in$start, gene.in$end, gene.in$genesymbol)

### get gene-set matrix
cnv.matrix <- getCNVGSMatrix(cnvs, genes, gsNFLD)
# sample.not.in <- info$WGS_ManuID[!info$WGS_ManuID %in% cnv.matrix$sample]
# sample.not.in <- data.frame("sample" = sample.not.in)
# sample.not.in[names(cnv.matrix)[-1]] <- 0
# cnv.matrix <- rbind(cnv.matrix, sample.not.in)

writeLines(as.character(unique(cnv.matrix$sample)), "samples.in.small.cnv.analysis.tsv")

### merge gene-set matrix with phenotype table
cnv.matrix <- merge(cnv.matrix, info[, c("WGS_ManuID", "Sex..CRV", "Platform", "Dysmorphology.classification",
                                         "Grouped_Dysmorphology", "Final.Dysmorphology.scores", "pc1", "pc2", "pc3")],
                    by.x = "sample", by.y = "WGS_ManuID", all.x = T)

### encode label variables used to build different regression models
cnv.matrix$Grouped_Dysmorphology[cnv.matrix$Grouped_Dysmorphology == "Equivocal"] <- "Non-Essential"
cnv.matrix$MultiClass <- ifelse(cnv.matrix$Dysmorphology.classification == "essential", 0, 1)
cnv.matrix$MultiClass <- ifelse(cnv.matrix$Dysmorphology.classification == "complex", 2, cnv.matrix$MultiClass)

cnv.matrix$ComplexVsEssential <- ifelse(cnv.matrix$MultiClass == 1, NA, 0)
cnv.matrix$ComplexVsEssential <- ifelse(cnv.matrix$MultiClass == 2, 1, cnv.matrix$ComplexVsEssential)

cnv.matrix$MultiClass <- factor(cnv.matrix$MultiClass, levels = 0:2)

cnv.matrix$GroupedDysmorphology <- ifelse(cnv.matrix$Grouped_Dysmorphology == "Non-Essential", 1, 0)

cnv.matrix$DysmorphologyScore <- factor(cnv.matrix$Final.Dysmorphology.scores, levels = sort(unique(as.numeric(cnv.matrix$Final.Dysmorphology.scores))))

write.table(cnv.matrix, "../cnv.smaller.coding.matrix.tsv", sep="\t", row.names=F, quote=F, col.names=T)
### separate gene-sets into 3 collections for multiple testing correction purpose
brainExp.gs <- grep("Expr|blueModule|Ilmn_BM|Bspan|FMR1|Synapse|Neurof", names(gsNFLD))
phenotype.gs <- grep("Ph", names(gsNFLD))

### define covariates
covariates <- c("Sex..CRV", "Platform", "pc1", "pc2", "pc3")

### run burden test for each of 3 gene-set collections
perm.values <- data.frame()

for(test in c("GroupedDysmorphology", "MultiClass")){
  brainExp.test.out <- CNVBurdenTest(cnv.matrix[!is.na(cnv.matrix[, test]),], gsNFLD[brainExp.gs], test, covariates, nperm = 1000)
  phenotype.test.out <- CNVBurdenTest(cnv.matrix[!is.na(cnv.matrix[, test]),], gsNFLD[phenotype.gs], test, covariates, nperm = 1000)
  coeff.out <- CNVBurdenTest(cnv.matrix[!is.na(cnv.matrix[, test]),], gsNFLD, test, covariates, 
                             permutation = F, correctGlobalBurden = F)$Test
  
  burden.test.out <- rbind(brainExp.test.out$Test, phenotype.test.out$Test)
  burden.test.out <- merge(burden.test.out[, -c(3:5)], coeff.out[, 3:5], by = "row.names", all.x = T)
  burden.test.out <- burden.test.out[, c(2:3, 6:8, 4:5)]
  perm.temp <- rbind(brainExp.test.out$Permutation.Test, phenotype.test.out$Permutation.Test)
  perm.temp$geneset <- gsub("[0-9]", "", rownames(perm.temp))
  perm.values <- rbind(perm.values, perm.temp)
  
  write.table(burden.test.out, sprintf("../%s.smaller.gsburden.tsv", test), sep="\t", row.names=F, quote=F, col.names=T)
}
perm.values$geneset <- rep(perm.temp$geneset, 4)
write.table(perm.values, "../perm.smaller.pvalues.tsv", sep="\t", row.names=F, quote=F, col.names=T)

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

perm.pvalues <- read.delim("../perm.smaller.pvalues.tsv", stringsAsFactors = F)
for(test in c("GroupedDysmorphology", "MultiClass")){
  dt <- read.delim(sprintf("../%s.smaller.gsburden.tsv", test), stringsAsFactors = F)
  
  brainExp.test <- dt[dt$geneset %in% names(gsNFLD)[brainExp.gs], ]
  phenotype.test <- dt[dt$geneset %in% names(gsNFLD)[phenotype.gs], ]
  
  perm.pvalues$geneset <- gsub("_deletion|_duplication", "", perm.pvalues$geneset)
  brain.perm <- perm.pvalues[perm.pvalues$geneset %in% names(gsNFLD)[brainExp.gs], ] 
  phenotype.perm <- perm.pvalues[perm.pvalues$geneset %in% names(gsNFLD)[phenotype.gs], ] 
  
  brainExp.test$permFDR <- sapply(brainExp.test$pvalue, getPermPvalue, brainExp.test$pvalue, brain.perm$pvalue)
  phenotype.test$permFDR <- sapply(phenotype.test$pvalue, getPermPvalue, phenotype.test$pvalue, phenotype.perm$pvalue)
  
  brainExp.test$permFDR <- sapply(1:nrow(brainExp.test), getPermFDR, brainExp.test)
  phenotype.test$permFDR <- sapply(1:nrow(phenotype.test), getPermFDR, phenotype.test)
  
  
  dt.out <- rbind(brainExp.test, phenotype.test)
  dt.out$coeff.lower[is.na(dt.out$coeff.lower) & !is.na(dt.out$coeff.upper)] <- -Inf
  dt.out$coeff.upper[!is.na(dt.out$coeff.lower) & is.na(dt.out$coeff.upper)] <- Inf
  
  write.table(dt.out, sprintf("../%s.smaller.gsburden.withPerm.tsv", test), sep="\t", row.names=F, quote=F, col.names=T)
}


### prepare data for plot
del <- data.frame()
dup <- data.frame()
for(test in c("GroupedDysmorphology", "MultiClass")){
  burden.test.out <- read.delim(sprintf("../%s.smaller.gsburden.withPerm.tsv", test), stringsAsFactors = F)
  test <- ifelse(test == "GroupedDysmorphology", "Two subtypes", "Three subtypes")
  
  burden.test.out$Test <- test
  burden.test.out$FDRRange <- "NotSignificant"
  burden.test.out$FDRRange[burden.test.out$permFDR < 0.20 & burden.test.out$pvalue < 0.05] <- "< 20%"
  burden.test.out$FDRRange[burden.test.out$permFDR < 0.10 & burden.test.out$pvalue < 0.05] <- "< 10%"
  burden.test.out$FDRRange[burden.test.out$permFDR < 0.05 & burden.test.out$pvalue < 0.05] <- "< 5%"
  burden.test.out$FDRRange <- factor(burden.test.out$FDRRange, levels = c("NotSignificant", "< 20%", "< 10%", "< 5%"))
  
  del <- rbind(del, burden.test.out[burden.test.out$type == "deletion", ])
  dup <- rbind(dup, burden.test.out[burden.test.out$type == "duplication", ])
}

del.significant <- unique(del$geneset[del$permFDR < 0.2])
dup.significant <- unique(dup$geneset[dup$permFDR < 0.2])
### plot deletion result
del$Test <- factor(del$Test, levels = c("Two subtypes", "Three subtypes"))
library(ggplot2)
ggplot(del, aes(x = geneset, y = coefficient, fill = FDRRange)) +
  
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = coeff.lower, ymax = coeff.upper, group = Test), 
                position = position_dodge(width = 0.7), width = 0) +
  geom_point(aes(shape = Test), size = 3, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("NotSignificant" = "white",
                               "< 20%" = "yellow",
                               "< 10%" = "orange",
                               "< 5%" = "red")) +
  scale_shape_manual(values = c(21, 23)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust =1)) + #coord_cartesian(ylim = c(-0.5, 1.5)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) + ggtitle("Deletions")

ggsave("../gs.smaller.deletion.burden.pdf", width = 12, height = 5)
dup$Test <- factor(dup$Test, levels = c("Two subtypes", "Three subtypes"))

### plot duplication result
ggplot(dup, aes(x = geneset, y = coefficient, fill = FDRRange)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = coeff.lower, ymax = coeff.upper, group = Test), 
                position = position_dodge(width = 0.7), width = 0) +
  geom_point(aes(shape = Test), size = 3, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = c("NotSignificant" = "white",
                               "< 20%" = "yellow",
                               "< 10%" = "orange",
                               "< 5%" = "red")) +
  scale_shape_manual(values = c(21, 23)) + coord_cartesian(ylim = c(-3, 1.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust =1)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) + ggtitle("Duplications")


ggsave("../gs.smaller.duplication.burden.pdf", width = 12, height = 5)


#### subset plot ####
# library(xlsx)
# gsname <- read.xlsx("../2019.03.14-geneset_CNVexonic_results_326samples.xlsx", sheetIndex = 4)
# # gsname <- gsname[, c(1, 3:4)]
# # write.table(gsname, "../gsName.tsv", sep="\t", row.names=F, quote=F, col.names = T)
# # gsname <- read.delim("../gsName.tsv", stringsAsFactors = F)
# 
# subset <- gsname$Original_name
# subset <- subset[subset %in% union(del.significant, dup.significant)]
# del.subset <- del[del$geneset %in% subset, ]
# dup.subset <- dup[dup$geneset %in% subset, ]
# del.subset <- merge(del.subset, gsname, by.x = "geneset", by.y = "Original_name", all.x = T)
# dup.subset <- merge(dup.subset, gsname, by.x = "geneset", by.y = "Original_name", all.x = T)
# del.subset$Ada_name <- factor(del.subset$Ada_name, levels = rev(gsname$Ada_name[c(1:14, 16, 15, 17:18)]))
# dup.subset$Ada_name <- factor(dup.subset$Ada_name, levels = rev(gsname$Ada_name[c(1:14, 16, 15, 17:18)]))
# 
# 
# del.s <- ggplot(del.subset, aes(x = Ada_name, y = coefficient, fill = FDRRange)) +
#   geom_hline(yintercept = 0) + xlab("Gene set") +
#   geom_errorbar(aes(ymin = coeff.lower, ymax = coeff.upper, group = Test), 
#                 position = position_dodge(width = 0.5), width = 0) +
#   geom_point(aes(shape = Test), size = 3, position = position_dodge(width = 0.5)) +
#   scale_fill_manual(values = c("NotSignificant" = "white",
#                                "< 20%" = "yellow",
#                                "< 10%" = "orange",
#                                "< 5%" = "red")) +
#   scale_shape_manual(values = c("GroupedDysmorphology" = 21,
#                                 "ComplexVsEssential" = 22,
#                                 "MultiClass" = 23,
#                                 "DysmorphologyScore" = 24)) +
#   theme_bw() + coord_flip(ylim = c(-1.2, 5)) + 
#   theme(axis.text.x = element_text(angle = 90, hjust =1),
#         legend.position = "none") +
#   guides(fill = guide_legend(override.aes=list(shape=21))) + ggtitle("Deletions")
# 
# ### plot duplication result
# dup.s <- ggplot(dup.subset, aes(x = Ada_name, y = coefficient, fill = FDRRange)) +
#   geom_hline(yintercept = 0) + xlab("Gene set") +
#   geom_errorbar(aes(ymin = coeff.lower, ymax = coeff.upper, group = Test), 
#                 position = position_dodge(width = 0.5), width = 0) +
#   geom_point(aes(shape = Test), size = 3, position = position_dodge(width = 0.5)) +
#   scale_fill_manual(values = c("NotSignificant" = "white",
#                                "< 20%" = "yellow",
#                                "< 10%" = "orange",
#                                "< 5%" = "red"), drop = FALSE) +
#   scale_shape_manual(values = c("GroupedDysmorphology" = 21,
#                                 "ComplexVsEssential" = 22,
#                                 "MultiClass" = 23,
#                                 "DysmorphologyScore" = 24)) + 
#   theme_bw() + coord_flip(ylim = c(-3.1, 3.1)) + 
#   theme(axis.text.x = element_text(angle = 90, hjust =1), axis.text.y = element_blank(), axis.title.y = element_blank()) +
#   guides(fill = guide_legend(override.aes=list(shape=21))) + ggtitle("Duplications")
# 
# library(cowplot)
# plot_grid(del.s, dup.s, rel_widths = c(1, 1))
# ggsave("../gs.smaller.subset.burden.pdf", width = 8.5, height = 9)

# ### Loci test for genes in the significant gene-sets
# gsburden <- rbind(read.delim("../ComplexVsEssential.smaller.gsburden.withPerm.tsv", stringsAsFactors = F),
#                   read.delim("../DysmorphologyScore.smaller.gsburden.withPerm.tsv", stringsAsFactors = F),
#                   read.delim("../GroupedDysmorphology.smaller.gsburden.withPerm.tsv", stringsAsFactors = F),
#                   read.delim("../MultiClass.smaller.gsburden.withPerm.tsv", stringsAsFactors = F))
# gsburden <- gsburden[gsburden$pvalue < 0.05 & gsburden$permFDR < 0.1, ]
# 
# loci.test <- CNVLociTest(cnvs, cnv.matrix, genes, "MultiClass", covariates, gsNFLD[unique(gsburden$geneset)],
#                          nsubject = 2, nperm = 1000)
# 
# grepInGeneset <- function(geneset, grep.pattern){
#   return(length(grep(grep.pattern, geneset)) > 0)
# }
# 
# for(i in 1:nrow(loci.test)){
#   gsTMP <- gsNFLD[unique(gsburden$geneset[gsburden$type == loci.test$type[i]])]
#   genesets  <- names(gsTMP)[sapply(gsTMP, grepInGeneset, paste(strsplit(loci.test$enzid[i], ",")[[1]], collapse = "|"))]
#   loci.test$genesets[i] <- paste(genesets, collapse = ";")
# }
# 
# loci.test <- loci.test[loci.test$genesets != "", ]
# write.table(loci.test, "../smaller.CNV.genes.tsv", sep="\t", row.names=F, quote=F, col.names=T)


####################################
######## global burden test ########

global.test <- data.frame()
for(test in c("GroupedDysmorphology", "ComplexVsEssential", "MultiClass", "DysmorphologyScore")){
  tmp.test <- CNVGlobalTest(cnv.matrix, test, covariates)
  tmp.test$Test <- test
  global.test <- rbind(global.test, tmp.test)
}

global.test$PvalueRange <- "NotSignificant"
global.test$PvalueRange[global.test$pvalue < 0.1] <- "< 0.1"
global.test$PvalueRange[global.test$pvalue < 0.05]  <- "< 0.05"
global.test$PvalueRange[global.test$pvalue < 0.01]  <- "< 0.01"
global.test$PvalueRange <- factor(global.test$PvalueRange, levels = c("NotSignificant", "< 0.1", "< 0.05", "< 0.01"))
write.table(global.test, "../global.smaller.burden.tsv", sep="\t", row.names=F, quote=F, col.names=T)

global.test <- read.delim("../global.smaller.burden.tsv", stringsAsFactors = F)
global.test <- global.test[global.test$Test %in% c("GroupedDysmorphology", "MultiClass"), ]
global.test$Test <- ifelse(global.test$Test == "GroupedDysmorphology", "Two subtypes", "Three subtypes")
global.test$Test <- factor(global.test$Test, levels = c("Two subtypes", "Three subtypes"))
ggplot(global.test[global.test$global == "gene_count", ], aes(x = type, y = coefficient, fill = PvalueRange)) +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = lowerbound, ymax = upperbound, group = Test), 
                position = position_dodge(width = 0.5), width = 0) +
  geom_point(aes(shape = Test), size = 3, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("NotSignificant" = "white",
                               "< 0.1" = "yellow",
                               "< 0.05" = "orange",
                               "< 0.01" = "red")) +
  scale_shape_manual(values = c(21, 23)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust =1)) +
  guides(fill = guide_legend(override.aes=list(shape=21))) + ggtitle("Global Burden")

ggsave("../global.smaller.burden.pdf", width = 6, height = 4)

