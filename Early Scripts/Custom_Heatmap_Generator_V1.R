###This Script will take a normalized counts file and output file generated using DESeq2 and create a lists of up/downregualted genes and create heatmap plots#

#Load Librarys#
library(ggplot2)
library(reshape2)
library(DESeq2)
library(reshape2)
library(pheatmap)

#Read In Differential expression Data#
statsDE <- read.csv("DEG_Table_Sorted_ptcG4_7AGFP_D8.csv", row.names= 1)

#Read In Expression file#
normExp <- read.csv("ptcG4_7AGFP_D8_normalized_counts.csv", row.names = 1)
#Remove N/A P values#
statsDE <- na.omit(statsDE)

#data filtering#

#Remove all rows of statsDE with a Pvalue of less than X and rowsum of less than 50.
moreThanFifty <- subset(statsDE, baseMean >50)
filteredStatsDE <- subset(moreThanFifty, padj <0.05)


#Pull out all rows with log fold change of >2#
sigFoldGenesUp <- subset(filteredStatsDE, log2FoldChange > 1.4)
write.csv(sigFoldGenesUp, "Significant_Upregulated_Genes_ptcG4_v_7AGFP_SG_D8.csv")

#Pull out all genes with a log fold change of <2#
sigFoldGenesDown <- subset(filteredStatsDE, log2FoldChange < -1.4)
write.csv(sigFoldGenesDown, "Significant_Downregulated_Genes_ptcG4_v_7AGFP_WD_D8.csv")

#Merge the above two datasets to create list of all DEG's#
SigFoldAllGenes <- rbind(sigFoldGenesUp, sigFoldGenesDown)

##Heat map Generation##

#Define Column Names to keep#
keep <- c("Row.names", "ptcG4_SG_D6_S36", "ptcG4_SG_D6_S37", "ptcG4_SG_D6_S38", "Fer12OG_SG_D6_S39", "Fer12OG_SG_D6_S40", "Fer12OG_SG_D6_S41", "Fer12OG_SG_D8_S44", "Fer12OG_SG_D8_S45", "Fer12OG_SG_D8_S46")

#Heatmap For Left Bias genes Only#
heatTableUp <- merge(x = sigFoldGenesUp, y = normExp, by = 0)
heatTableUp <- heatTableUp[order(heatTableUp$padj),]
keptUp <- heatTableUp[, names(heatTableUp) %in% keep]
newKeptUp <- keptUp[, -1]
rownames(newKeptUp) <- keptUp[, 1]
pheatmap(log2(newKeptUp+1), show_rownames = F)

#Heatmap for Right Bias genes only#
heatTableDown <- merge(x = sigFoldGenesDown, y = normExp, by = 0)
heatTableDown <- heatTableDown[order(heatTableDown$padj),]
keptDown <- heatTableDown[, names(heatTableDown) %in% keep]
newKeptDown <- keptDown[, -1]
rownames(newKeptDown) <- keptDown[, 1]
pheatmap(log2(newKeptDown+1), show_rownames = F)

#Heatmap for all DEG's#
heatTableAll <- merge(x = SigFoldAllGenes, y = normExp, by = 0)
heatTableAll <- heatTableAll[order(heatTableAll$padj),]
keptAll <- heatTableAll[,names(heatTableAll) %in% keep]
newKeptAll <- keptAll[, -1]
rownames(newKeptAll) <- keptAll[, 1]
pheatmap(log2(newKeptAll+1), show_rownames = F)

##Visualizing Expression differences in the top 15 up regulated DEG's#
topX <- newKeptUp[1:15,]
topXm <- melt(as.matrix(topX))
topXm <- topXm[order(topXm$value),]
colnames <- c("gene_ID", "lib", "read_count", "Day")
names(topXm) <- c("gene_ID", "lib", "read_count")
topXm$tissue <- ifelse(grepl("ptc", topXm$lib), "CD6", ifelse(grepl("D8", topXm$lib), "FD8", "FD6"))
geneExp <- ggplot(topXm, aes(x = tissue, y= log2(read_count+1), color = tissue)) +
  geom_point() +
  facet_grid(~gene_ID)
print(geneExp)

##Visualising Expression differences in the top 15 down regulated DEG's##
topX <- newKeptDown[1:15,]
topXm <- melt(as.matrix(topX))
topXm <- topXm[order(topXm$value),]
colnames <- c("gene_ID", "lib", "read_count", "Day")
names(topXm) <- c("gene_ID", "lib", "read_count")
topXm$tissue <- ifelse(grepl("ptc", topXm$lib), "CD6", ifelse(grepl("D8", topXm$lib), "FD8", "FD6") )
geneExp <- ggplot(topXm, aes (x = tissue, y = log2(read_count+1), color = tissue)) +
  geom_point() +
  facet_grid(~gene_ID)
print(geneExp)


