library(ggplot2)
library(reshape)
library(DESeq2)
library(reshape2)
library(pheatmap)

#Expression file#
normExp <- read.csv("Fer12OGSG_vs_ptcG4yw_D6_normalized_counts.csv", row.names = 1)

#DE Results File#
statsDE <- read.csv("DEG_Table_Sorted_Fer12OGSG_vs_ptcG4yw_D6.csv", row.names = 1)

#Remove N/A P values#
statsDE <- na.omit(statsDE)

#Is the data significant?#
statsDE$sig <- ifelse(statsDE$padj <= 0.001, "sig (p = <0.001)", "not")

#plotMA 
plotMA <- ggplot(statsDE, aes(x = log10(baseMean), y = log2FoldChange, color = sig)) + 
  geom_point() + 
  theme(legend.position = "none")
print(plotMA)                

#volcano
volcano <- ggplot(statsDE, aes(x = log2FoldChange, y =-log10(padj), color = sig)) + geom_point()
print(volcano)

##HEatmap##

#Getting top X rows of most differential expressed genes.# Merge pulls out same gene names from expression data 
top <- statsDE[1:15,]
top <- merge(top, normExp, by = 0)
keep <-c("Row.names", "rasyki_D5_2_S11","rasyki_D5_3_S12", "FRT82B_d5_S1", "FRT82B_d5_S2", "FRT82B_d5_S3")

#Look through keep and sift matrix to remove any column names that aren't specified above#
top <- top[,names(top) %in% keep]

#Fancy magic I dont understand# make gene ID's as row names.
newtop <- top[, -1]
rownames(newtop) <- top[,1]

#plot heatmap
pheatmap(log2(newtop+1), show_rownames = F)

#Analyse Specific Genes#
topX <- newtop[1:15,]
topXm <- melt(as.matrix(topX))
topXm <- topXm[order(topXm$value),]
colnames(topXm) <- c("gene_ID", "lib", "pvalue", "Genotype")
names(topXm) <- c("gene_ID", "lib", "exp")
topXm$tissue <- ifelse(grepl("ras", topXm$lib), "7A", "ptcG4")

geneExp <- ggplot(topXm, aes(x = tissue, y = log2(exp +1), color = tissue)) +
  geom_point() +
  facet_grid(~gene_ID)
print(geneExp)

