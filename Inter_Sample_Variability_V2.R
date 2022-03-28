#Install The required softwares using BioCManager:#
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", "RColorBrewer", "pheatmap", "matrixStats", "ggplot2", "reshape", "DeSeq2")

#InStall DESeq2#
BiocManager::install("reshape2")

#Load DESeq2 Library#
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("matrixStats")
library("ggplot2")
library("reshape2")

#load count matrix#
dat <- read.table("ptcG4_7AGFP_D8_readcounts.txt", header = T, row.names = 1 )
info <- read.table("ptcG4_7AGFP_D8_info.txt", header = T, sep = '\t')

#Create DESeq Datset using the two reference files outlined above:#
dds <- DESeqDataSetFromMatrix(dat, info, ~Genotype)

#Remove lolly expressed genes from the data set that have a total row count of below 10#
keep <- rowSums(counts(dds)) >= 50
dds <- dds[keep,]

#Run Deseq normalization. Factors in difference in library sizes and estimates dispersion levels.#
ddsDE <- DESeq(dds)

#Extracting Normalized Read Counts#
normcounts <- counts(ddsDE, normalized = T )

#Export normalized read counts to a .cvs file for later use#
write.csv(normcounts, "ptcG4_7AGFP_D8_normalized_counts.csv")

#Visualising DeSeq Results. if corrected p value is less than 0.5 define this as a differential expressed gene. 
res <-results(ddsDE, alpha = 0.05)
summary(res)

#Sorting and extracting Deferentially expressed genes into a csv file sorted by p value.
resOrdered <- res[order(res$padj),]
write.csv(resOrdered, "DEG_Table_Sorted_ptcG4_7AGFP_D8.csv")

#Plotting DEGs with mean expression along X and log fold change across the Y
plotMA(ddsDE)

###Graphs to Assess Variability between data###

#PCA Plot, vst performs a logarithmic transformation by log2.
vsd <- vst(dds, blind = F)
plotPCA(vsd , intgroup= c("X"))

#Heatmap plot#
sampleDist <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDist)
rownames(sampleDistMatrix) <- vsd$Genotype_Day_Organ
colors <- colorRampPalette( rev(brewer.pal(9, "YlOrRd")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDist,
         clustering_distance_cols = sampleDist,
         col = colors)

