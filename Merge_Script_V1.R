##Script To Merge Databases produced using DEG_finder V2 - Will identify mutuallly differential genes between two different comparisons. It will help with narrowing down relevent hits#

library(ggplot2)
library(reshape2)
library(DESeq2)
library(pheatmap)

#Read In upregualted Genes#
upreg1 <- read.csv("Significant_Upregulated_Genes_ptcG4_v_7AGFP_SG_D8.csv" , row.names = 1)
upreg2 <- read.csv("Significant_Upregulated_Genes_ptcG4_v_Fer12_SG_D8.csv" , row.names = 1)

#Read In Downregulated Genes#
downreg1 <- read.csv("Significant_Downregulated_Genes_ptcG4_v_7AGFP_SG_D8.csv", row.names = 1)
downreg2 <- read.csv("Significant_Downregulated_Genes_ptcG4_v_Fer12_SG_D8.csv", row.names = 1)

#Merge Upregualted  Data#
upregMerged <- merge(x = upreg2, y = upreg1, by = 0)
write.csv(upregMerged, "Fer12_7AGFP_SG_D8_Merged_Upregulated_Genes.csv")

#Merge Downregualted Data#
downregMerged <- merge(x = downreg1, y = downreg2, by = 0)
write.csv(downregMerged, "Fer12_7AGFP_SG_D8_Merged_Downregulated_Genes.csv")
