##Script To Merge Databases produced using DEG_finder V2 - Will identify mutuallly differential genes between two different comparisons. It will help with narrowing down relevent hits#

###LIBRARIES###
library(ggplot2)
library(reshape2)
library(DESeq2)
library(pheatmap)

###VARIABLES TO BE SET###
upregFile1 <- "<UPREGULATED_GENES_GENOTYPE_1.csv>"
upregFile2 <- "<UPREGULATED_GENES_GENOTYPE_2.csv>"
downRefFile1 <- "<DOWNREGULATED_GENES_GENOTYPE_1.csv">
downRegFile2 <- "<DOWNREGULATED_GENES_GENOTYPE_2.csv">
upregOutput <- "<MDEG_UPREGULATED_OUTPUT_FILENAME">
downregOutput <- "<MDEG_DOWNREGULATED_OUTPUT_FILENAME>"

###READ IN UPREGULATED GENES###
upreg1 <- read.csv(upregFile1 , row.names = 1)
upreg2 <- read.csv(upregFile2 , row.names = 1)

###READ IN DOWNREGULATED GENES###
downreg1 <- read.csv(downregFile1, row.names = 1)
downreg2 <- read.csv(downregFile2, row.names = 1)

###MERGE UPREGULATED DATA###
upregMerged <- merge(x = upreg2, y = upreg1, by = 0)
write.csv(upregMerged, upregOutput)

###MERGE DOWNREGULATED DATA###
downregMerged <- merge(x = downreg1, y = downreg2, by = 0)
write.csv(downregMerged, downregOutput)
