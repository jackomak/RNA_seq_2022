##This script will generate normalised counts data for a featurecounts rawcounts file and in combination with a DESEQ2 information matrix, summarise log2 expresion levels
##for particular genes of interest which can be specifed. The script also gives the option to remove unwanted genotypes. ##

###LIBRARIES###
library(ggplot2)
library(reshape2)
library(DESeq2)
library(pheatmap)
library(dplyr)

###VARIABLES TO SET###
rawData <- read.csv("<RAWCOUNTS_FILE>", header = T, row.names = 1) # <- Rawcounts Matrix File Name.
info <- read.csv("<INFO_FILE>", header = T, row.names = 1) # <- Info Matrix File Name.
experimentalGroup <- "<EXPERIMENT_NAME>" # <- What does the rawcounts file contain? - eg. "All_WD_samples"
experimentalCondition <- "<EXPERIMENTAL_CONDITION>" # <- Which experimental condition to test on <INFO_FILE> eg. GenotypeDay

###SELECT GENE/S FOR ANALYSIS###
genesForAnalysis <- list("FBgn0030411", "FBgn0040309", "FBgn0027087") ## <- Paste Gene List Here

###SELECT DATASETS TO ANALYSE ### (Row headers on normalised counts file)
GenotypesForAnalysis <- list("PtcG4_D5", "RasYki_D6", "RasYki_D8", "Fer12OG_D6", "Fer12OG_D8", "Fer12WT_D6", "ImpL2i_D6","ImpL2i_D8") ## <- Select Genotypes to analyse here

###FILENAME GENERATOR###
normalisedMatrixOutputFilename <- paste0(experimentalGroup,"_Normalised_Counts.csv") # <- May need to specify output file name for normalized counts matrix.

###RUN DESEQ2###
dds <- DESeqDataSetFromMatrix(rawData, info, ~ GenotypeDay)
ddsDE <- DESeq(dds)
normcounts <- counts(ddsDE, normalized = T )
write.csv(normcounts, normalisedMatrixOutputFilename)
normTable <- read.csv(normalisedMatrixOutputFilename, row.names = 1)

###PULL NORMALISED COUNT DATA TO WANTED GENE LIST###
wantedGenesTable <- normTable[paste(genesForAnalysis), ]
wantedGenesTable <- melt(as.matrix(wantedGenesTable))

###RENAME COLUMNS###
colnames(wantedGenesTable) <- c("Gene_ID", "Genotype", "Normcount")
ggplotOrder <- c("PtcG4_D6", "RasYki_D5", "RasYki_D8", "Fer12OG_D6", "Fer12OG_D8", "Fer12WT_D6", "ImpL2i_D6", "ImpL2_D8")
###Reorder Melted table###
wantedGenesTable$Genotype <- ifelse(grepl("PtcG4_D6", wantedGenesTable$Genotype), "PtcG4_D6",
                             ifelse(grepl("RasYki_D5", wantedGenesTable$Genotype), "RasYki_D5", 
                             ifelse(grepl("RasYki_D8", wantedGenesTable$Genotype), "RasYki_D8",
                             ifelse(grepl("Fer12OG_D6", wantedGenesTable$Genotype), "Fer12OG_D6",
                             ifelse(grepl("Fer12OG_D8", wantedGenesTable$Genotype), "Fer12OG_D8",
                             ifelse(grepl("Fer12WT_D6", wantedGenesTable$Genotype), "Fer12WT_D6",
                             ifelse(grepl("ImpL2i_D6", wantedGenesTable$Genotype), "ImpL2i_D6", "ImpL2i_D8" )))))))

#VISUALISE LOG2(NORMALISED_COUNTS) ACROSS GENOTYPES WITH GGPLOT)
log2NormCounts <- ggplot(wantedGenesTable, aes(x = Genotype, y = log2(Normcount+1), color = Genotype)) +
  geom_point() +
  facet_grid(~Gene_ID) +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none")
print(log2NormCounts)



