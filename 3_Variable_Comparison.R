#This Script can filter any two sample-sets from a rawcount master sheet and compare them each other - eg control (PtcG4 vs Fer12_D6). The script can also complete a three way comparison useful for running
#control samples against tumor samples at multiple time points.. eg, PtcG4 Vs RasYki_D5 Vs RasYki_D8. The script can then take this
#data and create graphs to visualize both log2foldchange against the control group, log2 fold normalized counts or rawcounts against the control.
library(ggplot2)
library(reshape2)
library(DESeq2)
library(pheatmap)
library(dplyr)

###VARIABLES TO SET###
rawdata <- read.csv("All_WD_Samples_rawcounts.csv", header = T, row.names = 1) # <- Rawcounts Matrix File Name
info <- read.csv("All_WD_Samples_info.csv", header = T, row.names = 1) # <- Info Matrix File Name
experimentalGroup <- "All_WD_Samples" # <- What does the rawcounts file contain - eg. "All_WD_samples"
group1 <- "PtcG4_D6" # Column Name of Rawcounts matrix file for the 2 way comparison
group2 <- "7AGFP_D8" # Column Name of Rawcounts matrix file for the 2 way comparison
experimentalGroup <- "All_WD_Samples" # <- What does the rawcounts file contain - eg. "All_WD_samples"


#Filename Generator#
normalisedMatrixOutputFilename <- paste0(experimentalGroup,"_Normalised_Counts.csv") # <- Specify output file name for normalised counts matrix
deseqResultsFilename <- paste0(group1,"_V_",group2,"_DESEQ_Results.csv")
sigFoldUpFilename  <- paste0(group1,"_V_",group2,"_Significant_Upregulated_Genes.csv")
sigFoldDownFilename <- paste0(group1,"_V_",group2,"_Significant_Downregulated_Genes.csv")
allSigGenesFilename <- paste0(group1,"_V_",group2,"_All_Significant_Genes.csv")


###RUN DESEQ2###
#Create DESeq data set using the two reference files outlined above:#
dds <- DESeqDataSetFromMatrix(rawdata, info, ~GenotypeDay)

#Remove lowly expressed genes from the data set that have a total row count of below X
keep <- rowSums(counts(dds)) >= 50
dds <- dds[keep,]

#Run Deseq normalization. Factors in difference in library sizes and estimates dispersion levels.#
ddsDE <- DESeq(dds)

#Extract All Normalised Counts#
normcounts <- counts(ddsDE, normalized = T )

#Export normalized read counts to a .cvs file for later use#
write.csv(normcounts, normalisedMatrixOutputFilename)

#Run DESEQ2 on normalized counts file#
resultsNames(ddsDE)
res <- results(ddsDE, contrast=c("GenotypeDay", "PtcG4_D6", "7AGFP_D8"))
statsDE <- res[order(res$padj),]
write.csv(statsDE, deseqResultsFilename)


###DIFFERENTIAL GENE IDENTIFIER (For 2 Variable Analysis)###
#Remove all genes with a total row count of < X (50 Recommended)
moreThanFifty <- subset(statsDE, baseMean >50)

#Remove all genes with a padj (ajusted pvalue) of less that X (.05 recommended)
filteredStatsDE <- subset(moreThanFifty, padj <0.05)

#Pull Out all genes with a log fold change of > X into a list and export to CSV#
sigFoldGenesUp <- subset(filteredStatsDE, log2FoldChange > 1.4)
write.csv(sigFoldGenesUp, sigFoldUpFilename)

#Pull out all genes with a log fold change of < X into a list and export to a CSV#
sigFoldGenesDown <- subset(filteredStatsDE, log2FoldChange < -1.4)
write.csv(sigFoldGenesDown, sigFoldDownFilename)

#Write all Differentially expressed genes to a csv file#
sigFoldAllGenes <- rbind(sigFoldGenesUp, sigFoldGenesDown)
write.csv(sigFoldAllGenes, allSigGenesFilename)



##Three Way Data Comparison by visualizing log2fold changes against a control.
##SET VARIABLES##
controlSample <- "PtcG4_D6"
liveSample1 <- "7AGFP_D5"
liveSample2 <- "7AGFP_D8"

#Run DESeq2 using 3+ way comparison test (LRT)#
dds <- DESeqDataSetFromMatrix(rawdata, info, ~GenotypeDay)
ddsTimeSeries <- DESeq(dds, test = "LRT" , reduced = ~1)

#Write results of Control V Group 1 and Control V Group 2 to variables res1Table & res2Table
res1 <- results(ddsTimeSeries, contrast= c("GenotypeDay", controlSample, liveSample1))
write.csv(res1, "tempFile1.csv")
res1Table <- read.csv("tempFile1.csv")
file.remove("tempFile1.csv")
res2 <- results(ddsTimeSeries, contrast = c("GenotypeDay", controlSample, liveSample2))
write.csv(res2, "tempFile2.csv")
res2Table <- read.csv("tempFile2.csv")
file.remove("tempFile2.csv")

#Merge Data sets# 
resMerged <- merge(x = res1Table, y= res2Table, by = 0)

#Clean Merged Table#
resMergedClean <- subset(resMerged, select = c(X.x, log2FoldChange.x, log2FoldChange.y))
row.names(resMergedClean) <- resMergedClean$X.x
resMergedClean <- subset(resMergedClean, select = -c(X.x))
colnames(resMergedClean) <- c(paste0(liveSample1,"_log2FoldChange"), paste0(liveSample2,"_log2FoldChange"))
resMergedFinal <- na.omit(resMergedClean)
write.csv(resMergedFinal, paste0(liveSample1,"_&_",liveSample2,"_V_",controlSample,".csv"))


##SELECT GENE/S FOR ANALYSIS##
genesForAnalysis <- list("FBgn0261529", "FBgn0261581", "FBgn0261797") ## <- Paste Gene List Here

#Pull row.names data for all genes specified above#
wantedGenesTable <- resMergedFinal[paste0(genesForAnalysis), ]
#Melt Matrix#
wantedGenesTable <- melt(as.matrix(wantedGenesTable))
#Rename Columns#
colnames(wantedGenesTable) <- c("Gene_ID", "Genotype", "Log2FoldChange") ## <- May need to edit wanted column names.
#Rename X.axis labels#
wantedGenesTable$Genotype <- ifelse(grepl("D5", wantedGenesTable$Genotype), "7AD5","7AD8") ##May Need to edit wanted X axis label names.
#GG plot##
log2FoldChangePlot <- ggplot(wantedGenesTable, aes(x =Genotype, y =Log2FoldChange, color= Genotype)) +
  geom_point() +
  facet_grid(~Gene_ID) +
  theme(axis.text.x = element_text(angle = 90))
print(log2FoldChangePlot)


##Graph Generator to compare Log2FoldChange. (Normalized Count) Levels across selected genotypes##

##SELECT GENE/S FOR ANALYSIS##
genesForAnalysis <- list("FBgn0030411", "FBgn0040309", "FBgn0027087") ## <- Paste Gene List Here
##SELECT DATASETS TO ANALYSE (Row headers on normalised counts file)
GenotypesForAnalysis <- list("PtcG4_D5", "RasYki_D6", "RasYki_D8", "Fer12OG_D6", "Fer12OG_D8", "Fer12WT_D6", "ImpL2i_D6","ImpL2i_D8") ## <- Select Genotypes to analyse here

#Convert Normcounts large dataset to matrix#
normTable <- read.csv(normalisedMatrixOutputFilename, row.names = 1)
#Pull Data For selected genotypes across to wanted genes##
wantedGenesTable <- normTable[paste(genesForAnalysis), ]
#Melt Matrix#
wantedGenesTable <- melt(as.matrix(wantedGenesTable))

#Rename Columns
colnames(wantedGenesTable) <- c("Gene_ID", "Genotype", "Normcount")
ggplotOrder <- c("PtcG4_D6", "RasYki_D5", "RasYki_D8", "Fer12OG_D6", "Fer12OG_D8", "Fer12WT_D6", "ImpL2i_D6", "ImpL2_D8")
#Reorder Melted table#
wantedGenesTable$Genotype <- ifelse(grepl("PtcG4_D6", wantedGenesTable$Genotype), "PtcG4_D6",
                             ifelse(grepl("RasYki_D5", wantedGenesTable$Genotype), "RasYki_D5", 
                             ifelse(grepl("RasYki_D8", wantedGenesTable$Genotype), "RasYki_D8",
                             ifelse(grepl("Fer12OG_D6", wantedGenesTable$Genotype), "Fer12OG_D6",
                             ifelse(grepl("Fer12OG_D8", wantedGenesTable$Genotype), "Fer12OG_D8",
                             ifelse(grepl("Fer12WT_D6", wantedGenesTable$Genotype), "Fer12WT_D6",
                             ifelse(grepl("ImpL2i_D6", wantedGenesTable$Genotype), "ImpL2i_D6", "ImpL2i_D8" )))))))
#GGplot#
log2NormCounts <- ggplot(wantedGenesTable, aes(x = Genotype, y = log2(Normcount+1), color = Genotype)) +
  geom_point() +
  facet_grid(~Gene_ID) +
  theme(axis.text.x = element_text(angle = 90), legend.position = "none")
print(log2NormCounts)

