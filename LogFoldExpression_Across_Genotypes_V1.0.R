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

