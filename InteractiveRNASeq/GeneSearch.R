#! R

#GeneSearch.R - An R script to search and generate heatamaps for a list of specified genes from the RNA seq datasets developed by the hamaratoglu lab 2021. Author Jack Bruton.

#Libraries#
library(ComplexHeatmap)
library(tidyverse)
library(readxl)

#VARIABLES TO SET#
geneList <- c("FBgn0015399", "FBgn0033395", "FBgn0026562") #A list of genes to generate a heatmap for.
tissuesForHeatmap <- c("Wingdisc", "Salivarygland", "Brain") #The tissues to include in the heatmap.
genotypesForAnalysis <- c("PtcG4_D6", "Yw_D5", "RasYki_D5", "RasYki_D8", "Fer12OG_D6", "Fer12OG_D8", "Fer12WT_D6","ImpL2D6", "ImpL2D8") #The genotypes to include in the heatmap.


#Read in raw log fold change database
wdLfcTable <- read_excel("All_Tissues_LFC_Database.xlsx", sheet = 1)
wdLfcTable <- column_to_rownames(wdLfcTable, var = "GeneID")

sgLfcTable <- read_excel("All_Tissues_LFC_Database.xlsx", sheet = 2)
sgLfcTable <- column_to_rownames(sgLfcTable, var = "GeneID")

bLfcTable <- read_excel("All_Tissues_LFC_Database.xlsx", sheet = 3)
bLfcTable <- column_to_rownames(bLfcTable, var = "GeneID")

#Filter datasets to only include rows that user has specified in "genelist" variable.
wdLfcTable <- wdLfcTable[rownames(wdLfcTable) %in% geneList, ]
sgLfcTable <- sgLfcTable[rownames(sgLfcTable) %in% geneList, ]
bLfcTable <- bLfcTable[rownames(bLfcTable) %in% geneList, ]

#Remove unwanted tissue datasets#
