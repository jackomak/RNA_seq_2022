#! R

#GeneSearch.R - An R script to search and generate heatamaps for a list of specified genes from the RNA seq datasets developed by the hamaratoglu lab 2021. Author Jack Bruton.

#Libraries#
library(ComplexHeatmap)
library(tidyverse)
library(readxl)

#Read in raw log fold change database

wdLfcTable <- read_excel("All_Tissues_LFC_Database.xlsx", sheet = 1)
sgLfcTable <- read_excel("All_Tissues_LFC_Database.xlsx", sheet = 2)
bLfcTable <- read_excel("All_Tissues_LFC_Database.xlsx", sheet = 3)


#VARIABLES TO SET#
geneList <- c("gene1", "gene2", "geneN") #A list of genes to generate a heatmap for.
tissuesForHeatmap <- c("Wingdisc", "Salivarygland", "Brain") #The tissues to include in the heatmap.
genotypesForAnalysis <- c("PtcG4_D6", "Yw_D5", "RasYki_D5", "RasYki_D8", "Fer12OG_D6", "Fer12OG_D8", "Fer12WT_D6","ImpL2D6", "ImpL2D8") #The genotypes to include in the heatmap.

