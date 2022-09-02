#! R

#GeneSearch.R - An R script to search and generate heatamaps for a list of specified genes from the RNA seq datasets developed by the hamaratoglu lab 2021. Author Jack Bruton.

#Libraries#
library(ComplexHeatmap)
library(tidyverse)


#Read in raw log fold change database
rawDataTable <- read.csv("GeneSearch_LFC_Database", header = TRUE, sep = ",", row.names = 1)

#VARIABLES TO SET#
geneList <- c("gene1", "gene2", "geneN") #A list of genes to generate a heatmap for.
tissuesForHeatmap <- c("Wingdisc", "Salivarygland", "Brain")
genotypesForAnalysis <- c("PtcG4_D6", "Yw_D5", "RasYki_D5", "RasYki_D8", "Fer12OG_D6", "Fer12OG_D8", "Fer12WT_D6","ImpL2D6", "ImpL2D8")