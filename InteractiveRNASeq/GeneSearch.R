#! R

#GeneSearch.R - An R script to search and generate heatamaps for a list of specified genes from the RNA seq datasets developed by the hamaratoglu lab 2021. Author Jack Bruton.

#Libraries#
library(ComplexHeatmap)
library(tidyverse)
library(readxl)
library(circlize)

#Set raw data file and gene Name database#
rawData <- "ALL_Tissues_LFC_Database.xlsx"
geneNames <- read_excel(rawData, sheet = 4)

#VARIABLES TO SET#
geneList <- c("FBgn0015399", "FBgn0033395", "FBgn0026562") #A list of genes to generate a heatmap for.
tissuesForHeatmap <- c("Wingdisc", "Salivarygland", "Brain") #The tissues to include in the heatmap.
genotypesForAnalysis <- c("PtcG4_D6", "Yw_D5", "RasYki_D5", "RasYki_D8", "Fer12OG_D6", "Fer12OG_D8", "Fer12WT_D6", "ImpL2i_D6", "ImpL2i_D8") #The genotypes to include in the heatmap.
convertGeneIds <- TRUE # If TRUE geneIDs are converted to primary gene names in final heatmap.

#Read in raw log fold change databases for each tissue.
wdLfcTable <- read_excel(rawData, sheet = 1)
wdLfcTable <- left_join(wdLfcTable, geneNames, by = "GeneID")
wdLfcTable <- column_to_rownames(wdLfcTable, var = "GeneID")

sgLfcTable <- read_excel(rawData, sheet = 2)
sgLfcTable <- left_join(sgLfcTable, geneNames, by = "GeneID")
sgLfcTable <- column_to_rownames(sgLfcTable, var = "GeneID")

bLfcTable <- read_excel(rawData, sheet = 3)
bLfcTable <- left_join(bLfcTable, geneNames, by = "GeneID")
bLfcTable <- column_to_rownames(bLfcTable, var = "GeneID")

#Filter datasets to only include rows that user has specified in "genelist" variable.
wdLfcTable <- wdLfcTable[rownames(wdLfcTable) %in% geneList, ]
sgLfcTable <- sgLfcTable[rownames(sgLfcTable) %in% geneList, ]
bLfcTable <- bLfcTable[rownames(bLfcTable) %in% geneList, ]

#Remove unwanted tissue datasets#
genotypesForAnalysis <- append(genotypesForAnalysis, c("LFC", "Gene Name"))
wdLfcTable <- wdLfcTable[, colnames(wdLfcTable) %in% genotypesForAnalysis]
sgLfcTable <- sgLfcTable[, colnames(sgLfcTable) %in% genotypesForAnalysis]
bLfcTable <- bLfcTable[, colnames(bLfcTable) %in% genotypesForAnalysis]

#Select Which tissues to add to heatmap#
Wingdisc <- wdLfcTable
Salivarygland <- sgLfcTable
Brain <- bLfcTable

#Create Heatmap Color Scale.
colorScale <- colorRamp2(c(-3,0,3), c("blue", "white", "red"))

##HEATMAP GENERATOR LOOP##
heatmapList <- list()
for (tissue in tissuesForHeatmap) {
  
  #Create normalised count level annotation.
  ncValues <- get(as.name(tissue))[, colnames(get(as.name(tissue))) == "LFC"]
  rowAnnotation <- rowAnnotation(LFC = anno_barplot(ncValues), border = TRUE)
  
  #Create Core Heatmap database.
  coreHeatmap <- get(as.name(tissue))[, colnames(get(as.name(tissue))) != "LFC"]
 
  if (convertGeneIds == TRUE){
    coreHeatmap <- rownames_to_column(coreHeatmap, var = "GeneID")
    coreHeatmap <- column_to_rownames(coreHeatmap, var = "Gene Name")
    coreHeatmap$GeneID = NULL
  } else {
    coreHeatmap$`Gene Name` = NULL
  }
  
  #Plot heatmap.
  heatmap <- Heatmap(as.matrix(coreHeatmap),
    col = colorScale,
    column_title = tissue,
    row_title = "Gene",
    row_title_gp = gpar(col = "white"),
    cluster_columns = FALSE,
    column_title_side = "top",
    row_title_side = "right",
    name = "Log2Fold Change",
    row_gap = unit(1, "mm"),
    border = TRUE,
    right_annotation = rowAnnotation,
    column_names_rot = 90)
  
heatmapList <- append(heatmap, heatmapList)
}

#Plot Heat map Based on Tissues Provided.
HeatmapListLength <- length(heatmapList)

if(HeatmapListLength == 1) {
  plot(heatmapList[[1]])
} else if (HeatmapListLength == 2) {
  plot(heatmapList[[2]] + heatmapList[[1]])
} else {
  plot(heatmapList[[3]] + heatmapList[[2]] + heatmapList[[1]])
}

