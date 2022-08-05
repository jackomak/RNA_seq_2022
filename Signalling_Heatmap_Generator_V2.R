#Load Libraries#
library(ComplexHeatmap)
library(circlize)

#Load log2FoldChange Counts#
mastersheet <- read.csv("Signalling_Pathway_Analysis_Mastersheet.csv", header = TRUE, row.names = 1)

#Create Core Heatmap Database#
coreHeatmap <- mastersheet
coreHeatmap[,c(1,7,8,9,10,11,12,13,14)] <- NULL 
coreHeatmap <- as.matrix(coreHeatmap)
colnames(coreHeatmap) <- c("RasYki_D5", "RasYki_D8", "Fer_D6", "Fer_D8", "FerWT_D6")

#Generate Colour Scale#
colorScale <- colorRamp2(c(-3,0,3), c("blue", "white", "red"))

#Generate Clusters#
clusters <- mastersheet
clusters[,c(1,2,3,4,5,6,7,8,10,11,12,13,14)] <- NULL
LA <- rowAnnotation(Group = clusters$Pathway, col = list(Group = c("EGFR/Ras" = "Red",
                                                                 "Dpp" = "Green",
                                                                 "Hedgehog" = "Blue",
                                                                 "Wingless" = "#FFFF00",
                                                                 "Hippo" = "#FFC0CB",
                                                                 "JNK" = "#00FFFF",
                                                                 "Jak-Stat" = "Black",
                                                                 "Notch" = "Orange",
                                                                 "Src" = "Purple")), border = TRUE)
#Normalised Count Annotation#
normalisedCountValues <- mastersheet
normalisedCountValues[,c(1,2,3,4,5,6,7,8,9,11,12,13,14)] <- NULL 
RA <- rowAnnotation(LCNV = anno_barplot(normalisedCountValues), border = TRUE)

#Generate Heatmap#
signallingHeatmap <- Heatmap(coreHeatmap,
                           col = colorScale,
                           column_title = "Salivary Gland",
                           row_title = "Gene",
                           row_title_gp = gpar(col = "white"),
                           cluster_columns = FALSE,
                           column_title_side = "top",
                           row_title_side = "right",
                           name = "Log2Fold Change",
                           row_order = 1:88,
                           row_gap = unit(1, "mm"),
                           border = TRUE,
                           left_annotation = LA,
                           right_annotation = RA,
                           column_names_rot = 90,
                           row_split = clusters$Pathway)
plot(signallingHeatmap)

#Create Brain Heatmap
#Create Brain Core Heatmap Data#
brainHeatmap <- mastersheet
brainHeatmap[,c(1,2,3,4,5,6,9,10,11,12,13,14)] <- NULL
brainHeatmap <- as.matrix(brainHeatmap)
colnames(brainHeatmap) <- c("RasYki_D5", "RasYki_D8")

#Brain Normalised Count Annotation#
brainNormalisedCountValues <- mastersheet
brainNormalisedCountValues[,c(1,2,3,4,5,6,7,8,9,10,12,13,14)] <- NULL
BRA <- rowAnnotation(LCNV = anno_barplot(brainNormalisedCountValues), border = TRUE)

#Generate Brain Heatmap#
brainSignallingHeatmap <- Heatmap(brainHeatmap,
                                col = colorScale,
                                column_title = "Brain",
                                cluster_columns = FALSE,
                                column_title_side = "top",
                                row_title_side = "right",
                                row_title = "Gene ID",
                                name = "Log2FoldChange",
                                show_heatmap_legend = FALSE,
                                row_order = 1:88,
                                row_gap = unit(1, "mm"),
                                border = TRUE,
                                right_annotation = BRA,
                                column_names_rot = 90,
                                row_split = clusters$Group)

#Concatenate salivary gland and brain heatmaps#
plot(brainSignallingHeatmap)
finalHeatmap <- signallingHeatmap + brainSignallingHeatmap
plot(finalHeatmap)

