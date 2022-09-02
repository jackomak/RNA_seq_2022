#Load Libraries#
library(ComplexHeatmap)
library(circlize)

#Load log2FoldChange Counts#
mastersheet <- read.csv("Ecdysone_Gene_List_SG.csv", header = TRUE, row.names = 1)

#Create Core Heatmap Database#
coreHeatmap <- mastersheet
coreHeatmap[,c(1,2,7,8,9,10)] <- NULL 
coreHeatmap <- as.matrix(coreHeatmap)
colnames(coreHeatmap) <- c("RasYki_D5", "RasYki_D8", "Ferritin_D6", "Ferritin_D8")

#Generate Colour Scale#
colorScale <- colorRamp2(c(-3,0,3), c("blue", "white", "red"))

#Generate Clusters#
clusters <- mastersheet
clusters[,c(1,3,4,5,6,7,8,9,10)] <- NULL
LA <- rowAnnotation(Group = clusters$Group, col = list(Group = c("AA Transporter" = "Red",
                                                                 "Clotting Factor" = "Green",
                                                                 "Ecdysone-Induced" = "Blue",
                                                                 "Ecdysone Related" = "#FFFF00",
                                                                 "Edysone Sythesis" = "#FFC0CB",
                                                                 "Nutritional Response" = "#00FFFF")), border = TRUE)
#Normalised Count Annotation#
normalisedCountValues <- mastersheet
normalisedCountValues[,c(1,2,3,4,5,6,7,8,10)] <- NULL 
RA <- rowAnnotation(LCNV = anno_barplot(normalisedCountValues), border = TRUE)
                                    
#Generate Heatmap#
ecdysoneHeatmap <- Heatmap(coreHeatmap,
                           col = colorScale,
                           column_title = "Salivary Gland",
                           row_title = "Gene ID",
                           row_title_gp = gpar(col = "white"),
                           cluster_columns = FALSE,
                           column_title_side = "top",
                           row_title_side = "right",
                           name = "Log2Fold Change",
                           row_order = 1:50,
                           row_gap = unit(1, "mm"),
                           border = TRUE,
                           left_annotation = LA,
                           right_annotation = RA,
                           column_names_rot = 90,
                           row_split = clusters$Group)
plot(ecdysoneHeatmap)

#Create Brain Comparison Heatmap#
#Create core brain heatmap data#
brainHeatmap <- mastersheet
brainHeatmap[,c(1,2,3,4,5,6,9,10)] <- NULL
brainHeatmap <- as.matrix(brainHeatmap)
colnames(brainHeatmap) <- c("RasYki_D5", "RasYki_D8")

#Brain Normalised Count Annotation#
brainNormalisedCountValues <- mastersheet
brainNormalisedCountValues[,c(1,2,3,4,5,6,7,8,9)] <- NULL
BRA <- rowAnnotation(LCNV = anno_barplot(brainNormalisedCountValues), border = TRUE)

#Generate Brain Heatmap#
brainEcdysoneHeatmap <- Heatmap(brainHeatmap,
                                col = colorScale,
                                column_title = "Brain",
                                cluster_columns = FALSE,
                                column_title_side = "top",
                                row_title_side = "right",
                                row_title = "Gene ID",
                                name = "Log2FoldChange",
                                show_heatmap_legend = FALSE,
                                row_order = 1:50,
                                row_gap = unit(1, "mm"),
                                border = TRUE,
                                right_annotation = BRA,
                                column_names_rot = 90,
                                row_split = clusters$Group)
plot(brainEcdysoneHeatmap)
finalHeatmap <- ecdysoneHeatmap + brainEcdysoneHeatmap
plot(finalHeatmap)

