#! - Rstudio / Shiny

#app.R - An interactive application to visualize log fold change in gene expression derived from data obtained from the 2022
#RNA-Seq cachexia project samples.

#Shiny packages ----
library(shiny)
library(shinythemes)

#core packages for data visualization ----
library(ComplexHeatmap)
library(tidyverse)
library(readxl)
library(circlize)

#Set core variable lists ----
genotypesForAnalysisNames <- list("RasYki Day 5", "RasYki Day 8", "Feritin (overgrown) Day 6", "Feritin (overgrown) Day 8",
                                 "Feritin (WT looking) Day 6", "ImpL2 RNAi Day 6", "ImpL2 RNAi Day 8")
genotypesForAnalysisIDs <- list("RasYki_D5", "RasYki_D8", "Fer12OG_D6", "Fer12OG_D8", "Fer12WT_D6", "ImpL2i_D6", "ImpL2i_D8")
initallySelectedGenotypes <- list("RasYki_D5", "RasYki_D8", "Fer12OG_D6", "Fer12OG_D8", "Fer12WT_D6")
tissueNames <- list("Wing disc", "Salivary Gland", "Brain")
tissueValues <- list("Wingdisc", "Salivarygland", "Brain")
initallySelectedTissues <- list("Wingdisc", "Salivarygland", "Brain")

#Define UI ----
ui <- fluidPage(
  theme = shinytheme("united"),
  titlePanel("RNA-Seq Heatmap Generator"),
  sidebarLayout(
    sidebarPanel(h2("Configuration Tab:"),
                 checkboxGroupInput(inputId = "genotypeSelector",
                                    label = "Select genotypes for analysis:",
                                    choiceNames = genotypesForAnalysisNames,
                                    choiceValues = genotypesForAnalysisIDs,
                                    selected = initallySelectedGenotypes),
                 
                  checkboxGroupInput(inputId = "tissueSelector",
                                     label = "Select tissues for analysis:",
                                     choiceNames = tissueNames,
                                     choiceValues = tissueValues,
                                     selected = initallySelectedTissues),
                  textAreaInput(inputId = "geneList", 
                               label = "Enter Gene List as flybase ID's:", 
                               placeholder = "Fbgn0000123...",
                               value = NULL,
                               height = 200,
                               cols = 1)),
    mainPanel(
      h1("RNA-seq Heatmap:"),
      checkboxInput("convertGeneIDs", label = "Convert Flybase IDs to gene names.", value = FALSE),
      sliderInput("annotationWidth", label = "Annotation Width (cm):", value = 3, min = 0, max = 5, step = 0.5),
      sliderInput("scaleMinMax", label = "Select Min/Max values for colour scale:", value = c(-3,3), min = -10, max = 10),
      plotOutput(outputId = "mainHeatmap", height = "1000")
      )
  ))


#Define server logic ----
server <- function(input, output) {
  output$mainHeatmap <- renderPlot({
    
    #Set path to raw data folder and load in geneNames conversion sheet.
    rawData <- "ALL_Tissues_LFC_Database.xlsx"
    geneNames <- read_excel(rawData, sheet = 4)
    
    #Search User input for genes using REGEX.
    regFilter <- regex("FBGN\\d\\d\\d\\d\\d\\d\\d", ignore_case = TRUE, )
    geneList <- as.list(str_extract_all(input$geneList, regFilter))[[1]]
    #Select which tissues to use in heatmap based on user input.
    tissuesForHeatmap <- as.list(input$tissueSelector)
    #Select which genotypes to use in heatmap based on user input.
    genotypesForAnalysis <- as.list(input$genotypeSelector)
    #Select whether heatmap should be built with gene names or gene symbols.
    convertGeneIds <- input$convertGeneIDs[1]
    
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
    genotypesForAnalysis <- append(genotypesForAnalysis, c("LFC", "GeneName"))
    wdLfcTable <- wdLfcTable[, colnames(wdLfcTable) %in% genotypesForAnalysis]
    sgLfcTable <- sgLfcTable[, colnames(sgLfcTable) %in% genotypesForAnalysis]
    bLfcTable <- bLfcTable[, colnames(bLfcTable) %in% genotypesForAnalysis]
    
    #Select Which tissues to add to heatmap#
    Wingdisc <- wdLfcTable
    Salivarygland <- sgLfcTable
    Brain <- bLfcTable
    
    #Create Heatmap Color Scale.
    colorScale <- colorRamp2(c(input$scaleMinMax[1],0,input$scaleMinMax[2]), c("blue", "white", "red"))
    lgd <- Legend(col_fun = colorScale, title = "Log 2 Fold-Change")
    
    ##HEATMAP GENERATOR LOOP##
    heatmapList <- list()
    for (tissue in tissuesForHeatmap) {
      
      #Create normalised count level annotation.
      ncValues <- get(as.name(tissue))[, colnames(get(as.name(tissue))) == "LFC"]
      rowAnnotation <- rowAnnotation(LFC = anno_barplot(ncValues, width = unit(input$annotationWidth, "cm")), border = TRUE)
      
      #Create Core Heat map database.
      coreHeatmap <- get(as.name(tissue))[, colnames(get(as.name(tissue))) != "LFC"]
      
      if (convertGeneIds == TRUE){
        coreHeatmap <- rownames_to_column(coreHeatmap, var = "GeneID")
        coreHeatmap <- column_to_rownames(coreHeatmap, var = "GeneName")
        coreHeatmap$GeneID = NULL
      } else {
        coreHeatmap$GeneName = NULL
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
                         column_names_rot = 90,
                         heatmap_legend_param = list(title = "Log-2 Fold Change",
                                                     direction = "vertical",
                                                     title_position = "leftcenter-rot",
                                                     legend_height = unit(10, "cm")))
      heatmapList <- append(heatmap, heatmapList)
    }
    
    #Plot Heat map Based on Tissues Provided.
    HeatmapListLength <- length(heatmapList)
    
    if(HeatmapListLength == 1) {
      draw(heatmapList[[1]], heatmap_legend_side = "topleft")
    } else if (HeatmapListLength == 2) {
      draw(heatmapList[[2]] + heatmapList[[1]], heatmap_legend_side = "topleft")
    } else {
      draw(heatmapList[[3]] + heatmapList[[2]] + heatmapList[[1]])
    }
    
  })
}

#Run app ----
shinyApp(ui = ui, server = server)
