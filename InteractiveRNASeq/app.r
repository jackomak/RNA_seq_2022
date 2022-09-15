#Shiny packages ----
library(shiny)
library(shinythemes)

#core packages for data visualization ----
library(ComplexHeatmap)
library(tidyverse)
library(readxl)
library(circlize)

#Define UI ----
ui <- fluidPage(
  theme = shinytheme("united"),
  titlePanel("RNA-Seq Heatmap Generator"),
  sidebarLayout(
    sidebarPanel(h2("Configuration Tab:"),
                 checkboxGroupInput(inputId = "genotypeSelector",
                                    label = "Select genotypes for analysis:",
                                    choiceNames = list("RasYki Day 5",
                                                       "RasYki Day 8",
                                                       "Feritin (overgrown) Day 6",
                                                       "Feritin (overgrown) Day 8",
                                                       "Feritin (WT looking) Day 6",
                                                       "ImpL2 RNAi Day 6",
                                                       "ImpL2 RNAi Day 8"),
                                    choiceValues = list("RasYki_D5",
                                                        "RasYki_D8",
                                                        "Fer12OG_D6",
                                                        "Fer12OG_D8",
                                                        "Fer12WT_D6",
                                                        "ImpL2i_D6",
                                                        "ImpL2i_D8"),
                                    selected = list("RasYki_D5",
                                                    "RasYki_D8",
                                                    "Fer12OG_D6",
                                                    "Fer12OG_D8",
                                                    "Fer12WT_D6")),
                  checkboxGroupInput(inputId = "tissueSelector",
                                     label = "Select tissues for analysis:",
                                     choiceNames = list("Wing disc",
                                                        "Salivary Gland",
                                                        "Brain"),
                                     choiceValues = list("Wingdisc",
                                                         "Salivarygland",
                                                         "Brain"),
                                     selected = list("Wingdisc",
                                                     "Salivarygland",
                                                     "Brain")),
                  textAreaInput(inputId = "geneList", 
                               label = "Enter Gene List as flybase ID's:", 
                               placeholder = "Fbgn0000123...",
                               value = NULL,
                               height = 200,
                               cols = 1)),
    mainPanel(
      h1("RNA-seq Heatmap:"),
      checkboxInput("convertGeneIDs", label = "Convert Flybase IDs to gene names.", value = FALSE),
      plotOutput(outputId = "mainHeatmap", height = "1000")),
  ))


#Define server logic ----
server <- function(input, output) {
  output$mainHeatmap <- renderPlot({
    
    #Set path to raw data folder and load in geneNames conversion sheet.
    rawData <- "ALL_Tissues_LFC_Database.xlsx"
    geneNames <- read_excel(rawData, sheet = 4)
    
    #Search User input for genes using REGEX.
    regFilter <- regex("FBgn\\d\\d\\d\\d\\d\\d\\d", ignore_case = TRUE, )
    geneList <- as.list(str_extract_all(input$geneList, regFilter))[[1]]
    #Select which tissues to use in heatmap based on user input.
    tissuesForHeatmap <- as.list(input$tissueSelector)
    #Select which genotypes to use in heatmap based on user input.
    genotypesForAnalysis <- as.list(input$genotypeSelector)
    #Select whether heatmap should be built with gene names or gene symbols.
    convertGeneIds <- input$convertGeneIDs
    
    #Read in raw data from external excel file.
    wdLfcTable <- read_excel(rawData, sheet = 1)
    wdLfcTable <- column_to_rownames(wdLfcTable, var = "GeneID")
    
    sgLfcTable <- read_excel(rawData, sheet = 2)
    sgLfcTable <- column_to_rownames(sgLfcTable, var = "GeneID")
    
    bLfcTable <- read_excel(rawData, sheet = 3)
    bLfcTable <- column_to_rownames(bLfcTable, var = "GeneID")
    
    #Filter data sets to only include rows that user has specified in "genelist" variable.
    wdLfcTable <- wdLfcTable[rownames(wdLfcTable) %in% geneList, ]
    sgLfcTable <- sgLfcTable[rownames(sgLfcTable) %in% geneList, ]
    bLfcTable <- bLfcTable[rownames(bLfcTable) %in% geneList, ]
    
    #Remove unwanted tissue datasets#
    genotypesForAnalysis <- append(genotypesForAnalysis, "LFC")
    wdLfcTable <- wdLfcTable[, colnames(wdLfcTable) %in% genotypesForAnalysis]
    sgLfcTable <- sgLfcTable[, colnames(sgLfcTable) %in% genotypesForAnalysis]
    bLfcTable <- bLfcTable[, colnames(bLfcTable) %in% genotypesForAnalysis]
    
    #Select Which tissues to add to heatmap#
    Wingdisc <- wdLfcTable
    Salivarygland <- sgLfcTable
    Brain <- bLfcTable
    
    #Create heatmap Color Scale.
    colorScale <- colorRamp2(c(-3,0,3), c("blue", "white", "red"))
    
    #Loop to create heatmap for each tissue.
    heatmapList <- list()
    for (tissue in tissuesForHeatmap) {
      
      #Create normalised count level annotation.
      ncValues <- get(as.name(tissue))[, colnames(get(as.name(tissue))) == "LFC"]
      rowAnnotation <- rowAnnotation(LFC = anno_barplot(ncValues), border = TRUE)
      
      #Create core heatmap database.
      coreHeatmap <- get(as.name(tissue))[, colnames(get(as.name(tissue))) != "LFC"]
      
      #Change gene IDs to gene names if user selects.
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
    
    #Plot heatmap Based on tissues provided.
    HeatmapListLength <- length(heatmapList)
    
    if(HeatmapListLength == 1) {
      plot(heatmapList[[1]])
    } else if (HeatmapListLength == 2) {
      plot(heatmapList[[2]] + heatmapList[[1]])
    } else {
      plot(heatmapList[[3]] + heatmapList[[2]] + heatmapList[[1]])
    }
    
  })
  
}

#Run app ----
shinyApp(ui = ui, server = server)
