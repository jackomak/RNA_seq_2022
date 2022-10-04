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
rawData <- "ALL_Tissues_LFC_Database.xlsx"

genotypesForAnalysisNames <- list("RasYki (D5)","RasYki (D8)","Feritin (D6)",
                                  "Feritin (D8)", "Feritin WT looking (D6)",
                                  "ImpL2 RNAi (D6)", "ImpL2 RNAi (D8)", "Wts (D6)",
                                  "Wts (D8)", "Cic Wts (D6)", "Cic Wts (D8)",
                                  "hh> RasScrib (D5-9)", "eyFLP > RasScrib (D5-9)")

genotypesForAnalysisIDs <- list("RasYki_D5","RasYki_D8","Fer12OG_D6",
                                "Fer12OG_D8","Fer12WT_D6","ImpL2i_D6",
                                "ImpL2i_D8", "WtsD6","WtsD9","CicWtsD6",
                                "CicWtsD9", "hh_RasScribD5_9", "_eyFLP_RasScrib_D5_9") 

initallySelectedGenotypes <- list("RasYki_D5", "RasYki_D8", "Fer12OG_D6", "Fer12OG_D8")

tissueNames <- list("FH 2022 Wing disc", "FH 2022 Salivary Gland", 
                    "FH 2022 Brain", "MA CicWts Wing Disc",
                    "MA RasScrib Wing Disc", "MA RasScrib Eye Disc")

tissueValues <- list("Wingdisc", "Salivarygland", "Brain", "Mardelle_CicWts_WD",
                     "Mardelle_hhRasScrib_WD", "Mardelle_EyFLPRasScrib_ED")

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
                 
                 p("Please note - All LFC annotation data and information on some gene ID's is unavailible for the MA 2016 dataset."),
                 
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
    geneNames <- read_excel(rawData, sheet = 7)
    
    #Search User input for genes using REGEX.
    regFilter <- regex("FBGN\\d\\d\\d\\d\\d\\d\\d", ignore_case = TRUE, )
    geneList <- as.list(str_extract_all(input$geneList, regFilter))[[1]]
    #Select which tissues to use in heatmap based on user input.
    tissuesForHeatmap <- as.list(input$tissueSelector)
    #Select which genotypes to use in heatmap based on user input.
    genotypesForAnalysis <- as.list(input$genotypeSelector)
    #Select whether heatmap should be built with gene names or gene symbols.
    convertGeneIds <- input$convertGeneIDs[1]
    
    #Read in raw log fold change databases for each tissue. # Create function to do this.
    formatLFCTable <- function(sheetnumber, inputVar) {
      inputVar <- read_excel(rawData, sheet = sheetnumber)
      inputVar <- left_join(inputVar, geneNames, by = "GeneID")
      inputVar <- column_to_rownames(inputVar, var = "GeneID")
      inputVar <- inputVar[rownames(inputVar) %in% geneList, ]
    }
    
    ####Could reduce redundancy by using sheet name instead of sheet number.
    wdLfcTable <- formatLFCTable(sheetnumber = 1, inputVar =  wdLfcTable)
    sgLfcTable <- formatLFCTable(sheetnumber = 2, inputVar = sgLfcTable)
    bLfcTable <- formatLFCTable(sheetnumber =  3, inputVar =  bLfcTable)
    MaWdLfcTable <- formatLFCTable(sheetnumber =  4, inputVar = MaWdLfcTable)
    MaHhWdLfcTable <- formatLFCTable(sheetnumber = 5, inputVar = MaHhWdLfcTable)
    MaEyFLPLfcTable <- formatLFCTable(sheetnumber = 6, inputVar = MaEyFLPLfcTable)
    
    #Remove unwanted tissue datasets#
    genotypesForAnalysis <- append(genotypesForAnalysis, c("LFC", "GeneName"))
    listOfGeneotypesToFilter <- c(wdLfcTable, sgLfcTable, bLfcTable, MaWdLfcTable, MaHhWdLfcTable, MaEyFLPLfcTable)
    
    #Function to remove unwanted tissue/genotypes from the dataset.
    genotypeFilter <- function(inputVar) {
      inputVar <<- inputVar[, colnames(inputVar) %in% genotypesForAnalysis]
    }
    
    #Run Function through each genotype - is there a way to de recurse this? loop through list an output to global env?
    wdLfcTable <- genotypeFilter(wdLfcTable)
    sgLfcTable <- genotypeFilter(sgLfcTable)
    bLfcTable <- genotypeFilter(bLfcTable)
    MaWdLfcTable <- genotypeFilter(MaWdLfcTable)
    MaHhWdLfcTable <- genotypeFilter(MaHhWdLfcTable)
    MaEyFLPLfcTable <- genotypeFilter(MaEyFLPLfcTable)
    
    #Select Which tissues to add to heatmap#
    Wingdisc <- wdLfcTable
    Salivarygland <- sgLfcTable
    Brain <- bLfcTable
    Mardelle_CicWts_WD <- MaWdLfcTable
    Mardelle_hhRasScrib_WD <- MaHhWdLfcTable
    Mardelle_EyFLPRasScrib_ED <- MaEyFLPLfcTable
    
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
    heatmapBuilder <- function(HeatmapListLength, heatmapList) {
    
    if(HeatmapListLength == 1) {
      draw(heatmapList[[1]], heatmap_legend_side = "right")
    } else if (HeatmapListLength == 2) {
      draw(heatmapList[[2]] + heatmapList[[1]], heatmap_legend_side = "right")
    } else if (HeatmapListLength == 3) {
      draw(heatmapList[[3]] + heatmapList[[2]] + heatmapList[[1]], heatmap_legend_side = "right")
    } else {
      draw(heatmapList[[4]] + heatmapList[[3]] + heatmapList[[2]] + heatmapList[[1]], heatmap_legend_side = "right")
    }}
    heatmapBuilder(HeatmapListLength, heatmapList)
    
  })
}

#Run app ----
shinyApp(ui = ui, server = server)
