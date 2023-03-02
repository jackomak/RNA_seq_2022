#Libraries
library(shiny)
library(readxl)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(gridExtra)
library(reshape2)

#Source String objects for UI:
source("./R/HeatmapStrings.R", local = TRUE)

# Design page UI ----

HeatmapUI <- function(id){
  ns <- NS(id)
  tagList(
    titlePanel(title = "Heatmap Visualisation Tool:"),
    sidebarLayout(sidebarPanel(h3("Mandatory Settings:"),
                               
                               #Checkbox group to select which genotypes to include in the heatmap.
                               checkboxGroupInput(inputId = ns("genotypeSelector"),
                                                  label = "Select genotypes for analysis:",
                                                  choiceNames = genotypesForAnalysisNames,
                                                  choiceValues = genotypesForAnalysisIDs,
                                                  selected = initallySelectedGenotypes),
                               
                               #Checkbox to select which datasets to include in the heatmap.
                               checkboxGroupInput(inputId = ns("tissueSelector"),
                                                  label = "Select tissues for analysis:",
                                                  choiceNames = tissueNames,
                                                  choiceValues = tissueValues,
                                                  selected = initallySelectedTissues),
                               
                               #Box for user to copy Flybase geneID's they want to generate a heatmap for.
                               textAreaInput(inputId = ns("geneList"), 
                                             label = "Enter Gene List as flybase ID's:", 
                                             placeholder = "Fbgn0000123...",
                                             value = NULL,
                                             height = 200,
                                             cols = 1),
                               
                               #Optional Paramaters.
                               h4("Optional parameters for heatmap aesthetics:"),
                               #Checkbox to convert gene symbols to gene names on heatmap.
                               checkboxInput(inputId = ns("convertGeneIDs"), 
                                             label = "Convert Flybase IDs to gene names.",
                                             value = TRUE),
                               sliderInput(inputId = ns("annotationWidth"),
                                           label = "Annotation Width (cm):", 
                                           value = 1.5, min = 0, max = 5, step = 0.5),
                               #Slider to add parameters to edit the boundaries of the color scale bar on the heatmap.
                               sliderInput(inputId = ns("scaleMinMax"),
                                           label = "Select Min/Max values for colour scale:",
                                           value = c(-3,3), min = -10, max = 10)),
                  
                  #Main Panel display for core heatmap.
                  mainPanel(h2("Heatmap for selected genes:"),
                            downloadButton(outputId = ns("downloadHeatmap"), label = "Download Heatmap as png"),
                            plotOutput(outputId = ns("mainHeatmap"), height = 1000))
                            ))
  
} 



# Define server logic ----

HeatmapServer <- function(id) {
  moduleServer(id, function(input, output, session){
    
    output$mainHeatmap <- renderPlot({
      
      #Only run if user adds gene ID list.
      req(input$geneList)
      
      #Get Gene names for gene ID conversion option.
      geneNames <- read_excel(rawdataPath, sheet = 4)
      
      #Search User input for genes using REGEX. (Regular expression search.)
      regFilter <- regex("FBGN\\d\\d\\d\\d\\d\\d\\d", ignore_case = TRUE, )
      geneList <- as.list(str_extract_all(input$geneList, regFilter))[[1]]
      
      #Select which tissues to use in heatmap based on user input.
      tissuesForHeatmap <- as.list(input$tissueSelector)
      
      #Select which genotypes to use in heatmap based on user input.
      genotypesForAnalysis <- as.list(input$genotypeSelector)
      
      #Select whether heatmap should be built with gene names or gene symbols.
      convertGeneIds <- input$convertGeneIDs[1]
    
      #Function to get gene names and replace geneIDs as rownames for core datasets.
      formatLFCTable <- function(sheetnumber, inputVar) {
        inputVar <- read_excel(rawdataPath, sheet = sheetnumber) %>%
        left_join(geneNames, by = "GeneID") %>%
        column_to_rownames(var = "GeneID")
        
        inputVar <- inputVar[rownames(inputVar) %in% geneList,]
      }
        
      #Read in Core Data for each organ.
      wdLfcTable <- formatLFCTable(sheetnumber = 1, inputVar =  wdLfcTable)
      sgLfcTable <- formatLFCTable(sheetnumber = 2, inputVar = sgLfcTable)
      bLfcTable <- formatLFCTable(sheetnumber =  3, inputVar =  bLfcTable)
      
      #Remove unwanted tissue datasets.
      genotypesForAnalysis <- append(genotypesForAnalysis, c("wt_NCL", "GeneName"))
      listOfGeneotypesToFilter <- c(wdLfcTable, sgLfcTable, bLfcTable)
      
      #Function to remove unwanted tissue/genotypes from the dataset.
      genotypeFilter <- function(inputVar) {
        inputVar <<- inputVar[, colnames(inputVar) %in% genotypesForAnalysis]}
      
      #Run Function through each tissue dataset - is there a way to compact this? - mapply() maybe?
      wdLfcTable <- genotypeFilter(wdLfcTable)
      sgLfcTable <- genotypeFilter(sgLfcTable)
      bLfcTable <- genotypeFilter(bLfcTable)
      
      #Select Which tissues to add to heatmap.
      Wingdisc <- wdLfcTable
      Salivarygland <- sgLfcTable
      Brain <- bLfcTable
      
      #Create Heatmap Color Scale.
      colorScale <- colorRamp2(c(input$scaleMinMax[1],0,input$scaleMinMax[2]), c("#1691D3", "#FFFFFF", "#E20F0F"))
      lgd <- Legend(col_fun = colorScale, title = "Log 2 Fold-Change")
      
      ##HEATMAP GENERATOR LOOP## -- Generates a heatmap for each tissue above and stores them all in a list variable.
      heatmapList <- list()
      for (tissue in tissuesForHeatmap) {
        
        #Create normalised count level annotation (LNC).
        ncValues <- get(as.name(tissue))[, colnames(get(as.name(tissue))) == "wt_NCL"]
        rowAnnotation <- rowAnnotation(wt_NCL = anno_barplot(ncValues, width = unit(input$annotationWidth, "cm")), border = FALSE)
        
        #Create Core Heat map database.
        coreHeatmap <- get(as.name(tissue))[, colnames(get(as.name(tissue))) != "wt_NCL"]
        
        #If user selects gene names, replace rownames with gene names.
        if (convertGeneIds == TRUE){
          coreHeatmap <- rownames_to_column(coreHeatmap, var = "GeneID")
          coreHeatmap <- column_to_rownames(coreHeatmap, var = "GeneName")
          coreHeatmap$GeneID = NULL
          
          #Else, move on to build heatmap with flybase geneID's as rownames.
        } else {
          coreHeatmap$GeneName = NULL
        }
        
        #Plot heatmap/s.
        heatmap <- Heatmap(as.matrix(coreHeatmap),                                        
                           col = colorScale,                                              
                           right_annotation = rowAnnotation,                              
                           column_title = tissue,                                         
                           cluster_columns = FALSE,                                       
                           column_title_side = "top",                                    
                           column_names_rot = ,                                         
                           row_title = "Gene",                                            
                           row_title_gp = gpar(col = "white"),                            
                           row_title_side = "right",                                      
                           row_gap = unit(1, "mm"),                                       
                           name = "Log2Fold Change",                                     
                           border = TRUE,                                                 
                           heatmap_legend_param = list(title = "Log-2 Fold Change",       
                                                       direction = "vertical",
                                                       title_position = "leftcenter-rot",
                                                       legend_height = unit(10, "cm")))
        #Add each generated heatmap to a list.
        heatmapList <- append(heatmap, heatmapList)
      }
      
      #Create figure using each heatmap stored in the list using the heatmapBuilder() custom function.
      HeatmapListLength <- length(heatmapList)
      heatmapBuilder <- function(HeatmapListLength, heatmapList) {
        
        if(HeatmapListLength == 1) {
          draw(heatmapList[[1]], heatmap_legend_side = "right")
        } else if (HeatmapListLength == 2) {
          draw(heatmapList[[2]] + heatmapList[[1]], heatmap_legend_side = "right")
        } else {
          draw(heatmapList[[3]] + heatmapList[[2]] + heatmapList[[1]], heatmap_legend_side = "right")
        }}
      
      #Plot Final Heatmap
      finalHeatmap <<- heatmapBuilder(HeatmapListLength, heatmapList)
      
      })
    
    output$downloadHeatmap <- downloadHandler(
      
      filename = function() {paste("Heatmap_Plot_",Sys.Date(),".png", sep = "")},
      content = function(file) {
        png(file)
        print(finalHeatmap)
        dev.off()
      })
  
    })
  
}

    
    
  