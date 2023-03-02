#UI and backend for individual gene count viewer tab.

#Libraries
library(ggplot2)
library(shiny)
library(dqshiny)
library(readxl)
library(tidyverse)
library(ggpubr)

#Source heatmap_string file - contains hardcoded function needed to reset genotype IDs during plot generation.
source("./R/HeatmapStrings.R")

#Read Raw data.
geneIDs <- read_excel("./Data/All_Tissue_Normcounts.xlsx", sheet = 1)



#Design UI ----
IndividualGeneUI <- function(id){
  ns <- NS(id)
  tagList(
    titlePanel("Individual Gene Explorer:"),
    sidebarLayout(sidebarPanel(h3("Explorer Settings:"),
                               autocomplete_input(id = ns("geneToView"), 
                                                  options = as.list(geneIDs$...1),
                                                  label = "Choose a Gene:",
                                                  placeholder = "FBgn0031085",
                                                  max_options = 10),
                               checkboxInput(inputId = ns("convertToLogTransformed"),
                                             label = "Convert Y axis to log format:",
                                             value = FALSE),
                               checkboxInput(inputId = ns("normalisePlotAxes"),
                                             label = "Normalize Plot Axes:",
                                             value = TRUE)
                               
  ),
  
  mainPanel(h1("Normalised count comparisons across tissues:"),
            plotOutput(outputId = ns("IndividualGeneCountGraph"), height = 750))
  
   
      ))
  
}
  




IndividualGeneServer <- function(id) {
  moduleServer(id, function(input, output, session){
      
               
               output$IndividualGeneCountGraph <- renderPlot({
                 
                 #Require a gene to be selected.
                 req(input$geneToView)
                 
                 #Create empty list to store genecounts.
                 geneCountList <- list()
                 
                 #Create MaxValue for y lim assessment placeholder (0).
                 maxNormalisedCountVal <- 0
                 
                 #For each Tissue/Sheet in normalised count file. Extract the normalised counts for the gene on interest.
                   for (i in 1:3){
                   
                   #Read in gene data for each tissue and index by geneID.
                   tissueInfo<- read_excel("./Data/All_Tissue_Normcounts.xlsx", sheet = i)
                   colnames(tissueInfo)[1] <- "geneID"
                   tissueInfo <- column_to_rownames(tissueInfo, "geneID")
                   
                   #Grab row user is interested in
                   geneInfo <- melt(as.matrix(tissueInfo[input$geneToView, ]))
                   colnames(geneInfo) <- c("Gene_ID", "Genotype", "Normcount")
                   
                   #If user selects, log transform normcounts.
                   if (input$convertToLogTransformed == TRUE){
                    geneInfo$Normcount <- log(geneInfo$Normcount+1 , base = 2)}
          
                   
                   #Format genotype column cells into readable format.
                   geneInfo$Genotype <- formatIndividualGenecells(geneInfo)
                   
                   #Append tissue gene count to list.
                   geneCountList <- append(geneCountList, values = list(geneInfo))
                   
                   #If user chooses to normalise plot scales. Find the value to set plot Y lims to.
                   if (input$normalisePlotAxes == TRUE) {
                     maxval <- max(geneInfo$Normcount, na.rm = TRUE)
                     
                      if (maxval > maxNormalisedCountVal){
                        maxNormalisedCountVal <- maxval
                      }
                    
                 }
                 
                 
                 #Generate facet grid plot for each organ and store as list.
                 facetPlotList <- list()
                 facetPlotTitleList <- c("Wing Disc", "Salivary Gland", "Brain")
                 counter <- 0
                 
                 for (tissue in geneCountList){
                   
                   #Get Core plot data.
                   coreData <- as.data.frame(tissue)
                   
                   #Get correct fact levels for each organ.
                   counter <- counter + 1
                   
                   if (counter == 3){plotLevels <- brainLevels}
                   else {plotLevels <- WDSGLevels}
                   
                   #Get correct plot title.
                   plotTitle <- facetPlotTitleList[1]
                   facetPlotTitleList <- facetPlotTitleList[-1]
                   
                   #Get Correct Y Label 
                   if (input$convertToLogTransformed == TRUE){
                     ylabel <- "Log2(Normalised Count +1)"}
                   else { ylabel <- " Normalised Count"}
                   
                   facetPlot <- ggplot(coreData, aes(x = Genotype, y = Normcount, color = Genotype)) +
                     geom_point(size= 4) +
                     facet_grid(~Gene_ID) +
                     ggtitle(plotTitle) +
                     scale_x_discrete(limits = plotLevels) +
                     scale_y_continuous(limits = if(input$normalisePlotAxes == TRUE){c(0, maxNormalisedCountVal)}) +
                     ylab(ylabel) +
                     theme(axis.text.x = element_text(size = 18, hjust = 1, angle = 45),
                           axis.text.y = element_text(size = 18), 
                           axis.title = element_text(size = 18),
                           strip.text.x = element_text(size = 18),
                           plot.title = element_text(size = 18, hjust = 0.5),
                           legend.position = "none")
                     
                 
                    #Append plot to list of plots for ggpubr.
                    facetPlotList <- c(facetPlotList, list(facetPlot))
                 

                    }
                 }
                
                #Append all elements into a single, readable figure.
                ggarrange(facetPlotList[[1]], facetPlotList[[2]], facetPlotList[[3]], ncol = 3)
                
               })}
               
              )}
   