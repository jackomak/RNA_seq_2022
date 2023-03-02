#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#Libraries
library(shiny)
library(shinythemes)
library(DESeq2)
library(rsconnect)
library(BiocManager)

#Some option needed for deployment purposes.
option <- options(repos = BiocManager::repositories())

# Define UI for application ----
ui <- fluidPage(theme = shinytheme("flatly"), 
                navbarPage(title = "FHL RNA-seq", position = "static-top",
                tabPanel(title = "Normalisation & DEG Anaysis", NormalisationUI(id = "Normalisation")),
                tabPanel(title = "Heatmap Visualiser", HeatmapUI(id = "Heatmap")),
                tabPanel(title = "Individual Gene Analysis", IndividualGeneUI(id = "IndividualGene")),
                tabPanel(title = "Gene Enrichment"))
                )

# Define server logic ----
server <- function(input, output, session) {
  
  #Logic for DEG Expression tab.
  NormalisationServer(id = "Normalisation")
  
  #Logic for Heat map visualization tab.
  HeatmapServer(id = "Heatmap")
  
  #Logic for Individual gene visualisation tool.
  IndividualGeneServer(id = "IndividualGene")
}

# Run the application 
shinyApp(ui = ui, server = server)
