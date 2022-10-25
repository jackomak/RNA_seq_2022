#! - Rstudio / Shiny

#app.R - An interactive application to visualize log fold change in gene expression derived from data obtained from the 2022
#RNA-Seq cachexia project samples.

#Shiny packages ----
library(shiny)
library(shinythemes)
library(ggplot2)


#core packages for data visualization ----
library(ComplexHeatmap)
library(tidyverse)
library(readxl)
library(circlize)
library(gridExtra)
library(reshape2)


#Set core variable lists ----
rawData <- "ALL_Tissues_LFC_Database.xlsx"

#List of names to describe the availible genotypes for the heatmap in the GUI.
genotypesForAnalysisNames <- list("RasYki (D5)","RasYki (D8)","Feritin (D6)",
                                  "Feritin (D8)", "Feritin WT looking (D6)",
                                  "ImpL2 RNAi (D6)", "ImpL2 RNAi (D8)", "Wts (D6)",
                                  "Wts (D8)", "Cic Wts (D6)", "Cic Wts (D8)",
                                  "hh> RasScrib (D5-9)", "eyFLP > RasScrib (D5-9)")

#List of genotypes stored in the raw data excel file that describe column headings for each genotype.
genotypesForAnalysisIDs <- list("RasYki_D5","RasYki_D8","Fer12OG_D6",
                                "Fer12OG_D8","Fer12WT_D6","ImpL2i_D6",
                                "ImpL2i_D8", "WtsD6","WtsD9","CicWtsD6",
                                "CicWtsD9", "hh_RasScribD5_9", "_eyFLP_RasScrib_D5_9") 

#Select which genotypes to initially mark as checked on the GUI.
initallySelectedGenotypes <- list("RasYki_D5", "RasYki_D8", "Fer12OG_D6", "Fer12OG_D8")

#List of names to describe which heatmaps can be built and visualized on the GUI.
tissueNames <- list("FH 2022 Wing disc", "FH 2022 Salivary Gland", 
                    "FH 2022 Brain", "MA CicWts Wing Disc",
                    "MA RasScrib Wing Disc", "MA RasScrib Eye Disc")

#List of sheet names in the raw data excel file availible for heatmap creation.
tissueValues <- list("Wingdisc", "Salivarygland", "Brain", "Mardelle_CicWts_WD",
                     "Mardelle_hhRasScrib_WD", "Mardelle_EyFLPRasScrib_ED")

#Select which datasets to intially mark as checked by the GUI.
initallySelectedTissues <- list("Wingdisc", "Salivarygland", "Brain")
cat("Created UI Variables successfully.")


#Define UI ----
ui <- fluidPage(
  navbarPage("FH Lab RNA-seq Visualiser:", tabPanel("Heatmap Generator",
  

  #GUI for heatmap generator ----                                                 
  
  #Add Title to the GUI.
  titlePanel("RNA-Seq Heatmap Generator"),
  
  #Create Sidebar containing all of the following reactive elements:
  sidebarLayout(
    sidebarPanel(h2("Configuration Tab:"),
                 
                 #Checkbox to select which genotypes to add to the heatmap.
                 checkboxGroupInput(inputId = "genotypeSelector",
                                    label = "Select genotypes for analysis:",
                                    choiceNames = genotypesForAnalysisNames,
                                    choiceValues = genotypesForAnalysisIDs,
                                    selected = initallySelectedGenotypes),
                 
                 #Checkbox to select which datasets to include in the heatmap.
                  checkboxGroupInput(inputId = "tissueSelector",
                                     label = "Select tissues for analysis:",
                                     choiceNames = tissueNames,
                                     choiceValues = tissueValues,
                                     selected = initallySelectedTissues),
                 #Descriptor.
                 p("Please note - All LFC annotation data and information on some gene ID's is unavailible for the MA 2016 dataset."),
                 
                 #Box for user to copy Flybase geneID's they want to generate a heatmap for.
                 textAreaInput(inputId = "geneList", 
                               label = "Enter Gene List as flybase ID's:", 
                               placeholder = "Fbgn0000123...",
                               value = NULL,
                               height = 200,
                               cols = 1)),
    
    #What to add reactive elements to the main panel of the GUI.             
    mainPanel(
      h1("RNA-seq Heatmap:"),
      #Checkbox to convert default Flybase gene ID's to gene symbols/names.
      checkboxInput("convertGeneIDs", label = "Convert Flybase IDs to gene names.", value = FALSE),
      #Slider to add parameters to edit the width of the LFC annotation columns on the heatmap.
      sliderInput("annotationWidth", label = "Annotation Width (cm):", value = 3, min = 0, max = 5, step = 0.5),
      #Slider to add parameters to edit the boundaries of the color scale bar on the heatmap.
      sliderInput("scaleMinMax", label = "Select Min/Max values for colour scale:", value = c(-3,3), min = -10, max = 10),
      #Add the final heatmap image to the main panel. 
      plotOutput(outputId = "mainHeatmap", height = "1000")
      ))),
  
  
  #GUI For Gene Count Visualizer ----
  tabPanel("Gene Count Visualiser", 
           
           #Add title.
           titlePanel("Individual Gene Count Viewer"),
           
           #Create sidebar containing all of the following reactive elements.
           sidebarLayout(sidebarPanel(
             h2("Configuration Tab:"),
             
             #Checkbox to select for which Tissues to include.
             
             #Checkbox to select which gene to visualize.
             textAreaInput(inputId = "GeneVis_GeneName", 
                           label = "Enter a FLybase ID to visualise - 1 gene only.", 
                           placeholder = "Fbgn0000123",
                           value = NULL,
                           height = 40,
                           cols = 1),
             #Checkbox to convert gene Id's to gene names.
             
             #Checkbox to convert to a log scale.
             
             #Search By gene Name rather than symbol.
           ),
           
          #Design Main panel - Graph.
          mainPanel(h1("Gene Count Visualiser."), 
                    plotOutput(outputId = "GC_visualiser", height = "750"))),
      )),
  
  theme = shinytheme("united"))
  


#Define server logic ----

server <- function(input, output) {
  
  #Heatmap Generator logic ----
  output$mainHeatmap <- renderPlot({
    
    #Set path to raw data folder and load in geneNames conversion sheet.
    geneNames <- read_excel(rawData, sheet = 7)
    
    #Search User input for genes using REGEX. (Regular expression search.)
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
    genotypesForAnalysis <- append(genotypesForAnalysis, c("wt_NCL", "GeneName"))
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
      heatmap <- Heatmap(as.matrix(coreHeatmap),                                        #Select the core data to build heatmap with.
                         col = colorScale,                                              #Add the colorscale to the heatmap.
                         right_annotation = rowAnnotation,                              #Add the row annotation column to the heatmap (LNC).
                         column_title = tissue,                                         #Add column descriptors to the heatmap
                         cluster_columns = FALSE,                                       #Do/do not cluster columns based on genotype.
                         column_title_side = "top",                                     #Choose side to add column titles to.
                         column_names_rot = 90,                                         #Rotate columns in degrees.
                         row_title = "Gene",                                            #Row descriptor.
                         row_title_gp = gpar(col = "white"),                            #Change color of row title.
                         row_title_side = "right",                                      #Choose side to add row title.
                         row_gap = unit(1, "mm"),                                       #Change size of gaps between rows of heatmap.
                         name = "Log2Fold Change",                                      #Title of heatmap.
                         border = TRUE,                                                 #Add border to heatmap cells.
                         heatmap_legend_param = list(title = "Log-2 Fold Change",       #Change parameters of heatmap legend
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
    } else if (HeatmapListLength == 3) {
      draw(heatmapList[[3]] + heatmapList[[2]] + heatmapList[[1]], heatmap_legend_side = "right")
    } else {
      draw(heatmapList[[4]] + heatmapList[[3]] + heatmapList[[2]] + heatmapList[[1]], heatmap_legend_side = "right")
    }}
    heatmapBuilder(HeatmapListLength, heatmapList)
    
  })
  
  # GC Visualizer logic ----
  output$GC_visualiser <- renderPlot({
   
  #Set empty variables.
  normalisedCountsList<- list()
  geneCountVisualsList <- list()
  valuesForScale <- c()
 
  #Set color scale for each genotype.
  geneotypeColors <- c("PtcG4_D6" = "#00BA38", "Yw_D5" = "#00BA38",
                      "RasYki_D5    " = "#F8766D", "RasYki_D8" = "#E68613",
                      "Fer12OG_D6" = "#00B8E7", "Fer12OG_D8" = "#00A9FF", "Fer12WT_D6" = "#8494FF",
                      "ImpL2i_D6" = "#ED68ED", "ImpL2i_D8" = "#FF61CC")
 
  #Select gene to analyse using flybase ID. Can add functionality later to convert this to gene Name.
  GeneID <- input$GeneVis_GeneName
 
  #Logic to assign the max normcount value from each of the three datasets for the gene of interest for y scale of graphs.
  for (i in c("Wing Disc", "Salivary Gland", "Brain")){
   
    #Read in databases from normalised count excel file. Each tissue dataset is stored on a seperate sheet.
    normalisedCountsDatabase <- read_excel("All_Tissue_Normcounts.xlsx", sheet = i) %>% column_to_rownames("...1")
   
    #Append the database to a list.
    normalisedCountsList[[i]] <- normalisedCountsDatabase
   
    #Grab Info for the gene of interest.
    geneInfo <- melt(as.matrix(normalisedCountsDatabase[GeneID,]))
    valuesForScale <- append(valuesForScale, geneInfo$value)
  }
 
  #Logic to generate matrix tables for each of the three tissue types and generate graphs for each. Outputs a list contianing all of the figures.
  for (i in c("Wing Disc", "Salivary Gland", "Brain")){
   
    #Grab the relevent database generated by the previous for loop.
    normalisedCountsDatabase <- normalisedCountsList[[i]]
   
    #Transform the matrix so graph can be generated using different genotypes.
    geneInfo <- melt(as.matrix(normalisedCountsDatabase[GeneID,]))
   
    #Rename genotypes.
    geneInfo$Genotype <- ifelse(grepl("PtcG4_D6", geneInfo$Var2), "PtcG4_D6",
                          ifelse(grepl("Yw_D5", geneInfo$Var2), "Yw_D5",
                          ifelse(grepl("RasYki_D5", geneInfo$Var2), "RasYki_D5    ", 
                          ifelse(grepl("RasYki_D8", geneInfo$Var2), "RasYki_D8",
                          ifelse(grepl("Fer12OG_D6", geneInfo$Var2), "Fer12OG_D6",
                          ifelse(grepl("Fer12OG_D8", geneInfo$Var2), "Fer12OG_D8",
                          ifelse(grepl("Fer12WT_D6", geneInfo$Var2), "Fer12WT_D6",
                          ifelse(grepl("ImpL2i_D6", geneInfo$Var2), "ImpL2i_D6", "ImpL2i_D8"))))))))
   
    #Reorder Levels of genotype factor.
    if (i == "Brain"){
      geneInfo$Genotype <- factor(geneInfo$Genotype, levels = c("Yw_D5", "RasYki_D5    ", "RasYki_D8", "ImpL2i_D6", "ImpL2i_D8"))
    }   
    else {
      geneInfo$Genotype <- factor(geneInfo$Genotype, levels = c("PtcG4_D6", "RasYki_D5    ", "RasYki_D8", "Fer12OG_D6", "Fer12OG_D8", "Fer12WT_D6", "ImpL2i_D6", "ImpL2i_D8"))
    }
   
    #Generate Graph for each tissue type and append to geneCountVisualsList.
    geneCountsVisual <- ggplot(geneInfo, aes(x = Genotype, y = value, color = Genotype)) +
      geom_point(size = 4) +
      facet_grid(~Var1) +
      theme(axis.text.x = element_text(angle = 90, size = 15, hjust = 0), 
            axis.text.y = element_text(angle = 0, size = 15),
            plot.title = element_text(hjust = 0.5, size = 20),
            axis.title = element_text(size = 15),
            legend.position = "none",
            ) +
      ggtitle(paste0(i)) +
      ylim(0, max(valuesForScale)) +
      ylab("Normalised Count") +
      scale_color_manual(values = geneotypeColors)
   
    geneCountVisualsList[[i]] <- geneCountsVisual
    
  }
    
    #Combine graphs generated from each tissue dataset to generate final figure.
    final_figure <- grid.arrange(geneCountVisualsList[[1]], geneCountVisualsList[[2]], geneCountVisualsList[[3]], nrow = 1)
    grid.draw(final_figure)

    })

  }


#Run app ----
shinyApp(ui = ui, server = server)
