#Libraries ----
library(shiny)
library(shinybusy)
library(DESeq2)


# Design page UI ----
NormalisationUI <- function(id){
  ns <- NS(id)
  
  tagList(
    titlePanel(title = "Normalisation and DEG Analysis:"),
    sidebarLayout(sidebarPanel(h3("Settings for Normalisation:"),
                               
                               #Add option to input raw count matrix as CSV file.
                               fileInput(inputId = ns("rawcounts_file"),
                                         label = "Select raw count matrix:",
                                         multiple = FALSE,
                                         accept = ".csv",
                                         buttonLabel = "Browse",
                                         placeholder = "No file selected."),
                               
                               #Add option to input sample info file as CSC file.
                               fileInput(inputId = ns("sampleinfo_file"),
                                         label = "Select sample information file:",
                                         multiple = FALSE,
                                         accept = ".csv",
                                         buttonLabel = "Browse",
                                         placeholder = "No file selected."),
                              
                                #Add option for user to change minimum count threshold
                               numericInput(inputId = ns("minCountCutoff"),
                                            label = "Minimum average expression cutoff:",
                                            min = 1,
                                            value = 30,
                                            width = "150px"),
                               
                               #Add mandatory settings for Differential gene expression.
                               h3("Settings for DEG Analysis:"),
                               numericInput(inputId = ns("pvalueThreshold"),
                                            label = ("P-value threshhold."),
                                            min = 0.000001,
                                            max = 1,
                                            value = 0.05,
                                            width = "150px"),
                               numericInput(inputId = ns("logFoldChangeThreshold"),
                                            label = "Log2Fold Change threshold (As positive number).",
                                            min = 0.001,
                                            max = 100,
                                            value = 1.4,
                                            width = "150px")
                               ),
                  
                  
                  
                  mainPanel(tabsetPanel(type = "tabs",
                                        tabPanel("Pre-normalisation Checks",
                                                 #Verify files:
                                                 h2(textOutput(ns("pre_normalisation_checker"))),
                                                 #Uploaded header of rawcounts file.
                                                 h3("Raw Counts File:"),
                                                 tableOutput(ns("rawcounts_file_table")),
                                                 #Upload sample Information file.
                                                 h3("Sample Information File:"),
                                                 tableOutput(ns("sampleinfo_file_table"))),
                                        
                                        

                                        #Create normalised count matrix using DESEQ2.
                                        tabPanel("Normalised Counts",
                                                 h3("Once Pre-normalisation check are complete:"),
                                                 selectInput(inputId = ns("experimental_Condition"), 
                                                             label = "Please select your experimental condtion:", choices = c(""), selected = ""),
                                                 actionButton(ns("generate_Normalised_Counts"), "Generate Normalised Count Table"),
                                                 tableOutput(outputId = ns("normalisedCountTable")),
                                                 add_busy_bar(color = "#FF0000"),
                                                 downloadButton(outputId = ns("downloadNormcounts"), label = "Download Normalised Count Matrix")),
                                        
                                        

                                        #Visualisation and Quality control of normalised counts.
                                        tabPanel("Normalisation Quality Control",
                                                 h3("Post Normalisation Quality Control:"),
                                                 p("Please make sure a normalised counts matrix has been generated
                                                   under the Normalised counts Tab prior to generating plots."),
                                                 h4("Generate PCA Plot:"),
                                                 actionButton(ns("Generate_PCA_Plot"), label = "Generate PCA Plot"),
                                                 plotOutput(outputId = ns("PCA_Plot")),
                                                 downloadButton(outputId = ns("downloadPCA"), label = "Save Plot")),
                                        tabPanel("Differential Gene Expression Analysis",
                                                 h3("Mandatory Inputs:"),
                                                 selectInput(inputId = ns("DEGControlGroup"),
                                                             label = "Please choose control group for comparison.",
                                                             choices = c("NONE FOUND"), selected = "NONE FOUND"),
                                                 p("Select a control group to use as the baseline for identifying 
                                                   differentially expressed genes"),
                                                 selectInput(inputId = ns("DEGTestGroup"), 
                                                             label = "Please choose test group for comparison.",
                                                             choices = "NONE FOUND", selected = "NONE FOUND"),
                                                 p("Select a group to identify differentially expressed genes for against the above control."),
                                                 selectInput(inputId = ns("DEGGeneSet"),
                                                             label = "Choose gene set to identify:",
                                                             choices = c("Upregulated Genes", "Downregulated Genes", "All DEG's"),
                                                             selected = ("All DEG's")),
                                                 actionButton(ns("generateDEGs"), label = "Generate DEG List"),
                                                 h3("Differentially Expressed Genes based on your criteria:"),
                                                 downloadButton(outputId = ns("downloadDEGs"), label = "Download as CSV"),
                                                 tableOutput(outputId = ns("DEGList"))
                                                 )
                                                 )
                                        )))
  
}

# Define server logic ----

NormalisationServer <- function(id) {
  moduleServer(id, function(input, output, session){
    
    #Output Top5 lines of users selected rawcount file.
    output$rawcounts_file_table <- renderTable({
      
      #Dont try to load file unless user adds it.
      req(input$rawcounts_file)
      
      #Output top 10 rows of dataframe to tab.
      head(read.csv(file = input$rawcounts_file$datapath, row.names = 1))
      
    }, rownames = TRUE, striped = TRUE, bordered = TRUE, spacing = "s")
    
    #Output sample information table logic.
    output$sampleinfo_file_table <- renderTable({
      
      #Dont try to load file unless user adds it.
      req(input$sampleinfo_file)
      
      #Output file as dataframe to tab.
      read.csv(file = input$sampleinfo_file$datapath, row.names = 1)
      
    }, rownames = TRUE, striped = TRUE, bordered = TRUE, spacing = "s")
    
    #If both files are uploaded, check to make sure both files are compatible for DEG analysis.
    output$pre_normalisation_checker <- renderText({
      
      req(input$rawcounts_file)
      req(input$sampleinfo_file)
      
      #Get colnames and rownames of relative files.
      rawcountColnames <- as.list(colnames(read.csv(file = input$rawcounts_file$datapath, row.names = 1)))
      sampleInfoRownames <- as.list(rownames(read.csv(file = input$sampleinfo_file$datapath, row.names = 1)))
      
      #Check they match.
      result = identical(rawcountColnames, sampleInfoRownames)
      
      #Inform user of outcome.
      if (result == TRUE){
        
      #Populate experimental condition widget for normalisation step.
      getExperimentalValues <- as.list(colnames(read.csv(file = input$sampleinfo_file$datapath, row.names = 1)))
      updateSelectInput(inputId = "experimental_Condition", choices = getExperimentalValues) 
      successMessage <- "Files compatible!"}
      
      else {
        errorMessage <- "Files incompatible, please check col/row names match."
      }})
  
    #Generate normalised count table if user specifies.
    observeEvent(input$generate_Normalised_Counts, {
      
      output$normalisedCountTable <<- renderTable({
        
        #Read in Files.
        rawcounts <- read.csv(file = input$rawcounts_file$datapath, header = TRUE, row.names = 1)
        sampleInfo <- read.csv(file = input$sampleinfo_file$datapath, header = TRUE, row.names = 1)
        
        #Create Deseq2 matrix datatype using user specified files.
        dds <<- DESeqDataSetFromMatrix(countData = rawcounts, colData = sampleInfo, design = as.formula(paste0("~",input$experimental_Condition)))
        
        #Populate Control and test group options for differential gene expression analysis.
        getExperimentalConditionValues <- as.list(sampleInfo[,paste0(input$experimental_Condition)])
        updateSelectInput(inputId = "DEGControlGroup", choices = getExperimentalConditionValues)
        updateSelectInput(inputId = "DEGTestGroup", choices = getExperimentalConditionValues)

        #Remove genes with a count below user specified low count threshold. 
        keep <- rowSums(counts(dds)) >= input$minCountCutoff
        ddsDE <<- DESeq(dds[keep,]) 
        
        #Generate normalised count matrix (save as higher order variable for access by other outputs).
        normcounts <<- as.data.frame(counts(ddsDE, normalized = TRUE))
        
        #Print first 100 rows of normalised count matrix.
        normcounts[1:100,]
        
        }, rownames = TRUE, bordered = TRUE, striped = TRUE) 

    })
    
    #Download Logic for downloading normcounts.
    output$downloadNormcounts <- downloadHandler(
      filename = function() {paste("Normalised_Count_Matrix_", Sys.Date(), ".csv", sep = "")},
      content = function(file){write.csv(normcounts, file)
      })
      
      
    #PCA Plot logic.
    observeEvent(input$Generate_PCA_Plot, {
      
      output$PCA_Plot <- renderPlot({
        vsd <- vst(dds, blind = F)
        PCA <<- plotPCA(vsd, intgroup = "GenotypeDay")
        plotPCA(vsd, intgroup = "GenotypeDay")
      })
    })
    
    
    #Download PCA Plot logic.
    output$downloadPCA <- downloadHandler(
      filename = function() {paste("PCA_Plot",Sys.Date(),".png", sep = "")},
      content = function(file) {
        png(file)
        print(PCA)
        dev.off()
    })
    
    #Generate Differential gene expression list.
    observeEvent(input$generateDEGs, {
      
      output$DEGList <- renderTable({
      
      #Get DEseq2 results for two way genotype comparison to identify DEG's.
      res <- results(ddsDE, contrast= c(paste0(input$experimental_Condition), 
                                        paste0(input$DEGTestGroup),
                                        paste0(input$DEGControlGroup)))
      statsDE <- res[order(res$padj),]
      filteredRes <- subset(statsDE, padj < input$pvalueThreshold)
      sigFoldGenesUp <- subset(filteredRes, log2FoldChange > input$logFoldChangeThreshold)
      sigFoldGenesDown <- subset(filteredRes, log2FoldChange < -input$logFoldChangeThreshold)

      if (input$DEGGeneSet == "Upregulated Genes"){finalDEGList <<- sigFoldGenesUp}
      if (input$DEGGeneSet == "Downregulated Genes"){finalDEGList <<- sigFoldGenesDown}
      else {finalDEGList <<- rbind(sigFoldGenesUp, sigFoldGenesDown)}
      
      finalDEGList
      
      }, rownames = TRUE, bordered = TRUE, striped = TRUE)
      
    })
    
    #Download handler for final list of differentially expressed genes.
    output$downloadDEGs <- downloadHandler(
                           filename = function() {paste0(input$DEGControlGroup,"_V_",input$DEGTestGroup,"_",Sys.Date(),".csv", sep = "")},
                           content = function(file){write.csv(finalDEGList, file)
    })
    
  })

}


