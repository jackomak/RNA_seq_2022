library(shiny)

#Define UI ----
ui <- fluidPage(
  theme = "Install Shiny theme and select package", # Complete
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
                                                        "Fer12OGD6",
                                                        "Fer12OGD8",
                                                        "Fer12WTD6",
                                                        "ImpL2iD6",
                                                        "ImpL2iD8")),
                  checkboxGroupInput(inputId = "tissueSelector",
                                     label = "Select tissues for analysis:",
                                     choiceNames = list("Wing disc",
                                                        "Salivary Gland",
                                                        "Brain"),
                                     choiceValues = list("Wingdisc",
                                                         "Salivarygland",
                                                         "Brain")),
                 actionButton(inputId = "Generate",
                              label = "Generate")),
    mainPanel(
      h1("I am a header"),
      h3("I am a smaller header"))
  ))


#Define server logic ----
server <- function(input, output) {
  
  
  
  
  

}

#Run app ----
shinyApp(ui = ui, server = server)