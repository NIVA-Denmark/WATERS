library(shiny)
library(DT)

# Define UI for random distribution application 
shinyUI(
  navbarPage(id = "inTabset",
    windowTitle="WATERS Status Assessment Tool",
    title=div(img(src="waters_2.gif") ), 
    #title= "WATERS Assessment", 
    tabPanel("Home",value="tabHome",
             h4("Introduction"),p("Explanation of the tool...")
    ),
    tabPanel("Get Data",value="tabData",
             navlistPanel(
               widths=c(2,10),well=F,
               tabPanel(
                 "New",
                 fluidRow(
                   column(3,
                          p(strong("New Assessment")),
                          #p("Select waterbodies and download data from database(s)."),
                          p("Select waterbodies"),
                          uiOutput("selectWaterDistrict"),
                          uiOutput("selectPeriod"),
                          uiOutput("selectWaterBodies"),
                          uiOutput("dataButton"),
                          uiOutput("nrows")
                          #actionButton("goButton", "Get data")
                          #conditionalPanel(condition = "input.waterbody != ''",)
                   ),
                   column(9, DT::dataTableOutput('x1'))
                 )
               ),
               tabPanel("Open",
                        fileInput('file1', 'Load existing assessment',
                                  accept=c('text/csv', 
                                           'text/comma-separated-values,text/plain', 
                                           '.csv'))  
               ),
               tabPanel("Save",
                        p(strong("Save the current assessment")),
                        downloadButton('downloadData', label="Save") 
               )
             )
             
    ),
    tabPanel("Options",value="tabOptions",
             navlistPanel(
               widths=c(2,10),well=F,
               tabPanel("Monte Carlo",
                        "Options for Monte Carlo simulations.",
                        numericInput("n",
                                     label = "Number of simulations", 
                                     value = 100)   
               ),
               tabPanel("Uncertainty",
                        "Uncertainty Library"),
               tabPanel("Aggregation",
                        "Aggregation methods")
             )
    ),
    
    tabPanel(
      "Results",value="tabResults",
      navlistPanel(
        widths=c(2,10),well=F,
        tabPanel("Overall",
                 "Overall Results",
                 uiOutput("resTableOverall")),
        tabPanel("Quality Element",
                 "QE Results",
                 uiOutput("resTableQE")),
        tabPanel("Indicator",
                 "Indicator Data",
                 uiOutput("resTableInd"))
      ),
      
      fluidRow(column(width=12,
                      downloadButton('downloadReport', label="Download report")
      ))
    )
  ))


