library(shiny)
library(DT)
library(PCRedux)
library(RDML)

options(shiny.maxRequestSize=10*1024^2)

source("functions.R", local = TRUE)

shinyServer(function(input, output) {
  
  raw_input <- reactive({
    withProgress({
      if (!is.null(input[["amp_file"]])) {
        input_curves <- switch(
          tolower(tools::file_ext(input[["amp_file"]][["datapath"]])),
          csv = read.csv(input[["amp_file"]][["datapath"]]),
          rdml = {
            rdml <- RDML$new(input[["amp_file"]][["datapath"]])
            rdml$GetFData()
          },
          xls = readxl::read_excel(input[["amp_file"]][["datapath"]],
                                   col_types = "numeric"),
          xlsx = readxl::read_excel(input[["amp_file"]][["datapath"]],
                                    col_types = "numeric"),
        )
      }
      input[["use_example"]]
      isolate({
        if (!is.null(input[["use_example"]]))
          if(input[["use_example"]] > 0)
            input_curves <- read.csv("example_dat.csv")
      })},
      message = "Reading file")
    
    if(exists("input_curves")) {
      # here we should put some tests of the input, like to maximum amount of curves
      # the first column should be named cycle
      colnames(input_curves)[1] <- "cycle"
      input_curves
    } else {
      NULL
    }
  })
  
  processed_input <- reactive({
    withProgress(shiny_encu(raw_input()),
                 message = "Processing data")
  })
  
  
  output[["raw_input_dt"]] <- DT::renderDataTable({
    my_DT(raw_input())
  })
  
  output[["processed_input_dt"]] <- DT::renderDataTable({
    formatRound(my_DT(processed_input()), 
                c(2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 12L, 13L, 14L, 15L, 16L, 
                  20L, 21L, 22L, 23L, 24L, 25L, 26L, 27L, 28L, 29L, 30L, 31L, 32L, 
                  33L, 34L, 35L, 36L, 37L, 45L, 46L, 47L, 48L, 49L),
                digits = 2)
  })
  
  
  output$dynamic_tabset <- renderUI({
    if(is.null(raw_input())) {
      tabsetPanel(
        tabPanel(title = "Data input",
                 fluidRow(
                   column(width = 5, 
                          fileInput('amp_file', 
                                    "Submit amplification data (.csv, .xls, .xlsx, .rdml):")
                   ),
                   column(width = 5,
                          p("Test the predPCR application with an example amplification data set"),
                          actionButton("use_example", "Load example")
                   ),
                   tags$style(type='text/css', "#use_example { width:100%; margin-top: -5px;}")
                 )
        ),
        tabPanel("About", includeMarkdown("readme.md"))
      )
    } else {
      tabsetPanel(
        tabPanel(title = "PCRedux data",
                 DT::dataTableOutput("processed_input_dt", width = "100%"),
                 tags$p(HTML("<button id=\"new_query\" type=\"button\" class=\"btn btn-default action-button\" onclick=\"javascript:history.go(0)\" style=\"width: 100%; margin: 20px;\"\">Start a new query</button>"))),
        tabPanel(title = "Input data",
                 DT::dataTableOutput("raw_input_dt")),
        tabPanel("About", includeMarkdown("readme.md"))
      )
    }
  })
})