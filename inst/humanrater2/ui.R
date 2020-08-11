library(shiny)
library(shinycssloaders)

shinyUI(fluidPage(
  

  titlePanel("Human Rater 2.0"),
  

  sidebarLayout(
    sidebarPanel(
      fileInput("amp_file", "Data file"),
      actionButton("next", "Next chunk of data"),
      actionButton("previous", "Previous chunk of data"),
      sliderInput("n_splits", "Number of runs in a split (how many runs are on a single plot):",
                  min = 50, max = 200, value = 100, step = 10)
    ),
    

    mainPanel(
      textOutput("chosen_split"),
      verbatimTextOutput("selected_info"),
      withSpinner(plotOutput("amp_plot", click = "plot_click"))
    )
  )
))
