library(shiny)

shinyUI(
  fluidPage(
    theme = "theme.css",
    tags$div(
      class = "site-backbone",
      tags$div(
        class = "logo-panel",
        img(src = "logo.png", class = "logo")
      )
    ),
    tags$style(HTML(".shiny-input-container:not(.shiny-input-container-inline) {
                        width: 100%;
                     }

                     pre{
                        background: white;
                     }

                     .shiny-notification {
                        height: 100px;
                        width: 800px;
                        position: fixed;
                        top: calc(50% - 50px);;
                        left: calc(50% - 400px);;
                    }
           ")),
    title = "PCRedux-app",
    
    headerPanel(""),
    
    # sidebarLayout(
    #   sidebarPanel(includeMarkdown("readme.md")
    #   ),
    
    fluidPage(uiOutput("dynamic_tabset"))
  )
)
