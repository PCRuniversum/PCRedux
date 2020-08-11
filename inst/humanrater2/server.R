library(shiny)
library(dplyr)
library(reshape2)
library(ggplot2)

raw_dat <- data.table::fread("../results/vermeulen_qPCR.csv", data.table = FALSE) 
cyc <- raw_dat[[1]]
dat <- raw_dat[, -1]





# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  all_splits <- reactive({
    split(1L:ncol(dat), ceiling(seq_along(1L:ncol(dat))/input[["n_splits"]]))
  })
  
  runs <- reactiveValues(
    selected = character(),
    split = 1
  )
  
  observeEvent(input[["next"]], {
    runs[["split"]] <- runs[["split"]] + 1
  })
  
  observeEvent(input[["previous"]], {
    runs[["split"]] <- runs[["split"]] - 1
  })
  
  observeEvent(input[["plot_click"]], {
      selected_runs <- as.character(nearPoints(plot_dat(), input[["plot_click"]], 
                                               xvar = "cycle", yvar = "value")[, "run"])
      
      if(length(selected_runs) > 0){
        remove <- selected_runs %in% runs[["selected"]]
        runs[["selected"]] <- runs[["selected"]][!runs[["selected"]] %in% selected_runs[remove]]
        runs[["selected"]] <- c(runs[["selected"]], selected_runs[!remove])
      }
  })
  
  plot_dat <- reactive({
    cbind(cycle = cyc, dat[all_splits()[[runs[["split"]]]]]) %>% 
      melt(id.vars = "cycle", variable.name = "run") %>% 
      mutate(selected = run %in% runs[["selected"]])
  })
  
  
  output[["amp_plot"]] <- renderPlot({
    dat <- plot_dat() 
    
    ggplot(dat, aes(x = cycle, y = value, color = run, group = run, linetype = selected)) +
      geom_line(size = 1.5) +
      geom_point(size = 3, show.legend = FALSE) +
      scale_linetype_manual(values = c("solid", "dashed")) +
      theme_bw(base_size = 17) +
      guides(color = "none") +
      facet_wrap(~ selected, ncol = 1)
    
  }, height = 800)
  

  output[["selected_info"]] <- renderPrint({
    dput(as.character(runs[["selected"]]))
  })
  
  output[["chosen_split"]] <- renderText({
    paste0("Chosen split: ", runs[["split"]], " (out of ", last(names(all_splits())), " splits).")
  })
  
})
