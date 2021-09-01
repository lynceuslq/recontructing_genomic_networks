#!/bin/Rscript

library(shiny)
library(ggplot2)
library(patchwork)
library(easyGgplot2)
library(shinyjs)

load("comparingdiamondstats.RData")

ui <- fluidPage(
  titlePanel("Plotting conservation status of protein clusters"),
  
  sidebarLayout(
    sidebarPanel(
      
      selectInput("pars", 
                  h4("parameters tp compare"), 
                  choices =  list("identity" = 1,
                            "mismatches" = 2,
                            "scores" = 3),
                  selected = 1
      ),
      
      checkboxInput("showby", h4("grid by max number of hits"), value = T),
      

      checkboxGroupInput("setacc", 
                         h4("please select test datasets from:"), 
                         choices = unique(comstats$testset),
                         selected = unique(comstats$testset)[1:3]
      )
    ),
    
    
    mainPanel(
      dataTableOutput("pc_tab"),
      plotOutput("distPlot")
    )
  )
)

server <- function(input, output,session) {
  
  output$pc_tab <- renderDataTable({
    
    subset(comstats, comstats$testset == input$setacc)
    
  },options = list(pageLength = 5))
  
  output$distPlot <- renderPlot({
    p <- as.numeric(input$pars)
    print(p)
    plottab <- subset(comstats, comstats$testset == input$setacc)
    max_hits <- plottab$kvalue
    set <- plottab$testset
    
    if(p == 1) {
      avrg_id <- plottab$avrg_id_tp
      tp <- data.frame(max_hits = max_hits, 
                       set=set, 
                       avrg_id = avrg_id,
                       treatment = c("TP"))
      avrg_id <- plottab$avrg_id_fp
      fp <- data.frame(max_hits = max_hits, 
                       set=set, 
                       avrg_id = avrg_id,
                       treatment = c("FP"))
      plotme <- rbind(tp, fp)
      
      if(input$showby == T){
        pp <- ggplot(plotme, aes(x=max_hits, y=avrg_id, fill=treatment)) +
        geom_boxplot() +
        facet_wrap(~set) + 
        geom_boxplot(alpha = 0.75)
      }
      if(input$showby == F){
        pp <- ggplot(plotme, aes(x=set, y=avrg_id, fill=treatment)) +
          geom_boxplot() +
          facet_wrap(~max_hits) + 
          geom_boxplot(alpha = 0.75)
      }
    }
    
    if(p == 2) {
      avrg_mm <- plottab$avrg_mm_tp
      tp <- data.frame(max_hits = max_hits, 
                       set=set, 
                       avrg_mm = avrg_mm,
                       treatment = c("TP"))
      avrg_mm <- plottab$avrg_mm_fp
      fp <- data.frame(max_hits = max_hits, 
                       set=set, 
                       avrg_mm = avrg_mm,
                       treatment = c("FP"))
      plotme <- rbind(tp, fp)
      
      if(input$showby == T){
        pp <- ggplot(plotme, aes(x=max_hits, y=avrg_mm, fill=treatment)) +
          geom_boxplot() +
          facet_wrap(~set) + 
          geom_boxplot(alpha = 0.75)
      }
      if(input$showby == F){
        pp <- ggplot(plotme, aes(x=set, y=avrg_mm, fill=treatment)) +
          geom_boxplot() +
          facet_wrap(~max_hits) + 
          geom_boxplot(alpha = 0.75)
      }
    }
    
    if(p == 3) {
      avrg_score <- plottab$avrg_score_tp
      tp <- data.frame(max_hits = max_hits, 
                       set=set, 
                       avrg_score = avrg_score,
                       treatment = c("TP"))
      avrg_score <- plottab$avrg_score_fp
      fp <- data.frame(max_hits = max_hits, 
                       set=set, 
                       avrg_score = avrg_score,
                       treatment = c("FP"))
      plotme <- rbind(tp, fp)
      
      if(input$showby == T){
        pp <- ggplot(plotme, aes(x=max_hits, y=avrg_score, fill=treatment)) +
          geom_boxplot() +
          facet_wrap(~set) + 
          geom_boxplot(alpha = 0.75)
      }
      if(input$showby == F){
        pp <- ggplot(plotme, aes(x=set, y=avrg_score, fill=treatment)) +
          geom_boxplot() +
          facet_wrap(~max_hits) + 
          geom_boxplot(alpha = 0.75)
      }
      
    }
    
    pp
    
  }, height = 600, width = 1000)
  
  
  observeEvent(input$refresh, {
    session$invalidate
  })
  
}

shinyApp(ui = ui, server = server)

