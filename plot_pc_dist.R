library(shiny)
library(ggplot2)
library(easyGgplot2)
library(ggpubr)
library(RColorBrewer)
library(shinyjs)


ui <- fluidPage(
  titlePanel("Plotting conservation status of protein clusters"),
  
  sidebarLayout(
    sidebarPanel(
      
      textInput("text", h4("input genome cluster:"), 
                value = "cluster11"),
      
      selectInput("prots", h4("please select from:"), choices = unique(ws2$V1))
      
    ),
    
    
    mainPanel(
      plotOutput("distPlot"),
      dataTableOutput("pc_tab")
    )
  )
)

server <- function(input, output,session) {
  
  output$pc_tab <- renderDataTable({
    pc <- as.character(input$prots)
    pc <- na.omit(pc)
    
    subset(pc_dist_tab, pc_dist_tab$PC == pc)
    
  },options = list(pageLength = 10))
  
  output$distPlot <- renderPlot({
    
    t <- as.character(input$text)
    t <- na.omit(t)
    
    pc <- as.character(input$prots)
    pc <- na.omit(pc)

    plotme_tmp <- data.frame()
    plotme <- data.frame()
    
    conservation <- as.numeric(strsplit(subset(ws2, ws2$V1 == pc)$V3, ",")[[1]])
    window <- 1:length(conservation)
    plotme_tmp <- data.frame(window, conservation)
    plotme_tmp$ws <- 2
    plotme <- rbind(plotme, plotme_tmp)
    
    conservation <- as.numeric(strsplit(subset(ws4, ws4$V1 == pc)$V3, ",")[[1]])
    window <- 1:length(conservation)
    plotme_tmp <- data.frame(window, conservation)
    plotme_tmp$ws <- 4
    plotme <- rbind(plotme, plotme_tmp)
    ggplot(data = plotme, aes(x = window, y=conservation, group=ws)) + geom_line(aes(color=ws)) + scale_y_continuous()
    
  }, height = 400, width = 600)
  
  observeEvent(input$refresh, {
    session$invalidate
  })
  
  observeEvent(input$text, {
    t <- as.character(input$text)
    f <- subset(ws2, ws2$V2 == t)$V1
    updateSelectInput(session, inputId = "prots", label = paste("Select PCs from VC", t, sep=": "), choices = f)
  })

  
}

shinyApp(ui = ui, server = server)

