library(shiny)
library(ggplot2)
library(patchwork)
library(easyGgplot2)
library(shinyjs)

load("plot_pc_dist.RData")

ui <- fluidPage(
  titlePanel("Plotting conservation status of protein clusters"),
  
  sidebarLayout(
    sidebarPanel(
      
      textInput("text", h4("input genome cluster:"), 
                value = "cluster11"),
      
      checkboxInput("bypc", h4("table showing only selected PC:"), value = F),
      
      checkboxInput("byth", h4("graph showing only PC which pass the threshoulds below:"), value = F),
      
      sliderInput("numcl", "maximum number of clustered proteins:", 
                  value = 5, min = 1, max = 10),
      
      sliderInput("maxd", "maximum pairwise distance of clustered proteins:", 
                  value = 0.5, min = 0, max = 1),
      
      sliderInput("mead", "mean pairwise distance of clustered proteins:", 
                  value = 0.1, min = 0, max = 1),
      
      sliderInput("cons", "threshould for conservation status (window size 4):", 
                  value = 2, min = 0, max = 4),
      
      sliderInput("perws", "percetage of windows beyond the threshould:", 
                  value = 80, min = 0, max = 100),
      
      #selectInput("prots", h4("please select from:"), choices = unique(ws2$V1))
      
      checkboxGroupInput("prots", h4("please select from:"), choices = unique(ws2$V1))
    ),
    
    
    mainPanel(
      dataTableOutput("pc_tab"),
      plotOutput("distPlot")
    )
  )
)

server <- function(input, output,session) {

  output$pc_tab <- renderDataTable({
    print(input$prots)
    
    if(input$bypc == FALSE) {
      t <- as.character(input$text)
      t <- na.omit(t)
      
      perc <- input$perws / 100
      
      showtmp <- subset(pc_dist_tab, pc_dist_tab$VC == t & 
                          pc_dist_tab$clustered_num <= input$numcl & 
                          pc_dist_tab$mean_dist <= input$mead & 
                          pc_dist_tab$max_dist <= input$maxd)
      
      
      pp <- showtmp$PC
      print(pp)
      
      con <- c()
      for(i in 1:length(pp)) {
        
        pc <- pp[i]
        
        truew <- sum(as.numeric(strsplit(subset(ws2, ws2$V1 == pc)$V3, ",")[[1]]) >= input$cons)
        con[i] <- truew / length(strsplit(subset(ws2, ws2$V1 == pc)$V3, ",")[[1]])
      }
      ppx <- data.frame(pc=pp, con=con)
      pc <- subset(ppx, ppx$con >= perc)$pc
      print(pc)
      show <- showtmp[match(pc, showtmp$PC),]
    }
    
    
    if(input$bypc == TRUE) {
      show <- pc_dist_tab[match(input$prots, pc_dist_tab$PC),]
    }
    
    show
    
  },options = list(pageLength = 5))
  
  output$distPlot <- renderPlot({
    
    t <- as.character(input$text)
    t <- na.omit(t)
    
    if(input$byth == FALSE) {
      
      pp <- as.list(input$prots)
      
    }
    
    if(input$byth == TRUE) {
      t <- as.character(input$text)
      t <- na.omit(t)
      
      perc <- input$perws / 100
      
      showtmp <- subset(pc_dist_tab, pc_dist_tab$VC == t & 
                          pc_dist_tab$clustered_num <= input$numcl & 
                          pc_dist_tab$mean_dist <= input$mead & 
                          pc_dist_tab$max_dist <= input$maxd)
      
      
      pp <- showtmp$PC
      print(pp)
      
      con <- c()
      for(i in 1:length(pp)) {
        
        pc <- pp[i]
        
        truew <- sum(as.numeric(strsplit(subset(ws2, ws2$V1 == pc)$V3, ",")[[1]]) >= input$cons)
        con[i] <- truew / length(strsplit(subset(ws4, ws4$V1 == pc)$V3, ",")[[1]])
      }
      ppx <- data.frame(pc=pp, con=con)
      pp <- as.list(subset(ppx, ppx$con >= perc)$pc)
      
    }
    
    
    dat <- pc_dist_tab[match(pp,pc_dist_tab$PC),]
    
    plotf <- ggplot(data = dat, aes(x=PC, y=jaccard)) + 
      geom_bar(stat="identity", fill="steelblue")
  
    
    for(i in 1:length(pp)) {
      
      pc <- pp[i]
      
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
      
      plot <- ggplot(data = plotme, aes(x = window, y=conservation, group=ws)) + 
        geom_line(aes(color=ws)) + scale_y_continuous() + ggtitle(pc)
      
      plotf <- plotf + plot
    }
    
    plotf + plot_layout(ncol = 3)
    
  }, height = 800, width = 1000)
  
  
  
  observeEvent(input$refresh, {
    session$invalidate
  })
  
  observeEvent(input$text, {
    t <- as.character(input$text)
    f <- subset(ws2, ws2$V2 == t)$V1
    updateCheckboxGroupInput(session, inputId = "prots", label = paste("Select PCs from VC", t, sep=": "), choices = f, selected = f[1])
  })
  
  
}

shinyApp(ui = ui, server = server)

