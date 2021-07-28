library(shiny)
library(ggplot2)
library(easyGgplot2)
library(ggpubr)
library(plot3D)
library(RColorBrewer)
library(shinyjs)

load('vcpc.new.RData')

ui <- fluidPage(
  titlePanel("Ploting proportions of shared genes against genome clusters"),
  
  sidebarLayout(
    sidebarPanel(
      
      textInput("text", h4("input genome cluster:"), 
                value = "cluster11"),
      sliderInput("min_prop_genome", "minimum proportion of genome to visualize:", 
                  value = 50, min = 0, max = 100),
      sliderInput("min_prop_proteins", "minimum proportion of proteins to visualize:", 
                  value = 50, min = 0, max = 100),
      sliderInput("dist", "threshold for jaccard similarity:", 
                  value = 0.5, min = 0, max = 1),
      numericInput("downj", "down_limit_of_Jaccard_index_(%)", 
                   value = 60, min = 0, max = 100),
      numericInput("upj", "up_limit_of_Jaccard_index_(%)", 
                   value = 80, min = 0, max = 100),
      
      numericInput("inter", "interval_of_Jaccard_index_(%)", 
                   value = 5, min = 1, max = 10),


      
    ),
    
    
    mainPanel(
      plotOutput("distPlot"),
      dataTableOutput("PC2VC")
    )
  )
)

server <- function(input, output,session) {
  
  output$PC2VC <- renderDataTable({
    
    bv <- subset(pairs, pairs$VC == input$text & pairs$prop_VC >= input$min_prop_genome/100 & pairs$prop_PC >= input$min_prop_proteins/100 & pairs$jaccard >= input$dist)
    bv[order(bv$jaccard, decreasing = TRUE),]
  },options = list(pageLength = 10))
  
  
  output$distPlot <- renderPlot({
    
    t <- as.character(input$text)
    t <- na.omit(t)
    
    w <- as.numeric(input$min_prop_genome) 
    
    u <- as.numeric(input$upj)
    d <- as.numeric(input$downj)
    v <- as.numeric(input$inter)
    
    p <- seq(100-w) / 100 + w/100
    
    m <- c()
    
    for(i in 1:length(p)){
      prop <- length(subset(pairs, pairs$VC == t & pairs$prop_VC >= p[i])$PC) / length(subset(pairs, pairs$VC == t)$PC)
      m <- c(m, prop)
    }
    
    plotme <- data.frame(prop_shared_genes=m, prop_genomes=p)
    
    bars <- ggplot(data=plotme, aes(x=prop_genomes, y=prop_shared_genes)) + geom_bar(stat="identity", fill="steelblue")
    
    
    plotme2 <- subset(pairs, pairs$VC == t)
    
    dots <- ggplot2.scatterplot(data=plotme2, x="prop_VC", y="prop_PC", color = "blue")  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    x <- sep.l <- plotme2$prop_VC
    y <- pet.l <- plotme2$prop_PC
    z <- sep.w <- plotme2$jaccard
    
    jdist <- as.ggplot(~scatter3D(x,y,z, bty = "g", theta = 20, phi = 0, xlab = "prop_VC",
                                  ylab ="prop_PC", zlab = "jaccard_index"
                                  ))
    
    countpcs <- data.frame()
    
    prop <- seq(d,u,v) / 100
    
    for(m in 1:length(prop)) {
      
      tmpp <- subset(pairs, pairs$gm_num_VC >= 10 & pairs$jaccard >= prop[m], select = c("PC", "VC"))
      
      countp <- data.frame(VC=unique(tmpp$VC))
      
      for(i in 1:length(countp$VC)){
        countp$sig_PC[i] <- length(subset(tmpp, tmpp$VC == countp$VC[i])$PC)
      }
      
      countp$jaccard <- paste(">", prop[m], sep="=")
      
      countpcs <- rbind(countpcs, countp)
      
    }
    
    box <- ggplot(countpcs, aes(x=jaccard, y=sig_PC, fill = jaccard )) + 
      geom_boxplot(alpha = 0.75) + 
      geom_jitter() + 
      scale_fill_brewer(palette = "Set3") +
      scale_y_continuous(name = "Number of PC per VC pass the threshold", 
                         limits=c(0, 200)) +
      scale_x_discrete(name = "Threshold for Jaccard index")
    
    
    ggarrange(bars, dots, jdist, box, widths = 2, heights = 4, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)
  }, height = 400, width = 800)
  
  observeEvent(input$refresh, {
    session$invalidate
  })
  
}

shinyApp(ui = ui, server = server)

