#!/bin/Rscript

library(shiny)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(ggpmisc)
library(easyGgplot2)
library(shinyjs)

#load("new.comparingdiamondstats.RData")

ui <- fluidPage(
  titlePanel("Plotting results from viral cluster profiling"),
  
  sidebarLayout(
    sidebarPanel(
      
      selectInput("showabd", 
                  h4("showing which plots below"), 
                  choices =  list("proportions of true abundance of test datasets" = 1,
                                  "diamond parameter selection" = 2,
                                  "profiling results and comparisons" = 3,
                                  "comparing abundance stats" = 4),
                  selected = 4
      ),
      
      selectInput("pars", 
                  h4("parameters tp compare"), 
                  choices =  list("identity" = 1,
                                  "mismatches" = 2,
                                  "scores" = 3),
                  selected = 1
      ),
      
      checkboxInput("showby", 
                    h4("grid by max number of hits"), 
                    value = T),
      
      textInput("min", 
                h4("min abundance to present:"), 
                value = "0.0001"),
      
      selectInput("est", 
                  h4("relative abundance esimated by"), 
                  choices =  list("negative log likelihood" = 1,
                                  "mean" = 2),
                  selected = 1
      ),
      checkboxGroupInput("tocompare", 
                  h4("comparing stats from different methods of estimation"), 
                  choices =  c("prop_fp_cl",
                                  "prop_fn_cl",
                                  "prop_num_tp_cl",
                                  "prop_num_fp_cl",
                                  "prop_num_fn_cl"),
                  selected = c("prop_fp_cl", "prop_num_fp_cl")
      ),
      sliderInput("yl", "up limit of y axis:", 
                  value = 0.05, min = 0, max = 1),
      checkboxGroupInput("selectk", 
                         h4("please select models from:"), 
                         choices = unique(prof$kvalue),
                         selected = unique(prof$kvalue)[1]
      ),
      checkboxGroupInput("setacc", 
                         h4("please select test datasets from:"), 
                         choices = unique(prof$testset),
                         selected = unique(prof$testset)[seq(1,30,5)]
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
    if(input$showabd == 3){
      mv <- as.numeric(input$min)
      if(input$est == 1) {
        tab <- subset(prof, prof$testset ==input$setacc && prof$prop_nll_abd_hm >= mv )
      }
      if(input$est == 2) {
        tab <- subset(prof, prof$testset ==input$setacc && prof$prop_sum_abd_hm_adj >= mv )
      }
    }
    if(input$showabd == 2) {
      tab <- subset(comstats, comstats$testset == input$setacc)
    }
    if(input$showabd == 1) {
      tab <- subset(true_abd, true_abd$testset == input$setacc)
    }
    if(input$showabd == 4){
      tab <- subset(sumall, sumall$testset == input$setacc)
      
    }
    tab
  },options = list(pageLength = 5))
  
  output$distPlot <- renderPlot({
    p <- as.numeric(input$pars)
    print(p)
    plottab <- subset(prof, prof$testset == input$setacc)
    max_hits <- plottab$kvalue
    set <- plottab$testset
    mv <- as.numeric(input$min)
    
    setlist <- unique(plottab$testset)
    
    if(input$showabd == 1){
      pp <- ggplot() +
        theme_void() +
        geom_text(aes(0,0,label='True relative abundance in test datasets')) +
        xlab(NULL)
      for(i in 1:length(setlist)){
        plotme <- subset(true_abd,true_abd$testset == setlist[i])
        bacprop <- sum(subset(plotme, plotme$bac != "NA" )$abundance)
        propclusteredphage <- sum(subset(plotme, plotme$cl != "NA" )$abundance)
        propunclusteredphage <- 1 - propclusteredphage - bacprop
        relative_abd <- c(propclusteredphage, propunclusteredphage, bacprop)
        overall_comp <- c("clustered_phages", "unclustered_phages", "bacteria")
        data <- data.frame(overall_comp = overall_comp, relative_abd = relative_abd)
        
        data$ymax <- cumsum(data$relative_abd)
        data$ymin <- c(0, head(data$ymax, n=-1))
        data$labelPosition <- (data$ymax + data$ymin) / 2
        data$label <- paste0(data$overall_comp, "\n value: ", round(data$relative_abd * 100, digits = 4), "%")
        
        mp <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=overall_comp)) +
          geom_rect() +
          geom_text( x=1.5, aes(y=labelPosition, label=label, color=overall_comp), size=3) + # x here controls label position (inner / outer)
          scale_fill_brewer(palette=3) +
          scale_color_brewer(palette=3) +
          coord_polar(theta="y") +
          xlim(c(-1, 4)) +
          theme_void() +
          theme(legend.position = "none") 
        
        
        cllist <- unique(subset(plotme, plotme$cl != "NA" )$cl)
        phabd <- c()
        for(k in 1:length(cllist)){
          phabd[k] <- sum(subset(plotme, plotme$cl == cllist[k])$abundance)
        }
        
        data <- data.frame(clustered_phage_comp = cllist, relative_abd = phabd)
        data<- data[order(data$relative_abd,decreasing = T),]
        #np <- ggplot(data=data, aes(x=clustered_phage_comp, y=relative_abd)) +
        #  geom_bar(stat="identity", width=0.5) + 
        #  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        np <- ggplot(data, aes(x=relative_abd)) +
          geom_histogram(color="black", fill="white")
        
        pp <- pp + ggarrange(mp,np, widths = 2, heights = 1, 
                             labels = c(paste("A",setlist[i],sep=":"), 
                                        paste("B",setlist[i],sep=":")), 
                             ncol = 2, nrow = 1)
        
      }
      
    }
    
    if(input$showabd == 2){
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
    }
    
    if(input$showabd == 3){
      pp <- ggplot() +
        theme_void() +
        geom_text(aes(0,0,label='comparing relative abundance calculated by the pipelines')) +
        xlab(NULL)
      
      for(i in 1:length(setlist)){
        prof1 <- subset(prof,prof$prop_nll_abd_hm >= mv)
        prof1 <- subset(prof, prof$kvalue == input$selectk)
        print(length(prof1$prop_nll_abd_hm))
        plotme <- subset(genomic, genomic$testset ==  setlist[i])
        if(length(plotme$cluster) > 0) {
          plotty <- data.frame(cluster=plotme$cluster, 
                               adj_ref_sum=plotme$adj_ref_sum, 
                               values=plotme$prop_sum_abd_gm, 
                               maxhits=c("genomic"), 
                               predictor=c("genomes"))
        }else{
          plotty <- data.frame()
        }
        plotme <- subset(prof1, prof1$testset == setlist[i])
        
        if(input$est ==1) {
          plotty2 <- data.frame(cluster=plotme$cluster, 
                                adj_ref_sum=plotme$adj_ref_sum, 
                                values=plotme$prop_nll_abd_hm, 
                                maxhits=plotme$kvalue, 
                                predictor=c("hallmarks"))
          
          
          plotty <- rbind(plotty, plotty2)
          
        }
        
        if(input$est == 2) {
          plotty2 <- data.frame(cluster=plotme$cluster, 
                                adj_ref_sum=plotme$adj_ref_sum, 
                                values=plotme$prop_sum_abd_hm_adj, 
                                maxhits=plotme$kvalue, 
                                predictor=c("hallmarks"))
          
          
          plotty <- rbind(plotty, plotty2)
        }
        
        
        formula <- y ~ x 
        pm <- ggplot(plotty, aes(x=adj_ref_sum,  y= values, color = predictor)) +
          geom_point(alpha = 0.3) +
          facet_wrap(~maxhits) +
          geom_smooth(method = "lm", formula = formula, se = F)  +
          stat_poly_eq(aes(label = paste(..rr.label..)), 
                       label.x.npc = "right", label.y.npc = 0.15,
                       formula = formula, parse = TRUE, size = 3) +
          stat_fit_glance(method = 'lm',
                          method.args = list(formula = formula),
                          geom = 'text',
                          aes(label = paste("P = ", signif(..p.value.., digits = 4), sep = "")),
                          label.x = 'right', label.y = 0.025, size = 3) +
          ylim(0, input$yl) +
          xlim(0, 0.05) +
          ggtitle(setlist[i])
        
        
        pp <- pp + pm
      }
    }
    
    if(input$showabd == 4){
      plotme <- data.frame()
    
      for(i in 1:length(setlist)) {
      plottmp <- subset(sumall, sumall$testset == setlist[i])
      print(setlist[i])
      plotme <- rbind(plotme, plottmp)
      }
      print(input$tocompare)
      chosefrom =  input$tocompare
      plotme1 <- data.frame()

      for(i in 1:length(chosefrom)) {
        print(chosefrom[i])
        plotme2 <- data.frame(
        kvalue = plotme$kvalue,
        testset = plotme$testset,
        stat_type = chosefrom[i],
        stat = plotme[,match(chosefrom[i], colnames(plotme))]
        
      )
        plotme1 <- rbind(plotme1, plotme2)
      }
      pp <- ggplot(plotme1, aes(x=kvalue, y=stat, fill=kvalue)) +
        geom_boxplot() +
        facet_wrap(~stat_type) +
        geom_boxplot(alpha = 0.75) +
        stat_compare_means(method = "anova")
        
    }
    # klist <- unique(plotme$kvalue)
    
    # for(k in 1:length(klist)) {
    #  print(paste(setlist[i], klist[k]))
    #   plotme <-subset(prof, prof$kvalue == klist[k])
    #  tt <- lm(prop_nll_abd_hm ~ adj_ref_sum, data=plotme)
    #   summary(tt)
    #   print(tab_model(tt))
    
    pp
    
  }, height = 1200, width = 1200)
  
  
  observeEvent(input$refresh, {
    session$invalidate
  })
  
}

shinyApp(ui = ui, server = server)

