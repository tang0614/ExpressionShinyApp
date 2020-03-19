library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(shiny)
library(shinyWidgets)
library(shinythemes)
library(shinycssloaders)
library(igraph)
library(ggraph)
library(scales)
library(DT)
library(treemap)
library(highcharter)
library(purrr)
library(stringr)
library(fuzzyjoin)
library(lexicon)
library(visNetwork)
library(httr)
library(rdrop2)
library(lubridate)
library(ggvis)
library(rsconnect)
library(plotly)
library(d3heatmap)
library(reshape2)
if (FALSE) {
  library(RSQLite)
  library(dbplyr)
}
source("helpers.R")
#CSS
button_color_css <- "
#DivCompClear, #FinderClear, #EnterTimes{
/* Change the background color of the update button
to blue. */
background: DodgerBlue;

/* Change the text size to 15 pixels. */
font-size: 15px;
}"

#Margin for all plots
mar = list(b=50)

filter_data <- function(tissue_input, gene_input, m, gene_col,tissue_col) {
  
  list_1<-m[ , c(gene_col)]
  sub_list1<- list_1 %in% gene_input
  
  list_2<-m[ , c(tissue_col)]
  sub_list2<- list_2 %in% tissue_input
  
  if ((is.null(gene_input))&&(is.null(tissue_input)) ){
    req(tissue_input)
    req(gene_input)
  } else
    if ((is.null(gene_input)) && (!is.null(tissue_input))) {
      m <- subset(m, sub_list2)
    } else
      
      if ((!is.null(gene_input)) && (is.null(tissue_input))){
        m <-subset(m, sub_list1)
      }else
        
        if ((!is.null(gene_input)) && (!is.null(tissue_input))){
          n <-subset(m, sub_list1)
          list_3<-n[ , c(tissue_col)]
          sub_list3<- list_3 %in% tissue_input
          
          m <- subset(n, sub_list3)
        }
  
  return(m)
  
}      



hearmap_matrix<-function(datainput, a, b, y_axis_name){
  m<-datainput
  
  if(TRUE){
    n<- m%>%
      group_by_(a,b)%>% 
      summarise(Median=median(log_yaxis_value))
  }else{
    n<- m%>%
      group_by_(a,b)%>% 
      summarise(Median=median(yaxis_value))
  }
  
  n <- n[is.finite(n$Median), ]
  df <-as.data.frame.matrix(xtabs(Median~., n))
  matrix <- data.matrix(df)
  matrix 
  
}

draw_heatmap<- function(heatmapInput, heatmapInput_log_button,input_gene){
  a <- list(
    showticklabels = TRUE,
    tickangle = -45
  )
  
  
  matrix<-heatmapInput
  x<-c(colnames(matrix))
  y<-c(rownames(matrix))
  
  
  p <-plot_ly(z = matrix, type = "heatmap", 
              colorscale= "Hot",
              x =colnames(matrix),
              y =rownames(matrix),
              reversescale=TRUE)
  
  
  if(heatmapInput_log_button){
    p <- layout(p,
                title = paste0('Log2 Median TPM Heatmap'), 
                margin=mar,xaxis = a)
    
    if(input_gene>30){
      
      p <- layout(p,height = 1000)
      
    }else if(input_gene>10){
      p <- layout(p,height = 800)
      
    }else{
      p <- layout(p,height = 500)
      
    }
  }else{
    p <- layout(p,
                title = paste0("Median TPM Heatmap"), 
                margin=mar,xaxis = a)
    
    if(input_gene>30){
      
      p <- layout(p,height = 1000)
      
    }else if(input_gene>10){
      p <- layout(p,height = 900)
      
    }else{
      p <- layout(p,height = 500)
      
    }
    
  }
  
  p
  
}


sortedData_1 <- function(dataInput){
  
  df <-dataInput  %>% group_by(cc) %>% summarize(n=n(),log_m=median(log_yaxis_value),m=median(yaxis_value))
  
  df['cc'] <- data.frame(lapply(df['cc'], as.character), stringsAsFactors=FALSE)
  dataInput['cc'] <- data.frame(lapply(dataInput['cc'], as.character), stringsAsFactors=FALSE)
  
  df2<-merge(x = dataInput, y = df, by = "cc", all.x = TRUE)
  df2
  
}


#groupby sorting
sortedData_2 <- function(dataInput, gene_col,tissue_col,control_gene,control_tissue){
  
  
  if((control_gene)&&(control_tissue)){
    lvls <- dataInput %>%
      group_by(cc) %>%
      summarise(m = median(yaxis_value))%>%
      arrange(m) %>%
      pull(cc)
  } else if(!(control_gene)&!(control_tissue)){
    lvls <- dataInput %>%
      group_by(cc) %>%
      summarise(m = median(yaxis_value))%>%
      arrange(m)%>%
      pull(cc)
  }else if(control_tissue){
    lvls <- dataInput %>%
      group_by_(tissue_col,'cc')%>%
      summarise(m = median(yaxis_value))%>%
      arrange(!!rlang::sym(tissue_col),m)%>%
      pull(cc)
    
  } else if(control_gene){
    lvls <- dataInput %>%
      group_by_(gene_col,'cc')%>%
      summarise(m = median(yaxis_value))%>%
      arrange(!!rlang::sym(gene_col),m)%>%
      pull(cc)
  }
  return(lvls)
  
}



draw_scatter <- function(datainput,logarithmic_button,rank_button,sortedData_2,sortedData_1){
  
  ax = list(showticklabels = TRUE,tickangle = -45)
  
  if(length(unique(datainput$cc)) <= 45){
    
    
    
    if(logarithmic_button){
      
      if(rank_button){
        
        p <- sortedData_1 %>%
          plot_ly(x = ~factor(cc,sortedData_2), 
                  y = ~log_yaxis_value,
                  height = 500,
                  type = 'violin',
                  scalemode = 'count',
                  
                  points='all', 
                  jitter=0.01,  
                  alpha = 0.1,
                  pointpos=0,
                  
                  
                  box = list(
                    visible = T
                  ),
                  meanline = list(
                    visible = T
                  ),
                  color = ~GeneName, 
                  colors ='Set2',
                  hoverinfo = "text",
                  hovertext = paste("Gene :", sortedData_1[,'GeneName'],
                                    "<br>Tissue :", sortedData_1[,'Tissue'],
                                    "<br>Median :", sortedData_1[,'log_m'],
                                    "<br>TPM :", sortedData_1[,'log_yaxis_value'],
                                    "<br>Sample size :", sortedData_1[,'n'])
          ) 
        
        
      }
      
      else{
        
        p <- sortedData_1 %>%
          plot_ly(x = ~cc, 
                  y = ~log_yaxis_value,
                  height = 500,
                  type = 'violin',
                  scalemode = 'count',
                  
                  points='all', 
                  jitter=0.01,  
                  alpha = 0.1,
                  pointpos=0,
                  
                  
                  box = list(
                    visible = T
                  ),
                  meanline = list(
                    visible = T
                  ),
                  color = ~GeneName, 
                  colors ='Set2',
                  hoverinfo = "text",
                  hovertext = paste("Gene :", sortedData_1[,'GeneName'],
                                    "<br>Tissue :", sortedData_1[,'Tissue'],
                                    "<br>Median :", sortedData_1[,'log_m'],
                                    "<br>TPM :", sortedData_1[,'log_yaxis_value'],
                                    "<br>Sample size :", sortedData_1[,'n']))
        
      }
      
      
      
      p<-layout(p,
                title = paste0("Log2 GTEx Violin Plot"), 
                xaxis = ax,
                yaxis = list(title =paste0("Log2 TPM")),
                margin=mar)
      
      
      
      
      
    }else{
      if(rank_button){
        
        p <- sortedData_1 %>%
          plot_ly(x = ~factor(cc,sortedData_2), 
                  y = ~yaxis_value,
                  height = 500,
                  type = 'violin',
                  scalemode = 'count',
                  
                  points='all', 
                  jitter=0.01,  
                  alpha = 0.1,
                  pointpos=0,
                  
                  
                  box = list(
                    visible = T
                  ),
                  meanline = list(
                    visible = T
                  ),
                  color = ~GeneName, 
                  colors ='Set2',
                  hoverinfo = "text",
                  hovertext = paste("Gene :", sortedData_1[,'GeneName'],
                                    "<br>Tissue :", sortedData_1[,'Tissue'],
                                    "<br>Median :", sortedData_1[,'m'],
                                    "<br>TPM :", sortedData_1[,'yaxis_value'],
                                    "<br>Sample size :", sortedData_1[,'n']))
        
      }
      
      else{
        p <- sortedData_1 %>%
          plot_ly(x = ~cc, 
                  y = ~yaxis_value,
                  height = 500,
                  type = 'violin',
                  scalemode = 'count',
                  
                  points='all', 
                  jitter=0.01,  
                  alpha = 0.1,
                  pointpos=0,
                  
                  
                  box = list(
                    visible = T
                  ),
                  meanline = list(
                    visible = T
                  ),
                  color = ~GeneName, 
                  colors ='Set2',
                  hoverinfo = "text",
                  hovertext = paste("Gene :", sortedData_1[,'GeneName'],
                                    "<br>Tissue :", sortedData_1[,'Tissue'],
                                    "<br>Median :", sortedData_1[,'m'],
                                    "<br>TPM :", sortedData_1[,'yaxis_value'],
                                    "<br>Sample size :", sortedData_1[,'n']))
        
        
        
      }
      
      p<-layout(p,
                title = paste0('GTEx Violin Plot'), 
                xaxis = ax,
                yaxis = list(title ='TPM'),
                margin=mar)
      
      
    }
    
  }
  
  
}