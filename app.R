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
library(gtools)
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
source("help_function.R")


#Upload files

#Upload  data
rna_tissue= read.csv('autophagy_gtex_outerjoin.csv', header = TRUE,stringsAsFactors=FALSE)
rna_tissue_hpa= read.csv('autophagy_hpa_outerjoin.csv', header = TRUE,stringsAsFactors=FALSE)

ccle = read.csv('ccle_autophogy_tpm.csv', header = TRUE,stringsAsFactors=FALSE)


#Create variables for dropdown input 
yaxis_vars <- c('TPM'='GETx_TPM')
yaxis_vars_hpa <- c('TPM'='hpa_TPM')
yaxis_vars_ccle <- c('TPM'='TPM')


#CSS
button_color_css <- "
#DivCompClear, #FinderClear, #EnterTimes{
/* Change the background color of the update button
to blue. */
background: DodgerBlue;

/* Change the text size to 15 pixels. */
font-size: 15px;
}"




## Shiny uses the function fluidPage to create a display that automatically adjusts to the dimensions of your user’s browser window. 
# Define UI ----

ui <- fluidPage(
      #Setup navigation bar
      navbarPage("Autophagy Gene Expression", theme = shinytheme("lumen"),
      
              #First Page
              tabPanel("Browse Tissue and Gene", fluid = TRUE, icon = icon("globe-americas"),tags$style(button_color_css),
                    
                      
                    # Sidebar layout with a input and output definitions
                    sidebarLayout(
                      
                    
                            #Input Panel
                            sidebarPanel(
                                  #Input Panel Title
                                  titlePanel("Visualize Gene Expression by Gene and Tissue"),
                              
                                  hr(),
                                  helpText("Choose from GTEx Dataset and HPA Dataset"),
                                  
                                  
                                  fluidRow(
                                          
                                    
                                          column(8, 
                                                
                                                   selectInput("gene", label = "Gene Symbol", 
                                                               choices=mixedsort(as.vector(unique(rna_tissue$GeneName))), selected=c('ATG10','ULK1'),multiple=TRUE
                                                                 )
                                                  
                                          ),
                                          br(),
                                          br(),
                                          column(4,prettyCheckbox("gene_selectall", "Select All"))
                                  ),
                                  
                                  
                                  br(),
                                  fluidRow(
                                          column(8,
                                                 selectInput("tissue", label = "Tissue Name", 
                                                                choices=mixedsort(as.vector(unique(rna_tissue$Tissue))), selected=c('liver','lung'),multiple=TRUE
                                                               )
                                               
                                                 
                                         ),
                                         
                                         br(),
                                         br(),
                                         column(4,prettyCheckbox("tissue_selectall", "Select All"))
                                  ),
                                  
                                  fluidRow(column(11,
                                        helpText("Tip: You can select all genes/tissues by clicking \"Select ALL\". 
                                         You can also select all tissues by only entering gene symbols; and
                                         select all genes by only entering tissue symbols.")
                                  )),
                             
                                 
                                  br(),
                                  hr(),
                                  fluidRow(
                                    
                                    column(10,
                                           span(textOutput("alertOne"),style="color:red")
                                    )
                                  ),
                                  
                                  br(),
                                  fluidRow(
                                    
                                    column(10,
                                           span(textOutput("alertTwo"),style="color:red")
                                    )
                                  ),
                                  br(),
                                  hr(),
                                  tags$p(span("Large graphs (e.g., selecting lots of genes and tissues) may take a few seconds to render.", style = "color:red")),
                                  tags$p(HTML("<b>Transcripts Per Million (TPM)</b> is a normalization method for RNA-seq, should be read as: for every 1,000,000 RNA molecules in the RNA-seq sample, x came from this gene.")),
                                  #tags$p(HTML("<b>pTPM</b> is calculated from scaling a sum of 1 million TPM to compensate for the non-coding transcripts that had been previously removed.")),
                                  tags$p(HTML("<b>GTEx Description</b> : Dataset comes from the 8th release of <a href=\"https://gtexportal.org/home/documentationPage\">Genotype-Tissue Expression</a>
                                  , which is a
                                  project that build a comprehensive public resource to study tissue-specific 
                                  gene expression and regulation. Samples were collected from 54 non-diseased tissue sites across 
                                  nearly 1000 individuals, primarily for molecular assays including WGS, WES, and RNA-Seq.")),
                                  
                                  tags$p(HTML("<b>HPA Description</b> : Dataset comes from the 19th release of the <a href=\"https://www.proteinatlas.org/about\">Human Protein Atlas</a>, which is a Swedish-based program initiated in 2003 with the aim to map 
                                  all the human proteins in cells, tissues and organs using integration of various omics technologies, 
                                  including antibody-based imaging, mass spectrometry-based proteomics, transcriptomics and systems biology. 
                                  All the data in the knowledge resource is open access to allow scientists both in 
                                              academia and industry to freely access the data for exploration of the human proteome. "))
                            ),
                                      
                    
                    #Output     
                    mainPanel(
                           
                            
                            fluidRow(
                                    column(11,offset = 1,
                                            helpText("Notes: NaN values (missing values) contained in a matrix are represented as white colors.")
                                    )
                            ),
                            br(),
                          
                            fluidRow(
                           
                                    column(2, 
                                           offset = 1,
                                           prettyCheckbox("logarithmicY_heat", "Log2 of Y", FALSE,
                                                          shape = "round", 
                                                          bigger=TRUE,
                                                          outline =TRUE, 
                                                          animation = "smooth")
                                    ),
                                    
                                    
                            ),
                            
                            br(),
                       
                            fluidRow(column(4,offset=4,h3('Median TPM Heatmap'))),
                            uiOutput("heatmap_ui"),
                            uiOutput("heatmap_hpa_ui"),
                     
                           
                            br(),
                            hr(),
                            
                            column(12,
                                   
                                   fluidRow(
                                     column(10,offset = 1,
                                            helpText("Notes: Sum of selected gene and tissue should be lower than 45 for voilin plot.")
                                     )
                                   ),
                                    
                                   fluidRow(
                                      column(5, offset = 1,
                                                   prettyCheckbox("rank", "Median sort all samples (click either one of below to sort by tissues or genes):", FALSE,
                                                                 shape = "round", 
                                                                 bigger=TRUE,
                                                                 outline =TRUE, 
                                                                 animation = "smooth")
                                            ),
                                      column(2,offset = 7,
                                             prettyCheckbox("control_tissue", "(1) By tissue", FALSE,
                                                            shape = "curve", 
                                                            bigger=FALSE,
                                                            animation = "smooth")
                                      ),
                                      column(2,
                                             prettyCheckbox("control_gene", "(2) By gene", FALSE, 
                                                            shape = "curve", 
                                                            bigger=FALSE,
                                                            animation = "smooth")
                                      )
            
                                  ),
                                  
                                 
                                  
                                  
                                
                                fluidRow(
                                
                                  column(3, offset =1 ,
                                         prettyCheckbox("logarithmicY_box", "Log2 of Y", TRUE,
                                                        shape = "round", 
                                                        bigger=TRUE,
                                                        outline =TRUE, 
                                                        animation = "smooth")
                                         
                                  )
                                      
                                
                                ),
                                
                                
                                fluidRow(column(4,offset=4,h3('TPM Voilin Plot')))
                          ),
                           
                          
                          br(),
                          br(),
                          br(),
                  
                          fluidRow(
                            
                                  column(11, withSpinner(plotlyOutput("scatter"))
                                         )
                          ),
                          
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
            
                          
                          fluidRow(
                            
                                  column(11, withSpinner(plotlyOutput("scatter_hpa")))
                          )
                           
        
                      )
              )
      ),
              
      tabPanel("Browse Cell Line and Gene", fluid = TRUE, icon = icon("globe-americas"),
               tags$style(button_color_css),
               
                         
                # Sidebar layout with a input and output definitions
                sidebarLayout(
                  
                  
                        #Input Panel
                        sidebarPanel(
                          
                                #Input Panel Title
                                titlePanel("Visualize Gene Expression by Gene and Cell-Line"),
                               
                               
                                hr(),
                                helpText("Choose from CCLE Dataset"),
                                
                                
                                fluidRow(
                                  
                                  
                                  column(8, 
                                         
                                         selectInput("gene_ccle", label = "Gene Symbol", 
                                                     choices=mixedsort(as.vector(sort(unique(ccle$GeneName)))), selected=c('ATG10','ULK1'),multiple=TRUE
                                         )
                                         
                                  ),
                                  br(),
                                  br(),
                                  column(4,prettyCheckbox("gene_selectall_cellline", "Select All"))
                                ),
                                
                                
                                br(),
                                fluidRow(
                                  column(8,
                                         selectInput("cell_line", label = "Cell line", 
                                                     choices=mixedsort(as.vector(sort(unique(ccle$Name)))), selected=c('22Rvl1','59M'),multiple=TRUE
                                         )
                                         
                                         
                                  ),
                                  
                                  br(),
                                  br(),
                                  column(4,prettyCheckbox("cellline_selectall", "Select All"))
                                ),
                                fluidRow(column(11,
                                                helpText("Tip: You can select all genes/cell-lines by clicking \"Select ALL\". 
                                         You can also select all tissues by only entering gene symbols; and
                                         select all genes by only entering tissue symbols.")
                                )),
                                
                                
                                
                     
                                         
                                br(),
                                hr(),
                                fluidRow(
                                  
                                  column(10,
                                         span(textOutput("alertThree"),style="color:red")
                                  )
                                ),
                                
                                tags$p(HTML("<b>CCLE Description</b> : Dataset comes from the latest release of <a href=\"https://portals.broadinstitute.org/ccle/about\">Cancer Cell Line Encyclopedia</a>
                                  , which is a
                                  project of collaboration between the Broad Institute, 
                                            and the Novartis Institutes for Biomedical Research
                                            and its Genomics Institute of the Novartis Research Foundation. It 
                                            conducts a detailed genetic and pharmacologic characterization of a large panel of human cancer models, 
                                            to develop integrated computational analyses that link distinct pharmacologic vulnerabilities to genomic patterns 
                                            and to translate cell line integrative genomics into cancer patient stratification. The CCLE provides public access to genomic data, analysis and visualization for over 1100 cell lines."))
                                
                                
                        ),
                        
                        #Output Panel
                        mainPanel(
                              fluidRow(column(11,
                                              helpText("Tip: You can select all genes/tissues by clicking \"Select ALL\". 
                                             You can also select all tissues by only entering gene symbols; and
                                             select all genes by only entering tissue symbols.")
                              )),
                                  
                              fluidRow(
                                
                                column(2, 
                                       offset = 1,
                                       prettyCheckbox("logarithmicY_heat_ccle", "Log2 of Y", FALSE,
                                                      shape = "round", 
                                                      bigger=TRUE,
                                                      outline =TRUE, 
                                                      animation = "smooth")
                                ),
                                column(4, 
                                       offset = 4,
                                       prettyCheckbox("fix_heatmap_ccle", "Fix current heatmap", FALSE,
                                                      shape = "round", 
                                                      bigger=TRUE,
                                                      outline =TRUE, 
                                                      animation = "smooth")
                                )
                              ),
                        
                          br(),
                          
                          h3('TPM Heatmap'),
                          uiOutput("heatmap_ccle_ui"),
                          
                          
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          hr(),
                          
                          column(12,
                                 
                                 fluidRow(helpText("Tip: (1) click submit every time after change any input. (2) zoom in to see individual points. (3) sum of gene and tissue should be lower than 45 for voilin plot.")),
                                 br(),
                                 br(),
                                 
                                 fluidRow(
                                   
                                   column(3, offset = 1,
                                          prettyCheckbox("logarithmicY_box_ccle", "Log2 of Y", FALSE,
                                                         shape = "round", 
                                                         bigger=TRUE,
                                                         outline =TRUE, 
                                                         animation = "smooth")
                                          
                                   ),
                                   column(3, 
                                          prettyCheckbox("tick_ccle", "Show Ticks", FALSE,
                                                         shape = "round", 
                                                         bigger=TRUE,
                                                         outline =TRUE, 
                                                         animation = "smooth")
                                   ),
                                   column(4, 
                                          prettyCheckbox("rank_ccle", "Median Sort (check one of below to group before median sort):", FALSE,
                                                         shape = "round", 
                                                         bigger=TRUE,
                                                         outline =TRUE, 
                                                         animation = "smooth")
                                   )
                                   
                                 ),
                                 
                                 
                                 
                                 
                                 fluidRow(
                                   
                                   column(3, offset = 7,
                                          prettyCheckbox("control_tissue_ccle", "(1) Cell-line sorting", FALSE,
                                                         shape = "curve", 
                                                         bigger=FALSE,
                                                         animation = "smooth")
                                   )
                                   
                                   
                                 ),
                                 
                                 fluidRow(
                                   
                                   column(3, offset = 7,
                                          prettyCheckbox("control_gene_ccle", "(2) Gene sorting", FALSE, 
                                                         shape = "curve", 
                                                         bigger=FALSE,
                                                         animation = "smooth")
                                   )
                                 )
                          ),
                          
                          br(),
                          br(),
                          br(),
                          
                          fluidRow(
                            
                            column(11, withSpinner(plotlyOutput("scatter_ccle"))
                            )
                          )
                          
                    

                          
                        )
                )
                            
               

        
      )
                
      
      
))

      


# Define server logic ----how to assenble input into output

server <- function(input, output,session) {
  

  #select button
  observe({
    
          if (input$gene_selectall){
                  updateSelectInput(session, "gene", selected = as.vector(unique(rna_tissue$GeneName)))
          }else{
                  updateSelectInput(session, "gene", selected = c('ATG10','ULK1'))
          }
    
          if (input$tissue_selectall){
            updateSelectInput(session, "tissue", selected = as.vector(unique(rna_tissue$Tissue)))
          }else{
            updateSelectInput(session, "tissue", selected = c('liver','lung'))
          }
    
          if (input$gene_selectall_cellline){
            updateSelectInput(session, "gene_ccle", selected = as.vector(unique(ccle$GeneName)))
          }else{
            updateSelectInput(session, "gene_ccle", selected = c('ATG10','ULK1'))
          }
          
          if (input$cellline_selectall){
            updateSelectInput(session, "cell_line", selected = as.vector(unique(ccle$Name)))
          }else{
            updateSelectInput(session, "cell_line", selected = c('59M','A101D'))
          }
    
    
          
  })
 
  
  #GTEx data input
  dataInput <- reactive ({
    #dataInput <- eventReactive (input$heat_button,{
    m<-rna_tissue
    m$yaxis_value <-m$gtex_TPM
    m$log_yaxis_value <-log2(m$yaxis_value)
    m<-filter_data(input$tissue,input$gene,m,'GeneName','Tissue')
    m$cc <- interaction(m$Tissue, m$GeneName)
    m

    
    
  })
  
  #Datainput hpa
  dataInput_hpa <- reactive({
    #dataInput_hpa <- eventReactive(input$heat_button, {
    m<-rna_tissue_hpa
    m$yaxis_value <-m$hpa_TPM
    m$log_yaxis_value <-log2(m$yaxis_value)
    m<-filter_data(input$tissue_hpa,input$gene_hpa,m,'GeneName','Tissue')
    m$cc <- interaction(m$Tissue, m$GeneName)
    m
  })
  
  
  #ccle data input
  dataInput_ccle <- reactive({
    m<-ccle
    m$yaxis_value <-m$TPM
    m$log_yaxis_value <-log2(m$yaxis_value) 
    m<-filter_data(input$cell_line,input$gene_ccle,m,'GeneName','Name')
    m$cc <- interaction(m$Name, m$GeneName)
    m
  })
  
  
  

    
  heatmapInput <-reactive({
    
    n<- dataInput()%>%
      group_by(Tissue,GeneName)
      
    if(input$logarithmicY_heat){n<- n%>%summarize(m=median(log_yaxis_value))}
    else if (!input$logarithmicY_heat){n<- n%>% summarise(m=median(yaxis_value))}
    
    n[is.na(n)] <- 0
    n <- n[is.finite(n$m), ]
    df <-as.data.frame.matrix(xtabs(m~., n))
    data.matrix(df)
  
    
  })
   
  heatmapInput_hpa <-reactive({

    n<- dataInput_hpa()%>%
      group_by(Tissue,GeneName)
   
    if(input$logarithmicY_heat){n<- n%>%summarize(m=median(log_yaxis_value))}
    else if (!input$logarithmicY_heat){n<- n%>% summarise(m=median(yaxis_value))}
    
    n[is.na(n)] <- 0
    n <- n[is.finite(n$m), ]
    df <-as.data.frame.matrix(xtabs(m~., n))
    data.matrix(df)
    
   
  })
  
  heatmapInput_ccle <-reactive({
   
    n<-dataInput_ccle()%>%
      group_by(GeneName,Name)
    
    
    if(input$logarithmicY_heat_ccle){n<- n%>%summarize(m=median(log_yaxis_value))}
    else if (!input$logarithmicY_heat_ccle){n<- n%>% summarise(m=median(yaxis_value))}
    
    n[is.na(n)] <- 0
    n <- n[is.finite(n$m), ]
    df <-as.data.frame.matrix(xtabs(m~., n))
    data.matrix(df)
    
 
  })
  
  #heatmap for tissue scale setup
  boundary <- reactive({
    matrix_gtex<-heatmapInput()
    matrix_hpa<-heatmapInput_hpa()
    
    matrix_gtex_max<-max(matrix_gtex)
    matrix_hpa_max<-max(matrix_hpa)
    
    if(matrix_gtex_max>matrix_hpa_max){
      return(matrix_gtex_max)
    }else{
      return(matrix_hpa_max)
    }
    
  })
  
  #heatmap for celline scale setup
  boundary_2 <- reactive({
    matrix_gtex<-heatmapInput_ccle()
    matrix_hpa<-heatmapInput_ccle()
    
    matrix_gtex_max<-max(matrix_gtex)
    matrix_hpa_max<-max(matrix_hpa)
 
    
    if(matrix_gtex_max>matrix_hpa_max){
      return(matrix_gtex_max)
    }else{
      return(matrix_hpa_max)
    }
    
  })
  
  # draw heatmap
  output$heatmap_gtex  <- renderPlotly({
    
    
    a <- list(
      showticklabels = TRUE,
      tickangle = -45
    )
    
    
    matrix<-heatmapInput()
    p <-plot_ly(z = matrix, type = "heatmap", 
                colorscale= "Hot",
                x =colnames(matrix),
                y =rownames(matrix),
                reversescale=TRUE)
    
    
    if(input$logarithmicY_heat){
      
      p <- layout(p,
                  title = paste0('GTEX'), 
                  margin=mar,xaxis = a)%>%
        colorbar(limits = c(0, boundary()))
      
      if(nrow(matrix)>20){
        
        p <- layout(p,height = 800)
        
      }else if(nrow(matrix)>10){
        p <- layout(p,height = 500)
        
      }else{
        p <- layout(p,height = 300)
        
      }
      
      
    }else{
      p <- layout(p,
                  title = paste0("GTEX"), 
                  margin=mar,xaxis = a)%>%
        colorbar(limits = c(0, boundary()))
      if(nrow(matrix)>20){
        
        p <- layout(p,height = 800)
        
      }else if(nrow(matrix)>10){
        p <- layout(p,height = 500)
        
      }else{
        p <- layout(p,height = 300)
        
      }
      
    }
    

      
})
  
  output$heatmap_hpa  <- renderPlotly({
    
    a <- list(
      showticklabels = TRUE,
      tickangle = -45
    )
    
    
    matrix<-heatmapInput_hpa()
    p <-plot_ly(z = matrix, type = "heatmap", 
                colorscale= "Hot",
                x =colnames(matrix),
                y =rownames(matrix),
                reversescale=TRUE)
    
    
    if(input$logarithmicY_heat){
      
      p <- layout(p,
                  title = paste0('HPA'), 
                  margin=mar,xaxis = a)%>%
        colorbar(limits = c(0, boundary()))
      
      if(nrow(matrix)>20){
        
        p <- layout(p,height = 800)
        
      }else if(nrow(matrix)>10){
        p <- layout(p,height = 500)
        
      }else{
        p <- layout(p,height = 300)
        
      }
      
      
    }else{
      p <- layout(p,
                  title = paste0("HPA"), 
                  margin=mar,xaxis = a)%>%
        colorbar(limits = c(0, boundary()))
      if(nrow(matrix)>20){
        
        p <- layout(p,height = 800)
        
      }else if(nrow(matrix)>10){
        p <- layout(p,height = 500)
        
      }else{
        p <- layout(p,height = 300)
        
      }
      
    }
    
    
})
  
output$heatmap_ccle  <- renderPlotly({

    
    a <- list(
      showticklabels = TRUE,
      tickangle = -45
    )
    
 
    matrix<-heatmapInput_ccle()
    p <-plot_ly(z = matrix, type = "heatmap", 
                colorscale= "Hot",
                x =colnames(matrix),
                y =rownames(matrix),
                reversescale=TRUE)
    
    
    if(input$logarithmicY_heat_ccle){
     
      p <- layout(p,
                  title = paste0('GTEX'), 
                  margin=mar,xaxis = a)%>%
        colorbar(limits = c(0, boundary_2()))
      
      if(nrow(matrix)>20){
        
        p <- layout(p,height = 800)
        
      }else if(nrow(matrix)>10){
        p <- layout(p,height = 500)
        
      }else{
        p <- layout(p,height = 300)
        
      }
      
      
    }else{
      p <- layout(p,
                  title = paste0("CCLE"), 
                  margin=mar,xaxis = a)%>%
        colorbar(limits = c(0, boundary_2()))
      if(nrow(matrix)>20){
        
        p <- layout(p,height = 800)
        
      }else if(nrow(matrix)>10){
        p <- layout(p,height = 500)
        
      }else{
        p <- layout(p,height = 300)
        
      }
      
    }
    
  })

  
  #scatter plot limit alert
  output$alertOne <- renderText({if (length(unique(dataInput()$cc)) > 45){print("You have exceeded max limit of input for GTEx voilin plot")}})
  output$alertTwo <- renderText({if (length(unique(dataInput_hpa()$cc)) > 45){print("You have exceeded max limit of input for HPA voilin plot")}})
  output$alertThree <- renderText({if (length(unique(dataInput_hpa()$cc)) > 45){print("You have exceeded max limit of input for plot")}})
  
  
  #scatter plot data preprocessing
  gtex_new_data <- reactive ({
    dataInput<-dataInput()
    dataInput[is.na(dataInput)] = 0
    sortedData_1(dataInput)
    
    })
  hpa_new_data <- reactive ({
    dataInput<-dataInput_hpa()
    dataInput[is.na(dataInput)] = 0
    sortedData_1(dataInput)
    })
  ccle_new_data <- reactive ({
    dataInput<-dataInput_ccle()
    dataInput[is.na(dataInput)] = 0
    sortedData_1(dataInput)
  
  })
 

  lvlsData <- reactive ({sortedData_2(gtex_new_data(), 'GeneName','Tissue',input$control_gene,input$control_tissue)})
  lvlsData_hpa <- reactive ({sortedData_2(hpa_new_data(), 'GeneName','Tissue',input$control_gene,input$control_tissue)})
  lvlsData_ccle <- reactive ({sortedData_2(ccle_new_data(), 'GeneName','Name',input$control_gene_ccle,input$control_tissue_ccle)})
  
  #scatter plot boundary limit
  boundary_scatter <- reactive({
    
    gtex_df<-gtex_new_data()
    hpa_df<-hpa_new_data()
    gtex_max<-max(gtex_df$yaxis_value)
    hpa_max<-max(hpa_df$yaxis_value)
    
    if(gtex_max>hpa_max){
      return(gtex_max)
    }else{
      return(hpa_max)
    }
    
  })
  boundary_scatter_log <- reactive({
   
    gtex_df<-gtex_new_data()
    hpa_df<-hpa_new_data()
 
    log_gtex_max<-max(gtex_df$log_yaxis_value)
    log_hpa_max<-max(hpa_df$log_yaxis_value)
    
    if(log_gtex_max>log_hpa_max){
      return(log_gtex_max+2)
    }else{
      return(log_hpa_max+2)
    }
    
  })
  
  
  output$scatter <- renderPlotly({
    
    if(length(unique(dataInput()$cc)) <= 45){
      
        ax = list(showticklabels = TRUE,
              tickangle = -45,title = "")
    

          p<-draw_scatter(dataInput(),input$logarithmicY_box,input$rank,lvlsData(),gtex_new_data())
          
          if(input$logarithmicY_box){
                  p<-layout(p,
                            title = paste0("GTEx"), 
                            xaxis = ax,
                            yaxis = list(title =paste0("TPM"),
                                         range = c(-5, boundary_scatter_log())),
                            margin=mar)
          }else{
            
            p<-layout(p,
                      title = paste0("GTEx"), 
                      xaxis = ax,
                      yaxis = list(title =paste0("TPM"),
                                   range = c(0, boundary_scatter())),
                      margin=mar)
            
          }
    }
    
  })
  
  output$scatter_hpa <- renderPlotly({
    if(length(unique(dataInput_hpa()$cc)) <= 45){
      
          ax = list(showticklabels = TRUE,
              tickangle = -45,title = "")
    
          p<- draw_scatter(dataInput_hpa(),input$logarithmicY_box,input$rank,lvlsData_hpa(),hpa_new_data())
          
          if(input$logarithmicY_box){
            
            p<-layout(p,
                      title = paste0("HPA"), 
                      xaxis = ax,
                      yaxis = list(title =paste0("TPM"),
                                   range = c(-5, boundary_scatter_log())),
                      margin=mar)
          }else{
            
            p<-layout(p,
                      title = paste0("HPA"), 
                      xaxis = ax,
                      yaxis = list(title =paste0("TPM"),
                                   range = c(0, boundary_scatter())),
                      margin=mar)
      
    }
    }
    
  })
  output$scatter_ccle <- renderPlotly({draw_scatter(dataInput_ccle(),input$logarithmicY_box,input$rank,lvlsData_ccle(),ccle_new_data())})
  
  
  output$heatmap_ui <- renderUI({
    
    df<-heatmapInput()
    nTissue= nrow(df)
    nGene = ncol(df)

    if((nGene<=10)&&(nTissue>10) ){
      column(6,withSpinner(plotlyOutput("heatmap_gtex")),style='padding-bottom:500px;')
      
    }else if((nGene<=10)&&(nTissue<=10)){
      column(6,withSpinner(plotlyOutput("heatmap_gtex")))
    }
    else{
      column(12,withSpinner(plotlyOutput("heatmap_gtex")))
      
    }
    
  })
  
  output$heatmap_hpa_ui <- renderUI({
    
    df<-heatmapInput_hpa()
    
    nTissue= nrow(df)
    nGene = ncol(df)
  
    nTissue= nrow(df)
    nGene = ncol(df)
    if((nGene<=10)&&(nTissue>10) ){
      column(6,withSpinner(plotlyOutput("heatmap_hpa")),style='padding-bottom:500px;')
      
    }else if((nGene<=10)&&(nTissue<=10)){
      column(6,withSpinner(plotlyOutput("heatmap_hpa")))
    }
    else{
      column(12,withSpinner(plotlyOutput("heatmap_hpa")))
      
    }
  })
  
  output$heatmap_ccle_ui <- renderUI({
  
    df<-heatmapInput_ccle()
    nTissue= nrow(df)
    nGene = ncol(df)
    if((nGene<=10)&&(nTissue>10) ){
      column(6,withSpinner(plotlyOutput("heatmap_ccle")),style='padding-bottom:500px;')
      
    }else if((nGene<=10)&&(nTissue<=10)){
      column(6,withSpinner(plotlyOutput("heatmap_ccle")))
    }
    else{
      column(12,withSpinner(plotlyOutput("heatmap_ccle")))
      
    }
    
  })
  
  
  
  
}

# Run the app ----This line needs to be the last line in your file.
shinyApp(ui = ui, server = server)

