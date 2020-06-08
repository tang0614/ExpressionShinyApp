library(dplyr)
library(tidyr)
library(tibble)
library(gtools)
library(shiny)
library(shinyWidgets)
library(shinythemes)
library(shinycssloaders)
library(scales)
library(DT)
library(treemap)
library(gtools)
library(highcharter)
library(stringr)
library(lexicon)
library(visNetwork)
library(httr)
library(rdrop2)
library(lubridate)
library(ggvis)
library(rsconnect)
library(plotly)
library(reshape2)
library(webshot)
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
hpa_cellline= read.csv('cellline_hpa_autophagy.csv', header = TRUE,stringsAsFactors=FALSE)


#Create variables for dropdown input 
yaxis_vars <- c('TPM'='TPM')
yaxis_vars_hpa <- c('TPM'='hpa_TPM')
yaxis_vars_ccle <- c('TPM'='TPM')

default_gene = c('ULK1','ATG14','ATG101')
default_tissue = c('LUNG','SKIN')

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
  
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
      ),
      #Setup navigation bar
      navbarPage("Autophagy Gene Expression", theme = shinytheme("spacelab"),
              #First Page
              tabPanel("Browse Tissue and Gene", fluid = TRUE, icon = icon("globe-americas"),tags$style(button_color_css),
                    # Sidebar layout with a input and output definitions
                    sidebarLayout(
                            #Input Panel
                            sidebarPanel(
                                  #Input Panel Title
                                  titlePanel("Visualize Gene Expression by Gene and Tissue"),
                                  hr(),
                                  helpText("Choose from GTEx, HPA Dataset"),
                                 fluidRow(
                                          column(8, 
                                                 selectInput("gene", label = "Gene Symbol", 
                                                               choices=mixedsort(as.vector(unique(rna_tissue$GeneName))), selected=default_gene,multiple=TRUE)
                                                  
                                          ),
                                          br(),
                                          br(),
                                          column(4,prettyCheckbox("gene_selectall", "Select All"))
                                  ),
                                  
                                  
                                  br(),
                                  fluidRow(
                                          column(8,
                                                 selectInput("tissue", label = "Tissue Name", 
                                                                choices=mixedsort(as.vector(unique(rna_tissue$Tissue))), selected=default_tissue,multiple=TRUE)
         
                                         ),
                                         br(),
                                         br(),
                                         column(4,prettyCheckbox("tissue_selectall", "Select All"))
                                  ),
  
                                  br(), 
                                  br(),  
                                  hr(),
                                  br(),
                                  fluidRow(
                                    column(10,span(textOutput("alertOne"),style="color:red"))
                                  ),  
                                  
                                  br(),
                                  fluidRow(
                                    column(10,span(textOutput("alertTwo"),style="color:red"))
                                  ),  
                                 
                                      
                                  br(),
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
                                     helpText("For Heatmap: NaN values (missing values) are represented as empty boxes."),
                                     helpText("“To generate new violin plots plots without changing the heatmaps, click ‘Keep All Heatmaps”")
                              )
                            ),
  
                            br(),
                            fluidRow(
                                    column(2, offset = 1,
                                           prettyCheckbox("logarithmicY_heat", "Log2 of Y", FALSE,
                                                          shape = "round", 
                                                          bigger=TRUE,
                                                          outline =TRUE, 
                                                          animation = "smooth")
                                    ),
                                    column(2,offset = 4,
                                           prettyCheckbox("keep_heat", "Keep All Heatmaps")
                                    )
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
                                            helpText("For Violin Plot: Sum of selected gene and tissue should be lower than 45 for violin plot."),
                                            helpText("Click \"Median sort with tissue-types\" and \"Median sort within genes\" will sort all selected samples by median TPM."),
                                            helpText("To generate new heatmap plots without changing the violin plots, click ‘Keep All Violin Plots”")
                                            
                                     )
                                   ),
                                    
                                    fluidRow(
                                            column(1, 
                                                   prettyCheckbox("rank", "", FALSE,
                                                                  shape = "round", 
                                                                  bigger=FALSE)
                                            ),
                                            column(4,
                                                   
                                                   prettyCheckbox("control_tissue", "Median sort with tissue-types", FALSE,
                                                                  shape = "curve", 
                                                                  bigger=FALSE,
                                                                  animation = "smooth")
                                                   
                                            ),
                                            column(4,
                                                   prettyCheckbox("control_gene", "Median sort within genes", FALSE, 
                                                                  shape = "curve", 
                                                                  bigger=FALSE,
                                                                  animation = "smooth")
                                            )
                                            
                                      
                                  ),
    
                                fluidRow(
                                          column(4, offset = 1,
                                                 prettyCheckbox("logarithmicY_box", "Log2 of Y", FALSE,
                                                                shape = "round", 
                                                                bigger=TRUE,
                                                                outline =TRUE, 
                                                                animation = "smooth")
                                          ),
                                          
                                          column(2, 
                                                 prettyCheckbox("keep_plot", "Keep All ViolinPlots")
                                          )
                                ),

                                fluidRow(column(4,offset=4,h3('TPM violin Plot')))
                                
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
                          fluidRow(
                            
                                  column(11, withSpinner(plotlyOutput("scatter_hpa")))
                          )

                      )
              )
      ),
      #Second Page
      tabPanel("Browse Cell-line and Gene", fluid = TRUE, icon = icon("globe-americas"),tags$style(button_color_css),
  
               # Sidebar layout with a input and output definitions
               sidebarLayout(
  
                 #Input Panel
                 sidebarPanel(
                   #Input Panel Title
                   titlePanel("Visualize Gene Expression by Gene and cell-line"),
                   hr(),
                   helpText("Choose from CCLE, HPA Dataset"),
              
                   fluidRow(

                     column(8, 
                            
                            selectInput("gene_cellline", label = "Gene Symbol", 
                                        choices=mixedsort(as.vector(unique(ccle$GeneName))), selected=default_gene,multiple=TRUE
                            )
                            
                     ),
                     br(),
                     br(),
                     column(4,prettyCheckbox("gene_selectall_cellline", "Select All"))
                   ),
  
                   br(),
                   fluidRow(
                     column(8,
                            selectInput("cellline_ccle", label = "Cell-line (CCLE)", 
                                        choices=mixedsort(as.vector(unique(ccle$cellline))), selected=c('A101D','A172','59M','697'),multiple=TRUE
                            )
                     ),
                     br(),
                     br(),
                     column(4,prettyCheckbox("cellline_ccle_selectall", "Select All"))
                   ),

                   fluidRow(
                     column(8,
                            selectInput("cellline_hpa", label = "Cell-line (HPA)", 
                                        choices=mixedsort(as.vector(unique(hpa_cellline$cellline))), selected=c('AF22','A-431','Daudi'),multiple=TRUE
                            )
                            
                            
                    ),
                    br(),
                    br(),
                    column(4,prettyCheckbox("cellline_hpa_selectall", "Select All"))
                   ),
                   br(), 
                   br(),  
                   hr(),
                   br(),
                   fluidRow(
                     column(10,span(textOutput("alertThree"),style="color:red"))
                   ),  
                   
                   br(),
                   fluidRow(
                     column(10,span(textOutput("alertFour"),style="color:red"))
                   ),  
                   
                   br(),
                   tags$p(span("Large graphs (e.g., selecting lots of genes and celllines) may take a few seconds to render.", style = "color:red")),
                   tags$p(HTML("<b>Transcripts Per Million (TPM)</b> is a normalization method for RNA-seq, should be read as: for every 1,000,000 RNA molecules in the RNA-seq sample, x came from this gene.")),
                   #tags$p(HTML("<b>pTPM</b> is calculated from scaling a sum of 1 million TPM to compensate for the non-coding transcripts that had been previously removed.")),
                   tags$p(HTML("<b>CCLE Description</b> : Dataset comes from the latest release <a href=\"https://portals.broadinstitute.org/ccle/about\">CCLE (Cancer Cell Line Encyclopedia) </a>project, which
                   is a collaboration between the Broad Institute, and the Novartis Institutes for Biomedical
                   Research and its Genomics Institute of the Novartis Research Foundation to conduct a detailed genetic 
                   and pharmacologic characterization of a large panel of human cancer models, to develop integrated computational analyses 
                   that link distinct pharmacologic vulnerabilities to genomic patterns and to translate cell line integrative genomics into cancer 
                   patient stratification. The CCLE provides public access to genomic data, analysis and visualization for over 1100 cell lines.")),
                   
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
                            helpText("For Heatmap: NaN values (missing values) are represented as empty boxes."),
                            helpText("To generate new bar plots plots without changing the heatmaps, click ‘Keep All Heatmaps”")
                     )
                   ),

                   br(),
                   
                   fluidRow(
                     
                     column(2, 
                            offset = 1,
                            prettyCheckbox("logarithmicY_heat_celline", "Log2 of Y", FALSE,
                                           shape = "round", 
                                           bigger=TRUE,
                                           outline =TRUE, 
                                           animation = "smooth")
                     ),
                     column(4, offset = 3,
                            prettyCheckbox("keep_heat_cell", "Keep All Heatmaps")
                     )
                     
                   ),
                   
                   
                   br(),
                   
                   fluidRow(column(4,offset=4,h3('TPM Heatmap'))),
                   
                   uiOutput("heatmap_ccle_ui"),
                   uiOutput("heatmap_hpac_ui"),
                   
                   
                   br(),
                   hr(),
                   
                   column(12,
                          fluidRow(
                            column(10,offset = 1,
                                   helpText("For Bar Plot: Sum of selected gene and tissue should be lower than 45 for bar plot."),
                                   helpText("Click \"Median sort with celllines\" and \"Median sort within genes\" will sort all selected samples by median TPM."),
                                   helpText("To generate new heatmaps without changing the bar plots, click ‘Keep All Bar Plots”")
                                   
                            )
                          ),
                          
                            fluidRow(
                            
                            
                            column(1, 
                                   prettyCheckbox("rank_ccle", "", FALSE,
                                                  shape = "round", 
                                                  bigger=FALSE)
                            ),
                            column(4,
                                   prettyCheckbox("control_cellline_ccle", "Median sort with celllines", FALSE,
                                                  shape = "curve", 
                                                  bigger=FALSE,
                                                  animation = "smooth")
                            ),
                            column(4,
                                   prettyCheckbox("control_gene_ccle", "Median sort within genes", FALSE, 
                                                  shape = "curve", 
                                                  bigger=FALSE,
                                                  animation = "smooth")
                            )
                            
                            
                          ),
                          
                          fluidRow(column(4, offset = 1,
                                          prettyCheckbox("logarithmicY_box_ccle", "Log2 of Y", FALSE,
                                                         shape = "round", 
                                                         bigger=TRUE,
                                                         outline =TRUE, 
                                                         animation = "smooth")),
                                   column(4, 
                                          prettyCheckbox("keep_plot_cell", "Keep All Bar Plots")
                                   )
                          ),
                          fluidRow(column(4,offset=4,h3('TPM Bar Plot'))
                                   )
                     
                         
                   ),
            
                   fluidRow(column(11, withSpinner(plotlyOutput("scatter_ccle")))),  
                   
                   br(),
                   br(),
                   br(),
                   br(),
                   br(),
                   br(),
                 fluidRow(column(11, withSpinner(plotlyOutput("scatter_hpac")))
      
                   
                 )
        )
      )
                   

)))

      


# Define server logic ----how to assenble input into output

server <- function(input, output,session) {
  
 
  #select button
  observe({
    if (input$tissue_selectall){
      updateSelectInput(session, "tissue", selected = as.vector(unique(rna_tissue$Tissue)))
    }else{
      updateSelectInput(session, "tissue", selected =  default_tissue)
      
    }
    
  })
  
  observe({
    if (input$gene_selectall){
      updateSelectInput(session, "gene", selected = as.vector(unique(rna_tissue$GeneName)))
    }else{
      updateSelectInput(session, "gene", selected =  default_gene)
    }
    
  })
  
  observe({
    
    if (input$gene_selectall_cellline){
      updateSelectInput(session, "gene_cellline", selected = as.vector(unique(ccle$GeneName)))
    }else{
      updateSelectInput(session, "gene_cellline", selected = default_gene)
    }
   
    
  })
  
  
  observe({
    if (input$cellline_ccle_selectall){
      updateSelectInput(session, "cellline_ccle", selected = as.vector(unique(ccle$cellline)))
    }else{
      updateSelectInput(session, "cellline_ccle", selected = c('A101D','A172','59M','697'))
    }
  })
  

  observe({
          if (input$cellline_hpa_selectall){
            updateSelectInput(session, "cellline_hpa", selected = as.vector(unique(hpa_cellline$cellline)))
          }else{
            updateSelectInput(session, "cellline_hpa", selected = c('AF22','A-431','Daudi'))
          }
  })
  observe({
    if(input$control_gene|input$control_tissue){
      updateCheckboxInput(session, "rank", value = TRUE)
    }
    
  })
  
  
  observe({
      if(input$control_gene_ccle|input$control_cellline_ccle){
        updateCheckboxInput(session, "rank_ccle", value = TRUE)
      }
    
 
    
})
  
 

  #GTEx data input
  dataInput <- reactive({
    m<-get_data(rna_tissue,'gtex_TPM')
    if(input$keep_plot){
       m<-filter_data(isolate(input$tissue),isolate(input$gene),m,'GeneName','Tissue')
    }else{
       m<-filter_data(input$tissue,input$gene,m,'GeneName','Tissue')
    }
    m
  })
  dataInput_heat <- reactive({
    m<-get_data(rna_tissue,'gtex_TPM')
    if(input$keep_heat){
      m<-filter_data(isolate(input$tissue),isolate(input$gene),m,'GeneName','Tissue')
    }else{
      m<-filter_data(input$tissue,input$gene,m,'GeneName','Tissue')
    }
    m
  })
  
  #Datainput hpa
  dataInput_hpa <-reactive({
    m<-get_data(rna_tissue_hpa,'hpa_TPM')
    if(input$keep_plot){
      m<-filter_data(isolate(input$tissue),isolate(input$gene),m,'GeneName','Tissue')
    }else{
      m<-filter_data(input$tissue,input$gene,m,'GeneName','Tissue')
    }
    m
  })
  
  dataInput_hpa_heat <-reactive({
    m<-get_data(rna_tissue_hpa,'hpa_TPM')
    if(input$keep_heat){
      m<-filter_data(isolate(input$tissue),isolate(input$gene),m,'GeneName','Tissue')
    }else{
      m<-filter_data(input$tissue,input$gene,m,'GeneName','Tissue')
    }
    m
  })
  
  #Datainput ccle
  dataInput_ccle <-reactive({
    m<-get_data(ccle,'TPM')
    if(input$keep_plot_cell){
      m<-filter_data(isolate(input$cellline_ccle),isolate(input$gene_cellline),m,'GeneName','cellline')
    }else{
      m<-filter_data(input$cellline_ccle,input$gene_cellline,m,'GeneName','cellline')
    }
    m
  })
  
  dataInput_ccle_heat <-reactive({
    m<-get_data(ccle,'TPM')
    if(input$keep_heat_cell){
      m<-filter_data(isolate(input$cellline_ccle),isolate(input$gene_cellline),m,'GeneName','cellline')
    }else{
      m<-filter_data(input$cellline_ccle,input$gene_cellline,m,'GeneName','cellline')
    }
    m
  })
  
  #Datainput hpa cellline
  dataInput_hpac <-reactive({
    m<-get_data(hpa_cellline,'TPM')
    if(input$keep_plot_cell){
      m<-filter_data(isolate(input$cellline_hpa),isolate(input$gene_cellline),m,'GeneName','cellline')
    }else{
      m<-filter_data(input$cellline_hpa,input$gene_cellline,m,'GeneName','cellline')
    }
    m
  })
  
  dataInput_hpac_heat <-reactive({
    m<-get_data(hpa_cellline,'TPM')
    if(input$keep_heat_cell){
      m<-filter_data(isolate(input$cellline_hpa),isolate(input$gene_cellline),m,'GeneName','cellline')
    }else{
      m<-filter_data(input$cellline_hpa,input$gene_cellline,m,'GeneName','cellline')
    }
    m
  })
  

  
   
  
#Heatmap Input

  heatmapInput <-reactive({heatmap_input(dataInput_heat(),'Tissue','GeneName',input$logarithmicY_heat)})
  heatmapInput_hpa <-reactive({heatmap_input(dataInput_hpa_heat(),'Tissue','GeneName',input$logarithmicY_heat)})
  heatmapInput_ccle <-reactive({heatmap_input(dataInput_ccle_heat(),'cellline','GeneName',input$logarithmicY_heat_celline)})
  heatmapInput_hpac <-reactive({heatmap_input(dataInput_hpac_heat(),'cellline','GeneName',input$logarithmicY_heat_celline)})

  # draw heatmap
  .heatmap_gtex<-reactive({
    df<-heatmapInput()
    
    df$z_hover <-df$m
    df$z_hover[df$z_hover<=-100]  <- "Missing" 
    
    a <- list(
      showticklabels = TRUE,
      tickangle = -45
    )
    
    p <-plot_ly(data = df, type = "heatmap", 
                colorscale= "Hot",
                x =df$GeneName,
                y =df$Tissue,
                z=df$m,
                reversescale=TRUE,
                hoverinfo = "text",
                hovertext = paste(
                  "<br>Gene :", df$GeneName,
                  "<br>Tissue :", df$Tissue,
                  "<br>Median TPM :", df$z_hover))
    
    p <- layout(p,
                title = paste0('GTEx'), 
                margin=mar,xaxis = a)%>%
      colorbar(limits = c(-10, boundary()))
    
    l<-length(unique(df$Tissue))
    if(l>20){
      
      p <- layout(p,height = 800)
      
    }else if(l>10){
      p <- layout(p,height = 500)
      
    }else{
      p <- layout(p,height = 300)
      
    }
    
    
    
  })
  output$heatmap_gtex <-renderPlotly({
    .heatmap_gtex()
  })
  
  output$heatmap_hpa  <- renderPlotly({
    
    a <- list(
      showticklabels = TRUE,
      tickangle = -45
    )
   
    
    df<-heatmapInput_hpa()
    df$z_hover <-df$m
    df$z_hover[df$z_hover<=-100]  <- "Missing" 
    
    p <-plot_ly(data = df, type = "heatmap", 
                colorscale= "Hot",
                x =df$GeneName,
                y =df$Tissue,
                z=df$m,
                reversescale=TRUE,
                hoverinfo = "text",
                hovertext = paste(
                  "<br>Gene :", df$GeneName,
                  "<br>Tissue :", df$Tissue,
                  "<br>Median TPM :", df$z_hover))
    
    p <- layout(p,
                title = paste0('HPA'), 
                margin=mar,xaxis = a)%>%
      colorbar(limits = c(-10, boundary()))
    l<-length(unique(df$Tissue))
    
    if(l>20){
      
      p <- layout(p,height = 800)
      
    }else if(l>10){
      p <- layout(p,height = 500)
      
    }else{
      p <- layout(p,height = 300)
      
    }
    
    

  })
  
  output$heatmap_ccle  <- renderPlotly({
  
    a <- list(
      showticklabels = TRUE,
      tickangle = -45
    )
    b<- list(
      showticklabels = TRUE,
      tickangle = 45
    )
    
    df<-heatmapInput_ccle()
    df$z_hover <-df$m
    df$z_hover[df$z_hover<=-100]  <- "Missing" 

    p <-plot_ly(data = df, type = "heatmap", 
                colorscale= "Hot",
                x =df$GeneName,
                y =df$cellline,
                z=df$m,
                reversescale=TRUE,
                hoverinfo = "text",
                hovertext = paste(
                  "<br>Gene :", df$GeneName,
                  "<br>Cellline :", df$cellline,
                  "<br>TPM :", df$z_hover))
    
    p <- layout(p,
                title = paste0('CCLE'), 
                margin=mar,xaxis = a)%>%
      colorbar(limits = c(-10, boundary_cell()))
    
    l<-length(unique(df$cellline))
    if(l>20){
      
      p <- layout(p,height = 800)
      
    }else if(l>10){
      p <- layout(p,height = 500)
      
    }else{
      p <- layout(p,height = 300)
      
    }
      
      
   
  })
  
  output$heatmap_hpac  <- renderPlotly({
    
    a <- list(
      showticklabels = TRUE,
      tickangle = -45
    )
    
    b<- list(
      showticklabels = TRUE,
      tickangle = 45
    )

    df<-heatmapInput_hpac()
    df$z_hover <-df$m
    df$z_hover[df$z_hover<=-100]  <- "Missing" 
    
    p <-plot_ly(data = df, type = "heatmap", 
                colorscale= "Hot",
                x =df$GeneName,
                y =df$cellline,
                z=df$m,
                reversescale=TRUE,
                hoverinfo = "text",
                hovertext = paste(
                  "<br>Gene :", df$GeneName,
                  "<br>Cellline :", df$cellline,
                  "<br>TPM :", df$z_hover))
    
    p <- layout(p,
                title = paste0('HPA'), 
                margin=mar,xaxis = a)%>%
      colorbar(limits = c(-10, boundary_cell()))
    l<-length(unique(df$cellline))
    if(l>20){
      
      p <- layout(p,height = 800)
      
    }else if(l>10){
      p <- layout(p,height = 500)
      
    }else{
      p <- layout(p,height = 300)
      
    }
  })
  
  
  
  
  #scatter plot limit alert
  output$alertOne <- renderText({if (length(unique(dataInput()$cc)) > 45){print("You have exceeded max limit of input for GTEx violin plot")}})
  output$alertTwo <- renderText({if (length(unique(dataInput_hpa()$cc)) > 45){print("You have exceeded max limit of input for HPA violin plot")}})
  output$alertThree <- renderText({if (length(unique(dataInput_ccle()$cc)) > 45){print("You have exceeded max limit of input for CCLE bar plot")}})
  output$alertFour <- renderText({if (length(unique(dataInput_hpac()$cc)) > 45){print("You have exceeded max limit of input for HPA bar plot")}})
  
  
  
  #scatter plot data preprocessing
  gtex_new_data <-reactive({
    dataInput<-dataInput()
    dataInput[is.na(dataInput)] = 0
    sortedData_1(dataInput)
    
  })
  hpa_new_data <-reactive({
    dataInput<-dataInput_hpa()
    dataInput[is.na(dataInput)] = 0
    sortedData_1(dataInput)
  })
  
  ccle_new_data <-reactive({
    
    dataInput<-dataInput_ccle()
    dataInput[is.na(dataInput)] = 0
    sortedData_1(dataInput)
  })
  
  hpac_new_data <- reactive({
    dataInput<-dataInput_hpac()
    dataInput[is.na(dataInput)] = 0
    sortedData_1(dataInput)
 
  })
  
  #boundary for tissue
  boundary <- reactive({get_boundary(heatmapInput(),heatmapInput_hpa(),'m')})
  boundary_cell <- reactive({get_boundary(heatmapInput_ccle(),dataInput_ccle_heat(),'m')})
  boundary_scatter <- reactive({get_boundary(gtex_new_data(),hpa_new_data(),'yaxis_value')})
  boundary_scatter_log <- reactive({get_boundary(gtex_new_data(),hpa_new_data(),'log_yaxis_value')})
  boundary_scatter_cell <- reactive({get_boundary(ccle_new_data(),hpac_new_data(),'yaxis_value')})
  boundary_scatter_log_cell <- reactive({get_boundary(ccle_new_data(),hpac_new_data(),'log_yaxis_value')})
  
  
  scatter_ccle_var <- reactive({
    d<-ccle_new_data()
    
    sub_list1<- d$cellline %in% input$celline
    m <-subset(d, sub_list1)
    as.data.frame(m)
    
  })
  scatter_hpa_var <- reactive({
    d<-hpac_new_data()
    
    sub_list1<- d$cellline %in% input$celline_hpa
    m <-subset(d, sub_list1)
    as.data.frame(m)
    
  })
  
  
  
  lvlsData <- reactive ({sortedData_2(gtex_new_data(), 'GeneName','Tissue',input$control_gene,input$control_tissue)})
  lvlsData_hpa <- reactive ({sortedData_2(hpa_new_data(), 'GeneName','Tissue',input$control_gene,input$control_tissue)})
  lvlsData_ccle <- reactive ({sortedData_2(ccle_new_data(), 'GeneName','cellline',input$control_gene_ccle,input$control_cellline_ccle)})
  lvlsData_hpac <- reactive ({sortedData_2(hpac_new_data(), 'GeneName','cellline',input$control_gene_ccle,input$control_cellline_ccle)})
  
  
  
  output$cellline_ui <- renderUI({
    column(3,offset = 7,
            selectInput("celline", label = "CCLE Cell-line", 
                        choices=as.vector(unique(dataInput_ccle()$cellline)), selected='',multiple=TRUE)
    )
  })
  
  output$cellline_ui_hpa <- renderUI({
    column(3,offset = 7,
           selectInput("celline_hpa", label = "HPA Cell-line", 
                       choices=as.vector(unique(dataInput_hpac()$cellline)), selected='',multiple=TRUE)
    )
  })
  

  output$scatter <- renderPlotly({
   
    ax = list(showticklabels = TRUE,
              tickangle = -45,title = "")
   
    if(length(unique(dataInput()$cc)) <= 45){
      
          p<-draw_scatter(dataInput(),input$logarithmicY_box,input$rank,lvlsData(),gtex_new_data())
          if(input$logarithmicY_box){
              p<-layout(p,
                        title = paste0("GTEx"), 
                        xaxis = ax,
                        yaxis = list(title =paste0("Log2 TPM"),
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
    
    ax = list(showticklabels = TRUE,
              tickangle = -45,title = "")
 
    if(length(unique(dataInput_hpa()$cc)) <= 45){
      
    
      p<- draw_scatter(dataInput_hpa(),input$logarithmicY_box,input$rank,lvlsData_hpa(),hpa_new_data())
      
      if(input$logarithmicY_box){
        
        p<-layout(p,
                  title = paste0("HPA"), 
                  xaxis = ax,
                  yaxis = list(title =paste0("Log2 TPM"),
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
  
  output$scatter_ccle <- renderPlotly({
    if(length(unique(dataInput_ccle()$cc)) <= 45){
    
    ax = list(showticklabels = TRUE,
              tickangle = -45,title = "")
    
    p <- draw_bar(input$logarithmicY_box_ccle,input$rank_ccle,lvlsData_ccle(),ccle_new_data(),'Site_Primary')
           
          if(input$logarithmicY_box_ccle){
                p<-p%>%layout(p,
                              title = paste0("CCLE"), 
                              xaxis = ax,
                              yaxis = list(title =paste0("Log2 TPM"),range = c(0, boundary_scatter_log_cell())),
                              margin=mar)
          }else{
                
             p<-p%>%layout(p,
                              title = paste0("CCLE"), 
                              xaxis = ax,
                              yaxis = list(title =paste0("TPM"),range = c(0, boundary_scatter_cell())),
                              margin=mar)
                
          }
  }

  })
  
  output$scatter_hpac <- renderPlotly({
  if(length(unique(dataInput_hpac()$cc)) <= 45){
   ax = list(showticklabels = TRUE,
              tickangle = -45,title = "")
   
    p <- draw_bar(input$logarithmicY_box_ccle,input$rank_ccle,lvlsData_hpac(),hpac_new_data(),'Organ')

      if(input$logarithmicY_box_ccle){
        p<-p%>%layout(p,
                      title = paste0("HPA"), 
                      xaxis = ax,
                      yaxis = list(title =paste0("Log2 TPM"),range = c(0, boundary_scatter_log_cell())),
                      margin=mar)
      }else{
        
        p<-p%>%layout(p,
                      title = paste0("HPA"), 
                      xaxis = ax,
                      yaxis = list(title =paste0("TPM"),range = c(0, boundary_scatter_cell())),
                      margin=mar)
        
      }
  }
  })
  
  
  
  
  output$heatmap_ui <- renderUI({
    
    df<-heatmapInput()
    nTissue=length(unique(df$Tissue))
    nGene<-length(unique(df$GeneName))
    
    if((nGene>=16)&&(nTissue>10) ){
      column(12,withSpinner(plotlyOutput("heatmap_gtex")),style='padding-bottom:500px;')
      
    }else 
      if((nGene<=16)&&(nTissue>10) ){
        column(5,withSpinner(plotlyOutput("heatmap_gtex")),style='padding-bottom:500px')
        
      }else if((nGene<=16)&&(nTissue<=10)){
        column(5,withSpinner(plotlyOutput("heatmap_gtex")))
      }
    else{
      column(12,withSpinner(plotlyOutput("heatmap_gtex")))
      
    }
    
  })
  
  output$heatmap_hpa_ui <- renderUI({
    
    df<-heatmapInput_hpa()
    
    nTissue=length(unique(df$Tissue))
    nGene<-length(unique(df$GeneName))
    if((nGene>=16)&&(nTissue>10) ){
      column(12,withSpinner(plotlyOutput("heatmap_hpa")),style='padding-bottom:500px;')
      
    }else 
      if((nGene<=16)&&(nTissue>10) ){
        column(5,withSpinner(plotlyOutput("heatmap_hpa")),style='padding-bottom:500px')
        
      }else if((nGene<=16)&&(nTissue<=10)){
        column(5,withSpinner(plotlyOutput("heatmap_hpa")))
      }
    else{
      column(12,withSpinner(plotlyOutput("heatmap_hpa")))
      
    }
  })
  
  output$heatmap_ccle_ui <- renderUI({
    
    df<-heatmapInput_ccle()
    
    nTissue=length(unique(df$cellline))
    nGene<-length(unique(df$GeneName))
    if((nGene>=16)&&(nTissue>10) ){
      column(12,withSpinner(plotlyOutput("heatmap_ccle")),style='padding-bottom:500px;')
      
    }else 
    if((nGene<=16)&&(nTissue>10) ){
      column(5,withSpinner(plotlyOutput("heatmap_ccle")),style='padding-bottom:500px')
      
    }else if((nGene<=16)&&(nTissue<=10)){
      column(5,withSpinner(plotlyOutput("heatmap_ccle")))
    }
    else{
      column(12,withSpinner(plotlyOutput("heatmap_ccle")))
      
    }
  })
  
  output$heatmap_hpac_ui <- renderUI({
    
    df<-heatmapInput_hpac()
    
    nTissue=length(unique(df$cellline))
    nGene<-length(unique(df$GeneName))
    
    if((nGene>=16)&&(nTissue>10) ){
      column(12,withSpinner(plotlyOutput("heatmap_hpac")),style='padding-bottom:500px;')
      
    }else 
    if((nGene<=16)&&(nTissue>10) ){
      column(6,withSpinner(plotlyOutput("heatmap_hpac")),style='padding-bottom:500px')
      
    }else if((nGene<=16)&&(nTissue<=10)){
      column(6,withSpinner(plotlyOutput("heatmap_hpac")))
    }
    else{
      column(12,withSpinner(plotlyOutput("heatmap_hpac")))
      
    }
  })
  

  
}
# Run the app ----This line needs to be the last line in your file.
shinyApp(ui = ui, server = server)