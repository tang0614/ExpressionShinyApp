#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


library(dplyr)
library(tidyr)
library(tibble)
library(gtools)
library(shiny)
library(shinyWidgets)
library(shinythemes)
library(scales)
library(gtools)
library(visNetwork)
library(rsconnect)
library(networkD3)
library(magrittr)

if (FALSE) {
    library(RSQLite)
    library(dbplyr)
}
library(shiny)

df <- read.csv("corr_rank_within70.txt", header = TRUE,stringsAsFactors = FALSE)
#default_gene = c('ATG14','ATG101','ULK1','WIPI2')


# Define UI for dataset viewer app ----
ui <- fluidPage(
  
    # App title ----
    titlePanel("Gene Network Analysis from Co-essentiality Dependencies"),

    # Sidebar layout with input and output definitions ----
    sidebarLayout(
        
        # Sidebar panel for inputs ----
        sidebarPanel(
            
            #Input Panel Title
            titlePanel("Input Panel"),
            hr(),
            helpText("Select the sign of pearson correlation."),
           
            fluidRow(
                column(11, 
                       radioButtons("corr", "Correlation:",
                                    c("Positive" = "positive",
                                      "Negative" = "negative",
                                      "All" = "all"),
                                    selected = "positive"),
                       
                )
            ),
            
            
 
            fluidRow(
              column(11, 
                     helpText("Enter capitalized gene symbol separatd by space/new lines."),
                     textAreaInput("caption", "Gene Symbol", "ATG14 ATG101 ULK1")
              ),
            
            ),

            fluidRow(
              column(11,
                     helpText("Select the top n correlated genes for each selected gene above."),
                     sliderInput("rank", "Rank:",
                                 min = 1, max = 140,
                                 value = 10)
                     
                     
              )
            ),
            

            fluidRow(
              column(11,
                     helpText("Select the numeric value indicating either the strength of the node repulsion (negative value) or attraction (positive value)."),
                     sliderInput("charge", "Charge:",
                                 min = -500, max = 5,
                                 value = -300)
                     
                     
              )
            ),
            
            helpText("Reflesh and submit the value you",
                     "entered above."),
            fluidRow(
              column(11,
                     submitButton("Submit", icon("refresh"))
                     
              )
            ),
           
            
            br(),
            hr(),
            br(),
            
            fluidRow(
                column(11, 
               
                tags$p(HTML("<b>Dataset Description</b> :This <a href=\"https://depmap.org/portal/download/\">DepMap release</a> contains data from CRISPR knockout screens from project Achilles which contains the results of genome-scale CRISPR knockout screens for 18,333 genes in 739 cell lines. These networks enable people to discover novel clusters of genes that might share functional
                            associations and connections with wach other based on their essentiality profiles from screen data across multiple cell lines.")),
                
            ))
            
            
            
     
            
        ),
        
        mainPanel(
            
            fluidRow(
               
                column(11,offset = 1,
                       h3('Network Visualization'),
                       helpText("The width of links indicates the strengh of gene correlations. Thick line represents greater positive/negative correlations."),
                       helpText("Use mouse to zoom in/out to change the size of the networks. ")
                )
            ),
  
            fluidRow(
               
  
               column(3, offset = 8,
                         uiOutput("render_ui_genes"),
                         downloadButton("downloadNet", "Download network below")
                         #conditionalPanel("output.show",  uiOutput("render_ui_genes"))
                )
            ),
            
            fluidRow(
              
              column(11,align="center",
                     tabPanel("Force Network", forceNetworkOutput("force"))
            ),
              

            
            ),
            
            
           
            
            
            br(),
            br(),
            
            br(),
            hr(),
            


            fluidRow(
                column(4,offset = 1,
                       
                       h3('Summary Table'),
                       helpText("Genes that are showing in the plot."),
                       downloadButton("downloadData", "Download Table")
                       
                ),
                
            
                
                column(4,offset = 3,
                       
                       h3('Additional Info'),
                       helpText("Genes that are not showing in the plot."),
                       downloadButton("downloadMissingData", "Download Missing Genes")
                )
              
                
            ),
            br(),
            br(),
            
           
           
            fluidRow(
                column(4, offset = 1,
                       tableOutput('table1')
                       
                ),
             
                column(4,offset = 3,
                       br(),
                       verbatimTextOutput("value", placeholder = TRUE)
                      
                ),
                
                
              
            ),
            
          
    
            
        )
        
       
    )
)

# Define server logic to summarize and view selected dataset ----
server <- function(input, output,session) {
  
  
   

}

# Create Shiny app ----
shinyApp(ui, server)