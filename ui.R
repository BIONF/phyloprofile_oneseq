library(shiny)
library(shinyBS)
library(DT)

shinyUI(fluidPage(
  
  # Application title
  titlePanel("Phylogenetic profile"),
  
  # Sidebar for input file
  sidebarLayout(
    sidebarPanel(
#      textOutput("testOutput"),    ### use for testing output ###
      selectInput("rankSelect", label = "Select taxonomy rank",
                  choices = list("1_genus" = "1_genus", "2_family" = "2_family", "3_order" = "3_order", "4_class" = "4_class",
                                 "5_phylum" = "5_phylum", "6_kingdom" = "6_kingdom", "7_superkingdom" = "7_superkingdom"), 
                  selected = 1),
      uiOutput("select"),     
      fileInput("file1","Upload input matrix file: ")
    ),
    
    # Main page
    mainPanel(
      tabsetPanel(
#        tabPanel ("Distribution",plotOutput("plot1")),
        tabPanel ("Present/absent profile",uiOutput("plot.ui"),
                  bsModal("modalBS", "Detailed plot", "go", size = "large",
                          uiOutput("detailPlot.ui"),
                          numericInput("detailedHeight","plot_height(px)",min=100,max=1600,step=50,value=100,width=100)
                          )
                ),
        tabPanel ("Data",dataTableOutput("dis"),
                  downloadButton('downloadData', 'Download filtered data'))
      )
    )
  ),
  
  ### show click info
  absolutePanel(
    bottom = 5, left = 30,
    fixed = TRUE,
#    wellPanel(
      h5("Point's info:"),
      verbatimTextOutput("pointInfo"),
      actionButton("go", "Plot"),
#    ),
    style = "opacity: 0.80"
  ),

  ### size of profile plot
  absolutePanel(
    bottom = 100, left = 30, width = 200,
    fixed = TRUE,
#    sliderInput("zoom", "Height", min=1, max=10, value=1, step = 0.5),
    numericInput("width","Width(px)",min=600,max=3200,step=50,value=600,width=100),
    style = "opacity: 0.92"
  ),
  absolutePanel(
    bottom = 100, left = 150, width = 200,
    fixed = TRUE,
    numericInput("height","Height(px)",min=600,max=1600,step=50,value=600,width=100),
    style = "opacity: 0.92"
  ),
  
  ### options for sub-plotting the data
  absolutePanel(
    bottom = 190, left = 30, width = 150,
    fixed = TRUE,
    numericInput("number","# rows for profile",min=1,max=1600,step=10,value=30,width=150),
    style = "opacity: 0.92"
  ),
  absolutePanel(
    bottom = 185, left = 115, width = 100,
    fixed = TRUE,
    textOutput("totalRows"),
    style = "opacity: 0.92"
  ),
  absolutePanel(
    bottom = 190, left = 190, width = 150,
    fixed = TRUE,
    numericInput("stIndex","start at:",min=1,max=1600,value=1,width=100),
    style = "opacity: 0.92"
  ),

  ### slider bar for FAS filter
  absolutePanel(
    bottom = 245, left = 30, width = 150,
    fixed = TRUE,
    sliderInput("fas",
                "FAS cutoff: ",
                min = 0,
                max = 1,
                step = 0.05,
                value = 0.0,
                width = 200),
    style = "opacity: 0.92"
  ),
  
  ### slider bar for percentage filter
  absolutePanel(
    bottom = 245, left = 185, width = 150,
    fixed = TRUE,
    sliderInput("percent",
                "Percentage of present species:",
                min = 0,
                max = 1,
                step = 0.05,
                value = 0.0,
                width = 200),
    style = "opacity: 0.92"
  )
))
