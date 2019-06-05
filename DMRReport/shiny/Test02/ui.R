
####### Loading Data
 


###


shinyUI(pageWithSidebar(

  #  Application title
  headerPanel("宫颈内膜癌差异甲基化(DMR)报告"),
  sidebarPanel(
    selectInput("disease name", "Choose the disease name:", 
                choices = c("BRCA", "CESC", "CRC")),
	 selectInput("dataset", "Choose a dataset:", 
                choices = c("rock", "pressure", "cars")),
    sliderInput("adjCancer","癌旁Beta不超过:", min=0, max=1, value=0.1,step=0.05),
    sliderInput("Cancer", "癌症Beta不小于:",   min= 0, max=1, value = 0.3, step= 0.05),
    sliderInput("dltB", "dltBeta不小于:", min = 0.1, max = 1, value = 0.25,step=0.05)
   ),

  # Show a table summarizing the values entered
  mainPanel(
    tabsetPanel(
        tabPanel("Sample Information", plotOutput("plot")), 
        tabPanel("DMR", verbatimTextOutput("summary")), 
        tabPanel("plot", tableOutput("table"))
  )
))
###########################################################################################