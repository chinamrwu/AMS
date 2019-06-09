
####### Loading Data
 


###


shinyUI(pageWithSidebar(
  
  #  Application title
  headerPanel("宫颈内膜癌差异甲基化(DMR)报告"),
  sidebarPanel(width=3,
    tags$style(HTML('hr {border-top: 2px solid black;width:100%}')),
    selectInput("cancerName", "Choose the disease name:", 
                choices = c('结直肠癌', 
					             '宫颈癌', 
									 '食管癌',
									 '子宫内膜癌'
									 )
	 ),
	 selectInput("site", "选择一个病灶位置:", 
                choices = c("rock", "pressure", "cars")),
	 #hr(),
    sliderInput("adjCancer","癌旁Beta不超过:", min=0, max=0.5, value=0.1,step=0.05),
    sliderInput("Cancer", "癌症Beta不小于:",   min= 0, max=1, value = 0.3, step= 0.05),
    sliderInput("dltB", "dltBeta不小于:", min = 0.2, max = 1, value = 0.25,step=0.05),
	 tags$hr(style="border-color: red;"),
	 tabPanel('test',plotOutput("rocPlot"))
   ),

  # Show a table summarizing the values entered
  #mainPanel(
  #  tabsetPanel(
  #      tabPanel("Sample Information", tableOutput("sampleInf")), 
  #      tabPanel("DMR", dataTableOutput("DMR")), 
  #  )
  #)
  mainPanel(
    tabsetPanel(
        tabPanel("Sample Information", tableOutput("sampleInf")), 
        tabPanel("DMR",  
					fluidRow(dataTableOutput("DMR1")),
					tags$hr(style="border-color: blue;"),
					fluidRow(
		           plotOutput('rocTest')
					)
		  )
    )
  )

))

###########################################################################################