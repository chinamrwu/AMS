
####### Loading Data
 


###


shinyUI(pageWithSidebar(
  #  Application title
  headerPanel('Differential Methylation Regions on Genome of '),
  sidebarPanel(width=3,
    #tags$style(HTML('hr {border-top: 2px solid black;width:100%}')),
    sliderInput("adjCancer","癌旁Beta不超过:", min=0,     max=0.5,   value=0.25, step=0.05),
    sliderInput("dltB", "dltBeta不小于:",      min = 0.1,  max =1,   value = 0.2,step=0.05),
	 tabPanel('test',plotOutput("rocPlot")),
 	 tabPanel('buttons',
	    fluidRow(
		   downloadButton("rocButton", "download ROC "),
			downloadButton("areaButton", "download Area Plot"),
			downloadButton("DMRButton", "download DMR ")
		 )
		
	)),

  mainPanel(fluidRow(dataTableOutput("DMR1")),
					tags$hr(style="border-color: blue;"),
					fluidRow(plotOutput('areaPlot'))
		  )

))

###########################################################################################