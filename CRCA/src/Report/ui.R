shinyUI(fluidPage(
   fluidRow(  
	   column(3,
				selectInput("site", "选择一个病灶位置:", 
				 choices = c(
						"overall",
						"Rectum, NOS" ,
						"Sigmoid colon",
						"Hepatic flexure of colon",
						"Cecum",
						"Descending colon",
						"Ascending colon",
						"Colon, NOS" ,             
						"Transverse colon",
						"Rectosigmoid junction"   
				 )),

				sliderInput("adjCancer","癌旁Beta不超过:", min=0, max=0.5, value=0.25,step=0.05),
				sliderInput("dltB", "dltBeta不小于:", min = 0.2, max = 1, value = 0.2,step=0.05),
				radioButtons("plotType","display plot", c("ROC" = "ROC", "HeatMap" = "HeatMap","BoxPlot"="BoxPlot"), inline=T)
		),
		column(9,dataTableOutput('DMR'),style='background-color:burlywood')
	),
	fluidRow(
    



   )




))