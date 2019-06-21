shinyUI(fluidPage(
   fluidRow(  
	   column(3,
		      selectInput("cancerCode", "癌症名称:", 
                choices = c(
						 '结直肠癌'   ='CRCA',
						 '胃 癌'      ='STAD',
						 '子宫内膜癌' ='UCEC',
						 '宫颈癌'     ='CECA',
						 '食管癌'     ='ESCA'

				)),
				selectInput("site", "选择一个病灶位置:",choices="overall"
				 ),

				sliderInput("adjCancer","癌旁Beta不超过:", min=0, max=1, value=0.2,step=0.1),
				sliderInput("dltB", "dltBeta不小于:", min = 0.2, max = 1, value = 0.1,step=0.1),
				fluidRow(
				   column(7,radioButtons("plotType","display plot:", c("ROC" = "ROC", "HeatMap" = "HeatMap"), inline=T)),
				   column(5,downloadButton("downloadLeft", "download it"),inline=T)
				),
				fluidRow(
				   column(7,radioButtons("plotType1",NULL, c("AreaPlot"="areaPlot","BoxPlot"="boxPlot"), inline=T)),
				   column(5,downloadButton("downloadRight", "download it"),style="float:left")
				)
		),
		column(9,tabsetPanel(
        tabPanel("DMR",dataTableOutput('DMR')),
		  tabPanel("Sample Information", tableOutput('sampleInf'))))
	),
	fluidRow(
	   column(4,plotOutput('leftPlot')),
		column(8,plotOutput('rightPlot'))
   ),
	hr(),
    HTML("<footer style='width:100%;margin-left:auto;margin-right:auto;text-align:center'>
            <span>Copyright<span  style='font-family:Microsoft yahei'> &copy; </span>2019-2022.Wuhan Ammunition 
				Life Technology Lt.,Company All rights reserved.</span>
          </footer>"
	)

))