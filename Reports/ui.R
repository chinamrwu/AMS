shinyUI(fluidPage(
   tags$head(
      #tags$style(type="text/css", "label{ display: table-cell; text-align: center; vertical-align: middle; } .form-group { display: table-row;}")
    ),
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
		  tabPanel("Sample Information", tableOutput('sampleInf')),
		  tabPanel("DMR Scanner",
		     fluidRow(
             column(3,selectInput('chr','Choromosome:',
						 choices=c(
                      'chr1' = 'chr1','chr2' = 'chr2','chr3' = 'chr3','chr4' = 'chr4','chr5' = 'chr5',
                      'chr6' = 'chr6','chr7' = 'chr7','chr8' = 'chr8','chr9' = 'chr9','chr10' = 'chr10',
							 'chr11' = 'chr11','chr12' = 'chr12','chr13' = 'chr13','chr14' = 'chr14','chr15' = 'chr15',
							 'chr16' = 'chr16','chr17' = 'chr17','chr18' = 'chr18','chr19' = 'chr19','chr20' = 'chr20',
							 'chr21' = 'chr21'))),
				 column(3,numericInput('chrStart','Start',1,min=1)),
				 column(3,numericInput('chrEnd','End',1,min=1)),
				 column(3,actionButton('scanButton','Goooo'))
		     )##DMRSCanner fluidRow
		  )
		  
		  
		  )
		  )
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