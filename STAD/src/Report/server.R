library(shiny)
library(data.table)
library(pROC)
library(sqldf)
library(DT)
library(tidyr)
library(ggplot2)


# Define server logic for slider examples


loadData <- function(){
   print("loading dataset  ......")
   mat450 <- fread('F:/projects/STAD/output/matDMR_20190614.txt',sep="\t",header=T,stringsAsFactors=F)
   mat450 <- as.data.frame(mat450)
   DMR <- read.csv('F:/projects/STAD/output/STAD_DMR_20190612.csv',header=T,stringsAsFactors=F)
	list('DMR'=DMR,'mat'=mat450)
}
shinyServer(function(input, output,session) {
   dat <- reactiveValues();

   dataSet <- reactive({
		 obj  <- loadData()
		 dat$DMR <- obj$DMR
		 dat$mat <- obj$mat
	})
   observe({
      dataSet();
		minV <- min(dat$DMR$betaN);
		maxV <- max(dat$DMR$betaN)
		updateSliderInput(session, 'adjCancer', min =minV, max =maxV,value=(maxV+minV)/2)
		
		minV <- min(dat$DMR$dltBeta);
		maxV <- max(dat$DMR$dltBeta)
		updateSliderInput(session, 'dltB', min =minV, max =maxV,value=(maxV+minV)/2)
	})
	##########################
output$DMR1 <- renderDataTable({
		obj <- dat$DMR
		dat$dataTable <- obj[obj$betaN < as.numeric(input$adjCancer) & obj$dltBeta >= as.numeric(input$dltB),]
		dat$dataTable[,-12]  %>% DT::datatable(options=list(scrollY = '250px',lengthMenu = c(10, 20, 50,100), pageLength = 10),
		selection='single') %>% 
		formatRound(columns=c('betaN', 'betaC', 'dltBeta', 'senesitivity', 'specificity', 'AUC'), digits=4) %>%
		formatStyle(c('dltBeta','betaN'), target = 'cell',
          backgroundColor = 'yellow') %>%
		formatStyle(0, target= 'row',lineHeight='50%')
      
	},server=T)
   ##################################################
areaPlot <- function( ){
        selectedRow <- input$DMR1_rows_selected;
		  if(!is.null(selectedRow)) {
				indx <- as.integer(selectedRow[1])
				selected <- dat$dataTable[indx,];
				probes   <- strsplit(selected[1,12],";")[[1]]
				mat <- dat$mat[,c('label',probes)];
            
				normalIndex <- which(mat$label=='normal')
				cancerIndex <- which(mat$label=='cancer')
				index <- c(normalIndex,cancerIndex)
				mat <- mat[c(normalIndex,cancerIndex),]
				geneSymbol <- selected[1,1]
				L <- dim(mat)[1]
				L1 <- length(probes)
				mat1 <- c()
            for(probe in probes) {
					 mat1 <- rbind(mat1,data.frame(
					 'X'=1:L,
					 'Beta'=as.numeric(mat[,probe]),
					 'probeId'=rep(probe,L),
					  'sampleType'=dat$mat$label[index])
				 
				)}
				#mat1$probeId <- factor(1:length(probes),levels =probes) 
				p1 <- ggplot(mat1,aes(x=X,y=Beta,fill=sampleType))+geom_area()+facet_grid(probeId ~ .)+
            ggtitle(sprintf("%s:%s",geneSymbol,paste0(selected[1,2:5],collapse="_")))+
				scale_x_discrete(expand = c(0.01,0))+
				theme(legend.position="none",
				    panel.grid.major = element_blank(),
	             panel.grid.minor = element_blank(),
	             panel.border = element_blank(),
	            axis.line.x = element_line(color="black", size = 0.25),
	            axis.line.y = element_line(color="black", size = 0.25),
	            plot.title   = element_text(size=16),
	          panel.background = element_blank())
				 return(p1)
		}
}
   ###################################################
rocPlot <- function(){
        selectedRow <- input$DMR1_rows_selected;
       
		  if(!is.null(selectedRow)) {
				indx <- as.integer(selectedRow[1])
				selected <- dat$dataTable[indx,];
				probes   <- strsplit(selected[1,12],";")[[1]]
				geneSymbol <- selected[1,1]
				mat <- dat$mat[,probes];
				
				Bta  <- as.numeric(apply(mat,1,function(v){mean(v,na.rm=T)}))
				lbls <- dat$mat$label
				objROC <- roc('controls'=Bta[lbls=='normal'],'cases'=Bta[lbls=='cancer'])
				plot(objROC,print.auc=TRUE,print.thres="best",col='blue',legacy.axes = TRUE,main=geneSymbol,print.auc.cex=2.0,print.thres.cex=1.5)
		}
	}
	#####################################################
output$areaPlot <- renderPlot({
areaPlot()},height=400)
output$rocPlot <- renderPlot(rocPlot())

#######################################
output$rocButton <- downloadHandler(
    filename = function() { 
	   selectedRow <- input$DMR1_rows_selected;
		fname <- NULL
		if(!is.null(selectedRow)){
			indx <- as.integer(selectedRow[1])
			selected <- dat$dataTable[indx,];
			fname <- paste0(selected[1,1:4],collapse='_')
			fname <- paste0('ROC_',fname,'.svg')
		}
		return(fname)
	 },
    content = function(file) {
      svg(file)
		rocPlot()
      dev.off()
    }
)
  #############################################
output$areaButton <- downloadHandler(
    filename = function() { 
	   selectedRow <- input$DMR1_rows_selected;
		fname <- NULL
		if(!is.null(selectedRow)){
			indx <- as.integer(selectedRow[1])
			selected <- dat$dataTable[indx,];
			fname <- paste0(selected[1,1:4],collapse='_')
			fname <- paste0('AreaPlot_',fname,'.svg')
		}
		return(fname)
	 },
    content = function(fname) {
      svg(fname,width=15,height=3)
		print(areaPlot())
      dev.off()
    }
  )
  #############################################################
  output$DMRButton <- downloadHandler(
    filename = function() { 
       paste0('DMR_',Sys.Date(),'_.csv')
	 },
    content = function(fname) {
	   obj <- dat$dataTable[as.integer(input$DMR1_rows_all),];
      write.csv(obj,fname)
    }
  )



  ####################################################
 


})