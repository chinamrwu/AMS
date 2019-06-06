library(shiny)
library(data.table)
library(pROC)
library(sqldf)


# Define server logic for slider examples

dataDir <- "F:/projects/UCEC/output/Report/"
loadData <- function(cancerCode){
   print(sprintf("loading dataset for %s ......",cancerCode))
   fnames <- list.files(dataDir)
	dmrFiles <- fnames[grepl('DMR_',fnames)]
   dmrs <- list()
	indx <- as.integer(sapply(dmrFiles,function(v){strsplit(v,"_")[[1]][2]}))
	dmrFiles <- dmrFiles[indx]
	dmrNames <- as.character(sapply(dmrFiles,function(v){a <- strsplit(v,"_")[[1]][3]; strsplit(a,"\\.")[[1]][1]}))
   dmrs <- sapply(dmrFiles,simplify=F,function(v){
      tmp <- read.csv(paste0(dataDir,v),stringsAsFactors=F,header=T);
		tmp
	})
   names(dmrs) <- dmrNames
	sampleInf <- read.csv(paste0(dataDir,'sampleInf.csv'),header=T,stringsAsFactors=F)
	list('DMR'=dmrs,'sampleInf'=sampleInf)
}

shinyServer(function(input, output,session) {
  
   dataSet <- reactive({
	    cn   <- as.character(input$cancerName);
		 obj  <- NULL
	    ifelse("结直肠癌" ==cn,   obj <- loadData('CRCA') ,
		 ifelse("宫颈癌"   ==cn,   obj <- loadData('CESA'),
		 ifelse("食管癌"   ==cn,   obj <- loadData('ESCA'),obj <- loadData('UCEC'))));
       return(obj)
	})
   observe({
      obj1   <- dataSet();
		updateSelectInput(session, "site", choices = names(obj1$DMR))
	})
	##########################
	output$sampleInf <- renderTable({
       obj <- dataSet();
		 obj$sampleInf
	})
	output$DMR <- renderDataTable({
      site <- as.character(input$site)
		obj <- dataSet();
		obj <- obj$DMR[[site]]
		obj
	})

})