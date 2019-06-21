library(shiny)
library(data.table)
library(pROC)
library(sqldf)
library(DT)
library(tidyr)
library(ggplot2)


# Define server logic for slider examples

dataDir <- "F:/projects/CRCA/output/Report/"
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
	fpath <- sprintf('F:/projects/allData/TCGA/%s_450k.txt',cancerCode);
	print(fpath)
	mat <- fread(sprintf('F:/projects/allData/TCGA/%s_450k.txt',cancerCode),sep="\t",header=T,stringsAsFactors=F)
	mat <- as.data.frame(mat)
	list('DMR'=dmrs,'sampleInf'=sampleInf,'mat'=mat)
}
shinyServer(function(input, output,session) {
   dat <- reactiveValues();
   dataSet <- reactive({
	    cn   <- as.character(input$cancerName);
		 print("loading all the data");
		 obj  <- NULL
		 ifelse("结直肠癌"     ==cn, obj <- loadData('CRCA'),
       ifelse("宫颈内膜癌"   ==cn, obj <- loadData('UCEC'),
		 ifelse("胃癌"         ==cn, obj <- loadData('STAD')
		 ifelse("宫颈癌"       ==cn, obj <- loadData('CESA'),
		 ifelse("食管癌"       ==cn, obj <- loadData('ESCA'),
		 obj <- loadData('CRCA'))))));
		 dat$DMR <- obj$DMR
		 dat$sampleInf <- obj$sampleInf
		 dat$mat <- obj$mat
       #return(obj)
	})
   observe({
      dataSet();
		updateSelectInput(session, "site", choices = names(dat$DMR))
	})
	##########################
	output$sampleInf <- renderTable({
       #obj <- dataSet();
		 dat$sampleInf
	})
	output$DMR1 <- renderDataTable({
      site <- as.character(input$site)
		#obj <- dataSet();
		obj <- dat$DMR[[site]]
		print(sprintf('%s - dim=%d %d ',site,dim(obj)[1],dim(obj)[2]))
		dat$dataTable <- obj[obj$betaN < as.numeric(input$adjCancer) & obj$dltBeta >= as.numeric(input$dltB),]
      
		dat$dataTable[,-c(1,12)]  %>% DT::datatable(options=list(scrollY = '320px',lengthMenu = c(10, 20, 50,100), pageLength = 10),selection='single') %>% 
		formatRound(columns=c('betaN', 'betaC', 'dltBeta', 'senesitivity', 'specificity', 'AUC'), digits=4)
      
	},server=T)

	output$areaPlot <- renderPlot(
      {
		 selectedRow <- input$DMR1_rows_selected;
		  if(!is.null(selectedRow)) {
				indx <- as.integer(selectedRow[1])
				site <- as.character(input$site)
				selected <- dat$dataTable[indx,];
				probes   <- strsplit(selected[1,12],";")[[1]]
				mat <- dat$mat[dat$mat$probeId %in% probes,];
            rownames(mat) <- mat$probeId
				mat <- mat[probes,]
				clnames <- colnames(mat)
				normalSamples <- colnames(mat)[grepl('-11A-',clnames)]
				cancerSamples <- colnames(mat)[grepl('-01A-',clnames)]
				mat <- mat[,c(normalSamples,cancerSamples)]
				geneSymbol <- selected[1,1]
            L <- dim(mat)[2]
				mat1 <- c()
            for(v in rownames(mat)) {
					 mat1 <- rbind(mat1,data.frame(
					 'X'=1:L,
					 'Beta'=as.numeric(mat[v,]),
					 'probeId'=rep(v,L),
					  'sampleType'=c(rep('normal',length(normalSamples)),rep('cancer',length(cancerSamples)))
					 ))
				}
				#mat1$probeId <- factor(1:length(probes),levels =probes) 
				ggplot(mat1,aes(x=X,y=Beta,fill=sampleType))+geom_area()+facet_grid(probeId ~ .)+
				scale_x_discrete(expand = c(0.01,0))+
				theme(legend.position="none",
				    panel.grid.major = element_blank(),
	             panel.grid.minor = element_blank(),
	             panel.border = element_blank(),
	            axis.line.x = element_line(color="black", size = 0.25),
	            axis.line.y = element_line(color="black", size = 0.25),
	            plot.title   = element_text(size=16),
	          panel.background = element_blank())
				

				
		   
		}
      },height=400
	)
   
	output$rocPlot <- renderPlot({
      selectedRow <- input$DMR1_rows_selected;
		  if(!is.null(selectedRow)) {
				indx <- as.integer(selectedRow[1])
				site <- as.character(input$site)
				selected <- dat$dataTable[indx,];
				probes   <- strsplit(selected[1,12],";")[[1]]
				geneSymbol <- selected[1,1]
				print(geneSymbol)
				mat <- dat$mat[dat$mat$probeId %in% probes,];
            rownames(mat) <- mat$probeId
				mat <- mat[probes,-1]
				
				Bta <- apply(mat,2,function(v){mean(v,na.rm=T)})
				lbls <- as.character(sapply(colnames(mat),function(v){a <- strsplit(v,"-")[[1]][4];ifelse(a=='01A','cancer','normal')}))
				objROC <- roc('controls'=Bta[lbls=='normal'],'cases'=Bta[lbls=='cancer'])
				plot(objROC,print.auc=TRUE,print.thres="best",col='blue',legacy.axes = TRUE,main=geneSymbol,print.auc.cex=2.0,print.thres.cex=1.5)
		}
   }
	)
})