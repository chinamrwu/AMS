library(shiny)
library(data.table)
library(pROC)
library(sqldf)
library(DT)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Define server logic for slider examples
loadData <- function(cancerCode){
   print(sprintf("loading dataset  from  %s......",cancerCode))
   mat450 <- fread(sprintf('%s/%s_reportMat.txt',cancerCode,cancerCode),sep="\t",header=T,stringsAsFactors=F)
   mat450 <- as.data.frame(mat450)
    
	fnames <- list.files(cancerCode)
	fnames <- fnames[grepl('DMR',fnames)]
   tmp <- sapply(fnames,function(v){
	   a <- strsplit(v,"_")[[1]];
      b <- strsplit(a[3],"\\.")[[1]][1];
		c(a[2],b)
	})
	indx <- order(as.integer(tmp[1,]));
	siteNames <- as.character(tmp[2,])
	siteNames <- siteNames[indx]

	DMR <- list()
   for(nms in fnames){
       tmp <- read.csv(sprintf("%s/%s",cancerCode,nms),header=T,stringsAsFactors=F)
		 DMR[[length(DMR)+1]] <- tmp
	}
   names(DMR) <- siteNames

	sampleInf <- read.csv(sprintf("%s/SampleInf.csv",cancerCode),header=T,stringsAsFactors=F,sep="\t")
	list('DMR'=DMR,'mat'=mat450,'sampleInf'=sampleInf)
}
shinyServer(function(input, output,session) {
  dat <- reactiveValues();
  observe({
      obj  <- loadData(input$cancerCode)
		 dat$DMR <- obj$DMR
		 dat$mat <- obj$mat
		 dat$sampleInf <- obj$sampleInf
		 updateSelectInput(session, "site", choices = names(dat$DMR))
	})
#################################################################
output$sampleInf <-renderDataTable({
		obj <- dat$sampleInf
		rownames(obj) <- 1:dim(obj)[1]
		obj  %>% DT::datatable(options = list(scrollY = '500px',dom = 't',pageLength = dim(obj)[1],
		columnDefs = list(list(className = 'dt-center', targets = 0:dim(obj)[2]))
		),selection='none') %>%
		 formatStyle(0, target= 'row',lineHeight='60%')
	},server=T)
#############################################################
areaPlot <- function(selectedRow){
        #selectedRow <- input$DMR1_rows_selected;
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
				theme(legend.position="bottom",
				      legend.title=element_blank(),
						legend.text=element_text(size=16),
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
rocPlot <- function(selectedRow){
        #selectedRow <- input$DMR1_rows_selected;
       
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
heatmapPlot <- function(selectedRow){
   if(!is.null(selectedRow)) {
				indx <- as.integer(selectedRow[1])
				selected <- dat$dataTable[indx,];
				probes   <- strsplit(selected[1,12],";")[[1]]
				geneSymbol <- selected[1,1]
				mat <- dat$mat[,probes];
				Label <- dat$mat$label
            
				mat[is.na(mat)] <- 0;
				mat <- data.frame(t(mat));

          
				ann_col <- data.frame('sampleType'=Label,row.names = colnames(mat)) 
				ann_col$sampleType <- factor(ann_col$sampleType)
				names(ann_col) <- " "

				p1 <- pheatmap(mat, color = brewer.pal(11,"RdYlBu")[11:1], fontsize_col = 8,
				annotation_col = ann_col,fontsize=15,
				cluster_rows = F, cluster_cols = T,show_rownames = F, show_colnames = F, 
				main = geneSymbol,silent=T)
				return(p1)
	 }
}
###############################
probeBoxPlot <- function(selectedRow){
  if(!is.null(selectedRow)) {
			indx <- as.integer(selectedRow[1])
			selected <- dat$dataTable[indx,];
			probes   <- strsplit(selected[1,12],";")[[1]]
			strTitle <- paste0(selected[1,1:4],collapse="_")
			mat   <- dat$mat[,c('label',probes)];
			matX  <- c();
			L <- dim(mat)[1];
			for(probe in probes){
            matX <- rbind(matX,data.frame('Beta'=mat[,probe],'label'=mat[,'label'],'probeId'=rep(probe,L),stringsAsFactors=F))
			}

			p1 <- ggplot(matX,aes(x=probeId, y=Beta, color = label,fill=label))+geom_boxplot()+
			ggtitle(strTitle)+
			stat_boxplot(geom ='errorbar')+
			#scale_x_discrete(labels=probes)+
			theme(
			axis.text.x = element_text(size=18,color="black"),
			axis.text.y = element_text(size=18,color="black"),
			axis.line.x = element_line(color="black", size = 0.25),
			axis.line.y = element_line(color="black", size = 0.25),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			#panel.border = element_blank(),
			panel.background = element_blank(),
			plot.title = element_text(hjust = 0.5,size=18))+
			scale_colour_manual(values=c('normal'='black','cancer'='black'))
			return(p1)
	 }
	 
}
probeViolin  <- function(selectedRow){
}
#####################################################
output$DMR <- renderDataTable({
      site <- as.character(input$site)
		obj <- dat$DMR[[site]]
      obj <- obj[order(obj$dltBeta,decreasing=T),]
		rownames(obj) <- 1:dim(obj)[1]

		dat$dataTable <- obj[obj$betaN < as.numeric(input$adjCancer) & obj$dltBeta >= as.numeric(input$dltB),]
		dat$dataTable[,-12]  %>% DT::datatable(options=list(scrollY = '300px',lengthMenu = c(50,100), pageLength = 50),
		selection='single') %>% 
		formatRound(columns=c('betaN', 'betaC', 'dltBeta', 'senesitivity', 'specificity', 'AUC'), digits=4) %>%
		formatStyle(c('dltBeta','betaN'), target = 'cell',
          backgroundColor = 'yellow') %>%
		formatStyle(0, target= 'row',lineHeight='50%')      
	},server=T)
output$rightPlot <- renderPlot(
 {
   plotType    <- input$plotType1
	p1 <- NULL;
	if('areaPlot'==plotType){     p1 <- areaPlot(input$DMR_rows_selected)}
	else if('boxPlot' ==plotType){p1 <- probeBoxPlot(input$DMR_rows_selected)}
	else if('violin'  ==plotType){p1 <-  probeViolin(input$DMR_rows_selected)}
	else{ p1 <- NULL}
	p1
   }
)
output$leftPlot <- renderPlot(
  {
    plotType <- input$plotType;
	 if(plotType=='ROC'){
       rocPlot(input$DMR_rows_selected)
	 }else{ heatmapPlot(input$DMR_rows_selected)
	 }
  }
)
####################################################
 
output$downloadLeft <- downloadHandler(
    filename = function() { 
	   selectedRow <- input$DMR_rows_selected;
		#print(selectedRow)
		plotType <- input$plotType
		cancerCode <- input$cancerCode
		#print(plotType)
		fname <- NULL
		if(!is.null(selectedRow)){
			indx <- as.integer(selectedRow[1])
			selected <- dat$dataTable[indx,];
			fname <- paste0(selected[1,1:4],collapse='_')
			fname <- paste0(c(cancerCode,plotType,fname,'.svg'),collapse="_")
			#print(fname)
		}
		return(fname)
	 },
    content = function(filename) {
      svg(filename)
		if(grepl('roc',tolower(input$plotType))){rocPlot(input$DMR_rows_selected)}
		else if(grepl('heatmap',tolower(input$plotType))) { 
		 print(heatmapPlot(input$DMR_rows_selected))
		}
      dev.off()
    }
)
##############################################################
output$downloadRight <- downloadHandler(
    filename = function() { 
	   selectedRow <- input$DMR_rows_selected;
		plotType <- input$plotType1
		cancerCode <- input$cancerCode
		fname <- NULL
		if(!is.null(selectedRow)){
			indx <- as.integer(selectedRow[1])
			selected <- dat$dataTable[indx,];
			fname <- paste0(selected[1,1:4],collapse='_')
			fname <- paste0(c(cancerCode,plotType,fname,'.svg'),collapse="_")
		}
		return(fname)
	 },
    content = function(filename) {
      svg(filename)
		if(grepl('areaplot',tolower(input$plotType1))){print(areaPlot(input$DMR_rows_selected))}
		else if(grepl('boxplot',tolower(input$plotType1))) { print(probeBoxPlot(input$DMR_rows_selected))}
      dev.off()
    }
)
##################################################################
scannerClick <- 



})