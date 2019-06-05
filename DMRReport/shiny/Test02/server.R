library(shiny)
library(data.table)
library(pROC)
library(sqldf)
library(

# Define server logic for slider examples

loadData <- function(){
   dataDir <- 'F:/projects/UCEC/output/Report'
   fnames <- list.list('F:/projects/


}

shinyServer(function(input, output) {

  # Reactive expression to compose a data frame containing all of the values
  sliderValues <- reactive({
    # Compose data frame
  data.frame(
      Name = c("Integer", 
               "Decimal",
               "Range",
               "Custom Format",
               "Animation"),
      Value = as.character(c(input$integer, 
                             input$decimal,
                             paste(input$range, collapse=' '),
                             input$format,
                             input$animation)), 
      stringsAsFactors=FALSE)
  }) 

  # Show the values using an HTML table
  output$values <- renderTable({
    sliderValues()
  })
})