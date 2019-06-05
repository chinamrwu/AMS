shinyUI(pageWithSidebar(

  #  Application title
  headerPanel("宫颈内膜癌差异甲基化(DMR)报告"),

  # Sidebar with sliders that demonstrate various available options
  sidebarPanel(
    # Simple integer interval
    sliderInput("adjCancer", "癌旁Beta不超过:", 
                min=0, max=1, value=0.1,step=0.05),

    # Decimal interval with step value
    sliderInput("Cancer", "癌症Beta不小于:", 
                min = 0, max = 1, value = 0.3, step= 0.05),

    # Specification of range within an interval
    sliderInput("dltB", "dltBeta不小于:",
                min = 0.1, max = 1, step=0.05)
   ),

  # Show a table summarizing the values entered
  mainPanel(
    tableOutput("values")
  )
))
#######################################################################