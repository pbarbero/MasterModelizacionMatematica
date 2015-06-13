
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)

shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Suavizado Tipo Núcleo"),
  
  # Sidebar with a slider input for number of bins
  sidebarPanel(
    numericInput(label = "sample size",inputId = "ssize", 
                 value = 150,min = 10, max=500),
    sliderInput("bw",
                "Select the bandwidth:",
                min = 0,
                max = 1,
                value = 0.5),
  
    selectInput("kernel", label="Select Kernel", 
                choices=c("Uniforme"="box", "Normal"="normal")),
    checkboxInput("scatter",label = "Show a scatterplot", value=TRUE),
    checkboxInput("truereg", label="Show the true regression function", value=FALSE)
  ),

  # Show a plot of the generated distribution
  mainPanel(
        plotOutput("ksPlot")
  )
))
