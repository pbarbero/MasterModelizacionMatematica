
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
# 
# http://www.rstudio.com/shiny/
#

library(shiny)



shinyServer(function(input, output) {
   
  dataset<-reactive({
    
    x<-runif(input$ssize)
    y<- .4+ 3*x+ 2*exp(-(x-.4)^2/.05) +.09*rnorm(input$ssize)
    data.frame(x=x,y=y)
  })
  
  output$ksPlot <- renderPlot({
        
    # draw the scatterplot and the kernel smoothing line  with the specified bandwidth
    plot(x = dataset()$x, y = dataset()$y, type="n")
    if(input$scatter) points(x = dataset()$x, y = dataset()$y)
    if(input$truereg) 
      {xg<-seq(0,1,length.out = 200)
       yg<-.4+3*xg+ 2*exp(-(xg-.4)^2/0.05)
      lines(yg~xg,  lwd=3,col=3)
         
    }
    lines(ksmooth(x = dataset()$x, y = dataset()$y, 
                  kernel = input$kernel,
                  bandwidth = input$bw,
                  range.x = c(0,1),
            n.points = 200),col=2,lwd=.5)
    
  })
  
})
