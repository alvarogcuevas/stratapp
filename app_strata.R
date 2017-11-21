library(shiny)
library(ggplot2)
library(plotly)
library(partitions)
library(combinat)
library(parallel)

# The code has to have three main parts: 
#   -UI (how the user sees the app)
#   -server(the logic statements R needs to produce the output)
#   -call to the function shinyApp() (to launch the app)

RMSE_fn <- function(nvector) {
  
  nvalues <- array(nvector, c(2, 2, 2))
  
  ###### Bias: Naive effect estimator
  
  py2_t2 <- sum(nvalues[2, , 2]) / sum(nvalues[2, ,])
  py2_t1 <- sum(nvalues[1, , 2]) / sum(nvalues[1, ,])
  
  b.naive_effect <- py2_t2 - py2_t1
  
  ##### Bias: True effect
  
  py2_t2x1 <- nvalues[2, 1, 2] / sum(nvalues[2, 1,])
  py2_t1x1 <- nvalues[1, 1, 2] / sum(nvalues[1, 1,])
  px1 <- sum(nvalues[, 1,]) / sum(nvalues)
  
  py2_t2x2 <- nvalues[2, 2, 2] / sum(nvalues[2, 2,])
  py2_t1s1 <- nvalues[1, 2, 2] / sum(nvalues[1, 2,])
  
  b.true_effect <- (py2_t2x1 - py2_t1x1) * px1 + (py2_t2x2 - py2_t1s1) * (1 - px1)
  
  bias <- b.naive_effect - b.true_effect
  
  ##### Variance: Naive effect estimator
  
  v.naive_effect <-
    (py2_t2 * (1 - py2_t2)) / sum(nvalues[2, ,]) + (py2_t1 * (1 - py2_t1)) / sum(nvalues[1, ,])
  
  ##### Variance: True effect
  
  v.true_effect <-
    ((py2_t2x1 * (1 - py2_t2x1)) / sum(nvalues[2, 1,]) + (py2_t1x1 * (1 - py2_t1x1)) /
       sum(nvalues[1, 1,])) * px1 + ((py2_t2x2 * (1 - py2_t2x2)) / sum(nvalues[2, 2,]) +
                                       (py2_t1s1 * (1 - py2_t1s1)) / sum(nvalues[1, 1,])) * (1 - px1)
  
  variance <- v.naive_effect + v.true_effect
  
  ############## Both ################################
  
  data.frame(
    bias = bias,
    sd = sqrt(variance),
    rmse = sqrt(bias ^ 2 + variance)
  )
  
}

Perm_fn <- function(nnn) {
  nijk <- restrictedparts(n = nnn,
                          m = 8,
                          include.zero = FALSE)
  
  cl<-makeCluster(8)
  clusterExport(cl,c("nijk", "permn"),envi=environment())
  
  n_p <- parApply(cl = cl,X = nijk,MARGIN = 2,FUN = function(x) {
    permut <- permn(x)
    permut2 <- do.call(rbind.data.frame, permut)
    colnames(permut2) <- c("n000","n100","n010","n110","n001","n101","n011","n111")
    unique(permut2)
  }
  )
  
  stopCluster(cl)
  
  do.call(rbind.data.frame, n_p)
}

############################## UI ########################

ui <- fluidPage(
  
  h2("TITLE"),
  
  tabsetPanel(
    
    ##### FIRST TAB 
    tabPanel("n_txy input",
             h2("n_txy as an input"),
             sidebarLayout(
               sidebarPanel(
                 sliderInput("n000", 
                             label = "n_000",
                             min=1,max=10,value = 5),
                 sliderInput("n100", 
                             label = "n_100",
                             min=1,max=10,value = 5),
                 sliderInput("n010", 
                             label = "n_010",
                             min=1,max=10,value = 5),
                 sliderInput("n110", 
                             label = "n_110",
                             min=1,max=10,value = 5),
                 sliderInput("n001", 
                             label = "n_001",
                             min=1,max=10,value = 5),
                 sliderInput("n101", 
                             label = "n_101",
                             min=1,max=10,value = 5),
                 sliderInput("n011", 
                             label = "n_011",
                             min=1,max=10,value = 5),
                 sliderInput("n111", 
                             label = "n_111",
                             min=1,max=10,value = 5)
               ),
               mainPanel(
                 verbatimTextOutput("graph_header"),
                 plotlyOutput("dotmap")
               )
             )
    ),
    
    ##### SECOND TAB 
    tabPanel("N input",
             h2("N as an input"),
             sidebarLayout(
               sidebarPanel(
                 sliderInput("nn", 
                             label = "N",
                             min=1,max=20,value = 10)
               ),
               mainPanel(
                 verbatimTextOutput("no_cores"),
                 plotlyOutput("heatmap_N")
               )
             )
    )
    )
)


#################################### SERVER #############################

# server is always a function of two list-type arguments: input and output

server <- function(input, output) {
  
  output$heatmap_N<-renderPlotly({
    #Generate the dataset
    dataset <- Perm_fn(input$nn)
    perm_lis<-apply(X = dataset,MARGIN = 1,FUN = RMSE_fn)
    perm_lis2<-do.call(rbind.data.frame,perm_lis)
    df<-cbind(dataset,perm_lis2)

    ll <- ggplot(data = df, aes(x = bias, y = sd,n000=n000,n100=n100,n010=n010,n110=n110,n001=n001,n101=n101,n011=n011,n111=n111)) +
      geom_point(aes(colour=rmse))+
      scale_colour_gradient(low = "green",high="red")
    ggplotly(ll)
  })
  
  output$dotmap<-renderPlotly({
    values<-c(input$n000,input$n100,input$n010,input$n110,input$n001,input$n101,input$n011,input$n111)
    df<-as.data.frame(c(values,RMSE_fn(values)))
    colnames(df)<-c("n000","n100","n010","n110","n001","n101","n011","n111","bias","sd","rmse")
    
    ll<-ggplot()+
      geom_point(data=df,aes(x = bias, y = sd,z=rmse))+
      xlim(-1,1)+
      ylim(0,0.6)
    ggplotly(ll)
  })
  
  output$graph_header<-renderText({
    values<-c(input$n000,input$n100,input$n010,input$n110,input$n001,input$n101,input$n011,input$n111)
    paste("The sample size is N=",sum(values))
  })
  
  output$no_cores<-renderText({
    paste("The number of available cores is",detectCores())
  })
  
}



############################## CALL SHINYAPP TO RENDER APP ############################################

shinyApp(ui = ui, server = server)
