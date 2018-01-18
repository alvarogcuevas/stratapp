library(shiny)
library(ggplot2)
library(plotly)
library(partitions)
library(combinat)
library(parallel)

########################################
# Shiny app for MSE analysis of strata #
########################################

# Every Shiny app has to have three main parts: 
#   -UI (how the user sees the app)
#   -server(the logic statements R needs to produce the output)
#   -call to the function shinyApp() (to launch the app)

###############################################################
# Additional functions #
########################
 
# RMSE_fn does the following:
#   -gets a vector of n_{ijk} values, 
#   -places them in a three-dimensional matrix and 
#   -calculates the bias and variance components of the rmse

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

# Perm_fn does the following:
#   -gets an integer,
#   -obtains all the partitions of said integer,
#   -in a parallel setting, it obtains all the permutations of each partition
#   -deletes cases where n_{ijk}=0 would make a problem in the denominator when calculating the RMSE

Perm_fn <- function(nnn) {
  nijk <- restrictedparts(n = nnn,
                          m = 8,
                          include.zero = T)
  
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, c("nijk", "permn"), envi = environment())
  
  n_p <- parApply(
    cl = cl,
    X = nijk,
    MARGIN = 2,
    FUN = function(x) {
      permut <- permn(x)
      permut2 <- do.call(rbind.data.frame, permut)
      colnames(permut2) <-
        c("n000",
          "n100",
          "n010",
          "n110",
          "n001",
          "n101",
          "n011",
          "n111")
      unique(permut2)
    }
  )
  
  stopCluster(cl)
  
  before_forbidden <- do.call(rbind.data.frame, n_p)
  subset(
    before_forbidden,
    !(n100 == 0 & n101 == 0) &
      !(n111 == 0 & n110 == 0) &
      !(n001 == 0 & n000 == 0) &
      !(n011 == 0 & n010 == 0)
  )
}

############################## UI ########################

ui <- fluidPage(
  
  h2("TITLE"),
  
  tabsetPanel(
    
    ##### FIRST TAB 
    # -Enter all the n_{ijk} to obtain a single dot
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
                             min=1,max=10,value = 5),
                 actionButton("go", "Go")
               ),
               mainPanel(
                 verbatimTextOutput("graph_header"),
                 plotlyOutput("dotmap")
               )
             )
    ),
    
    ##### SECOND TAB 
    # -Enter a single integer to obtain a heat map with all the permutations 
    #  of the partions of that integer
    tabPanel("N input",
             h2("N as an input"),
             sidebarLayout(
               sidebarPanel(
                 sliderInput("nn", 
                             label = "N",
                             min=8,max=20,value = 10),
                 actionButton("go2", "Go")
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
  
  #The dataset used in output$heatmap_n will be updated each time a new integer is entered
  reac_dataset<-eventReactive(input$go2,{
    Perm_fn(input$nn)
  })
  
  #This is a heatmap of rmse for the permutations of partitions of a given integer
  output$heatmap_N<-renderPlotly({
    
    #Generate the dataset
    dataset <- reac_dataset()
    
    #Compute bias and variance for each permutation
    perm_lis<-apply(X = dataset,MARGIN = 1,FUN = RMSE_fn)
    
    #Turn the list into a data.frame
    perm_lis2<-do.call(rbind.data.frame,perm_lis)
    
    #Attach each permutation to its (bias,variance)
    df<-cbind(dataset,perm_lis2)
    
    #Plot the interactive heatmap
    ll <- ggplot(data = df, aes(x = bias, y = sd,n000=n000,n100=n100,n010=n010,n110=n110,n001=n001,n101=n101,n011=n011,n111=n111)) +
      geom_point(aes(colour=rmse))+
      scale_colour_gradient(low = "green",high="red")+
      geom_density2d()
    
    ggplotly(ll)
  })
  
  #The dataset used in output$dotmap is updated each time one of the n_{ijk} is changed
  reac_nijk<-eventReactive(input$go,{
    c(input$n000,input$n100,input$n010,input$n110,input$n001,input$n101,input$n011,input$n111)
  })
  
  #This is a single point plot for a single case of n_{ijk}
  output$dotmap<-renderPlotly({
    
    #Update the dataset
    values<-reac_nijk()
    
    #Calculate the (bias,variance) for this specific case
    df<-as.data.frame(c(values,RMSE_fn(values)))
    colnames(df)<-c("n000","n100","n010","n110","n001","n101","n011","n111","bias","sd","rmse")
    
    #Generate the interactive plot
    ll<-ggplot()+
      geom_point(data=df,aes(x = bias, y = sd,z=rmse))+
      xlim(-1,1)+
      ylim(0,0.6)
    ggplotly(ll)
  })
  
  
  #For the single case, this shows the sum of n_{ijk}
  output$graph_header<-renderText({
    values<-c(input$n000,input$n100,input$n010,input$n110,input$n001,input$n101,input$n011,input$n111)
    paste("The sample size is N=",sum(values))
  })
  
  #Small detail that shows how many available cores for parallel computing
  output$no_cores<-renderText({
    paste("The number of available cores is",detectCores())
  })
  
}



############################## CALL SHINYAPP TO RENDER APP ############################################

shinyApp(ui = ui, server = server)
