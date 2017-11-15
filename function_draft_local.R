library(partitions)
library(combinat)
library(ggplot2)
library(plotly)
library(akima)


Perm_fn <- function(nnn) {
  nijk <- restrictedparts(n = nnn,
                          m = 8,
                          include.zero = FALSE)
  nijk_diff <- apply(X = nijk, MARGIN = 2, FUN = permn)
  nijk_mat <- lapply(nijk_diff, FUN = rbind.data.frame)
  n_p <- apply(
    X = nijk,
    MARGIN = 2,
    FUN = function(x) {
      permut <- permn(x)
      permut2 <- do.call(rbind.data.frame, permut)
      colnames(permut2) <- c("n000","n100","n010","n110","n001","n101","n011","n111")
      unique(permut2)
    }
  )
  do.call(rbind.data.frame, n_p)
}

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

start_time<-Sys.time()
prueba<-Perm_fn(20)
lista<-apply(X = prueba,MARGIN = 1,FUN = RMSE_fn)
lista2<-do.call(rbind.data.frame, lista)
df<-cbind(prueba,lista2)
finish_time<-Sys.time()-start_time
# With N=12: Time difference of 15.1287 secs
# With N=20: Time difference of 4.791182 mins

ll <- ggplot(data = df, aes(x = bias, y = sd,n000=n000,n100=n100,n010=n010,n110=n110,n001=n001,n101=n101,n011=n011,n111=n111)) + 
  geom_point(aes(colour=rmse))+ 
  scale_colour_gradient(low = "green",high="red")

ggplotly(ll)



values<-rep(5,8)
df<-as.data.frame(c(values,RMSE_fn(values)))
colnames(df)<-c("n000","n100","n010","n110","n001","n101","n011","n111","bias","sd","rmse")

x_axis<-seq(from=df$bias-0.2,to=df$bias+0.2,by=0.01)
y_axis<-seq(from=df$sd-0.2,to=df$sd+0.2,by=0.01)
df_background<-expand.grid(bias=x_axis,sd=y_axis)
df_background$rmse<-sqrt(df_background$bias^2+df_background$sd^2)

ll<-ggplot()+
  geom_raster(data=df_background,aes(x=bias,y=sd,fill=rmse))+
  scale_fill_gradientn(colours = rev(heat.colors(5)))+
geom_point(data=df,aes(x = bias, y = sd,z=rmse))
ggplotly(ll)

mm<-ggplot()+
  geom_bar(data = df,aes(rmse))
mm



