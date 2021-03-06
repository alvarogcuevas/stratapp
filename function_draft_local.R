library(partitions)
library(combinat)
library(ggplot2)
library(plotly)
library(akima)
library(parallel)

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

start_time <- Sys.time()
prueba <- Perm_fn(15)
lista <- apply(X = prueba, MARGIN = 1, FUN = RMSE_fn)
lista2 <- do.call(rbind.data.frame, lista)
df <- cbind(prueba, lista2)
finish_time <- Sys.time() - start_time
# Without parallel:
# With N=12: Time difference of 15.1287 secs
# With N=20: Time difference of 4.791182 mins
# With parallel:
# N=12: Time difference of 1.949439 secs
# N=20: Time difference of 21.85425 secs

ll <- ggplot(data = df, aes(x = bias, y = sd,n000=n000,n100=n100,n010=n010,n110=n110,n001=n001,n101=n101,n011=n011,n111=n111)) + 
  geom_point(aes(colour=rmse))+ 
  scale_colour_gradient(low = "green",high="red")

ggplotly(ll)


mm<- ll+geom_density2d()
mm

values <- rep(5, 8)
df <- as.data.frame(c(values, RMSE_fn(values)))
colnames(df) <-
  c("n000",
    "n100",
    "n010",
    "n110",
    "n001",
    "n101",
    "n011",
    "n111",
    "bias",
    "sd",
    "rmse")

x_axis <- seq(from = df$bias - 0.2,
              to = df$bias + 0.2,
              by = 0.01)
y_axis <- seq(from = df$sd - 0.2,
              to = df$sd + 0.2,
              by = 0.01)
df_background <- expand.grid(bias = x_axis, sd = y_axis)
df_background$rmse <- sqrt(df_background$bias ^ 2 + df_background$sd ^ 2)

ll <- ggplot() +
  geom_raster(data = df_background, aes(x = bias, y = sd, fill = rmse)) +
  scale_fill_gradientn(colours = rev(heat.colors(5))) +
  geom_point(data = df, aes(x = bias, y = sd, z = rmse))
ggplotly(ll)

mm <- ggplot() +
  geom_bar(data = df, aes(rmse))
mm

# Double checking that the numbers add up when excluding the forbidden cases:

prueba<-Perm_fn(9)
#6435
#11440

prueba2<-subset(prueba,(n100==0 & n101==0))
#1287
#2002

prueba3<-subset(prueba,!(n100==0 & n101==0))
#6435-1287=5148
#11440-2002=9438

prueba4<-subset(prueba,(n111==0 & n110==0))
#1287
#2002

prueba5<-subset(prueba,(n100==0 & n101==0 & n111==0 & n110==0))
#165
#220

prueba6<-subset(prueba,!(n100==0 & n101==0) & !(n111==0 & n110==0))
#4026
#7656




