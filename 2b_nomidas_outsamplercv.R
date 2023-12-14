library(tidyverse)
library(doParallel)
library(foreach)
library(invgamma)
library(statmod)
library(MASS)
library(mvtnorm)

setwd("C:/Users/esthe/Downloads/Acads/FYP/MyCode")
sapply(paste0("./functions/", list.files("./functions")), source)
load("data.Rdata")
load("agg.Rdata")

load("2108rollingcv.Rdata")

#set up
cl <- makeCluster(3)
registerDoParallel(cl)
cat("starting code")
cat("/n")

nmidasrollingcv <- function(row, y, z.low){
  train.nmidasbl <- list()
  train.nmidasbal <- list()
  train.nmidasar <- list()
  
  y.train <- as.matrix(y[1:row-1,1])
  y.test <- as.matrix(y[row,1])
  z.low.train <- z.low[1:row-1,1:12]
  z.low.test <- t(as.matrix(z.low[row,1:12]))
  
  for(h in 1:12){
    train.nmidasbl[[h]] <- nomidas.bl(h, y.train, z.low.train, iterations=5000)
    train.nmidasbal[[h]] <- nomidas.bal(h, y.train, z.low.train, iterations=5000)
    train.nmidasar[[h]] <- nomidas.ar(h, y.train, z.low.train, iterations=5000)
  }
  cat("timepoint:", row, "progress: 50%")
  cat("\n")
  
  summary.nmidasbl <- summarisenmidas(train.nmidasbl, 1000)
  summary.nmidasbal <- summarisenmidas(train.nmidasbal, 1000)
  summary.nmidasar <- summarisenmidas(train.nmidasar, 1000)
  cat("timepoint:", row, "progress: 60%")
  cat("\n")
  
  pred.nmidasbl <- predictnmidas(summary.nmidasbl, z.low.test)
  pred.nmidasbal <- predictnmidas(summary.nmidasbal, z.low.test)
  pred.nmidasar <- predictnmidas(summary.nmidasar, z.low.test)
  cat("timepoint:", row, "progress: 70%")
  cat("\n")
  
  dens.nmidasbl <- predictdensnmidas(summary.nmidasbl, z.low.test)  
  dens.nmidasbal <- predictdensnmidas(summary.nmidasbal, z.low.test)
  dens.nmidasar <- predictdensnmidas(summary.nmidasar, z.low.test)
  cat("timepoint:", row, "progress: 95%")
  cat("\n")
  
  nmidaspred <- list()
  for(m in 1:3){
    nmidaspred[[m]] <- list()
  }
  for(h in 1:12){
    nmidaspred[[1]][[h]] <- list(pred.nmidasar[[h]], dens.nmidasar[[h]])
    nmidaspred[[2]][[h]] <- list(pred.nmidasbal[[h]], dens.nmidasbal[[h]])
    nmidaspred[[3]][[h]] <- list(pred.nmidasbl[[h]], dens.nmidasbl[[h]])
  }
  
  save(nmidaspred, file=paste0("./rollingcv/", "nmidas", row, ".Rdata"))
}

#split the data
split <- round(1/2*nrow(y))
foreach(row=split:nrow(y), .packages=c('invgamma', 'statmod', 'MASS', 'mvtnorm', 'tidyverse')) %dopar% {
  nmidasrollingcv(row=row, y=y, z.low=z.low)
}

agg.nmidasrollingcv <- function(row, y, z.agg){
  train.nmidasagg <- list()
  
  y.train <- as.matrix(y[1:row-1,1])
  y.test <- as.matrix(y[row,1])
  z.low.train <- z.agg[1:row-1,]
  z.low.test <- t(as.matrix(z.agg[row,]))
  
  for(h in 1:12){
    train.nmidasagg[[h]] <- nomidas.ar(h, y.train, z.low.train, iterations=5000)
  }
  cat("timepoint:", row, "progress: 50%")
  cat("\n")
  
  summary.nmidasagg <- summarisenmidas(train.nmidasagg, 1000)
  cat("timepoint:", row, "progress: 60%")
  cat("\n")
  
  pred.nmidasagg <- predictnmidas(summary.nmidasagg, z.low.test)
  cat("timepoint:", row, "progress: 70%")
  cat("\n")

  dens.nmidasagg <- predictdensnmidas(summary.nmidasagg, z.low.test)
  cat("timepoint:", row, "progress: 95%")
  cat("\n")
  
  nmidasagg <- list()
  for(h in 1:12){
    nmidasagg[[h]] <- list(pred.nmidasagg[[h]], dens.nmidasagg[[h]])
  }
  
  save(nmidasagg, file=paste0("./rollingcv/", "nmidasagg", row, ".Rdata"))
}

#split the data
split <- round(1/2*nrow(y))
foreach(row=split:nrow(y), .packages=c('invgamma', 'statmod', 'MASS', 'mvtnorm', 'tidyverse')) %dopar% {
  agg.nmidasrollingcv(row=row, y=y, z.agg=z.agg)
}

stopCluster(cl)



  
  
















