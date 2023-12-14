library(tidyverse)
library(tseries)
library(car)
setwd("C:/Users/esthe/Downloads/Acads/FYP/MyCode")
setwd("/Users/estherchoo/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/FYP/MyCode")
sapply(paste0("./functions/", list.files("./functions/")), source)

#data starts y2000 w9 20 feb

###DATA CLEAN##
weather <- read.csv("../SingaporeHourly.csv")

#multicollinearity
newweather <- weather[seq(from=1,by=6,to=186264),]
cor(newweather[2:5])

#25.76	0.1031	30.05	88.01
#8329 8160

data <- read.csv("../sgdengue.csv")
data <- data[,-c(171:339, 507:675, 843:1011)]

#y and z.low -> low freq
y <- as.vector(data[,3])

#check stationary with augmented dicker-fuller test
#p value is smaller thn 0.01, series is stationary
adf.test(y) #DF val -5.021, lag order 10, pval<0.01 

scaleddata <- scale(data)
y <- as.vector(scaleddata[,3])

z.low <- data.frame(rep(NA, length(y)))
for(i in 1:12){
  z.low <- cbind(z.low, dplyr::lag(y,i))
  colnames(z.low)[i+1] <- paste0("ylag", i)
}
z.low <- as.matrix(z.low[13:1059,])
z.low <- z.low[,-1]
y <- as.matrix(y[13:1059])

#high freq vars
z.high.raw <- as.matrix(scaleddata[13:1059, 4:ncol(scaleddata)])
ah <- z.high.raw[,1:167]
rh <- z.high.raw[,168:334]
tp <- z.high.raw[,335:501]
at <- z.high.raw[,502:668]
z.high <- list(ah, rh, tp, at)

save(list=c("z.low","z.high", "y"), file="C:/Users/ASUS/Downloads/Acads/FYP/MyCode/data.Rdata")

#all aggregated (normal regression)
wea <- weather[8161:186264,]
wea <- wea[-1]
wea <- scale(wea)
z.agg <- z.low
for(i in 1:4){
  wk <- data.frame((colMeans(matrix(wea[,i], nrow=168))))
  colnames(wk) <- paste0(colnames(wea)[i], "_lag", 1)
  for(j in 1:5){
    wk <- data.frame(cbind(wk, dplyr::lag(wk[1],j)))
    colnames(wk)[j+1] <- paste0(colnames(wea)[i], "_lag", j+1)
  }
  wk <- wk[-c(1:12),]
  wk <- wk[1:1047,]
  z.agg <- cbind(z.agg, wk)
}
z.agg <- as.matrix(z.agg)

save(z.agg, file="./agg.Rdata")



