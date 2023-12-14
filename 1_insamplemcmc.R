library(tidyverse)
library(invgamma)
library(statmod)
library(MASS)
library(mvtnorm)

folder <- "C:/Users/esthe/OneDrive - Nanyang Technological University/FYP/MyCode"
#folder <- "/Users/estherchoo/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/FYP/MyCode"
setwd(folder)
sapply(paste0("./functions/", list.files("./functions")), source)
load(file="data.Rdata")
load(file="agg.Rdata")

#RUN
out.mcmcar1 <- list()
out.mcmcar2 <- list()
out.mcmcar3 <- list()

out.mcmcbl1 <- list()
out.mcmcbl2 <- list()
out.mcmcbl3 <- list()

out.mcmcbal1 <- list()
out.mcmcbal2 <- list()
out.mcmcbal3 <- list()

out.mcmcbgl1 <- list()
out.mcmcbgl2 <- list()
out.mcmcbgl3 <- list()

out.nmidasar <- list()
out.nmidasbl <- list()
out.nmidasbal <- list()
out.nmidasagg <- list()


for(h in 1:12){
  out.mcmcar1[[h]] <- mcmc.ar(h, y, z.low, z.high, iterations=5000, shape=1)
  out.mcmcar2[[h]] <- mcmc.ar(h, y, z.low, z.high, iterations=5000, shape=2)
  out.mcmcar3[[h]] <- mcmc.ar(h, y, z.low, z.high, iterations=5000, shape=3)
}
 

for(h in 1:12){
  print(h)
  out.mcmcbl1[[h]] <- mcmc.bl(h, y, z.low, z.high, iterations=5000, shape=1)
  out.mcmcbl2[[h]] <- mcmc.bl(h, y, z.low, z.high, iterations=5000, shape=2)
  out.mcmcbl3[[h]] <- mcmc.bl(h, y, z.low, z.high, iterations=5000, shape=3)
}


for(h in 1:12){
  print(h)
  out.mcmcbal1[[h]] <- mcmc.bal(h, y, z.low, z.high, iterations=5000, shape=1)
  out.mcmcbal2[[h]] <- mcmc.bal(h, y, z.low, z.high, iterations=5000, shape=2)
  out.mcmcbal3[[h]] <- mcmc.bal(h, y, z.low, z.high, iterations=5000, shape=3)
}

for(h in 1:12){
  print(h)
  out.mcmcbgl1[[h]] <- mcmc.bgl(h, y, z.low, z.high, iterations=5000, shape=1)
  out.mcmcbgl2[[h]] <- mcmc.bgl(h, y, z.low, z.high, iterations=5000, shape=2)
  out.mcmcbgl3[[h]] <- mcmc.bgl(h, y, z.low, z.high, iterations=5000, shape=3)
}

for(h in 1:12){
  out.nmidasar[[h]] <- nomidas.ar(h, y, z.low, 5000)
  out.nmidasbl[[h]] <- nomidas.bl(h, y, z.low, 5000)
  out.nmidasbal[[h]] <- nomidas.bal(h, y, z.low, 5000)
  out.nmidasagg[[h]] <- nomidas.ar(h, y, z.agg, 5000)
}


save(list=c("out.mcmcar1", "out.mcmcar2", "out.mcmcar3", "out.mcmcbl1", "out.mcmcbl2", "out.mcmcbl3", 
            "out.mcmcbal1", "out.mcmcbal2", "out.mcmcbal3", "out.mcmcbgl1", "out.mcmcbgl2", "out.mcmcbgl3",
            "out.nmidasar", "out.nmidasbl", "out.nmidasbal", "out.nmidasagg"),
     file="./out/1_insampleresults.Rdata")





    









