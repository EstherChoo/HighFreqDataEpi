#combine all separate timept files
#[ar/bal/bgl/bl][timept].Rdata
#"ar1", "ar2", "ar3", "bal1", "bal2", "bal3", "bgl1", "bgl2", "bgl3", "bl1", "bl2", "bl3", "nmidasar", "nmidasbal", "nmidasbl", "nmidasagg"
#obj 1: pt pred - list of 15 model, each model list of 20 horizons, each horizon is 524x1
#obj 2: dens pred - list of 15 model, each model list of 20 horizons, each horizon is 524x4000

setwd("C:/Users/esthe/Downloads/Acads/FYP/MyCode")
setwd("C:/Users/esthe/OneDrive - Nanyang Technological University/FYP/MyCode")

#loading files
files <- list.files("./rollingcv")
start <- 524
end <- 1047

#set up structure of output data
out.ptpred <- list()
out.denspred <- list()
for(i in 1:16){
  out.ptpred[[i]] <- list()
  out.denspred[[i]] <- list()
}
for(i in 1:16){
  for(j in 1:12){
    out.ptpred[[i]][[j]] <- matrix(nrow=524, ncol=1)
    out.denspred[[i]][[j]] <- matrix(nrow=524, ncol=4000)
  }
}
names(out.ptpred) <- c("ar1", "ar2", "ar3", "bal1", "bal2", "bal3", "bgl1", "bgl2", "bgl3", "bl1",
                       "bl2", "bl3", "nmidasar", "nmidasbal", "nmidasbl", "nmidasagg")
names(out.denspred) <- names(out.ptpred)

#TO CHANGE (ar-0, bal-3, bgl-6, bl-9, nmidas-12, nmidasagg-15)
model <- "bl" #change 
ind <- 9 #change acc to model index!!!

for(i in start:end){
  n <- i-523
  load(file=paste0("./rollingcv/", model, i, ".Rdata"))
  assign("lst", get(paste0(model, "pred", i)))
  for(mod in 1:3){
    for(hor in 1:12){
      out.ptpred[[(mod+ind)]][[hor]][n, 1] <- lst[[mod]][[hor]][[1]]
      out.denspred[[(mod+ind)]][[hor]][n,] <- lst[[mod]][[hor]][[2]]
    }
  }
  rm(list=paste0(model, "pred", i))
}

model <- "nmidas"
ind <- 12

for(i in start:end){
  n <- i-523
  load(file=paste0("./rollingcv/", model, i, ".Rdata"))
  assign("lst", get(paste0(model, "pred")))
  for(mod in 1:3){
    for(hor in 1:12){
      out.ptpred[[(mod+ind)]][[hor]][n, 1] <- lst[[mod]][[hor]][[1]]
      out.denspred[[(mod+ind)]][[hor]][n,] <- lst[[mod]][[hor]][[2]]
    }
  }
  rm(list=paste0(model, "pred"))
}

for(i in 16:16){
  out.ptpred[[i]] <- list()
  out.denspred[[i]] <- list()
}

for(i in 16:16){
  for(j in 1:20){
    out.ptpred[[i]][[j]] <- matrix(nrow=524, ncol=1)
    out.denspred[[i]][[j]] <- matrix(nrow=524, ncol=4000)
  }
}

model <- "nmidasagg"
ind <- 15

mod <- 1
for(i in start:end){
  n <- i-523
  load(file=paste0("./rollingcv/", model, i, ".Rdata"))
  assign("lst", get(paste0(model)))
  for(hor in 1:12){
    out.ptpred[[(mod+ind)]][[hor]][n, 1] <- lst[[hor]][[1]]
    out.denspred[[(mod+ind)]][[hor]][n,] <- lst[[hor]][[2]]
  }
}

save(out.ptpred, file="./out/4b_outptpred.Rdata")
save(out.denspred, file="./out/4b_outdenspred.Rdata")

