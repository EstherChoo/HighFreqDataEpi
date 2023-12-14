library(tidyverse)
library(forecast)
library(EnvStats)
library(scoringRules)

folder <- "C:/Users/esthe/OneDrive - Nanyang Technological University/FYP/MyCode"
folder <- "/Users/estherchoo/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/FYP/MyCode"
setwd(folder)

load("./data.Rdata")
load("./out/4b_outptpred.Rdata")
load("./out/4b_outdenspred.Rdata")

#1. mape 2. rmae 3. dm test 4. pit+ks 5. logscore 6. crps
#1-3: ar, 4-6: bal, 7-9: bgl, 10-12:bl, 13-16: nmidas

model.names <- c("ar1", "ar2", "ar3", "bal1", "bal2", "bal3", "bgl1", "bgl2", "bgl3", "bl1",
                "bl2", "bl3", "nmidasar", "nmidasbal", "nmidasbl", "nmidasagg")

unscaled.y <- as.matrix(read.csv("../sgdengue.csv")[3])

descale <- function(y, val){
  center <- mean(y)
  scale <- sd(y)
  return(val*scale + center)
}

og.y <- unscaled.y[13:1059,]
testy <- y[524:nrow(y),]

calc.err <- function(ypred, yactual){ #func to calculate residuals
  h <- length(ypred)
  res <- list()
  mape <- list()
  for(i in 1:h){
      realyactual <- descale(og.y, yactual[(1+i):length(yactual)])
      realypred <- descale(og.y, ypred[[i]][1:(length(ypred[[i]])-i)])
      res[[i]] <- abs(ypred[[i]][1:(length(ypred[[i]])-i)] - yactual[(1+i):length(yactual)])
      mape[[i]] <- mean(abs((realyactual-realypred)/realyactual))
      mae[[i]] <- mean(abs(realyactual-realypred))
  }
  return(list(res, mape, mae))
}


#residuals and mean abs error
i <- 1
error <- list()
mape <- list()
mae <- list()
for(model in out.ptpred){
  result <- calc.err(model, testy)
  error[[i]] <- result[[1]]
  mape[[i]] <- result[[2]]
  mae[[i]] <- result[[3]]
  
  i <- i+1
}
mapedf <- sapply(mape, c)
colnames(mapedf) <- model.names
write.csv(mapedf, "./out/tables/mape.csv")

#relative mean abs error
rmae <- list()
for(h in 1:20){
  mat <- matrix(0, nrow=16, ncol=16) #total 16 models
  for(i in 1:length(out.ptpred)){
    for(j in i:length(out.ptpred)){
      mat[i,j] <- mae[[i]][[h]]/mae[[j]][[h]] #val < 1 => model i (row) better than model j (col)
      mat[j,i] <- mae[[j]][[h]]/mae[[i]][[h]]
    }
  }
  rmae[[h]] <- mat
}

newdf <- data.frame(t(model.names))

for(i in 1:length(rmae)){
  newdf <- rbind(newdf, rep(i,16))
  colnames(rmae[[i]]) <- names(newdf)
  newdf <- rbind(newdf, rmae[[i]])
}

colnames(newdf) <- model.names
newdf <- newdf[-1,]
write.csv(newdf, "./out/tables/rmae.csv")


#diebold mariano test
dm <- list()
dmp <- list()
for(h in 1:20){
  dmmat <- matrix(0, nrow=16, ncol=16)
  dmmatp <- matrix(0, nrow=16, ncol=16)
  for(i in 1:length(out.ptpred)){
    for(j in 1:length(out.ptpred)){
      if(i==j){
        next
      }
      dmmat[i,j] <- dm.test(error[[i]][[h]][1:504], error[[j]][[h]][1:504], h=1, alternative="g")$statistic #is model in col better than model in row
      dmmatp[i,j] <- dm.test(error[[i]][[h]][1:504], error[[j]][[h]][1:504], h=1, alternative="g")$p.value
    }
  }
  dm[[h]] <- dmmat
  dmp[[h]] <- dmmatp
}

#dm test statistic
dmdf <- data.frame(t(model.names))

for(i in 1:length(dm)){
  dmdf <- rbind(dmdf, rep(i,16))
  colnames(dm[[i]]) <- names(dmdf)
  dmdf <- rbind(dmdf, dm[[i]])
}

colnames(dmdf) <- model.names
dmdf <- dmdf[-1,]
write.csv(dmdf, "./out/tables/dieboldmariano.csv")

#dm test p value
dmpdf <- data.frame(t(model.names))

for(i in 1:length(dmp)){
  dmpdf <- rbind(dmpdf, rep(i,16))
  colnames(dmp[[i]]) <- names(dmpdf)
  dmpdf <- rbind(dmpdf, dmp[[i]])
}

colnames(dmpdf) <- model.names
dmpdf <- dmpdf[-1,]
write.csv(dmpdf, "./out/tables/dieboldmarianopvalues.csv")

  
##DENSITY PERFORMANCE
#pit
pit <- function(mat, actualy, h, unif){ #takes in matrix of timepts x density - each col is a density for 1 timepoint
  pitval <- c()
  for(n in 1:ncol(mat)){
    pitval <- c(pitval, pemp(actualy[n+h-1], mat[,n]))
  }
  return(ks.test(pitval, unif))
}

pit.ks <- matrix(0, nrow=20, ncol=16)
pit.ks.p <- matrix(0, nrow=20, ncol=16)
i <- 1
unif <- runif(1000)
for(file in out.denspred){
  for(h in 1:20){
    ks <- pit(file[[h]], testy, h=1, unif) #get pit values
    pit.ks[h,i] <- ks$statistic #check pit values against uniform using ks
    pit.ks.p[h,i] <- ks$p.value
  }
  i <- i + 1
}

colnames(pit.ks) <- model.names
colnames(pit.ks.p) <- model.names

#add signif to pit.ks.p
symp <- symnum(pit.ks.p, corr = FALSE,
               cutpoints = c(0,  .001,.01,.05, .1, 1),
               symbols = c("***","**","*over","."," "))

for(row in 1:20){
  for(col in 1:16){
    pit.ks.p[row, col] <- paste0("$", pit.ks.p[row,col], "^{", symp[row, (col-1)], "}$")
  }
}

write.csv(pit.ks, "./out/tables/pit.csv")
write.csv(pit.ks.p, "./out/tables/pitpvalues.csv")

#log score and crps
ls.mat <- matrix(0, nrow=20, ncol=16)
crps.mat <- matrix(0, nrow=20, ncol=16)
i <- 1
for(i in 1:16){
  file <- out.denspred[[i]]
  for(h in 1:20){
    mat <- file[[h]]
    mat <- mat[-c((nrow(mat)-h+1):nrow(mat)),]
    yactual <- testy[(h+1):length(testy)]
    ls <- logs_sample(yactual, mat)
    ls <- mean(ls[which(is.finite(ls))])
    crps <- mean(crps_sample(yactual, mat))
    ls.mat[h,i] <- ls
    crps.mat[h,i] <- crps
  }
}

colnames(ls.mat) <- model.names
colnames(crps.mat) <- model.names
write.csv(ls.mat, "./out/tables/logscore.csv")
write.csv(crps.mat, "./out/tables/crps.csv")

save(list=c("mapedf", "rmae", "dmdf", "dmpdf", "ls.mat", "crps.mat", "pit.ks", "pit.ks.p"), file="./out/outsamplemeasures.Rdata")

# rmsfe <- list()
# sum.ar1sfe <- rep(0,20)
# sum.ar2sfe <- rep(0,20)
# sum.ar3sfe <- rep(0,20)
# 
# ###MAE###
# for(arfile in arfiles){
#   ind <- parse_number(arfile)
#   ar <- load(file=paste0(folder, arfile))
#   assign("ar", get(ar))
#   rm(list=ls(pattern="(arpred)"))
#   
#   ar1pred <- sapply(ar[[1]], "[[", 1)
#   ar1sfe <- fe(ar1pred, y[(ind+1):(ind+20)])
#   
#   ar2pred <- sapply(ar[[2]], "[[", 1)
#   ar2sfe <- fe(ar2pred, y[(ind+1):(ind+20)])
#   
#   ar3pred <- sapply(ar[[3]], "[[", 1)
#   ar3sfe <- fe(ar3pred, y[(ind+1):(ind+20)])
#   
#   sum.ar1sfe <- sum.ar1sfe + ar1sfe
#   sum.ar2sfe <- sum.ar2sfe + ar2sfe
#   sum.ar3sfe <- sum.ar3sfe + ar3sfe
# }
# rmsfe[[1]] <- sum.ar1sfe/length(arfiles)
# rmsfe[[2]] <- sum.ar2sfe/length(arfiles)
# rmsfe[[3]] <- sum.ar3sfe/length(arfiles)
# 
# ##############
# blfiles <- grep("bl", files, value = TRUE)
# 
# rmsfe <- list()
# sum.ar1sfe <- rep(0,20)
# sum.ar2sfe <- rep(0,20)
# sum.ar3sfe <- rep(0,20)
# 
# ###MSFE###
# for (blfile in arfiles){
#   ind <- parse_number(blfile)
#   ar <- load(file=paste0(folder, arfile))
#   assign("ar", get(ar))
#   rm(list=ls(pattern="(arpred)"))
#   
#   ar1pred <- sapply(ar[[1]], "[[", 1)
#   ar1sfe <- fe(ar1pred, y[(ind+1):(ind+20)])
#   
#   ar2pred <- sapply(ar[[2]], "[[", 1)
#   ar2sfe <- fe(ar2pred, y[(ind+1):(ind+20)])
#   
#   ar3pred <- sapply(ar[[3]], "[[", 1)
#   ar3sfe <- fe(ar3pred, y[(ind+1):(ind+20)])
#   
#   sum.ar1sfe <- sum.ar1sfe + ar1sfe
#   sum.ar2sfe <- sum.ar2sfe + ar2sfe
#   sum.ar3sfe <- sum.ar3sfe + ar3sfe
# }
# rmsfe[[1]] <- sum.ar1sfe/length(arfiles)
# rmsfe[[2]] <- sum.ar2sfe/length(arfiles)
# rmsfe[[3]] <- sum.ar3sfe/length(arfiles)
# 
# fe <- function(ypred, yactual){
#   if(length(ypred)!=length(yactual)){
#     cat("Pred and Actual must be same lengths")
#     return(NULL)
#   }
#   for(i in 1:length(ypred)){
#     ypred[[i]] <- abs(ypred[[i]] - yactual[[i]])
#   }
#   return(unlist(ypred))
# }
