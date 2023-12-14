library(ggplot2)
library(cowplot)
library(GGally)
library(lubridate)
library(ggprism)
library(tidyverse)

folder <- "C:/Users/esthe/OneDrive - Nanyang Technological University/FYP/MyCode"
folder <- "/Users/estherchoo/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/FYP/MyCode"
setwd(folder)

#best model based on aic,dic,waic is BAL3

#trace plot
load("./out/4a_traces.Rdata")

traceplotter <- function(model, h, modelname){
  plotlist <- list()
  if(ncol(model[[h]])==26){
    colnames(model[[h]])[1] <- "intercept"
  }
  
  for(coef in 1:ncol(model[[h]])){
    plotlist[[coef]] <- local({
      coef <- coef
      ggplot() +
      geom_line(aes(x=1:nrow(model[[h]]), y=model[[h]][,coef])) +
      xlab("Iterations") +
      ylab(colnames(model[[h]])[coef]) +
      theme_bw() + 
      theme(axis.line = element_line(color='black'),
              plot.background = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())
    })
    title <- ggdraw() + 
      draw_label(
        paste0("Traceplot of ", modelname),
        fontface = 'bold',
        x = 0,
        hjust = 0
      ) +
      theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 7)
      )
  }
  return(plot_grid(title, plotlist=plotlist, ncol=3))
}

traceplotter1 <- function(model, h, modelname){ #use for nmidasagg
  plotlist1 <- list()
  plotlist2 <- list()
  
  colnames(model[[h]]) <- str_replace_all(colnames(model[[h]]), "t2m_C", "at")
  colnames(model[[h]]) <- str_replace_all(colnames(model[[h]]), "tp_mm", "tp")
  colnames(model[[h]]) <- str_replace_all(colnames(model[[h]]), "ah_gm3", "ah")
  colnames(model[[h]]) <- str_replace_all(colnames(model[[h]]), "rh_p", "rh")
  
  if(ncol(model[[h]])==26){
    colnames(model[[h]])[1] <- "intercept"
  }
  
  for(coef in 1:26){
    plotlist1[[coef]] <- local({
      coef <- coef
      ggplot() +
        geom_line(aes(x=1:nrow(model[[h]]), y=model[[h]][,coef])) +
        xlab("Iterations") +
        ylab(colnames(model[[h]])[coef]) +
        theme_bw() + 
        theme(axis.line = element_line(color='black'),
              plot.background = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())
    })
    title <- ggdraw() + 
      draw_label(
        paste0("Traceplot of ", modelname),
        fontface = 'bold',
        x = 0,
        hjust = 0
      ) +
      theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 7)
      )
  }
  
  for(coef in 27:ncol(model[[h]])){
    plotlist2[[coef-26]] <- local({
      coef <- coef
      ggplot() +
        geom_line(aes(x=1:nrow(model[[h]]), y=model[[h]][,coef])) +
        xlab("Iterations") +
        ylab(colnames(model[[h]])[coef]) +
        theme_bw() + 
        theme(axis.line = element_line(color='black'),
              plot.background = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())
    })
    title <- ggdraw() + 
      draw_label(
        paste0("Traceplot of ", modelname),
        fontface = 'bold',
        x = 0,
        hjust = 0
      ) +
      theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 7)
      )
  }
  
  return(list(plot_grid(title, plotlist=plotlist1, ncol=3), plot_grid(title, plotlist=plotlist2, ncol=3, nrow=9)))
}

tp_ar1 <- traceplotter(ar1trace, 1, "AR-(1)")
tp_ar2 <- traceplotter(ar2trace, 1, "AR-(2)")
tp_ar3 <- traceplotter(ar3trace, 1, "AR-(3)")

tp_bl1 <- traceplotter(bl1trace, 1, "BL-(1)")
tp_bl2 <- traceplotter(bl2trace, 1, "BL-(2)")
tp_bl3 <- traceplotter(bl3trace, 1, "BL-(3)")

tp_bal1 <- traceplotter(bal1trace, 1, "BAL-(1)")
tp_bal2 <- traceplotter(bal2trace, 1, "BAL-(2)")
tp_bal3 <- traceplotter(bal3trace, 1, "BAL-(3)")

tp_bgl1 <- traceplotter(bgl1trace, 1, "BGL-(1)")
tp_bgl2 <- traceplotter(bgl2trace, 1, "BGL-(2)")
tp_bgl3 <- traceplotter(bgl3trace, 1, "BGL-(3)")

tp_nmagg <- traceplotter1(nmagg.trace, 1, "NMIDAS-AGG")
tp_nmar <- traceplotter(nmar.trace, 1, "NMIDAS-AR")
tp_nmbal <- traceplotter(nmbaltrace, 1, "NMIDAS-BAL")
tp_nmbl <- traceplotter(nmbl.trace, 1, "NMIDAS-BL")

pdf(file="./plots/tp_nmagg1.pdf", width=8.3 , height=11.7)
tp_nmagg[[1]]
dev.off() 

pdf(file="./plots/tp_nmagg2.pdf", width=8.3 , height=11.7)
tp_nmagg[[2]]
dev.off()

for(file in ls(pattern="(tp)")){
  if(file=="tp_nmagg"){next}
  pdf(file=paste0("./plots/", file,".pdf"), width=8.3 , height=11.7)
  print(get(file))
  dev.off() 
}


#fitted vs actual values
load("./out/3_insamplepoint.Rdata")
load("./unscaleddata.Rdata")

#descale y and predictions
unscaled.y <- as.matrix(read.csv("../sgdengue.csv")[3])

descale <- function(y, val){
  center <- mean(y)
  scale <- sd(y)
  return(val*scale + center)
}

og.y <- unscaled.y[13:1059,]

date <- seq(from=as.Date("2000-02-27"), length.out=1047, by="weeks")
yearbreaks <- date[match(unique(year(date)), year(date))]
monthbreaks <- as.Date(c())
for(year in unique(year(date))){
  subsetdates <- date[year(date)==year]
  monthbreaks <- c(monthbreaks, subsetdates[match(unique(month(subsetdates)), month(subsetdates))])
}


plot1.fitted1 <- ggplot() +
  geom_point(aes(x=date, y=og.y, colour="Actual")) + 
  geom_line(aes(x=date[2:1047], y=descale(unscaled.y, in.pred.bal2[[1]][1:1046]), colour="Fitted")) +
  scale_colour_manual("", breaks=c("Actual", "Fitted"), values=c("black", "firebrick")) +
  ylab("Case Counts") +
  xlab("Year") +
  labs(title="A. Fitted vs Actual Case Counts (1 Week Ahead)") +
  theme(legend.position = "top") +
  theme_bw() + 
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title=element_text(face="bold"),
        axis.ticks.length = unit(5, "pt")) +
  theme(axis.text.x=element_text(hjust=-0.7)) +
  scale_x_date(guide="prism_minor", minor_breaks=monthbreaks, breaks=yearbreaks,labels=sprintf("%02d", seq(0, 20, 1)) )

plot1.fitted2 <- ggplot() +
  geom_point(aes(x=date, y=og.y, colour="Actual")) + 
  geom_line(aes(x=date[6:1047], y=descale(unscaled.y,in.pred.bal2[[5]][1:1042]), colour="Fitted")) +
  scale_colour_manual("", breaks=c("Actual", "Fitted"), values=c("black", "firebrick")) +
  ylab("Case Counts") +
  xlab("Year") +
  labs(title="B. Fitted vs Actual Case Counts (5 Weeks Ahead)") +
  theme(legend.position = "top") +
  theme_bw() + 
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title=element_text(face="bold"),
        axis.ticks.length = unit(5, "pt")) +
  theme(axis.text.x=element_text(hjust=-0.7)) +
  scale_x_date(guide="prism_minor", breaks=yearbreaks, minor_breaks=monthbreaks, labels=sprintf("%02d", seq(0, 20, 1)))

plot1.fitted3 <- ggplot() +
  geom_point(aes(x=date, y=og.y, colour="Actual")) + 
  geom_line(aes(x=date[11:1047], y=descale(unscaled.y,in.pred.bal2[[10]][1:1037]), colour="Fitted")) +
  scale_colour_manual("", breaks=c("Actual", "Fitted"), values=c("black", "firebrick")) +
  ylab("Case Counts") +
  xlab("Year") +
  labs(title="C. Fitted vs Actual Case Counts (10 Weeks Ahead)") +
  theme(legend.position = "top") +
  theme_bw() + 
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title=element_text(face="bold"),
        axis.ticks.length = unit(5, "pt")) +
  theme(axis.text.x=element_text(hjust=-0.7)) +
  scale_x_date(guide="prism_minor", breaks=yearbreaks, minor_breaks=monthbreaks, labels=sprintf("%02d", seq(0, 20, 1)))

pdf(file="./plots/fig2-FittedvsActual.pdf", width=8.27, height=11.69)
plot_grid(plot1.fitted1, plot1.fitted2, plot1.fitted3, ncol=1)
dev.off()

#case counts
casets <- ggplot() +
  geom_line(aes(x=date, y=unscaled.y[13:1059,])) +
  theme_bw() + 
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks.length = unit(5, "pt"),
        axis.text.x=element_text(hjust=-0.4)) +
  labs(x="Year", y="Dengue Case Counts") +
  scale_x_date(guide="prism_minor", breaks=yearbreaks, minor_breaks=monthbreaks, labels=sprintf("%02d", seq(0, 20, 1)))

pdf(file="./plots/fig1-cases.pdf", width=7.29, height=3.69)
casets
dev.off()


#shape of weights - best model
load("./out/3_mcmcsummary.Rdata")
load("./unscaleddata.Rdata")
full_clim <- read.csv("./out/tables/descaledclimcoef.csv")
sapply(paste0("./functions/", list.files("./functions/")), source)
the <- list()
the.low <- list()
the.high <- list()
for(i in 1:4){
  the[[i]] <- sapply(in.summary.bal3[[1]][["theta"]][[i]], "[[", "mean")
}
for(i in 1:4){
  the.low[[i]] <- sapply(lapply(in.summary.bal3[[1]][["theta"]][[i]], "[[", "quantile"), "[[", "2.5%")
}
for(i in 1:4){
  the.high[[i]] <- sapply(lapply(in.summary.bal3[[1]][["theta"]][[i]], "[[", "quantile"), "[[", "97.5%")
}
weight <- list()
weight.high <- list()
weight.low <- list()

for(i in 1:4){
  weight[[i]] <- midasWeights(the, 3)[[i]] * full_clim$estimate[20+i]
  weight.low[[i]] <- midasWeights(the.low, 3)[[i]] * full_clim$estimate[20+i]
  weight.high[[i]] <- midasWeights(the.high, 3)[[i]] * full_clim$estimate[20+i]
}


x <- seq(from=6, by=6, length.out=167)
plot1.w1 <- ggplot()+
  geom_line(aes(x=x, y=weight[[1]], colour="Posterior Mean", linetype="Posterior Mean")) +
  geom_line(aes(x=x, y=weight.low[[1]], colour="2.5% LB", linetype="2.5% LB")) +
  geom_line(aes(x=x, y=weight.high[[1]], colour="97.5% UB", linetype="97.5% UB")) +
  scale_color_manual(values=c("blue", "red", "black")) +
  scale_linetype_manual(values=c(2, 2, 1)) +
  labs(x="Hours", y=expression(Absolute ~ Humidity ~ (gm^-3)), colour=NULL, linetype=NULL) +
  theme_bw() + 
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  theme(legend.position="none") +
  #theme(legend.position = c(0.75,0.75)) +
  scale_x_continuous(sec.axis=sec_axis(trans=~./24, name="Days"))
  
plot1.w2 <- ggplot() +
  geom_line(aes(x=x, y=weight[[2]], colour="Posterior Mean", linetype="Posterior Mean")) +
  geom_line(aes(x=x, y=weight.low[[2]], colour="2.5% LB", linetype="2.5% LB")) +
  geom_line(aes(x=x, y=weight.high[[2]], colour="97.5% UB", linetype="97.5% UB")) +
  scale_color_manual(values=c("blue", "red", "black")) +
  scale_linetype_manual(values=c(2, 2, 1)) +
  labs(x="Hours", y="Relative Humidty (%)", colour=NULL, linetype=NULL) +
  theme_bw() + 
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  theme(legend.position="none") +
  #theme(legend.position = c(0.75,0.75))+
  scale_x_continuous(sec.axis=sec_axis(trans=~./24, name="Days"))

plot1.w3 <- ggplot()+
  geom_line(aes(x=x, y=weight[[3]], colour="Posterior Mean", linetype="Posterior Mean")) +
  geom_line(aes(x=x, y=weight.low[[3]], colour="2.5% LB", linetype="2.5% LB")) +
  geom_line(aes(x=x, y=weight.high[[3]], colour="97.5% UB", linetype="97.5% UB")) +
  scale_color_manual(values=c("blue", "red", "black")) +
  scale_linetype_manual(values=c(2, 2, 1)) +
  labs(x="Hours", y="Total Precipitation (mm)", colour=NULL, linetype=NULL) +
  theme_bw() + 
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  theme(legend.position="none") +
  #theme(legend.position = c(0.75,0.75)) +
  scale_x_continuous(sec.axis=sec_axis(trans=~./24, name="Days"))

plot1.w4 <- ggplot()+
  geom_line(aes(x=x, y=weight[[4]], colour="Posterior Mean", linetype="Posterior Mean")) +
  geom_line(aes(x=x, y=weight.low[[4]], colour="2.5% LB", linetype="2.5% LB")) +
  geom_line(aes(x=x, y=weight.high[[4]], colour="97.5% UB", linetype="97.5% UB")) +
  scale_color_manual(values=c("blue", "red", "black")) +
  scale_linetype_manual(values=c(2, 2, 1)) +
  labs(x="Hours", y="Temperature (\u00B0C)", colour=NULL, linetype=NULL) +
  theme_bw() + 
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  #theme(legend.position = c(0.75,0.75)) +
  scale_x_continuous(sec.axis=sec_axis(trans=~./24, name="Days"))

title <- ggdraw() + 
  draw_label(
    "Contribution of Climate to 1-Week Ahead Dengue Case Counts over Time",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
weightlegend <- get_legend(plot1.w4)
weightplot <- plot_grid(plot1.w1, plot1.w2, plot1.w3, plot1.w4 +
                        theme(legend.position="none"), labels="AUTO")
weightplot2 <- plot_grid(title, weightplot, ncol=1, rel_heights = c(0.1, 1))

pdf(file="./plots/MIDASWeights.pdf", width=9.3, height=7.63)
plot_grid(weightplot2, legend=weightlegend, rel_widths = c(3, .5))
dev.off()

#which day contribute the most to dengue case counts
which.max(midasWeights(the, 3)[[1]])*6/24
which.max(midasWeights(the, 3)[[2]])*6/24
which.max(midasWeights(the, 3)[[3]])*6/24
which.max(midasWeights(the, 3)[[4]])*6/24

#coefficients
load("./out/3_mcmcsummary.Rdata")
coef.df <- function(model){
  lst <- model[[1]][["phi"]]
  if(length(lst)==17){lst <- lst[2:17]} #remove intercept
  names(lst) <- c("Cases Lag 1", "Cases Lag 2", "Cases Lag 3", "Cases Lag 4", "Cases Lag 5",
                    "Cases Lag 6", "Cases Lag 7", "Cases Lag 8", "Cases Lag 9", "Cases Lag 10", "Cases Lag 11",
                    "Cases Lag 12", "Absolute Humidity", "Relative Humidity", "Total Precipitation", "Air Temperature")
  
  conf <- lapply(lst, "[[", "quantile")
  df <- data.frame(conf.low = sapply(conf, "[[", "2.5%"),
                 conf.high = sapply(conf, "[[", "97.5%"),
                 estimate = sapply(lst, "[[", "mean"),
                 term = names(lst))
  return(df)
}

coef.df1 <- function(model){
  lst <- model[[1]][["phi"]]
  if(length(lst)==17){names(lst)[1] <-"Intercept"}
  conf <- lapply(lst, "[[", "quantile")
  df <- data.frame(conf.low = sapply(conf, "[[", "2.5%"),
                   conf.high = sapply(conf, "[[", "97.5%"),
                   estimate = sapply(lst, "[[", "mean"),
                   term = names(lst))
  return(df)
}

modelnames <- c("AR-(1)", "AR-(2)", "AR-(3)", "BAL-(1)", "BAL-(2)", "BAL-(3)", "BGL-(1)", "BGL-(2)", "BGL-(3)",
                "BL-(1)", "BL-(2)", "BL-(3)", "NMIDAS-AGG", "NMIDAS-AR", "NMIDAS-BAL", "NMIDAS-BL")


#coef table
coef <- data.frame()
i <- 1
for(model in ls(pattern="in.summary")){
  print(model)
  df <- coef.df1(get(model))
  df <- cbind(df, "model" = modelnames[i])
  i <- i + 1
  coef <- rbind(coef, df)
}
coef$term <- str_replace_all(coef$term, "ylag", "Cases Lag ")
coef$term <- str_replace_all(coef$term, "t2m_C", "at")
coef$term <- str_replace_all(coef$term, "tp_mm", "tp")
coef$term <- str_replace_all(coef$term, "ah_gm3", "ah")
coef$term <- str_replace_all(coef$term, "rh_p", "rh")
coef$term <- str_replace_all(coef$term, "at", "Air Temperature")
coef$term <- str_replace_all(coef$term, "tp", "Total Precipitation")
coef$term <- str_replace_all(coef$term, "ah", "Absolute Humidity")
coef$term <- str_replace_all(coef$term, "rh", "Relative Humidity")
colnames(coef) <- c("2.5% Quantile", "97.5% Quantile", "Estimate", "Variable", "Model")
write.csv(coef, file="./out/tables/coef.csv")

###coef plots
full_cases <- data.frame()
full_clim <- data.frame()
i <- 1
for(model in ls(pattern="in.summary")){
  print(model)
  if(grepl("nm", model)){next}
  df <- coef.df(get(model))
  df <- cbind(df, "model" = modelnames[i])
  i <- i + 1
  #full_cases <- rbind(full_cases, df[1:12,])
  full_clim <- rbind(full_clim, df[13:16,])
  
}

#unscale MIDAS climate coefs
data <- read.csv("../sgdengue.csv")
at <- read.csv("../AT.csv")
at <- at[-1]
data <- data[,-c(171:339, 507:675, 843:1011)]
data <- cbind(data, at)

z.high.raw <- as.matrix(data[13:1059, 4:ncol(data)])
ah <- z.high.raw[,1:167]
rh <- z.high.raw[,168:334]
tp <- z.high.raw[,335:501]
at <- z.high.raw[,502:668]
z.high <- list(ah, rh, tp, at)

scaleddata <- scale(data)
z.high.raw <- as.matrix(scaleddata[13:1059, 4:ncol(scaleddata)])
ah <- z.high.raw[,1:167]
rh <- z.high.raw[,168:334]
tp <- z.high.raw[,335:501]
at <- z.high.raw[,502:668]
z.high.scaled <- list(ah, rh, tp, at)

r <- 1
for(model in ls(pattern="in.summary")){
  if(grepl("nm",model)){next}
  shape <- as.numeric(substr(model, start=nchar(model), stop=nchar(model)))
  model <- get(model)
  unscaledclim <- list()
  for(i in 1:4){
    unscaledclim[[i]] <- midaspoly(sapply(model[[1]][["theta"]][[i]], "[[", "mean"), z.high[[i]], shape)
  }
  scaledclim <- list()
  for(i in 1:4){
    scaledclim[[i]] <- midaspoly(sapply(model[[1]][["theta"]][[i]], "[[", "mean"), z.high.scaled[[i]], shape)
  }
  climateintep <- c()
  for(i in 1:4){
    climateintep[i] <- ((1 - mean(scaledclim[[i]]))/sd(scaledclim[[i]])) * sd(unscaledclim[[i]]) + mean(unscaledclim[[i]])
  }
  y <- as.vector(data[,3])
  for(i in 1:4){
    full_clim[(r+i-1),1:3] <- full_clim[(r+i-1),1:3] * sd(y)  + mean(y)
    if(i==3){full_clim[(r+i-1),1:3] <- full_clim[(r+i-1),1:3]/100}
    full_clim[(r+i-1),1:3] <- full_clim[(r+i-1),1:3] / climateintep[i]
  }
  r <- r+4
}

full_clim <- full_clim %>%
  mutate(model = factor(model, levels=unique(model))) %>%
  mutate(term = factor(term, levels=unique(term)))

colours <- c("royalblue1", "royalblue2", "royalblue3", "mediumorchid1", "mediumorchid2", "mediumorchid3", "seagreen2",
                         "seagreen3", "seagreen4", "orange1", "orange2", "orange3")

#clim coefs
strip_data_midas <- data.frame(unique(full_clim$term)) %>%
  mutate(y_position = 1:nrow(.),
         ymin = y_position - 0.5,
         ymax = y_position + 0.5,
         xmin = -5,
         xmax = 5,
         fill = c(rep(c("a", "b"), length.out=4))
  ) %>%
  pivot_longer(cols=c(xmin, xmax), names_to="min_max", values_to="x")

ymin <- c()
ymax <- c()
full <- full_clim
n <- 4
for(i in seq(1, n, 2)){
  if(i==1){
    ymin <- c(ymin, as.numeric(full$term[[i]])-0.6)
  } else {
    ymin <- c(ymin, as.numeric(full$term[[i]])-0.5)
  }
  
  if(i==n){
    ymax <- c(ymax, as.numeric(full$term[[i]])+0.6)
  } else {
    ymax <- c(ymax, as.numeric(full$term[[i]])+0.5)
  }
  
}
full_clim$line <- as.character(as.integer(full_clim$model == "BAL-(3)"))
full_clim$model <- as.factor(full_clim$model)
linemap <- c("AR-(1)" = "dashed" , "AR-(2)" = "dashed" ,"AR-(3)" = "dashed" , "BAL-(1)" = "dashed" ,"BAL-(2)" = "dashed" ,
              "BAL-(3)" = "solid" , "BGL-(1) = dashed" , "BGL-(2)" = "dashed" ,
              "BGL-(3)" = "dashed" , "BL-(1)" = "dashed" ,  "BL-(2)" = "dashed" , "BL-(3)" = "dashed")
linemap <- c(rep("dashed",5), "solid", rep("dashed",6))

stripe_midas <- data.frame(ymin=ymin, ymax=ymax, xmin=-Inf, xmax=Inf)

plot1.coef2 <- ggplot(full_clim) +
  aes(y=term, x=estimate, color=model, linetype=model, group=model) + 
  geom_rect(xmin=stripe_midas$xmin[1], ymin=stripe_midas$ymin[1], xmax=stripe_midas$xmax[1], ymax=stripe_midas$ymax[1], fill="grey90", inherit.aes=F)+
  geom_rect(xmin=stripe_midas$xmin[2], ymin=stripe_midas$ymin[2], xmax=stripe_midas$xmax[2], ymax=stripe_midas$ymax[2], fill="grey90", inherit.aes=F)+
  geom_point(aes(fill="Posterior Mean"), position = position_dodge(width=.75)) + 
  geom_errorbarh(aes(xmin=conf.low, xmax=conf.high),position=position_dodge(width=.75), height=0) + 
  labs(x="Estimate", y="", colour="Models", title="Estimated Change in 1-Week Ahead Cases", fill=NULL, linetype="Models") + 
  geom_vline(xintercept=0, linetype="dashed") +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_color_manual(name="Models", values=colours) +
  scale_linetype_manual(name="Models", values=c(rep("dashed",5), "solid", rep("dashed",6))) +
  scale_fill_manual(values=NA) #+
  #guides(linetype=guide_legend(override.aes=list(linetype=c(rep("dashed",5), "solid", rep("dashed",6)))))

pdf(file="./plots/fig2-climcoef.pdf", width=5.93, height=5.07)
plot1.coef2
dev.off()



#underfit/overfit
load("./out/3_insamplepoint.Rdata")
load("./unscaleddata.Rdata")

plot1.fit1 <- ggplot() +
  geom_point(aes(x=og.y[2:1047],y=descale(og.y, in.pred.bal2[[1]][1:1046]), colour="red"), show.legend=F) +
  geom_abline()+
  xlab("Actual Case Counts") +
  ylab("Fitted Case Counts") +
  labs(title="A. 1 Weeks Ahead") +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title=element_text(face="bold")) +
  scale_x_continuous(breaks=seq(0,900,100)) +
  scale_y_continuous(breaks=seq(0,900,100)) +
  annotate("text", x=100, y=800, label="Overpredict", size=4) +
  annotate("text", x=750, y=100, label="Underpredict", size=4)

plot1.fit2 <- ggplot() +
  geom_point(aes(x=og.y[6:1047],y=descale(og.y, in.pred.bal2[[5]][1:1042]), colour="red"), show.legend=F) +
  geom_abline()+
  xlab("Actual Case Counts") +
  ylab("Fitted Case Counts") +
  labs(title="B. 5 Weeks Ahead") +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title=element_text(face="bold")) +
  scale_x_continuous(breaks=seq(0,900,100)) +
  scale_y_continuous(breaks=seq(0,900,100), limits=c(0,900)) +
  annotate("text", x=100, y=800, label="Overpredict", size=4) +
  annotate("text", x=750, y=100, label="Underpredict", size=4)

plot1.fit3 <- ggplot() +
  geom_point(aes(x=og.y[11:1047],y=descale(og.y, in.pred.bal2[[10]][1:1037]), colour="red"), show.legend=F) +
  geom_abline()+
  xlab("Actual Case Counts") +
  ylab("Fitted Case Counts") +
  labs(title="C. 10 Weeks Ahead") +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title=element_text(face="bold")) +
  scale_x_continuous(breaks=seq(0,900,100)) +
  scale_y_continuous(breaks=seq(0,900,100),limits=c(0,900)) +
  annotate("text", x=100, y=800, label="Overpredict", size=4) +
  annotate("text", x=750, y=100, label="Underpredict", size=4)
  
  
pdf(file="./plots/underoverpredict.pdf", width=12.46, height=4.47)
plot_grid(plot1.fit1, plot1.fit2, plot1.fit3, nrow=1)
dev.off()

###out of sample plots
load(file="./out/4c_outsamplemeasures.Rdata")

#MAPE
mapedf <- mapedf[1:12,]
mae1 <- gather(data.frame(mapedf), key=Model, value=MAPE, ar1:nmidasagg)
mae1$MAPE <- as.numeric(mae1$MAPE) *100
mae1 <- cbind(mae1, c(rep("Midas", 144), rep("No Midas", 48)))
mae1 <- cbind(V1=c(1:12), mae1)
colnames(mae1)[4] <- "type"

colours <- c("royalblue1", "royalblue2", "royalblue3", "mediumorchid1", "mediumorchid2", "mediumorchid3", "seagreen2",
                         "seagreen3", "seagreen4", "orange1", "orange2", "orange3", "skyblue", "brown1", "brown2", "brown3")

plot2.mae <- ggplot(data=mae1, aes(x=V1, y=MAPE, group=Model, colour=factor(Model))) +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none") +
  geom_line(linewidth=1) +
  geom_point(aes(shape=type), size=3) +
  scale_colour_manual(values=colours) +
  scale_shape_manual(values=c(1, 4)) +
  labs(colour="Model", shape="") +
  xlab("Forecast Horizon (Weeks)") +
  ylab("Mean Absolute Percentage Error (%)") +
  scale_x_continuous(guide="prism_minor", breaks=c(1,4,8,12), minor_breaks=c(1:12), labels=c(1,4,8,12))


#PIT
pit.ks <- pit.ks[1:12,]
pit.ks.p <- pit.ks.p[1:12,]
pit1 <- gather(data.frame(pit.ks), key=Model, value=Value, ar1:nmidasagg)
pit1 <- cbind(pit1, c(rep("Midas", 144),rep("No Midas", 48)))
pit1 <- cbind(c(1:12), pit1)
colnames(pit1)[4] <- "type"
colnames(pit1)[1] <- "V1"

plot2.pit <- ggplot(data=pit1, aes(x=V1, y=Value, group=Model, colour=factor(Model))) +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none") +
  geom_line(linewidth=1) +
  geom_point(aes(shape=type), size=3) +
  scale_colour_manual(values=colours) +
  scale_shape_manual(values=c(1, 4)) +
  labs(colour="Model", shape="") +
  xlab("Forecast Horizon (Weeks)") +
  ylab("Kolmogorov-Smirnov Test Statistic") +
  scale_x_continuous(guide="prism_minor", breaks=c(1,4,8,12), minor_breaks=c(1:12), labels=c(1,4,8,12))


#Logscore
ls.mat <- ls.mat[1:12,]
ls1 <- gather(data.frame(ls.mat), key=Model, value=`Logarithmic Score`, ar1:nmidasagg)
ls1 <- cbind(ls1, c(rep("Midas", 144), rep("No Midas", 48)))
ls1 <- cbind(c(1:12), ls1)
colnames(ls1)[4] <- "type"
colnames(ls1)[1] <- "V1"

plot2.ls <- ggplot(data=ls1, aes(x=V1, y=`Logarithmic Score`, group=Model, colour=factor(Model))) +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none") +
  geom_line(linewidth=1) +
  geom_point(aes(shape=type), size=3) +
  scale_colour_manual(values=colours) +
  scale_shape_manual(values=c(1, 4)) +
  labs(colour="Model", shape="") +
  xlab("Forecast Horizon (Weeks)") +
  ylab("Logarithmic Score") +
  scale_x_continuous(guide="prism_minor", breaks=c(1,4,8,12), minor_breaks=c(1:12), labels=c(1,4,8,12))


#CRPS
crps.mat <- crps.mat[1:12,]
crps1 <- gather(data.frame(crps.mat), key=Model, value=`Continuous Ranked Probability Score`, ar1:nmidasagg)
crps1 <- cbind(crps1, c(rep("Midas", 144), rep("No Midas", 48)))
crps1 <- cbind(c(1:12), crps1)
colnames(crps1)[4] <- "type"
colnames(crps1)[1] <- "V1"
crpsmodels <- c("AR-(1)", "AR-(2)", "AR-(3)", "BAL-(1)", "BAL-(2)", "BAL-(3)", "BGL-(1)", "BGL-(2)", "BGL-(3)",
                "BL-(1)", "BL-(2)", "BL-(3)", "NMIDAS-AR", "NMIDAS-BAL", "NMIDAS-BL", "NMIDAS-AGG")

plot2.crps <- ggplot(data=crps1, aes(x=V1, y=`Continuous Ranked Probability Score`, group=Model, colour=factor(Model))) +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  geom_line(linewidth=1) +
  geom_point(aes(shape=type), size=3) +
  scale_colour_manual(values=colours, labels=modelnames) +
  #scale_colour_discrete(labels=crpsmodels) +
  scale_shape_manual(values=c(1, 4)) +
  labs(colour="Model", shape="") +
  xlab("Forecast Horizon (Weeks)") +
  ylab("Continuous Ranked Probability Score") +
  scale_x_continuous(guide="prism_minor", breaks=c(1,4,8,12), minor_breaks=c(1:12), labels=c(1,4,8,12))


saveplot <- plot_grid(plot2.mae, plot2.pit, plot2.ls, plot2.crps +
          theme(legend.position="none"), labels="AUTO")
savelegend <- get_legend(plot2.crps) 

pdf(file="./plots/outsampleperf.pdf", width=12, height=8.5)
plot_grid(saveplot, legend=savelegend, rel_widths = c(3, .5))
dev.off()

tiff(file="./plots/outsampleperf.tif", width=2250, height=2200, res=300, compression="lzw")
plot_grid(saveplot, legend=savelegend, rel_widths = c(3, .45))
dev.off()

#DM test (1 step ahead)
#darker colour -> lower pval -> x is better than y
dmdat <- dmpdf[2:17,] #change
diag(dmdat) <- NA
dmdatlong <- gather(dmdat, key=x, value=pvalue)
dmdatlong <- cbind(dmdatlong, y=rep(colnames(dmdat), 16))
dmdatlong[which(dmdatlong$pvalue < 0.05), "signif"] <- "*"
dmdatlong$pvalue <- as.numeric(dmdatlong$pvalue)

plot2.dm <- ggplot(data=dmdatlong, aes(x=x, y=y, fill=pvalue)) +
  geom_tile() +
  theme_bw() +
  labs(x="Base Model", y="Competing Model", title="Diebold-Mariano Test") +
  geom_text(aes(label=signif, colour="p-value < 0.05"), key_glyph="point") +
  scale_colour_manual(name="", values="brown1") +
  scale_fill_continuous(breaks=c(0.05, 0.25, 0.5, 0.75, 0.99))+
  annotate("rect", xmin=c(12.5, 0.5), xmax=c(16.5, 16.5), ymin=c(0.5,12.5), ymax=c(16.5,16.5), colour="coral", fill="transparent", linewidth=1) +
  theme(plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title=element_text(face="bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

pdf(file="./plots/dieboldmariano.pdf", width=8.98, height=7.03)
plot2.dm
dev.off()

tiff(file="./plots/dieboldmariano.tif", width=1800, height=1600, res=300, compression="lzw")
plot2.dm
dev.off()

plots <- c(ls(pattern="(plot1)"), ls(pattern="(plot2)"))
save(plots,file="./out/plots.Rdata")


