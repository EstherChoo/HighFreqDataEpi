#folder <- "C:/Users/ASUS/Downloads/Acads/FYP/MyCode/"
folder <- "/Users/estherchoo/Library/CloudStorage/OneDrive-NanyangTechnologicalUniversity/FYP/MyCode"
setwd(folder)
sapply(paste0("./functions/", list.files("./functions")), source)
load(file="./out/1_insampleresults.Rdata")
load(file="./data.Rdata")
load("agg.Rdata")

#####IN-SAMPLE#####

###summarise results from mcmc
in.summary.ar1 <- summariseMCMC(out.mcmcar1, 1000)
in.summary.ar2 <- summariseMCMC(out.mcmcar2, 1000)
in.summary.ar3 <- summariseMCMC(out.mcmcar3, 1000)

in.summary.bl1 <- summariseMCMC(out.mcmcbl1, 1000)
in.summary.bl2 <- summariseMCMC(out.mcmcbl2, 1000)
in.summary.bl3 <- summariseMCMC(out.mcmcbl3, 1000)

in.summary.bal1 <- summariseMCMC(out.mcmcbal1, 1000)
in.summary.bal2 <- summariseMCMC(out.mcmcbal2, 1000)
in.summary.bal3 <- summariseMCMC(out.mcmcbal3, 1000)

in.summary.bgl1 <- summariseMCMC(out.mcmcbgl1, 1000)
in.summary.bgl2 <- summariseMCMC(out.mcmcbgl2, 1000)
in.summary.bgl3 <- summariseMCMC(out.mcmcbgl3, 1000)

in.summary.nmar <- summarisenmidas(out.nmidasar, 1000)
in.summary.nmbl <- summarisenmidas(out.nmidasbl, 1000)
in.summary.nmbal <- summarisenmidas(out.nmidasbal, 1000)
in.summary.nmagg <- summarisenmidas(out.nmidasagg, 1000)

#needed for traceplots
save(list=grep("in.summary", ls(), value=T),
  file="./out/3_mcmcsummary.Rdata")

###use summary to get point forecasts
in.pred.ar1 <- predictMCMC(in.summary.ar1, z.high, z.low, 1)
in.pred.ar2 <- predictMCMC(in.summary.ar2, z.high, z.low, 2)
in.pred.ar3 <- predictMCMC(in.summary.ar3, z.high, z.low, 3)

in.pred.bl1 <- predictMCMC(in.summary.bl1, z.high, z.low, 1)
in.pred.bl2 <- predictMCMC(in.summary.bl2, z.high, z.low, 2)
in.pred.bl3 <- predictMCMC(in.summary.bl3, z.high, z.low, 3)

in.pred.bal1 <- predictMCMC(in.summary.bal1, z.high, z.low, 1)
in.pred.bal2 <- predictMCMC(in.summary.bal2, z.high, z.low, 2)
in.pred.bal3 <- predictMCMC(in.summary.bal3, z.high, z.low, 3)

in.pred.bgl1 <- predictMCMC(in.summary.bgl1, z.high, z.low, 1)
in.pred.bgl2 <- predictMCMC(in.summary.bgl2, z.high, z.low, 2)
in.pred.bgl3 <- predictMCMC(in.summary.bgl3, z.high, z.low, 3)

in.pred.nmar <- predictnmidas(in.summary.nmar, z.low)
in.pred.nmbl <- predictnmidas(in.summary.nmbl, z.low)
in.pred.nmbal <- predictnmidas(in.summary.nmbal, z.low)
in.pred.nmagg <- predictnmidas(in.summary.nmagg, z.agg)

save(list=grep("in.pred", ls(), value=T),
     file="./out/3_insamplepoint.Rdata")

###use summary to get density forecasts
in.dens.ar1 <- predictdensMCMC(in.summary.ar1, z.high, z.low, 1)
in.dens.ar2 <- predictdensMCMC(in.summary.ar2, z.high, z.low, 2)
in.dens.ar3 <- predictdensMCMC(in.summary.ar3, z.high, z.low, 3)

in.dens.bl1 <- predictdensMCMC(in.summary.bl1, z.high, z.low, 1)
in.dens.bl2 <- predictdensMCMC(in.summary.bl2, z.high, z.low, 2)
in.dens.bl3 <- predictdensMCMC(in.summary.bl3, z.high, z.low, 3)

in.dens.bal1 <- predictdensMCMC(in.summary.bal1, z.high, z.low, 1)
in.dens.bal2 <- predictdensMCMC(in.summary.bal2, z.high, z.low, 2)
in.dens.bal3 <- predictdensMCMC(in.summary.bal3, z.high, z.low, 3)

in.dens.bgl1 <- predictdensMCMC(in.summary.bgl1, z.high, z.low, 1)
in.dens.bgl2 <- predictdensMCMC(in.summary.bgl2, z.high, z.low, 2)
in.dens.bgl3 <- predictdensMCMC(in.summary.bgl3, z.high, z.low, 3)

in.dens.nmar <- predictdensnmidas(in.summary.nmar, z.low)
in.dens.nmbl <- predictdensnmidas(in.summary.nmbl, z.low)
in.dens.nmbal <- predictdensnmidas(in.summary.nmbal, z.low)
in.dens.nmagg <- predictdensnmidas(in.summary.nmagg, z.agg)

save(list=grep("in.dens", ls(), value=T),
     file="./out/3_insampledensity.Rdata")












  