library(MIAmaxent)

?readData
sp1po <- readData("D:\\Rpackage\\testdata\\sommerfeltia\\Sp1.csv", contEV="D:\\Rpackage\\testdata\\sommerfeltia\\EV_continuous\\")

?deriveVars
DVs <- deriveVars(sp1po)

?selectDVforEV
seldvs <- selectDVforEV(sp1po, DVs[[1]], alpha = 0.4)

?selectEV
selevs <- selectEV(sp1po, seldvs[[1]], alpha=0.4, interaction = TRUE)
# choose lambdas file for round 5 model 1.