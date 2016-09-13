library(MIAmaxent)

?readData
toydata_sp1po <- readData("D:\\Rpackage\\testdata\\sommerfeltia\\Sp1.csv",
  contEV="D:\\Rpackage\\testdata\\sommerfeltia\\EV_continuous\\")

?deriveVars
toydata_dvs <- deriveVars(toydata_sp1po)

?selectDVforEV
toydata_seldvs <- selectDVforEV(toydata_sp1po, toydata_dvs[[1]], alpha = 0.4)

?selectEV
toydata_selevs <- selectEV(toydata_sp1po, toydata_seldvs[[1]], alpha=0.4,
  interaction = TRUE)
# choose lambdas file for round 5 model 1.