library(MIAmaxent)

?readData
toydata_sp1po <- readData("Z:\\Rpackage\\MIAmaxent\\inst\\extdata\\sommerfeltia\\Sp1.csv",
  contEV="Z:\\Rpackage\\MIAmaxent\\inst\\extdata\\sommerfeltia\\EV_continuous\\")

?deriveVars
toydata_dvs <- deriveVars(toydata_sp1po)

?selectDVforEV
toydata_seldvs <- selectDVforEV(toydata_dvs$dvdata, alpha = 0.4)

?selectEV
toydata_selevs <- selectEV(toydata_seldvs$dvdata, alpha=0.4,
  interaction = TRUE)
