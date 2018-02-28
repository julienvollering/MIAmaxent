### Remember to install source MIAmaxent before testing

### Returned values per v.0.5.0.9000


### Grassland

data <- readData(occurrence = system.file("extdata", "occurrence_PO.csv", package = "MIAmaxent"),
                        contEV = system.file("extdata", "EV_continuous", package = "MIAmaxent"),
                        catEV = system.file("extdata", "EV_categorical", package = "MIAmaxent"),
                        maxbkg = 20000)
str(data)

plotFOP(data, "pca1") #[1] 0.001485
plotFOP(data, 2, span=1, intervals=20, ranging=T) #[1] 0.8639161

derived <- deriveVars(data)
length(derived$transformations) #[1] 117
transformations <- derived$transformations
temp <- deriveVars(data, transformtype = c("L", "M", "D", "HF", "HR"), allsplines = T)
length(temp$transformations) #[1] 361

selDVs <- selectDVforEV(derived$dvdata)
sum(sapply(selDVs$dvdata[-1], length)) #[1] 23
temp <- selectDVforEV(derived$dvdata, alpha = 1E-9, test="F")
sum(sapply(temp$dvdata[-1], length)) #[1] 7

selEVs <- selectEV(selDVs$dvdata)
sum(sapply(selEVs$dvdata[-1], length)) #[1] 20
names(selEVs$dvdata) #[1] "RV"        "pr.bygall" "geoberg"   "lcucor1"   "tertpi09"  "pca1"      "geolmja1"  "teraspif"
model <- selEVs$selectedmodel
temp <- selectEV(selDVs$dvdata, alpha = 1E-9, interaction = TRUE, formula = "~geoberg", test="F")
sum(sapply(temp$dvdata[-1], length)) #[1] 13
names(temp$dvdata) #[1] "RV"        "geoberg"   "pr.bygall" "lcucor1"   "tertpi09"

model2 <- chooseModel(selDVs$dvdata, formula = "~ pr.bygall + geoberg + lcucor1 + tertpi09 + pca1 + geolmja1 + teraspif")
all.equal(model$betas, model2$betas) #[1] TRUE

preds <- projectModel(model, transformations, data)
mean(preds$output$PRO) #[1] 1
var(preds$output$PRO) #[1] 0.8483016
EVstack <- raster::stack(c(
  list.files(system.file("extdata", "EV_continuous", package="MIAmaxent"), full.names=TRUE),
  list.files(system.file("extdata", "EV_categorical", package="MIAmaxent"), full.names=TRUE)))
preds <- projectModel(model, transformations, EVstack)
temp <- projectModel(model, transformations, data, clamping=T, raw=T, rescale=T)
sum(temp$output$raw) #[1] 1
var(temp$output$raw) #[1] 2.776624e-09

plotResp(model, transformations, "tertpi09")
plotResp(model, transformations, "tertpi09", logscale=T, lty=5, cex.main=1.5)
plotResp(model, transformations, "pr.tilany") #The 'EV' specified cannot be found in the model

plotResp2(model, transformations, "tertpi09")
plotResp2(model, transformations, "tertpi09", logscale=T, lty=5, cex.main=1.5)
plotResp2(model, transformations, "pr.tilany") #The 'EV' specified cannot be found in the model

trainauc <- testAUC(model, transformations, data) #Warning message: The test data consist of 1/NA only, so NA is treated as absence. Be aware of implications for the interpretation of the AUC value.
trainauc #[1] 0.7010111
datapa <- readData(occurrence = system.file("extdata", "occurrence_PA.csv", package = "MIAmaxent"),
                 contEV = system.file("extdata", "EV_continuous", package = "MIAmaxent"),
                 catEV = system.file("extdata", "EV_categorical", package = "MIAmaxent"),
                 PA=TRUE)
testauc <- testAUC(model, transformations, datapa)
testauc #[1] 0.7891876

rm(list=ls())


### Bradypus

?maxnet::bradypus # version 0.1.2
data("bradypus", package="maxnet")
data <- bradypus
names(data) <- make.names(names(data), allow_ = F)
str(data)
tail(data)

plotFOP(data, "h_dem") #[1] 23.76923
plotFOP(data, "h_dem", span=1, intervals=20, ranging=T) #[1] 0.02092482

derived <- deriveVars(data) #Warning messages: glm.fit: algorithm did not converge (x3)
length(derived$transformations) #[1] 125
transformations <- derived$transformations
temp <- deriveVars(data, transformtype = c("L", "M", "D", "HF", "HR"), allsplines = T)
length(temp$transformations) #[1] 586

selDVs <- selectDVforEV(derived$dvdata)
sum(sapply(selDVs$dvdata[-1], length)) #[1] 19
temp <- selectDVforEV(derived$dvdata, alpha = 1E-9, test="F")
sum(sapply(temp$dvdata[-1], length)) #[1] 6

selEVs <- selectEV(selDVs$dvdata)
sum(sapply(selEVs$dvdata[-1], length)) #[1] 7
names(selEVs$dvdata) #[1] "RV"          "pre6190.l10" "tmn6190.ann" "ecoreg"      "pre6190.l7"
model <- selEVs$selectedmodel
temp <- selectEV(selDVs$dvdata, alpha = 0.2, interaction = TRUE, formula = "~ecoreg", test="F")
sum(sapply(temp$dvdata[-1], length)) #[1] 8
names(temp$dvdata) #[1] "RV"          "ecoreg"      "pre6190.l10" "tmn6190.ann" "pre6190.l7"  "h.dem"

model2 <- chooseModel(selDVs$dvdata, formula = "~ pre6190.l10 + tmn6190.ann + ecoreg + pre6190.l7")
all.equal(model$betas, model2$betas) #[1] TRUE

preds <- projectModel(model, transformations, data)
mean(preds$output$PRO) #[1] 1
var(preds$output$PRO) #[1] 2.151611
temp <- projectModel(model, transformations, data, clamping=T, raw=T, rescale=T)
sum(temp$output$raw) #[1] 1
var(temp$output$raw) #[1] 1.727569e-06

plotResp(model, transformations, "pre6190.l10")
plotResp(model, transformations, "pre6190.l10", logscale=T, lty=5, cex.main=1.5)
plotResp(model, transformations, "h.dem") #The 'EV' specified cannot be found in the model

plotResp2(model, transformations, "pre6190.l10")
plotResp2(model, transformations, "pre6190.l10", logscale=T, lty=5, cex.main=1.5)
plotResp2(model, transformations, "h.dem") #The 'EV' specified cannot be found in the model

trainauc <- testAUC(model, transformations, data)
trainauc #[1] 0.8747716

rm(list=ls())


### Anguilla

?dismo::Anguilla_train # version 1.1-4
data(Anguilla_train, package = "dismo")
data <- Anguilla_train
str(data)
data <- data[, -1]
plot(data$DSDam)
data$DSDam <- as.factor(data$DSDam)
data[which(apply(data, 1, function(x) {any(is.na(x))})), ]
plotFOP(data[complete.cases(data), ], "LocSed")
data <- data[, -match("LocSed", names(data))]

plotFOP(data, "SegSumT") #[1] 19.4
plotFOP(data, "SegSumT", span=1, intervals=20, ranging=T) #[1] 0.9824561

derived <- deriveVars(data, algorithm = "LR") #Warning messages: glm.fit: algorithm did not converge (x1)
length(derived$transformations) #[1] 73
transformations <- derived$transformations
temp <- deriveVars(data, transformtype = c("L", "M", "D", "HF", "HR"), allsplines = T, algorithm = "LR")
length(temp$transformations) #[1] 406

selDVs <- selectDVforEV(derived$dvdata, algorithm = "LR")
sum(sapply(selDVs$dvdata[-1], length)) #[1] 12
temp <- selectDVforEV(derived$dvdata, alpha = 1E-9, test="F", algorithm = "LR")
sum(sapply(temp$dvdata[-1], length)) #[1] 8

selEVs <- selectEV(selDVs$dvdata, algorithm = "LR")
sum(sapply(selEVs$dvdata[-1], length)) #[1] 6
names(selEVs$dvdata) #[1] "RV"         "SegSumT"    "USNative"   "SegTSeas"   "USRainDays" "DSMaxSlope"
model <- selEVs$selectedmodel
temp <- selectEV(selDVs$dvdata, alpha = 0.2, interaction = TRUE, formula = "~USNative + SegTSeas",
                 test="F", algorithm = "LR")
sum(sapply(temp$dvdata[-1], length)) #[1] 7
names(temp$dvdata) #[1] "RV"         "USNative"   "SegTSeas"   "SegSumT"    "USRainDays" "DSMaxSlope" "DSDist"
length(temp$selectedmodel$betas) #[1] 15

model2 <- chooseModel(selDVs$dvdata, formula = "~ SegSumT + USNative + SegTSeas + USRainDays + DSMaxSlope",
                      algorithm = "LR")
all.equal(model$betas, model2$betas) #[1] TRUE

preds <- projectModel(model, transformations, data)
mean(preds$output$response) #[1] 0.202
var(preds$output$response) #[1] 0.04298421
temp <- projectModel(model, transformations, data, clamping=T, raw=T, rescale=T)
mean(temp$output$response) #[1] 0.202
var(temp$output$response) #[1] 0.04298421

plotResp(model, transformations, "SegSumT")
plotResp(model, transformations, "SegSumT", logscale=T, lty=5, cex.main=1.5)
plotResp(model, transformations, "DSDist") #The 'EV' specified cannot be found in the model

plotResp2(model, transformations, "SegSumT")
plotResp2(model, transformations, "SegSumT", logscale=T, lty=5, cex.main=1.5)
plotResp2(model, transformations, "DSDist") #The 'EV' specified cannot be found in the model

trainauc <- testAUC(model, transformations, data)
trainauc #[1] 0.8479677

rm(list=ls())
