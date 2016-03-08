# Hello, world!
#
# This is a function named 'modelFromLambdas' which reproduces a model
# produced by maxent.jar in R, from the model's lambda file.

# Details:

modelFromLambdas <- function(file) { # file is the pathway to the lambdas file of the desired Maxent model
  lambdas <- read.csv(file, header = FALSE)
  dvrows <- lambdas[1:(nrow(lambdas)-4),]
  m <- nrow(dvrows)
  thetas <- dvrows[,2]
  xmins <- dvrows[,3]
  xmaxs <- dvrows[,4]
  linPredNorm <- lambdas[nrow(lambdas)-3,2]
  densNorm <- lambdas[nrow(lambdas)-2,2]
  function(X) { # X must be an array (matrix or df) with m columns, where column names match variable names in lambda file.
    matchorder <- match(as.character(dvrows[,1]), colnames(X))

    if (ncol(X) != m) {
      stop("Input must have as many columns as there are variables in the model", call. = F)
    }

    if (length((matchorder[!is.na(matchorder)])) < m) {
      stop("Input column names must match the names of the variables in the model", call. = F)
    }

    orderedX <- as.matrix(X[,matchorder])
    thetaX <- matrix(nrow = nrow(orderedX), ncol = m)
    for (j in 1:nrow(orderedX)) {
      for (i in 1:m) {
        thetaX[j,i] <- thetas[i]*((orderedX[j,i]-xmins[i])/(xmaxs[i]-xmins[i]))
      }
    }
    rawoutput <- apply(thetaX, 1, function(x) ((exp(sum(x)-linPredNorm))/densNorm))
    projection <- cbind(orderedX, rawoutput)
    return(projection)
  }
}




# below is for development stage testing only

mymodel_1dv <- modelFromLambdas("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3A_All\\M3A_Al1\\EV01\\GP01\\BIO_10_ZSk\\RV01Al1E01P01.lambdas")
values_1dv <- as.data.frame(read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3A_All\\M3A_Al1\\EV01\\GP01\\BIO_10_ZSk.csv")[,4])
colnames(values_1dv) <- "BIO_10_ZSk"
head(values_1dv)
Rpreds_1dv <- mymodel_1dv(values_1dv)[,"rawoutput"]
Jarpreds_1dv <- read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3A_All\\M3A_Al1\\EV01\\GP01\\BIO_10_ZSk\\RV01Al1E01P01.csv")[,3]
plot (Jarpreds_1dv, Rpreds_1dv)
abline(a=0, b=1, col="red")

mymodel_2dv <- modelFromLambdas("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3A_All\\M3A_Al1\\EV01\\GP02\\EV01_GP02_DV1_08\\RV01Al1E01P02.lambdas")
values_2dv <- read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3A_All\\M3A_Al1\\EV01\\GP02\\EV01_GP02_DV1_08.csv")[,4:5]
head(mymodel_2dv(values_2dv))
Rpreds_2dv <- mymodel_2dv(values_2dv)[,"rawoutput"]
Jarpreds_2dv <- read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3A_All\\M3A_Al1\\EV01\\GP02\\EV01_GP02_DV1_08\\RV01Al1E01P02.csv")[,3]
plot (Jarpreds_2dv, Rpreds_2dv)
abline(a=0, b=1, col="red")

mymodel_3dv <- modelFromLambdas("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3A_All\\M3A_Al1\\EV01\\GP03\\EV01_GP03_DV1_06\\RV01Al1E01P03.lambdas")
values_3dv <- read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3A_All\\M3A_Al1\\EV01\\GP03\\EV01_GP03_DV1_06.csv")[,4:6]
values_3dv <- as.matrix(values_3dv)
head(mymodel_3dv(values_2dv))
Rpreds_3dv <- mymodel_3dv(values_3dv)[,"rawoutput"]
Jarpreds_3dv <- read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3A_All\\M3A_Al1\\EV01\\GP03\\EV01_GP03_DV1_06\\RV01Al1E01P03.csv")[,3]
plot (Jarpreds_3dv, Rpreds_3dv)
abline(a=0, b=1, col="red")

mymodel_4dv <- modelFromLambdas("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3B_All\\M3B_aNI_RV_Al1\\GP02\\GrPr02pv1\\aNI_RV_Al1_P02.lambdas")
samplevalues_4dv <- read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3B_All\\M3B_aNI_RV_Al1\\GP02\\GrPr02pv1_aNI_RV_Al1_P02.csv")[,4:7]
Rpreds_4dv <- mymodel_4dv(samplevalues_4dv)[,"rawoutput"]
Jarpreds_4dv <- read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3B_All\\M3B_aNI_RV_Al1\\GP02\\GrPr02pv1\\aNI_RV_Al1_P02_samplePredictions.csv")[,4]
plot (Jarpreds_4dv, Rpreds_4dv)
abline(a=0, b=1, col="red")

mymodel_45dv <- modelFromLambdas("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3C_All\\M3C_aWI_RV_Al1\\GP03\\GrPr03pv4\\aWI_RV_Al1P03.lambdas")
samplevalues_45dv <- read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3C_All\\M3C_aWI_RV_Al1\\GP03\\GrPr03pv4_aWI_RV_Al1P03.csv")[,4:48]
Rpreds_45dv <- mymodel_45dv(samplevalues_45dv)[,"rawoutput"]
Jarpreds_45dv <- read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3C_All\\M3C_aWI_RV_Al1\\GP03\\GrPr03pv4\\aWI_RV_Al1P03_samplePredictions.csv")[,4]
plot (Jarpreds_45dv, Rpreds_45dv)
abline(a=0, b=1, col="red")

mymodel_24dv <- modelFromLambdas("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3C_All\\M3C_aWI_RV_Al2\\GP02\\GrPr02pv10\\aWI_RV_Al2P02.lambdas")
samplevalues_24dv <- read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3C_All\\M3C_aWI_RV_Al2\\GP02\\GrPr02pv10_aWI_RV_Al2P02.csv")[,4:27]
Rpreds_24dv <- mymodel_24dv(samplevalues_24dv)[,"rawoutput"]
Jarpreds_24dv <- read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3C_All\\M3C_aWI_RV_Al2\\GP02\\GrPr02pv10\\aWI_RV_Al2P02_samplePredictions.csv")[,4]
plot (Jarpreds_24dv, Rpreds_24dv)
abline(a=0, b=1, col="red")
