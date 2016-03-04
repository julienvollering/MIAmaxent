# Hello, world!
#
# This is a function named 'modelFromLambdas' which reproduces a model
# produced by maxent.jar in R, from the model's lambda file.

# under construction


modelFromLambdas <- function(file) { # file is the pathway to the lambdas file of the desired Maxent model
  lambdas <- read.csv(file, header = FALSE)
  dvrows <- lambdas[1:(nrow(lambdas)-4),]
  m <- nrow(dvrows)
  thetas <- dvrows[,2]
  xmins <- dvrows[,3]
  xmaxs <- dvrows[,4]
  linPredNorm <- lambdas[nrow(lambdas)-3,2]
  densNorm <- lambdas[nrow(lambdas)-2,2]
  function(dvvalues) { # dvvalues must be a vector of length m
    summationvector <- numeric(length = m)
    for (i in 1:m) {
      summationvector[i] <- thetas[i]*((dvvalues[i]-xmins[i])/(xmaxs[i]-xmins[i]))
    }
    rawoutput <- ((exp(sum(summationvector)-linPredNorm))/densNorm)
    return(rawoutput)
  }
}


# below is for development stage testing only

mymodel_1dv <- modelFromLambdas("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3A_All\\M3A_Al1\\EV01\\GP01\\BIO_10_ZSk\\RV01Al1E01P01.lambdas")
samplevalues_1dv <- read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3A_All\\M3A_Al1\\EV01\\GP01\\BIO_10_ZSk_RV01Al1E01P01.csv")[,4]
Rpreds_1dv <- sapply(samplevalues_1dv, mymodel_1dv)
Jarpreds_1dv <- read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3A_All\\M3A_Al1\\EV01\\GP01\\BIO_10_ZSk\\RV01Al1E01P01_samplePredictions.csv")[,4]
plot (Jarpreds_1dv, Rpreds_1dv)
abline(a=0, b=1, col="red")

mymodel_2dv <- modelFromLambdas("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3A_All\\M3A_Al1\\EV01\\GP02\\EV01_GP02_DV1_08\\RV01Al1E01P02.lambdas")
samplevalues_2dv <- read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3A_All\\M3A_Al1\\EV01\\GP02\\EV01_GP02_DV1_08_RV01Al1E01P02.csv")[,4:5]
samplevalues_2dv <- as.matrix(samplevalues_2dv) # dvvalues need to be in vectors
Rpreds_2dv <- apply(samplevalues_2dv, 1, mymodel_2dv)
Jarpreds_2dv <- read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3A_All\\M3A_Al1\\EV01\\GP02\\EV01_GP02_DV1_08\\RV01Al1E01P02_samplePredictions.csv")[,4]
plot (Jarpreds_2dv, Rpreds_2dv)
abline(a=0, b=1, col="red")

mymodel_3dv <- modelFromLambdas("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3A_All\\M3A_Al1\\EV01\\GP03\\EV01_GP03_DV1_06\\RV01Al1E01P03.lambdas")
samplevalues_3dv <- read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3A_All\\M3A_Al1\\EV01\\GP03\\EV01_GP03_DV1_06_RV01Al1E01P03.csv")[,4:6]
samplevalues_3dv <- as.matrix(samplevalues_3dv) # dvvalues need to be in vectors
Rpreds_3dv <- apply(samplevalues_3dv, 1, mymodel_3dv)
Jarpreds_3dv <- read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3A_All\\M3A_Al1\\EV01\\GP03\\EV01_GP03_DV1_06\\RV01Al1E01P03_samplePredictions.csv")[,4]
plot (Jarpreds_3dv, Rpreds_3dv)
abline(a=0, b=1, col="red")

mymodel_45dv <- modelFromLambdas("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3C_All\\M3C_aWI_RV_Al1\\GP03\\GrPr03pv4\\aWI_RV_Al1P03.lambdas")
samplevalues_45dv <- read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3C_All\\M3C_aWI_RV_Al1\\GP03\\GrPr03pv4_aWI_RV_Al1P03.csv")[,4:48]
samplevalues_45dv <- as.matrix(samplevalues_45dv) # dvvalues need to be in vectors
Rpreds_45dv <- apply(samplevalues_45dv, 1, mymodel_45dv)
Jarpreds_45dv <- read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3C_All\\M3C_aWI_RV_Al1\\GP03\\GrPr03pv4\\aWI_RV_Al1P03_samplePredictions.csv")[,4]
plot (Jarpreds_45dv, Rpreds_45dv)
abline(a=0, b=1, col="red")
