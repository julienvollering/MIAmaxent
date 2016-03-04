# Hello, world!
#
# This is a function named 'MaxentModel' which reproduces a model
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
  function(dvvalues) {
    summationvector <- numeric(length = m)
    for (i in 1:m) {
      summationvector[i] <- thetas[i]*((dvvalues[i]-xmins[i])/(xmaxs[i]-xmins[i]))
    }
    rawoutput <- ((exp(sum(summationvector)-linPredNorm))/densNorm)
    return(rawoutput)
  }
}

# below is for development stage only
mymodel <- MaxentModel("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3A_All\\M3A_Al1\\EV01\\GP02\\EV01_GP02_DV1_02\\RV01Al1E01P02.lambdas")
sample_dvvalues <- read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3A_All\\M3A_Al1\\EV01\\GP02\\EV01_GP02_DV1_02_RV01Al1E01P02.csv")
c(sample_dvvalues[1,4],sample_dvvalues[1,5]
csample_dvvalues <- t(unname(apply(sample_dvvalues[,4:5], 1, c)))
csample_dvvalues[1032,]
as.vector(csample_dvvalues[1032,])
Rpreds <- numeric(length=nrow(csample_dvvalues))
for (i in 1:nrow(csample_dvvalues)) {
  Rpreds[i] <- mymodel(csample_dvvalues[i,])
}

Rpreds <- sapply(csample_dvvalues, mymodel, simplify=T)
Jarpreds <- read.csv("D:\\Rpackage\\MIATtest\\MIATtest2\\AM\\M3A_All\\M3A_Al1\\EV01\\GP02\\EV01_GP02_DV1_02\\RV01Al1E01P02_samplePredictions.csv")[,4]
plot (Jarpreds, Rpreds)
abline(a=0, b=1)
