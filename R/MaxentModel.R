# Hello, world!
#
# This is a function named 'MaxentModel' which reproduces a model
# produced by maxent.jar in R, from the model's lambda file.

# under construction

MaxentModel <- function(file) { # file is the pathway to the lambda file of the desired Maxent model
  lambdas <- read.csv(file, header = FALSE)
  dvrows <- lambdas[1:(nrow(lambdas)-4),]
  m <- nrow(dvrows)
  thetas <- dvrows[,2]
  xmins <- dvrows[,3]
  xmaxs <- dvrows[,4]
  linPredNorm <- lambdas[nrow(lambdas)-3,2]
  densNorm <- lambdas[nrow(lambdas)-2,2]
  Maxentpredict <- function(dvvalues) { # dvvalues is a vector containing values for all variables included in the model
    summationvector <- numeric(length = m)
    for (i in 1:m) {
      summationvector[i] <- thetas[i]*((dvvalues[i]-xmins[i])/(xmaxs[i]-xmins[i]))
    }
    rawoutput <- ((exp(sum(summationvector)-linPredNorm))/densNorm)
    return(rawoutput)
  }
  message(paste("The function Maxentpredict() now calculates raw output for the model specified by:", "\n",
    file, "\n", "The dvvalues argument must be a vector of length ", m, sep=""))
}
