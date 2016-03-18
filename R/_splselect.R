#' selects a subset of spline dvs based on Maxent FTVE
#'
#' @param rv Vector of response variable values
#' @param dv Dataframe of spline dvs to be selected from (HF, HR, or Th)
#' @param dir Directory to which Maxent runs are
#'   written
#' @param jarpath Pathway to maxent.jar
#'

.splselect <- function(rv, dv, dir, jarpath) {
  for (i in 1:ncol(dv))
    dvname <- colnames(dv)[i]
    df <- data.frame("RV" = rv, "X" = -9999, "Y" = -9999, dv[,i])
    colnames(df)[4] <- dvname

    dvdir <- paste(dir, "\\", dvname, sep = "")
    dir.create(dvdir)

    samplesdf <- na.omit(df)
    environlayersdf <- df

    csvfiles <- paste(dvdir, c("\\samples.csv", "\\environlayers.csv"), sep="")
    write.csv(samplesdf, csvfiles[1], row.names = F)
    write.csv(environlayersdf, csvfiles[2], row.names = F)

    commandsettings <- " threads=8 pictures=FALSE plots=FALSE autofeature=FALSE
    beta_categorical=0 beta_hinge=0 beta_lqp=0 beta_threshold=0 betamultiplier=0
    doclamp=FALSE extrapolate=FALSE jackknife=FALSE
    outputformat=raw outputgrids=FALSE perspeciesresults=FALSE
    hinge=FALSE product=FALSE quadratic=FALSE threshold=FALSE randomtestpoints=0
    removeduplicates=FALSE responsecurves=FALSE responsecurvesexponent=FALSE
    writeplotdata=FALSE warnings=FALSE skipifexists=FALSE redoifexists=TRUE\n"
    command <- paste("java -Xmx512m -jar ",
                     "\"", jarpath, "\"",
                     commandsettings,
                     " samplesfile=","\"", csvfiles[1], "\"", "\n",
                     " outputdirectory=", "\"", dvdir, "\\", "\"", "\n",
                     " environmentallayers=", "\"", csvfiles[2], "\"", "\n",
                     " visible=FALSE",
                     " autorun", sep="")
    command <- gsub("\\\\","/", command)
    system(paste(command), wait=T)
}