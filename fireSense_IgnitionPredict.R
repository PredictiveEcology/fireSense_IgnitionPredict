defineModule(sim, list(
  name = "fireSense_IgnitionPredict",
  description = "Predict rates of fire frequency from a model fitted using the
                 fireSense_IgnitionFit module. These can be used to feed the
                 ignition component of a landscape fire model (e.g fireSense).",
  keywords = c("fire frequency", "additive property", "poisson", "negative binomial", "fireSense"),
  authors = c(
    person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = c("aut")),
    person("Eliot", "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = "aut"),
    person("Ian", "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = c("ctb"))
  ),
  childModules = character(),
  version = list(SpaDES.core = "0.1.0", fireSense_IgnitionPredict = "0.2.0"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_IgnitionPredict.Rmd"),
  reqdPkgs = list("magrittr", "terra",
                  "PredictiveEcology/fireSenseUtils@development (>=0.0.5.9090)"),
  loadOrder = list(after = "fireSense_dataPrepPredict"),
  parameters = bindrows(
    defineParameter("ignitionFit_Predict_Package", "character", "glmmTMB", NA, NA,
                    desc = paste("The package used to fit the ignitionFit model.",
                                 "It wil be loaded using Require.")),
    defineParameter(".runInitialTime", "numeric", start(sim), NA, NA,
                    desc = "when to start this module? By default, the start
                            time of the simulation."),
    defineParameter(".runInterval", "numeric", 1, NA, NA,
                    desc = "optional. Interval between two runs of this module,
                            expressed in units of simulation time. By default, 1 year."),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    desc = "optional. When to start saving output to a file."),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    desc = "optional. Interval between save events."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    desc = paste("Should this entire module be run with caching activated?",
                                 "This is generally intended for data-type modules,",
                                 "where stochasticity and time are not relevant"))
  ),
  inputObjects = bindrows(
    expectsInput("fireSense_EscapeFitted", "fireSense_EscapeFit",
                 desc = "An object of class `fireSense_EscapeFit` created with the `fireSense_IgnitionFit` module.",
                 sourceURL = NA),
    expectsInput("fireSense_IgnitionFitted", "fireSense_IgnitionFit",
                 desc = "An object of class `fireSense_IgnitionFit` created with the `fireSense_IgnitionFit` module.",
                 sourceURL = NA),
    expectsInput("fireSense_IgnitionAndEscapeCovariates", c("data.table", "SpatRaser"),
                 desc = paste("An object of class `SpatRaster` (named according to variables)",
                              "or `data.frame`/`data.table` with prediction variables.",
                              "If a `data.frame`/`data.table`, then a",
                              "column named 'pixelID' needs to be supplied")),
    expectsInput("flammableRTM", "SpatRaster", sourceURL = NA,
                 desc = "RTM without ice/rocks/urban/water. Flammable map with 0 and 1."),
  ),
  outputObjects = bindrows(
    createsOutput("fireSense_IgnitionProb", "SpatRaster",
                  desc = "a raster layer of ignition probabilities"),
    createsOutput("ignitionsAndEscapes", "data.table",
                 "A data.table containing pixelID (referencing flammableRTM), ignition, and escape")
  )
))

doEvent.fireSense_IgnitionPredict = function(sim, eventTime, eventType, debug = FALSE) {
  moduleName <- currentModule(sim)

  switch(
    eventType,
    init = {

      Require(P(sim)$ignitionFit_Predict_Package)

      sim <- scheduleEvent(sim, eventTime = P(sim)$.runInitialTime, moduleName, "run", eventPriority = 5.11)

      if (!is.na(P(sim)$.saveInitialTime))
        sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, moduleName, "save", .last())
    },
    run = {
      sim <- IgnitionPredictRun(sim)
      if (!is.na(P(sim)$.runInterval))
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.runInterval, moduleName, "run", eventPriority = 5.11)
    },
    save = {
      sim <- IgnitionPredictSave(sim)

      if (!is.na(P(sim)$.saveInterval))
        sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, moduleName, "save", .last())
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

IgnitionPredictRun <- function(sim) {

  ## checks
  if (is.null(sim$fireSense_IgnitionFitted$lambdaRescaleFactor)) {
    sim$fireSense_IgnitionFitted$lambdaRescaleFactor <- 1
  }

  fireSense_IgnitionCovariates <- sim$fireSense_IgnitionAndEscapeCovariates

  nonNaPixels <- sim$fireSense_IgnitionAndEscapeCovariates$pixelID
  dataForPredict <- na.omit(fireSense_IgnitionCovariates)

  if (!is.null(sim$fireSense_IgnitionFitted$rescales)) {
    dataForPredict <- rescaleVarsByMagnitude(dataForPredict,
                                             sim$fireSense_IgnitionFitted$rescales)
  }

  sim$fireSense_IgnitionAndEscapeCovariates[, igProb := predictIgnition(model = sim$fireSense_IgnitionFitted$model,
                                                                        dataForPredict,
                                                                        rescaleFactor = 1,
                                                                        sim$fireSense_IgnitionFitted$lambdaRescaleFactor)]
  # Create outputs

  ## Ignite - accounting for spatial resolution of models
  igDisAggFactor <- sim$fireSense_IgnitionFitted$fittingRes/c(res(sim$flammableRTM)[1])

  igs <- as.data.table(sim$fireSense_IgnitionPredicted, cells = TRUE)
  igs[, ignited := rpois(n = length(ignitionProbs),
                         prob = ignitionProbs)]
  igRas <- rast(sim$fireSense_IgnitionPredicted, vals = igs$cells)
  igRas <- disagg(igRas, fact= igDisAggFactor)
  igRas[sim$flammableRTM[] != 1] <- NA
  igRas <- as.data.table(igRas, cells = TRUE)
  setnames(igRas, new = c("pixelID", "chunkyPixels"))
  igs <- igs[igRas, on = c("cells" = "chunkyPixels")]

  #now, randomly sample
  #use unique because of multiple pixelID per chunkyPixel
  drawList <- unique(igs[ignited > 0, .(chunkyPixels, ignited)])
  pixelID_igs <- sapply(1:nrow(drawList),
                        FUN = function(n, drawList = drawList, igs = igs){
                          sampledChunk <- drawList[n,]
                          pixelID_pop <- igs[chunkyPixel == sampledChunk$chunkyPixel]
                          drawn <- sample(pixelID_pop$pixelID, size = sampledChunk$ignited)
                          return(drawn)
                        })

  ignited <- sample(ignited) # Randomize order




  return(invisible(sim))
}

IgnitionPredictSave <- function(sim) {
  timeUnit <- timeunit(sim)
  currentTime <- time(sim, timeUnit)

  writeRaster(
    sim$fireSense_IgnitionPredicted,
    filename = file.path(paths(sim)$out, paste0("fireSense_IgnitionPredicted_", timeUnit, currentTime, ".tif"))
  )

  return(invisible(sim))
}

.inputObjects <- function(sim) {
  # cacheTags <- c(currentModule(sim), "otherFunctions:.inputObjects")
  # dPath <- asPath(inputPath(sim), 1)
  # message(currentModule(sim), ": using dataPath '", dPath, "'.")

  return(invisible(sim))
}
