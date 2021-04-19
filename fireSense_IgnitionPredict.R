defineModule(sim, list(
  name = "fireSense_IgnitionPredict",
  description = "Predict rates of fire frequency from a model fitted using the
                 fireSense_IgnitionFit module. These can be used to feed the
                 ignition component of a landscape fire model (e.g fireSense).",
  keywords = c("fire frequency", "additive property", "poisson", "negative binomial", "fireSense"),
  authors = c(person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = c("aut", "cre"))),
  childModules = character(),
  version = list(SpaDES.core = "0.1.0", fireSense_IgnitionPredict = "0.2.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_IgnitionPredict.Rmd"),
  reqdPkgs = list("magrittr", "raster",
                  "PredictiveEcology/fireSenseUtils@development (>=0.0.4.9080)"),
  parameters = bindrows(
    defineParameter(name = ".runInitialTime", class = "numeric", default = start(sim),
                    desc = "when to start this module? By default, the start
                            time of the simulation."),
    defineParameter(name = ".runInterval", class = "numeric", default = 1,
                    desc = "optional. Interval between two runs of this module,
                            expressed in units of simulation time. By default, 1 year."),
    defineParameter(name = ".saveInitialTime", class = "numeric", default = NA,
                    desc = "optional. When to start saving output to a file."),
    defineParameter(name = ".saveInterval", class = "numeric", default = NA,
                    desc = "optional. Interval between save events."),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant")
  ),
  inputObjects = bindrows(
    expectsInput(objectName = "covMinMax_ignition",
                 objectClass = "data.table",
                 desc = "Table of the original ranges (min and max) of covariates"),
    expectsInput(objectName = "fireSense_IgnitionFitted", objectClass = "fireSense_IgnitionFit",
                 desc = "An object of class fireSense_IgnitionFit created with the fireSense_IgnitionFit module."),
    expectsInput(objectName = "fireSense_IgnitionAndEscapeCovariates", objectClass = c("data.frame", "data.table", "RasterStack"),
                 desc = paste("An object of class RasterStack (named according to variables)",
                              "or data.frame/data.table with prediction variables. If a data.frame/data.table, then a",
                              "column named 'pixelID' needs to be supplied when outputAsRaster is TRUE"),
                 sourceURL = NA_character_),
    expectsInput(objectName = "flammableRTM", objectClass = "RasterLayer",
                 desc = paste("OPTIONAL. A raster with values of 1 for every flammable pixel, required if",
                              "class(fireSense_IgnitionAndEscapeCovariates) == 'data.table'"))
  ),
  outputObjects = bindrows(
    createsOutput(objectName = "fireSense_IgnitionPredicted", objectClass = "RasterLayer",
                  desc = "a raster layer of ignition probabilities"),
    createsOutput(objectName = "fireSense_IgnitionPredictedVec", objectClass = "numeric",
                  desc = paste("a named numeric vector ignition probabilities, with names",
                               "corresponding to non-NA pixels in 'fireSense_IgnitionPredicted'",
                               "and 'flammableRTM'"))
  )
))

doEvent.fireSense_IgnitionPredict = function(sim, eventTime, eventType, debug = FALSE) {
  moduleName <- currentModule(sim)

  switch(
    eventType,
    init = {
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
    stop("sim$fireSense_IgnitionFitted$lambdaRescaleFactor must be non-NULL and > 0")
  }

  isRasterStack <- is(sim$fireSense_IgnitionAndEscapeCovariates,  "RasterStack")
  covsUsed <- rownames(attr(terms(sim$fireSense_IgnitionFitted$formula[-2]), "factors"))
  covsUsed <- grep("pw", covsUsed, invert = TRUE, value = TRUE)

  if (isRasterStack) {
    fireSense_IgnitionCovariates <- as.data.table(sim$fireSense_IgnitionAndEscapeCovariates[[covsUsed]][])
    rasterTemplate <- raster(sim$fireSense_IgnitionAndEscapeCovariates[[1]])
    nonNaPixels <- which(rowSums(is.na(fireSense_IgnitionCovariates)) != NCOL(fireSense_IgnitionCovariates))
    fireSense_IgnitionCovariates <- fireSense_IgnitionCovariates[nonNaPixels]
  } else {
    fireSense_IgnitionCovariates <-
      if (!is(sim$fireSense_IgnitionAndEscapeCovariates, "data.table")) {
        as.data.table(sim$fireSense_IgnitionAndEscapeCovariates)
      } else {
        sim$fireSense_IgnitionAndEscapeCovariates
      }

    fireSense_IgnitionCovariates <- fireSense_IgnitionCovariates[, ..covsUsed]
    ## checks
    if (is.null(sim$flammableRTM)) {
      stop("'fireSense_IgnitionAndEscapeCovariates' is a table. Please supply 'flammableRTM'")
    }
    if (!"pixelID" %in% colnames(sim$fireSense_IgnitionAndEscapeCovariates)) {
      stop("fireSense_IgnitionAndEscapeCovariates must have a 'pixelID' column")
    }
    nonNaPixels <- sim$fireSense_IgnitionAndEscapeCovariates$pixelID
    rasterTemplate <- sim$flammableRTM
  }

  rescaleFactor <- (raster::res(rasterTemplate)[1]/sim$fireSense_IgnitionFitted$fittingRes)^2

  if (!is.null(sim$fireSense_IgnitionFitted$rescales)) {
    rescaledLayers <- names(sim$fireSense_IgnitionFitted$rescales)

    # # rescale the relevant values
    rescaledVals <-
      Map(r = fireSense_IgnitionCovariates[ , ..rescaledLayers],
          cmm = sim$covMinMax_ignition[, ..rescaledLayers],
          rescaler = sim$fireSense_IgnitionFitted$rescales[rescaledLayers],
          function(r, cmm, rescaler) {
            if (grepl("rescale", rescaler)) {
              rescaleKnown2(r[], 0, 1, min(cmm), max(cmm))
            } else {
              eval(parse(text = rescaleFun), env = fireSense_IgnitionCovariates)
            }
          })
    # # update original object
    fireSense_IgnitionCovariates[, eval(rescaledLayers) := rescaledVals]
  }

  if (!is.null(sim$fireSense_IgnitionFitted$knots)) {
    knots <- as.list(sim$fireSense_IgnitionFitted$knots)
    dataForPredict <- data.frame(fireSense_IgnitionCovariates[], knots)
  } else {
    dataForPredict <- data.frame(fireSense_IgnitionCovariates[])
  }
  dataForPredict <- na.omit(dataForPredict[])
  mu <- predictIgnition(sim$fireSense_IgnitionFitted[["formula"]][-2],
                        dataForPredict,
                        sim$fireSense_IgnitionFitted$coef,
                        rescaleFactor,
                        sim$fireSense_IgnitionFitted$lambdaRescaleFactor,
                        sim$fireSense_IgnitionFitted$family$linkinv)
  # Create outputs
  sim$fireSense_IgnitionPredicted <- raster(rasterTemplate)
  sim$fireSense_IgnitionPredicted[nonNaPixels] <- mu
  sim$fireSense_IgnitionPredictedVec <- mu
  names(sim$fireSense_IgnitionPredictedVec) <- as.character(nonNaPixels)

  return(invisible(sim))
}

IgnitionPredictSave <- function(sim) {
  timeUnit <- timeunit(sim)
  currentTime <- time(sim, timeUnit)

  raster::writeRaster(
    sim$fireSense_IgnitionPredicted,
    filename = file.path(paths(sim)$out, paste0("fireSense_IgnitionPredicted_", timeUnit, currentTime, ".tif"))
  )

  return(invisible(sim))
}

.inputObjects <- function(sim) {
  # cacheTags <- c(currentModule(sim), "function:.inputObjects")
  # dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  # message(currentModule(sim), ": using dataPath '", dPath, "'.")

  return(invisible(sim))
}
