defineModule(sim, list(
  name = "fireSense_IgnitionPredict",
  description = paste(
    "Predicts annual ignition and escape probabilities from the models fitted by",
    "fireSense_IgnitionFit and fireSense_EscapeFit, and draws the pixels that ignite and escape."),
  keywords = c("fire frequency", "additive property", "poisson", "negative binomial", "fireSense"),
  authors = c(
    person("Eliot", "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person("Ian", "Eddy", email = "ian.eddy@nrcan-rncan.gc.ca", role = "aut"),
    person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = "aut"),
    person("Alex M", "Chubaty", email = "achubaty@for-cast.ca", role = "ctb")
  ),
  childModules = character(),
  version = list(SpaDES.core = "0.1.0", fireSense_IgnitionPredict = "1.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_IgnitionPredict.Rmd"),
  reqdPkgs = list(
    "magrittr", "terra",
    "PredictiveEcology/fireSenseUtils@development (>=0.1.0)"
  ),
  loadOrder = list(after = "fireSense_dataPrepPredict"),
  parameters = bindrows(
    defineParameter("modelAlgorithm", "character", "xgboost", NA, NA,
                    paste("Algorithm used to fit the models; only `xgboost` is supported.",
                          "Must agree with the value in the other fireSense modules.")),
    defineParameter("rescaleVars", "logical", default = TRUE,
                    desc = paste("Rescale the covariates before predicting? With `xgboost` they are standardized",
                                 "with `scale()`. Must agree with the value in the other fireSense modules.")),
    defineParameter(".runInitialTime", "numeric", start(sim), NA, NA,
                    desc = "Time of the first prediction."
    ),
    defineParameter(".runInterval", "numeric", 1, NA, NA,
                    desc = "Interval between predictions, in years. `NA` predicts once."
    ),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    desc = "Time of the first `save` event. `NA` means never."
    ),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    desc = paste("Interval between `save` events.",
                                 "If not `NA`, the ignition probability raster is also plotted each year.")
    ),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    desc = paste(
                      "Should this entire module be run with caching activated?",
                      "This is generally intended for data-type modules,",
                      "where stochasticity and time are not relevant"
                    )
    )
  ),
  inputObjects = bindrows(
    expectsInput("fireSense_EscapeFitted", "fireSense_EscapeFit",
                 desc = "Fitted escape models (`$modelList$model`, one per fold), from `fireSense_EscapeFit`.",
                 sourceURL = NA
    ),
    expectsInput("fireSense_IgnitionFitted", "fireSense_IgnitionFit",
                 desc = paste("Fitted ignition models (`$modelList$model`, one per fold) and `$modelList$fittingRes`,",
                              "from `fireSense_IgnitionFit`."),
                 sourceURL = NA
    ),
    expectsInput("fireSense_igAndEscapePred_Covariates", "data.table",
                 desc = paste(
                   "This year's covariates, from `fireSense_dataPrepPredict`.",
                   "`pixelID` is the cell index of `ignitionFitRTM`."
                 )
    ),
    expectsInput("flammableRTM", "SpatRaster",
                 sourceURL = NA,
                 desc = "Binary raster, 1 where the pixel is flammable."
    ),
  ),
  outputObjects = bindrows(
    createsOutput("fireSense_IgAndEscapeProbRas", "SpatRaster",
                  desc = paste("Two layers, `ignitionProb` (expected ignitions per pixel) and `escapeProb`,",
                               "at the resolution of `ignitionFitRTM`.")
    ),
    createsOutput(
      "ignitionsAndEscapes", "data.table",
      paste("One row per ignited pixel, in random order: `pixelID` (cell index of `flammableRTM`),",
            "and `igProb`, `ignitions`, `escapeProb`, `escapes` of the coarse pixel it was drawn from.")
    )
  )
))

#' Event dispatcher
#'
#' Events: `init`, `run` (predict, repeated every `.runInterval`), `save`.
#'
#' @param sim A `simList`.
#' @param eventTime Time of the event.
#' @param eventType Name of the event.
#' @param debug Not used.
#'
#' @return The `simList`, invisibly.
doEvent.fireSense_IgnitionPredict <- function(sim, eventTime, eventType, debug = FALSE) {
  moduleName <- currentModule(sim)

  switch(eventType,
         init = {
           sim <- scheduleEvent(sim, eventTime = P(sim)$.runInitialTime, moduleName, "run")

           if (!is.na(P(sim)$.saveInitialTime)) {
             sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, moduleName, "save", .last())
           }
         },
         run = {
           sim <- IgnitionPredictRun(sim)
           if (!is.na(P(sim)$.runInterval)) {
             sim <- scheduleEvent(sim, time(sim) + P(sim)$.runInterval, moduleName, "run")
           }
         },
         save = {
           sim <- IgnitionPredictSave(sim)

           if (!is.na(P(sim)$.saveInterval)) {
             sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, moduleName, "save", .last())
           }
         },
         warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                       "' in module '", current(sim)[1, "moduleName", with = FALSE], "'",
                       sep = ""
         ))
  )
  return(invisible(sim))
}

#' Predict and draw this year's ignitions and escapes
#'
#' Averages the per-fold predictions of the ignition models, draws ignitions per coarse pixel
#' (Poisson), then does the same for escapes (binomial, given ignition). Each ignition is placed in
#' a random flammable `flammableRTM` pixel inside its coarse pixel.
#'
#' @param sim A `simList`.
#'
#' @return The `simList`, invisibly, with `fireSense_IgAndEscapeProbRas` and `ignitionsAndEscapes`.
IgnitionPredictRun <- function(sim) {
  ## checks
  if (is.null(sim$fireSense_IgnitionFitted$lambdaRescaleFactor)) {
    sim$fireSense_IgnitionFitted$lambdaRescaleFactor <- 1
  }

  igCov <- copy(sim$fireSense_igAndEscapePred_Covariates)

  igCov <- na.omit(igCov)
  
  modelAlgorithm <- paramCheckOtherMods(sim, "modelAlgorithm")
  rescaleVars <- paramCheckOtherMods(sim, "rescaleVars")
  data <- prepareCovariatesOuter(igCov,
                                 algorithm = modelAlgorithm, 
                                 rescaleVars = rescaleVars, useCache = FALSE)
  pixelId <- igCov$pixelID
  igCov <- data$covariates
  set(igCov, NULL, "pixelID", pixelId)
  

  # From here, it makes predictions from each KFold, then averages them to get the probabilities for
  #   each cell; same for Escape, which is conditional on having ignited.
  modsHere <- sim$fireSense_IgnitionFitted$modelList$model
  modsOnly <- modsHere[grep("Fold", names(modsHere))]
  if (!is.null(modsOnly)) {
    predsIgnList <- Map(model = modsOnly, function(model) {
      predict(model, newdata = igCov)
    })

    predsIgnsMat <- do.call(cbind, predsIgnList)
    predsIgns <- rowMeans(predsIgnsMat)
    igns <- rpois(NROW(predsIgns), lambda = predsIgns) # can't use rtweedie because don't know the dispersion parameter

    modsEscHere <- sim$fireSense_EscapeFitted$modelList$model
    modsEscOnly <- modsEscHere[grep("Fold", names(modsEscHere))]
    whHasIgns <- which(igns > 0)
    
    predsEscList <- Map(model = modsEscOnly,
                        function(model) {
                          predsEsc <- predict(model, newdata = igCov[whHasIgns])
                          pmax(pmin(1, predsEsc), 0) # the escape model is tweedie, so can go above 1, below 0 rarely
                        })
    predsEscsMat <- do.call(cbind, predsEscList)
    predsEscs <- rowMeans(predsEscsMat)
    escs <- rbinom(NROW(predsEscs), size = igns[whHasIgns], prob = predsEscs)

    set(igCov, NULL, c("igProb", "ignitions"), list(predsIgns, igns))
    set(igCov, whHasIgns, c("escapeProb", "escapes"), list(predsEscs, escs))

    igDisAggFactor <- ceiling(sim$fireSense_IgnitionFitted$modelList$fittingRes / c(res(sim$flammableRTM)[1]))
  } else {
    stop("Not tested anymore")
  }
  # Ignite - accounting for spatial resolution of models

  # disaggregate the coarse raster to size of flammableRTM
  # draw ignitions/escapes from smaller pixel size
  igRas <- rast(sim$ignitionFitRTM)
  igRas[igCov$pixelID] <- igCov$pixelID

  # can't use disagg because it may not be round pixels
  igRas <- postProcess(igRas, to = sim$flammableRTM, method = "near")
  igRas[sim$flammableRTM[] != 1] <- NA

  igRasDT <- as.data.table(igRas, cells = TRUE)
  setnames(igRasDT, new = c("pixelID", "chunkyPixels"))
  # nomatch = NULL --> removes pixels that aren't in igCov, which are non flammable
  igs <- igRasDT[igCov, on = c("chunkyPixels" = "pixelID"), nomatch = NULL]
  # randomly sample
  # use unique because of multiple pixelID per chunkyPixel
  drawTable <- unique(igs[ignitions > 0, .(chunkyPixels, ignitions, escapes)])
  pixelID_igs <- lapply(drawTable$chunkyPixels,
                        FUN = function(n, draw = drawTable, ig = igs) {
                          sampledChunk <- draw[chunkyPixels == n]
                          pixelID_pop <- igs[chunkyPixels == sampledChunk$chunkyPixels]
                          drawn <- sample(pixelID_pop$pixelID, size = sampledChunk$ignitions)
                          igPixels <- try(ig[pixelID %in% drawn, .(pixelID, igProb, ignitions, escapeProb, escapes)])
                          return(igPixels)
                        }
  ) |>
    rbindlist() #will this fail with zero ignitions?

  #for posterior analysis
  igProbRas <- rast(sim$ignitionFitRTM)
  igProbRas[igCov$pixelID] <- igCov$igProb
  escapeProbRas <- igProbRas
  escapeProbRas[igCov$pixelID] <- igCov$escapeProb
  sim$fireSense_IgAndEscapeProbRas <- c(igProbRas, escapeProbRas)
  names(sim$fireSense_IgAndEscapeProbRas) <- c("ignitionProb", "escapeProb")
  
  if (!is.na(P(sim)$.saveInterval)) {
    #TODO: implement proper plotting control
    Plots(igProbRas, types = Par$.plots, filename = paste0("IgnitionProbability_yr", time(sim)))
  }
  randomOrder <- sample(1:nrow(pixelID_igs))
  sim$ignitionsAndEscapes <- pixelID_igs[randomOrder,] #Randomize order

  return(invisible(sim))
}

#' Write the predicted raster to `outputPath`
#'
#' @param sim A `simList`.
#'
#' @return The `simList`, invisibly.
IgnitionPredictSave <- function(sim) {
  timeUnit <- timeunit(sim)
  currentTime <- time(sim, timeUnit)

  writeRaster(
    sim$fireSense_IgnitionPredicted,
    filename = file.path(paths(sim)$out, paste0("fireSense_IgnitionPredicted_", timeUnit, currentTime, ".tif"))
  )

  return(invisible(sim))
}

#' Supply default inputs
#'
#' There are none.
#'
#' @param sim A `simList`.
#'
#' @return The `simList`, invisibly.
.inputObjects <- function(sim) {
  return(invisible(sim))
}
