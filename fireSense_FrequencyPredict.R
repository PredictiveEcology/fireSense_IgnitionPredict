# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireSense_FrequencyPredict",
  description = "Predict fire frequencies from a model fitted using fireSense_FrequencyFit.
                 These can be used to feed the ignition component of a fire simulation 
                 model (e.g fireSense).",
  keywords = c("fire frequency", "additive property", "poisson", "negative binomial", "fireSense"),
  authors = c(person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = c("aut", "cre"))),
  childModules = character(),
  version = numeric_version("1.2.0.9004"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = NA_character_, # e.g., "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_FrequencyPredict.Rmd"),
  reqdPkgs = list("magrittr", "raster"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", default, min, max, "parameter description")),
    defineParameter("sf", "numeric", 1, 
      desc = "numeric. Rescale predicted rates of fire count at any given temporal and spatial
              resolutions by a factor sf = new_res / old_res. sf describes the ratio between the
              scale of data aggregation at which the statistical model was fitted to the scale at
              which predictions should be made."),
    defineParameter(name = "data", class = "character", default = NULL,
      desc = "optional. A character vector indicating the names of objects present in the simList
              environment, in which to look for variables with which to predict. Objects can be
              data.frames or named lists of RasterLayers. However, objects of different classes
              cannot be mixed. For example, variables cannot be searched simultaneously within an
              object of class data.frame and within an object of class RasterLayer. If omitted, 
              or if variables are not found in the data objects, variables are searched in the
              simList environment."),
    defineParameter(name = "mapping", class = "character", default = NULL,
      desc = "optional. Named character vector to map variable names in the formula to those in
              data objects. Names of unmapped variables are used directly to look for variables
              in data objects or in the simList environment."),
    defineParameter(name = "initialRunTime", class = "numeric", default = NA, 
      desc = "optional. Simulation time at which to start this module. If omitted, start at start(simList)."),
    defineParameter(name = "intervalRunModule", class = "numeric", default = NA, 
      desc = "optional. Interval in simulation time units between two runs of this module.")
  ),
  inputObjects = data.frame(
    objectName = "fireSense_FrequencyFitted",
    objectClass = "fireSense_FrequencyFit",
    sourceURL = "",
    other = NA_character_,
    stringsAsFactors = FALSE
  ),
  outputObjects = data.frame(
    objectName = NA_character_,
    objectClass = NA_character_,
    other = NA_character_,
    stringsAsFactors = FALSE
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.fireSense_FrequencyPredict = function(sim, eventTime, eventType, debug = FALSE) {
  if (eventType == "init") {
    sim <- sim$fireSense_FrequencyPredictInit(sim)

  } else if (eventType == "run") {
    
    sim <- sim$fireSense_FrequencyPredictRun(sim)

  } else if (eventType == "save") {
    # ! ----- EDIT BELOW ----- ! #
    # do stuff for this event
    
    # e.g., call your custom functions/methods here
    # you can define your own methods below this `doEvent` function
    
    # schedule future event(s)
    
    # e.g.,
    # sim <- scheduleEvent(sim, time(sim) + increment, "fireSense_FrequencyPredict", "save")
    
    # ! ----- STOP EDITING ----- ! #
    
  } else {
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  }
  invisible(sim)
}

## event functions
#   - follow the naming convention `modulenameEventtype()`;
#   - `modulenameInit()` function is required for initiliazation;
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
fireSense_FrequencyPredictInit <- function(sim) {

  sim <- scheduleEvent(sim, eventTime = if (is.na(p(sim)$initialRunTime)) start(sim) else p(sim)$initialRunTime, "fireSense_FrequencyPredict", "run")
  sim
  
}


fireSense_FrequencyPredictRun <- function(sim) {

  ## Toolbox: set of functions used internally by the module
    ## Raster predict function
      fireSense_FrequencyPredictRaster <- function(model, data, sim) {

        model %>%
          model.matrix(c(data, sim$fireSense_FrequencyFitted$knots)) %>%
          `%*%` (sim$fireSense_FrequencyFitted$coef) %>%
          drop %>% sim$fireSense_FrequencyFitted$family$linkinv(.)
        
      }
      
    ## Handling piecewise terms in a formula
      pw <- function(v, k) pmax(v - k, 0)

  envData <- new.env(parent = envir(sim))
  on.exit(rm(envData))
  list2env(as.list(envir(sim)), envir = envData)
  envData$pw <- pw
  
  if (!is.null(p(sim)$data))
    lapply(p(sim)$data, function(x, envData) if (is.list(sim[[x]])) list2env(sim[[x]], envir = envData), envData = envData)

  terms <- delete.response(terms.formula(sim$fireSense_FrequencyFitted$formula))
  
  ## Mapping variables names to data
  if (!is.null(p(sim)$mapping)) {
    
    for (i in 1:length(p(sim)$mapping)) {
      
      attr(terms, "term.labels") <- gsub(pattern = names(p(sim)$mapping[i]),
                                         replacement = p(sim)$mapping[i], x = attr(terms, "term.labels"))

    }
    
  }
  
  formula <- reformulate(attr(terms, "term.labels"), intercept = attr(terms, "intercept"))
  allxy <- all.vars(formula)

  if (!is.null(sim$fireSense_FrequencyFitted$knots)) {
    
    list2env(as.list(sim$fireSense_FrequencyFitted$knots), envir = envData)
    kNames <- names(sim$fireSense_FrequencyFitted$knots)
    allxy <- allxy[!allxy %in% kNames]
    
  } else {
    
    kNames <- NULL
    
  }

  if (all(unlist(lapply(allxy, function(x) is.vector(envData[[x]]))))) {
    
    sim$fireSense_FrequencyPredict <- formula %>%
      model.matrix(envData) %>%
      `%*%` (sim$fireSense_FrequencyFitted$coef) %>%
      drop %>% sim$fireSense_FrequencyFitted$family$linkinv(.)

  } else if (all(unlist(lapply(allxy, function(x) is(envData[[x]], "RasterLayer"))))) {

    sim$fireSense_FrequencyPredict <- mget(allxy, envir = envData, inherits = FALSE) %>%
      stack %>% predict(model = formula, fun = fireSense_FrequencyPredictRaster, na.rm = TRUE, sim = sim)

  } else {
    
    varsExist <- allxy %in% ls(envData)
    varsClass <- unlist(lapply(allxy, function(x) is.data.frame(envData[[x]]) || is(envData[[x]], "RasterLayer")))

    if (any(!varsExist)) {
      stop(paste0("fireSense_FrequencyPredict> Variable '", allxy[which(!varsExist)[1L]], "' not found."))
    } else if (any(varsClass)) {
      stop("fireSense_FrequencyPredict> Data objects are not of the same class (e.g. data.frames).")
    } else {
      stop(paste0("fireSense_FrequencyPredict> Variable '", allxy[which(!varsClass)[1L]], "' is not of a data.frame or a RasterLayer."))
    }
  }
  
  if (!is.na(p(sim)$intervalRunModule))
    sim <- scheduleEvent(sim, time(sim) + p(sim)$intervalRunModule, "fireSense_FrequencyPredict", "run")
  
  sim
  
}


### template for save events
fireSense_FrequencyPredictSave <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  sim <- saveFiles(sim)
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}
