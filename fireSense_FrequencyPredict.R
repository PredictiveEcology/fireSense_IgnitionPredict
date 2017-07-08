# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireSense_FrequencyPredict",
  description = "Predict rates of fire frequency from a model fitted using the
                 fireSense_FrequencyFit module. These can be used to feed the
                 ignition component of a landscape fire model (e.g fireSense).",
  keywords = c("fire frequency", "additive property", "poisson", "negative binomial", "fireSense"),
  authors = c(person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = c("aut", "cre"))),
  childModules = character(),
  version = numeric_version("0.1.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = NA_character_, # e.g., "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_FrequencyPredict.Rmd"),
  reqdPkgs = list("magrittr", "raster"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", default, min, max, "parameter description")),
    defineParameter("f", "numeric", 1, 
                    desc = "Rescale predicted rates of fire counts at any given 
                            temporal and spatial resolutions by a factor 
                            `f = new_res / old_res`. `f` is the ratio between 
                            the aggregation scale of data to which the
                            statistical model has been fitted and the scale at
                            which predictions are to be made."),
    defineParameter(name = "data", class = "character", default = "dataFireSense_FrequencyPredict",
                    desc = "a character vector indicating the names of objects 
                            in the `simList` environment in which to look for 
                            variables in the model. `data` objects should be
                            data.frames or named lists of RasterLayers. If
                            omitted, or if variables are not found in `data` 
                            objects, variables are searched in  the `simList`
                            environment."),
    defineParameter(name = "mapping", class = "character, list", default = NULL,
                    desc = "optional named vector or list of character strings
                            mapping one or more variables in the model formula
                            to those in data objects."),
    defineParameter(name = "initialRunTime", class = "numeric", default = start(sim),
                    desc = "when to start this module? By default, the start 
                            time of the simulation."),
    defineParameter(name = "intervalRunModule", class = "numeric", default = NA, 
                    desc = "optional. Interval between two runs of this module,
                            expressed in units of simulation time.")
  ),
  inputObjects = rbind(
    expectsInput(
      objectName = "fireSense_FrequencyFitted",
      objectClass = "fireSense_FrequencyFit",
      sourceURL = NA_character_,
      desc = "An object of class fireSense_FrequencyFit created with the fireSense_FrequencyFit module."
    ),
    expectsInput(
      objectName = "dataFireSense_FrequencyPredict",
      objectClass = "data.frame, raster",
      sourceURL = NA_character_,
      desc = "One or more data.frames or named lists of RasterLayers in which to look for variables with which to predict."
    )
  ),
  outputObjects = createsOutput(
    objectName = NA_character_,
    objectClass = NA_character_,
    desc = "An object whose class depends on that of the inputs, could be a raster or a vector of type numeric."
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

  sim <- scheduleEvent(sim, eventTime = P(sim)$initialRunTime, current(sim)$moduleName, "run")
  sim
  
}


fireSense_FrequencyPredictRun <- function(sim) {

  moduleName <- current(sim)$moduleName
  
  ## Toolbox: set of functions used internally by the module
    ## Raster predict function
      fireSense_FrequencyPredictRaster <- function(model, data, sim) {
  
        model %>%
          model.matrix(c(data, sim[[P(sim)$model]]$knots)) %>%
          `%*%` (sim[[P(sim)$model]]$coef) %>%
          drop %>% sim[[P(sim)$model]]$family$linkinv(.) %>%
          `*` (P(sim)$f)
        
      }
      
    ## Handling piecewise terms in a formula
    pw <- function(v, k) pmax(v - k, 0)

  envData <- new.env(parent = envir(sim))
  on.exit(rm(envData))

  # Load data in the container
  list2env(as.list(envir(sim)), envir = envData)
  
  lapply(P(sim)$data, function(x, envData) {
    
    if (!is.null(sim[[x]])) {
      
      if (is.list(sim[[x]]) && !is.null(names(sim[[x]]))) {
        
        list2env(sim[[x]], envir = envData)
        
      } else stop(paste0(moduleName, "> '", x, "' is not a data.frame or a named list."))
      
    }
    
  }, envData = envData)
  
  # Define pw() within the data container
  envData$pw <- pw

  terms <- delete.response(terms.formula(sim[[P(sim)$model]]$formula))

  ## Mapping variables names to data
  if (!is.null(P(sim)$mapping)) {
    
    for (i in 1:length(P(sim)$mapping)) {
      
      attr(terms, "term.labels") <- gsub(
        pattern = names(P(sim)$mapping[i]),
        replacement = P(sim)$mapping[[i]],
        x = attr(terms, "term.labels")
      )

    }
    
  }

  formula <- reformulate(attr(terms, "term.labels"), intercept = attr(terms, "intercept"))
  allxy <- all.vars(formula)

  if (!is.null(sim[[P(sim)$model]]$knots)) {
    
    list2env(as.list(sim[[P(sim)$model]]$knots), envir = envData)
    kNames <- names(sim[[P(sim)$model]]$knots)
    allxy <- allxy[!allxy %in% kNames]
    
  } else {
    
    kNames <- NULL
    
  }

  if (all(unlist(lapply(allxy, function(x) is.vector(envData[[x]]))))) {
    
    sim$fireSense_FrequencyPredict <- (formula %>%
      model.matrix(envData) %>%
      `%*%` (sim[[P(sim)$model]]$coef) %>%
      drop %>% sim[[P(sim)$model]]$family$linkinv(.)) %>%
      `*` (P(sim)$f)

  } else if (all(unlist(lapply(allxy, function(x) is(envData[[x]], "RasterLayer"))))) {

    sim$fireSense_FrequencyPredict <- mget(allxy, envir = envData, inherits = FALSE) %>%
      stack %>% predict(model = formula, fun = fireSense_FrequencyPredictRaster, na.rm = TRUE, sim = sim)

  } else {
    
    exist <- allxy %in% ls(envData)
    class <- unlist(lapply(allxy, function(x) is.data.frame(envData[[x]]) || is(envData[[x]], "RasterLayer")))

    if (any(!exist)) {
      stop(paste0(moduleName, "> Variable '", allxy[which(!exist)[1L]], "' not found."))
    } else if (any(class)) {
      stop(paste0(moduleName, "> Data objects are not of the same class (e.g. data.frames)."))
    } else {
      stop(paste0(moduleName, "> Variable '", allxy[which(!class)[1L]], "' does not match a data.frame's column, a list component, or a RasterLayer."))
    }
  }
  
  if (!is.na(P(sim)$intervalRunModule))
    sim <- scheduleEvent(sim, time(sim) + P(sim)$intervalRunModule, moduleName, "run")
  
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
