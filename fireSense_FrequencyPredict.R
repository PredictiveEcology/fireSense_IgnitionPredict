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
  version = list(SpaDES.core = "0.1.0", fireSense_FrequencyPredict = "0.1.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = NA_character_, # e.g., "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_FrequencyPredict.Rmd"),
  reqdPkgs = list("magrittr", "raster"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", default, min, max, "parameter description")),
    defineParameter(name = "modelName", class = "character", 
                    default = "fireSense_FrequencyFitted",
                    desc = "name of the object of class fireSense_FrequencyFit
                            describing the statistical model used for
                            predictions."),
    defineParameter(name = "data", class = "character",
                    default = "dataFireSense_FrequencyPredict",
                    desc = "a character vector indicating the names of objects 
                            in the `simList` environment in which to look for 
                            variables present in the model formula. `data`
                            objects can be data.frames, RasterStacks or
                            RasterLayers. However, data.frames cannot be mixed
                            with objects of other classes. If variables are not 
                            found in `data` objects, they are searched in the
                            `simList` environment."),
    defineParameter(name = "mapping", class = "character, list", default = NULL,
                    desc = "optional named vector or list of character strings
                            mapping one or more variables in the model formula
                            to those in `data` objects."),
    defineParameter("f", "numeric", 1, 
                    desc = "rescale predicted rates of fire counts at any given 
                            temporal and spatial resolutions by a factor 
                            `f = new_res / old_res`. `f` is the ratio between 
                            the data aggregation scale used for model fitting
                            and the scale at which predictions are to be made."),
    defineParameter(name = ".runInitialTime", class = "numeric",
                    default = start(sim),
                    desc = "when to start this module? By default, the start 
                            time of the simulation."),
    defineParameter(name = ".runInterval", class = "numeric", default = NA, 
                    desc = "optional. Interval between two runs of this module,
                            expressed in units of simulation time."),
    defineParameter(name = ".saveInitialTime", class = "numeric", default = NA, 
                    desc = "optional. When to start saving output to a file."),
    defineParameter(name = ".saveInterval", class = "numeric", default = NA, 
                    desc = "optional. Interval between save events."),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant")
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
      objectClass = "data.frame, RasterLayer, RasterStack",
      sourceURL = NA_character_,
      desc = "One or more objects of class data.frame, RasterLayer or RasterStack in which to look for variables with which to predict."
    )
  ),
  outputObjects = createsOutput(
    objectName = "frequencyPredicted",
    objectClass = NA_character_,
    desc = "An object whose class depends on that of the inputs, could be a RasterLayer or a vector of type numeric."
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.fireSense_FrequencyPredict = function(sim, eventTime, eventType, debug = FALSE)
{
  switch(
    eventType,
    init = { sim <- frequencyPredictInit(sim) },
    run = { sim <- frequencyPredictRun(sim) },
    save = { sim <- frequencyPredictSave(sim) },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  invisible(sim)
}

## event functions
#   - follow the naming convention `modulenameEventtype()`;
#   - `modulenameInit()` function is required for initiliazation;
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
frequencyPredictInit <- function(sim) 
{
  moduleName <- current(sim)$moduleName
  
  sim <- scheduleEvent(sim, eventTime = P(sim)$.runInitialTime, moduleName, "run")
  
  if (!is.na(P(sim)$.saveInitialTime))
    sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, moduleName, "save", .last())
  
  invisible(sim)
}


frequencyPredictRun <- function(sim) 
{
  stopifnot(is(sim[[P(sim)$modelName]], "fireSense_FrequencyFit"))
  
  moduleName <- current(sim)$moduleName
  currentTime <- time(sim, timeunit(sim))
  endTime <- end(sim, timeunit(sim))
  
  ## Toolbox: set of functions used internally by the module
    ## Raster predict function
      fireSense_FrequencyPredictRaster <- function(model, data, sim)
      {
        model %>%
          model.matrix(c(data, sim[[P(sim)$modelName]]$knots)) %>%
          `%*%` (sim[[P(sim)$modelName]]$coef) %>%
          drop %>% sim[[P(sim)$modelName]]$family$linkinv(.) %>%
          `*` (P(sim)$f)
      }
      
    ## Handling piecewise terms in a formula
    pw <- function(v, k) pmax(v - k, 0)
  
  # Load inputs in the data container
  list2env(as.list(envir(sim)), envir = mod)
  
  for (x in P(sim)$data) 
  {
    if (!is.null(sim[[x]])) 
    {
      if (is.data.frame(sim[[x]])) 
      {
        list2env(sim[[x]], envir = mod)
      } 
      else if (is(sim[[x]], "RasterStack"))
      {
        list2env(setNames(unstack(sim[[x]]), names(sim[[x]])), envir = mod)
      } 
      else if (is(sim[[x]], "RasterLayer"))
      {
        next
      } 
      else stop(moduleName, "> '", x, "' is not a data.frame, a RasterLayer or a RasterStack.")
    }
  }
  
  # Define pw() within the data container
  mod$pw <- pw

  terms <- delete.response(terms.formula(sim[[P(sim)$modelName]]$formula))

  ## Mapping variables names to data
  if (!is.null(P(sim)$mapping)) 
  {
    for (i in 1:length(P(sim)$mapping))
    {
      attr(terms, "term.labels") %<>% gsub(
        pattern = names(P(sim)$mapping[i]),
        replacement = P(sim)$mapping[[i]],
        x = .
      )
    }
  }

  formula <- reformulate(attr(terms, "term.labels"), intercept = attr(terms, "intercept"))
  allxy <- all.vars(formula)

  if (!is.null(sim[[P(sim)$modelName]]$knots)) 
  {
    list2env(as.list(sim[[P(sim)$modelName]]$knots), envir = mod)
    kNames <- names(sim[[P(sim)$modelName]]$knots)
    allxy <- allxy[!allxy %in% kNames]
  } 
  else kNames <- NULL

  if (all(unlist(lapply(allxy, function(x) is.vector(mod[[x]])))))
  {
    sim$frequencyPredicted <- (
      formula %>%
        model.matrix(mod) %>%
        `%*%` (sim[[P(sim)$modelName]]$coef) %>%
        drop %>% sim[[P(sim)$modelName]]$family$linkinv(.)
    ) %>% `*` (P(sim)$f)
    
  } 
  else if (all(unlist(lapply(allxy, function(x) is(mod[[x]], "RasterLayer"))))) 
  {
    sim$frequencyPredicted <- mget(allxy, envir = mod, inherits = FALSE) %>%
        stack %>% predict(model = formula, fun = fireSense_FrequencyPredictRaster, na.rm = TRUE, sim = sim)
  } 
  else 
  {
    missing <- !allxy %in% ls(mod, all.names = TRUE)
    
    if (s <- sum(missing))
      stop(moduleName, "> '", allxy[missing][1L], "'",
           if (s > 1) paste0(" (and ", s-1L, " other", if (s>2) "s", ")"),
           " not found in data objects.")

    badClass <- unlist(lapply(allxy, function(x) is.vector(mod[[x]]) || is(mod[[x]], "RasterLayer")))
    
    if (any(badClass))
    {
      stop(moduleName, "> Data objects of class 'data.frame' cannot be mixed with objects of other classes.")
    } 
    else
    {
      stop(moduleName, "> '", paste(allxy[which(!badClass)], collapse = "', '"),
           "' does not match a data.frame's column, a RasterLayer or a RasterStack's layer.")
    }
  }
  
  if (!is.na(P(sim)$.runInterval))
    sim <- scheduleEvent(sim, currentTime + P(sim)$.runInterval, moduleName, "run")
  
  invisible(sim)
}


frequencyPredictSave <- function(sim)
{
  moduleName <- current(sim)$moduleName
  timeUnit <- timeunit(sim)
  currentTime <- time(sim, timeUnit)
  
  saveRDS(
    sim$frequencyPredicted, 
    file = file.path(paths(sim)$out, paste0("fireSense_FrequencyPredicted_", timeUnit, currentTime, ".rds"))
  )
  
  if (!is.na(P(sim)$.saveInterval))
    sim <- scheduleEvent(sim, currentTime + P(sim)$.saveInterval, moduleName, "save", .last())
  
  invisible(sim)
}
