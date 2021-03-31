defineModule(sim, list(
  name = "fireSense_IgnitionPredict",
  description = "Predict rates of fire frequency from a model fitted using the
                 fireSense_IgnitionFit module. These can be used to feed the
                 ignition component of a landscape fire model (e.g fireSense).",
  keywords = c("fire frequency", "additive property", "poisson", "negative binomial", "fireSense"),
  authors = c(person("Jean", "Marchal", email = "jean.d.marchal@gmail.com", role = c("aut", "cre"))),
  childModules = character(),
  version = list(SpaDES.core = "0.1.0", fireSense_IgnitionPredict = "0.1.0"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fireSense_IgnitionPredict.Rmd"),
  reqdPkgs = list("magrittr", "raster"),
  parameters = rbind(
    # defineParameter("rescaleFactor", "numeric", (250 / 10000)^2,
    #                 desc = paste("rescale predicted rates of fire counts at any given temporal and spatial",
    #                              "resolutions by a factor `rescaleFactor = new_res / old_res`.",
    #                              "`rescaleFactor` is the ratio between the data aggregation scale used",
    #                              "for model fitting and the scale at which predictions are to be made")),
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
  inputObjects = rbind(
    expectsInput(objectName = "fireSense_IgnitionFitted", objectClass = "fireSense_IgnitionFit",
                 desc = "An object of class fireSense_IgnitionFit created with the fireSense_IgnitionFit module."),
    expectsInput(objectName = "fireSense_IgnitionAndEscapeCovariates", objectClass = "data.frame",
                 desc = "An object of class RasterStack or data.frame with prediction variables"),
    expectsInput(objectName = "flammableRTM", objectClass = "RasterLayer",
                 desc = "a raster with values of 1 for every flammable pixel"),
    expectsInput(objectName = "rescaleFactor", objectClass = "numeric",
                 desc = paste("rescale predicted rates of fire counts at any given temporal and spatial",
                              "resolutions by a factor `rescaleFactor = new_res / old_res`.",
                              "`rescaleFactor` is the ratio between the data aggregation scale used",
                              "for model fitting and the scale at which predictions are to be made",
                              "If not provided, defaults to (250 / 10000)^2"))
  ),
  outputObjects = rbind(
    createsOutput(objectName = "fireSense_IgnitionProbRaster", objectClass = "RasterLayer",
                  desc = "a raster layer of ignition probabilities")
  )
))

doEvent.fireSense_IgnitionPredict = function(sim, eventTime, eventType, debug = FALSE)
{
  moduleName <- current(sim)$moduleName

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
  invisible(sim)
}

IgnitionPredictRun <- function(sim) {

  moduleName <- currentModule(sim)

  fireSense_IgnitionCovariates <- copy(sim$fireSense_IgnitionAndEscapeCovariates)
  #TODO: IE wrote this - please review it
  if (!is.null(sim$fireSense_IgnitionFitted$rescales)) {
    for (name in names(sim$fireSense_IgnitionFitted$rescales)) {
      op <- eval(parse(text = sim$fireSense_IgnitionFitted$rescales[[name]]),
                 env = fireSense_IgnitionCovariates)
      fireSense_IgnitionCovariates[, eval(name) := op]
    }
  }

  ## Toolbox: set of functions used internally by IgnitionPredictRun
  IgnitionPredictRaster <- function(model, data, sim) {
    model %>%
      model.matrix(c(data, sim$knots)) %>%
      `%*%` (sim$fireSense_IgnitionFitted$coef) %>%
      drop %>% sim$fireSense_IgnitionFitted$family$linkinv(.) %>%
      `*` (sim$rescaleFactor)
  }

  ## Handling piecewise terms in a formula
  pw <- function(v, k) pmax(v - k, 0)

  # Load inputs in the data container
  #TODO: We can get rid of this 'environment data container' construct
  mod_env <- new.env(parent = baseenv()) #get access to base R

  # Define pw() within the data container
  mod_env$pw <- pw

  terms <- delete.response(terms.formula(sim$fireSense_IgnitionFitted[["formula"]]))

  formula_fire <- reformulate(attr(terms, "term.labels"), intercept = attr(terms, "intercept"))
  allxy <- all.vars(formula_fire)

  #add knots
  #TODO: I'm not 100% sure this works with a linear model only
  if (!is.null(sim$fireSense_IgnitionFitted$knots)){
    for (knot in names(sim$fireSense_IgnitionFitted$knots)) {
      fireSense_IgnitionCovariates[, eval(knot) := sim$fireSense_IgnitionFitted$knots[knot]]
    }
    kNames <- names(sim$fireSense_IgnitionFitted$knots)
    allxy <- allxy[!allxy %in% kNames]
  }

  list2env(fireSense_IgnitionCovariates, env = mod_env)

  if (all(unlist(lapply(allxy, function(x) is.vector(mod_env[[x]]))))) {

    sim$fireSense_IgnitionPredicted <- (
      formula_fire %>%
        model.matrix(mod_env) %>%
        `%*%`(sim$fireSense_IgnitionFitted$coef) %>%
        drop %>% sim$fireSense_IgnitionFitted$family$linkinv(.)
    ) %>% `*`(sim$rescaleFactor)
  } else if (all(unlist(lapply(allxy, function(x) is(mod_env[[x]], "RasterLayer"))))) {
    covList <- mget(allxy, envir = mod_env, inherits = FALSE)
    tryCatch({
      raster::stack(covList)
    }, error = function(e){
      stop("At least one of the covariate rasters does not align with the others. Please debug your inputs.
              Consider using a function like reproducible::postProcess on your layers to make sure these align.")
    })
    sim$fireSense_IgnitionPredicted <- raster::stack(covList) %>%
      raster::predict(model = formula_fire, fun = IgnitionPredictRaster, na.rm = TRUE, sim = sim)
  } else {
    missing <- !allxy %in% ls(mod_env, all.names = TRUE)

    if (s <- sum(missing))
      stop(
        moduleName, "> '", allxy[missing][1L], "'",
        if (s > 1) paste0(" (and ", s - 1L, " other", if (s > 2) "s", ")"),
        " not found in data objects."
      )

    badClass <- unlist(lapply(allxy, function(x) is.vector(mod_env[[x]]) || is(mod_env[[x]], "RasterLayer")))

    if (any(badClass)) {
      stop(moduleName, "> Data objects of class 'data.frame' cannot be mixed with objects of other classes.")
    } else {
      stop(moduleName, "> '", paste(allxy[which(!badClass)], collapse = "', '"),
           "' does not match a data.frame's column, a RasterLayer or a layer from a RasterStack or RasterBrick.")
    }
  }

  #force ignitionPredicted a raster
  if (!class(sim$fireSense_IgnitionPredicted) == "RasterLayer"){
    IgnitionRas <- raster(sim$flammableRTM)
    IgnitionRas[fireSense_IgnitionCovariates$pixelID] <- sim$fireSense_IgnitionPredicted
    sim$fireSense_IgnitionPredicted <- IgnitionRas
  }

  invisible(sim)
}

IgnitionPredictSave <- function(sim) {
  moduleName <- current(sim)$moduleName
  timeUnit <- timeunit(sim)
  currentTime <- time(sim, timeUnit)

  raster::writeRaster(
    sim$fireSense_IgnitionPredicted,
    filename = file.path(paths(sim)$out, paste0("fireSense_IgnitionPredicted_", timeUnit, currentTime, ".tif"))
  )

  invisible(sim)
}

.inputObjects <- function(sim) {
  # cacheTags <- c(currentModule(sim), "function:.inputObjects")
  # dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  # message(currentModule(sim), ": using dataPath '", dPath, "'.")

  if (!suppliedElsewhere("rescaleFactor", sim)) {
    sim$rescaleFactor <- (250 / 10000)^2
  }

  return(invisible(sim))
}
