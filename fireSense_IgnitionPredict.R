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
  parameters = bindrows(
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
                              "class(fireSense_IgnitionAndEscapeCovariates) == 'data.table'")),
    expectsInput(objectName = "rescaleFactor", objectClass = "numeric",
                 desc = paste("rescale predicted rates of fire counts at any given temporal and spatial",
                              "resolutions by a factor `rescaleFactor = new_res / old_res`.",
                              "`rescaleFactor` is the ratio between the data aggregation scale used",
                              "for model fitting and the scale at which predictions are to be made",
                              "If not provided, defaults to (250 / 10000)^2"))
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
  return(invisible(sim))
}

IgnitionPredictRun <- function(sim) {

  moduleName <- currentModule(sim)

  if (class(sim$fireSense_IgnitionAndEscapeCovariates) == "RasterStack") {
    fireSense_IgnitionCovariates <- sim$fireSense_IgnitionAndEscapeCovariates
  } else {
    fireSense_IgnitionCovariates <- copy(setDT(sim$fireSense_IgnitionAndEscapeCovariates))

    ## check
    if (is.null(sim$flammableRTM)) {
      stop("'fireSense_IgnitionAndEscapeCovariates' is a table. Please supply 'flammableRTM'")
    }
  }

  #TODO: IE wrote this - please review it
  ## TODO: Ceres added more - please review AGAIN
  if (!is.null(sim$fireSense_IgnitionFitted$rescales)) {
    ## make data frame if raster stack has been applied - easier and less repetitive coding
    if (class(fireSense_IgnitionCovariates) == "RasterStack") {
      fireSense_IgnitionCovariatesSc <- raster::as.data.frame(fireSense_IgnitionCovariates)
      fireSense_IgnitionCovariatesSc <- copy(setDT(fireSense_IgnitionCovariatesSc))
    } else {
      fireSense_IgnitionCovariatesSc <- copy(fireSense_IgnitionCovariates)
    }

    for (cn in names(sim$fireSense_IgnitionFitted$rescales)) {
      rescaleFun <- sim$fireSense_IgnitionFitted$rescales[[cn]]
      if (grepl("rescale", rescaleFun)) {
        set(
          fireSense_IgnitionCovariatesSc, NULL, cn,
          rescaleKnown2(x = fireSense_IgnitionCovariatesSc[[cn]],
                        minNew = 0,
                        maxNew = 1,
                        minOrig = min(sim$covMinMax_ignition[[cn]]),
                        maxOrig = max(sim$covMinMax_ignition[[cn]]))
        )
      } else {
        op <- eval(parse(text = rescaleFun),
                   env = fireSense_IgnitionCovariatesSc)
        fireSense_IgnitionCovariatesSc[, eval(cn) := op]   ## this will fail if raster stack is provided!
      }
    }

    if (class(fireSense_IgnitionCovariates) == "RasterStack") {
      fireSense_IgnitionCovariates <- sapply(names(fireSense_IgnitionCovariates), FUN = function(var, stk, DT) {
        ras <- stk[[var]]
        ras[] <- DT[[var]]
        ras
      }, stk = fireSense_IgnitionCovariates, DT = fireSense_IgnitionCovariatesSc,
      simplify = FALSE, USE.NAMES = TRUE)
      fireSense_IgnitionCovariates <- stack(fireSense_IgnitionCovariates)
    } else {
      fireSense_IgnitionCovariates <- fireSense_IgnitionCovariatesSc
    }
    rm(fireSense_IgnitionCovariatesSc)
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
    if (class(fireSense_IgnitionCovariates) == "RasterStack") {
      ## make a template mask-type raster with 1s where there are no NAs
      rasTemp <- sum(!is.na(fireSense_IgnitionCovariates))
      rasTemp[rasTemp[] == 0] <- NA
      rasTemp[!is.na(rasTemp)] <- 1

      ## make new raster stack with knot values
      rasStk <- sapply(names(sim$fireSense_IgnitionFitted$knots), FUN = function(knotVar, ras, knots) {
        ras[!is.na(ras)] <- rep(knots[[knotVar]], sum(!is.na(ras[])))
        ras
      }, ras = rasTemp, knots = sim$fireSense_IgnitionFitted$knots,
      simplify = FALSE, USE.NAMES = TRUE) %>%
        raster::stack(.)

      fireSense_IgnitionCovariates <- stack(fireSense_IgnitionCovariates, rasStk)
      rm(rasStk, rasTemp); gc()
    } else {
      for (knot in names(sim$fireSense_IgnitionFitted$knots)) {
        fireSense_IgnitionCovariates[, eval(knot) := sim$fireSense_IgnitionFitted$knots[knot]]
      }
    }

    kNames <- names(sim$fireSense_IgnitionFitted$knots)
    # allxy <- allxy[!allxy %in% kNames]   ## Ceres: not needed, but more testing may be necessary
  }

  if (class(fireSense_IgnitionCovariates) == "RasterStack") {
    list2env(setNames(unstack(fireSense_IgnitionCovariates), names(fireSense_IgnitionCovariates)),
             env = mod_env)
  } else {
    list2env(fireSense_IgnitionCovariates, env = mod_env)
  }

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
    ## knots are not supplied by the user. so exclude them from these checks
    allxyNoK <- allxy[!allxy %in% kNames]
    missing <- !allxyNoK %in% ls(mod_env, all.names = TRUE)

    if (s <- sum(missing)) {
      stop(
        moduleName, "> '", allxyNoK[missing][1L], "'",
        if (s > 1) paste0(" (and ", s - 1L, " other", if (s > 2) "s", ")"),
        " not found in data objects."
      )
    }

    badClass <- unlist(lapply(allxyNoK, function(x) is.vector(mod_env[[x]]) || is(mod_env[[x]], "RasterLayer")))

    if (any(badClass)) {
      stop(moduleName, "> Data objects of class 'data.frame' cannot be mixed with objects of other classes.")
    } else {
      stop(moduleName, "> '", paste(allxyNoK[which(!badClass)], collapse = "', '"),
           "' does not match a data.frame's column, a RasterLayer or a layer from a RasterStack or RasterBrick.")
    }
  }

  ## output raster AND table
  ## TODO: add more assertions - do pixelIDs correspond to raster cells? are there NAs to be accounted for
  ## compareRaster(flammableRTM, ...)
  if (class(sim$fireSense_IgnitionPredicted) != "RasterLayer") {
    sim$fireSense_IgnitionPredictedVec <- sim$fireSense_IgnitionPredicted
    ## checks
    if (!"pixelID" %in% names(fireSense_IgnitionCovariates)) {
      stop("fireSense_IgnitionAndEscapeCovariates must have a 'pixelID' column")
    }

    IgnitionRas <- raster(sim$flammableRTM)
    IgnitionRas[fireSense_IgnitionCovariates$pixelID] <- sim$fireSense_IgnitionPredictedVec
    ## "re-write" object as a raster
    sim$fireSense_IgnitionPredicted <- IgnitionRas
  } else {
    fireSense_IgnitionPredictedVec <- raster::as.data.frame(sim$fireSense_IgnitionPredicted, na.rm = TRUE)
    sim$fireSense_IgnitionPredictedVec <- fireSense_IgnitionPredictedVec[["layer"]]
    names(sim$fireSense_IgnitionPredictedVec) <- rownames(fireSense_IgnitionPredictedVec)
  }

  return(invisible(sim))
}

IgnitionPredictSave <- function(sim) {
  moduleName <- current(sim)$moduleName
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

  if (!suppliedElsewhere("rescaleFactor", sim)) {
    sim$rescaleFactor <- (250 / 10000)^2
  }

  return(invisible(sim))
}
