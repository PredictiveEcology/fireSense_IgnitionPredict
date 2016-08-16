library(SpaDES)

## data.frame
  # Poisson
    # mySim <- simInit(
    #   times = list(start = 1, end = 1, timeunit = "year"),
    #   modules = list("fireSense_FrequencyPredict"),
    #   paths = list(modulePath = " # replace with empty string instead"),
    #   params = list(fireSense_FrequencyPredict = list(data = "dataFireSense_Frequency")),
    #   inputs = data.frame(
    #     files = c("Z:/fireSense_FrequencyFitted_P.rds", "Z:/dataFireSense_Frequency.rds"),
    #     functions = c("readRDS", "readRDS"),
    #     package = c("base", "base"),
    #     objectName = c("fireSense_FrequencyFitted", NA),
    #     stringsAsFactors = FALSE)
    # )

  ## Negative binomial
    # mySim <- simInit(
    #   times = list(start = 1, end = 1, timeunit = "year"),
    #   modules = list("fireSense_FrequencyPredict"),
    #   paths = list(modulePath = " # replace with empty string instead"),
    #   params = list(fireSense_FrequencyPredict = list(data = "dataFireSense_Frequency")),
    #   inputs = data.frame(
    #     files = c("Z:/fireSense_FrequencyFitted_NB.rds", "Z:/dataFireSense_Frequency.rds"),
    #     functions = c("readRDS", "readRDS"),
    #     package = c("base", "base"),
    #     objectName = c("fireSense_FrequencyFitted", NA),
    #     stringsAsFactors = FALSE)
    # )

## RasterLayer
  ## 1 var
    # mySim <- simInit(
    #   times = list(start = 1, end = 1, timeunit = "year"),
    #   modules = list("fireSense_FrequencyPredict"),
    #   paths = list(modulePath = " # replace with empty string instead"),
    #   params = list(fireSense_FrequencyPredict = list(mapping = c(MDC_JUL = "MDC_JUN"))), ## An example of mapping
    #   inputs = data.frame(
    #     files = c("Z:/fireSense_FrequencyFitted_P.rds", "Z:/MDC_JUN.tif"),
    #     functions = c("readRDS", "raster"),
    #     package = c("base", "raster"),
    #     objectName = c("fireSense_FrequencyFitted", NA),
    #     stringsAsFactors = FALSE)
    # )

  ## 6 var
    mySim <- simInit(
      times = list(start = 1, end = 1, timeunit = "year"),
      modules = list("fireSense_FrequencyPredict"),
      paths = list(modulePath = " # replace with empty string instead"),
      inputs = data.frame(
        files = c("Z:/fireSense_FrequencyFitted_NB.rds", "Z:/MDC_JUN.tif", "Z:/HW.tif", "Z:/HW.tif", "Z:/HW.tif", "Z:/HW.tif"),
        functions = c("readRDS", "raster", "raster", "raster", "raster", "raster"),
        package = c("base", "raster", "raster", "raster", "raster", "raster"),
        objectName = c("fireSense_FrequencyFitted", "MDC_JUL", "HW", "CN", "D", "O"),
        stringsAsFactors = FALSE)
    )

spades(mySim, debug = FALSE)

# str(mySim$fireSense_FrequencyPredict)
# x11(); Plot(mySim$fireSense_FrequencyPredict)