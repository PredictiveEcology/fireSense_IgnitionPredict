inputsNoIgnition <- function() toyIgInputs(list(constant(0)), list(constant(0)))

test_that("run repeats every .runInterval from .runInitialTime", {
  sim <- toyIgRun(inputsNoIgnition(), times = list(start = 2001, end = 2003))
  expect_equal(evOf(SpaDES.core::completed(sim), "run")$eventTime, 2001:2003)
  expect_equal(evOf(SpaDES.core::events(sim), "run")$eventTime, 2004)
  expect_equal(evOf(SpaDES.core::completed(sim), "init")$eventTime, 2001)
  ## one call to the ignition model per run event
  expect_length(toyCalls("ign1"), 3L)
})

test_that(".runInitialTime and .runInterval are honoured", {
  sim <- toyIgRun(inputsNoIgnition(), params = list(.runInitialTime = 2002, .runInterval = 2),
                  times = list(start = 2001, end = 2005))
  expect_equal(evOf(SpaDES.core::completed(sim), "run")$eventTime, c(2002, 2004))
  expect_equal(evOf(SpaDES.core::events(sim), "run")$eventTime, 2006)
})

test_that(".runInterval = NA predicts once", {
  sim <- toyIgRun(inputsNoIgnition(), params = list(.runInterval = NA_real_),
                  times = list(start = 2001, end = 2003))
  expect_equal(evOf(SpaDES.core::completed(sim), "run")$eventTime, 2001)
  expect_identical(nrow(evOf(SpaDES.core::events(sim), "run")), 0L)
})

test_that("no save event by default; .saveInitialTime schedules one, last in its year", {
  sim <- toyIgRun(inputsNoIgnition(), times = list(start = 2001, end = 2002))
  expect_identical(nrow(evOf(SpaDES.core::completed(sim), "save")), 0L)
  expect_identical(nrow(evOf(SpaDES.core::events(sim), "save")), 0L)
  ## scheduled beyond end(sim) so that the (broken) save event itself never runs
  sim <- toyIgRun(inputsNoIgnition(), params = list(.saveInitialTime = 2010),
                  times = list(start = 2001, end = 2001))
  sv <- evOf(SpaDES.core::events(sim), "save")
  expect_equal(sv$eventTime, 2010)
  expect_equal(sv$eventPriority, SpaDES.core::.last())
})

test_that("nothing is predicted before the first run event", {
  sim <- toyIgRun(inputsNoIgnition(), params = list(.runInitialTime = 2005),
                  times = list(start = 2001, end = 2002))
  expect_null(sim$ignitionsAndEscapes)
  expect_null(sim$fireSense_IgAndEscapeProbRas)
})
