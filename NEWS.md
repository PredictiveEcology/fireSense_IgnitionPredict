# fireSense_IgnitionPredict (development version)

- Several fitted ELFs in one study area: with `fireSense_IgnitionFittedList` and `fireSense_EscapeFittedList` (one fit per ELF, named by `ELFind`) and `rasterToMatchLargeELF`, each ELF's models predict the coarse pixels of that ELF, with its own `scaleData`. The Poisson and binomial draws are still made once over all pixels, so one ELF gives exactly the previous result.

# fireSense_IgnitionPredict (development version)

## Bug fixes

- Covariates are now standardized with the center and scale stored by the fit (`fireSense_IgnitionFitted$scaleData`, `fireSense_EscapeFitted$scaleData`). They were standardized with each year's own mean and sd, so every year looked average to the models. With `rescaleVars = TRUE`, a fitted object without `scaleData` is now an error.

## Cleanup

- The `save` event now does nothing and says so. It used to fail, because it wrote `sim$fireSense_IgnitionPredicted`, which the module no longer creates. `IgnitionPredictSave()` is removed.

# fireSense_IgnitionPredict 1.0.0

First release from `development` since `master` was last updated (2021-03-11). Full history: https://github.com/PredictiveEcology/fireSense_IgnitionPredict/compare/5a242fa...v1.0.0

## Breaking changes

- Removed input `dataFireSense_IgnitionPredict` (data.frame, RasterLayer, RasterStack).
- Removed output `fireSense_IgnitionPredicted`.
- Removed parameters: `data`, `mapping`, `modelObjName`, `rescaleFactor`.

## New features

- New inputs: `fireSense_EscapeFitted`, `fireSense_igAndEscapePred_Covariates`, `flammableRTM`.
- New outputs: `fireSense_IgAndEscapeProbRas`, `ignitionsAndEscapes`.
- New parameters: `modelAlgorithm`, `rescaleVars`.

## Dependencies

- No longer depends on `raster`.
- Now depends on `fireSenseUtils`, `terra`.

## Testing

- testthat suite and CI (`testthat-module`), including a snapshot of the module's inputs, outputs and parameters in `tests/testthat/test-metadata.R`.
