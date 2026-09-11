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
