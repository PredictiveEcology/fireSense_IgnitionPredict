---
title: "fireSense_IgnitionPredict Manual"
subtitle: "v.1.0.0.9000"
date: "Last updated: 2026-09-21"
output:
  bookdown::html_document2:
    toc: true
    toc_float: true
    theme: sandstone
    number_sections: false
    df_print: paged
    keep_md: yes
editor_options:
  chunk_output_type: console
bibliography: citations/references_fireSense_IgnitionPredict.bib
link-citations: true
always_allow_html: true
---

# fireSense_IgnitionPredict Module

<!-- the following are text references used in captions for LaTeX compatibility -->
(ref:fireSense-IgnitionPredict) *fireSense_IgnitionPredict*



[![made-with-Markdown](figures/markdownBadge.png)](https://commonmark.org)

<!-- if knitting to pdf remember to add the pandoc_args: ["--extract-media", "."] option to yml in order to get the badge images -->

#### Authors:

Eliot McIntire <eliot.mcintire@nrcan-rncan.gc.ca> [aut, cre], Ian Eddy <ian.eddy@nrcan-rncan.gc.ca> [aut], Jean Marchal <jean.d.marchal@gmail.com> [aut], Alex M Chubaty <achubaty@for-cast.ca> [ctb]
<!-- ideally separate authors with new lines, '\n' not working -->

## Module Overview

### Module summary

Predicts fire frequency or rates of fire counts using a model fitted with the `fireSense_IgnitionFit` module.
Use them to feed the ignition component of a landscape fire model (e.g fireSense [@Marchal:2017a; @Marchal:2017b; @Marchal:2019]).

### Module inputs and parameters

Table \@ref(tab:moduleInputs-fireSense-IgnitionPredict) shows the full list of module inputs.

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleInputs-fireSense-IgnitionPredict)(\#tab:moduleInputs-fireSense-IgnitionPredict)List of (ref:fireSense-IgnitionPredict) input objects and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
   <th style="text-align:left;"> sourceURL </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> fireSense_EscapeFitted </td>
   <td style="text-align:left;"> fireSense_EscapeFit </td>
   <td style="text-align:left;"> An object of class `fireSense_EscapeFit` created with the `fireSense_IgnitionFit` module. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_IgnitionFitted </td>
   <td style="text-align:left;"> fireSense_IgnitionFit </td>
   <td style="text-align:left;"> An object of class `fireSense_IgnitionFit` created with the `fireSense_IgnitionFit` module. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_igAndEscapePred_Covariates </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> A `data.table` with prediction variables and a column named 'pixelID' </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> flammableRTM </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> RTM without ice/rocks/urban/water. Flammable map with 0 and 1. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
</tbody>
</table>

Summary of user-visible parameters (Table \@ref(tab:moduleParams-fireSense-IgnitionPredict))


<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleParams-fireSense-IgnitionPredict)(\#tab:moduleParams-fireSense-IgnitionPredict)List of (ref:fireSense-IgnitionPredict) parameters and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> paramName </th>
   <th style="text-align:left;"> paramClass </th>
   <th style="text-align:left;"> default </th>
   <th style="text-align:left;"> min </th>
   <th style="text-align:left;"> max </th>
   <th style="text-align:left;"> paramDesc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> modelAlgorithm </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> xgboost </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Can be `xgboost`, `glmmtmb`, `glm.nb`, `glmmadaptive`, `glm`; only `xgboost` is supported currently </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rescaleVars </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Attempt to rescale variables? If `rescalers` is defined, use it to rescale variables as `var / rescalers['var']`. Otherwise, `scale()` will be used to rescale variables to `[0,1]`, if they are not already within this range. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .runInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> when to start this module? By default, the start time of the simulation. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .runInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> optional. Interval between two runs of this moduleexpressed in units of simulation time. By default, 1 year. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .saveInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> optional. When to start saving output to a file. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .saveInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> optional. Interval between save events. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useCache </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant </td>
  </tr>
</tbody>
</table>

### Events
<!-- TODO -->
- Module initialization;
- Make predictions;

### Plotting
<!-- TODO -->
Write what is plotted.

### Saving
<!-- TODO -->
Write what is saved.

### Module outputs

Description of the module outputs (Table \@ref(tab:moduleOutputs-fireSense-IgnitionPredict)).

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleOutputs-fireSense-IgnitionPredict)(\#tab:moduleOutputs-fireSense-IgnitionPredict)List of (ref:fireSense-IgnitionPredict) outputs and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> fireSense_IgAndEscapeProbRas </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> a raster layer of the annual ignition and escape probabilities </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ignitionsAndEscapes </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> A data.table containing pixelID (referencing flammableRTM), ignitions, escapes, and their associated probabilities </td>
  </tr>
</tbody>
</table>

### Links to other modules

<!-- TODO: add links to other fireSense modules -->
Predictions made with this module can be used to feed the ignition component of a landscape fire model (e.g fireSense).

### Getting help

- <https://github.com/PredictiveEcology/fireSense_IgnitionPredict/issues>

## References

<!-- autogenerated from bibligraphy -->
