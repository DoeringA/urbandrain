[![R-CMD-check](https://github.com/mrustl/urbandrain/workflows/R-CMD-check/badge.svg)](https://github.com/mrustl/urbandrain/actions?query=workflow%3AR-CMD-check)
[![pkgdown](https://github.com/mrustl/urbandrain/workflows/pkgdown/badge.svg)](https://github.com/mrustl/urbandrain/actions?query=workflow%3Apkgdown)

This R package provides functions to automatically generate a stormwater drainage network. The aim is to keep existing urban structures and rebuild the pipe system as close to reality as possible. Since public drainage pipes are arranged in the road cross-section as gravity driven open channels, we setup the drainage network below the existing streets and define flow directions following the surface topology. Required input data are publicly available. The result is a SWMM input file including a complete drainage network model with discharging surfaces, junctions, conduits and outfalls. 


## Installation

You can install the package from github with:

```{r gh-installation, eval = FALSE}
# install.packages("remotes")
remotes::install_github("mrustl/urbandrain")
```

## Workflows

Checkout the [workflows](articles/index.html), what can be done 
with this R package.

For further details on the input data please view vignette: [How to extract and preprocess input data from open source](articles/How_to_extract_and_preprocess_input_data_from_open_source.html). For further information on the usage of the functions provided by the package please view vignette: [How to built a drainage network model using urbandrain](articles/How_to_built_a_drainage_network_model_using_urbandrain.html)
