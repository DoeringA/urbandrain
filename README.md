[![R-CMD-check](https://github.com/mrustl/urbandrain/workflows/R-CMD-check/badge.svg)](https://github.com/mrustl/urbandrain/actions?query=workflow%3AR-CMD-check)
[![pkgdown](https://github.com/mrustl/urbandrain/workflows/pkgdown/badge.svg)](https://github.com/mrustl/urbandrain/actions?query=workflow%3Apkgdown)


# urbandrain

This R package provides functions to automatically generate a stormwater
drainage network. The aim is to keep existing urban structures and
rebuild the pipe system as close to reality as possible. Since public
drainage pipes are arranged in the road cross-section as gravity driven
open channels, we setup the drainage network below the existing streets
and define flow directions following the surface topology. Required
input data are publicly available. The result is a SWMM input file
including a complete drainage network model with discharging surfaces,
junctions, conduits and outfalls.


## Documentation 

For further details checkout the documentation website:

Release: [https://mrustl.de/urbandrain](https://mrustl.de/urbandrain)

Development: [https://mrustl.de/urbandrain/dev](https://mrustl.de/urbandrain/dev)