cchunts -- Cross-cultural hunting life history analysis R package
==========

This R package contains the data and code necessary to repeat and extend the Koster et al life history analysis of human foraging.

## Installation

Install from github directly:
```
library(devtools)
devtools::install_github("rmcelreath/cchunts")
library(cchunts)
```

## Data

The data sets made available in the package are documented:
```
?cchunts
```
See also ``?make_joint``.


## Simulation and validation

The package provides code to simulate foraging trips, both for exploring study design and validating the analytical code. The simulation approach is documented:
```
?sim_foragers
```
See also ``?sim_multi``.


## Data analysis

To repeat the data analysis from the paper, see the examples at:
```
?run_cchunt_model
```

