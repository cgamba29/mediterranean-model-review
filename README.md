# Mediterranean shortfin mako spatio-temporal model

This repository contains the processed model inputs and R code required to fit
the spatio-temporal models described in the associated manuscript.

The analysis integrates shortfin mako (*Isurus oxyrinchus*) occurrence records
from the Mediterranean Sea and accounts for source-specific variation in the
observation process.

## Repository contents

- `01_fit_effort_corrected_model.R`  
  Fits the main effort-corrected thinned log-Gaussian Cox process model.

- `02_fit_no_effort_model.R`  
  Fits the corresponding model without observation-effort correction, used as
  a sensitivity analysis.

- `model_inputs.RData`  
  Contains the processed occurrence data, environmental covariates, spatial
  mesh and barrier objects, and source-specific observation-effort surfaces
  required to fit the models.

## Software

The models were fitted in R using the `INLA` and `inlabru` packages.

## Reproducibility

To fit the main model:

```r
source("01_fit_effort_corrected_model.R")
```

To fit the model without observation-effort correction:
```r
source("02_fit_no_effort_model.R")
```
Both scripts load the processed inputs from:

```r
model_inputs.RData
```





This repository has been prepared for anonymous peer review.
