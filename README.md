# Mediterranean shortfin mako spatio-temporal model

This repository contains the processed model inputs and R code required to fit
the spatio-temporal models described in the associated manuscript.

The analysis integrates shortfin mako (*Isurus oxyrinchus*) occurrence records
from the Mediterranean Sea and accounts for source-specific variation in the
observation process.

## Repository contents

- `01_fit_effort_corrected_model.R`  
Fits the main effort-corrected thinned log-Gaussian Cox process model, computes the annual relative abundance index, and estimates the relative change between 2017 and 2022 from joint posterior draws.

- `02_fit_no_effort_model.R`  
  Fits the corresponding model without observation-effort correction, used as
  a sensitivity analysis.

- `model_inputs_1.RData`  
  Contains environmental and spatial raster inputs (`depth.r`, `bpi.r`, `sstc_s.mesh`, `sstc_res.mesh`).

- `model_inputs_2.RData`  
  Contains source-specific observation-effort raster surfaces (`eff.SAR.r`, `eff.SP.r`).

- `model_inputs_3.RData`  
  Contains additional spatial and model inputs (`mk.MED`, `mk.SP`, `mesh_medit`, `triBarrier`, `aoi.sf`, `sstc_t`).

The model inputs are provided in three separate files to accommodate GitHub file upload size limitations.

## Software

The models were fitted in R using the `INLA` and `inlabru` packages.

## Reproducibility
All three `.RData` files are required to reproduce the analyses.

The input datasets can be loaded into R using:

```r
load("model_inputs_1.RData")
load("model_inputs_2.RData")
load("model_inputs_3.RData")
```

Raster objects were stored using `terra::wrap()` to facilitate serialization. They can be restored using `terra::unwrap()` before model fitting.

To fit the main model:

```r
source("01_fit_effort_corrected_model.R")
```

To fit the model without observation-effort correction:
```r
source("02_fit_no_effort_model.R")
```








This repository has been prepared for anonymous peer review.
