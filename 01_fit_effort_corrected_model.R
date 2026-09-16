
# ==============================================================================
# Effort-corrected thinned LGCP model for Mediterranean shortfin mako
#
# Associated study:
# "Integrating opportunistic records to model spatio-temporal abundance of a
# Critically Endangered pelagic shark in the Mediterranean Sea"
#
# This script fits the effort-corrected model used in the manuscript.
# All covariates and observation-effort surfaces are supplied as processed
# model inputs. 
# ==============================================================================


## ---- Packages ----------------------------------------------------------------

library(INLA)
library(inlabru)
library(fmesher)
library(sf)
library(sp)
library(terra)


## ---- Load processed model inputs ---------------------------------------------

# The repository should contain: model_inputs.RData
#
# Required objects:
#   mk.MED          MEDLEM occurrences used for model fitting
#   mk.SP           sharkPulse occurrences used for model fitting
#   mesh_medit      spatial mesh
#   triBarrier      indices of mesh triangles classified as terrestrial barriers
#   aoi.sf          Mediterranean study domain
#   depth.r         depth covariate projected to the model domain
#   bpi.r           Bathymetric Position Index covariate
#   sstc_s.mesh     spatial SST component
#   sstc_res.mesh   annual residual spatio-temporal SST component
#   sstc_t          temporal SST component (2017--2022)
#   eff.SAR.r       annual MEDLEM observation-effort proxy
#   eff.SP.r        annual sharkPulse observation-effort proxy

load("data/model_inputs.RData")


## ---- Basic input checks -------------------------------------------------------

required_objects <- c(
  "mk.MED", "mk.SP", "mesh_medit", "triBarrier", "aoi.sf",
  "depth.r", "bpi.r", "sstc_s.mesh", "sstc_res.mesh", "sstc_t",
  "eff.SAR.r", "eff.SP.r"
)


# The fitted model uses six annual layers corresponding to 2017--2022.
stopifnot(
  terra::nlyr(sstc_res.mesh) >= 6,
  terra::nlyr(eff.SAR.r) >= 6,
  terra::nlyr(eff.SP.r) >= 6,
  length(sstc_t) == 6
)


## ---- Barrier SPDE model -------------------------------------------------------

barrier.model <- INLA::inla.barrier.pcmatern(
  mesh_medit,
  barrier.triangles = triBarrier,
  prior.range = c(800, 0.7),
  prior.sigma = c(1, 0.1)
)


## ---- Temporal component -------------------------------------------------------

# One temporal mesh node per modelled year (2017-2022).
knots <- seq(1, 6, 1)

mesh1D <- fmesher::fm_mesh_1d(
  knots,
  boundary = "free"
)

spline.t <- INLA::inla.spde2.pcmatern(
  mesh1D,
  prior.range = c(5, 0.8),
  prior.sigma = c(1, 0.1)
)


## ---- Source-specific detection functions -------------------------------------

qexppnorm <- function(x, rate) {
  qexp(
    pnorm(x, lower.tail = FALSE, log.p = TRUE),
    rate = rate,
    log.p = TRUE,
    lower.tail = FALSE
  )
}

# MEDLEM: SAR-derived apparent fishing activity.
log_detect_year_MED <- function(xy, t1, t2, Year1) {

  dens <- inlabru::eval_spatial(
    data = eff.SAR.r,
    where = xy,
    layer = Year1
  )

  dens[is.na(dens)] <- min(dens, na.rm = TRUE)

  pnorm(
    dens / qexppnorm(t1, rate = 0.5) + t2,
    log.p = TRUE
  )
}

# sharkPulse: species-independent sharkPulse observation intensity.
log_detect_year_SP <- function(xy, t1, t2, Year1) {

  dens <- inlabru::eval_spatial(
    data = eff.SP.r,
    where = xy,
    layer = Year1
  )

  dens[is.na(dens)] <- min(dens, na.rm = TRUE)

  pnorm(
    dens / qexppnorm(t1, rate = 0.5) + t2,
    log.p = TRUE
  )
}


## ---- Model components ---------------------------------------------------------

# Environmental covariates are fitted as linear fixed effects.
# Year_cov is the separate latent temporal component.

cmp <- ~ 0 +
  intercept1(
    1,
    mean.linear = 0,
    prec.linear = 2
  ) +
  SPDE(
    geometry,
    model = barrier.model,
    mapper = bru_mapper(mesh_medit)
  ) +
  depth_s(
    depth.r$depth,
    mean.linear = 0.5,
    prec.linear = 2
  ) +
  sstc_s(
    sstc_s.mesh$sstc_s.r,
    mean.linear = 0,
    prec.linear = 5
  ) +
  sstc_res(
    eval_spatial(sstc_res.mesh, .data., layer = Year1),
    model = "linear",
    mean.linear = 1,
    prec.linear = 2
  ) +
  sstc_temp(
    sstc_t,
    model = "linear",
    mean.linear = 1,
    prec.linear = 5
  ) +
  bpi(
    bpi.r$bpi,
    mean.linear = 1,
    prec.linear = 5
  ) +
  Year_cov(
    Year1,
    model = spline.t
  ) +
  effort_scale_MED(
    1,
    prec.linear = 1
  ) +
  effort_shift_MED(
    1,
    prec.linear = 5
  ) +
  effort_scale_SP(
    1,
    prec.linear = 1
  ) +
  effort_shift_SP(
    1,
    prec.linear = 5
  )


## ---- Source-specific likelihoods ---------------------------------------------

formula_MED <- geometry + Year1 ~
  intercept1 +
  SPDE +
  depth_s +
  sstc_s +
  sstc_res +
  sstc_temp +
  bpi +
  Year_cov +
  log_detect_year_MED(
    geometry,
    effort_scale_MED,
    effort_shift_MED,
    Year1
  )

formula_SP <- geometry + Year1 ~
  intercept1 +
  SPDE +
  depth_s +
  sstc_s +
  sstc_res +
  sstc_temp +
  bpi +
  Year_cov +
  log_detect_year_SP(
    geometry,
    effort_scale_SP,
    effort_shift_SP,
    Year1
  )

# Preserve the sampler/domain construction used in the analysis.
boundary <- as(aoi.sf, "Spatial")
boundary_domain <- sf::st_as_sf(boundary)

lik_MED <- bru_obs(
  "cp",
  formula = formula_MED,
  samplers = boundary_domain,
  domain = list(
    geometry = mesh_medit,
    Year1 = 1:6
  ),
  data = mk.MED,
  used = bru_used(
    effect_exclude = c(
      "effort_scale_SP",
      "effort_shift_SP"
    )
  )
)

lik_SP <- bru_obs(
  "cp",
  formula = formula_SP,
  samplers = boundary_domain,
  domain = list(
    geometry = mesh_medit,
    Year1 = 1:6
  ),
  data = mk.SP,
  used = bru_used(
    effect_exclude = c(
      "effort_scale_MED",
      "effort_shift_MED"
    )
  )
)


## ---- Fit model ----------------------------------------------------------------

fit_effort_corrected <- bru(
  components = cmp,
  lik_MED,
  lik_SP,
  options = list(
    verbose = FALSE,
    inla.mode = "experimental",
    bru_max_iter = 35,
    control.inla = list(
      int.strategy = "eb"
    )
  )
)


## ---- Model summary ------------------------------------------------------------

summary(fit_effort_corrected)


## ---- Optional: save fitted model ---------------------------------------------

# Uncomment if a serialized fitted-model object is required.
#
# dir.create("output", showWarnings = FALSE)
# save(
#   fit_effort_corrected,
#   file = "output/fit_effort_corrected.RData"
# )
