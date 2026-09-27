# inlaws 0.0.0.9000

* Added `cmpls()` with separate log-mean and log-dispersion predictors for
  Conway-Maxwell-Poisson counts, REML and NCV support in `gam()`, and
  simulation, CDF, and quantile callbacks.

* Added `pigls()` for Poisson inverse Gaussian counts with separate log-mean
  and log-dispersion predictors, REML and NCV fitting in `gam()`, and
  simulation, CDF, and quantile callbacks.

* Added `pig()` for Poisson inverse Gaussian count models with fixed or
  estimated global dispersion, REML/ML/NCV support, ordinary and discrete
  `bam()` fitting, and simulation, CDF, and quantile callbacks.

* Added `cgammals()` for censored gamma models with separate log-mean and
  log-dispersion predictors, REML and NCV fitting, and distribution callbacks.

* Renamed `gp()` to `gpoisson()` to distinguish it from `gpd()`. Update
  existing calls to the new name; the old constructor is no longer exported.

* Renamed `ocat_link()` to `cumulative_link()`. Update existing calls to the
  new name; the old constructor is no longer exported.

* Migrated 14 additional mgcv family constructors, their distribution callbacks,
  mathematical notes, and regression tests from mgcvUtils development branches.
* Included the family-specific `predict_ordbeta()` helper.

* Replaced the statmod dependency with internally generated and cached
  Gauss-Legendre quadrature rules for the censored gamma and log-normal families.
