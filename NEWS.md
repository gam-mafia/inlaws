# inlaws 0.0.0.9000

## New features

* Added `betals()` for beta responses with separate logit-mean and
  log-precision smooth predictors, REML and NCV fitting in `gam()`,
  simulation, CDF, quantile, and diagnostic callbacks.

- Add `log_logistic()` and `log_logisticls()` for arithmetic mean and
  log-response logistic scale, using log/logit links. Includes distribution
  callbacks, analytic fourth-order derivatives, and tested REML/NCV routes;
  the single-predictor family also supports ML and discrete `bam()`.

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

* Migrated 14 additional mgcv family constructors, their distribution callbacks,
  mathematical notes, and regression tests from mgcvUtils development branches.
  Native response predictions for `zinb()` and `zanb()` include the complete
  response mean, variance, and zero probability after the original three
  component parameters, with full-covariance delta-method standard errors.
* Included the family-specific `predict_ordbeta()` helper.

* All family constructors now supply native `gratia_diagnostics()` factories
  without requiring gratia. They bind fitted parameters to response moments,
  count left limits, ordinal and ordered-beta masses, censoring intervals, and
  beta-binomial count functions. Multivariate families provide component moments
  and marginal CDFs with joint simulation. Heavy-tail families report undefined
  moments explicitly.
* Censored constructors accept a known `observe` function for simulated
  diagnostic references. `censoring_process()` implements shared scalar limits;
  interval PIT and quantile checks remain available when the rule is unknown.

* Added `cgen_gammalss()` for censored generalized gamma mean, scale and
  transformed-shape models, with exact, left-, right- and interval-censored
  responses, REML and NCV/QNCV in `gam()`, distribution callbacks, and native
  gratia diagnostics.

* Added `zoibetals()` with configurable zero, one, or both endpoint
  masses, separate interior mean, precision and membership predictors, EFS
  REML fitting, named component/moment predictions, distribution callbacks,
  and native randomized quantile diagnostics.

* Added `gen_gammalss()` for generalized gamma mean, scale and transformed
  shape models. Log links model the arithmetic mean, Prentice scale, and
  `r = 1 + sigma * Q`, ensuring a finite mean. Includes REML, NCV and QNCV
  in `gam()`, distribution callbacks, and native gratia diagnostics.
  Constant scale and shape predictors are supported with both REML optimizers.

## User visible changes

* `ordbeta()` response predictions now return a five-column matrix ordered
  `mu`, `p0`, `pI`, `p1`, `mean`, following `zoibetals()` conventions. These
  contain the conditional interior mean, zero/interior/one probabilities, and
  complete response mean. Code expecting the previous mean vector should use
  column 5 or `predict_ordbeta(fit)`. The helper selects standard predictions
  and restores column names, with probabilities now named `p0`, `pI`, `p1`.
  Delta-method standard errors condition on fitted precision and active
  cutpoints, including with `unconditional = TRUE`. Link predictions and
  the interior means stored in `fitted.values` retain their representation.

* Renamed `gp()` to `gpoisson()` to distinguish it from `gpd()`. Update
  existing calls to the new name; the old constructor is no longer exported.

* Renamed `ocat_link()` to `cumulative_link()`. Update existing calls to the
  new name; the old constructor is no longer exported.

* Replaced the statmod dependency with internally generated and cached
  Gauss-Legendre quadrature rules for the censored gamma and log-normal families.
