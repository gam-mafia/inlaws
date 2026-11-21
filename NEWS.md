# inlaws 0.0.0.9000

## New features

- Add `log_logistic()` and `log_logisticls()` for arithmetic mean and
  log-response logistic scale, using log/logit links and analytic fourth-order
  derivatives. Both support REML, NCV, and QNCV in `gam()`. The single-predictor
  family also supports ML and ordinary and discrete `bam()`.

- Add `cmpls()` and `pigls()` for Conway-Maxwell-Poisson and Poisson inverse
  Gaussian counts, respectively. Both use separate log-mean and log-dispersion
  predictors, support REML, NCV, and QNCV in `gam()`, and provide simulation,
  CDF, and quantile callbacks.

- Add `pig()` with fixed or estimated global dispersion, distribution callbacks,
  REML, ML, NCV, and QNCV in `gam()`, and ordinary and discrete `bam()` fitting.

- Add `cgammals()` with separate log-mean and log-dispersion predictors, REML,
  NCV, and QNCV fitting in `gam()`, and distribution callbacks.

- Migrate 14 additional mgcv family constructors, their distribution callbacks,
  mathematical notes, and regression tests from mgcvUtils development branches.
  Add NCV and QNCV for `dirichlet()`, `gpd()`, `nbls()`, and `ordbeta()` in
  `gam()`, and QNCV for `betabinomial()` and `clognormalls()`. Add ordinary
  and discrete `bam()` fitting for `ordbeta()` and verify additional fitting
  routes for `cgamma()`, `clognormal()`, `cmp()`, and `gpoisson()`. Include
  the family-specific `predict_ordbeta()` helper.

- All family constructors now supply native `gratia_diagnostics()` factories
  without requiring gratia. They bind fitted parameters to response moments,
  count left limits, ordinal and ordered-beta masses, censoring intervals, and
  beta-binomial count functions. Multivariate families provide component moments
  and marginal CDFs with joint simulation. Heavy-tail families report undefined
  moments explicitly.

- Censored constructors accept a known `observe` function for simulated
  diagnostic references. `censoring_process()` implements shared scalar limits.
  Interval PIT and quantile checks remain available when the rule is unknown.

- Accept `scaled.pearson` residuals for `nbls()`, `zinb()`, and `zanb()`.
  The `cumulative_link()` CDF accepts the `logp` argument used by other families.

## User visible changes

- Rename `gp()` to `gpoisson()` to distinguish it from `gpd()`. Update
  existing calls to the new name.

- Rename `ocat_link()` to `cumulative_link()`. Update existing calls to the
  new name.

- Replace the statmod dependency with internally generated and cached
  Gauss-Legendre quadrature rules for the censored gamma and log-normal families.

- Interpret `clognormal(theta = 0)` as estimation from the default start,
  matching the other scalar-dispersion constructors.

## Bug fixes

- Give `gpoisson()` positive working information at zero counts in ordinary
  `bam()`, while retaining its observed likelihood derivatives.

- Isolate estimated family state across direct fits and use the final nuisance
  parameters in likelihood and deviance reporting. Preserve mgcv's method flags
  during preinitialization and reject nonunit likelihood scales in `gam()` for
  the fixed-scale extended families.

- Accumulate shared predictor contributions to the Dirichlet and GPD REML
  Hessian and keep singleton predictor blocks intact in BFGS. Compute grouped
  NCV and QNCV gradients using the actual prediction rows and offsets.

- Preserve category probabilities in one-row cumulative-link predictions and
  stabilize derivatives for near-certain cloglog and loglog categories. Use the
  exact GPD scale link for fitted values and prediction standard errors.

- Handle censoring in gamma residuals and omit zero-weight rows from censored
  gamma derivative calculations. Honor the formula intercept in the CMP and
  scalar censored families' null deviance. Improve small-shape Dirichlet and
  high-shape censored gamma derivatives.

- Stabilize beta-binomial diagnostic probabilities at high precision. Preserve
  missing rows in zero-modified negative binomial quantiles and draws, and
  stabilize the zero-inflated negative binomial likelihood at large logits.
