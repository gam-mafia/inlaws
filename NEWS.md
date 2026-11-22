# inlaws 0.0.0.9000

## New features

* Added `bgevlss()` for blended GEV responses with separate location
  quantile, log-spread and log-shape predictors. Includes analytic
  fourth-order derivatives, REML Newton/BFGS/EFS and NCV/QNCV BFGS in
  `gam()`, shape-constrained REML/EFS in `scasm()`, distribution callbacks
  and native gratia diagnostics.
  Constrained fits report AIC and log likelihood at the returned coefficients.

* Added `egpdlss()` for extended generalized Pareto model 1 with separate
  scale, tail shape, and bulk shape predictors. Includes distribution and
  gratia diagnostic callbacks, REML with Newton, BFGS and EFS in `gam()`,
  NCV/QNCV with BFGS, and REML with EFS for constrained fits in `scasm()`.

* Added `weibulls()` for separate Weibull arithmetic-mean and shape
  predictors in `gam()` (REML, NCV and QNCV) and `scasm()` (REML/EFS),
  and `weibull()` for a mean predictor with estimated or fixed global shape
  in `gam()` and ordinary or discrete `bam()`. Both use the exact positive
  response likelihood, analytic fourth-order derivatives, and native
  response-distribution and gratia diagnostic callbacks.

* Added `hyper_poissonls()` with exact log-mean and log-dispersion predictors
  for hyper-Poisson counts, shared terms and offsets, REML Newton/BFGS/EFS,
  NCV/QNCV BFGS, distribution callbacks with stable log tails and CDF-step
  quantiles, and native gratia diagnostics.

* Added `scatlss()` for Student-t location, scale and degrees-of-freedom
  predictors, stable density and deviance evaluation, analytic fourth-order
  derivatives, REML Newton/BFGS and NCV/QNCV BFGS in `gam()`, distribution
  callbacks and native diagnostics.
  NCV/QNCV are tested on mgcv 1.9-4 with signed-Hessian deletion solves
  and local score/jackknife compatibility measures; the extended-family
  response-subsetting restriction does not apply to this numeric response.
  Uses the shared structure guard for missing-row neighbourhood preprocessing,
  retaining valid whole-group omissions.

* Added `betals()` for beta responses with separate logit-mean and
  log-precision smooth predictors, REML and NCV fitting in `gam()`,
  simulation, CDF, quantile, and diagnostic callbacks.

- Add `log_logistic()` and `log_logisticls()` for arithmetic mean and
  log-response logistic scale, using log/logit links. Includes distribution
  callbacks, analytic fourth-order derivatives, and REML, NCV and QNCV in `gam()`.
  The single-predictor family also supports ML and discrete `bam()`.

* Added `cmpls()` with separate log-mean and log-dispersion predictors for
  Conway-Maxwell-Poisson counts, REML, NCV, and QNCV support in `gam()`, and
  simulation, CDF, and quantile callbacks.

* Added `pigls()` for Poisson inverse Gaussian counts with separate log-mean
  and log-dispersion predictors, REML, NCV, and QNCV fitting in `gam()`, and
  simulation, CDF, and quantile callbacks.

* Added `pig()` for Poisson inverse Gaussian count models with fixed or
  estimated global dispersion, REML/ML/NCV/QNCV support, ordinary and discrete
  `bam()` fitting, and simulation, CDF, and quantile callbacks.

* Added `cgammals()` for censored gamma models with separate log-mean and
  log-dispersion predictors, REML, NCV, and QNCV fitting, and distribution callbacks.

* Migrated 14 additional mgcv family constructors, their distribution callbacks,
  mathematical notes, and regression tests from mgcvUtils development branches.
  Extended-family fitting for `cgamma()`, `clognormal()`, `cmp()` and
  `gpoisson()` now has regression coverage for GAM REML/ML, ordinary
  BAM REML/ML/fREML and discrete BAM fREML/NCV, with fixed or estimated
  dispersion, offsets and case weights.
  Native response predictions for `zinb()` and `zanb()` include the complete
  response mean, variance, and zero probability after the original three
  component parameters, with full-covariance delta-method standard errors.
  NCV and QNCV are available for `dirichlet()`, `gpdlss()`, `nbls()`,
  and `ordbeta()`, and QNCV is available for `betabinomial()` and
  `clognormalls()`. `ordbeta()` supports ordinary and discrete `bam()`.
  The `predict_ordbeta()` helper exposes the response components.
  `dirmult()` has stable analytic derivatives through fourth order for REML
  Newton and BFGS, and NCV/QNCV with BFGS.

* Accept `scaled.pearson` residuals for `nbls()`, `zinb()`, and `zanb()`.
  The `cumulative_link()` CDF accepts the `logp` argument used by other families.

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

* `gpdlss()` is the canonical name for the generalized Pareto family.
  `gpd()` remains an exported alias. Previously fitted models retain
  their family identifier and callbacks.

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

* Interpret `clognormal(theta = 0)` as estimation from the default start,
  matching the other scalar-dispersion constructors.

## Bug fixes

* Isolate estimated nuisance state across direct fits, including for the
  censored, ordinal, and count families. Report likelihood and deviance
  at the final parameters in discrete BAM. Reject nonunit likelihood scales
  in `gam()` for fixed-scale extended families.

* Stabilize generalized Poisson BAM working information and preserve the
  NCV settings when initializing its family.

* Reject GAM NCV and QNCV for `cgamma()`, `clognormal()`, `cmp()` and
  `gpoisson()` under mgcv 1.9-4 for independent upstream numerical NCV bugs;
  the censored families also lose response censoring information during
  held-out evaluation. Later versions require reassessment. Discrete BAM NCV
  remains available subject to the shared neighbourhood qualification below.

* Guard explicit NCV neighbourhoods corrupted by mgcv missing-row preprocessing
  in all NCV-capable general families, including QNCV and `nei` during REML.
  Correctly adjusted whole-group omissions and default leave-one-out after
  omission remain available. Document the complete-case/rebuild workaround
  and the interception limits for extended families, BAM and incomplete
  neighbourhood lists, whose upstream validation can precede family hooks.

* Accumulate shared predictor contributions to the Dirichlet and GPD REML
  Hessian and retain singleton predictor blocks in BFGS. Compute grouped
  NCV and QNCV gradients using the actual prediction rows and offsets.

* Preserve category probabilities in one-row cumulative-link predictions
  and stabilize derivatives for near-certain cloglog and loglog categories.
  Keep GPD's fitted scales and prediction standard errors on the exact log
  link. Stabilize ordinary ordered-beta BAM working steps and support
  response predictions from its discrete fits.

* Handle censoring in gamma residuals and omit zero-weight rows from
  censored gamma derivative calculations. Honor the formula intercept in
  the CMP and scalar censored families' null deviance. Improve small-shape
  Dirichlet and high-shape censored gamma derivatives.

* Stabilize beta-binomial diagnostic probabilities at high precision.
  Preserve missing rows in zero-modified negative binomial quantiles and
  draws and stabilize the zero-inflated likelihood at large logits.
