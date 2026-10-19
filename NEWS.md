# inlaws 0.0.0.9000

## New features

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
  Extended-family fitting for `cgamma()`, `clognormal()`, `cmp()` and
  `gpoisson()` now has regression coverage for GAM REML/ML, ordinary
  BAM REML/ML/fREML and discrete BAM fREML/NCV, with fixed or estimated
  dispersion, offsets and case weights.
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

## Bug fixes

* Isolate mutable state for `cgamma()`, `clognormal()`, `cmp()` and
  `gpoisson()` before BAM fitting,
  preserving constructor objects and earlier fits. Report deviance and
  likelihood statistics using fitted dispersion, including discrete BAM.

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
