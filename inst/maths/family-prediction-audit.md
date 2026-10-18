# Native family prediction audit

Reviewed 2026-10-18 against main at 33b8fb1 and installed mgcv 1.9-4.
Integration review also includes main's ordered-beta update at c2b8af4.
The package requires mgcv >= 1.9-3. This audit concerns the native
`family$predict` interface used by `predict.gam(type = "response")`.
It does not introduce prediction wrappers or change fitting parameterizations.

## Inventory and decision rule

Inspection of every constructor in `R/` found **14** families without a
custom callback before this change: all 13 originally proposed candidates,
plus `betals()`. `cgen_gammalss()` inherits a callback from `gen_gammalss()`;
it is not an additional missing callback.

Add a callback when the existing response predictions leave a material gap
in interpreting the response distribution, and a joint transformation closes
that gap with correctly propagated uncertainty. Do not expand every family's
response matrix merely because additional distribution functionals exist.
In particular, a variance formula involving two predictors is evidence that
its SE needs covariance, but does not alone establish a need to replace a
useful, documented mean/dispersion prediction interface.

**Decision:** add callbacks to `zinb()` and `zanb()`. Retain the default for
the other 12 families. The former return a component mean that can differ
substantially from the response expectation, while the retained families
already expose the complete arithmetic mean or a documented mean proportion.

## mgcv's contract

Inspected installed `mgcv::predict.gam`, `mgcv:::gam.fit5`, and the callback
implementations in this package. These are the relevant behaviors:

* Without a custom callback, response prediction inverse-links each predictor
  independently. For predictor j its SE is
  `abs(linfo[[j]]$mu.eta(eta_j)) * sqrt(X_j Vb[j,j] X_j')` (with j denoting
  a coefficient index block). For a single predictor the family-level
  `linkinv` and `mu.eta` play the same role. This is correct for each
  component parameter; it does not provide SEs for a function combining them.
* The multi-predictor callback receives `family, se, eta, y, X, beta, off, Vb`.
  It returns `list(fit = ...)`, or `list(fit = ..., se.fit = ...)` when SEs
  are requested. In the prediction path, mgcv supplies X with attribute
  `lpi`, the coefficient vector, fixed offsets and the selected covariance.
  The callback can also evaluate supplied eta without X for values alone.
* Each predictor uses its own `lpi` coefficient indices; these can overlap.
  Derivatives for shared coefficients must be **added**, rather than assigned
  or propagated independently. Then contract the complete coefficient
  Jacobian with Vb, including cross-predictor covariance.
* `predict.gam` selects the smoothing-corrected covariance before calling
  the callback for `unconditional = TRUE`, when available. The callback
  must use the supplied Vb rather than fetch or cache a covariance itself.
* mgcv assembles predictions in blocks. Response output can have a different
  number of columns from `nlp`; its block assembler reallocates accordingly.
  In 1.9-4 it drops callback column names. The new output therefore has a
  documented fixed column order; no wrapper is added to repair names.
* `type = "link"`, `"lpmatrix"`, `"terms"`, and `"iterms"` retain mgcv's
  ordinary paths. `gam.fit5` inverse-links predictors into fitted values
  independently of response prediction. The callback must not change their
  dimensions or the input parameterization of simulation/CDF/quantile tools.
* A callback can receive observed y, but prediction on covariates alone must
  work. Neither response values nor training-row state is needed here.
  Formula offsets, terms excluded by mgcv, NA restoration, new rows and block
  sizes are handled through the supplied design and offsets.

All SEs below are first-order **estimation** SEs. They are neither predictive
standard deviations nor an integration over coefficient or random-effect
uncertainty. Unless requested and available, smoothing uncertainty is omitted.
For extended families, global estimated nuisance parameters do not become
regression coefficients in Vb, so adding a callback cannot by itself supply
their full estimation uncertainty.

## Per-family audit (before the additions)

Notation: mu is the family's arithmetic mean, phi dispersion or precision as
specified, sigma the standard deviation of log(Y), theta NB size, and p the
count-component probability. SE rules in this table describe existing default
response output, not the new six-column mixture interface.

| Family | Current response output and SE | Interpretation and useful joint transformations | Decision and reason |
|:--|:--|:--|:--|
| `betabinomial()` | `(mu, phi)`; slopes `mu*(1-mu)` and `phi` times respective link SEs. | mu is a mean success proportion, phi beta precision. For known trial total m, count mean is m*mu; its SE is m*SE(mu), needing no joint transformation. Count variance `m*mu*(1-mu)*(m+phi)/(1+phi)` needs both parameters and their covariance for its SE. | Retain. The proportion interface is explicit and supports predictions without trial totals. There is no universal newdata trial total; reusing fitted totals would silently misalign new rows. Appending precision or inverse links would add nothing. |
| `betals()` (additional discovery) | `(mu, phi)`; slopes `mu*(1-mu)` and `phi`. | Complete beta mean and precision. Variance `mu*(1-mu)/(1+phi)` combines predictors and needs covariance for an SE. | Retain. The response mean is already correct and the standard beta mean/precision interface is documented. An optional variance SE alone does not justify widening this interface. |
| `cgamma()` | mu, with `mu*SE(eta)`. | Arithmetic mean of the underlying uncensored gamma; variance `theta*mu^2`. Global theta is fixed or estimated outside predictor coefficients. Recorded censoring limits/events do not define a common numerical observed-response mean. | Retain. No bias correction is needed: the predictor already models the arithmetic mean. Adding theta would expose a global nuisance parameter, with incomplete uncertainty. |
| `cgammals()` | `(mu, phi)`; slopes mu and phi. | Uncensored gamma mean and dispersion. Variance `phi*mu^2` and shape/scale transformations are available algebraically; a variance or gamma-scale SE needs covariance. | Retain. The uncensored mean is correct. A callback cannot infer a future observation process from censoring intervals. An extra variance column would be optional distribution enrichment, rather than a correction to response semantics. |
| `clognormal()` | mu, with `mu*SE(eta)`. | Arithmetic mean of the uncensored lognormal, already parameterized as `meanlog = log(mu)-sigma^2/2`. Variance `mu^2*expm1(sigma^2)` and median depend on global sigma. | Retain. Multiplying mu by `exp(sigma^2/2)` would be wrong. Adding sigma/median would condition on a globally estimated nuisance parameter rather than improve its uncertainty. |
| `clognormalls()` | `(mu, sigma)`; slopes mu and sigma. | Uncensored arithmetic mean and log-response SD. Median `mu*exp(-sigma^2/2)` and variance `mu^2*expm1(sigma^2)` need both predictors and covariance for SEs. | Retain. The current response mean is already arithmetic, not a median needing conversion. Censoring-process moments require a specified future observation rule; `observe` is a diagnostic simulation rule, not an additional regression parameter. |
| `cmp()` | mu, with log-link delta SE. | Exact arithmetic count mean, not the conventional CMP intensity lambda. Variance and lambda require numerical distribution calculations conditional on global nu. | Retain. Numerical mean inversion is already in the likelihood, so no mean correction is needed. Adding nu/lambda would expose nuisance/distribution parameters without the missing joint nuisance covariance. |
| `gpoisson()` | mu, with log-link delta SE. | GP-1 arithmetic count mean. Variance `mu*phi` uses global phi. | Retain. A single predictor already represents the response expectation. A global dispersion column or duplicate inverse link brings no benefit. |
| `log_logistic()` | mu, with log-link delta SE. | Arithmetic response mean, with global log-response logistic scale `0 < s < 1`. Median `mu*sin(pi*s)/(pi*s)` differs from mu; variance exists only for s < 1/2. | Retain. There is no missing mean correction; additional functionals condition on global s. A generic always-finite variance/SE would be misleading. |
| `log_logisticls()` | `(mu, s)`; slopes mu and `s*(1-s)`. | Arithmetic mean and varying log-response logistic scale. Median combines both predictors; variance and its SE require covariance and exist only for s < 1/2. | Retain. The mean is already complete. Do not introduce a partially undefined variance column or expand the interface merely to expose another optional functional. |
| `nbls()` | `(mu, theta)`; mean-link slope mu (log), 1 (identity), or `2*eta` (sqrt); size slope theta. SEs multiply absolute slopes by link SEs. | Complete NB arithmetic mean and size. Variance `mu+mu^2/theta` and zero probability combine both predictors, so their SEs require covariance. | Retain. Unlike a mixture, its first response column is the complete count mean and theta completes a standard NB parameterization. Additional variance/zero-probability SEs are possible but do not resolve a semantic gap in the response prediction. |
| `pig()` | mu, with log-link delta SE. | Poisson inverse Gaussian arithmetic mean; variance `mu+theta*mu^2`, with global dispersion theta. | Retain. No response mean conversion or missing multi-predictor uncertainty; adding theta would expose a global nuisance parameter. |
| `zinb()` | `(mu, p, theta)`; slopes mu, `p*(1-p)`, theta. | Component mean, membership probability and NB size. Complete mean is p*mu, variance is `p*(mu+mu^2/theta)+p*(1-p)*mu^2`, actual zero probability `1-p+p*q0`. p is not P(Y>0). All three complete quantities involve multiple predictors and their covariance. | Add. Users otherwise have only component predictions for an intrinsically mixed response, with no native complete-moment SEs. Keep the first three columns, append complete mean, variance and p0, and supply full-covariance SEs. |
| `zanb()` | `(mu, p, theta)` with the same parameter-wise SE rules as ZINB. | mu is the **untruncated** NB mean; p=P(Y>0). Complete mean is `p*mu/(1-q0)`, not p*mu. Mean and variance combine three predictors; p0=1-p uses just the membership predictor. | Add. Truncation makes the component interpretation especially easy to misuse. A joint analytic transformation supplies complete moments and their SEs; p0 is the actual zero mass and has a simpler derivative. |

The retained defaults already account for coefficient correlations **within**
each parameter's coefficient block. Their limitation is that there is no
returned joint functional whose cross-parameter covariance could be propagated;
the existing parameter SEs themselves are not erroneous.

## Critical review of existing callbacks

* `zoibetals()` is a useful precedent for preserving component parameters
  alongside complete mixture probabilities and moments, and for accumulating
  overlapping coefficient blocks. Its seven-column shape and softmax logits
  are specific to inflated beta. NB mixtures have a different membership
  interpretation and hurdle normalization, so neither its number of outputs
  nor its derivative formulas can be copied. Its wrapper repairs names lost
  by mgcv; this audit deliberately uses only the native callback interface.
* `ordbeta()` converts a logit-linked interior beta mean into the complete
  mean, which is a clear semantic improvement. Main's subsequent update at
  c2b8af4 returns that mean alongside the interior mean and endpoint/interior
  probabilities, rather than returning only the complete mean. It still
  conditions on fitted global precision/cutpoints for its SE. That is not a
  template for uncertainty in the NB families, where size and membership are
  regression coefficients and their covariance must enter response moments.
  Preserving the NB component columns also retains their existing uses.
* `dirichlet()` needs a joint softmax: identity inverse links alone would
  expose log ratios, not proportions. Its Jacobian correctly combines mean
  predictors; concentration has zero derivative for mean proportions. This
  is stronger justification than merely attaching column names to inverse links.
* `dirmult()` similarly returns category probabilities plus concentration,
  using full covariance and adding shared-coefficient contributions. It does
  not fabricate newdata count totals. Its fitted offsets and category naming
  state are family-specific; the new NB callback needs no captured training
  state or observed responses.
* The remaining existing callbacks are `cumulative_link()`, `gpd()`,
  `cmpls()`, `pigls()`, `gen_gammalss()` and the inherited `cgen_gammalss()`
  callback. Their presence is not itself a reason to add callbacks elsewhere.
  In particular, callbacks that just inverse-link and name parameter columns
  do not establish a benefit under this audit's criteria.

## Implemented interface and validation

`zinb()` and `zanb()` now supply `.zinb_predict` from `R/zinb-predict.R`.
Response columns are **mu, p, theta, mean, variance, p0**. Fits and link
predictions remain three-column component quantities. Distribution callbacks
consume those three columns; code using all response columns as parameters
should explicitly select columns 1:3. Callback results are named, while
`predict.gam` may discard these names. Existing first-three-column values and
SEs are preserved. This is an appended native capability, documented in the
existing unreleased new-feature entry, the family help, and compatibility
vignette. No `predict_fam*` API or family-specific prediction wrapper is added.

`inst/maths/zinb-zanb.md` records the moment formulas, Jacobian contraction,
and small-mean/near-Poisson calculations. Tests in `tests/zinb-prediction.R`
use independent NB PMF sums and coefficient finite differences, not just the
same moment helper twice. They cover both families, full correlated and
smoothing-corrected covariance, real and synthetic overlapping coefficient
blocks, offsets in all predictors, covariates-only newdata, row reordering,
block sizes including one, one-row output, missing rows and `na.exclude`,
excluded smooths, series boundaries, extreme logits, tiny means (including
ratio underflow), and the Poisson limit. Existing model tests check unchanged
fitted parameters, distribution utilities and residuals against the new
response moments. The audit also guards the inventory and retained-default
decisions in the prediction regression test.

Completed validation on macOS Tahoe 26.6.1, R 4.6.1 (aarch64), mgcv 1.9-4:

* `R CMD INSTALL` into an isolated temporary library: passed.
* Standalone `tests/zinb-prediction.R` and `tests/zinb-models.R`: passed.
* Help regenerated with roxygen2 8.1.1; `man/zinb.Rd` matches that output.
* `R CMD build`: passed, including the compatibility vignette. The first
  sandboxed attempt could not write Quarto's Sass cache; rebuilding with
  access to the standard cache succeeded without source changes.
* `R CMD check inlaws_0.0.0.9000.tar.gz`: **Status: OK**, zero errors,
  warnings or notes. All **50** test scripts passed; examples, vignette
  rebuilding and PDF manual compilation passed. No validation-skipping
  build or check options were used.
* `git diff --check`: passed. Release-note additions were reviewed against
  the repository's change classification rules and amend the existing
  unreleased feature entry under New features.

Only the installed mgcv 1.9-4 was tested; this does not establish runtime
coverage of the minimum supported 1.9-3. Nor does a successful native
prediction test validate arbitrary unrepresentable predictor extremes;
existing family numerical and identifiability limitations continue to apply.
