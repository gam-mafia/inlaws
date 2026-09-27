# Validation of the pure-R CMP reference family

Validated on 2026-09-25 using R 4.6.1 (aarch64-apple-darwin23) and mgcv 1.9-4.
The reproducible base-R regression suite is `tests/cmp.R`.

## Numerical coverage

- Mean inversion, probabilities, log normalisers and variances: the Cartesian
  product of means 1e-6, 0.01, 0.2, 1, 2.5, 10, 100 and nu 0.1, 0.3, 1, 3,
  10, 30, compared with independent direct sums on extended support.
- All required first through fourth likelihood derivatives, including mixed
  mean/log-dispersion terms: means 0.2, 0.8, 2.5, 10 at nu 0.3, 1, 3.
  Neighbouring-value checks of derivatives differed by at most about 2e-9
  after scaling by 1 + absolute derivative, using summation and inversion
  tolerances of 1e-12 and finite-difference steps of 1e-5.
- Expected curvature derivatives: independently averaged over the response
  PMF, including at the Poisson point. Saturated-likelihood derivatives and
  the deviance/likelihood identity are checked separately.
- These are tested grids, not promises that arbitrary parameters outside
  them are supported, or that every boundary within their envelope is
  well-conditioned. Near-degenerate distributions can have unstable high
  derivatives even when their probabilities are accurately computed.

## Model checks

Fixed nu=1 reproduces Poisson GAM fits, smoothing parameters, predictions,
log likelihood, AIC, deviance, null deviance and deviance residuals, including
exposure offsets, zero weights and non-unit weights. The different fitting
paths have small optimisation differences, so model-level checks allow
relative tolerances up to 2e-5 rather than machine precision.

For intercept-only models the conditional mean estimate is the sample mean.
The independent reference evaluator in the tests inverts that mean using
`uniroot()` and direct sums, and minimises the restricted objective

    -ell + 0.5 log(n * mean(y)^2 / Var(Y)).

The omitted Laplace constant is independent of nu. This checks dispersion
estimation against an independently implemented REML objective for both
underdispersed and overdispersed samples. Other tests cover smooth-model
recovery, supplied dispersion starting values, simulation, quantiles, CDFs,
invalid inputs, numerical limits and reuse of the same family object.

An additional independent glmmTMB comparison used seed 9, 120 independent
uniform covariates x, mean exp(0.3+x), and nu=0.5. Unpenalised ML fits with
`gam(y~x, family=cmp(), method="ML")` and
`glmmTMB(y~x, family=compois())` differed by:

- maximum absolute coefficient difference: 4.1e-7;
- absolute nu difference (CMP nu versus 1/sigma(glmmTMB fit)): 1.2e-6;
- absolute log-likelihood difference: 4.6e-11.

The installed glmmTMB emitted a TMB build-version mismatch warning. This
comparison is supplementary; the regression suite uses independent R sums
and optimisation and has no dependency on glmmTMB.

## Representative runtime

A 100-observation estimated-dispersion REML GAM took approximately 3.7 seconds
on the validation machine. Reproduction:

```r
library(inlaws)
set.seed(24)
x <- runif(100)
off <- log(runif(100, 0.5, 2))
mu <- exp(0.3 + sin(6*x) + off)
y <- cmp(2)$rd(mu)
system.time(b <- gam(y ~ s(x, k=6) + offset(off),
                     family=cmp(), method="REML"))
b$family$getTheta(TRUE) # about 2.039; total EDF about 5.463
```

Runtime depends on support width, sample size, iterations and hardware. This
is a correctness-oriented R reference, not a large-data performance claim.

## Known limitations

Only REML is the supported smoothing-selection interface. A default-start
intercept-only ML fit on an underdispersed sample reached a trial nu near 23
with negative individual observed curvature, then failed in mgcv's
ML-specific BLAS/LAPACK path (DGESDD error -10), despite finite family
values and derivatives. The REML fit and ML fits started closer to the
solution succeeded. No workaround that changes the likelihood or its
curvature is used. ML and other optimisation methods require further
validation before being advertised.

bam(), varying dispersion, additional links, asymptotic approximations and
compiled acceleration are deferred. Summation limits fail explicitly.
CDFs and quantiles use the normalised finite support, so extreme-tail
accuracy is limited by the summation tolerance. All-zero positive-weight
responses have no finite log-mean fit and are rejected.
