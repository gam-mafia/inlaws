# Censored gamma mean and dispersion family

Write `u = log(mu)`, `t = log(phi)`, `k = exp(-t)`, and
`x = k*y/mu`. The gamma density on the original response scale is

```
log f(y) = k*log(x) - x - lgamma(k) - log(y).
```

An exact response contributes this density. A censored response contributes
`log(F(hi)-F(lo))`, with the same encoding and support as `cgamma()`:
left censoring is `(threshold, -Inf)`, right censoring is `(threshold, Inf)`,
and a finite ordered pair is an interval. Equal entries are exact. Zero lower
bounds are valid; `[0, Inf]` contributes zero log likelihood and derivatives.
Case weights multiply these contributions, without altering the gamma law.

## Derivatives and numerical evaluation

The evaluator propagates bivariate Taylor coefficients in `(u,t)` through total
order four using the polynomial arithmetic also used by the lognormal family.
The coefficient of `du^a dt^b` is the corresponding derivative divided by
`a! b!`. Derivatives of `lgamma` use the polygamma functions.

For the lower incomplete gamma ratio, when `x < k+1`, differentiate

```
P(k,x) = exp(k*log(x)-x-lgamma(k))/k * sum_j term_j
term_0 = 1
term_j = term_(j-1)*x/(k+j).
```

For `x >= k+1`, differentiate the upper-tail continued fraction using the
modified Lentz recurrence, as in `cgamma()`. Its prefactor is
`exp(k*log(x)-x-lgamma(k))`. Convergence must hold for every retained Taylor
coefficient, including at integer shapes where a scalar coefficient can vanish
while its derivatives remain nonzero. Nonconvergence is an explicit error;
there is no finite-difference fallback in fitting.

Use base R `pgamma(log.p=TRUE)` values before complementing or subtracting tails.
Interval probabilities are subtracted in the smaller tail using `expm1` and
Taylor arithmetic. This works directly with log-probability derivatives, avoiding
subtraction of large boundary-density derivative polynomials in the far tail.
For narrow positive intervals (`width < 1e-4*lo` and `k*width/mu < 1`), integrate
the density jet using 16-point Gauss-Legendre quadrature, factoring out the
full midpoint log-density jet to prevent both underflow and cancellation of
high-order derivatives. Bounds zero and infinity are handled
as constant zero/one probabilities rather than differentiating infinite jets.

The likelihood packs derivatives in mgcv's symmetric tensor order and uses
`gamlss.gH` for coefficient contractions. NCV passes row indices through mgcv's
subsetting machinery and restores censoring and offsets in its likelihood
callback. This also preserves reordered neighbourhoods and NCV gamma adjustments.

## Diagnostics and initialization

The saturated mean for an exact observation is `y`. For a positive finite
interval it is `(hi-lo)/log(hi/lo)`, obtained by equating the scaled boundary
densities. One-sided events have supremum probability one; intervals starting
at zero have the same limiting optimum. Deviance fixes each fitted dispersion
and is twice the weighted saturated-minus-fitted log likelihood. Null deviance
optimizes a constant mean intercept, retaining the mean offset and row-specific
fitted dispersions. These are conditional mean diagnostics.

Initialization uses censor-based pseudo-responses (midpoints for finite
intervals, half-thresholds for left censoring), a weighted log-response regression,
and the weighted mean of squared relative residuals for dispersion. Only the
initial dispersion is limited to `[0.05,5]`; the fitted parameter is unbounded.

## Validation

`tests/cgammals-derivatives.R` compares all mixed derivatives through order four
with numerical derivatives, and pure dispersion derivatives with independent
integrals of differentiated gamma densities. `tests/cgammals.R` checks coefficient
contractions, aligned mgcv gamma fits, independent likelihood optimization,
REML, NCV and its score derivatives, and prediction/distribution interfaces.
