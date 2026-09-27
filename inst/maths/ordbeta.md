# Ordered beta family: likelihood and derivatives

This implementation was independently derived for inlaws. It does not
copy mgcv family source. The model is from Kubinec (2023), Political Analysis
31:519–536, https://doi.org/10.1017/pan.2022.20.

Let `eta = logit(mu)`, `phi > 0`, and `c1 < c2`. The probabilities of zero,
an interior response, and one are respectively

    p0 = logistic(c1 - eta)
    pI = logistic(eta - c1) - logistic(eta - c2)
    p1 = logistic(eta - c2).

An interior response has beta density with shapes `a = mu*phi` and
`b = (1-mu)*phi`. Its log likelihood is

    log(pI) + lgamma(phi) - lgamma(a) - lgamma(b)
            + (a-1)*log(y) + (b-1)*log1p(-y).

The endpoint log likelihoods are `log(p0)` and `log(p1)`. Observation
weights multiply log likelihood contributions; zero weights are excluded
before evaluation, including unsupported endpoints in a reduced model.
The optimization deviance is negative twice log likelihood and the `ls`
terms are zero. The extra GLM scale is fixed at one.

## Parameterization and exact reductions

Full-model optimization uses `theta = (log(phi), c1, log(c2-c1))`. Without
zeros, `c1 = -Inf` and only `(log(phi), c2)` are estimated. Without ones,
`c2 = Inf` and only `(log(phi), c1)` are estimated. Without either endpoint
the beta likelihood estimates just `log(phi)`. A dataset with no interior
observations cannot identify precision and is rejected.

These reductions impose zero population probability on an absent endpoint.
They are boundary likelihood models, not estimates that an unobserved
endpoint can never occur in the underlying population. Prediction uncertainty
does not account for selecting the reduced model.

For gap `g = exp(theta[3])`, stable interior log probability is

    -softplus(c1-eta) - softplus(eta-c2) + log(1-exp(-g)).

`softplus(x) = max(x,0) + log1p(exp(-abs(x)))`. The last term uses
`log(-expm1(-g))` for small `g` and `log1p(-exp(-g))` otherwise. Missing
cutpoints remove the corresponding softplus term and the gap term.

For differentiation with respect to mu, use a further algebraic rewrite to
avoid cancelling derivatives of logit(mu). Define

    D(c) = logistic(c)*(1-mu) + logistic(-c)*mu.

Then `log(p0) = -softplus(-c1) + log(1-mu) - log(D(c1))`, and
`log(p1) = -softplus(c2) + log(mu) - log(D(c2))`. The full interior term is

    log(mu) + log(1-mu) - log(D(c1)) - log(D(c2))
      - softplus(-c2) - softplus(c1) + log(1-exp(-g)).

This form keeps finite endpoint derivatives accurate when mu approaches
the endpoint carrying almost all the mass. Logistic complements are
evaluated separately instead of subtracting rounded probabilities from one.

## Analytic Taylor differentiation

`ordbeta-derivatives.R` implements forward analytic differentiation by
truncated multivariate Taylor arithmetic, vectorized over observations.
For variables `(mu, theta)`, coefficient `T[alpha]` is

    partial^alpha(f) / product(factorial(alpha)).

Multiplication uses polynomial convolution. For a unary function `h`,
composition uses the finite Taylor expansion around the constant term:

    h(T) = sum_{k=0}^4 h^(k)(T[0]) * (T-T[0])^k / k!.

Terms have total degree at most four and nuisance-parameter degree at most
two. This downward-closed set includes every derivative requested by
`gam.fit4`, without computing unused higher nuisance derivatives. Lower
`level` requests use smaller sets. This is analytic differentiation using
known derivatives, not finite differences or a statistical approximation.
Finite differences are used only in validation.

The elementary derivatives are those of exp, log, lgamma (via polygamma),
logistic, softplus, and log(1-exp(-x)). In particular, writing
`h = 1/expm1(x)`, the first four derivatives of the last function are

    h
    -h*(1+h)
    h*(1+h)*(1+2*h)
    -h*(1+h)*(1+6*h+6*h^2).

Actual derivatives are recovered by multiplying Taylor coefficients by
their multi-index factorial. Nuisance Hessians are packed in mgcv's
row-wise upper-triangle order `(1,1), (1,2), ..., (2,2), ...`.

## Expected curvature

The categorical score with respect to eta is `-(1-p0)` at zero, `1-p1`
at one, and `p0-p1` for an interior response. Thus its Fisher information is

    Ic = p0*(1-p0)^2 + p1*(1-p1)^2 + pI*(p0-p1)^2.

The expected curvature of negative twice log likelihood with respect to mu is

    EDmu2 = 2 * (Ic / (mu^2*(1-mu)^2)
                   + pI*phi^2*(trigamma(a)+trigamma(b))).

The beta score has conditional expectation zero, eliminating its cross
term with the categorical score. First-order Taylor differentiation of this
formula supplies `EDmu3` and `EDmu2th`. Observed higher derivatives come
from differentiating the log likelihood itself.

## Response interpretation and diagnostics

The complete mean and variance are

    m = p1 + pI*mu
    v = p1 + pI*(mu^2 + mu*(1-mu)/(phi+1)) - m^2.

The derivatives with respect to eta are

    p0' = -p0*(1-p0)
    p1' =  p1*(1-p1)
    pI' =  pI*(p0-p1)
    m'  = p1' + pI'*mu + pI*mu*(1-mu).

Prediction standard errors multiply the linear-predictor SE by the absolute
appropriate derivative. Family parameters are treated as fixed, even when
the coefficient covariance includes mgcv's smoothing-parameter correction.
`fitted.values` remains the interior mean, consistent with the internal
logit link; response predictions and response/Pearson residuals use `m`.

The saturated likelihood holds nuisance parameters fixed and maximizes each
interior contribution over eta by scalar optimization. Endpoint suprema
are zero. Deviance residuals are signed square roots of twice the weighted
log likelihood difference. Null deviance optimizes a common intercept with
the fitted nuisance parameters and original offsets fixed (or retains only
the offset when there is no intercept).

For `0 <= q < 1`, the CDF is `p0 + pI*pbeta(q,a,b)`; at one it jumps to one.
Quantiles account for both atoms. Random generation selects the component
with a single uniform draw, then samples a beta value only for interiors.

## Integration and validation

Preinitialization uses mgcv's `needG` contract to access weights and returns
a newly constructed family. Each fit therefore owns its nuisance-parameter
closure; the setup object is not retained. This is an internal mgcv contract
and needs rechecking when supported mgcv versions change.

Base-R tests compare every derivative field with centered finite differences,
check Fisher information by integration, compare an independent MLE objective
and the beta reduction, recover simulated smooths, and check distribution,
prediction, weighting, and family-state behavior. No fitting-time numerical
differentiation or additional automatic-differentiation dependency is used.
