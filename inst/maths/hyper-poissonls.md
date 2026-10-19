# Hyper-Poisson mean and dispersion predictors

`hyper_poissonls()` uses `eta1 = log(mu)` and `eta2 = log(sigma)`.
The mass on the nonnegative integers is

    p(y) = lambda^y Gamma(sigma) / {Gamma(sigma+y) Z(lambda,sigma)},
    Z(lambda,sigma) = sum_{j>=0} lambda^j Gamma(sigma)/Gamma(sigma+j).

This is the HYPERPO2 parameterization in Hernandez-Barajas et al. (2026),
<https://journal.r-project.org/articles/RJ-2026-004/>. Here lambda is an
internal rate, obtained by inversion, and is never a third regression parameter.
The recurrence derived directly from the mass is

    (sigma+y) p(y+1) = lambda p(y).

In particular, the rate and dispersion in this recurrence must not be
interchanged. Summing gives

    mu = lambda - (sigma-1)*(1-p(0)),
    variance = lambda + (lambda-sigma+1)*mu - mu^2.

The implementation computes centered variance, avoiding cancellation in the
second identity. The identity is an independent validation check.

## Summation and inversion

Let `a = log(lambda)`. Relative weights are accumulated from a mode using
`w(j+1)/w(j) = exp(a)/(sigma+j)`, without subtracting large log-Gamma values.
Both omitted tails are bounded using the monotone recurrence, multiplied by
`abs(j-mode)^r`, for `r = 0,...,5`. Beyond either endpoint the ratio of these
terms decreases, so a geometric series bounds the remainder. The target is
`sum_tol * sum(w*(1+abs(j-mode))^r)`. These envelopes also bound all the
Taylor coefficients needed below: each log-dispersion derivative of a single
recurrence step is bounded, and polynomial products through total degree five
grow at most as the fifth power of distance from the mode. The bounds control
absolute derivative accuracy at the corresponding moment scale, not relative
error in a derivative that vanishes.

In the natural rate direction,

    d E(Y) / d a = Var(Y) > 0.

The expectation increases from zero to infinity, so every positive mu and
sigma has a unique finite rate. A safeguarded Newton solve uses the variance
as Jacobian, tracks lower and upper brackets, and bisects or expands when the
Newton step leaves the bracket. It stops at relative mean error `mean_tol`.
The defaults are `sum_tol = mean_tol = 1e-10`, `max_terms = 100000`, and
`max_iter = 100`. Exhaustion errors identify the failing operation. Likelihood
only trial steps outside numerical limits return negative infinity so mgcv can
backtrack; accepted derivative evaluations and public callbacks fail explicitly.
The mode must be representable below `2^50` for reliable integer arithmetic.

At fixed mu, sigma below one is underdispersed and sigma above one is
overdispersed. Sigma one is exactly Poisson. As sigma tends to infinity,
lambda/sigma tends to mu/(1+mu), giving a geometric distribution with variance
mu*(1+mu). As sigma tends to zero, the limit is Bernoulli for mu <= 1, or
1 + Poisson(mu-1) for mu > 1. Thus the family cannot represent arbitrary
underdispersion or arbitrarily large variance at fixed mean. These are limits,
not fitted bounds; finite dispersion estimates need not exist at a boundary.

## Analytic derivatives

Write `b_y(t) = log Gamma(exp(t)) - log Gamma(exp(t)+y)` and choose a mode m.
Centering at the fitted mean gives the normalizer variation

    G(da,dt) = log E[exp{(Y-mu)*da +
                         b_Y(t+dt)-b_m(t+dt)-b_Y(t)+b_m(t)}].

Then `G_a(0,0)=0`, `G_aa(0,0)=Var(Y)`, and the new mean is
`mu + G_a(da,dt)`. Recurrence-step derivatives of `-log(sigma+k)` are
polynomials in `u=sigma/(sigma+k)` and `v=k/(sigma+k)`:

    -u,
    -u*v,
    -u*v*(1-2*u),
    -u*v*(1-6*u+6*u^2),
    -u*v*(1-14*u+36*u^2-24*u^3).

Sum these steps relative to the mode. Center each resulting coefficient by
its expectation before exponentiating, and restore that expectation after
taking the logarithm. The residual moment-generating polynomial then starts
at degree two, so its log through degree five is `M-M^2/2`. This produces G
without cancellation of large polygamma terms. Distant observation-to-mode
differences use log-Gamma ratio and polygamma formulas; a Stirling correction
stabilizes the log-Gamma ratio at large arguments.

Substitute `dmu=mu*(exp(h1)-1)` and `dt=h2`. Four formal Newton updates with
the constant Jacobian Var(Y) solve `G_a(da,dt)=dmu` through degree four.
The likelihood variation is

    (y-mu)*da + b_y(t+dt)-b_m(t+dt)-b_y(t)+b_m(t) - G(da,dt).

Taylor coefficients include multi-index factorial divisors, removed before
packing the symmetric derivative arrays for `mgcv::gamlss.gH()`. Normalizer
derivatives through degree five supply the implicit mean derivatives needed
for fourth-order likelihood derivatives. No numerical differentiation is used
in fitting, and no Poisson shortcut removes dispersion derivatives at sigma one.
Weights multiply every coefficient; zero-weight rows are neutralized before
parameter transformation. Shared coefficient contractions reuse the package's
fourth-order trace workaround, and intercept-only predictors retain matrix
dimensions in the trace-only third-order route.

## REML and cross-validation

REML supports outer Newton (fourth derivatives), outer BFGS (third derivatives),
and EFS (second derivatives). EFS uses gamma one and requires all smoothing
parameters free and unlinked, or all fixed, because its mgcv driver does not
apply reduced smoothing-parameter maps. A preinitialization hook supplies
`offxset` for mgcv 1.9-3's misspelled EFS offset argument. NCV/QNCV use BFGS.

The NCV callback follows the existing signed-Hessian machinery in scatlss.
It uses mgcv's conjugate-gradient deletion solve with all curvature signs
retained; CG iteration counts are not reported as Cholesky failures. Exhausted
deletion solves and nonfinite predictions fail explicitly. Negative derivative
requests use local direct coefficient jackknife solves; shared score/Hessian
entries accumulate every predictor contribution.

With full-fit contributions ell0 and predictor changes Delta, scores are

    NCV  = -sum{gamma*ell(eta_cv) + (1-gamma)*ell0},
    QNCV = -sum{ell0 + gamma*(ell1^T Delta + .5 Delta^T ell2 Delta)}.

Scores, smoothing gradients, and row-gradient cross-products Vg are assembled
on actual prediction indices, including repetitions and subsets. Both offsets
are retained. The shared missing-row neighbourhood guard rejects inconsistent
upstream preprocessing and permits valid whole-group omissions; it does not
reconstruct neighbourhoods. On mgcv 1.9-4, fully unpenalized NCV/QNCV models
must disable jackknife covariance because of its empty-penalty reparameterization
defect. Default leave-one-out for penalized models is supported.

## Distribution and diagnostic callbacks

Fitted and response prediction matrices contain mu, sigma, in that order.
Delta-method standard errors describe parameter estimation, using mgcv's selected
coefficient covariance, including shared terms. Additional likelihood scale is
one. `bam()`, `gamm()` and ML smoothing selection are unsupported.

CDFs directly sum the requested small tail with a geometric remainder bound;
large tails use stable log complements. They can extend beyond normalization
support. Quantiles bracket and search using log-tail comparisons, preserving
tiny upper-tail probabilities, and random draws invert uniforms. The zero-mean
callback limit is a point mass at zero. Missing inputs propagate and prior
weights do not modify an individual response distribution.

Poisson probability callbacks use `ppois()` and `qpois()` directly, preserving
their discrete-boundary conventions without changing dispersion derivatives in
the likelihood. General quantile comparisons allow eight machine epsilons at
the target log-probability scale so complement rounding at an exact CDF step
does not advance a quantile by one.

The native gratia factory binds mean, centered variance, CDF, left CDF,
quantiles, simulation and residual functions to the supplied parameter rows.
It explicitly identifies a discrete response. Left limits use
`F(ceiling(q)-1)`, including between integers and at zero. Probability and
simulation diagnostics require no gratia runtime dependency for fitting.

Deviance saturates mu at y while holding sigma fixed; at y zero the saturated
log mass is zero. Null deviance optimizes a constant log-mean plus its offset,
holding fitted sigma fixed, or uses the offset alone without a mean intercept.
Pearson residuals use centered response variance and square-root case weights.
These are conditional mean-fit diagnostics, not joint saturation over sigma.

## Validation

The hyper-Poisson tests use independent wide-support sums and `uniroot()`;
Poisson reductions include log tails far outside normalization support.
Every derivative order is checked with differences of lower partials and
independent tensor-product differences. Coefficient, sandwich and REML
contractions include shared terms. NCV/QNCV references assemble signed deletion
Hessians directly and compare scores, gradients and jackknife perturbations;
fixed-basis deletion refits assess the one-step approximation. Model tests cover
the supported optimizers, prediction, missingness, weights and serialization,
plus optional gratia quantile, QQ, worm and appraisal integration.
