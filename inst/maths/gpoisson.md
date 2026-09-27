# GP-1 likelihood and derivatives

The implemented domain is mu > 0 and phi >= 1. Write t = log(phi - 1),
c = phi^(-1/2), b = 1-c, z = c*mu + b*y, and h = logistic(t). Then

    log p = log(mu) + log(c) + (y-1)*log(z) - z - lgamma(y+1).

Evaluation uses `dpois(y, z, log=TRUE) + log(mu) + log(c) - log(z)`
for positive counts, retaining R's stable Poisson calculation. For zero,
log p = -c*mu exactly. At phi=1 the distribution is exactly Poisson.

## Analytic derivatives

The derivatives of c with respect to t are

    c'  = -c*h/2
    c'' = c*(3*h^2/4 - h/2).

Put u=c/z, v=c'*(mu-y)/z. The derivatives used by `Dd` follow from

    L_mu   = 1/mu + (y-1)*u - c
    L_mu^k = (-1)^(k-1)*(k-1)! * (mu^(-k) + (y-1)*u^k), k >= 2
    u_t    = c'*y/z^2
    u_tt   = c''*y/z^2 - 2*u_t*v
    L_t    = -h/2 + (y-1)*v - c'*(mu-y)
    L_tt   = -h*(1-h)/2 + (y-1)*(c''*(mu-y)/z-v^2) - c''*(mu-y).

Differentiating the mu derivatives in t gives all mixed derivatives through
mu order four and t order two (total order at most four). Zero counts are
handled separately to avoid cancellation. No finite differences are used
in fitting. Zero-weight observations are excluded before evaluation.

## Saturated likelihood

For positive y, maximizing in a=c*mu gives

    a^2 - y*c*a - b*y = 0
    a_sat = (y*c + sqrt((y*c)^2 + 4*b*y))/2.

At zero the supremum is zero, attained as mu tends to zero. If S is the
saturated log likelihood, its derivatives are obtained by the envelope theorem:

    S_t  = L_t(mu_sat, t)
    S_tt = L_tt(mu_sat, t) - L_mut(mu_sat, t)^2 / L_mumu(mu_sat, t).

The family deviance is 2*w*(S-L); `ls` returns the weighted sum of S and
its dispersion derivatives. The likelihood constants are retained for AIC.
The extra scale is fixed at one. Each new fit rebuilds the family closure.

Observed mu information estimates expected information in mgcv's `EDmu2`
slot, as permitted by the extended-family interface. Exact expected Poisson
information is used for fixed phi=1. Likelihood optimization uses the full
observed log-link curvature, which can be negative for individual observations.

Fixed phi=1 has no internal nuisance coordinate: `getTheta()` returns
`numeric(0)` because mgcv's C interface rejects an infinite coordinate even
when fixed. `getTheta(TRUE)` returns one. Distribution helpers internally
use t=-Inf for the exact Poisson calculation.

## Distribution callbacks

The CDF sums actual probabilities through floor(q) in log-scaled blocks.
Quantiles sum until crossing the requested probability. Neither operation
renormalizes a finite support approximation. Both have a one-million-term
limit and fail explicitly if exhausted. Values at probability zero/one and
zero means are handled directly. Floating-point roundoff can affect extreme
tails near one.

Simulation uses a Poisson branching representation: draw Poisson(a)
immigrants, then independently generate Poisson(b) offspring per individual
until extinction. The sum over generations has the generalized Poisson
law. An explicit generation limit and integer-range checks prevent silent
nonconvergence or overflow. Weights do not change the sampling distribution.

## Upstream ML limitation

In mgcv 1.9-4, an entirely unpenalized ML model with negative observed
log-link curvature can enter `MLpenalty1` in `src/gdi.c` with zero-dimensional
penalty range space. Its negative-weight correction calls SVD on a matrix
with zero columns, producing a LAPACK error. REML avoids that ML-only path;
penalized ML models also avoid the zero-dimensional correction. This family
does not change likelihood derivatives to disguise the upstream issue.
