# Zero-inflated and hurdle NB regression

Both families use predictors u = log(mu), b = logit(p), v = log(theta).
The underlying NB has variance mu + mu^2/theta. Write f(y) for its PMF,
q0 = exp(-a), a = theta log1p(mu/theta), and z = I(y > 0).

## Likelihoods and interpretation

ZINB is a mixture with log likelihood

    y = 0: log(1-p + p*q0)
    y > 0: log(p) + log f(y).

Its probability p describes entering the NB component; it is not P(Y>0).
ZANB has log likelihood

    (1-z) log(1-p) + z {log(p) + log f(y) - log(1-q0)}.

Here p = P(Y>0). The binomial likelihood separates from the mean/size
likelihood. Every derivative mixing b with u or v is exactly zero. Sharing
coefficients or smoothing parameters can nevertheless couple coefficient fits.

All three predictors are regression coefficients under general.family REML.
A constant size formula does not reproduce the auxiliary-size marginal
likelihood used by mgcv::nb(). Neither size nor probabilities are artificially
bounded. Separation, absent zeros, positive counts all one, or a Poisson limit
can lead to unidentified parameters. A large number of zeros does not by itself
establish zero inflation. All-zero positive-weight samples are rejected.

## Analytic derivatives

The existing nbls natural-parameter derivatives and gamlss.etamu provide the
NB derivatives with respect to u and v through fourth order. Pure u derivatives
are evaluated directly using s = mu/(mu+theta), t = 1-s:

    l_u    = y - (theta+y)*s
    l_uu   = -(theta+y)*s*t
    l_uuu  = -(theta+y)*s*t*(1-2*s)
    l_uuuu = -(theta+y)*s*t*(1-6*s+6*s^2).

This avoids subtraction of order-one quantities for tiny means. Small-count
likelihood values use the finite-product identity

    log f(y) = sum_{k=0}^{y-1} log1p(k/theta) - lgamma(y+1)
               + y*u - (theta+y)*log1p(mu/theta).

For larger counts in the near-Poisson regime, values use the same twelve-term
convergent expansion as the nbls size derivatives (see nbls.md).

Derivatives of mixture and truncation transformations are assembled by finite
Taylor polynomial arithmetic in three variables, truncated at degree four.
Each coefficient is a derivative divided by the multi-index factorial.
Products use multi-index convolution. A scalar composition uses

    F(x0 + dx) = sum_{k=0}^4 F^(k)(x0) * dx^k/k!.

This is analytic forward differentiation, not numerical differencing or a
statistical approximation to the likelihood. No differentiation dependency is
required. Coefficients are converted back to mgcv's symmetric derivative
packing order before coefficient contractions. EFS requires degree two, BFGS
and NCV degree three, and full Newton REML degree four.

The logistic log probabilities use softplus and its stable derivatives. For
ZINB zeros, the identity

    log P(Y=0) = log(1-p) + softplus(b + log(q0))

avoids forming the probability mixture on the response scale.

For ZANB, directly subtracting log f(y) and log(1-q0) loses precision when
mu is small. Put r=mu/theta and h(a)=log((1-exp(-a))/a). Then

    log(1-q0) = u + log(log1p(r)/r) + h(a).

Subtract u from the NB jet before adding the correction functions. For small
r, evaluate log1p(r)/r by its alternating power series through degree twelve.
For small a, use

    h(a) = -a/2 + a^2/24 - a^4/2880 + a^6/181440
           - a^8/9676800 + a^10/479001600 + O(a^12).

Elsewhere use log1p/expm1 and analytic compositions. The series switch points
are tested on both sides and the tiny-mean y=1 curvature is tested explicitly.

Weights multiply every observation likelihood and derivative. Inactive rows
are replaced by harmless parameters before calculations, then contribute zero.
The shared-coefficient fourth-order trace workaround from nbls is generalized
to any number of predictors; nbls's public interface is unchanged.

## NCV

The family uses mgcv's compiled Newton deletion updates via gamlss.ncv, with
its ordinary unadjusted score path. The returned deleted-fit predictors and
their smoothing derivatives are reused to calculate the requested score.
This avoids upstream QNCV bookkeeping assumptions about offsets and prediction
row ordering. No mgcv namespace or installed source is modified.

Let eta be the full-fit predictor at a predicted row, eta_cv the approximate
omission predictor, delta = eta_cv - eta, and g, H, T the first three
observation log-likelihood derivatives. For QNCV the negative score contribution
is

    -l(eta) - gamma {g' delta + delta' H delta/2}.

For a smoothing direction, its derivative simplifies to

    -(1-gamma) g' d_eta
    -gamma (g + H delta)' d_eta_cv
    -gamma T[d_eta, delta, delta]/2.

For ordinary NCV, use -gamma*l(eta_cv)+(gamma-1)*l(eta). All scores and
score derivatives include observation weights. Offsets enter both eta and
eta_cv; smoothing derivatives exclude fixed offsets. Repeated, reordered,
and partial prediction indices are handled explicitly. The per-row gradients
supply the covariance ingredient Vg. Jackknife coefficient perturbations are
preserved from mgcv. Full Newton NCV is not implemented; mgcv selects BFGS.

## Prediction and distribution utilities

Fitted values have named columns mu, p, theta. Response predictions return the
same three parameters and parameter-wise delta-method SEs. mgcv may discard
column names when assembling prediction blocks. Link predictions retain the
three predictor scales. Generic mgcv prediction handles shared coefficients,
term offsets, and new data without a custom prediction override.

For ZINB, E(Y)=p*mu and Var(Y)=p*(mu+mu^2/theta)+p*(1-p)*mu^2.
For ZANB, put m=mu/(1-q0). Its positive-component variance is

    vplus = m * {1 + mu*(1+1/theta) - m},

so E(Y)=p*m and Var(Y)=p*vplus+p*(1-p)*m^2. These moments define
response and Pearson residuals, while parameter predictions remain unchanged.

CDFs and quantiles use NB upper-tail probabilities in log space. For ZINB,
P(Y>q)=p*P(N>q). For ZANB and q>=0,
P(Y>q)=p*P(N>q)/(1-q0). Inverse upper-tail quantiles avoid rounding
q0+(1-q0)*u to one when truncation is extreme. Simulation applies these
quantiles to uniforms; weights do not alter an individual response law.

Deviance holds fitted theta fixed. Zero observations have saturated log
likelihood zero. ZINB positive observations saturate at p=1, mu=y.
ZANB positives saturate at p=1 and the mu solving mu/(1-q0)=y;
y=1 has limiting saturated log likelihood zero as mu tends to zero.
Null deviance optimizes constant u and b with their offsets retained and fitted
theta held fixed. It is a conditional diagnostic, not a joint test of all
three predictors.

## Validation

Tests differentiate each analytic order to check the next, including mixed
partials, small means, extreme logits, near-Poisson sizes, and series switches.
They independently check coefficient/smoothing contractions, shared coefficients,
weighted scores, fixed-basis deletion refits, NCV/QNCV gradients and offsets,
probability normalization, moments, simulation, independent MLEs, and hurdle
component separation. Regression checks cover all fitting methods, prediction,
weights, utilities, and the pre-existing nbls family.

For mgcv 1.9-3, the family retains nbls's supported preinitialize hook that
supplies the misspelled offxset field needed by EFS in that release. A separate
upstream fixed-sp setup bug affects formulas whose later predictors have no
smooths. Use an unfixed gam(..., fit=FALSE) setup and then gam(G=G, sp=...)
in 1.9-3, or use 1.9-4. The tests cover the setup workaround on 1.9-3 without
patching mgcv.
