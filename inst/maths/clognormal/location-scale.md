# Censored lognormal mean and scale regression

`clognormalls()` uses two predictors, a = log(mu) and t = log(sigma).
The normal location of log(Y) is m = a - exp(2t)/2. Consequently E[Y] = mu
and Var[Y] = mu^2 expm1(sigma^2), even when sigma varies across observations.
The exact and censored likelihood contributions are those in README.md.
Weights multiply those contributions, and exact observations retain -log(y).

## Derivatives and fitting

The shared evaluator has a predictor mode in which its two Taylor variables
are (a,t), rather than (mu,t). It includes *all* monomials a^j t^k with
j+k <= 4, including t^3, a*t^3, and t^4, which the extended family does not
need. Evaluating m as a - exp(2t)/2 automatically includes its scale dependence.
Coefficients multiplied by j! k! give the joint log-likelihood derivatives.
For each order they are packed with the number of t derivatives increasing:
for example, the Hessian order is aa, at, tt. mgcv::gamlss.gH performs the
coefficient-space contractions needed for REML, using these derivatives
already on the predictor scale; no further link transformation is applied.

Each link's inverse is the exact exponential. No variance/scale floor or
likelihood clipping is imposed. Initialization uses rough logged censor
bounds, fits a regularized location regression, estimates a residual scale,
and shifts the mean start by sigma^2/2. These are starting values only. The
scale model should not be so flexible that individual exact observations
can be interpolated with sigma approaching zero. Coefficients shared between
predictors are explicitly rejected: supported formulas have separate
coefficient blocks, which may involve the same covariates.

## NCV and censoring state

The family delegates neighbourhood calculations to mgcv's gamlss.ncv helper
and provides per-observation predictor derivatives and a direct-eta likelihood
path. This helper subsets its response as a vector, which would discard the
censor attribute or flatten a two-column response. A local adapter passes
row indices as the response and restores both observed bounds before each
likelihood evaluation. It does not mutate the fitted family or cache training
state in it. The adapter passes no offsets to the helper and adds the appropriate
row offsets in its likelihood callback. This keeps held-out evaluations and
the full-data adjustment for gamma consistent; mgcv otherwise adds offsets
to only the held-out predictors. Reported CV predictors have offsets restored. Tests explicitly reorder neighbourhood prediction
rows and check the resulting NCV score against an independent likelihood.
The adapter uses the unexported gamlss.ncv interface available in the declared
mgcv dependency; tests exercise this integration rather than duplicating its
numerical implementation.

## Inference and constant-scale comparison

The scale intercept is a regression coefficient, not an extra family
parameter. REML integrates it along with other unpenalized coefficients;
therefore `~ 1` for scale does not promise equality to the extended-family
REML fit. Unpenalized likelihood fits do agree with clognormal's ML fits.
Prediction columns are mu and sigma, and standard errors use the delta
method on each predictor. Posterior-mean/epsilon corrections and population
averaging over random effects are not included.

Deviance keeps each fitted sigma fixed and uses the saturated-location
likelihood. Null deviance optimizes a constant log-mean with the fitted sigma
vector and mean offsets retained. It is a mean-fit diagnostic, not a comparison
against a model having constant mean *and* scale. The distribution components
operate on uncensored responses using the two-column parameter matrix.
