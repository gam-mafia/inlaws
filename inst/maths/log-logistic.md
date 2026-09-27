# Log-logistic arithmetic mean and log-response scale

Let log(Y) have a logistic distribution with location l and scale s, where
0 < s < 1. Put c(s) = log(sin(pi*s)/(pi*s)), and l = log(mu) + c(s).
The conventional log-logistic shape is 1/s and median is exp(l).
Then E(Y) = mu and, for s < 1/2,

    Var(Y) = mu^2 * {tan(pi*s)/(pi*s) - 1}.

For s >= 1/2 variance is infinite. The mean restriction excludes shape <= 1;
this is deliberately a mean model rather than an unrestricted median model.
The log-response standard deviation is pi*s/sqrt(3), not s.

## Likelihood and saturation

With z = (log(y) - log(mu) - c(s))/s, the original-response log density is

    ell = dlogis(z, log = TRUE) - log(s) - log(y).

At fixed s, saturation occurs at log(mu) = log(y) - c(s), or z = 0.
Thus ell_sat = -log(4) - log(s) - log(y). Deviance is twice the weighted
saturated-minus-fitted log likelihood, and its residual sign compares y to
the fitted median. The expected negative curvature with respect to mu is
1/(3*mu^2*s^2), even when response variance is infinite.

For theta = logit(s), the extended-family expected deviance curvature is

    EDmu2 = 2*w/(3*mu^2*s^2)
    EDmu3 = -2*EDmu2/mu
    EDmu2th = -2*(1-s)*EDmu2.

The additional mgcv likelihood scale is fixed at one. Replication weights
multiply row likelihoods, without changing the response distribution.

## Analytic derivatives and stability

The shared evaluator propagates bivariate Taylor coefficients through total
order four, using the existing polynomial convolution/composition helper.
Variables are (mu, theta) for the extended family and (log(mu), theta) for
the general family. All fitting derivatives are analytic; finite differences
are used only in tests. The extended family requires at most two theta
derivatives, while the general family requires every mixed fourth derivative.

For p = plogis(z), q = plogis(-z), derivatives of log logistic density are
q-p, -2*p*q, -2*p*q*(q-p), and -2*p*q*(1-6*p*q). Evaluating both tails
avoids cancellation from forming 1-p. Softplus uses
max(t,0) + log1p(exp(-abs(t))). Its analytic derivatives generate log(s),
s and 1-s without overflowing exp(theta).

For s > 1/2, c(s) = c(1-s) - theta. This retains accurate derivatives near
the mean-existence boundary. Near zero, log sinc uses its even power series
through degree 12, including derivatives; elsewhere use sinpi/cospi formulas.
Variance uses a small-argument series to avoid tan(x)/x - 1 cancellation.
Distribution callbacks use base R logistic CDF/quantile/random generators on
the log-response scale, supporting upper tails and log probabilities.

## Fitting contracts and validation

The extended family supports gam REML/ML/NCV, ordinary bam REML/ML/fREML,
and discrete bam fREML/NCV (the latter requires mgcv >= 1.9-4). The general
family supports gam REML/NCV with two separate coefficient blocks. Its NCV
adapter subsets offsets using observation indices and applies them exactly
once, including the full-data gamma adjustment.

The family is cloned before mgcv captures its scale state, so later fits do
not mutate earlier fits. No explicit scale argument is needed in bam: its
extended-family fitter fixes the additional scale at one, implicitly for
nb() and via family$scale = 1 for pig() and log_logistic(). Default and
explicit scale=1 calls give identical fits, covariance, family parameters,
smoothing criteria and AIC across the supported routes.

There is a separate mgcv 1.9-4 reporting inconsistency, also present for nb():
the outer bam function can set scale.estimated = TRUE despite the inner
fitter fixing scale at one. logLik.gam then adds one to the reported log
likelihood and its degrees of freedom. The tests retain explicit scale=1
for absolute likelihood checks and compare default fits with those reference
fits. This does not impose a caller requirement or justify refitting models.

Validation covers numerical mixed derivatives, independent likelihoods and
moments, coefficient contractions, NCV directional derivatives, fitting
routes, replication weights, offsets, predictions, and distribution tails.

Reference: [flexsurv log-logistic distribution documentation](https://search.r-project.org/CRAN/refmans/flexsurv/help/dllogis.html).
