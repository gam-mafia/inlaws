# Dirichlet mean/concentration family

For D = K + 1 response components, let eta contain K mean log-ratios and
one log-concentration predictor. The first response component is the reference:

    mu = softmax(0, eta[1:K])
    phi = exp(eta[D])
    alpha[j] = phi * mu[j]

For one observation, the unweighted log density is

    lgamma(phi) - sum(lgamma(alpha)) + sum((alpha - 1) * log(y)).

Observation weights multiply the entire density and every derivative. They
are likelihood weights, not multinomial trial counts. Zero-weight observations
are excluded before problematic arithmetic, as well as from initialization.

## Analytic derivatives

The implementation propagates finite multivariate Taylor polynomials of total
degree at most four. A column with multi-index a stores the derivative divided
by a! = product(factorial(a)). These are analytic derivatives, not finite
differences or numerical perturbations. Multiplication is the truncated
multi-index convolution; a univariate composition f(g) is

    sum(r = 0,...,4) f^(r)(g0) * (g - g0)^r / factorial(r).

Exponentials of the individual predictors have nonzero Taylor coefficients
only on the corresponding pure powers. Compute the softmax denominator,
compose its reciprocal, and multiply by the numerator and concentration.
Compose lgamma using digamma, trigamma, and higher polygamma derivatives.
The response logarithm contributes a linear term in alpha.

This automatically includes all mean/mean and mean/concentration cross
terms. It is necessary because mgcv's scalar link-transformation helper
cannot represent a coupled softmax transformation. Hold the rowwise
exponential stabilizing shift constant during differentiation: it cancels
between the numerator and denominator for every nearby predictor value.

Restore multi-index factorials before passing derivatives to gamlss.gH.
Nondecreasing predictor tuples enumerate columns in the same order as
trind.generator. The likelihood evaluator computes only the derivative orders
requested by mgcv. Fourth-order storage alone uses choose(D + 3, 4) columns;
the jet representation through fourth order uses choose(D + 4, 4) columns.
Convolutions also carry computational overhead, so this implementation targets
modest numbers of components.

## Initialization and numerical domain

A weighted pooled Dirichlet MLE provides initial log-ratios and concentration.
The pooled initialization optimization bounds log shapes to [-20, 20] solely
to produce finite starts. Moment concentration starts are bounded to [0.1, 1e4].
Neither bound constrains the fitted model. Joint penalized least squares
projects these starts onto all designs, accounting for shared coefficients
and offsets. Invalid likelihood trial points return -Inf so the fitter can
reduce its step; responses and fitted parameters are never silently clipped.

Ordinary floating-point lgamma/polygamma evaluation ultimately limits extreme
concentrations and component ratios. Tests exercise concentration from 0.001
to 10000 and small proportions; this is not a guarantee for arbitrary extremes.
For simulation, evaluate Gamma(alpha + 1) * U^(1/alpha) on the log scale to
avoid premature gamma underflow. Normalization can still produce exact zero
when a simulated proportion is below floating-point range. Such draws cannot
be passed directly back to the fitting family.

## Diagnostics and integration

The Dirichlet covariance is (diag(mu) - mu mu') / (phi + 1). On the tangent
space sum(y - mu) = 0, its quadratic distance is

    (phi + 1) * sum((y - mu)^2 / mu).

Pearson residuals are the square root of this distance times the observation
weight. They are nonnegative and are not normal residuals; use simulated QQ
comparisons. An unrestricted per-observation saturated likelihood is unbounded
as concentration increases around the observed composition. Therefore this
family does not supply likelihood-ratio deviance residuals or deviance explained.

Use qq.gam(fit, type = "pearson", rep = 100). The current gam.check() calls
k.check(), which requests default deviance residuals even when the user requests
Pearson residuals. It is consequently unsupported. No replacement k-index test
is claimed. Inspect componentwise response residuals against covariates and
compare models with increased basis dimensions when assessing basis adequacy.

predict.gam() drops column names returned by general-family prediction hooks.
The order always matches colnames(fit$y); simulation retains those names.
Also, summary.gam() currently computes its residual.df field with length(y)
for matrix responses. The fit's df.residual uses observation rows. Scale is
fixed at one and standard general-family asymptotic inference is used here;
this package does not patch mgcv's summary method.

## brms correspondence

brms::dirichlet(link = "logit", link_phi = "log", refcat = NULL) uses the same
first-component reference and shapes alpha = softmax(eta) * phi. brms expected
response predictions likewise return mean proportions. Our response-scale
standard errors use the full coefficient covariance and a delta approximation;
brms instead integrates predictions over posterior draws. Penalized REML and
Bayesian posterior summaries need not coincide.

Run inst/validation/dirichlet-brms.R for an optional comparison of density
values, conditional means, and generated Stan code. This requires brms but
neither Stan compilation nor sampling. See:
https://paulbuerkner.com/brms/reference/brmsfamily.html

## Fitting support

The family uses exported mgcv derivative helpers and needs no core patch.
REML is the supported smoothing criterion. NCV is not provided; mgcv may
warn and fall back to REML. Both bam paths reject general families in the
supported mgcv releases. Enabling bam requires a separate multiple-predictor
fitting-engine project, not a family flag or a change of class.
