# Censored generalized gamma mean, scale and transformed shape

The three predictors and finite-mean parameterization are identical to
`gen_gammalss()`: u = log(mu), t = log(sigma), v = log(r),
Q = (r - 1)/sigma. See `gen-gammalss.md` for the mean correction C and
Prentice location a = u - C, including its analytic continuation at Q = 0.
Censoring is non-informative conditional on covariates. Case weights multiply
observed log likelihoods without changing the underlying response law.

For exact observations use log f(y). A left-censored observation at h contributes
log F(h), a right-censored observation at l contributes log S(l), and a finite
interval (l, h) contributes log{F(h) - F(l)}. Zero lower bounds are allowed;
(0, infinity) contributes zero at every derivative order.

## Probability derivatives

Let T = log Y and z = (T - a)/sigma. Its log density is

    A(Q) - t - z^2 E(Q z),
    E(x) = (exp(x) - 1 - x)/x^2, E(0) = 1/2.

Both A and the mean correction use the same stable Taylor arithmetic as the
uncensored family. Integrating exp(log f_T) over fixed observed log bounds
integrates its full mixed parameter derivatives through order four. This
retains shape information at Q = 0 without evaluating singular incomplete-gamma
shape derivatives or replacing a neighbourhood of zero by a normal family.

For each one-sided probability, integrate the smaller tail and complement its
Taylor jet if necessary. A rational change of variable maps its infinite log
interval to (0, 1). Its scale is chosen using the boundary log-density slope
and the natural log-tail scale, and held fixed while differentiating.
A full log-density Taylor jet at the integration peak is factored out to
reduce cancellation in extreme-tail higher derivatives. Finite intervals
shorter than sigma on the log scale use direct integration;
others subtract the better-conditioned pair of log-tail jets. A separately
computed log interval width preserves narrow intervals even when adding that
width to a large log boundary would round it away.

Gauss-Legendre rules double from 16 to at most 256 nodes. Convergence is checked
on every log-probability Taylor coefficient, with target relative difference
2e-10 (scale 1 + absolute coefficient). The last two rules must agree within
1e-7; otherwise the likelihood wrapper rejects that numerical evaluation.
Far-tail nodes more than 650 log units below the scaled integrand are discarded
to avoid overflow in irrelevant exponential terms. Scalar one-sided
probabilities use the existing stable gamma/log-normal CDF implementation.

In super-exponential tails the integration width can be smaller than one ulp
of the log boundary. When Q*z > 30 and Q^2*exp(-Q*z) < 1e-12, use the
boundary expansion log f_T(boundary) - log(abs(d log f_T/dT)) with its full
Taylor jet. Its relative error is O(Q^2*exp(-Q*z)), below the quadrature
accuracy threshold. This preserves finite log probabilities and derivatives
where directly forming quadrature nodes would lose the tail width.

## Deviance and residuals

Saturate over mean only, keeping sigma and Q fixed (equivalently sigma and r).
For exact observations the saturated Prentice location is log y. One-sided
censoring has saturated probability one. For finite intervals with positive
lower bound l, set D = log(h/l), d = Q D/sigma. Equating the log-response
boundary densities gives

    a_sat = log(l) + (sigma/Q) log(expm1(d)/d).

At Q = 0 this is (log l + log h)/2. A power series for
log(expm1(d)/d)/d avoids cancellation at zero. The arithmetic mean at this
location is exp(a_sat + C). Deviance is twice the weighted difference between
saturated and fitted observed log likelihood. Signed residuals compare the
fitted location with a_sat; left and right censoring have negative and positive
sign respectively. Full-support intervals have zero deviance.

Null deviance retains fitted sigma, r and mean offsets and optimizes a single
mean intercept where present. Response/Pearson residuals require exact
observations; Pearson residuals additionally require r > 1/2.

## Fitting and diagnostics

Use the shared generalized-gamma coefficient tensor contractions and mgcv
1.9-4 NCV/QNCV corrections. NCV transports observation indices and restores
both censoring bounds inside every likelihood callback, including reordered
or repeated prediction rows. Supplied eta already contains offsets.

Native CDF, quantile and simulation callbacks describe the uncensored law.
The gratia diagnostics factory additionally supplies observed bounds, underlying
response moments and optional known censoring rules for simulated references.
The family remains a general.family, which mgcv 1.9-4 does not fit with bam.

For the REML/BFGS trace calculation, mgcv 1.9-4 drops matrix dimensions
when a predictor contains only an intercept. In that case the shared family
helper obtains full Hessian derivatives and contracts them with the inverse
penalized Hessian locally. Both generalized-gamma constructors use this route.
