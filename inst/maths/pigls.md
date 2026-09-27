# Poisson inverse Gaussian mean and dispersion predictors

## Model and likelihood

Use the same mixture as `pig()`: conditional on `Z`, counts are Poisson with
mean `mu*Z`, and `Z` is inverse Gaussian with mean one and variance `theta`.
Both parameters now have predictors:

```
mu_i = exp(eta_i1), theta_i = exp(eta_i2)
eta_ij = X_j[i, ] beta_j + offset_ij
E(Y_i) = mu_i, Var(Y_i) = mu_i + theta_i * mu_i^2.
```

The exact log probability is the half-integer Bessel recurrence documented
in [pig.md](pig.md). Constants are retained. Prior weights multiply each
observation's log likelihood and derivatives, including fractional weights;
zero-weight rows are replaced by harmless values before numerical evaluation.
No mutable family parameter environment stores fitted dispersion.

## Derivatives and REML

The shared recurrence propagates bivariate Taylor polynomials through total
order four. The natural-mean derivative path remains the default for `pig()`.
For `pigls()`, its mean polynomial is instead

```
mu(eta1 + d1) = mu * sum(d1^k/k!, k=0,...,4).
```

The dispersion polynomial is `theta*sum(d2^k/k!, k=0,...,4)` for both paths.
Thus every mixed derivative in `(eta1, eta2)`, including pure third/fourth
dispersion derivatives, comes from the same analytic recurrence. Coefficients
include inverse factorials; unpacking multiplies by `a!*b!` for derivative
order `(a,b)` and uses mgcv's lexicographic predictor-index order.

`mgcv::gamlss.gH` contracts these observation derivatives into coefficient
gradients, Hessians, higher derivatives and smoothing-parameter trace terms.
The sandwich calculation uses weighted observation scores. Infeasible trial
predictors return negative-infinite likelihood, allowing the optimizer to
backtrack. There is no fitted dispersion bound. Starting coefficients come
from weighted, lightly regularized regressions of log(count+1/6) and a
method-of-moments dispersion target; only the latter starting target is
restricted to a moderate positive interval.

The two coefficient blocks are separate. Reusing covariates in both formulae
is supported, but shared coefficients are rejected. REML is supported with
outer Newton and BFGS. A constant second predictor still represents a
regression coefficient, so its REML treatment differs from the extra family
parameter in `pig()`; equality of those REML estimates is not an invariant.
With fully parametric identifiable predictors, joint likelihood optimization
provides an independent check of the fitted coefficients.

## NCV

The likelihood supplies observation-level first, second and third predictor
derivatives, including for supplied held-out predictors whose model matrix
has no `lpi` attribute. Supplied predictors already contain offsets; they are
never added a second time.

The wrapper obtains neighbourhood updates from `mgcv:::gamlss.ncv` at gamma=1.
For prediction indices `ix=nei$d`, it then evaluates the actual PIG likelihood
at the returned held-out predictors, and reconstructs the score as

```
NCV = -gamma * sum(l_cv[ix]) - (1-gamma) * sum(l_full[ix]).
```

Both terms retain offsets, weights and the exact prediction-row order,
including repeated indices and prediction subsets of a different length from
the training data. The smoothing gradient applies the chain rule using
`deta.cv` for the first term and `X_j %*% db_j` for the second. `Vg` is the
cross-product of these individual score-gradient contributions. The returned
CV predictor attributes retain the original offsets. QNCV is explicitly
rejected. NCV uses BFGS smoothing selection.

## Predictions, residuals and callbacks

Response predictions contain mean and dispersion in that order. The family
prediction callback supplies both columns and delta-method standard errors;
mgcv may drop column labels when assembling prediction blocks. Fitted-value
columns are named `mu` and `theta`.

Deviance is `2*w*(l_sat(y,theta)-l(y,mu,theta))`, where the saturated mean
is found at each observation's fitted dispersion. For zero counts its value
is zero. Unlike Poisson, the saturated mean generally differs from y.
Null deviance optimizes a constant mean predictor while retaining dispersion
and mean offsets; without a mean intercept it uses the offset alone.
Response residuals are `y-mu`; Pearson residuals divide the weighted response
residuals by `sqrt(mu+theta*mu^2)`; deviance residuals use the sign of `y-mu`.

The simulation, CDF and quantile callbacks accept a two-column parameter
matrix and reuse the PIG numerical helpers with row-specific dispersion.
A single parameter row can be recycled over probability/quantile arguments.
The existing PIG boundary conventions, missing-value propagation, scalar or
equal-length recycling, and numerical limits apply. Simulation uses the
inverse Gaussian mixture directly. The helpers retain their original scalar
dispersion behavior for `pig()`.

## Validation

Tests cover each mixed predictor derivative through order four, coefficient
contractions and the sandwich matrix, independent likelihood maximization,
Newton/BFGS agreement, offsets, starts, weighted replication, missing rows,
serialization, parameter recovery, prediction covariance, conditional null
deviance and distribution moments. Vector-dispersion numerical results are
compared with scalar-dispersion calculations row by row.

NCV tests compare implicit and explicit leave-one-out fits, exercise grouped
and reordered neighbourhoods, and check the score directly from held-out
predictors. Smoothing gradients are compared with finite differences of
fully converged inner fits, including nonunit gamma and repeated prediction
rows. Exact neighbourhood refits check the approximation's predictors.
Existing `pig()` tests continue to validate its original numerical contracts.
