# Conway–Maxwell–Poisson mean and dispersion predictors

The family uses `eta1 = log(mu)` and `eta2 = log(nu)`, with

    log p(y) = y a - nu log(y!) - A(a, nu),
    A(a, nu) = log sum_{j>=0} exp(j a - nu log(j!)),
    A_a(a, nu) = mu.

Thus `a = log(lambda)` is implicit, not a third free parameter. The mean
parameterization, adaptive normalizer and centered cumulants are the same
as in `cmp()`. See the CMP numerical notes in this directory.

For a bivariate Taylor perturbation `(h1, h2)`, substitute
`dmu = mu * (exp(h1) - 1)` and `dnu = nu * (exp(h2) - 1)` in the implicit
mean equation. Four formal Newton updates with the constant Jacobian
`A_aa = Var(Y)` determine `da` through total degree four. Substitution in
the log likelihood gives all predictor derivatives through fourth order.
Centered joint moments of `(Y, -log(Y!))` through degree five supply the
normalizer derivatives needed for this inversion. Coefficients include
factorial divisors, removed when packing derivatives for `gamlss.gH`.

REML uses the coefficient score and Hessian and their third/fourth-order
contractions. Both predictors have separate coefficients. NCV returns
predictor derivatives through third order to mgcv's neighbourhood updater.
The family callback evaluates the likelihood on the updated predictors,
including both offsets. Gamma adjustment combines updated and full-fit
log likelihoods on the actual prediction indices, including repeated or
reordered indices. The derivative uses the corresponding predictor
sensitivities; this is tested against finite differences in log smoothing
parameters, as are the REML gradient and Hessian.

Prior weights multiply all log likelihood derivatives. Zero-weight rows
have zero contribution, including when their predictors overflow.
Dispersion and mean remain positive without fitted bounds. Extreme
parameters may exhaust the numerical summation or inversion controls.

Deviance is twice the weighted log likelihood difference from the
saturated mean at each fitted nu; it is not a joint saturation over nu.
Null deviance holds fitted nu and mean offsets fixed and fits an intercept
in the mean predictor (or uses zero if that predictor has no intercept).
