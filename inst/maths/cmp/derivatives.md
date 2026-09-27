# CMP reference implementation

The response is a nonnegative integer, with

    ell(y; a, nu) = y a - nu log(y!) - A(a, nu),
    A(a, nu) = log sum_j exp(j a - nu log(j!)),
    a = log(lambda), nu > 0.

The user models log(mu), with mu = A_a. The family optimises theta = log(nu).
The separate GLM likelihood scale is fixed at one. Likelihood weights multiply
individual contributions; they do not change the response distribution.

## Summation and tail bounds

Summation is centred at m = floor(exp(a/nu)), a mode. Weights are scaled by
the mass at m, and support expands in both directions. For each p + q <= 5,
we bound omitted sums of

    w_j |j-m|^p |log(j!)-log(m!)|^q.

For the upper tail beginning at k > m, k >= 3, use

    log(j!)-log(m!) <= (j-m) log(j),
    log(j) <= log(k) j/k,                  j >= k.

The weight ratio is at most exp(a)/(k+1)^nu. The envelope's successive ratio
is at most (1 + 1/(k-m))^(p+q) (1 + 1/k)^q. Their product r gives the bound

    w_k (k-m)^(p+q) log(k)^q / (1-r),      r < 1.

For the lower tail beginning at k < m, use

    log(m!)-log(j!) <= (m-j) log(m).

The corresponding ratio is exp(-a) k^nu (1 + 1/(m-k))^(p+q).
The envelope starts at w_k (m-k)^(p+q) log(m)^q. A tail consisting only of
zero has no further terms. If m=1, factorial-log differences below the mode
are exactly zero. A ratio >= 1 requires expansion, not acceptance.

Each combined tail bound must be less than sum_tol times the computed
absolute-moment sum. The zeroth moment checks normalisation. This conservative
criterion also controls the moment inputs, not just the probability mass.
Subsequent centring and cumulant cancellation can amplify relative error in
near-degenerate cases; the tolerance is not a promise of relative accuracy
for a cumulant arbitrarily close to zero.

max_terms limits evaluated support length, rather than the largest possible
count. Exceeding a limit or obtaining nonfinite derivatives is an explicit
error. There is no asymptotic approximation or silent truncation fallback.

## Mean inversion

Solve A_a(a,nu) = mu on the a scale. Its derivative is A_aa = Var(Y) > 0.
Maintain a bracket as evaluations become available and use safeguarded Newton
steps; fall back to bisection within a finite bracket. Convergence requires
absolute mean error <= mean_tol * mu. The mu=0 boundary is a point mass at
zero and is evaluated without finite-a inversion.

## Cumulants and analytic implicit differentiation

Let X=Y-E(Y) and T=-(log(Y!)-E(log(Y!))). Form the polynomial

    M(u,v) = sum_{2 <= i+j <= 5} E(X^i T^j) u^i v^j / (i! j!).

The degree-five expansion of log E exp(uX+vT) is M - M^2/2: M has no
constant or linear term. Thus the natural-parameter Taylor expansion is

    A(a+u,nu+v) = A(a,nu) + mu u - E(log(Y!)) v + M - M^2/2.

Fifth-order cumulants are needed because fourth-order differentiation of
the inverse mean map differentiates A_a, not A itself.

All polynomial coefficients include their factorial divisors. Work in formal
increments h of mu and t of theta, truncated at total degree four. Set

    v = nu (t + t^2/2 + t^3/6 + t^4/24).

Starting with u=0, apply four formal corrections

    u <- u + [h - {A_a(a+u,nu+v)-A_a(a,nu)}] / Var(Y).

Each correction determines one further homogeneous degree. This is analytic
Taylor arithmetic, not finite differencing. Compose the likelihood with
these polynomials and multiply coefficients by factorials when returning
ordinary partial derivatives to mgcv. Standard mgcv link machinery converts
mean-scale derivatives to log-mean-scale derivatives.

For example, ell_mu = (y-mu)/V and E(-2 ell_mumu) = 2/V.
All likelihood derivatives involving mu are affine in y, so substituting
mu for y gives expected observed curvature and higher derivatives. This
substitution is *not* used for pure dispersion derivatives involving log(y!).

## Saturation, deviance and parameter state

For fixed nu the saturated positive mean is y. At y=0 the saturated
log likelihood and all dispersion derivatives are zero. The deviance is

    D_i = 2 w_i [ell(y_i; mu=y_i, nu) - ell(y_i; mu_i, nu)].

Pure theta derivatives include derivatives of the saturated term; mixed
mean/theta derivatives do not. The ls callback returns the weighted saturated
log likelihood and its first two theta derivatives, including per-observation
first derivatives. Consequently D/2-ls equals the negative weighted log
likelihood. AIC is computed from that likelihood independently.

Caches are bounded to the most recent mean evaluation and saturated response
vector. They are local to each family. preinitialize clones the family before
fitting, retaining mgcv-added components while preventing reuse of a family
object from mutating an earlier fit's dispersion.

## References

Huang (2017), Mean-parametrized Conway-Maxwell-Poisson regression models for
dispersed counts, Statistical Modelling 17, 359–380.
https://doi.org/10.1177/1471082X17697749

Raim and Sellers (2022), COMPoissonReg: Usage, the Normalizing Constant, and
Other Computational Details, US Census Bureau RRC2022-01.
https://www.census.gov/library/working-papers/2022/adrm/RRC2022-01.html

The implementation uses original R Taylor arithmetic and tail bounds;
it does not copy mgcv or another package's numerical implementation.
