# Negative binomial mean and size regression

The distribution has mean mu > 0, size t = theta > 0, and variance
mu + mu^2/t. The predictors are g(mu) and log(t). The first link can be
log, identity, or square root; the latter two have positive predictor domains.
Every likelihood contribution and derivative below is multiplied by its prior
weight before coefficient derivatives are assembled.

## Natural-parameter derivatives

For an integer observation y >= 0, put z = mu + t. The log likelihood is

    l = lgamma(y+t) - lgamma(t) - lgamma(y+1)
        + t log(t) + y log(mu) - (y+t) log(z).

Likelihood values in the implementation come from R's `dnbinom(log = TRUE)`.
Write l[p,q] for p derivatives in mu and q in t, and (a)_q for the rising
factorial a(a+1)...(a+q-1). For p >= 1,

    l[p,0] = (-1)^(p-1) (p-1)! { y/mu^p - (y+t)/z^p }.

In particular, the score is evaluated as t/z * (y-mu)/mu. For p,q >= 1,

    l[p,q] = (-1)^(p+q) (p-1)! {
                 (p)_q (y-mu)/z^(p+q) + (p-1)_q/z^(p+q-1) }.

For the pure size derivatives, write

    Delta_q = psi^(q-1)(t+y) - psi^(q-1)(t).

Then

    l[0,1] = Delta_1 - log1p(mu/t) + (mu-y)/z

and, for q >= 2,

    l[0,q] = Delta_q + (-1)^q (q-2)! { 1/t^(q-1) - 1/z^(q-1) }
                       + (-1)^q (q-1)! (y-mu)/z^q.

These formulae give all 2, 3, 4, and 5 distinct derivatives of orders one
through four. Their columns are ordered with an increasing number of size
derivatives, matching `mgcv::trind.generator(2)`.

For small counts, polygamma differences are evaluated using the exact recurrence

    Delta_q = (-1)^(q-1) (q-1)! sum_{k=0}^{y-1} (t+k)^(-q).

This avoids subtracting two almost equal special-function evaluations and gives
exactly zero when y=0. Direct special-function differences are used for larger
counts outside the near-Poisson region.

## Near-Poisson evaluation

Direct evaluation of the pure size derivatives loses precision as t increases.
Expand the finite-product likelihood instead:

    l = log Pois(y; mu) + sum_{r>=1} C_r / t^r,
    C_r = (-1)^(r+1) { S_r/r - y mu^r/r + mu^(r+1)/(r+1) },
    S_r = sum_{k=0}^{y-1} k^r.

Consequently

    l[0,q] = sum_{r>=1} (-1)^q (r)_q C_r / t^(r+q).

Use twelve terms when t > 100 max(1,y,mu). The expansion ratios are then
less than 0.01, including for the finitely many polynomial factors introduced
by differentiation through order four. This is a numerical evaluation of the
same NB likelihood, not a switch to a Poisson family or a cap on size.

Compute scaled sums T_r = S_r/t^r, rather than large powers of y, using
the telescoping identity for (k+1)^(r+1) - k^(r+1):

    T_0 = y,
    T_r = { y(y/t)^r - sum_{j=0}^{r-1} choose(r+1,j) T_j/t^(r-j) }/(r+1).

This keeps the score, curvature, and higher derivatives nonzero and mutually
consistent close to the Poisson limit. As with other location/scale families,
identifiability of size itself is lost at that limit.

## Links and coefficient derivatives

`mgcv::gamlss.etamu()` applies the chain rule through fourth order, using the
inverse first link derivative and the second through fourth forward link
derivatives. `mgcv::gamlss.gH()` contracts these observation-level arrays with
the predictor model matrices and the derivatives of coefficients with respect
to log smoothing parameters. The log inverse links are exact exponentials,
without the small-value flooring used by `stats::make.link("log")`.

The likelihood callback follows mgcv's conventions:

| `deriv` | Returned information |
| --- | --- |
| 0 | Total and observation log likelihoods |
| 1 | Coefficient score and Hessian |
| 2 | First Hessian derivative traces using third derivatives |
| 3 | Full first Hessian derivatives |
| 4 | Second Hessian derivative traces using fourth derivatives |

EFS only needs the first two likelihood derivatives. Outer BFGS also needs
third derivatives; outer Newton needs all four. The finished family advertises
`available.derivs = 2` and retains all three optimizer choices.

In mgcv 1.9-3 and 1.9-4, the fourth-order trace accumulator in `gamlss.gH()`
overwrites overlapping coefficient columns for shared predictors. For that
case only, `nbls` computes the contraction explicitly. If V is the inverse
penalized Hessian, X_i is predictor i's matrix, and a_i,k and b_i,kl are
the first and second predictor derivatives in smoothing-parameter directions,
the contribution from predictor indices i,j is

    sum_n diag(X_i V_ij X_j')[n] * {
      sum_q l_ijq[n] b_q,kl[n] +
      sum_q sum_s l_ijqs[n] a_q,k[n] a_s,l[n] }.

Off-diagonal predictor pairs contribute twice. All other contractions use
mgcv's exported helper directly. Tests compare the shared and unshared
contractions with independent finite differences.

## Diagnostics and interpretation

Deviance compares l(y; mu,t) with l(y; y,t), holding each observation's size
fixed. Null deviance minimizes this conditional deviance over a constant mean
predictor while preserving mean offsets and the fitted size vector. Pearson
residuals are sqrt(weight) (y-mu)/sqrt(mu+mu^2/t); response residuals are y-mu.

Observation weights multiply log likelihoods. Sandwich filling uses the
cross-product of these weighted individual coefficient scores, so its weights
are squared. Simulation and quantile callbacks describe individual responses
and do not scale their distributions by prior weights.

With an intercept-only size predictor, coefficient estimation maximizes the
joint penalized likelihood, but REML integrates over that coefficient along
with the other regression coefficients. `mgcv::nb()` instead treats size as
an extra family parameter. Exact agreement between their REML smoothing
parameters or resulting fits is therefore not an acceptance criterion.

The derivative checks differentiate the preceding analytic order instead of
taking repeated fourth differences of likelihood values. Coefficient and
smoothing-parameter contractions, joint intercept-only maximum likelihood,
weight replication, and the three optimizers provide independent checks.

For compatibility with mgcv 1.9-3, the family's `preinitialize` hook also
supplies `G$offxset` as an alias of `G$offset`: that version's EFS dispatch
uses the misspelled name. This preserves offsets in fitting and fitted values
without changing mgcv or persisting training offsets in family callbacks.
