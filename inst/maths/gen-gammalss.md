# Generalized gamma mean, scale and transformed shape

The response is positive. The three regression predictors are

\[
 u=\log m,\quad t=\log\sigma,\quad v=\log r,\qquad
 r=1+\sigma Q,\quad Q=\operatorname{expm1}(v)e^{-t}.
\]

Here \(m=E[Y\mid x]\) is the arithmetic mean, not Prentice's log-response
location. Predictors model `mu`, `sigma`, `r` in that order, with separate
coefficient blocks and offsets. Penalties act on `u`, `t`, `v`. A constant
third predictor holds `r` constant; `Q` can vary with scale. All real predictor
triples define a finite-mean distribution mathematically, subject to finite
floating-point evaluation. There is no truncation at a fitted shape bound.

## Density and moments

For \(Q\ne0\), let \(k=Q^{-2}\), \(G\sim\mathrm{Gamma}(k,1)\), and

\[
 C(\sigma,Q)=\log\Gamma(k+\sigma/Q)-\log\Gamma(k)
               -(\sigma/Q)\log k,\qquad a=u-C(\sigma,Q).
\]

Then \(Y=\exp\{a+\sigma\log(G/k)/Q\}\) has mean \(m\).
The first moment requires \(k+\sigma/Q=kr>0\), exactly `r > 0`.
For negative Q higher moments need not
exist: in general the positive moment of order \(d\) exists iff
\(1+d\sigma Q>0\). In particular variance requires `r > 0.5`, and then

\[
 \mathrm{Var}(Y)=m^2\operatorname{expm1}
             \{C(2\sigma,Q)-2C(\sigma,Q)\}.
\]

With \(z=(\log y-a)/\sigma\), define

\[
 A(Q)=\log|Q|-\log\Gamma(k)+k\log k-k,\qquad
 E(x)=\frac{\exp(x)-1-x}{x^2}.
\]

The log density, including its original-response Jacobian, is

\[
 \ell=A(Q)-t-\log y-z^2E(Qz).
\]

At Q=0, `C = sigma^2/2`, `A = -log(2*pi)/2`, `E(0) = 1/2`:
this is exactly the log-normal with meanlog `log(mu)-sigma^2/2`.
Gamma occurs at `Q = sigma`; Weibull occurs at `Q = 1`.

## Stable derivatives

All mixed predictor derivatives through total order four are obtained by
truncated multivariate Taylor arithmetic, with coefficients equal to the
corresponding derivative divided by its multi-factorial. No finite differences
are used in fitting. Symmetric derivative tensors use mgcv's lexicographic
packing. Multiplication tables are cached, not observation values.

Near Q=0 we use Stirling's expansion, with
\(b_j=B_{2j}/\{2j(2j-1)\}\), through eight terms:

\[
 A=-\tfrac12\log(2\pi)-\sum_{j=1}^{8}b_j Q^{4j-2},
\]
\[
 C=\sigma^2 F(h)-\tfrac12v+
       \sum_{j=1}^{8}b_j Q^{4j-2}\operatorname{expm1}\{-(2j-1)v\},
 \quad h=\operatorname{expm1}(v),
\]
\[
 F(h)=\frac{(1+h)\log(1+h)-h}{h^2}
     =\sum_{j=0}^{\infty}\frac{(-h)^j}{(j+1)(j+2)}.
\]

The correction expansion requires both gamma arguments to be at least 20.
Otherwise analytic derivatives of log gamma are composed directly. `F` uses
61 series terms for `abs(h) <= 0.25`; `E` uses 33 terms for `abs(Q*z) <= 1`.
These series retain mixed shape derivatives even at exactly Q=0. They do not
replace a neighbourhood of zero by a shape-independent log-normal density.

## Probability callbacks and diagnostics

CDFs use gamma lower tails for positive Q and upper tails for negative Q;
quantiles reverse the corresponding transformation. Gamma probabilities and
quantiles are evaluated in their requested tail, including log probabilities.
At exactly Q=0 use the log-normal callbacks. For `abs(Q) <= 1e-5`, when
`abs(Q)*(1+abs(z))^3 <= 0.01`, a second-order expansion avoids roundoff in
huge gamma arguments. For the standardized log response it is

\[
 F_W(z)=\Phi(z)+\phi(z)\left[
 Q\frac{z^2+2}{6}-Q^2\frac{z^5+2z^3+6z}{72}\right]+O(Q^3),
\]
\[
 F_W^{-1}(p)=z-Q\frac{z^2+2}{6}+Q^2\frac{z^3+5z}{36}+O(Q^3),
 \quad z=\Phi^{-1}(p).
\]

Other tails use gamma functions. Samplers invert uniforms. The native callbacks
and diagnostic factory describe the same full response distribution, ignoring
case weights and external scale. Response residuals subtract the mean; Pearson
residuals require finite variance. The factory reports undefined variance
without disabling CDF, quantile or simulation diagnostics.

Deviance saturates only the mean, holding scale and transformed shape fixed;
its saturated log-response location is `log(y)`. Null deviance fits a constant
log-mean predictor, retaining mean offsets and fitted scale and shape. Without
a mean intercept it evaluates the zero coefficient predictor plus offsets.

## NCV and QNCV under mgcv 1.9-4

The family uses mgcv's compiled conjugate-gradient neighbourhood deletion
updates with the full signed Hessian. Its Cholesky route can produce different
updates for three predictors with negative cross-curvature; independent direct
matrix solves in the tests guard against this discrepancy. CG iteration counts
are distinct from the Cholesky failure count expected by `gam.fit5`. Reaching
the CG iteration limit is an explicit failure.

The 1.9-4 CG helper does not initialize its output buffer for negative derivative
requests (coefficient jackknife perturbations). Those requests are evaluated
locally by solving the signed penalized deletion Hessian. This also avoids a
native crash in mgcv's post-fit jackknife calculations. This path uses dense
coefficient matrices and is not a large-data performance guarantee.

NCV scores and derivatives are assembled locally, preserving offsets, grouped
neighbourhoods, reordered/repeated prediction rows, prediction subsets, and
`gamma`. QNCV uses the quadratic Taylor approximation to each held-out log
likelihood and its analytic smoothing-parameter gradient. These compatibility
measures target mgcv 1.9-4 and should be reassessed after upstream NCV fixes,
expected in later 1.9-5/1.9-6 releases; they are not version-based assumptions
about unreleased implementations. No changes to other families are made here.

## Validation

The tests check all mixed derivatives through order four against differences
of lower derivatives, independent gamma-change-of-variables densities, exact
log-normal/gamma/Weibull/exponential submodels, probability inversions in both
tails and on log scales, simulated response moments, and zero-weight rows.
Model tests cover optimized REML/NCV/QNCV, Newton and BFGS REML, offsets,
weights/replicated data, independent unpenalized optimization, predictions and
standard errors, omission and serialization. Fixed-penalty tests compare REML
score gradients/Hessians and NCV/QNCV gradients with numerical derivatives,
compare deletion predictions with independent matrix solves and actual refits,
and verify coefficient jackknife updates.

`bam()` rejects general families upstream: ordinary/discrete fREML and NCV
attempts are tested as unsupported. No mgcv functions are patched.
