# Student-t location, scale and degrees of freedom

`scatlss(min.df = 3)` defines three predictors, in this order:

\[
\mu=u,\qquad \sigma=\exp(t),\qquad \nu=b+\exp(v),\quad b=\texttt{min.df}>2.
\]

`min.df` must be a finite scalar. Each predictor has its own coefficient
block, smooths and offsets. `~ 1` estimates a constant parameter when no offset
is present. The response is an ordinary finite numeric vector, so subsetting
requires no response attributes or external row-indexed parameter state.

For \(Y=\mu+\sigma T_\nu\),

\[
E(Y)=\mu,\qquad \operatorname{Var}(Y)=\sigma^2\nu/(\nu-2).
\]

Scale is **not** the response standard deviation. The model is symmetric;
`shash()` instead supplies a different four-parameter sinh-arcsinh distribution
including skewness. `scat()` estimates or fixes global scale and df separately
from its location predictor. Matching parameters give matching likelihoods;
REML estimates need not match different parameter treatments.

## Likelihood and stable analytic derivatives

With \(z=(y-u)\exp(-t)\), the normalized log likelihood is

\[
\ell=\log\Gamma((\nu+1)/2)-\log\Gamma(\nu/2)
 -\tfrac12\log(\pi\nu)-t-\tfrac12(\nu+1)\log(1+z^2/\nu).
\]

Finite nonnegative likelihood weights multiply this expression and every
partial derivative. Zero-weight rows are neutralized before evaluation.
Let \(h=1/\nu\) and \(F(a)=\log(1+a)/a\), with \(F(0)=1\). Then the last term is
\(-\tfrac12(1+h)z^2F(z^2h)\), avoiding multiplication of a vanishing logarithm
by enormous df. `F` uses its 41-term Taylor series for `a <= 0.25` and
`log1p(a)/a` otherwise. Reciprocal df is evaluated as
`exp(-v)/(1 + b*exp(-v))` for positive `v`, and `1/(b + exp(v))` otherwise.

For df >= 40 the gamma-ratio normalization uses the Stirling ratio expansion

\[
-\tfrac12\log(2\pi)-\frac h4+\frac{h^3}{24}-\frac{h^5}{20}
+\frac{17h^7}{112}-\frac{31h^9}{36}+\frac{691h^{11}}{88}
-\frac{5461h^{13}}{52}+\frac{929569h^{15}}{480}.
\]

Below 40 it uses log-gamma and polygamma functions directly. The expansion
retains finite-df dependence; it is not a Gaussian substitution or df cap.
`inst/validation/scatlss.R` reproduces the coefficients by expanding the
Stirling corrections. Tests cover both sides of both numerical transitions.

The existing distribution-independent three-variable Taylor algebra produces
all analytic mixed predictor partials through total order four, including link
transformations. Coefficients are derivatives divided by multi-factorials.
They are repacked in mgcv's lexicographic symmetric-tensor order before
`gamlss.gH()` assembles coefficient gradients/Hessians and smoothing traces.
For trace-only requests with single-coefficient predictors, the full Hessian
derivative is contracted locally, avoiding mgcv 1.9-4's dimension dropping.
Fitting never uses numerical differentiation.

Location initialization uses weighted penalized regression; residual spread
initializes scale with a Student-t variance adjustment. Initial df is
`min.df + max(2, 5 - min.df)`. All three initial regressions subtract offsets.
Invalid or unrepresentable parameter trials return negative infinity for step control; failed
NCV solves or non-finite criteria/gradients produce explicit errors.
Extreme predictors and residuals can exceed floating-point range. Very large
df is weakly identified, and scale and tail weight can be confounded; singular
or poorly conditioned fits remain possible. No finite upper df bound is imposed.

## NCV, QNCV and mgcv 1.9-4

Let \(H=S-\ell_{\beta\beta}\) be the full penalized negative Hessian.
For deletion set \(A\), define

\[
H_{-A}=H+\sum_{i\in A}\ell_{i,\beta\beta},\qquad
s_A=\sum_{i\in A}\ell_{i,\beta},\qquad
\delta_A=H_{-A}^{-1}s_A,\quad \beta_{-A}^{(1)}=\beta-\delta_A.
\]

These are one-step deletion estimates, not exact refits. All Hessian entries
retain their signs: Student-t location curvature can be non-positive in the
negative log likelihood, and mixed curvature can have either sign. The
compiled mgcv conjugate-gradient route uses the signed Hessian, avoiding the
1.9-4 Cholesky deletion discrepancy already covered by generalized-gamma tests.
CG iteration counts are distinguished from Cholesky failure counts; exhausting
the CG iteration limit fails explicitly. Negative derivative requests for
coefficient jackknife output use local direct solves because mgcv 1.9-4's CG
helper leaves that output buffer uninitialized.

Scores and gradients are assembled locally from prediction rows in `nei$d`,
preserving their order and repetitions, offsets and weights. If \(\ell_i^0\)
is the full-data contribution and \(\ell_i^{cv}\) its deletion prediction,
NCV is \(-\gamma\sum_i\ell_i^{cv}-(1-\gamma)\sum_i\ell_i^0\).
QNCV replaces each held-out contribution by
\(\ell_i^0+\ell_{i,\eta}^T\Delta_i+
\tfrac12\Delta_i^T\ell_{i,\eta\eta}\Delta_i\), with
\(\Delta_i=\eta_i^{cv}-\eta_i\). Analytic smoothing gradients differentiate
both deletion predictions and full-data derivatives. `Vg` is the cross-product
of prediction-row gradient contributions. Jackknife output uses the same
signed deletion updates.

The response-attribute-loss issue guarded by the package's extended-family
helper does not apply here: the numeric response subsets directly. No
blanket extended-family guard is applied to `scatlss()`. End-to-end default and
grouped NCV/QNCV, and NCV requested through `nei` during REML, are regression tested
on installed mgcv 1.9-4. No installed mgcv code is patched and there is no public
bypass. A separate mgcv 1.9-4 `nanei()` defect can corrupt indices or cumulative
boundaries when missing observations adjust explicit neighbourhoods.
The constructor uses the shared `.inlaws_ncv_family()` installer, which
prepends `.inlaws_check_ncv_nei()` before its existing BFGS optimizer setup.
It checks adjusted structure on the installed version without presuming a
future fix, rejects unsafe NCV/QNCV and REML-with-`nei` requests, and permits
valid whole-group omissions. It does not reconstruct neighbourhoods. Filter
all relevant response, predictor, offset and weight variables after subset
selection, then rebuild `nei` against the resulting rows. Default leave-one-out
after omission, complete-case grouped fits and REML without `nei` remain
available. Incomplete lists may fail inside mgcv before the family hook.
For a supplied `G`, neighbourhoods must refer to its retained rows; model-frame
adjustment is not repeated. See `ncv-neighbourhood-audit.md` for hook limits.
Guard tests are rejection coverage, separate from successful end-to-end fits.

An additional mgcv 1.9-4 empty-penalty reparameterization defect prevents
jackknife covariance for fully unpenalized models. Such NCV/QNCV models can
evaluate their criterion using explicit `nei` with `jackknife = FALSE`.
The affected covariance request errors explicitly; ordinary fully unpenalized
REML models remain supported. There are no smoothing parameters to select
in these fully parametric fits.

These compatibility measures make no claim about untested future
mgcv releases. REML uses outer Newton or BFGS; NCV/QNCV use BFGS. Other
criteria, `bam()`, `gamm()` and shared predictor coefficients are unsupported.

## Predictions, diagnostics and validation

Fitted and response predictions contain `mu`, `sigma`, `nu`. Link predictions
contain `u`, `t`, `v`. Marginal delta-method standard errors use multipliers
`1`, `sigma`, `nu - min.df`; estimation errors are not predictive deviations.
Default inference conditions on smoothing parameters; `unconditional = TRUE`
requests mgcv's available smoothing correction. Prediction does not average
over coefficient or random-effect uncertainty.

Simulation, CDFs and quantiles use `rt`, `pt` and `qt`, in the requested tail
and log-probability convention. Callbacks ignore likelihood weights and external
scale. Diagnostics provide the complete Student-t mean and variance without
requiring gratia. Deviance holds scale and df fixed, saturating location at
`y`. Its log1p expression retains tiny deviances, including when `z^2/nu`
underflows; a logarithmic tail expression avoids squaring huge residuals.
Response residuals are `y-mu`, and Pearson residuals multiply by
`sqrt(weight/variance)`. Null deviance optimizes a constant location plus its
offset, retaining fitted scale and df; without an intercept it uses the offset.

Tests separately check analytic derivatives against independent tensor-product
finite differences of `stats::dt()`, differences of all lower analytic partials,
coefficient gradients/Hessians, REML score derivatives, and probability/moment
references. Model tests cover all supported optimizers, global and varying
parameters, offsets, weights, missing data, subsets and serialization. NCV
references explicitly assemble deleted Hessians and solve them without native
deletion updates; independently evaluated Student-t densities and numerical
score derivatives check NCV/QNCV and `Vg`. Exact deletion refits provide an
additional approximation assessment, not an equality requirement. Guarded or
skipped end-to-end routes must be reported separately from numerical tests.
