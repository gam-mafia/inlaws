# Cumulative-link ordered categorical models for mgcv

This document derives the likelihood and derivatives implemented by
`inlaws::cumulative_link()`. It follows the purpose of the package's censored
log-normal mathematical notes: make the parameterization, calculus and mapping
to the fitting code explicit. It is Markdown with LaTeX mathematics; no LyX
source or generated PDF is required.

## 1. Model, links and identification

Let $Y_i\in\{1,\ldots,K\}$ and introduce a latent variable

$$
Z_i=\mu_i+\epsilon_i,\qquad \mu_i=X_i\beta,
$$

where the error has a fixed, unit-scale CDF $F$. The model matrix can include
mgcv smooth bases and random-effect terms. Define

$$
-\infty=\alpha_0<\alpha_1<\cdots<\alpha_{K-1}<\alpha_K=\infty,
\qquad Y_i=j\ \Longleftrightarrow\ \alpha_{j-1}<Z_i\le\alpha_j.
$$

Thus

$$
P(Y_i\le j)=F(\alpha_j-\mu_i),\qquad
p_{ij}=F(\alpha_j-\mu_i)-F(\alpha_{j-1}-\mu_i).
$$

The word **location** matters: the Cauchy distribution has no mean, and the
extreme-value distributions below are not centered to have zero mean.

| Cumulative link $g=F^{-1}$ | $F(z)$ | $g(p)$ |
| --- | --- | --- |
| logit | $(1+e^{-z})^{-1}$ | $\log\{p/(1-p)\}$ |
| probit | $\Phi(z)$ | $\Phi^{-1}(p)$ |
| cloglog | $1-\exp(-e^z)$ | $\log\{-\log(1-p)\}$ |
| loglog | $\exp(-e^{-z})$ | $-\log\{-\log p\}$ |
| cauchit | $1/2+\arctan(z)/\pi$ | $\tan\{\pi(p-1/2)\}$ |

For logit, $\operatorname{logit}P(Y_i\le j)=\alpha_j-\mu_i$:
this is proportional odds, including when $\mu_i$ contains nonlinear smooths.
Other links share a location effect across thresholds but are not
proportional-odds models. Positive effects shift mass towards larger categories.
The `link` argument names the **cumulative** link; the returned mgcv family has
`link="identity"` because its single distribution parameter is $\mu_i$ itself.

Simultaneously adding a constant to all thresholds and $\mu$ leaves probabilities
unchanged. Following `mgcv::ocat()`, retain a model intercept and set

$$
\alpha_1=-1,\qquad
\alpha_j=-1+\sum_{r=1}^{j-1}e^{\theta_r},\quad 2\le j\le K-1.
$$

There are $K-2$ unconstrained log-increments. The latent scale is fixed; it is
not an extra dispersion parameter. For comparison with `MASS::polr()`, subtract
the fitted location intercept from every finite threshold. With $K=2$ there
are no free threshold increments and this is a binary regression model.

The public `theta` argument follows `ocat()`'s convention: positive increments
fix thresholds; a vector containing a negative entry supplies starting
increments via absolute values. Internally all increments are logged.
`getTheta()` returns log-increments and `getTheta(TRUE)` returns finite cut-points.

## 2. Likelihood and deviance

For one observation suppress its index and write

$$
a=\alpha_{y-1}-\mu,\quad b=\alpha_y-\mu,\quad
p=F(b)-F(a),\quad L=\log p.
$$

With a nonnegative case weight $w$, the log-likelihood contribution is $wL$.
The saturated categorical likelihood contribution is zero on the log scale, so

$$
D=-2wL.
$$

All derivatives below are derivatives of **unweighted $L$**. The implementation
multiplies them by $-2w$ to return mgcv's deviance derivatives. Zero-weight
contributions are explicitly zeroed. The smoothing penalty adds
$\tfrac12\sum_r\lambda_r\beta^\mathsf{T}S_r\beta$ to the negative
log likelihood; it does not change the observation-level formulas.

## 3. Endpoint and threshold derivatives

Let $F^{(m)}$ denote the $m$th derivative of the CDF. The endpoint derivative for
log-increment $\theta_r$ is

$$
A_{jr}=\frac{\partial\alpha_j}{\partial\theta_r}
=e^{\theta_r}\mathbf1\{r<j<K\},
\qquad
\frac{\partial^2\alpha_j}{\partial\theta_r\partial\theta_s}
=\mathbf1\{r=s\}A_{jr}.
$$

Both derivatives of infinite endpoints are defined as zero. Put
$B_r=A_{yr}$ and $A_r=A_{y-1,r}$. Since both endpoints have location derivative
$-1$, define the following ratios:

$$
r_m=\frac{\partial_\mu^m p}{p}
=\frac{(-1)^m\{F^{(m)}(b)-F^{(m)}(a)\}}p,
\quad m=1,\ldots,4,
$$

$$
s_{mr}=\frac{\partial_\mu^m\partial_{\theta_r}p}{p}
=\frac{(-1)^m\{F^{(m+1)}(b)B_r-F^{(m+1)}(a)A_r\}}p,
\quad m=0,\ldots,3,
$$

$$
\begin{split}
t_{mrs}
&=\frac{\partial_\mu^m\partial_{\theta_r}\partial_{\theta_s}p}{p}\\
&=\frac{(-1)^m}{p}\left[
 F^{(m+2)}(b)B_rB_s-F^{(m+2)}(a)A_rA_s
 +\mathbf1\{r=s\}\{F^{(m+1)}(b)B_r-F^{(m+1)}(a)A_r\}
\right],\quad m=0,1,2.
\end{split}
$$

These formulas include the second derivative of the exponential threshold
transform; omitting it gives incorrect diagonal threshold Hessian entries.

## 4. Derivatives of the log likelihood

Differentiating $L=\log p$ gives the location derivatives

$$
\begin{aligned}
L_\mu &= r_1,\\
L_{\mu\mu} &= r_2-r_1^2,\\
L_{\mu\mu\mu} &= r_3-3r_1r_2+2r_1^3,\\
L_{\mu\mu\mu\mu} &= r_4-4r_1r_3-3r_2^2+12r_1^2r_2-6r_1^4.
\end{aligned}
$$

For one threshold derivative,

$$
\begin{aligned}
L_{\theta_r} &= s_{0r},\\
L_{\mu\theta_r} &= s_{1r}-r_1s_{0r},\\
L_{\mu\mu\theta_r} &= s_{2r}-2r_1s_{1r}+(2r_1^2-r_2)s_{0r},\\
L_{\mu\mu\mu\theta_r} &= s_{3r}-3r_1s_{2r}
 +(6r_1^2-3r_2)s_{1r}
 +(-r_3+6r_1r_2-6r_1^3)s_{0r}.
\end{aligned}
$$

For two threshold derivatives, put
$C_{rs}=s_{1r}s_{0s}+s_{0r}s_{1s}$ and $P_{rs}=s_{0r}s_{0s}$. Then

$$
\begin{aligned}
L_{\theta_r\theta_s} &= t_{0rs}-P_{rs},\\
L_{\mu\theta_r\theta_s} &= t_{1rs}-r_1t_{0rs}-C_{rs}+2r_1P_{rs},\\
L_{\mu\mu\theta_r\theta_s} &= t_{2rs}-2r_1t_{1rs}
 +(2r_1^2-r_2)t_{0rs}\\
&\quad-s_{2r}s_{0s}-s_{0r}s_{2s}-2s_{1r}s_{1s}
 +4r_1C_{rs}+(2r_2-6r_1^2)P_{rs}.
\end{aligned}
$$

In `Dd()`, these become `Dmu`, `Dmu2`, `Dmu3`, `Dmu4`, `Dth`,
`Dmuth`, `Dmu2th`, `Dmu3th`, `Dth2`, `Dmuth2`, and `Dmu2th2` after
multiplication by $-2w$. Symmetric pairs are packed in mgcv's order
$(1,1),(1,2),\ldots,(1,q),(2,2),\ldots,(q,q)$, where $q=K-2$.
`level=0` supplies the working derivatives; levels 1 and 2 supply the higher
and mixed derivatives needed by the outer optimization. These are also supplied
for binary outcomes even though no threshold derivatives are required.

## 5. Distribution-specific derivatives

Write $f=F'$. Only $\log f$ and
$H=(1,f'/f,f''/f,f'''/f)$ are needed; the ratios above are formed using
$\exp(\log f-\log p)H$ at each endpoint.

For logit, let $u=F(z)$:

$$
H=(1,\ 1-2u,\ 1-6u+6u^2,\ 1-14u+36u^2-24u^3).
$$

For probit:

$$
H=(1,\ -z,\ z^2-1,\ 3z-z^3).
$$

For cloglog, let $t=e^z$, with $\log f=z-t$:

$$
H=(1,\ 1-t,\ 1-3t+t^2,\ 1-7t+6t^2-t^3).
$$

Loglog is its reflection: evaluate at $-z$ and reverse the signs of the odd
log-density-normalized derivatives $f'/f$ and $f'''/f$.

For cauchit, let $v=z/(1+z^2)$ and $r=1/(1+z^2)$:

$$
H=(1,\ -2v,\ 8v^2-2r,\ 24vr-48v^3).
$$

At infinite endpoints all CDF derivatives are zero.

### Fisher working curvature for cauchit

Unlike the other four densities, the Cauchy density is not log-concave.
The observed location curvature can be negative. Supply mgcv with positive
expected deviance curvature for its Fisher fallback, retaining the observed
likelihood derivatives above for the objective. With $u_j=\partial_\mu\log p_j$,

$$
I_D=2w\sum_j p_j u_j^2,
$$

$$
\partial_\mu I_D=2w\sum_j p_j
\{u_j^3+2u_j L_{\mu\mu,j}\},
$$

$$
\partial_{\theta_r} I_D=2w\sum_j p_j
\{L_{\theta_r,j}u_j^2+2u_j L_{\mu\theta_r,j}\}.
$$

These populate `EDmu2`, `EDmu3`, and `EDmu2th`. They require summation over all
categories, adding work to the cauchit fit. For the other links, the family
uses the nonnegative observed deviance curvature, as `ocat()` does.

## 6. Numerical evaluation

Direct subtraction of CDFs close to one is unreliable. For left/central
intervals compute

$$
\log p=\log F(b)+\log[-\operatorname{expm1}\{\log F(a)-\log F(b)\}],
$$

and for intervals with $F(a)>1/2$ use the survival equivalent

$$
\log p=\log S(a)+\log[-\operatorname{expm1}\{\log S(b)-\log S(a)\}],
\qquad S=1-F.
$$

R's log-CDF/log-survival functions are used for logistic, normal and Cauchy
errors. For cloglog, $\log S(z)=-e^z$ and
$\log F(z)=\log[-\operatorname{expm1}(-e^z)]$; for $z<-35$ the latter is
replaced by $z$, whose approximation error is below ordinary double precision
at that magnitude. Loglog follows by reflection.

For extreme-value links, derivatives expressed as differences of large powers
of $e^z$ suffer catastrophic cancellation even when the log probability is
accurate. If the opposite endpoint has log-probability ratio below -50, and
we are in the extreme tail (or the opposite endpoint is infinite), differentiate
the dominant endpoint log probability directly. These are $-e^a$ for cloglog
survival and $-e^{-b}$ for loglog CDF. Their successive derivatives are all
$-e^a$, or alternate signs for loglog. Threshold chain rules remain those in
Section 3. The finite-endpoint approximation discards terms below numerical
precision in this tail region; infinite-endpoint formulas are exact.

No probability floor is used to redefine the likelihood. This preserves tail
log likelihoods when probabilities themselves underflow. Finite precision
still limits extremely narrow intervals and extremely large predictors;
exponential overflow, coincident floating-point thresholds, or derivatives
outside the representable range cannot be repaired by this formula. Sparse
categories and separation can also cause estimation difficulties, as in other
ordinal models. Initialization uses empirical cumulative frequencies with a
half-observation per category and a minimum starting increment of 0.01.

## 7. Prediction, CDF, quantiles and simulation

`predict(fit, type="response")` returns the $n\times K$ matrix $p_{ij}$.
For fixed thresholds,

$$
\frac{\partial p_{ij}}{\partial\mu_i}=f(\alpha_{j-1}-\mu_i)-f(\alpha_j-\mu_i),
$$

so the returned delta-method standard error is the absolute derivative times
$\sqrt{X_i V_\beta X_i^\mathsf{T}}$. It **conditions on thresholds**; it does not
include their uncertainty or their covariance with coefficients. Discrete bam
prediction uses mgcv's compressed design-matrix operations.

For real $q$, the categorical CDF is zero below 1, one at or above $K$, and
$F(\alpha_{\lfloor q\rfloor}-\mu)$ otherwise. The quantile is

$$
Q(p\mid\mu)=\min\{j:P(Y\le j\mid\mu)\ge p\},\quad 0<p<1,
$$

with endpoint convention $Q(0)=1$, $Q(1)=K$. It is calculated by comparing
against cumulative probabilities, using strict comparisons so an exact boundary
belongs to the lower category. For log probabilities, comparisons stay on the
log scale; upper-tail input uses survival probabilities directly. Simulated
categories use $Q(U\mid\mu)$ for $U\sim\mathrm{Uniform}(0,1)$.

All three helpers accept `mu` as a latent location. Weights determine likelihood
contributions, not a different sampling distribution, and `scale` is fixed;
compatibility arguments `wt` and `scale` are ignored. Missing inputs propagate.
The family owns the threshold state through mgcv's `getTheta`/`putTheta`
closures, so these helpers use the thresholds fitted by gam or bam, including
after serialization. Response residuals are the observed category minus the
conditional median category. Deviance residuals use the sign of the latent
interval midpoint minus the location, as in `ocat()`.

## 8. Fitting methods and scope

| Engine | Supported smoothing selection |
| --- | --- |
| gam | REML, ML, NCV |
| bam, ordinary | REML, ML, fREML |
| bam, discrete | fREML, NCV (mgcv 1.9-4 or later) |

This is an `extended.family` with one predictor and global threshold parameters.
Gam optimizes thresholds within its outer criterion; bam updates them
conditionally within fitting iterations. Bam NCV optimizes the working weighted
linear model criterion, whereas gam NCV evaluates the ordinal model loss using
approximate leave-neighbourhood-out updates. Neither is exact repeated refitting.
Gam does not acquire bam's fast-REML algorithm by requesting `method="fREML"`.

A specific limitation observed with mgcv 1.9-4 is a native solver failure for
cauchit `gam(method="ML")` with **no penalized terms** and negative observation
curvature. Smooth-model ML and NCV are tested. For purely parametric cauchit ML,
use `bam(method="ML")` with tight iteration tolerances, or a dedicated optimizer.
The family does not replace the observed likelihood Hessian by a positive
surrogate to disguise this engine limitation.

Also, MASS 7.3-66's `polr.fit` truncates CDF arguments to $[-100,100]$, including
infinite category endpoints. Cauchy tail mass beyond those limits is appreciable,
so its cauchit fitting objective is not the exact likelihood used here. The
cauchit regression test instead uses direct, untruncated CDF differences in
an independent numerical optimizer. The other four links are compared to polr.

Case weights, offsets, fixed thresholds, and random effects in the location
predictor are compatible with this representation. There is no scale predictor,
threshold-specific predictor, or structured-threshold option. Those are separate
extensions. Per-observation threshold second-derivative storage is
$O\{n(K-2)^2\}$; many-level semiparametric transformation models need a different
computational design.

## 9. Verification and references

`tests/cumulative-link.R` checks central finite differences of the location, threshold
and mixed derivatives; comparison to mgcv's logistic derivatives; unpenalized
ML probabilities and thresholds against `MASS::polr()` for four links and an
independent exact-likelihood optimizer for cauchit; fitting routes; fixed
thresholds; distribution identities and boundaries; simulation; prediction
standard errors; binary reductions; weighting; and serialization. Numerical
agreement of small examples does not guarantee convergence for every dataset.

- [mgcv ordered categorical family](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/ocat.html).
- [mgcv family interfaces](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/family.mgcv.html).
- [mgcv bam fitting details](https://stat.ethz.ch/R-manual/R-devel/library/mgcv/html/bam.html).
- [MASS cumulative-link parameterization](https://stat.ethz.ch/R-manual/R-devel/library/MASS/html/polr.html).
- Wood, S. N., Pya, N. and Saefken, B. (2016). Smoothing parameter and model
  selection for general smooth models. *JASA* 111, 1548–1575.
  [doi:10.1080/01621459.2016.1180986](https://doi.org/10.1080/01621459.2016.1180986).
