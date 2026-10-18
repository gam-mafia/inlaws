# Configurable inflated beta

`zoibetals(inflation = "both")` is one constructor with `"zero"` and
`"one"` configurations. It complements `ordbeta()`: membership and the
conditional beta response have separate predictors. No GAMLSS implementation
code was used. The calculations below follow directly from the normalized
mixture and the beta density.

## Distribution and predictor order

The first formula models $u=\operatorname{logit}(\mu)$, the second
$v=\log\phi$. For both endpoints, the third and fourth model
$z=\log(p_0/p_I)$ and $o=\log(p_1/p_I)$. Thus

$$ (p_0,p_I,p_1) = (e^z,1,e^o)/(1+e^z+e^o). $$

For a single endpoint, the third predictor is its binary logit relative to
the interior; the disabled endpoint probability is exactly zero. All formulas
can contain their own offsets and smooth terms, or be intercept-only. Shared
coefficients follow mgcv's list-formula syntax. Coefficients are not otherwise
coupled across predictors.

With $a=\mu\phi$, $b=(1-\mu)\phi$,

$$ \ell_i = w_i\begin{cases}
\log p_{0i},&y_i=0,\\
\log p_{1i},&y_i=1,\\
\log p_{Ii}-\log B(a_i,b_i)+(a_i-1)\log y_i+(b_i-1)\log(1-y_i),&0<y_i<1.
\end{cases} $$

Integrating the continuous term gives $p_I$; adding the atoms gives one.
Prior weights multiply this likelihood and do not change individual response
moments or draws. Zero-weight rows contribute exactly zero, including when
predictors would overflow. Responses outside the configured support are
rejected before fitting, including in zero-weight rows. Missing rows are
handled by the model frame, without clipping or changing endpoint values.

## Analytic gradient and Hessian

For the interior beta log density $L$, let
$T=\log y-\log(1-y)-\psi(a)+\psi(b)$. With $\psi$ and $\psi_1$
the digamma and trigamma functions,

$$\begin{aligned}
L_\mu &= \phi T,\\
L_\phi &= \psi(\phi)-\mu\psi(a)-(1-\mu)\psi(b)
 +\mu\log y+(1-\mu)\log(1-y),\\
L_{\mu\mu} &= -\phi^2\{\psi_1(a)+\psi_1(b)\},\\
L_{\mu\phi} &= T+\phi\{-\mu\psi_1(a)+(1-\mu)\psi_1(b)\},\\
L_{\phi\phi} &= \psi_1(\phi)-\mu^2\psi_1(a)-(1-\mu)^2\psi_1(b).
\end{aligned}$$

Writing $s=\mu(1-\mu)$, the predictor derivatives are

$$\begin{aligned}
L_u &= sL_\mu,& L_v&=\phi L_\phi,\\
L_{uu}&=s^2L_{\mu\mu}+s(1-2\mu)L_\mu, &
L_{uv}&=s\phi L_{\mu\phi},\\
L_{vv}&=\phi^2L_{\phi\phi}+\phi L_\phi.
\end{aligned}$$

These derivatives are zero for endpoint observations. For active membership
logits indexed by $j,k$, the score is $1\{c_i=j\}-p_j$ and the Hessian is
$-p_j(1\{j=k\}-p_k)$. Membership and interior parameter mixed derivatives
are zero. Weights multiply every derivative. The code packs symmetric
second derivatives using `mgcv::trind.generator()` and maps them to
coefficients with `mgcv::gamlss.gH()`. Shared coefficient contributions add.

REML smoothing selection uses EFS and requires only these first and second
derivatives. `available.derivs = 0` selects EFS automatically. Higher-order
outer Newton/BFGS smoothing optimization and NCV are not implemented. EFS in
mgcv 1.9-3/1.9-4 does not apply reduced smoothing-parameter maps, so partly
fixed smoothing parameters and shared smoothing IDs are rejected. All-free
unlinked and all-fixed smoothing parameters are supported. The preinitialize
hook provides the `offxset` alias used by those versions' EFS caller, ensuring
that fitting, initialization and prediction see the same offsets.

## Moments, predictions, and diagnostics

$$M=p_1+p_I\mu,\qquad
V=p_1+p_I\{\mu^2+\mu(1-\mu)/(1+\phi)\}-M^2.$$

Variance is evaluated using the nonnegative within/between mixture expression
$p_0M^2+p_I\{\mu(1-\mu)/(1+\phi)+(\mu-M)^2\}+p_1(1-M)^2$.

`fitted()` and native distribution callbacks use link-scale predictor matrices.
Response predictions have columns `mu`, `phi`, `p0`, `pI`, `p1`, `mean`,
`variance` for every configuration. `predict_zoibetals()` restores names
that mgcv may discard and selects complete means, variances, interior
parameters, probabilities or all columns. Prediction standard errors use
$J V_\beta J^T$, including cross-predictor covariance. For probabilities,
$\partial p_k/\partial\eta_j=p_k(1\{k=j\}-p_j)$; the mean and variance
follow by differentiating the formulas above. Disabled probabilities and their
standard errors are zero. Standard errors condition on selected smoothing
parameters; EFS supplies no smoothing uncertainty correction, even if
`unconditional=TRUE` is requested. They describe estimation uncertainty,
not predictive variation or uncertainty from selecting the model.

For $0\le q<1$, $F(q)=p_0+p_I I_q(a,b)$; $F(q)=0$ below zero and $F(q)=1$
at and above one. The left limits are $F(0-)=0$ and $F(1-)=p_0+p_I$.
Survival probabilities are evaluated directly as $p_1+p_I\{1-I_q(a,b)\}$,
using log-sum-exp for stable log probabilities. Quantiles are zero through
$p_0$, beta quantiles of $(u-p_0)/p_I$ in the interior, and one from
$p_0+p_I$ upwards, with support limits at probabilities zero and one.
Simulation first selects the component, then draws from the beta distribution
only for interior selections. Extreme beta shapes can produce rounded
endpoint values through finite precision; observations are never clipped.

The diagnostics factory binds these CDFs, left limits, quantiles, simulations,
and complete moments. Randomized PIT uses a uniform draw between $F(y-)$ and
$F(y)$; normal quantile residuals transform this PIT. Response residuals are
$y-M$ and Pearson residuals are $\sqrt{w}(y-M)/\sqrt V$.
An unrestricted saturated beta likelihood diverges as $\phi\to\infty$ with
$\mu=y$, so deviance and null deviance are `NA` and deviance residuals are
explicitly unavailable. Use PIT, quantile, response or Pearson diagnostics.

## Identification and validation scope

Every requested component must appear in positive-weight observations.
An active endpoint with no observations has a membership intercept tending to
minus infinity; it is rejected rather than disabled. A configuration explicitly
disabling that endpoint defines a different model. At least two distinct
interior responses are required; otherwise interior mean and precision are
unidentified or the precision likelihood is unbounded. These are necessary
screening conditions, not sufficient identification guarantees. Separation,
sparse components, a rank-deficient design, or a mean predictor interpolating
interior observations can still prevent finite well-estimated parameters.
There are no artificial parameter bounds; nearly constant interior responses
can lead to very large precision.

The likelihood factorizes into membership and interior-beta likelihoods when
coefficients are separate. Unpenalized joint estimates can therefore be
checked against independently fitted component likelihoods. Separate REML
smoothing selections are not required to match joint REML selection: the
joint marginal criterion and effective sample-size terms differ. The regression
script `tests/zoibetals.R` checks finite differences of gradients,
Hessians and prediction transformations; normalization, integrated moments,
endpoint CDF jumps and quantiles; simulated moments and uniform randomized PIT;
all three fits with varying precision; independent unpenalized component fits;
weights, offsets, missing rows, absent and sparse components; full covariance
prediction errors; EFS selection and fixed penalties; shared coefficients; and
the native diagnostics contract. `inst/examples/zoibetals.R` demonstrates
independent membership and interior covariates alongside `ordbeta()`.
