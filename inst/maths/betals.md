# Beta mean and precision family

For independent observations strictly inside (0, 1), write

\[
Y_i\sim\operatorname{Beta}(a_i,b_i),\quad
 a_i=\mu_i\phi_i,\quad b_i=(1-\mu_i)\phi_i,
\qquad \operatorname{logit}(\mu_i)=\eta_{1i},\quad
 \log\phi_i=\eta_{2i}.
\]

The mean is mu and variance is mu (1 - mu) / (1 + phi). Case weights
multiply the entire log density, including its normalization:

\[
\ell_i=w_i\{(a_i-1)\log y_i+(b_i-1)\log(1-y_i)
 +\log\Gamma(\phi_i)-\log\Gamma(a_i)-\log\Gamma(b_i)\}.
\]

Zero-weight observations have zero likelihood and derivative contributions.
They are replaced with harmless interior responses and predictors before
arithmetic. Endpoints are not part of this likelihood and are rejected.

## Stable analytic derivatives

Let u = mu, v = 1 - mu, t = log(phi), A = logit(y) - eta1, and

\[
 B=u\log(y/u)+v\log\{(1-y)/v\},\qquad
 H(z)=\log\Gamma(e^z)-(e^z-1/2)z+e^z-\tfrac12\log(2\pi).
\]

Then the unweighted log density is exactly

\[
 e^t B+\tfrac12(t+\log u+\log v)-\log y-\log(1-y)
 -\tfrac12\log(2\pi)+H(t)-H(t+\log u)-H(t+\log v).
\]

The implementation evaluates the log density with `stats::dbeta()` and uses
this decomposition for derivatives. For p >= 1,

\[
 \partial_{\eta_1}^p B=u^{(p)}A-(p-1)u^{(p-1)},
\]

where the second term is zero for p = 1. Thus every mixed derivative of
`exp(t) * B` is `exp(t)` times the corresponding mean derivative of B,
regardless of how many t derivatives are taken. Logistic derivatives through
fourth order are explicit. Near y = u, with d = y - u, compute B by

\[
 B=-\sum_{k=2}^{10}\{u(-d/u)^k+v(d/v)^k\}/k,
\]

when abs(d) < 0.01 min(u, v), avoiding cancellation of linear terms.
The omitted terms are below double-precision accuracy relative to B.

For exp(z) < 16, derivatives of H use `lgamma(exp(z) + 1)` and polygamma
recurrences, avoiding singular polygamma terms at small shapes. At larger
shapes use the Stirling correction series through the inverse thirteenth
power; its derivatives are analytic powers of the same exponentials. This
avoids subtraction of large gamma derivatives in the precision direction.
Partial Bell polynomials compose H with log(u) + t and log(v) + t.

Derivatives of total order r are packed as eta1^r, eta1^(r-1) eta2, ...,
eta2^r, and contracted with the design matrices by `mgcv::gamlss.gH()`.
REML outer Newton requires up to fourth order; BFGS and NCV use the same
analytic evaluator. NCV passes row identifiers through the mgcv helper to
keep response and offset subsets aligned, including repeated deletion rows
and gamma adjustments. Sandwich information is the cross-product of weighted
observation scores.

## Conditional deviance

At fixed phi the saturated mean solves

\[
 \log\{y/(1-y)\}-\psi(\mu\phi)+\psi\{(1-\mu)\phi\}=0.
\]

The score decreases strictly in mu. Its unique solution lies between y and
1/2. The implementation solves for the logit of this mean. Deviance is twice
the weighted difference between saturated and fitted log densities, with
residual sign given by y - fitted mean. The null model maximizes likelihood
over a constant mean intercept with mean offsets and fitted precisions
retained. These measure conditional mean fit, holding precision fixed.

Unrestricted observation-specific precision would give an unbounded saturated
likelihood; it is deliberately held fixed in these diagnostics. Constant
precision specifies the same response distribution as `mgcv::betar()`, but
REML need not agree because of the different treatment of precision in the
restricted likelihood.

## Validation

`tests/betals-derivatives.R` checks all packed derivatives by differentiating
the preceding order, compares scores and Hessians with `mgcv::betar()`, tests
small and large shapes and near-boundary means, and verifies third- and
fourth-order coefficient contractions and sandwich information.
`tests/betals.R` covers separate smooths, REML Newton/BFGS, ordinary and grouped
NCV, offsets, NCV gamma adjustments and derivatives, replication weights,
direct unpenalized maximum likelihood, conditional deviance, response/link
predictions and standard errors, missing observations, serialization,
distribution callbacks, and invalid input.
