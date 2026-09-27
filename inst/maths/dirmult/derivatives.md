# Dirichlet-multinomial family

Each independent response row contains C = K + 1 integer counts y_j. Its total
N is conditioned upon. Let eta_1,...,eta_K be baseline logits, eta_C = log(phi),
p = softmax(0, eta_1,...,eta_K), and alpha_j = phi p_j.

The complete log probability is

    lgamma(N+1) - sum(lgamma(y_j+1))
    + lgamma(phi) - lgamma(phi+N)
    + sum(lgamma(alpha_j+y_j) - lgamma(alpha_j)).

Prior weights multiply this expression and both derivative orders once.
There is no independent-observation approximation across categories.

## Derivatives

For F(a,m) = log Gamma(a+m) - log Gamma(a), define

    A(a,m) = a F'(a,m)
    B(a,m) = a^2 F''(a,m).

Then dF/dlog(a) = A and d^2F/dlog(a)^2 = A+B. Write A_j, B_j for
(a,m) = (alpha_j,y_j), and A_0,B_0 for (phi,N). Let D_j = A_j+B_j,
sA = sum(A_j), sD = sum(D_j), and u = p[-1].

The score on the linear predictor scale is

    score_k = A_(k+1) - u_k sA       (k = 1,...,K)
    score_C = sA - A_0.

The symmetric observed Hessian is

    H_kl = u_k u_l (sD+sA) - u_k D_(l+1) - u_l D_(k+1)
           + I(k=l) (D_(k+1)-u_k sA)
    H_kC = D_(k+1) - u_k sD
    H_CC = sD - A_0 - B_0.

The term involving sA is the second derivative of log(p_j). Omitting it is an
incorrect chain rule and yields the wrong coefficient Hessian. `gamlss.gH()`
then assembles every coefficient block, including shared predictor coefficients.
Only the packed second-order index matrix is allocated; there are no third- or
fourth-order tensors. EFS needs only this score and Hessian.

## Stable evaluation

Compute R(a,m) = F(a,m) - m log(a). The log probability becomes

    log multinomial coefficient + sum(y_j log(p_j))
    + sum(R(alpha_j,y_j)) - R(phi,N).

This cancels the N log(phi) terms symbolically. For m <= 64 use the exact
finite sums

    R = sum_{t=1}^{m-1} log1p(t/a)
    A = sum_{t=0}^{m-1} a/(a+t)
    B = -sum_{t=0}^{m-1} (a/(a+t))^2.

Evaluate log1p(t/a) as softplus(log(t)-log(a)) and the ratios as logistic
functions. The t=0 contribution is handled separately. Thus underflow of
alpha does not invalidate zero counts or positive counts in rare categories.

For larger m with a > 1000 m, use four terms of the series in t/a. With
S_r = sum_{t=1}^{m-1} (t/a)^r, computed from power-sum formulas,

    R = S_1 - S_2/2 + S_3/3 - S_4/4
    A = m - S_1 + S_2 - S_3 + S_4
    B = -m + 2 S_1 - 3 S_2 + 4 S_3 - 5 S_4.

The neglected R remainder is bounded by m (m/a)^5 / 5. Otherwise use gamma
and polygamma differences; for a < 1 first apply Gamma(a+1) = a Gamma(a)
to remove the singular contribution analytically. At N=1 all concentration
derivatives vanish. As phi tends to infinity the multinomial log probability
and its probability-predictor derivatives remain.

## Fitting and interpretation

The last formula can be `~1` or contain smooth/parametric effects. Initialization
uses pseudocounts only for starting logits, a one-dimensional search for starting
concentration, and joint penalized least squares to respect shared coefficients.
The concentration search interval constrains initialization only, never fitting.

Residual deviance is -2 times log probability (a probability-one reference),
not a saturated Dirichlet-multinomial likelihood ratio. Null deviance fits a
weighted intercept-only composition/concentration model with the same offsets.
Pearson residuals are marginal, correlated category residuals. Response-scale
standard errors propagate the full coefficient covariance.

The package supplies a family to unmodified mgcv. It does not enable general
families in `bam()`. Future scalable fitting needs chunked accumulation of
likelihood, score and all cross-predictor Hessian blocks, then discrete-design
support, covariance/prediction integration, and comparisons against `gam()`.
Flattening counts into independent category observations is not equivalent.
