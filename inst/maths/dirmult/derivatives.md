# Dirichlet-multinomial family

Each independent response row contains C = K + 1 integer counts y_j. Its total
N is conditioned upon. Let eta_1,...,eta_K be baseline logits, eta_C = log(phi),
p = softmax(0, eta_1,...,eta_K), and alpha_j = phi p_j.

The complete log probability is

    lgamma(N+1) - sum(lgamma(y_j+1))
    + lgamma(phi) - lgamma(phi+N)
    + sum(lgamma(alpha_j+y_j) - lgamma(alpha_j)).

Prior weights multiply the log probability and each derivative once.
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

The term involving sA comes from the second derivative of log(p_j).
For higher orders let x_j = log(alpha_j), f_j(x_j) = F(exp(x_j), y_j)
and take predictor derivatives indexed by i_1, ..., i_r. The contribution
of category j is

    sum over partitions P of {1,...,r}:
      f_j^(|P|)(x_j) * product over blocks B in P of x_j,[i_b for b in B].

For order three the block sizes are 3, 2+1 (three partitions) and 1+1+1.
For order four they are 4, 3+1 (four), 2+2 (three), 2+1+1 (six) and
1+1+1+1. The same formula applies to repeated indices. First derivatives
of x_j with respect to a logit eta_k are I(j = k+1) - p_(k+1), and its
first concentration derivative is 1. For a block of order two or larger,
the derivative is zero if any index is the log-concentration predictor.
Otherwise it is the common derivative of the negative log-softmax
normalizer. Subtract f_0^(r)(log(phi)), with f_0(x) = F(exp(x), N), in the
all-concentration entry at each order.

Unique mixed derivatives are packed in mgcv's symmetric order. The
general-family coefficient transformation uses the third derivatives for
outer REML BFGS and neighbourhood cross-validation, and the third and
fourth derivatives for outer REML Newton. EFS uses the score and Hessian.

## Stable evaluation

Compute R(a,m) = F(a,m) - m log(a). The log probability becomes

    log multinomial coefficient + sum(y_j log(p_j))
    + sum(R(alpha_j,y_j)) - R(phi,N).

This cancels the N log(phi) terms symbolically. For m <= 64 use the exact
finite sums

    R = sum_{t=1}^{m-1} log1p(t/a)
    A = sum_{t=0}^{m-1} a/(a+t)
    B = -sum_{t=0}^{m-1} (a/(a+t))^2.

At small category shapes, separate the first log shape from the other
rising factors before summing the likelihood. For phi < 1 this gives
a leading term (the number of positive counts minus one) times log(phi),
instead of subtracting terms proportional to N log(phi).

For each t > 0 let q = a/(a+t) and v = q(1-q). The second through fourth
derivatives of F with respect to log(a) are the sums of v, v(1-2q) and
v(1-6v). Compute both q and its complement as logistic functions of signed
log ratios. Evaluate log1p(t/a) as softplus(log(t)-log(a)). The t=0 term is
exactly log(a) in F, so it contributes only to the first log derivative.
The formulas remain defined when exp(log(a)) underflows.

For larger counts at m/a < 0.05, expand log(1+t/a) through 16 terms. The
r-th log derivative multiplies series term k by (-k)^r. Normalized power
sums keep evaluation on the scale of m(m/a)^k. At small a, shift the gamma
ratio to Gamma(a+m)/Gamma(a+1) and account for log(a) analytically. Apply
the chain rule to the scaled polygamma differences: if T_k = a^k d^kF/da^k,
the second through fourth log derivatives are T1+T2, T1+3T2+T3, and
T1+7T2+6T3+T4. At intermediate and large a, a Stirling expansion in ratios
avoids subtracting large nearby gamma values. Its derivatives are
differentiated at the same order as the likelihood.

At N=1 all concentration derivatives vanish. As phi tends to infinity the
probability-predictor derivatives approach those of the multinomial, and
the concentration derivatives tend to zero. The first concentration
derivative uses differences of A-m near this limit to retain the small
correction to the multinomial likelihood.

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
