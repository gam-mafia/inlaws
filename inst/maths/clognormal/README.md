# Arithmetic-mean censored lognormal family

Let m = log(mu) - sigma^2/2, t = log(sigma), and z(v) = (log(v)-m)/sigma.
The predictor is eta = log(mu). Thus mu is E[Y], not exp(E[log(Y)]).
For non-informative censoring the individual log likelihood contributions are:

- Exact y: -log(y) - t - log(2*pi)/2 - z(y)^2/2.
- Left censoring at U: log Phi(z(U)).
- Right censoring at L: log Phi(-z(L)).
- Interval (L,U): log[Phi(z(U))-Phi(z(L))].

Case weights multiply these contributions. They do not modify sigma. The
Jacobian -log(y) belongs only to exact observations: interval probabilities
are invariant under the transformation.

For mgcv's extended-family interface, use D_i = 2(s_i - ell_i). Here s_i is the
supremum over the location m with sigma fixed: -log(y)-t-log(2*pi)/2 for exact
observations, zero for one-sided censoring, and log[2 Phi(d/(2 sigma))-1] for
an interval of log width d. The ls component returns sum(w_i s_i) and its
first two t derivatives, so sum(D_i)/2 - ls is exactly the negative log
likelihood on the original response scale.

The implementation propagates truncated bivariate Taylor polynomials in mu
and t. A coefficient of mu^a t^b is the corresponding derivative divided by
a! b!. Polynomial multiplication is convolution; exp and log are composed
using their fourth-order Taylor series. Only coefficients with a+b <= 4 and
b <= 2 are needed. This supplies all Dd derivatives analytically, including
mixed derivatives induced by m = log(mu)-exp(2t)/2, without differencing or
calling another family.

For censoring probabilities, use log CDF differences (or survival differences
in the upper tail). The first four CDF derivatives are phi(z) times
1, -z, z^2-1, and 3z-z^3. Scale those derivatives by the endpoint probability
before taking the Taylor logarithm, then take a stable log difference.
Below a standardized endpoint of -12, use the 16-term Mills asymptotic
series for log Phi to avoid cancellation in high-order tail derivatives. For small standardized intervals, use
12-point Gauss-Legendre integration over the log-response interval, factoring
out the midpoint log density. This avoids cancellation in CDF derivatives.
Compute the interval log width with log1p((U-L)/L), retaining widths smaller
than the floating-point spacing of log(L). Wide-ratio overflow falls back to
log(U)-log(L).

EDmu2 supplies the complete-data expected curvature 2w/(mu^2 sigma^2) as a
positive fallback if observed curvature is unsuitable. The likelihood and
observed derivatives remain those of the censored data; this is not a claim
that censoring preserves complete-data information.

The tests check each derivative slot against centered differences, compare
probabilities with R's lognormal functions, check likelihood optimization,
and compare fits against independently transformed Gaussian/censored-normal
fits. These reference families are used only in tests. Standard predict.gam
inference conditions on estimated sigma. No integration over coefficient
uncertainty or random effects is added to response means.
