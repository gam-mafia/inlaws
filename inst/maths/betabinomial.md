# Beta-binomial mean and precision family

For independent groups, `Y | P ~ Binomial(m, P)` and
`P ~ Beta(mu phi, (1-mu) phi)`. The two linear predictors are
`eta1 = logit(mu)` and `eta2 = log(phi)`. The count variance is
`m mu (1-mu) {1 + (m-1)/(1+phi)}`. The trials within a group have
intraclass correlation `1/(1+phi)`; this does not model dependence between
separate response rows.

`m=0` contributes no likelihood or derivatives. For `m=1` the likelihood is
Bernoulli and all derivatives involving precision are exactly zero. Precision
is unrestricted within its positive domain, so binomial/underdispersed data
can push its estimate to infinity. All-success/all-failure datasets and datasets
without a positive-weight group of at least two trials are rejected. Other
forms of separation and local nonidentifiability remain possible.

## Likelihood and analytic predictor derivatives

Write `a=mu phi`, `b=(1-mu) phi`, `f=m-y`. The log probability is

```
l = lchoose(m,y) + lbeta(y+a,f+b) - lbeta(a,b).
```

Direct beta-function subtraction is inaccurate near the binomial limit.
Instead define

```
F(t,n) = sum(k=1,...,n-1) log(1 + k exp(-t)),  F(t,0)=F(t,1)=0.
l = lchoose(m,y) + y log(mu) + f log(1-mu)
    + F(log(mu)+eta2,y) + F(log(1-mu)+eta2,f) - F(eta2,m).
```

For each summand set `r=k/(exp(t)+k)` and `s=1-r`, evaluated using the two
opposite logistic functions to avoid cancellation. Its first four derivatives
with respect to `t` are

```
-r, r*s, r*s*(r-s), r*s*(1-6*r*s).
```

`F` and its derivatives are evaluated together:

* For `n<=64`, use the finite recurrence, vectorized in batches of at most 256
  rows. Softplus evaluates the log term without overflowing.
* For larger `n` and `n/exp(t)<0.05`, use 16 terms of
  `F = sum(r>=1) (-1)^(r+1) exp(-r*t) sum(k^r)/r`. Each `t` derivative
  multiplies term `r` by `-r`. Normalized Faulhaber recurrences calculate the
  power sums without forming large powers of counts or precision. This retains
  small, nonzero precision derivatives near the binomial limit.
* Otherwise use gamma/polygamma differences shifted to `exp(t)+1`, omitting
  the singular zero-index term analytically. Stirling numbers convert the
  natural-scale derivatives to `t` derivatives. For shape at least 8, evaluate
  the log-gamma difference with `log1p(n/a)` and the Stirling correction through
  inverse power 13, avoiding subtraction of large log-gamma values.

For `u=log(mu)`, its derivatives with respect to `eta1` are

```
u1 = 1-mu
u2 = -mu*(1-mu)
u3 = -mu*(1-mu)*(1-2*mu)
u4 = -mu*(1-mu)*(1-6*mu*(1-mu)).
```

For `v=log(1-mu)`, only the first derivative differs: `v1=-mu`.
For `p>=1`, the mixed derivative of `F(u+eta2,n)` is

```
D_eta1^p D_eta2^q F = sum(k=1,...,p) B[p,k](u1,...) F^(q+k),
```

using partial Bell polynomials. For `p=0` it is simply `F^(q)`.
Add the corresponding `v` contribution, subtract the total-trials term when
`p=0`, and add `y*u[p]+f*v[p]` when `q=0`. Store symmetric derivatives in mgcv's
packed order, from all-mean to all-precision derivatives. This performs the
link transformation analytically, directly on the predictor scale; a second
call to `gamlss.etamu` is neither needed nor correct.

Multiply every log likelihood and derivative by its likelihood weight exactly
once. Pass derivatives to `gamlss.gH` for coefficient and smoothing-parameter
contractions. For fourth-order contractions with shared coefficients, sum
predictor-pair contributions explicitly: mgcv 1.9-3/1.9-4's expanded-design
trace code overwrites overlapping blocks. This correction is local to this
family and leaves the fitter and other families unchanged.

## Fitting and response representation

`preinitialize` runs after model-frame subsetting and missing-value removal.
For a two-column response it extracts counts and totals and retains prior
weights as likelihood multipliers. For a proportion response it reconstructs
integer counts using weights as totals and replaces likelihood weights with
ones. Integer tolerance is `min(1e-7,64*eps*max(1,abs(value)))`; totals above
`2^53` are rejected. No response rounding beyond this tolerance is allowed.

The fitted object stores success counts in `y`, totals in `family$trials`, and
likelihood weights in `prior.weights`. Each fit receives its own family copy.
There is no global or shared mutable training-data state. Means and precisions
are predicted without needing future trial totals.

Starting coefficients jointly regress adjusted logits and constant `log(10)`
precision on both predictor blocks, subtracting offsets and adding a weak
penalty for numerical regularization. This handles shared coefficients.
User-supplied starts take precedence. These starting adjustments do not alter
the likelihood or bound fitted parameters.

REML uses mgcv's Laplace criterion, with analytic fourth derivatives for outer
Newton or third derivatives for outer BFGS. There is no extra estimated scale.
True ML smoothing selection is not provided by the general-family fitter;
mgcv can silently convert an ML request to REML. `bam()` and `gamm()` are outside
this family's interface.

## NCV

The family supplies observation-level predictor derivatives through order
three, direct `eta` evaluation, and a wrapper around mgcv's internal
`gamlss.ncv` helper. The wrapper passes scalar row identifiers as the helper's
response. A local likelihood callback maps those identifiers back to counts
and totals, restores predictor indices on subset model matrices, and evaluates
the original likelihood with the supplied `eta`. This preserves repeated or
reordered prediction indices and varying trial totals. Supplied `eta` already
includes offsets, so the callback must not add offsets a second time. For
the full-data score used with `gamma != 1`, the helper instead supplies
`X beta` without offsets; reconstruct this predictor through the ordinary
likelihood path so that the score and its derivative agree.

The family supports NCV only. QNCV is rejected: the current mgcv helper's
QNCV derivative path does not correctly align general prediction indices
with full-data predictor derivatives. NCV's one-step deletion approximation
is not an exact refit. Leave out entire rows
(groups), or neighbourhoods of rows to address between-group dependence.
Weak precision information can make deletion approximations unreliable.

## Diagnostics

Count response residuals are `y-m*mu`; Pearson residuals divide by the count
standard deviation and multiply by the square root of the likelihood weight.
Inactive Pearson and deviance residuals are zero.

Conditional deviance holds precision fixed and maximizes each row likelihood
over its mean. Endpoint counts have saturated log probability zero. For an
interior count, the maximizing mean lies between `y/m` and `1/2`; use bounded
optimization on the corresponding logit interval. It is generally not `y/m`.
Null deviance fits a constant mean predictor with first-predictor offsets and
fitted precisions retained. These diagnostics do not compare both predictors
to a joint null model.

Warn if any informative row has a precision-predictor standard error greater
than 5, or the smallest eigenvalue of the unpenalized precision information
block is at most `1e-10 * max(1,max(abs(eigenvalues)))`. This is a heuristic
warning about weak or boundary identification, not a proof of regularity when
no warning occurs. Smoothing cannot identify an uninformed precision intercept.

## Validation

The three beta-binomial test scripts cover normalization and moments,
Bernoulli/binomial limits, all mixed derivatives through fourth order,
coefficient contractions including shared coefficients, and derivatives of the
REML/NCV criteria away from their optima. Integration checks cover both response
formats, row replication, offsets, missing data, zero-information rows,
independent unpenalized likelihood optimization, simulated predictor recovery,
Newton/BFGS agreement, and approximate NCV versus explicit leave-out refits.
