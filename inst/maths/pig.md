# Poisson inverse Gaussian: likelihood and derivatives

## Parameterization

Let `Z ~ IG(mean = 1, shape = 1/theta)` and `Y | Z ~ Poisson(mu * Z)`.
Then `E(Y) = mu` and `Var(Y) = mu + theta * mu^2`. The mean has a log link;
`t = log(theta)` is one global nuisance parameter. This is the PIG
parameterization of GAMLSS, with its `sigma` equal to our `theta` (not squared).
See the [GAMLSS PIG documentation](https://search.r-project.org/CRAN/refmans/gamlss.dist/html/PIG.html).

The mixing density is

```
f(z) = (2*pi*theta*z^3)^(-1/2) * exp(-(z-1)^2/(2*theta*z)).
```

Integrating the Poisson probability gives, with `s = sqrt(1 + 2*mu*theta)`,

```
p(y) = sqrt(2/(pi*theta)) * mu^y * s^(1/2-y) / y! *
       exp(1/theta) * K_(y-1/2)(s/theta).
```

## Stable recurrence

Use `K_(1/2)(x) = sqrt(pi/(2*x))*exp(-x)` and the Bessel recurrence.
Set `a = theta/s`, `h_0 = 1`, and for k >= 1 let

```
h_k = (2*k-1)*a + 1/h_(k-1).
log p(0) = -2*mu/(1+s).
log p(y) = log p(0) + y*(log(mu)-log(s)) - lgamma(y+1)
           + sum(log(h_k), k=1,...,y-1).
```

An empty sum is zero. The rationalized zero probability avoids subtracting
`1-s` near theta=0. Each recurrence step uses positive terms; no unscaled
Bessel function, factorial, or probability product is formed. The recurrence
converges continuously to Poisson as theta decreases, without model switching.
Its cost grows with the largest count and the number of observations still
active at each step. No recurrence state is shared across family instances.

## Analytic derivatives

The numerical core propagates bivariate Taylor polynomials through order four
in `(mu, t)`. A coefficient at `(i,j)` equals the derivative divided by `i!*j!`.
Multiplication is truncated convolution. Logarithms, reciprocals and square
roots use their series in the relative increment around the positive constant
coefficient. The theta polynomial is `exp(t)*sum(dt^j/j!,j=0,...,4)`.
Thus the likelihood and every derivative use the same finite recurrence;
there is no numerical differentiation in fitting.

`Dd` returns derivatives of `D = 2*w*(l_sat(t)-l(y,mu,t))`. This includes
mean derivatives through order four and mixed derivatives with up to two
log-dispersion derivatives as required by mgcv.

For y>0, bracket and solve `l_mu(y, mu_sat, t)=0` in log mean. For y=0 the
saturated mean is zero and `l_sat=0`. The envelope derivatives are

```
d l_sat/dt = l_t
 d^2 l_sat/dt^2 = l_tt - l_mut^2/l_mumu,
```

where the right-hand sides are evaluated at the saturated mean. `ls` returns
these weighted saturated likelihood terms. Consequently `ls - D/2` is the
actual weighted log likelihood, including its dispersion derivatives. A
single-entry cache of saturated results belongs to each family environment.

## Working information and fitting

Observed mean curvature is negative at zero counts. It cannot serve as
mgcv's positive expected-information fallback, including in `initial.spg`.
The implementation therefore provides variance-based working information:

```
V = mu + theta*mu^2
EDmu2 = 2*w/V
EDmu3 = -2*w*(1+2*theta*mu)/V^2
EDmu2th = -2*w*theta*mu^2/V^2.
```

These are a positive working approximation, not a claim that PIG is a
one-parameter exponential family in mu at fixed theta. All observed
likelihood/deviance derivatives remain exact. REML/ML optimization and NCV
use mgcv's extended-family machinery. Discrete bam NCV uses its working-model
criterion, so its optimum is not expected to equal gam NCV.

Before each fit, serialization clones family state while preserving fields
injected by mgcv (notably `qapprox`). A `fix.family.link.pig.family` S3 method
clones before bam saves its local and setup-family references, then delegates
to mgcv's extended-family link method. The preinitialization callback clones
only when this has not already happened. This matters because discrete bam
finalizes AIC using its original family reference; a later clone would leave
that reference with the starting dispersion. No mgcv namespace or calling
frame is modified. Postprocessing recomputes deviance at the final parameters.

Like nb(), pig() fixes the additional scale at one inside bam's fitter;
no explicit scale argument is needed for fitting. In mgcv 1.9-4 the outer
bam routine nevertheless marks scale as estimated under the default argument.
This makes logLik and its attached degrees of freedom each one too high.
Explicit scale=1 avoids that reporting inconsistency without changing
coefficients, smoothing parameters, dispersion, covariance, or AIC (the two
reporting errors cancel in AIC). Tests compare default and explicit-scale
calls for every supported bam route, with fixed and estimated dispersion.

Zero-weight rows are excluded before
numerical evaluation. Null deviance optimizes a constant log-mean predictor
with offsets and the fitted dispersion retained, or fixes that predictor at
zero if the model has no intercept.

## Distribution callbacks and boundaries

Simulation uses the Michael--Schucany--Haas inverse Gaussian construction:
for `v=theta*N(0,1)^2/2`, take the small root
`z=1/(1+v+sqrt(v)*sqrt(2+v))`, then replace it by its reciprocal with
probability `z/(1+z)`. Finally draw Poisson with mean `mu*z`.

CDF and quantile callbacks accumulate the actual recurrence probabilities in
log space. They never renormalize a truncated distribution. The hard limit
is one million terms; exhaustion is an explicit error. For positive means,
`qf(0)=0`, `qf(1)=Inf`; zero means give a point mass at zero. CDFs floor their
argument. Missing arguments propagate, and distribution arguments use scalar
or equal-length recycling. Prior weights do not change individual draws.

The all-zero positive-weight response has no finite log-mean estimate and is
rejected. Very large counts can be expensive, extreme tails have floating-point
limits, and near-Poisson dispersion may be weakly identified. There is no
fixed theta=0 constructor: theta=0 retains the package convention of requesting
estimation from the default starting value.

## Validation

The base-R tests check the likelihood by independent mixing-density integration
and scaled Bessel evaluation; normalization and moments; each analytic
Taylor derivative against a numerical derivative of the preceding order;
saturated envelope derivatives; weighted callback identities; and all
advertised fitting routes with fixed and estimated dispersion. They also
cover default/grouped NCV, integer-weight replication, offsets, template
isolation, conditional simulation, CDF/quantile inversion and Poisson limits.
