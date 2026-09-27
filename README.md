# inlaws

Additional distribution families for **mgcv**. This is a development package;
consult each family's help for parameterizations, limitations, and supported
fitting methods. Availability does not imply identical validation maturity.

| Response | Family constructors |
| --- | --- |
| Beta-binomial | `betabinomial()` |
| Censored gamma | `cgamma()` |
| Censored log-normal | `clognormal()`, `clognormalls()` |
| Conway-Maxwell-Poisson | `cmp()` |
| Cumulative-link ordinal | `ocat_link()` |
| Dirichlet | `dirichlet()` |
| Dirichlet-multinomial | `dirmult()` |
| Generalized Poisson | `gp()` |
| Negative binomial location and size | `nbls()` |
| Ordered beta | `ordbeta()` |
| Zero-inflated / hurdle negative binomial | `zinb()`, `zanb()` |
| Generalized Pareto | `gpd()` |

```r
library(inlaws)

set.seed(1)
d <- data.frame(x = runif(300))
d$y <- rnbinom(300, mu = exp(1 + d$x), size = 3)
fit <- mgcv::gam(list(y ~ s(x, k = 6), ~ 1),
                 data = d, family = nbls(), method = "REML")
predict(fit, type = "response")
fit$family$rd(fitted(fit))
fit$family$qf(0.95, fitted(fit))
```

The package contains families and their distribution callbacks, without general
simulation or quantile wrappers for fitted models. Callback parameter scales and
matrix-response support differ by family; see its documentation. Mathematical
notes are under `inst/maths`, and optional external comparisons under
`inst/validation`. The existing mgcvUtils `clognorm()` and smooth constructors
are outside this package's scope.

## Development

From the repository's parent directory:

```sh
R CMD build inlaws
R CMD check --no-manual inlaws_0.0.0.9000.tar.gz
```

The base-R regression tests live in `tests/` and run during `R CMD check`.
Regenerate documentation with `roxygen2::roxygenise("inlaws")` after editing
roxygen comments. See `MIGRATION.md` for provenance and original attribution.
