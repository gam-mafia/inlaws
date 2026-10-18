# inlaws <img src="man/figures/logo.png" align="right" height="160" alt="inlaws hex sticker with a family of three distribution characters" />

[![Lifecycle: experimental](man/figures/lifecycle-experimental.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R CMD check](https://github.com/gam-mafia/inlaws/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gam-mafia/inlaws/actions/workflows/R-CMD-check.yaml)
[![Test coverage](https://github.com/gam-mafia/inlaws/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/gam-mafia/inlaws/actions/workflows/test-coverage.yaml)
[![pkgdown](https://github.com/gam-mafia/inlaws/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/gam-mafia/inlaws/actions/workflows/pkgdown.yaml)

Documentation: <https://gam-mafia.github.io/inlaws/>.

Additional distribution families for **mgcv**.

**inlaws is highly experimental and largely untested.** All families are
currently marked as experimental. Interfaces, defaults, and behaviour may
change without notice, and fitted results may be incorrect. The existing
tests and checks provide only partial validation; they do not establish
reliability across datasets, model specifications, or fitting methods.
Independently validate results before relying on them.

Consult each family's help for parameterizations, limitations, and supported
fitting methods. Experimental markers will be removed from individual
families as their behaviour is confirmed and tested.

See the [family compatibility vignette](https://gam-mafia.github.io/inlaws/articles/family-compatibility.html)
for the current `gam()`/`bam()` and smoothing-selection support matrix, or run
`vignette("family-compatibility", package = "inlaws")` after installing with
vignettes. Building the vignette requires the Quarto CLI and the R package
**quarto**.

| Response | Family constructors |
| --- | --- |
| Beta mean and precision | `betals()` |
| Beta-binomial | `betabinomial()` |
| Censored gamma | `cgamma()`, `cgammals()` |
| Log-logistic mean and scale | `log_logistic()`, `log_logisticls()` |
| Generalized gamma mean, scale and shape | `gen_gammalss()` |
| Censored generalized gamma | `cgen_gammalss()` |
| Censored log-normal | `clognormal()`, `clognormalls()` |
| Conway-Maxwell-Poisson | `cmp()`, `cmpls()` |
| Cumulative-link ordinal | `cumulative_link()` |
| Dirichlet | `dirichlet()` |
| Dirichlet-multinomial | `dirmult()` |
| Generalized Poisson | `gpoisson()` |
| Poisson inverse Gaussian | `pig()`, `pigls()` |
| Negative binomial location and size | `nbls()` |
| Ordered beta | `ordbeta()` |
| Configurable zero/one inflated beta | `zoibetals()` |
| Zero-inflated / hurdle negative binomial | `zinb()`, `zanb()` |
| Generalized Pareto | `gpd()` |

```r
library("inlaws")

set.seed(1)
d <- data.frame(x = runif(300))
d$y <- rnbinom(300, mu = exp(1 + d$x), size = 3)
fit <- mgcv::gam(list(y ~ s(x, k = 6), ~ 1),
                 data = d, family = nbls(), method = "REML")
predict(fit, type = "response")
fit$family$rd(fitted(fit))
fit$family$qf(0.95, fitted(fit))
```

`zoibetals(inflation = "both")` complements `ordbeta()` with separate
predictors for the interior beta mean, precision, and endpoint membership.
Use `"zero"` or `"one"` for one active endpoint. Formula order is mean,
precision, then zero and/or one logits relative to the interior; `~ 1` is
allowed for every parameter. `predict_zoibetals()` returns named
component parameters, probabilities, or complete response moments. It supports
REML with EFS and randomized quantile diagnostics. See `?zoibetals` and
`inst/examples/zoibetals.R` for a reproducible comparison with `ordbeta()`.

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

GitHub Actions runs `R CMD check` on Linux, macOS, and Windows. The coverage
workflow runs the base-R tests with **covr**, reports the percentage in the run
summary, and uploads XML, CSV, and RDS coverage reports as the `coverage`
artifact. No external coverage service or token is required.

The secondary `brms-validation` workflow compares the Dirichlet likelihood,
conditional means, and generated Stan parameterization with **brms** on Linux.
It runs for changes to R code, package metadata, validation scripts, or its
workflow, and can also be started manually. It runs independently of the
ordinary checks, coverage, and site build, which do not require **brms**.
To run the comparison locally, install **brms** explicitly with
`install.packages("brms")`, install the current **inlaws** package, and run
`source("inst/validation/dirichlet-brms.R")`. No Stan compilation or sampling
is performed by this comparison.

Build the documentation locally with `pkgdown::build_site()` (requires
**pkgdown** and Quarto). The pkgdown workflow builds pull requests and uploads
a `pkgdown-preview` artifact; pushes to `main` build and deploy the site to
GitHub Pages. The repository's Pages publishing source must be **GitHub
Actions**. All workflows can also be run manually from the Actions tab.

## Fitted-model diagnostics

Every family supplies a `gratia_diagnostics(mu, wt, scale, y, eta)` function
for gratia's extensible diagnostic interface. It describes the complete response
distribution at fitted parameters, including response moments and count or
endpoint masses. Multivariate families retain response-column dependence in
joint draws. Constructing and fitting inlaws families does not require gratia.
See `?family_diagnostics` for moment conditions and componentwise checks.

Censored families use observed intervals for PIT and quantile residuals.
Traditional simulated QQ or worm references need the known observation rule:

```r
library("inlaws")
f <- cgamma(theta = 0.6, observe = censoring_process(lower = 0.4, upper = 4))
```

The limits in `censoring_process()` are scalars shared across observations, so
they remain aligned after subsetting and missing-row removal. Unknown censoring
rules cannot be recovered from realised exact observations. Public response
samplers continue to draw the uncensored response distribution.
