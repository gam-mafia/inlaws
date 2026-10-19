#' Censored generalized gamma mean, scale and shape family
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Model the conditional arithmetic mean, scale and transformed shape of a
#' generalized gamma response with non-informative censoring.
#' @param link Three links, all `"log"`, for mean, scale and transformed shape.
#' @param observe Optional function applying a known observation rule to
#'   simulated responses for diagnostic references. See [censoring_process()].
#' @return An mgcv `general.family` with distribution and diagnostic callbacks.
#' @details Supply three formulae to [mgcv::gam()], with the response in the
#' first and one-sided scale and shape formulae second and third. As in
#' [gen_gammalss()], the predictors are `log(mu)`, `log(sigma)`, and `log(r)`,
#' where `mu = E(Y | x)` is the arithmetic mean of the underlying response,
#' `sigma > 0`, and positive `r`, equal to `1 + sigma * Q`. Recover shape as
#' `Q = (r - 1) / sigma`. Use `~ 1` in both nuisance formulae for constant
#' sigma and Q. The mean always exists; variance exists only for \eqn{r>0.5}.
#' The log-normal limit has `r` equal to 1; the gamma submodel has `r` equal
#' to `1 + sigma^2`.
#'
#' Supply a positive numeric vector for exact responses, or a two-column
#' matrix as in [cgamma()]: `(y, y)` for an exact observation, `(limit, -Inf)`
#' for left censoring, `(limit, Inf)` for right censoring, and `(lower, upper)`
#' for finite intervals. All types can be mixed. Zero lower bounds and the
#' uninformative interval `(0, Inf)` are supported. Missing rows follow
#' `gam()`'s `na.action`. Finite nonnegative case weights multiply likelihood
#' contributions; at least one weight must be positive.
#'
#' Supports REML (outer Newton or BFGS) and NCV/QNCV (BFGS), including
#' neighbourhoods supplied through `nei` and offsets in all predictors.
#' Separate coefficient blocks are required. Shared coefficients, ML and
#' `gamm()` are outside the supported scope. mgcv currently rejects general
#' families in `bam()`, including discrete fitting. Censoring can weaken
#' identification of scale and shape; sufficient information about the
#' underlying distribution is needed to estimate all three predictors.
#'
#' Censored contributions are log probabilities, while exact observations
#' contribute log densities including their constants and response Jacobian.
#' Censored derivatives through order four integrate analytic density
#' derivatives using adaptive quadrature, preserving shape information at
#' `Q = 0`. This costs more than the uncensored family. Extreme predictors or
#' poorly identified models can exceed numerical accuracy or yield singular
#' fits. No additional fitted shape constraint is imposed.
#'
#' Fitted values and response predictions have columns `mu`, `sigma`, `r`;
#' link predictions give their logarithms. Standard errors use the delta
#' method. The `$cdf`, `$qf`, and `$rd` callbacks have the same arguments and
#' underlying, uncensored response distribution as [gen_gammalss()]. They
#' ignore weights and external scale. The `observe` argument supplies a
#' censoring process for diagnostic reference simulations, without changing
#' these distribution callbacks.
#'
#' Deviance maximizes each observation's likelihood over mean at fitted sigma
#' and r. Null deviance retains those nuisance parameters and mean offsets,
#' estimating a constant mean predictor when an intercept is present.
#' Response and Pearson residuals are `NA` for censored observations; Pearson
#' residuals are also `NA` when variance is infinite. Probability diagnostics
#' use the censoring bounds. Mathematical details are in
#' `inst/maths/cgen-gammalss.md`.
#' @seealso [gen_gammalss()], [cgamma()], [censoring_process()]
#' @template ncv-neighbourhoods
#' @export
#' @examples
#' set.seed(31)
#' dat <- data.frame(x = runif(200))
#' pars <- cbind(exp(0.5 + sin(6 * dat$x)), 0.5, 1.2)
#' y <- gen_gammalss()$rd(pars)
#' dat$y <- cbind(pmax(y, 0.5), ifelse(y < 0.5, -Inf, y))
#' fit <- mgcv::gam(list(y ~ s(x, k = 6), ~ 1, ~ 1), data = dat,
#'   family = cgen_gammalss(observe = censoring_process(lower = 0.5)),
#'   method = "REML")
#' pars_fit <- fitted(fit)
#' head(cbind(pars_fit, Q = (pars_fit[, "r"] - 1) / pars_fit[, "sigma"]))
cgen_gammalss <- function(link = list("log", "log", "log"), observe = NULL) {
  family <- gen_gammalss(link)
  family$family <- "cgen_gammalss"
  family$prepare.response <- .cgen_gamma_response
  family$preinitialize <- function(G)
    list(y = .cgen_gamma_response(G$y), offxset = G$offset)
  family$evaluate <- .cgen_gamma_eval
  family$ll <- .cgen_gamma_ll
  family$init <- .cgen_gamma_init
  family$ncv <- .cgen_gamma_ncv
  family$dev.resids <- .cgen_gamma_deviance
  family$residuals <- .cgen_gamma_residuals
  family$null.deviance <- .cgen_gamma_null_deviance
  .inlaws_ncv_family(.inlaws_diagnostic_family(family, "cgen_gammalss", observe))
}

.cgen_gamma_response <- function(y) {
  r <- .cg_response(y)
  structure(r$y, censor = r$censor)
}

.cgen_gamma_ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
  d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL, D = NULL, eta = NULL,
  ncv = FALSE, sandwich = FALSE) {
  y <- .cgen_gamma_response(y)
  zero <- rep_len(if (is.null(wt)) 1 else wt, length(y)) == 0
  # Neutralize both censoring bounds before the shared likelihood wrapper
  # replaces unused observations. Otherwise their interval could become invalid.
  attr(y, "censor")[zero] <- y[zero] <- 1
  .gen_gamma_ll(y, X, coef, wt, family, offset, deriv, d1b, d2b, Hp,
    rank, fh, D, eta, ncv, sandwich)
}

.cgen_gamma_init <- function(y, x, E, wt, offset, start) {
  r <- .cg_response(y); pseudo <- r$y
  left <- !r$exact & r$lo == 0 & is.finite(r$hi)
  interval <- !r$exact & r$lo > 0 & is.finite(r$hi)
  pseudo[left] <- r$hi[left]/2
  pseudo[interval] <- r$lo[interval]+(r$hi[interval]-r$lo[interval])/2
  pseudo[pseudo == 0] <- 1
  .gen_gamma_init(pseudo, x, E, wt, offset, start)
}

.cgen_gamma_ncv <- function(X, y, wt, nei, beta, family, llf, H = NULL,
  Hi = NULL, R = NULL, offset = NULL, dH = NULL, db = NULL, deriv = FALSE, nt = 1) {
  # mgcv subsets y without preserving attributes. Encode row indices and
  # restore both values and bounds in every held-out likelihood evaluation.
  response <- .cgen_gamma_response(y); bounds <- attr(response, "censor")
  likelihood <- family$ll
  family$ll <- function(y, X, coef, wt, family, ..., eta = NULL) {
    yy <- response[y]; attr(yy, "censor") <- bounds[y]
    likelihood(yy, X, coef, wt, family, ..., eta = eta)
  }
  .gen_gamma_ncv(X, seq_along(response), wt, nei, beta, family, llf,
    H, Hi, R, offset, dH, db, deriv, nt)
}

.cgen_gamma_deviance <- function(y, mu, wt) {
  r <- .cg_response(y); y <- r$y; wt <- rep_len(wt, length(y)); keep <- wt > 0
  out <- numeric(length(y))
  if (any(keep)) {
    yy <- r$y[keep]; attr(yy, "censor") <- r$censor[keep]
    sat <- .cgen_gamma_saturated(yy, mu[keep, , drop = FALSE])$loglik
    ell <- .cgen_gamma_eval(yy, log(mu[keep, , drop = FALSE]))$ell[, 1]
    out[keep] <- pmax(0, 2*wt[keep]*(sat-ell))
  }
  out
}

.cgen_gamma_residuals <- function(object, type = c("deviance", "pearson", "response", "scaled.pearson"), ...) {
  type <- match.arg(type); mu <- object$fitted.values; w <- object$prior.weights
  r <- .cg_response(object$y)
  if (type == "deviance") return(.cgen_gamma_saturated(object$y, mu)$sign *
    sqrt(.cgen_gamma_deviance(object$y, mu, w)))
  out <- .gen_gamma_residuals(object, type)
  out[!r$exact] <- NA_real_
  out
}

.cgen_gamma_null_deviance <- function(object, offset) {
  keep <- object$prior.weights > 0; r <- .cg_response(object$y)
  y <- r$y[keep]; attr(y, "censor") <- r$censor[keep]
  w <- object$prior.weights[keep]; mu <- object$fitted.values[keep, , drop = FALSE]
  sat <- .cgen_gamma_saturated(y, mu)$loglik; eta <- log(mu)
  off <- if (length(offset) && !is.null(offset[[1]]))
    rep_len(offset[[1]], length(keep))[keep] else numeric(length(y))
  objective <- function(a) {
    eta[, 1] <- a+off
    ell <- .cgen_gamma_eval(y, eta)$ell[, 1]
    ans <- sum(2*w*(sat-ell))
    if (is.finite(ans)) ans else .Machine$double.xmax/100
  }
  if (isTRUE(attr(object$pterms[[1]], "intercept") == 0)) return(max(0, objective(0)))
  initial <- stats::weighted.mean(eta[, 1]-off, w)
  max(0, stats::optim(initial, objective, method = "BFGS")$value)
}
