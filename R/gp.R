#' Generalized Poisson family for count responses
#'
#' A GP-1 extended family for [mgcv::gam()], with a log link for the mean
#' and one fixed or estimated global dispersion parameter.
#'
#' @param theta The variance-to-mean ratio phi. `NULL` or zero estimates phi
#'   starting at two; a positive value at least one fixes phi; a negative
#'   value less than minus one estimates phi starting at `abs(theta)`.
#' @param link Only `"log"` is supported.
#' @return An `extended.family` object for `gam()`, including `rd`, `cdf`
#'   and `qf` distribution callbacks.
#' @details
#' With `a = mu / sqrt(phi)` and `b = 1 - 1 / sqrt(phi)`, the probability
#' of a nonnegative integer y is `a * (a + b*y)^(y-1) * exp(-a-b*y) / y!`.
#' The mean is mu and the variance is `mu * phi`. Only overdispersion is
#' estimated: phi > 1. Fixed `theta = 1` gives the exact Poisson family.
#' Underdispersion and a dispersion linear predictor are not supported.
#'
#' Use `method = "REML"` or `"ML"`. The additional likelihood scale is
#' fixed at one. `getTheta(TRUE)` returns phi; `getTheta()` returns
#' `log(phi - 1)`, or `numeric(0)` for fixed Poisson (no nuisance parameter).
#' Estimation approaches
#' the Poisson boundary without automatically switching models; almost
#' Poisson or underdispersed data can give weak dispersion identification.
#' `bam()` and other smoothing-selection methods have not been validated.
#' In mgcv 1.9-4, ML can fail with a LAPACK error for entirely unpenalized
#' models when individual observations have negative log-link curvature.
#' This is an mgcv zero-column SVD limitation; use REML for such models.
#'
#' Responses must be finite nonnegative integers. At least one positive
#' count must have positive weight for a finite log-mean fit. Prior weights
#' multiply likelihood contributions; they do not alter the distribution
#' used for simulation or quantiles. Offsets have their usual log-mean
#' interpretation. Deviance uses the saturated mean at fixed phi, which
#' generally differs from y. The fitting callbacks use observed information
#' as an estimate of expected information (exact expected information for
#' fixed Poisson).
#'
#' Distribution callbacks condition on fitted parameters. The CDF is
#' right-continuous and accepts `logp = TRUE`. Quantiles at probabilities
#' zero and one are zero and infinity for positive means; a zero mean is
#' a point mass at zero. Simulation uses a Poisson branching process.
#' CDF/quantile summation is limited to one million terms, and simulation
#' to one million generations; exhausting a limit is an error. Probabilities
#' are not renormalized after truncation. Extreme tails remain subject to
#' floating-point accuracy.
#' @references Consul, P. C. and Jain, G. C. (1973). A generalization of the
#'   Poisson distribution. Technometrics, 15, 791--799.
#'   doi:10.1080/00401706.1973.10489112.
#' @export
#' @examples
#' set.seed(12)
#' dat <- data.frame(x = runif(200))
#' dat$y <- gp(theta = 2)$rd(exp(1 + sin(2 * pi * dat$x)))
#' fit <- mgcv::gam(y ~ s(x, k = 6), data = dat,
#'                  family = gp(), method = "REML")
#' fit$family$getTheta(TRUE)
#' predict(fit, type = "response")
gp <- function(theta = NULL, link = "log") {
  if (!identical(link, "log")) stop("gp supports only the log link")
  if (!is.null(theta) && (!is.numeric(theta) || length(theta) != 1L ||
      !is.finite(theta) || !(theta == 0 || theta >= 1 || theta < -1)))
    stop("theta must be NULL, zero, at least one, or less than minus one")
  fixed <- !is.null(theta) && theta >= 1
  initial <- if (is.null(theta) || theta == 0) 0 else log(abs(theta) - 1)
  .gp_family(initial, fixed)
}

.gp_family <- function(initial, fixed) {
  .Theta <- initial
  poisson <- fixed && initial == -Inf
  getTheta <- function(trans = FALSE) {
    if (trans) return(1 + exp(.Theta))
    if (poisson) numeric(0) else .Theta
  }
  putTheta <- function(theta) {
    if (!poisson) .Theta <<- theta
  }
  dev.resids <- function(y, mu, wt, theta = NULL) {
    if (!length(theta)) theta <- .Theta
    wt <- rep_len(wt, length(y)); active <- wt > 0
    ans <- numeric(length(y))
    if (theta == -Inf) {
      ans[active] <- stats::poisson()$dev.resids(y[active], mu[active], wt[active])
    } else if (any(active)) {
      ans[active] <- 2 * wt[active] * pmax(0,
        .gp_saturated(y[active], theta)[, 1] - .gp_logpmf(y[active], mu[active], theta))
    }
    ans
  }
  Dd <- function(y, mu, theta, wt, level = 0) {
    if (!length(theta)) theta <- .Theta
    wt <- rep_len(wt, length(y)); active <- wt > 0
    fields <- c("Dmu", "Dmu2", "EDmu2")
    if (level > 0) fields <- c(fields, "Dth", "Dmuth", "Dmu3", "Dmu2th", "EDmu3", "EDmu2th")
    if (level > 1) fields <- c(fields, "Dmu4", "Dth2", "Dmuth2", "Dmu2th2", "Dmu3th")
    ans <- stats::setNames(lapply(fields, function(x) numeric(length(y))), fields)
    if (!any(active)) return(ans)
    ya <- y[active]; ma <- mu[active]; wa <- wt[active]
    d <- .gp_lld(ya, ma, theta, level)
    map <- c(Dmu = "m", Dmu2 = "mm", Dth = "t", Dmuth = "mt",
             Dmu3 = "mmm", Dmu2th = "mmt", Dmu4 = "mmmm", Dth2 = "tt",
             Dmuth2 = "mtt", Dmu2th2 = "mmtt", Dmu3th = "mmmt")
    for (nm in intersect(names(map), fields)) ans[[nm]][active] <- -2 * wa * d[[map[[nm]]]]
    if (level > 0) {
      sat <- .gp_saturated(ya, theta)
      ans$Dth[active] <- ans$Dth[active] + 2 * wa * sat[, 2]
      if (level > 1) ans$Dth2[active] <- ans$Dth2[active] + 2 * wa * sat[, 3]
    }
    ans$EDmu2 <- ans$Dmu2
    if (level > 0) {
      ans$EDmu3 <- ans$Dmu3
      ans$EDmu2th <- ans$Dmu2th
    }
    if (theta == -Inf) {
      ans$EDmu2[active] <- 2 * wa / ma
      if (level > 0) ans$EDmu3[active] <- -4 * wa / ma^2
    }
    ans
  }
  ls <- function(y, w, theta, scale) {
    if (!length(theta)) theta <- .Theta
    if (length(scale) != 1L || !is.finite(scale) || scale != 1)
      stop("gp likelihood scale is fixed at one")
    w <- rep_len(w, length(y)); active <- w > 0
    s <- matrix(0, length(y), 3L)
    if (any(active)) s[active, ] <- .gp_saturated(y[active], theta) * w[active]
    if (poisson) return(list(ls = sum(s[, 1]), lsth1 = numeric(0),
      LSTH1 = matrix(0, length(y), 0L), lsth2 = matrix(0, 0L, 0L)))
    list(ls = sum(s[, 1]), lsth1 = sum(s[, 2]),
         LSTH1 = matrix(s[, 2], ncol = 1), lsth2 = sum(s[, 3]))
  }
  aic <- function(y, mu, theta = NULL, wt, dev) {
    if (!length(theta)) theta <- .Theta
    wt <- rep_len(wt, length(y)); active <- wt > 0
    -2 * sum(wt[active] * .gp_logpmf(y[active], mu[active], theta))
  }
  variance <- function(mu) mu * getTheta(TRUE)
  rd <- function(mu, wt = 1, scale = 1) .gp_random(mu, .Theta)
  cdf <- function(q, mu, wt = 1, scale = 1, logp = FALSE) {
    lp <- .gp_probability(q, mu, .Theta)
    if (logp) lp else exp(lp)
  }
  qf <- function(p, mu, wt = 1, scale = 1) .gp_probability(p, mu, .Theta, TRUE)
  initialize <- expression({
    if (!is.numeric(y) || is.matrix(y) || any(!is.finite(y)) ||
        any(y < 0 | y != floor(y))) stop("gp requires finite nonnegative integer responses")
    if (any(!is.finite(weights)) || any(weights < 0) || !is.finite(sum(weights)))
      stop("gp requires finite nonnegative weights")
    if (!any(weights > 0 & y > 0)) stop("gp needs a positive count with positive weight")
    mustart <- y + (y == 0) / 6
  })
  preinitialize <- function(y, family) {
    # Rebuild closures for each fit; never retain or modify the template.
    fresh <- .gp_family(initial, fixed)
    list(family = mgcv::fix.family.link(fresh))
  }
  postproc <- function(family, y, prior.weights, fitted, linear.predictors, offset, intercept) {
    obj <- function(b) sum(family$dev.resids(y, exp(b + offset), prior.weights))
    null <- if (!intercept) obj(0) else {
      centre <- stats::weighted.mean(linear.predictors - offset, prior.weights)
      width <- 1
      repeat {
        opt <- stats::optimize(obj, centre + c(-width, width), tol = 1e-9)
        if (abs(opt$minimum - centre) < .95 * width) break
        width <- width * 2
        if (width > 128) stop("gp null deviance optimization failed")
      }
      opt$objective
    }
    list(null.deviance = null)
  }
  lk <- stats::make.link("log")
  structure(list(family = "Generalized Poisson", link = "log",
    linkfun = lk$linkfun, linkinv = lk$linkinv, mu.eta = lk$mu.eta,
    valideta = lk$valideta, validmu = function(mu) all(is.finite(mu) & mu > 0),
    dev.resids = dev.resids, Dd = Dd, ls = ls, aic = aic,
    initialize = initialize, preinitialize = preinitialize, postproc = postproc,
    variance = variance, n.theta = as.integer(!fixed), ini.theta = if (poisson) numeric(0) else initial,
    getTheta = getTheta, putTheta = putTheta, scale = 1, rd = rd, cdf = cdf, qf = qf),
    class = c("extended.family", "family"))
}
