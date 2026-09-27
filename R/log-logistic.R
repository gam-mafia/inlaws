#' Log-logistic arithmetic mean and scale family
#'
#' A positive response family with a log-linked arithmetic mean and a global
#' logistic scale for the log response.
#'
#' @param theta Logistic scale of `log(Y)`, strictly between zero and one.
#'   `NULL` or zero estimates it starting at 0.25; a positive value fixes it;
#'   a negative value estimates it starting at `abs(theta)`.
#' @param link Only `"log"` is supported.
#' @details
#' Write `s` for the log-response logistic scale and `mu` for the arithmetic
#' mean. Then `log(Y)` has a logistic distribution with scale `s` and location
#' `log(mu) + log(sin(pi*s)/(pi*s))`. Conventional log-logistic shape is `1/s`;
#' conventional scale (the response median) is `mu*sin(pi*s)/(pi*s)`.
#' The restriction `0 < s < 1` ensures a finite arithmetic mean.
#' Variance is `mu^2 * (tan(pi*s)/(pi*s) - 1)` for `s < 0.5`, and infinite
#' otherwise. `s` is neither the response standard deviation nor the
#' standard deviation of `log(Y)` (the latter is `pi*s/sqrt(3)`).
#'
#' Scale is optimized on the logit scale: `getTheta()` returns `qlogis(s)`
#' and `getTheta(TRUE)` returns `s`. Weights multiply likelihood contributions;
#' they do not change the response distribution. Zero weights are allowed.
#' Responses must be finite positive numeric vectors; censoring is not supported.
#' Missing rows are handled by the fitter's `na.action`.
#'
#' Supported routes are `gam()` with REML, ML or NCV, ordinary `bam()` with
#' REML, ML or fREML, and `bam(discrete = TRUE)` with fREML or NCV. Discrete
#' NCV requires mgcv >= 1.9-4. Ordinary REML and ML do not use discretization;
#' `bam()` NCV requires discretization. NCV accepts grouped neighbourhoods
#' through `nei`. `gamm()` and unlisted smoothing criteria are not supported.
#'
#' Response predictions are plug-in conditional arithmetic means, with
#' delta-method standard errors. Deviance holds scale fixed and maximizes
#' each observation's likelihood over its mean. Deviance residuals are signed
#' relative to the fitted median. Pearson residuals are `NA` when variance
#' is infinite; response and deviance residuals remain available.
#'
#' Distribution callbacks are `rd(mu, wt = 1, scale = 1)`,
#' `qf(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)`, and
#' `cdf(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE)`.
#' They use the current family scale and the supplied response-scale mean.
#' `rd` draws one observation per mean. `wt` and external `scale` are ignored.
#' CDF and quantile support upper tails and log probabilities; `logp` in the
#' CDF follows mgcv's spelling. They do not draw uncertain model parameters.
#' @return An `extended.family` object for [mgcv::gam()] or [mgcv::bam()].
#' @seealso [log_logisticls()], [stats::Logistic]
#' @export
#' @examples
#' set.seed(21)
#' dat <- data.frame(x = runif(200))
#' dat$y <- log_logistic(0.25)$rd(exp(sin(6*dat$x)))
#' fit <- mgcv::gam(y ~ s(x, k = 7), data = dat,
#'                  family = log_logistic(), method = "REML")
#' fit$family$getTheta(TRUE)
#' fit$family$qf(0.95, fitted(fit))
log_logistic <- function(theta = NULL, link = "log") {
  if (!identical(link, "log")) stop('log_logistic supports only link = "log"')
  if (!is.null(theta) && (!is.numeric(theta) || length(theta) != 1L ||
      !is.finite(theta) || abs(theta) >= 1)) stop("theta must be NULL or a finite number strictly between -1 and 1")
  estimated <- is.null(theta) || theta <= 0
  .Theta <- if (is.null(theta) || theta == 0) stats::qlogis(0.25) else stats::qlogis(abs(theta))
  getTheta <- function(trans = FALSE) if (trans) stats::plogis(.Theta) else .Theta
  putTheta <- function(theta) {
    if (length(theta) != 1L || !is.finite(theta) ||
        stats::plogis(theta) <= 0 || stats::plogis(theta) >= 1)
      stop("invalid log-logistic logit scale")
    .Theta <<- theta
  }
  dev.resids <- function(y, mu, wt, theta = NULL) {
    if (!length(theta)) theta <- getTheta()
    wt <- rep_len(wt, length(y))
    y[wt == 0] <- mu[wt == 0] <- 1
    r <- .llog_eval(y, mu, theta)
    d <- -2*wt*(r$ell[, 1] - r$sat[, 1])
    d[wt == 0] <- 0
    pmax(d, 0)
  }
  Dd <- function(y, mu, theta, wt, level = 0) {
    if (!length(theta)) theta <- getTheta()
    wt <- rep_len(wt, length(y))
    y[wt == 0] <- mu[wt == 0] <- 1
    r <- .llog_eval(y, mu, theta, if (level == 0) 2L else if (level == 1) 3L else 4L)
    d <- -2*wt*(r$ell - r$sat)
    d[wt == 0, ] <- 0
    take <- function(a, b) d[, which(r$powers[, 1] == a & r$powers[, 2] == b)] * factorial(a)*factorial(b)
    out <- list(Dmu = take(1, 0), Dmu2 = take(2, 0),
      EDmu2 = 2*wt/(3*mu^2*stats::plogis(theta)^2))
    if (level > 0) out <- c(out, list(Dth = take(0, 1), Dmuth = take(1, 1),
      Dmu3 = take(3, 0), Dmu2th = take(2, 1),
      EDmu3 = -2*out$EDmu2/mu,
      EDmu2th = -2*stats::plogis(-theta)*out$EDmu2))
    if (level > 1) out <- c(out, list(Dmu4 = take(4, 0), Dth2 = take(0, 2),
      Dmuth2 = take(1, 2), Dmu2th2 = take(2, 2), Dmu3th = take(3, 1)))
    out
  }
  ls <- function(y, w, theta, scale) {
    if (!length(theta)) theta <- getTheta()
    if (length(scale) != 1L || !is.finite(scale) || scale != 1)
      stop("log_logistic likelihood scale is fixed at one")
    w <- rep_len(w, length(y))
    y[w == 0] <- 1
    r <- .llog_eval(y, rep(1, length(y)), theta, 2L)
    take <- function(b) {
      v <- w*r$sat[, which(r$powers[, 1] == 0 & r$powers[, 2] == b)]*factorial(b)
      v[w == 0] <- 0
      v
    }
    list(ls = sum(take(0)), lsth1 = sum(take(1)),
         LSTH1 = matrix(take(1), ncol = 1), lsth2 = matrix(sum(take(2)), 1, 1))
  }
  aic <- function(y, mu, theta = NULL, wt, dev) {
    if (!length(theta)) theta <- getTheta()
    wt <- rep_len(wt, length(y))
    y[wt == 0] <- mu[wt == 0] <- 1
    v <- -2*wt*.llog_eval(y, mu, theta)$ell[, 1]
    v[wt == 0] <- 0
    sum(v)
  }
  initialize <- expression({
    y <- family$prepare.response(y)
    if (any(!is.finite(weights) | weights < 0) || !any(weights > 0))
      stop("log_logistic requires finite non-negative weights with positive total")
    mustart <- as.numeric(y)
  })
  preinitialize <- function(y, family) {
    if (!isTRUE(family$.llog_isolated)) family <- unserialize(serialize(family, NULL))
    list(family = family)
  }
  postproc <- function(family, y, prior.weights, fitted, linear.predictors, offset, intercept) {
    null <- utils::getFromNamespace("find.null.dev", "mgcv")
    list(null.deviance = null(family, .llog_response(y), linear.predictors,
      offset, prior.weights),
      deviance = sum(family$dev.resids(y, fitted, prior.weights)), family = paste0("log_logistic(", round(family$getTheta(TRUE), 3), ")"))
  }
  residuals <- function(object, type = "deviance", ...) {
    type <- match.arg(type, c("deviance", "pearson", "scaled.pearson", "working", "response"))
    if (type == "working") return(object$residuals)
    y <- object$y
    mu <- object$fitted.values; w <- object$prior.weights
    if (type == "deviance") {
      loc <- log(mu) + .llog_correction(object$family$getTheta(TRUE))
      return(sign(log(y)-loc)*sqrt(pmax(0, object$family$dev.resids(y, mu, w))))
    }
    ans <- y - mu
    if (type != "response") {
      v <- object$family$variance(mu)
      ans <- ans*sqrt(w/v)
      ans[!is.finite(v)] <- NA_real_
    }
    ans
  }

  lk <- stats::make.link(link)
  structure(c(list(family = "log_logistic", link = link), lk[c("linkfun", "linkinv", "mu.eta", "valideta")],
    list(validmu = function(mu) all(is.finite(mu) & mu > 0),
      dev.resids = dev.resids, Dd = Dd, ls = ls, aic = aic,
      initialize = initialize, prepare.response = .llog_response,
      preinitialize = preinitialize, postproc = postproc, residuals = residuals,
      n.theta = as.integer(estimated), ini.theta = .Theta,
      getTheta = getTheta, putTheta = putTheta, no.r.sq = TRUE,
      scale = 1,
      variance = function(mu) .llog_variance(mu, getTheta(TRUE)),
      rd = function(mu, wt = 1, scale = 1) .llog_random(mu, getTheta(TRUE)),
      qf = function(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)
        .llog_quantile(p, mu, getTheta(TRUE), lower.tail, log.p),
      cdf = function(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE)
        .llog_cdf(q, mu, getTheta(TRUE), logp, lower.tail))),
    class = c("log_logistic.family", "extended.family", "family"))
}

# Clone before bam captures the family, as for pig.family.
#' @export
#' @importFrom mgcv fix.family.link
fix.family.link.log_logistic.family <- function(fam) {
  fam <- unserialize(serialize(fam, NULL))
  fam$.llog_isolated <- TRUE
  NextMethod()
}
