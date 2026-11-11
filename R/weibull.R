#' Weibull arithmetic mean and global shape family
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Model the arithmetic mean of a positive Weibull response with a log-linked
#' additive predictor and a single estimated or fixed shape.
#' @param shape Weibull shape. `NULL` or zero estimates shape starting at one.
#'   A positive value fixes shape. A negative value estimates shape starting
#'   at its absolute value.
#' @param link The mean link, `"log"`.
#' @return An `extended.family` object for [mgcv::gam()] or [mgcv::bam()].
#' @details For arithmetic mean `mu` and Weibull shape `k`, the conventional
#' scale is `lambda = mu / Gamma(1 + 1/k)`. The response variance is
#' `mu^2 * (Gamma(1 + 2/k) / Gamma(1 + 1/k)^2 - 1)`. At `k = 1` the
#' distribution is exponential with mean `mu`. Responses must be finite,
#' positive, uncensored numeric vectors. Weights multiply likelihood
#' contributions, with zero weights allowed. The additional likelihood scale
#' in mgcv is one. `getTheta()` returns log shape and `getTheta(TRUE)` returns
#' shape on its original scale.
#'
#' In `gam()`, REML supports outer Newton, BFGS, `nlm`, `optim` and EFS.
#' ML with positive weights supports outer Newton. On mgcv through 1.9-5, omit
#' zero-weight rows before using `gam()` ML because of an upstream
#' extended-family ML penalty issue. NCV and QNCV use outer BFGS,
#' including leave-one-out and grouped neighbourhoods. mgcv can retain
#' the requested optimizer in the fitted object when it uses BFGS for NCV.
#'
#' Ordinary `bam()` supports REML, ML and fREML. Discrete `bam()` supports
#' fREML and, with mgcv 1.9-4 or later, NCV. Discrete NCV uses mgcv's
#' neighbourhood criterion on the PIRLS working model. Set `scale = 1`
#' in `bam()` for correct scale and log-likelihood reporting on mgcv
#' 1.9-4 and 1.9-5. The default scale also fits the likelihood with scale one,
#' but mgcv reports an estimated scale and adds one to `logLik()`.
#' `gamm()` is outside the supported scope.
#'
#' Response predictions are plug-in conditional means with delta-method
#' standard errors. Deviance maximizes each observation's likelihood over its
#' mean with shape held fixed. Its residual sign compares the response to
#' `lambda`. Pearson residuals use the Weibull response variance and
#' square-root case weights.
#'
#' Distribution callbacks are `rd(mu, wt = 1, scale = 1)`,
#' `cdf(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE)`, and
#' `qf(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)`.
#' They use the fitted family shape and response-scale means. `rd` draws
#' one observation per mean. Weights and extra scale do not alter the response
#' distribution. Native gratia diagnostics use these same parameters.
#'
#' For a shape predictor use [weibulls()]. At matched parameters the two
#' families give identical likelihoods. Estimating shape as a family parameter
#' or as a regression coefficient can give different REML estimates.
#' @seealso [weibulls()], [stats::Weibull], [family_diagnostics]
#' @template ncv-neighbourhoods-extended
#' @export
#' @examples
#' set.seed(79)
#' d <- data.frame(x = runif(200))
#' d$y <- weibull(2)$rd(exp(sin(6 * d$x)))
#' b <- mgcv::gam(y ~ s(x, k = 6), data = d, family = weibull(), method = "REML")
#' b$family$getTheta(TRUE)
#' head(b$family$qf(0.95, fitted(b)))
#' fast <- mgcv::bam(y ~ s(x, k = 6), data = d, family = weibull(),
#'                   method = "fREML", discrete = TRUE, scale = 1)
weibull <- function(shape = NULL, link = "log") {
  if (!identical(link, "log")) stop('weibull requires the "log" link', call. = FALSE)
  if (!is.null(shape) && (!is.numeric(shape) || length(shape) != 1L ||
      !is.finite(shape)))
    stop("shape must be NULL or a finite numeric scalar", call. = FALSE)
  estimated <- is.null(shape) || shape <= 0
  .Theta <- if (is.null(shape) || shape == 0) 0 else log(abs(shape))
  getTheta <- function(trans = FALSE) if (trans) exp(.Theta) else .Theta
  putTheta <- function(theta) {
    if (length(theta) != 1L || !is.finite(exp(theta)) || exp(theta) <= 0)
      stop("Weibull shape is outside floating-point range", call. = FALSE)
    .Theta <<- theta
  }
  dev.resids <- function(y, mu, wt, theta = NULL) {
    if (!length(theta)) theta <- getTheta()
    .weibull_deviance(y, mu, exp(theta), wt)
  }
  Dd <- function(y, mu, theta, wt, level = 0) {
    if (!length(theta)) theta <- getTheta()
    .weibull_Dd(y, mu, theta, wt, level)
  }
  ls <- function(y, w, theta, scale) {
    if (length(scale) != 1L || !is.finite(scale) || scale != 1)
      stop("Weibull likelihood scale is fixed at one", call. = FALSE)
    if (!length(theta)) theta <- getTheta()
    w <- rep_len(w, length(y))
    keep <- w > 0
    list(ls = sum(w[keep] * (theta - log(y[keep]) - 1)), lsth1 = sum(w),
      LSTH1 = matrix(w, ncol = 1L), lsth2 = matrix(0, 1L, 1L))
  }
  aic <- function(y, mu, theta = NULL, wt, dev) {
    if (!length(theta)) theta <- getTheta()
    w <- rep_len(wt, length(y))
    keep <- w > 0
    -2 * sum(w[keep] * .weibull_eval(y[keep], NULL, theta, mu = mu[keep])[, 1L])
  }
  postproc <- function(family, y, prior.weights, fitted, linear.predictors, offset, intercept) {
    k <- family$getTheta(TRUE)
    list(null.deviance = .weibull_null_deviance(y, k, prior.weights, offset, intercept),
      deviance = sum(family$dev.resids(y, fitted, prior.weights)),
      family = paste0("weibull(", signif(k, 4L), ")"))
  }
  residuals <- function(object, type = c("deviance", "pearson", "scaled.pearson",
                                          "working", "response"), ...) {
    type <- match.arg(type)
    if (type == "working") return(object$residuals)
    .weibull_residuals(object$y, object$fitted.values, object$family$getTheta(TRUE),
                       object$prior.weights, type)
  }
  lk <- .weibull_log_link()
  family <- structure(c(list(family = "weibull", link = link),
    lk[c("linkfun", "linkinv", "mu.eta", "valideta")], list(
      validmu = function(mu) all(is.finite(mu) & mu > 0),
      dev.resids = dev.resids, Dd = Dd, ls = ls, aic = aic,
      prepare.response = .weibull_response, check.weights = .weibull_check_weights,
      initialize = expression({
        y <- family$prepare.response(y)
        family$check.weights(weights)
        mustart <- as.numeric(y)
      }),
      postproc = postproc, residuals = residuals, n.theta = as.integer(estimated),
      ini.theta = .Theta, getTheta = getTheta, putTheta = putTheta,
      no.r.sq = TRUE, scale = 1,
      variance = function(mu) exp(.weibull_log_variance(mu, getTheta(TRUE))),
      rd = function(mu, wt = 1, scale = 1) .weibull_random(mu, getTheta(TRUE)),
      cdf = function(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE)
        .weibull_probability(q, mu, getTheta(TRUE), lower.tail = lower.tail, log.p = logp),
      qf = function(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)
        .weibull_probability(p, mu, getTheta(TRUE), quantile = TRUE,
                             lower.tail = lower.tail, log.p = log.p))),
    class = c("extended.family", "family"))
  family <- .inlaws_extended_family(.inlaws_diagnostic_family(family, "weibull"))
  family$preinitialize <- .weibull_preinitialize
  family
}

.weibull_preinitialize <- function(y, family) {
  caller <- parent.frame()
  scale <- get0("scale", envir = caller, inherits = FALSE)
  if (length(scale) && scale > 0 && scale != 1)
    stop("Weibull likelihood scale is fixed at one", call. = FALSE)
  criterion <- get0("criterion", envir = caller, inherits = FALSE)
  weights <- attr(family, "G")$w
  attr(family, "G") <- NULL
  # The extended-family ML penalty takes a faulty zero-curvature path in
  # mgcv through 1.9-5. Ordinary bam applies ML to a working Gaussian fit.
  if (identical(criterion, "ML") && any(weights == 0) &&
      utils::packageVersion("mgcv") <= "1.9.5")
    stop("weibull gam() ML with zero weights requires omitting zero-weight rows (mgcv ML penalty limitation)",
         call. = FALSE)
  .inlaws_extended_preinitialize(y, family, guard_ncv = FALSE)
}
attr(.weibull_preinitialize, "needG") <- TRUE

# The same analytic partials, with d/dmu = (1/mu) d/dlog(mu).
.weibull_Dd <- function(y, mu, theta, wt, level) {
  w <- rep_len(wt, length(y))
  keep <- w > 0
  m <- mu[keep]
  order <- min(level + 2L, 4L)
  d <- 2 * (w[keep] * .weibull_eval(y[keep], NULL, theta, order, m, deviance = TRUE))
  full <- function(z, power = 0L) {
    # Divide successively so that mu^power itself need not be representable.
    for (j in seq_len(power)) z <- z / m
    out <- numeric(length(y))
    out[keep] <- z
    out
  }
  info <- exp(log(2) + log(w[keep]) + 2 * (theta - log(m)))
  out <- list(Dmu = full(d[, 2L], 1L), Dmu2 = full(d[, 4L] - d[, 2L], 2L),
    EDmu2 = full(info))
  if (level > 0L) out <- c(out, list(Dth = full(d[, 3L]), Dmuth = full(d[, 5L], 1L),
    Dmu3 = full(d[, 7L] - 3 * d[, 4L] + 2 * d[, 2L], 3L),
    Dmu2th = full(d[, 8L] - d[, 5L], 2L),
    EDmu3 = full(-2 * info, 1L), EDmu2th = full(2 * info)))
  if (level > 1L) out <- c(out, list(Dmu4 = full(
      d[, 11L] - 6 * d[, 7L] + 11 * d[, 4L] - 6 * d[, 2L], 4L),
    Dth2 = full(d[, 6L]), Dmuth2 = full(d[, 9L], 1L),
    Dmu2th2 = full(d[, 13L] - d[, 9L], 2L),
    Dmu3th = full(d[, 12L] - 3 * d[, 8L] + 2 * d[, 5L], 3L)))
  out
}
