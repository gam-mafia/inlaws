#' Poisson inverse Gaussian family for count responses
#'
#' An extended family with a log link for the mean and one fixed or estimated
#' global dispersion parameter, for [mgcv::gam()] and [mgcv::bam()].
#' @param theta Dispersion: `NULL` or zero estimates it starting at one;
#'   a positive finite value fixes it; a negative finite value estimates it
#'   starting at `abs(theta)`.
#' @param link Only `"log"` is supported.
#' @return An `extended.family` object, with `rd`, `cdf`, and `qf` callbacks.
#' @details
#' Conditional on Z, Y is Poisson with mean `mu * Z`, where Z is inverse
#' Gaussian with mean one and variance theta. Consequently the response mean
#' is mu and its variance is `mu + theta * mu^2`. There is one mean predictor;
#' dispersion is global. `getTheta(TRUE)` returns theta and `getTheta()`
#' returns log(theta).
#'
#' Supported smoothing methods are REML, ML, and NCV for `gam()`; REML, ML,
#' and fREML for ordinary `bam()`; and fREML and NCV for discrete `bam()`.
#' Discrete NCV requires mgcv >= 1.9-4 and uses a working-model criterion
#' which need not agree with `gam()` NCV. `gamm()` is not supported.
#' The fitting callbacks use variance-based positive working information
#' when mgcv requests expected information; likelihood derivatives are exact.
#' Dispersion can approach zero for Poisson or
#' underdispersed data and may then be weakly identified; the family does not
#' automatically switch to an exact Poisson model.
#'
#' Responses must be finite nonnegative integers, and weights finite and
#' nonnegative. At least one positive count must have positive weight.
#' Weights multiply likelihood contributions; they do not change individual
#' observation distributions. Offsets have the usual log-mean interpretation.
#' Deviance compares with the saturated mean at fixed dispersion, which
#' generally differs from y. Null deviance retains offsets and dispersion.
#'
#' The simulation (`rd(mu)`), CDF (`cdf(q, mu, logp = FALSE)`), and quantile
#' (`qf(p, mu)`) callbacks condition on
#' the family dispersion. Distribution arguments accept scalar recycling or
#' equal lengths. A zero mean is a point mass at zero. For positive means,
#' quantiles at zero and one are zero and infinity. The CDF is right-continuous.
#' CDF and quantile summation uses actual probabilities without truncated
#' renormalization, with a one-million-term limit; exhausting it is an error.
#' Extreme tails remain subject to floating-point accuracy; exceeding the
#' representable numerical range is an error. Simulation uses
#' inverse Gaussian mixing followed by Poisson sampling.
#'
#' Likelihood and analytic derivatives use a half-integer Bessel recurrence.
#' Computation grows with the largest observed count; very large counts may
#' be expensive. Mathematical details are in `inst/maths/pig.md` in the source.
#' @export
#' @examples
#' set.seed(12)
#' dat <- data.frame(x = runif(200))
#' dat$y <- pig(theta = 0.5)$rd(exp(1 + sin(2 * pi * dat$x)))
#' fit <- mgcv::gam(y ~ s(x, k = 6), data = dat,
#'                  family = pig(), method = "REML")
#' fit$family$getTheta(TRUE)
#' predict(fit, type = "response")
pig <- function(theta = NULL, link = "log") {
  if (!identical(link, "log")) stop("pig supports only the log link")
  if (!is.null(theta) && (!is.numeric(theta) || length(theta) != 1L ||
      !is.finite(theta))) stop("theta must be NULL or one finite number")
  fixed <- !is.null(theta) && theta > 0
  initial <- if (is.null(theta) || theta == 0) 0 else log(abs(theta))
  .pig_family(initial, fixed)
}

.pig_family <- function(initial, fixed) {
  .Theta <- initial
  sat_cache <- NULL
  saturated <- function(y, theta) {
    if (is.null(sat_cache) || !identical(sat_cache$y, y) ||
        !identical(sat_cache$theta, theta))
      sat_cache <<- list(y = y, theta = theta, value = .pig_saturated(y, theta))
    sat_cache$value
  }
  getTheta <- function(trans = FALSE) {
    if (trans) exp(.Theta) else .Theta
  }
  putTheta <- function(theta) {
    if (length(theta) != 1L || !is.finite(theta) ||
        !is.finite(exp(theta)) || exp(theta) <= 0) stop("invalid PIG log dispersion")
    .Theta <<- theta
  }
  dev.resids <- function(y, mu, wt, theta = NULL) {
    if (!length(theta)) theta <- .Theta
    wt <- rep_len(wt, length(y)); active <- wt > 0
    ans <- numeric(length(y))
    if (any(active)) {
      ans[active] <- 2 * wt[active] * pmax(0,
        saturated(y[active], theta)[, 1] - .pig_logpmf(y[active], mu[active], theta))
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
    d <- .pig_lld(ya, ma, theta, level)
    map <- c(Dmu = "m", Dmu2 = "mm", Dth = "t", Dmuth = "mt",
             Dmu3 = "mmm", Dmu2th = "mmt", Dmu4 = "mmmm", Dth2 = "tt",
             Dmuth2 = "mtt", Dmu2th2 = "mmtt", Dmu3th = "mmmt")
    for (nm in intersect(names(map), fields)) ans[[nm]][active] <- -2 * wa * d[[map[[nm]]]]
    if (level > 0) {
      sat <- saturated(ya, theta)
      ans$Dth[active] <- ans$Dth[active] + 2 * wa * sat[, 2]
      if (level > 1) ans$Dth2[active] <- ans$Dth2[active] + 2 * wa * sat[, 3]
    }
    # Positive variance-based working information is needed when observed
    # mean curvature is negative (in particular at every zero count).
    v <- ma + exp(theta) * ma^2
    ans$EDmu2[active] <- 2 * wa / v
    if (level > 0) {
      ans$EDmu3[active] <- -2 * wa * (1 + 2 * exp(theta) * ma) / v^2
      ans$EDmu2th[active] <- -2 * wa * exp(theta) * ma^2 / v^2
    }
    ans
  }
  ls <- function(y, w, theta, scale) {
    if (!length(theta)) theta <- .Theta
    if (length(scale) != 1L || !is.finite(scale) || scale != 1)
      stop("pig likelihood scale is fixed at one")
    w <- rep_len(w, length(y)); active <- w > 0
    s <- matrix(0, length(y), 3L)
    if (any(active)) s[active, ] <- saturated(y[active], theta) * w[active]
    list(ls = sum(s[, 1]), lsth1 = sum(s[, 2]),
         LSTH1 = matrix(s[, 2], ncol = 1), lsth2 = sum(s[, 3]))
  }
  aic <- function(y, mu, theta = NULL, wt, dev) {
    if (!length(theta)) theta <- .Theta
    wt <- rep_len(wt, length(y)); active <- wt > 0
    -2 * sum(wt[active] * .pig_logpmf(y[active], mu[active], theta))
  }
  variance <- function(mu) mu + getTheta(TRUE) * mu^2
  rd <- function(mu, wt = 1, scale = 1) .pig_random(mu, .Theta)
  cdf <- function(q, mu, wt = 1, scale = 1, logp = FALSE) {
    lp <- .pig_probability(q, mu, .Theta)
    if (logp) lp else exp(lp)
  }
  qf <- function(p, mu, wt = 1, scale = 1) .pig_probability(p, mu, .Theta, TRUE)
  initialize <- expression({
    if (!is.numeric(y) || is.matrix(y) || any(!is.finite(y)) ||
        any(y < 0 | y != floor(y))) stop("pig requires finite nonnegative integer responses")
    if (any(!is.finite(weights)) || any(weights < 0) || !is.finite(sum(weights)))
      stop("pig requires finite nonnegative weights")
    if (!any(weights > 0 & y > 0)) stop("pig needs a positive count with positive weight")
    mustart <- y + (y == 0) / 6
  })
  preinitialize <- function(y, family) {
    # gam calls this before fixing links; bam has already cloned the family
    # in fix.family.link. A second clone there would leave bam's final AIC
    # calculation referring to the original, unestimated dispersion.
    if (!isTRUE(family$.pig_isolated))
      family <- unserialize(serialize(family, NULL))
    list(family = family)
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
        if (width > 128) stop("pig null deviance optimization failed")
      }
      opt$objective
    }
    # bam may retain deviance from before its final dispersion update.
    list(null.deviance = null,
         deviance = sum(family$dev.resids(y, fitted, prior.weights)))
  }
  lk <- stats::make.link("log")
  structure(list(family = "Poisson inverse Gaussian", link = "log",
    linkfun = lk$linkfun, linkinv = lk$linkinv, mu.eta = lk$mu.eta,
    valideta = lk$valideta, validmu = function(mu) all(is.finite(mu) & mu > 0),
    dev.resids = dev.resids, Dd = Dd, ls = ls, aic = aic,
    initialize = initialize, preinitialize = preinitialize, postproc = postproc,
    variance = variance, n.theta = as.integer(!fixed), ini.theta = initial,
    getTheta = getTheta, putTheta = putTheta, scale = 1, rd = rd, cdf = cdf, qf = qf),
    class = c("pig.family", "extended.family", "family"))
}

# Clone at the public mgcv link-fixing dispatch, before bam retains separate
# local and G$family references. Preserve mgcv's NCV fields and delegate all
# actual link derivatives to its extended-family implementation.
#' @export
#' @importFrom mgcv fix.family.link
fix.family.link.pig.family <- function(fam) {
  fam <- unserialize(serialize(fam, NULL))
  fam$.pig_isolated <- TRUE
  NextMethod()
}
