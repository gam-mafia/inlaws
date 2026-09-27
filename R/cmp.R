#' Mean-parametrised Conway-Maxwell-Poisson family
#'
#' An extended family for [mgcv::gam()] with a log link for the exact mean
#' and a fixed or estimated global dispersion parameter.
#'
#' @param theta Dispersion nu: `NULL` or zero estimates nu starting at one;
#'   a positive value fixes nu; a negative value estimates nu starting at
#'   `abs(theta)`.
#' @param link Only `"log"` is supported.
#' @param control Named list with `sum_tol` (relative tail tolerance, 1e-10),
#'   `mean_tol` (relative mean inversion tolerance, 1e-10), `max_terms`
#'   (maximum summation terms, 100000), and `max_iter` (maximum mean inversion
#'   iterations, 100). Unknown controls are errors.
#' @details The probability at integer y >= 0 is proportional to
#'   lambda^y / (y!)^nu. The mean is modelled directly; lambda is obtained by
#'   numerical inversion and generally differs from the mean. nu = 1 is
#'   Poisson, nu > 1 is underdispersed, and 0 < nu < 1 is overdispersed.
#'   The boundary nu = 0 is excluded. The additional likelihood scale is one.
#'   Supplying another scale to `gam()` is an error. An all-zero response
#'   among observations with positive weights has no finite log-mean fit
#'   and is rejected.
#'
#'   Use `method = "REML"`. This pure-R reference implementation uses adaptive
#'   summation and analytic Taylor arithmetic for derivatives. Extreme means
#'   or dispersion can exhaust numerical limits, which produces an error.
#'   Large-data performance, `bam()`, and other smoothing-selection methods
#'   have not been validated. Near-deterministic data may not have a finite
#'   dispersion estimate.
#'
#'   Prior weights multiply log likelihood contributions; they do not change
#'   the distribution used by the `rd`, `qf`, and `cdf` callbacks. These
#'   callbacks condition on the fitted mean and dispersion. Probabilities and
#'   quantiles use the numerically truncated distribution; extreme-tail
#'   accuracy is limited by `sum_tol`. `qf(0, ...)` is zero and `qf(1, ...)`
#'   is infinity for positive means. A zero mean is a point mass at zero.
#'
#'   `getTheta(TRUE)` returns nu; `getTheta()` returns log(nu).
#' @return An object of classes `extended.family` and `family`, including
#'   variance, simulation (`rd`), quantile (`qf`) and distribution (`cdf`)
#'   functions.
#' @references Huang, A. (2017). Mean-parametrized Conway-Maxwell-Poisson
#'   regression models for dispersed counts. Statistical Modelling, 17,
#'   359--380. doi:10.1177/1471082X17697749.
#' @export
#' @examples
#' set.seed(17)
#' dat <- data.frame(x = runif(60))
#' dat$y <- rpois(60, exp(0.4 + sin(2 * pi * dat$x)))
#' fit <- mgcv::gam(y ~ s(x, k = 5), data = dat,
#'                  family = cmp(theta = 1), method = "REML")
#' predict(fit, type = "response")
cmp <- function(theta = NULL, link = "log", control = list()) {
  if (!is.character(link) || length(link) != 1L || is.na(link) || link != "log")
    stop("cmp supports only the log link", call. = FALSE)
  if (!is.null(theta) && (!is.numeric(theta) || length(theta) != 1L || !is.finite(theta)))
    stop("theta must be NULL or one finite number", call. = FALSE)
  control <- .cmp_control(control)
  fixed <- !is.null(theta) && theta > 0
  .Theta <- if (is.null(theta) || theta == 0) 0 else log(abs(theta))
  # Bounded caches belong to this family instance, never to a namespace/global.
  cache <- sat_cache <- NULL
  getTheta <- function(trans = FALSE) if (trans) exp(.Theta) else .Theta
  putTheta <- function(theta) {
    if (length(theta) != 1L || !is.finite(theta) || !is.finite(exp(theta)) || exp(theta) <= 0)
      stop("invalid CMP log dispersion", call. = FALSE)
    .Theta <<- as.numeric(theta)
  }
  evaluate <- function(mu, theta, jets = FALSE) {
    mu <- as.numeric(mu); nu <- exp(theta)
    if (is.null(cache) || !identical(cache$mu, mu) || !identical(cache$theta, theta)) {
      cache <<- list(mu = mu, theta = theta, z = .cmp_moments(mu, nu, control))
    }
    if (jets && is.null(cache$jet)) {
      if (any(mu <= 0)) stop("CMP fitting means must be positive", call. = FALSE)
      cache$jet <<- .cmp_lljet(mu, mu, nu, cache$z)
    }
    cache
  }
  saturated <- function(y, theta) {
    if (is.null(sat_cache) || !identical(sat_cache$y, y) || !identical(sat_cache$theta, theta)) {
      u <- unique(y); ans <- matrix(0, length(u), 3)
      pos <- u > 0
      if (any(pos)) {
        z <- .cmp_moments(u[pos], exp(theta), control)
        ll <- .cmp_lljet(u[pos], u[pos], exp(theta), z)$ll
        ans[pos, ] <- ll[, c(1, 3, 6), drop = FALSE]
        ans[pos, 3] <- 2 * ans[pos, 3]
      }
      sat_cache <<- list(y = y, theta = theta, ans = ans[match(y, u), , drop = FALSE])
    }
    sat_cache$ans
  }
  loglik <- function(y, mu, theta) {
    z <- evaluate(mu, theta)$z
    a <- vapply(z, `[[`, numeric(1), "a")
    ya <- y * a; ya[y == 0] <- 0
    ya - exp(theta) * lgamma(y + 1) - vapply(z, `[[`, numeric(1), "logZ")
  }
  dev.resids <- function(y, mu, wt, theta = NULL) {
    if (is.null(theta)) theta <- .Theta
    active <- rep_len(wt, length(y)) > 0
    out <- numeric(length(y))
    if (any(active)) out[active] <- 2 * rep_len(wt, length(y))[active] *
      pmax(0, saturated(y[active], theta)[, 1] - loglik(y[active], mu[active], theta))
    out
  }
  Dd <- function(y, mu, theta, wt, level = 0) {
    wt <- rep_len(wt, length(y))
    # Zero-weight observations contribute nothing, including at extreme means.
    active <- wt > 0
    fields <- c("Dmu", "Dmu2", "EDmu2")
    if (level > 0) fields <- c(fields, "Dth", "Dmuth", "Dmu3", "Dmu2th", "EDmu3", "EDmu2th")
    if (level > 1) fields <- c(fields, "Dmu4", "Dth2", "Dmuth2", "Dmu2th2", "Dmu3th")
    ans <- stats::setNames(lapply(fields, function(x) numeric(length(y))), fields)
    if (!any(active)) return(ans)
    ya <- y[active]; ma <- mu[active]; wa <- wt[active]
    e <- evaluate(ma, theta, TRUE)
    # The expected observed derivatives are evaluated by replacing y with mu.
    el <- e$jet$ll
    ll <- el + (ya - ma) * e$jet$da
    dl <- lgamma(ya + 1) - lgamma(ma + 1)
    # The jet at y=mu uses lgamma(mu+1), whereas E[log(Y!)] is needed
    # only for pure dispersion derivatives. All derivatives involving mu
    # are affine in y and unaffected by this term.
    for (k in 1:4) {
      j <- k * (k + 1) / 2 + k + 1
      ll[, j] <- ll[, j] - exp(theta) * dl / factorial(k)
    }
    sat <- if (level > 0) saturated(ya, theta) else NULL
    getd <- function(i, j, expected = FALSE) {
      col <- (i + j) * (i + j + 1) / 2 + j + 1
      -2 * wa * (if (expected) el[, col] else ll[, col]) * factorial(i) * factorial(j)
    }
    ans$Dmu[active] <- getd(1, 0)
    ans$Dmu2[active] <- getd(2, 0)
    ans$EDmu2[active] <- getd(2, 0, TRUE)
    if (level > 0) {
      ans$Dth[active] <- getd(0, 1) + 2 * wa * sat[, 2]
      ans$Dmuth[active] <- getd(1, 1)
      ans$Dmu3[active] <- getd(3, 0)
      ans$Dmu2th[active] <- getd(2, 1)
      ans$EDmu3[active] <- getd(3, 0, TRUE)
      ans$EDmu2th[active] <- getd(2, 1, TRUE)
    }
    if (level > 1) {
      ans$Dmu4[active] <- getd(4, 0)
      ans$Dth2[active] <- getd(0, 2) + 2 * wa * sat[, 3]
      ans$Dmuth2[active] <- getd(1, 2)
      ans$Dmu2th2[active] <- getd(2, 2)
      ans$Dmu3th[active] <- getd(3, 1)
    }
    ans
  }
  ls <- function(y, w, theta, scale) {
    if (length(scale) != 1L || !is.finite(scale) || scale != 1)
      stop("CMP likelihood scale is fixed at one", call. = FALSE)
    active <- rep_len(w, length(y)) > 0
    s <- matrix(0, length(y), 3)
    if (any(active)) s[active, ] <- saturated(y[active], theta) * rep_len(w, length(y))[active]
    list(ls = sum(s[, 1]), lsth1 = sum(s[, 2]),
         LSTH1 = matrix(s[, 2], ncol = 1), lsth2 = sum(s[, 3]))
  }
  aic <- function(y, mu, theta = NULL, wt, dev) {
    if (is.null(theta)) theta <- .Theta
    active <- rep_len(wt, length(y)) > 0
    -2 * sum(rep_len(wt, length(y))[active] * loglik(y[active], mu[active], theta))
  }
  variance <- function(mu) vapply(evaluate(mu, .Theta)$z, `[[`, numeric(1), "variance")
  qf <- function(p, mu, wt = 1, scale = 1) {
    v <- .cmp_recycle(p, mu); p <- v[[1]]; mu <- v[[2]]
    if (anyNA(p) || any(p < 0 | p > 1)) stop("CMP probabilities must lie in [0, 1]", call. = FALSE)
    z <- evaluate(mu, .Theta)$z
    vapply(seq_along(p), function(i) {
      if (p[i] == 0 || mu[i] == 0) return(0)
      if (p[i] == 1) return(Inf)
      cs <- cumsum(z[[i]]$p); cs[length(cs)] <- 1
      z[[i]]$y[which(cs >= p[i])[1]]
    }, numeric(1))
  }
  cdf <- function(q, mu, wt = 1, scale = 1, logp = FALSE) {
    v <- .cmp_recycle(q, mu); q <- v[[1]]; mu <- v[[2]]
    z <- evaluate(mu, .Theta)$z
    ans <- vapply(seq_along(q), function(i) {
      if (is.na(q[i])) return(NA_real_)
      min(1, sum(z[[i]]$p[z[[i]]$y <= floor(q[i])]))
    }, numeric(1))
    if (logp) log(ans) else ans
  }
  rd <- function(mu, wt = 1, scale = 1) qf(stats::runif(length(mu)), mu)
  initialize <- expression({
    if (any(!is.finite(y)) || any(y < 0 | y != floor(y)))
      stop("CMP response must contain nonnegative finite integers")
    if (any(!is.finite(weights)) || any(weights < 0))
      stop("CMP weights must be nonnegative and finite")
    if (!any(weights > 0 & y > 0))
      stop("CMP needs a positive count with positive weight for a finite log-mean fit")
    mustart <- y + (y == 0) / 6
  })
  # gam() mutates family parameter environments. Clone before each fit so
  # reusing a constructor object cannot change a previously fitted model.
  preinitialize <- function(y, family) {
    list(family = unserialize(serialize(family, NULL)))
  }
  postproc <- function(family, y, prior.weights, fitted, linear.predictors, offset, intercept) {
    # Convex in the intercept's log-lambda only for a constant mean. Use a
    # scalar mean-link optimisation to accommodate exposure offsets as well.
    obj <- function(b) sum(family$dev.resids(y, exp(b + offset), prior.weights))
    centre <- sum((linear.predictors - offset) * prior.weights) / sum(prior.weights)
    width <- 1
    repeat {
      opt <- stats::optimize(obj, c(centre - width, centre + width))
      if (abs(opt$minimum - centre) < 0.95 * width) break
      width <- width * 2
      if (width > 128) stop("CMP null deviance optimisation failed", call. = FALSE)
    }
    list(null.deviance = opt$objective)
  }
  lk <- stats::make.link("log")
  structure(list(family = "Conway-Maxwell-Poisson", link = "log",
    linkfun = lk$linkfun, linkinv = lk$linkinv, mu.eta = lk$mu.eta,
    valideta = lk$valideta, validmu = function(mu) all(is.finite(mu) & mu > 0),
    dev.resids = dev.resids, Dd = Dd, ls = ls, aic = aic,
    initialize = initialize, preinitialize = preinitialize,
    postproc = postproc, variance = variance,
    n.theta = as.integer(!fixed), ini.theta = .Theta, getTheta = getTheta,
    putTheta = putTheta, scale = 1, rd = rd, qf = qf, cdf = cdf),
    class = c("extended.family", "family"))
}
