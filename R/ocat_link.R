#' Cumulative-link family for ordered categorical responses
#'
#' A single-predictor ordinal family extending [mgcv::ocat()] to other latent
#' distributions. The cumulative link is selected by `link`; the returned
#' family's location link remains the identity.
#'
#' @param R Integer number of ordered categories, at least two.
#' @param link Cumulative link: `"logit"`, `"probit"`, `"cloglog"`,
#'   `"loglog"`, or `"cauchit"`.
#' @param theta Optional vector of `R - 2` threshold increments, following
#'   [mgcv::ocat()]: all positive fixes the increments; any negative value
#'   supplies starting increments through their absolute values. Zeros are
#'   not allowed. If omitted, thresholds are estimated. `R` may be omitted
#'   when `theta` is supplied.
#' @details
#' For category \eqn{j}, \eqn{P(Y \le j) = F(\alpha_j - \eta)}, where
#' \eqn{F} is the inverse cumulative link, \eqn{\eta} is the smooth location
#' predictor, and the first finite threshold is fixed at -1. Keep the formula
#' intercept. The latent scale is fixed at one. Responses must be integers
#' from 1 to `R`; convert ordered factors explicitly with `as.integer()`.
#'
#' `gam()` supports REML, ML and NCV; `bam()` supports REML, ML and fREML,
#' and (with a sufficiently recent mgcv) NCV with `discrete = TRUE`.
#' Discrete fitting supports fREML and NCV. As for `ocat()`, bam estimates
#' thresholds conditionally within its fitting iteration, rather than by
#' gam's outer smoothing criterion. Bam NCV uses the working linear model.
#' Cauchit likelihoods are not globally concave: use care with convergence.
#' In mgcv 1.9-4, cauchit `gam(method="ML")` fits with no penalized terms
#' can fail in the native solver; `bam(method="ML")` is an alternative
#' (tight convergence tolerances may be needed). Smooth-model ML is supported.
#' `MASS::polr()` clips latent CDF arguments at +/-100 when fitting, which
#' materially affects Cauchy tails, so its cauchit estimates can differ.
#'
#' `predict(..., type = "response")` returns a matrix of category
#' probabilities. Its standard errors condition on the fitted thresholds.
#' `fitted.values` and `predict(..., type = "link")` are latent locations.
#' A positive location effect shifts probability towards larger categories.
#' Only the logit choice is a proportional-odds model.
#'
#' The family includes `cdf(q, mu, wt = 1, scale = 1, lower.tail = TRUE,
#' log.p = FALSE)`, `qf(p, mu, wt = 1, scale = 1, lower.tail = TRUE,
#' log.p = FALSE)`, and `rd(mu, wt = 1, scale = 1)`.
#' Here `mu` is the latent location, not a category probability or mean.
#' Scalars recycle against vectors; other lengths must match. `cdf` gives
#' \eqn{P(Y \le q)} (using the floor of noninteger `q`), `qf` gives the
#' smallest category with cumulative probability at least `p`, and `rd`
#' draws one category per location. The endpoint quantiles are 1 and `R`.
#' Weights are likelihood weights and do not alter these distributions;
#' `wt` and `scale` are accepted for compatibility and ignored by the helpers.
#' Missing inputs propagate; invalid probabilities give NaN with a warning.
#'
#' The derivation, derivative conventions, identification, and numerical
#' details are in `system.file("maths", "ocat_link.md", package="inlaws")`.
#' This implementation is intended for modest numbers of categories:
#' threshold derivative storage grows quadratically in `R`.
#' @return An `extended.family` object for [mgcv::gam()] and [mgcv::bam()].
#' @seealso [mgcv::ocat()], [MASS::polr()]
#' @export
#' @examples
#' set.seed(12)
#' x <- runif(300, -2, 2)
#' y <- as.integer(cut(sin(x) + rnorm(300),
#'                     c(-Inf, -1, 0.5, 1.5, Inf)))
#' fit <- mgcv::gam(y ~ s(x, k = 8), family = ocat_link(R = 4, link = "probit"),
#'                  method = "REML")
#' head(predict(fit, type = "response"))
#' fit$family$getTheta(TRUE)
#' head(fit$family$cdf(2, fitted(fit)))
#' head(fit$family$qf(0.9, fitted(fit)))
#' head(fit$family$rd(fitted(fit)))
ocat_link <- function(R = NULL,
                      link = c("logit", "probit", "cloglog", "loglog", "cauchit"),
                      theta = NULL) {
  link <- match.arg(link)
  if (!is.null(theta)) {
    if (!is.numeric(theta) || !is.null(dim(theta)) ||
        any(!is.finite(theta)) || any(theta == 0))
      stop("theta must contain finite nonzero increments")
    if (!is.null(R) && !identical(as.numeric(R), length(theta) + 2))
      stop("R must equal length(theta) + 2")
    R <- length(theta) + 2L
  }
  if (length(R) != 1L || !is.numeric(R) || !is.finite(R) || R < 2 || R != floor(R))
    stop("R must be an integer of at least two")
  R <- as.integer(R)
  # Reuse mgcv's threshold state and extended-family infrastructure, not its
  # logistic likelihood. All distribution-specific components are replaced.
  fam <- mgcv::ocat(R = R, theta = theta)
  get_theta <- fam$getTheta
  dist <- .ocat_distribution(link)
  cuts <- function(theta) c(-Inf, -1, -1 + cumsum(exp(theta)), Inf)
  fam$cumulative.link <- link
  fam$family <- paste0("Cumulative link (", link, ")")
  fam$scale <- 1
  fam$dev.resids <- function(y, mu, wt, theta = NULL) {
    if (is.null(theta)) theta <- get_theta()
    a <- cuts(theta)
    out <- -2 * wt * .ocat_log_interval(a[y] - mu, a[y + 1L] - mu, dist)
    out[wt == 0] <- 0
    attr(out, "sign") <- sign((a[y] + a[y + 1L])/2 - mu)
    out
  }
  fam$Dd <- function(y, mu, theta, wt = NULL, level = 0) {
    if (is.null(wt)) wt <- rep(1, length(y))
    out <- .ocat_derivatives(y, mu, theta, wt, level, dist)
    if (link %in% c("cloglog", "loglog"))
      out <- .ocat_extreme(out, y, mu, theta, wt, level, link)
    if (link == "cauchit") {
      fisher <- .ocat_fisher(mu, theta, wt, level, dist)
      out[names(fisher)] <- fisher
    }
    out
  }
  fam$aic <- function(y, mu, theta = NULL, wt, dev) {
    if (is.null(theta)) theta <- get_theta()
    a <- cuts(theta)
    z <- wt * .ocat_log_interval(a[y] - mu, a[y + 1L] - mu, dist)
    z[wt == 0] <- 0
    -2 * sum(z)
  }
  fam$preinitialize <- function(y, family) {
    if (!is.numeric(y) || any(!is.finite(y)) || any(y != floor(y)) ||
        any(y < 1 | y > R)) stop("response must be integer category labels from 1 to R")
    if (R > 2L && family$n.theta > 0) {
      if (!is.null(theta)) return(list(Theta = log(abs(theta))))
      p <- cumsum(tabulate(y, nbins = R) + 0.5)/(length(y) + R/2)
      list(Theta = log(pmax(diff(dist$q(p[-R])), 0.01)))
    } else list()
  }
  fam$initialize <- expression({
    a <- c(family$getTheta(TRUE))
    a <- c(a[1L] - 1, a, a[length(a)] + 1)
    mustart <- (a[y] + a[y + 1L])/2
  })
  fam$postproc <- function(family, y, prior.weights, fitted, linear.predictors,
                           offset, intercept) {
    null_dev <- utils::getFromNamespace("find.null.dev", "mgcv")
    list(null.deviance = null_dev(family, y, eta = linear.predictors,
                                 offset = offset, weights = prior.weights),
         family = paste0("Cumulative link (", link, "; ",
                         paste(round(family$getTheta(TRUE), 3), collapse = ","), ")"))
  }
  fam$cdf <- function(q, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) {
    z <- .ocat_recycle(q, mu)
    q <- z[[1L]]; mu <- z[[2L]]
    a <- cuts(get_theta())
    index <- pmin(R, pmax(0, floor(q))) + 1L
    dist$p(a[index] - mu, lower.tail = lower.tail, log.p = log.p)
  }
  fam$qf <- function(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) {
    z <- .ocat_recycle(p, mu)
    p <- z[[1L]]; mu <- z[[2L]]
    bad <- if (log.p) p > 0 else p < 0 | p > 1
    if (any(bad, na.rm = TRUE)) warning("NaNs produced")
    p[which(bad)] <- NaN
    # Compare on the supplied probability scale, retaining exact boundaries
    # and using log CDF/survival probabilities for log-probability inputs.
    lp <- if (log.p) p else log(p)
    a <- get_theta(TRUE)
    ans <- rep(1L, length(p))
    for (aj in a) {
      v <- dist$p(aj - mu, lower.tail = lower.tail, log.p = log.p)
      ans <- ans + if (lower.tail) as.integer(v < p) else as.integer(v > p)
    }
    low <- if (lower.tail) lp == -Inf else lp == 0
    high <- if (lower.tail) lp == 0 else lp == -Inf
    ans[which(low)] <- 1L; ans[which(high)] <- R
    ans[is.na(p) | is.na(mu)] <- p[is.na(p) | is.na(mu)] + mu[is.na(p) | is.na(mu)]
    ans
  }
  qfun <- fam$qf
  fam$rd <- function(mu, wt = 1, scale = 1) qfun(stats::runif(length(mu)), mu)
  fam$predict <- function(family, se = FALSE, eta = NULL, y = NULL, X = NULL,
                          beta = NULL, off = NULL, Vb = NULL) {
    if (!is.null(eta)) return(list(family$qf(0.5, eta)))
    if (is.null(off)) off <- 0
    if (is.list(X)) {
      args <- list(X$Xd, beta, k = X$kd, ks = X$ks, ts = X$ts,
                   dt = X$dt, v = X$v, qc = X$qc, drop = X$drop)
      mu <- off + do.call(utils::getFromNamespace("Xbd", "mgcv"), args)
      if (se) {
        args[[2L]] <- Vb
        vse <- do.call(utils::getFromNamespace("diagXVXd", "mgcv"), args)
      }
    } else {
      mu <- off + drop(X %*% beta)
      if (se) vse <- rowSums((X %*% Vb) * X)
    }
    a <- cuts(family$getTheta())
    prob <- vapply(seq_len(R), function(j)
      exp(.ocat_log_interval(a[j] - mu, a[j + 1L] - mu, dist)), numeric(length(mu)))
    if (!se) return(list(prob))
    dp <- vapply(seq_len(R), function(j)
      exp(dist$d(a[j] - mu)$logd) - exp(dist$d(a[j + 1L] - mu)$logd), numeric(length(mu)))
    list(fit = prob, se.fit = abs(dp) * sqrt(pmax(0, vse)))
  }
  fam$residuals <- function(object, type = c("deviance", "working", "response")) {
    type <- match.arg(type)
    if (type == "working") return(object$residuals)
    if (type == "response") return(object$y - object$family$qf(0.5, object$fitted.values))
    d <- object$family$dev.resids(object$y, object$fitted.values, object$prior.weights)
    attr(d, "sign") * sqrt(pmax(0, d))
  }
  fam
}

.ocat_recycle <- function(x, mu) {
  vector_like <- function(z) is.numeric(z) &&
    (is.null(dim(z)) || length(dim(z)) == 1L || (is.matrix(z) && ncol(z) == 1L))
  if (!vector_like(x) || !vector_like(mu))
    stop("distribution arguments must be numeric vectors")
  x <- as.numeric(x); mu <- as.numeric(mu)
  if (any(is.infinite(mu))) stop("mu must be finite or missing")
  if (!length(x) || !length(mu)) return(list(numeric(), numeric()))
  n <- max(length(x), length(mu))
  if (!length(x) %in% c(1L, n) || !length(mu) %in% c(1L, n))
    stop("distribution arguments must have matching lengths or length one")
  list(rep_len(x, n), rep_len(mu, n))
}

# F, its quantile, and log f together with f'/f, f''/f, f'''/f.
.ocat_distribution <- function(link) {
  if (link == "logit") return(list(p = stats::plogis, q = stats::qlogis,
    d = function(x) {
      u <- stats::plogis(x)
      list(logd = stats::dlogis(x, log = TRUE),
           h = cbind(1, 1 - 2*u, 1 - 6*u + 6*u^2,
                     1 - 14*u + 36*u^2 - 24*u^3))
    }))
  if (link == "probit") return(list(p = stats::pnorm, q = stats::qnorm,
    d = function(x) {
      z <- x; z[is.infinite(z)] <- 0
      list(logd = stats::dnorm(x, log = TRUE), h = cbind(1, -z, z^2 - 1, 3*z - z^3))
    }))
  if (link == "cauchit") return(list(p = stats::pcauchy, q = stats::qcauchy,
    d = function(x) {
      z <- x; z[is.infinite(z)] <- 0
      # Reciprocal formulation avoids powers of large x.
      r <- 1/(1 + z^2); v <- z*r
      list(logd = stats::dcauchy(x, log = TRUE),
           h = cbind(1, -2*v, 8*v^2 - 2*r, 24*v*r - 48*v^3))
    }))
  # loglog is the reflection of cloglog. Stable log-CDFs in both tails.
  p <- function(q, lower.tail = TRUE, log.p = FALSE) {
    z <- if (link == "cloglog") q else -q
    t <- exp(z)
    lf <- log(-expm1(-t))
    lf[z < -35 & !is.na(z)] <- z[z < -35 & !is.na(z)]
    ls <- -t
    out <- if (xor(lower.tail, link == "loglog")) lf else ls
    if (log.p) out else exp(out)
  }
  q <- function(p) {
    if (link == "cloglog") log(-log1p(-p)) else -log(-log(p))
  }
  d <- function(x) {
    z <- if (link == "cloglog") x else -x
    t <- exp(z)
    ld <- z - t
    ld[is.infinite(z)] <- -Inf
    # Infinite t gives vanishing derivatives; avoid 0 * Inf downstream.
    t[!is.finite(t) & !is.na(t)] <- 0
    h <- cbind(1, 1-t, 1-3*t+t^2, 1-7*t+6*t^2-t^3)
    if (link == "loglog") h[, c(2L, 4L)] <- -h[, c(2L, 4L)]
    list(logd = ld, h = h)
  }
  list(p = p, q = q, d = d)
}

.ocat_log_interval <- function(a, b, dist) {
  la <- dist$p(a, log.p = TRUE); lb <- dist$p(b, log.p = TRUE)
  sa <- dist$p(a, lower.tail = FALSE, log.p = TRUE)
  sb <- dist$p(b, lower.tail = FALSE, log.p = TRUE)
  right <- la > log(0.5)
  hi <- ifelse(right, sa, lb); lo <- ifelse(right, sb, la)
  hi + log(-expm1(lo - hi))
}

# Analytic log-likelihood derivatives, using endpoint density/probability ratios.
.ocat_derivatives <- function(y, mu, theta, wt, level, dist) {
  n <- length(y); nt <- length(theta)
  wt <- rep_len(wt, n)
  a <- c(-Inf, -1, -1 + cumsum(exp(theta)), Inf)
  lo <- a[y] - mu; hi <- a[y + 1L] - mu
  lp <- .ocat_log_interval(lo, hi, dist)
  dl <- dist$d(lo); dh <- dist$d(hi)
  l <- dl$h * exp(dl$logd - lp); h <- dh$h * exp(dh$logd - lp)
  l[is.infinite(lo), ] <- 0; h[is.infinite(hi), ] <- 0
  r <- sweep(h - l, 2L, c(-1, 1, -1, 1), `*`)
  r1 <- r[, 1L]; r2 <- r[, 2L]; r3 <- r[, 3L]; r4 <- r[, 4L]
  mul <- function(x) { z <- -2 * wt * x; z[wt == 0] <- 0; z }
  out <- list(D = mul(lp), Dmu = mul(r1), Dmu2 = mul(r2 - r1^2))
  out$EDmu2 <- out$Dmu2
  if (level < 1L) return(out)
  out$Dmu3 <- mul(r3 - 3*r1*r2 + 2*r1^3)
  s <- array(0, c(n, nt, 4L))
  up <- down <- matrix(0, n, nt)
  for (j in seq_len(nt)) {
    up[, j] <- exp(theta[j]) * (y >= j + 1L & y < nt + 2L)
    down[, j] <- exp(theta[j]) * (y >= j + 2L)
    s[, j, ] <- sweep(h * up[, j] - l * down[, j], 2L, c(1, -1, 1, -1), `*`)
  }
  sm <- lapply(seq_len(4L), function(k) matrix(s[, , k], n, nt))
  s0 <- sm[[1L]]; s1 <- sm[[2L]]; s2 <- sm[[3L]]; s3 <- sm[[4L]]
  out$Dth <- mul(s0)
  out$Dmuth <- mul(s1 - r1*s0)
  out$Dmu2th <- mul(s2 - 2*r1*s1 + (2*r1^2-r2)*s0)
  out$EDmu2th <- out$Dmu2th
  if (level < 2L) return(out)
  out$Dmu4 <- mul(r4 - 4*r1*r3 - 3*r2^2 + 12*r1^2*r2 - 6*r1^4)
  out$Dmu3th <- mul(s3 - 3*r1*s2 + (6*r1^2-3*r2)*s1 +
                     (-r3+6*r1*r2-6*r1^3)*s0)
  out$Dth2 <- out$Dmuth2 <- out$Dmu2th2 <- matrix(0, n, nt*(nt+1L)/2L)
  col <- 0L
  for (j in seq_len(nt)) for (k in j:nt) {
    col <- col + 1L
    tt <- matrix(0, n, 3L)
    for (m in 0:2) {
      tt[, m+1L] <- (-1)^m * (h[, m+2L]*up[, j]*up[, k] -
        l[, m+2L]*down[, j]*down[, k] +
        (j == k)*(h[, m+1L]*up[, j] - l[, m+1L]*down[, j]))
    }
    t0 <- tt[, 1L]; t1 <- tt[, 2L]; t2 <- tt[, 3L]
    cross <- s1[, j]*s0[, k] + s0[, j]*s1[, k]
    prod <- s0[, j]*s0[, k]
    out$Dth2[, col] <- mul(t0 - prod)
    out$Dmuth2[, col] <- mul(t1 - r1*t0 - cross + 2*r1*prod)
    out$Dmu2th2[, col] <- mul(t2 - 2*r1*t1 + (2*r1^2-r2)*t0 -
      s2[, j]*s0[, k] - s0[, j]*s2[, k] - 2*s1[, j]*s1[, k] +
      4*r1*cross + (2*r2-6*r1^2)*prod)
  }
  out
}

# Cauchy interval log likelihoods need not be concave in location. Supply
# positive Fisher working curvature and its derivatives for mgcv's fallback.
.ocat_fisher <- function(mu, theta, wt, level, dist) {
  n <- length(mu); nt <- length(theta)
  info <- dmu <- numeric(n); dth <- matrix(0, n, nt)
  for (j in seq_len(nt + 2L)) {
    d <- .ocat_derivatives(rep(j, n), mu, theta, rep(1, n), min(level, 1L), dist)
    p <- exp(-d$D/2); score <- -d$Dmu/2
    info <- info + 2*p*score^2
    if (level > 0L) {
      dmu <- dmu + 2*p*(score^3 - score*d$Dmu2)
      dth <- dth + 2*p*(-d$Dth/2*score^2 - score*d$Dmuth)
    }
  }
  out <- list(EDmu2 = wt*info)
  if (level > 0L) { out$EDmu3 <- wt*dmu; out$EDmu2th <- wt*dth }
  out
}

# In extreme-value tails one endpoint can dominate to machine precision.
# Differentiate the endpoint log probability directly: expanding p derivatives
# and then subtracting powers of exp(z) would lose all lower-order terms.
.ocat_extreme <- function(out, y, mu, theta, wt, level, link) {
  nt <- length(theta); n <- length(y); wt <- rep_len(wt, n)
  a <- c(-Inf, -1, -1+cumsum(exp(theta)), Inf)
  lo <- a[y]-mu; hi <- a[y+1L]-mu
  if (link == "cloglog") {
    z <- lo; gap <- exp(lo)-exp(hi); endpoint <- y-1L
  } else {
    z <- -hi; gap <- exp(-hi)-exp(-lo); endpoint <- y
  }
  open_end <- if (link == "cloglog") is.infinite(hi) else is.infinite(lo)
  take <- which(is.finite(z) & gap < -50 & (z >= 0 | open_end))
  if (!length(take)) return(out)
  # G^(m)(endpoint-mu) = -exp(z) for cloglog;
  #                           (-1)^(m+1) exp(z) for loglog.
  deriv <- function(m) if (link == "cloglog") -exp(z[take]) else
    (-1)^(m+1L)*exp(z[take])
  weight <- -2*wt[take]
  for (m in seq_len(if (level < 1L) 2L else if (level < 2L) 3L else 4L)) {
    name <- if (m == 1L) "Dmu" else paste0("Dmu",m)
    out[[name]][take] <- weight * (-1)^m * deriv(m)
  }
  out$EDmu2[take] <- out$Dmu2[take]
  if (level < 1L) return(out)
  ends <- vapply(seq_len(nt), function(j)
    exp(theta[j])*(endpoint[take] > j & endpoint[take] < nt+2L), numeric(length(take)))
  ends <- matrix(ends, length(take), nt)
  names <- c("Dth","Dmuth","Dmu2th","Dmu3th")
  for (m in 0:(if (level < 2L) 2L else 3L))
    out[[names[m+1L]]][take, ] <- weight * (-1)^m * deriv(m+1L) * ends
  out$EDmu2th[take, ] <- out$Dmu2th[take, , drop=FALSE]
  if (level < 2L) return(out)
  names <- c("Dth2","Dmuth2","Dmu2th2")
  col <- 0L
  for (j in seq_len(nt)) for (k in j:nt) {
    col <- col+1L
    for (m in 0:2) out[[names[m+1L]]][take,col] <- weight*(-1)^m *
      (deriv(m+2L)*ends[,j]*ends[,k] + (j==k)*deriv(m+1L)*ends[,j])
  }
  out
}
