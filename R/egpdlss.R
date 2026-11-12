#' Extended generalized Pareto models
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A three-predictor family for strictly positive responses with location fixed
#' at zero. This is model 1 of Naveau et al. (2016), as implemented in
#' `evgam`: the CDF is \eqn{G(y)^\kappa}, where \eqn{G} is the generalized
#' Pareto CDF from [gpdlss()]. The density is
#' \eqn{\kappa g(y) G(y)^{\kappa-1}}. At \eqn{\kappa=1} both reduce to
#' the generalized Pareto distribution. For negative tail shape \eqn{\xi},
#' the upper endpoint is \eqn{-\sigma/\xi}.
#'
#' @param link Three links in scale, tail shape, bulk shape order. Scale
#'   and bulk shape use `"log"`. Tail shape uses `"logshift"` (the default,
#'   constraining \eqn{\xi>-0.5}) or `"identity"` (unrestricted).
#'
#' @details Supply three formulae in a list: response and scale first, then
#'   tail shape, then bulk shape. Response predictions return the natural
#'   parameters in order. The fitted matrix names these columns `scale`,
#'   `shape`, and `kappa`.
#'   Use strictly positive observations at positive weight because the density
#'   at zero can vanish or diverge depending on bulk shape. The probability
#'   callbacks include zero as the lower endpoint. Prior weights multiply
#'   log likelihoods.
#'
#'   With mgcv 1.9-5, [mgcv::gam()] supports REML with outer Newton, outer
#'   BFGS, or extended Fellner–Schall (EFS). For NCV and QNCV, mgcv selects
#'   outer BFGS internally regardless of the requested optimizer.
#'   `mgcv::scasm()` uses EFS for REML smoothing selection, including with
#'   constrained smooths. It exposes only this smoothing selection route.
#'   See the family compatibility vignette for details.
#'
#'   The family supplies `rd(mu, wt = NULL, scale = NULL)`,
#'   `qf(p, mu, wt = NULL, scale = NULL)` and
#'   `cdf(q, mu, wt = NULL, scale = NULL, log_p = FALSE)`.
#'   Supply `mu` as a numeric three-column matrix of natural parameters.
#'   `p` and `q` can be scalars or have one value per row. Simulation returns
#'   one draw per row. The compatibility arguments `wt` and `scale` do not
#'   affect these distributions.
#'   Set `lower.tail = FALSE` for survival probabilities or upper-tail
#'   quantiles. Quantiles also accept log probabilities via `log.p = TRUE`.
#'
#'   Response residuals require \eqn{\xi<1} and Pearson residuals require
#'   \eqn{\xi<1/2}. Residuals requiring unavailable moments are `NA` with
#'   a warning. Use PIT or quantile diagnostics when moments are unavailable.
#'   Deviance and deviance residuals are unavailable.
#'
#' @references
#' Naveau, P., Huser, R., Ribereau, P., and Hannart, A. (2016).
#' Modeling jointly low, moderate, and heavy rainfall intensities without a
#' threshold selection. *Water Resources Research*, 52, 2753–2769.
#' \doi{10.1002/2015WR018552}
#' @return An mgcv `general.family` object.
#' @template ncv-neighbourhoods
#' @export
#' @examples
#' set.seed(12)
#' d <- data.frame(x = runif(400))
#' fam <- egpdlss()
#' d$y <- fam$rd(cbind(exp(sin(2 * pi * d$x)), 0.2, 1.5))
#' fit <- mgcv::gam(list(y ~ s(x, k = 6), ~ 1, ~ 1), family = fam,
#'                  data = d, method = "REML")
#' head(predict(fit, type = "response"))
#' head(fit$family$qf(0.95, fitted(fit)))
egpdlss <- function(link = list("log", "logshift", "log")) {
  if (!(is.list(link) || is.character(link)) || length(link) != 3L ||
      any(lengths(link) != 1L) ||
      !identical(link[[1L]], "log") || !identical(link[[3L]], "log") ||
      !is.character(link[[2L]]) || is.na(link[[2L]]) ||
      !link[[2L]] %in% c("logshift", "identity")) {
    stop('egpdlss requires links "log", "logshift" or "identity", and "log"')
  }
  loglink <- .gpd_log_link()
  loglink$linkinv <- loglink$mu.eta <- exp
  family <- structure(list(
    family = "egpdlss", link = unlist(link), nlp = 3L,
    linfo = list(loglink, .gpd_shape_link(link[[2L]]), loglink),
    tri = mgcv::trind.generator(3), ll = .egpd_ll, ncv = .inlaws_ncv3,
    preinitialize = .egpd_preinitialize,
    initialize.coef = .gpd_start, check.nei = .inlaws_check_ncv_nei,
    initialize = expression({
      n <- rep(1, nobs)
      start <- family$initialize.coef(y, x, weights, offset, start, family)
    }),
    presetup = expression({
      # scasm has no nei formal and ignores method/optimizer in its dots.
      if (exists("nei", inherits = FALSE)) family$check.nei(nei, attr(mf, "na.action"))
    }),
    postproc = expression({
      object$deviance <- object$null.deviance <- NA_real_
      colnames(object$fitted.values) <- c("scale", "shape", "kappa")
      if (isTRUE(object$scale.known)) {
        # scasm.newton can stop with AIC evaluated at the last trial step.
        # Evaluate it at the returned coefficients used by logLik.gam().
        ll <- object$family$ll(object$y, NULL, NULL, object$prior.weights,
                               object$family, eta = object$linear.predictors)
        object$aic <- 2 * (sum(object$edf) - ll$l)
      }
    }),
    predict = .gpd_predict, residuals = .egpd_residuals,
    rd = .egpd_rd, qf = .egpd_qf, cdf = .egpd_cdf,
    d2link = 1, d3link = 1, d4link = 1, ls = 1, scale = 1,
    available.derivs = 2, discrete.ok = FALSE
  ), class = c("general.family", "extended.family", "family"))
  .inlaws_diagnostic_family(family, "egpdlss")
}

.egpd_preinitialize <- function(G) {
  # estimate.gam() resolves the actual optimizer before calling this hook.
  # scasm.fit() has no optimizer formal and handles linked/fixed penalties.
  opt <- get0("optimizer", envir = parent.frame(), inherits = FALSE)
  if (length(opt) && opt[1L] == "efs" && length(G$sp) &&
      ((!is.null(G$L) && !isTRUE(all.equal(G$L, diag(length(G$sp)),
                                          check.attributes = FALSE))) || any(G$lsp0 != 0)))
    stop("egpdlss EFS requires all smoothing parameters free and unlinked, or all fixed; use an outer optimizer")
  # The general-family EFS dispatch uses this spelling in mgcv 1.9-3.
  list(offxset = G$offset)
}

.egpd_ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                     d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                     D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
  if (is.null(eta)) eta <- .gpd_eta(X, coef, offset)
  if (is.null(wt)) wt <- rep(1, length(y))
  active <- wt > 0
  eta[!active, ] <- 0
  eta[!active, 2L] <- family$linfo[[2L]]$linkfun(0)
  y[!active] <- 1
  rho <- eta[, 1L]
  xi <- family$linfo[[2L]]$linkinv(eta[, 2L])
  sigma <- exp(rho)
  kappa <- exp(eta[, 3L])
  t <- y / sigma
  end <- -sigma / xi
  w <- 1 + xi * t
  # The scaled ratio can round onto the bound on an interior observation.
  near_end <- which(xi < 0 & w <= 0 & is.finite(end) & y < end)
  w[near_end] <- (end[near_end] - y[near_end]) / end[near_end]
  valid <- y > 0 & is.finite(t) & is.finite(sigma) & sigma > 0 &
    is.finite(xi) & w > 0 & (xi >= 0 | y < end) & is.finite(kappa) & kappa > 0
  if (family$link[2L] == "logshift") valid <- valid & xi > -0.5
  if (!all(valid)) {
    l0 <- numeric(length(y))
    l0[!valid] <- -Inf
    return(list(l = -Inf, l0 = l0))
  }
  order <- if (deriv == 0L) 0L else if (deriv == 1L) 2L else if (deriv < 4L) 3L else 4L
  ds <- .egpd_derivatives(y, rho, xi, eta[, 3L], order, w)
  if (!deriv) return(list(l = sum(wt * ds$l0), l0 = wt * ds$l0))
  .gpd_gH(ds, X, wt, xi, eta, family, order, deriv, d1b, d2b, fh, D, ncv, sandwich)
}

# Packed derivatives of l = l_GPD + log(kappa) + (kappa - 1) A(H),
# where A = log(1 - exp(-H)) and H is the GPD cumulative hazard.
.egpd_derivatives <- function(y, rho, xi, logk, order, w) {
  gpd <- .gpd_derivatives(y, rho, xi, order, hazard = TRUE, w = w)
  h <- gpd$H[[1L]]
  a <- .gpd_logcdf(h)
  tiny <- h < .Machine$double.xmin
  a[tiny] <- log(y[tiny]) - rho[tiny]
  # For large H, log(-A) is -H to working precision. This retains
  # representable products even when A alone underflows.
  log_a <- ifelse(h > 37, -h, log(-a))
  delta <- expm1(logk)
  ka <- -exp(logk + log_a)
  l0 <- gpd$l0 + logk - sign(delta) * exp(log(abs(delta)) + log_a)
  out <- list(l0 = l0)
  if (!order) return(out)

  # Apply the chain rule to A(H) using derivatives of H divided by H.
  # The compensating factors H^j A^(j) stay finite as H tends to zero.
  # Incorporate the kappa factor before exponentiating the first factor.
  b <- h / -expm1(-h)
  b[tiny] <- 1
  log_r <- log(b) - h
  r <- exp(log_r)
  # In this range even the largest finite derivative products vanish.
  far_tail <- h > 1e6
  b[far_tail] <- 1
  f2 <- -b
  f3 <- b * (b + r)
  f4 <- -b * (b^2 + 4 * b * r + r^2)
  hder <- function(...) {
    ind <- c(...)
    d <- gpd$H[[length(ind) + 1L]][, 1L + sum(ind == 2L)] / h
    # In the lower tail H ~ y / sigma; its shape derivatives vanish.
    d[tiny] <- if (any(ind == 2L)) 0 else (-1)^length(ind)
    d[far_tail] <- 0
    d
  }
  chain <- function(ind, bulk) {
    log_weight <- (if (bulk) logk else log(abs(delta))) + log_r
    direction <- if (bulk) 1 else sign(delta)
    # Multiply after combining magnitudes to retain very small derivatives.
    weighted <- function(v) direction * sign(v) * exp(log_weight + log(abs(v)))
    i <- ind[1L]; j <- ind[2L]; k <- ind[3L]; l <- ind[4L]
    h1 <- hder(i)
    if (length(ind) == 1L) return(weighted(h1))
    h2 <- hder(j)
    v <- hder(i, j) + f2 * h1 * h2
    if (length(ind) == 2L) return(weighted(v))
    h3 <- hder(k)
    v <- hder(i, j, k) + f2 * (hder(i, j) * h3 + hder(i, k) * h2 +
      hder(j, k) * h1) + f3 * h1 * h2 * h3
    if (length(ind) == 3L) return(weighted(v))
    h4 <- hder(l)
    weighted(hder(i, j, k, l) +
      f2 * (hder(i, j, k) * h4 + hder(i, j, l) * h3 +
        hder(i, k, l) * h2 + hder(j, k, l) * h1 +
        hder(i, j) * hder(k, l) + hder(i, k) * hder(j, l) + hder(i, l) * hder(j, k)) +
      f3 * (hder(i, j) * h3 * h4 + hder(i, k) * h2 * h4 + hder(i, l) * h2 * h3 +
        hder(j, k) * h1 * h4 + hder(j, l) * h1 * h3 + hder(k, l) * h1 * h2) +
      f4 * h1 * h2 * h3 * h4)
  }
  for (d in seq_len(order)) {
    ld <- matrix(0, length(y), choose(d + 2L, 2L))
    col <- 0L
    for (nr in d:0) for (nx in (d - nr):0) {
      nk <- d - nr - nx
      col <- col + 1L
      if (nk == d) {
        ld[, col] <- ka + as.numeric(d == 1L)
      } else {
        ld[, col] <- chain(rep(1:2, c(nr, nx)), bulk = nk > 0L)
        if (nk == 0L) ld[, col] <- ld[, col] + gpd[[paste0("l", d)]][, nx + 1L]
      }
    }
    out[[d + 1L]] <- ld
    names(out)[d + 1L] <- paste0("l", d)
  }
  out
}

.egpd_args <- function(x, mu) {
  if (!is.matrix(mu) || !is.numeric(mu) || ncol(mu) != 3L)
    stop("mu must be a numeric matrix with scale, shape and kappa columns")
  a <- .gpd_args(x, mu[, 1:2, drop = FALSE])
  a$kappa <- mu[, 3L]
  a$missing <- a$missing | is.na(a$kappa)
  a$bad <- !a$missing & (a$bad | !is.finite(a$kappa) | a$kappa <= 0)
  a
}

.egpd_cdf <- function(q, mu, wt = NULL, scale = NULL, log_p = FALSE, lower.tail = TRUE) {
  a <- .egpd_args(q, mu)
  if (any(a$bad)) warning("NaNs produced by invalid egpdlss parameters")
  q <- a$x
  good <- !a$missing & !a$bad
  end <- q > 0 & (q == Inf | (a$xi < 0 & q >= -a$sigma / a$xi))
  inside <- good & q > 0 & !end
  log_neg_lp <- rep(NA_real_, length(q))
  log_neg_lp[good & q <= 0] <- Inf
  log_neg_lp[good & end] <- -Inf
  h <- .gpd_phazard(q[inside], a$sigma[inside], a$xi[inside])
  lp <- .gpd_logcdf(h)
  tiny <- h < .Machine$double.xmin
  lp[tiny] <- log(q[inside][tiny]) - log(a$sigma[inside][tiny])
  log_neg_lp[inside] <- log(a$kappa[inside]) + ifelse(h > 37, -h, log(-lp))
  lp <- if (lower.tail) -exp(log_neg_lp) else .egpd_log_survival(log_neg_lp)
  lp[a$bad] <- NaN
  if (log_p) lp else exp(lp)
}

# If log(F) = -exp(z), return log(1 - F), including very small survival.
.egpd_log_survival <- function(z) {
  ifelse(z < -37, z, .gpd_logcdf(exp(z)))
}

.egpd_qf <- function(p, mu, wt = NULL, scale = NULL, log.p = FALSE, lower.tail = TRUE) {
  a <- .egpd_args(p, mu)
  p <- a$x
  outside <- if (log.p) p > 0 else p < 0 | p > 1
  bad <- a$bad | (!a$missing & outside)
  if (any(bad)) warning("NaNs produced by invalid egpdlss probabilities or parameters")
  good <- !a$missing & !bad
  p <- p[good]
  # z is log(-log(G)), for G = p_lower^(1/kappa).
  z <- if (lower.tail) {
    if (log.p) log(-p) else log(-log(p))
  } else {
    if (log.p) ifelse(p < -37, p, log(-.gpd_logcdf(-p))) else log(-log1p(-p))
  }
  z <- z - log(a$kappa[good])
  h <- -.egpd_log_survival(z)
  logh <- ifelse(z > log(37), -exp(z), log(h))
  q <- rep(NA_real_, length(good))
  q[bad] <- NaN
  sigma <- a$sigma[good]
  xi <- a$xi[good]
  value <- numeric(length(z))
  inside <- is.finite(z)
  value[inside] <- .gpd_qhazard(h[inside], sigma[inside], xi[inside], logh[inside])
  end <- z == -Inf
  value[end] <- ifelse(xi[end] < 0, -sigma[end] / xi[end], Inf)
  q[good] <- value
  q
}

.egpd_rd <- function(mu, wt = NULL, scale = NULL) {
  .egpd_qf(stats::runif(nrow(mu)), mu, wt, scale)
}

# Let E = -log(1 - U), where U has CDF u^kappa. Then
# E[exp(t*E)] = B(1-t, kappa)/B(1, kappa).
# For xi near zero, expand its log so the first divided difference and
# the variance retain their limits even at xi = 0.
.egpd_moments <- function(mu) {
  sigma <- mu[, 1L]
  xi <- mu[, 2L]
  kappa <- mu[, 3L]
  ok <- is.finite(sigma) & sigma > 0 & is.finite(xi) & is.finite(kappa) & kappa > 0
  mean <- variance <- rep(NA_real_, nrow(mu))
  mean_ok <- ok & xi < 1
  var_ok <- ok & xi < 0.5
  variance[ok & !var_ok] <- Inf
  gpd <- ok & kappa == 1
  i <- mean_ok & gpd
  mean[i] <- sigma[i] / (1 - xi[i])
  i <- var_ok & gpd
  variance[i] <- (sigma[i] / (1 - xi[i]))^2 / (1 - 2 * xi[i])

  i <- which(mean_ok & !gpd & abs(xi) > 1e-3)
  l1 <- .egpd_log_hazard_mgf(xi[i], kappa[i])
  log_diff <- ifelse(l1 >= 0, .egpd_log_expm1(pmax(0, l1)), log(-expm1(pmin(0, l1))))
  mean[i] <- exp(log(sigma[i]) + log_diff - log(abs(xi[i])))
  j <- var_ok[i]
  i <- i[j]
  l1 <- l1[j]
  l2 <- .egpd_log_hazard_mgf(2 * xi[i], kappa[i])
  variance[i] <- exp(2 * (log(sigma[i]) + l1 - log(abs(xi[i]))) +
                      .egpd_log_expm1(pmax(0, l2 - 2*l1)))

  i <- which(mean_ok & !gpd & abs(xi) <= 1e-3)
  if (length(i)) {
    x <- xi[i]
    k <- kappa[i]
    coefficient <- function(j) {
      diff <- psigamma(1, j - 1L) - psigamma(1 + k, j - 1L)
      small <- k < 1e-5
      diff[small] <- -k[small] * (psigamma(1, j) + k[small] / 2 * psigamma(1, j + 1L))
      (-1)^j * diff / factorial(j)
    }
    mean_scaled <- coefficient(1L)
    var_scaled <- numeric(length(i))
    for (j in 2:8) {
      cj <- coefficient(j)
      mean_scaled <- mean_scaled + cj * x^(j - 1L)
      var_scaled <- var_scaled + (2^j - 2) * cj * x^(j - 2L)
    }
    exprel <- function(z) {
      out <- rep(1, length(z))
      nz <- z != 0
      out[nz] <- expm1(z[nz]) / z[nz]
      out
    }
    l1 <- x * mean_scaled
    mean[i] <- sigma[i] * mean_scaled * exprel(l1)
    variance[i] <- exp(2 * (log(sigma[i]) + l1) + log(var_scaled) +
                        log(exprel(x^2 * var_scaled)))
  }
  list(mean = mean, variance = variance, mean_available = all(mean_ok),
       variance_available = all(var_ok))
}

.egpd_log_hazard_mgf <- function(t, kappa) {
  a <- 1 - t
  out <- lbeta(a, kappa) + log(kappa)
  # Expand log Gamma(1 + kappa) - log Gamma(a + kappa) + log Gamma(a)
  # when adding log(kappa) to log Beta would lose the result.
  small <- kappa < 1e-5 * pmin(1, a)
  b <- a[small]
  k <- kappa[small]
  out[small] <- k * (digamma(1) - digamma(b) +
    k/2 * (trigamma(1) - trigamma(b) + k/3 * (psigamma(1, 2) - psigamma(b, 2))))
  out
}

.egpd_log_expm1 <- function(x) {
  ifelse(x > 37, x, log(expm1(pmin(x, 37))))
}

.egpd_residuals <- function(object, type = c("deviance", "pearson", "response")) {
  type <- match.arg(type)
  if (type == "deviance")
    stop("egpdlss deviance residuals are unavailable; use PIT, response or Pearson residuals")
  mu <- object$fitted.values
  wt <- object$prior.weights
  if (is.null(wt)) wt <- rep(1, nrow(mu))
  moments <- .egpd_moments(mu)
  out <- object$y - moments$mean
  available <- mu[, 2L] < 1
  if (type == "pearson") {
    out <- sqrt(wt) * out / sqrt(moments$variance)
    available <- mu[, 2L] < 0.5
    out[wt == 0] <- 0
  }
  if (any(!available & wt > 0))
    warning("egpdlss residuals undefined because required moments do not exist")
  out[!available & (type == "response" | wt > 0)] <- NA_real_
  out
}
