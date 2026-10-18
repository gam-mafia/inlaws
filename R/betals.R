#' Beta mean and precision family
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Model the interior mean and precision of a beta response with separate
#' additive predictors. Unlike [mgcv::betar()], precision can vary with
#' covariates and smooths.
#' @param link Two links, which must be `"logit"` and `"log"`, respectively.
#' @return A `general.family` object for [mgcv::gam()].
#' @details
#' Supply a list of two formulae: the first includes the response and models
#' `logit(mu)`; the second is one-sided and models `log(phi)`. Both support
#' smooths and offsets. The distribution is `Beta(mu * phi, (1 - mu) * phi)`,
#' with mean `mu` and variance `mu * (1 - mu) / (1 + phi)`. Larger precision
#' means less variability. Use `~ 1` for constant precision.
#'
#' Responses must be finite numeric vectors strictly inside `(0, 1)`.
#' Endpoints are rejected, without truncation or adjustment; see [ordbeta()]
#' for a model including endpoint masses. Case weights are nonnegative
#' likelihood replication weights, not binomial trial totals. Zero weights
#' contribute nothing; at least one weight must be positive. Missing rows
#' follow `gam()`'s `na.action`.
#'
#' Supports `gam()` REML with outer Newton or BFGS optimization, and NCV
#' including grouped neighbourhoods and offsets. Analytic likelihood
#' derivatives through fourth order are supplied. ML, QNCV, `bam()` and
#' `gamm()` are not supported. mgcv can change an ML request to REML for a
#' general family. Separate coefficients are required for the two predictors;
#' additional formulae sharing coefficients are not supported. Precision is
#' not artificially bounded; near-degenerate data or an overly flexible
#' precision model can give boundary estimates and unreliable uncertainty.
#'
#' Fitted values and response predictions have two columns, `mu` and `phi`.
#' Link predictions contain `qlogis(mu)` and `log(phi)`. Response standard
#' errors use the delta method. Predictions condition on fitted coefficients
#' without averaging over coefficient uncertainty or random effects.
#' Constant precision gives the same likelihood as [mgcv::betar()], but
#' REML estimates can differ because precision is an unpenalized coefficient
#' here rather than a family parameter.
#'
#' Deviance maximizes each observation's likelihood over its mean with fitted
#' precision held fixed. The maximizing mean is generally not the response.
#' Deviance residuals use the sign of `y - mu`. Null deviance fits a constant
#' mean predictor retaining mean offsets and fitted precisions; these are
#' conditional mean-fit diagnostics. Log likelihood and AIC include density
#' constants.
#'
#' The distribution callbacks `rd(mu, wt = 1, scale = 1)`,
#' `qf(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE)` and
#' `cdf(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE)` take a
#' two-column response-scale matrix as `mu`. They ignore `wt` and external
#' `scale`. `rd` draws one response per row. Very small shapes can produce
#' simulated endpoints through floating-point rounding. Native
#' `gratia_diagnostics()` supplies response moments and distribution callbacks
#' without requiring gratia. Mathematical details are in
#' `inst/maths/betals.md` in the package source.
#' @seealso [mgcv::betar()], [ordbeta()], [betabinomial()]
#' @export
#' @examples
#' set.seed(42)
#' dat <- data.frame(x = runif(300), z = runif(300))
#' mu <- plogis(sin(2 * pi * dat$x))
#' phi <- exp(2 + dat$z)
#' dat$y <- rbeta(300, mu * phi, (1 - mu) * phi)
#' fit <- mgcv::gam(list(y ~ s(x, k = 6), ~ s(z, k = 5)),
#'   family = betals(), method = "REML", data = dat
#' )
#' head(predict(fit, type = "response"))
betals <- function(link = list("logit", "log")) {
  if (length(link) != 2L ||
    !all(vapply(link, function(x) {
      is.character(x) && length(x) == 1L &&
        !is.na(x)
    }, logical(1))) ||
    link[[1L]] != "logit" || link[[2L]] != "log") {
    stop('betals requires links "logit" and "log"')
  }
  linfo <- lapply(link, stats::make.link)
  linfo[[1L]]$linkinv <- stats::plogis
  linfo[[1L]]$mu.eta <- function(eta) stats::plogis(eta) * stats::plogis(-eta)
  linfo[[2L]]$linkinv <- linfo[[2L]]$mu.eta <- function(eta) exp(eta)
  initialize <- expression({
    y <- family$prepare.response(y)
    n <- rep.int(1, nobs)
    start <- family$init(y, x, E, weights, offset, start)
  })
  postproc <- expression({
    colnames(object$fitted.values) <- c("mu", "phi")
    object <- object$family$finish(object, G$offset)
  })
  family <- structure(
    list(
      family = "betals", link = unlist(link), nlp = 2L,
      linfo = linfo, tri = mgcv::trind.generator(2), ll = .betals_ll,
      ncv = .betals_ncv, initialize = initialize, init = .betals_init,
      prepare.response = .betals_response,
      preinitialize = function(G) {
        if (isTRUE(G$family$qapprox)) stop("betals supports NCV but not QNCV")
        list(y = .betals_response(G$y), offxset = G$offset)
      },
      postproc = postproc, residuals = .betals_residuals,
      finish = .betals_finish,
      sandwich = function(y, X, coef, wt, family, offset = NULL) {
        family$ll(y, X, coef, wt, family, offset, deriv = 1, sandwich = TRUE)$lbb
      },
      rd = function(mu, wt = 1, scale = 1) {
        p <- .betals_parameters(mu)
        stats::rbeta(nrow(mu), p$a, p$b)
      },
      qf = function(p, mu, wt = 1, scale = 1, lower.tail = TRUE, log.p = FALSE) {
        par <- .betals_parameters(mu)
        stats::qbeta(p, par$a, par$b, lower.tail = lower.tail, log.p = log.p)
      },
      cdf = function(q, mu, wt = 1, scale = 1, logp = FALSE, lower.tail = TRUE) {
        par <- .betals_parameters(mu)
        stats::pbeta(q, par$a, par$b, lower.tail = lower.tail, log.p = logp)
      },
      d2link = 1, d3link = 1, d4link = 1, ls = 1, scale = 1,
      no.r.sq = TRUE, available.derivs = 2L, discrete.ok = FALSE
    ),
    class = c("general.family", "extended.family", "family")
  )
  .inlaws_diagnostic_family(family, "betals")
}

.betals_ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                       d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                       D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
  if (is.list(X)) stop("betals does not support discrete model matrices")
  jj <- attr(X, "lpi")
  if (is.null(eta)) {
    eta <- matrix(0, nrow(X), 2L)
    for (k in 1:2) {
      eta[, k] <- drop(X[, jj[[k]], drop = FALSE] %*% coef[jj[[k]]])
      if (length(offset) >= k && !is.null(offset[[k]])) {
        eta[, k] <- eta[, k] + offset[[k]]
      }
    }
  }
  # Direct eta already includes offsets in gamlss.ncv; do not add them twice.
  wt <- rep_len(if (is.null(wt)) 1 else wt, nrow(eta))
  eta[wt == 0, ] <- 0
  u <- stats::plogis(eta[, 1])
  v <- stats::plogis(-eta[, 1])
  phi <- exp(eta[, 2])
  valid <- is.finite(eta[, 1]) & is.finite(eta[, 2]) &
    u > 0 & u < 1 & is.finite(phi) & phi > 0 &
    u * phi > 0 & v * phi > 0
  if (!all(valid)) {
    l0 <- numeric(nrow(eta))
    l0[!valid & wt > 0] <- -Inf
    return(list(l = -Inf, l0 = l0))
  }
  order <- if (deriv == 0) 0L else if (deriv == 1) 2L else if (deriv < 4) 3L else 4L
  y[wt == 0] <- 0.5
  de <- .betals_derivatives(y, eta, order)
  de <- lapply(de, function(x) {
    x <- x * wt
    if (is.matrix(x)) x[wt == 0, ] <- 0 else x[wt == 0] <- 0
    x
  })
  ret <- list()
  if (deriv > 0) {
    packed <- de[-1L]
    l3 <- if (order >= 3) packed[[3]] else 0
    l4 <- if (order >= 4) packed[[4]] else 0
    # NCV supplies eta for row subsets whose X may lack lpi. It only needs
    # observation derivatives, not coefficient-space contractions there.
    if (!is.null(jj)) {
      ret <- mgcv::gamlss.gH(X, jj, packed[[1]], packed[[2]], family$tri$i2,
        l3 = l3, i3 = family$tri$i3, l4 = l4, i4 = family$tri$i4,
        d1b = d1b, d2b = d2b, deriv = deriv - 1L, fh = fh, D = D,
        sandwich = sandwich
      )
    }
    if (ncv) {
      ret$l1 <- packed[[1]]
      ret$l2 <- packed[[2]]
      ret$l3 <- l3
      if (order >= 4) ret$l4 <- l4
    }
  }
  ret$l0 <- de$l0
  ret$l <- sum(ret$l0)
  ret
}

.betals_ncv <- function(X, y, wt, nei, beta, family, llf, H = NULL,
                        Hi = NULL, R = NULL, offset = NULL, dH = NULL,
                        db = NULL, deriv = FALSE, nt = 1) {
  # Pass row indices to preserve correct offset selection for NCV subsets.
  response <- .betals_response(y)
  likelihood <- family$ll
  fit_offset <- lapply(1:2, function(k) {
    if (length(offset) >= k && !is.null(offset[[k]])) {
      rep_len(offset[[k]], length(response))
    } else {
      numeric(length(response))
    }
  })
  family$ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                        d1b = 0, ..., eta = NULL) {
    yy <- response[y]
    off <- lapply(fit_offset, function(x) x[y])
    if (!is.null(eta)) for (k in 1:2) eta[, k] <- eta[, k] + off[[k]]
    likelihood(yy, X, coef, wt, family,
      offset = off, deriv = deriv,
      d1b = d1b, ..., eta = eta
    )
  }
  # Keep the helper's predictors offset-free and add offsets only in the
  # callback. Its gamma adjustment evaluates the full-data eta without
  # adding offsets, whereas its held-out eta normally includes them.
  ret <- utils::getFromNamespace("gamlss.ncv", "mgcv")(
    X, seq_along(response), wt, nei, beta, family, llf, H = H, Hi = Hi,
    R = R, offset = NULL, dH = dH, db = db, deriv = deriv, nt = nt)
  eta <- attr(ret$NCV, "eta.cv")
  if (!is.null(eta)) {
    for (k in 1:2) eta[, k] <- eta[, k] + fit_offset[[k]][nei$d]
    attr(ret$NCV, "eta.cv") <- eta
  }
  ret
}

.betals_init <- function(y, x, E, wt, offset, start) {
  if (is.list(x)) stop("betals requires gam(), not bam()")
  jj <- attr(x, "lpi")
  if (length(jj) != 2L) stop("betals requires two formulae")
  if (anyDuplicated(unlist(jj))) stop("betals does not support shared predictor coefficients")
  if (any(!is.finite(wt) | wt < 0) || !any(wt > 0)) {
    stop("betals requires finite non-negative weights with positive total")
  }
  if (!is.null(start)) {
    return(start)
  }
  keep <- wt > 0
  off <- lapply(1:2, function(k) {
    if (length(offset) >= k && !is.null(offset[[k]])) {
      rep_len(offset[[k]], length(y))
    } else {
      numeric(length(y))
    }
  })
  target <- stats::qlogis(y)
  solve_start <- function(k, target) {
    xx <- x[keep, jj[[k]], drop = FALSE] * sqrt(wt[keep])
    ee <- E[, jj[[k]], drop = FALSE]
    if (length(ee) && sum(ee^2) > 0) {
      ee <- ee * (0.01 * sqrt(sum(xx^2) / sum(ee^2)))
    }
    ans <- qr.coef(
      qr(rbind(xx, ee)),
      c((target[keep] - off[[k]][keep]) * sqrt(wt[keep]), rep(0, nrow(ee)))
    )
    ans[!is.finite(ans)] <- 0
    ans
  }
  start <- numeric(ncol(x))
  start[jj[[1]]] <- solve_start(1, target)
  loc <- drop(x[, jj[[1]], drop = FALSE] %*% start[jj[[1]]]) + off[[1]]
  m <- stats::plogis(loc[keep])
  variance <- stats::weighted.mean((y[keep] - m)^2, wt[keep])
  phi <- max(0.1, min(1e4, stats::weighted.mean(m * (1 - m), wt[keep]) / max(variance, 1e-8) - 1))
  start[jj[[2]]] <- solve_start(2, rep(log(phi), length(y)))
  start
}


.betals_response <- function(y) {
  if (!is.numeric(y) || !is.null(dim(y)) ||
    any(!is.finite(y) | y <= 0 | y >= 1)) {
    stop("betals requires a finite numeric response vector strictly between zero and one")
  }
  y
}

.betals_parameters <- function(mu) {
  if (!is.matrix(mu) || !is.numeric(mu) || ncol(mu) != 2L ||
    any(!is.finite(mu)) || any(mu[, 1] <= 0 | mu[, 1] >= 1) ||
    any(mu[, 2] <= 0)) {
    stop("mu must be a two-column matrix containing interior mean and positive precision")
  }
  list(a = mu[, 1] * mu[, 2], b = (1 - mu[, 1]) * mu[, 2])
}

# H(t) is the Stirling correction to log Gamma(exp(t)). Keeping the
# large-shape terms separate prevents cancellation in precision derivatives.
.betals_gamma_correction <- function(t, order) {
  a <- exp(t)
  out <- matrix(0, length(t), order + 1L)
  large <- a >= 16
  if (any(large)) {
    powers <- c(1, 3, 5, 7, 9, 11, 13)
    coefficients <- c(
      1 / 12, -1 / 360, 1 / 1260, -1 / 1680, 1 / 1188,
      -691 / 360360, 1 / 156
    )
    for (k in seq_along(powers)) {
      term <- coefficients[k] * exp(-powers[k] * t[large])
      for (q in 0:order) out[large, q + 1L] <- out[large, q + 1L] + (-powers[k])^q * term
    }
  }
  small <- !large
  if (any(small)) {
    x <- a[small]
    z <- t[small]
    out[small, 1] <- lgamma(x + 1) - z - (x - .5) * z + x - log(2 * pi) / 2
    stirling <- list(1, c(1, 1), c(1, 3, 1), c(1, 7, 6, 1))
    for (q in seq_len(order)) {
      g <- numeric(length(x))
      for (k in seq_len(q)) g <- g + stirling[[q]][k] * x^k * psigamma(x + 1, k - 1L)
      out[small, q + 1L] <- g - (q == 1L) / 2 - x * (z + q - 1)
    }
  }
  out
}

# Packed predictor derivatives: eta1^r, eta1^(r-1)*eta2, ..., eta2^r.
.betals_derivatives <- function(y, eta, order = 4L) {
  u <- stats::plogis(eta[, 1])
  v <- stats::plogis(-eta[, 1])
  phi <- exp(eta[, 2])
  out <- list(l0 = stats::dbeta(y, u * phi, v * phi, log = TRUE))
  if (!order) {
    return(out)
  }
  lu <- -pmax(-eta[, 1], 0) - log1p(exp(-abs(eta[, 1])))
  lv <- -pmax(eta[, 1], 0) - log1p(exp(-abs(eta[, 1])))
  du <- list(v, -u * v, -u * v * (v - u), -u * v * (1 - 6 * u * v))
  dv <- du
  dv[[1]] <- -u
  dm <- list(u * v, u * v * (v - u), u * v * (1 - 6 * u * v), u * v * (v - u) * (1 - 12 * u * v))
  # B = -KL(Bernoulli(mu) || Bernoulli(y)). Use a series near y=mu
  # to preserve the quadratic term when the linear log terms cancel.
  B <- u * (log(y) - lu) + v * (log1p(-y) - lv)
  delta <- y - u
  close <- abs(delta) < .01 * pmin(u, v)
  B[close] <- 0
  for (k in 2:10) {
    B[close] <- B[close] -
      (u[close] * (-delta[close] / u[close])^k + v[close] * (delta[close] / v[close])^k) / k
  }
  bd <- c(list(B), lapply(1:4, function(p) {
    dm[[p]] * (stats::qlogis(y) - eta[, 1]) - if (p > 1) (p - 1) * dm[[p - 1]] else 0
  }))
  ha <- .betals_gamma_correction(lu + eta[, 2], order)
  hb <- .betals_gamma_correction(lv + eta[, 2], order)
  hc <- .betals_gamma_correction(eta[, 2], order)
  bell <- function(d) {
    list(
      list(1), list(d[[1]]),
      list(d[[2]], d[[1]]^2),
      list(d[[3]], 3 * d[[1]] * d[[2]], d[[1]]^3),
      list(
        d[[4]], 4 * d[[1]] * d[[3]] + 3 * d[[2]]^2,
        6 * d[[1]]^2 * d[[2]], d[[1]]^4
      )
    )
  }
  bu <- bell(du)
  bv <- bell(dv)
  for (r in seq_len(order)) {
    dd <- matrix(0, length(y), r + 1L)
    for (q in 0:r) {
      p <- r - q
      value <- phi * bd[[p + 1L]]
      if (p == 0L) {
        value <- value + (q == 1L) / 2 + hc[, q + 1L] - ha[, q + 1L] - hb[, q + 1L]
      } else {
        for (k in seq_len(p)) {
          value <- value -
            bu[[p + 1L]][[k]] * ha[, q + k + 1L] - bv[[p + 1L]][[k]] * hb[, q + k + 1L]
        }
        if (q == 0L) value <- value + (du[[p]] + dv[[p]]) / 2
      }
      dd[, q + 1L] <- value
    }
    out[[paste0("l", r)]] <- dd
  }
  out
}

.betals_deviance <- function(y, mu, wt) {
  ans <- numeric(length(y))
  keep <- which(wt > 0)
  for (i in keep) {
    ends <- sort(c(stats::qlogis(y[i]), 0))
    logphi <- log(mu[i, 2])
    score <- function(a) .betals_derivatives(y[i], matrix(c(a, logphi), 1), 1L)$l1[1, 1]
    if (diff(ends) < 1e-12) {
      location <- mean(ends)
    } else {
      ss <- vapply(ends, score, numeric(1))
      location <- if (ss[1] * ss[2] < 0) {
        stats::uniroot(score, ends, tol = 1e-11)$root
      } else {
        ends[which.min(abs(ss))]
      }
    }
    sat <- .betals_derivatives(y[i], matrix(c(location, logphi), 1), 0L)$l0
    fit <- stats::dbeta(y[i], mu[i, 1] * mu[i, 2], (1 - mu[i, 1]) * mu[i, 2], log = TRUE)
    ans[i] <- max(0, 2 * wt[i] * (sat - fit))
  }
  ans
}

.betals_residuals <- function(object, type = c("deviance", "pearson", "response", "scaled.pearson"), ...) {
  type <- match.arg(type)
  mu <- object$fitted.values
  r <- object$y - mu[, 1]
  if (type == "deviance") {
    return(sign(r) * sqrt(object$beta.deviance))
  }
  if (type != "response") r <- r * sqrt(object$prior.weights * (1 + mu[, 2]) / (mu[, 1] * (1 - mu[, 1])))
  r
}

.betals_finish <- function(object, offset) {
  y <- object$y
  w <- object$prior.weights
  mu <- object$fitted.values
  object$beta.deviance <- .betals_deviance(y, mu, w)
  object$deviance <- sum(object$beta.deviance)
  keep <- w > 0
  off <- if (length(offset) && !is.null(offset[[1]])) rep_len(offset[[1]], length(y)) else numeric(length(y))
  objective <- function(a) {
    eta <- cbind(a + off[keep], log(mu[keep, 2]))
    -sum(w[keep] * .betals_derivatives(y[keep], eta, 0L)$l0)
  }
  initial <- stats::weighted.mean(stats::qlogis(mu[, 1]) - off, w)
  opt <- stats::optim(initial, objective, method = "BFGS", control = list(reltol = 1e-10))
  if (opt$convergence != 0L) warning("betals conditional null fit did not converge")
  fitted.loglik <- sum(w[keep] * stats::dbeta(y[keep], mu[keep, 1] * mu[keep, 2],
    (1 - mu[keep, 1]) * mu[keep, 2],
    log = TRUE
  ))
  object$null.deviance <- max(0, object$deviance + 2 * (fitted.loglik + opt$value))
  object
}
