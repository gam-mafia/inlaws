#' Zero-inflated and hurdle negative binomial families
#'
#' Three additive predictors model the underlying NB mean, count-component
#' probability, and NB size, in that order. The NB variance is
#' `mu + mu^2/theta` before mixing or truncation.
#' @param link Three links, currently only `list("log", "logit", "log")`.
#' @return A `general.family` for [mgcv::gam()].
#' @details Supply three formulae, with the response in the first. Use `~1`
#'   for constant parameters. Offsets and shared coefficients are supported.
#'   `zinb()` is a true mixture: the NB component can produce zeros. Its zero
#'   probability is `1-p+p*q0`, where `q0` is the NB zero probability.
#'   `zanb()` is a hurdle model: zero probability is `1-p`, and positive counts
#'   follow the zero-truncated NB. Here `mu` always refers to the untruncated NB.
#'
#'   Fitted values and response predictions are matrices with columns ordered `mu`,
#'   `p`, `theta`, including parameter-wise delta-method standard errors.
#'   Fitted-value columns are named; mgcv may omit response-prediction names.
#'   Overall response expectations are `p*mu` for ZINB and `p*mu/(1-q0)` for
#'   ZANB. Response residuals use these expectations, not the first column.
#'
#'   REML supports Newton, BFGS and EFS; NCV and QNCV use BFGS. EFS does not
#'   optimize NCV. `bam()` and `gamm()` are unsupported. In mgcv 1.9-3,
#'   fixed smoothing parameters with a later predictor having no smooths can
#'   trigger an upstream setup error: use `G <- gam(..., fit=FALSE)` followed
#'   by `gam(G=G, sp=...)`, or upgrade to mgcv 1.9-4. Weights multiply the
#'   likelihood and do not change individual simulated response distributions.
#'   The family supplies CDF, quantile and random generation callbacks.
#'
#'   Data must be finite nonnegative integers with finite nonnegative weights.
#'   All-zero positive-weight responses are rejected. No zeros, positive counts
#'   all equal to one, and the Poisson limit can cause boundary estimates.
#'   Flexible mean, size and mixture predictors can be weakly identified,
#'   especially for ZINB. Check convergence and uncertainty; size is not capped.
#'   Constant size is estimated as a coefficient, unlike [mgcv::nb()].
#'
#'   Deviance holds fitted size fixed and uses the observation-wise saturated
#'   likelihood. Null deviance fits constant mean and probability predictors,
#'   retaining offsets and fitted sizes; it is conditional on the size fit.
#'   Mathematical details are in `inst/maths/zinb-zanb.md`.
#' @export
#' @examples
#' set.seed(29)
#' d <- data.frame(x = runif(600), z = runif(600))
#' d$y <- rbinom(600, 1, plogis(1 + d$x)) *
#'   rnbinom(600, mu = exp(1 + sin(6 * d$x)), size = exp(1 + d$z))
#' b <- mgcv::gam(list(y ~ s(x, k = 6), ~x, ~z),
#'   data = d,
#'   family = zinb(), method = "REML"
#' )
#' head(predict(b, type = "response"))
#' head(b$family$rd(fitted(b)))
zinb <- function(link = list("log", "logit", "log")) .zinb_family(link, FALSE)

#' @rdname zinb
#' @export
zanb <- function(link = list("log", "logit", "log")) .zinb_family(link, TRUE)

.zinb_family <- function(link, hurdle) {
  if (length(link) != 3L || !identical(unname(unlist(link)), c("log", "logit", "log"))) {
    stop('links must be list("log", "logit", "log")')
  }
  linfo <- lapply(link, stats::make.link)
  for (j in c(1L, 3L)) linfo[[j]]$linkinv <- linfo[[j]]$mu.eta <- function(x) exp(x)
  linfo[[2L]]$linkinv <- stats::plogis
  linfo[[2L]]$mu.eta <- function(x) stats::plogis(x) * stats::plogis(-x)
  name <- if (hurdle) "zanb" else "zinb"
  structure(
    list(
      family = name, link = unlist(link), nlp = 3L, hurdle = hurdle,
      linfo = linfo, tri = mgcv::trind.generator(3), ll = .zinb_ll,
      preinitialize = function(G) list(offxset = G$offset),
      initialize = expression({
        n <- rep.int(1, nobs)
        start <- family$init(y, x, E, weights, offset, start, family)
      }),
      init = .zinb_initialize,
      postproc = expression({
        colnames(object$fitted.values) <- c("mu", "p", "theta")
        object$null.deviance <- object$family$null.deviance(object, G$offset)
      }),
      residuals = .zinb_residuals, ncv = .zinb_ncv,
      null.deviance = .zinb_null,
      sandwich = function(y, X, coef, wt, family, offset = NULL) {
        family$ll(y, X, coef, wt, family, offset = offset, deriv = 1, sandwich = TRUE)$lbb
      },
      rd = function(mu, wt = NULL, scale = NULL) .zinb_quantile(stats::runif(nrow(mu)), mu, hurdle),
      qf = function(p, mu, wt = NULL, scale = NULL) .zinb_quantile(p, mu, hurdle),
      cdf = function(q, mu, wt = NULL, scale = NULL, logp = FALSE) .zinb_cdf(q, mu, hurdle, logp),
      d2link = 1, d3link = 1, d4link = 1, ls = 1, available.derivs = 2L, discrete.ok = FALSE
    ),
    class = c("general.family", "extended.family", "family")
  )
}

# Multivariate Taylor coefficients (derivative / multi-index factorial).
# Polynomial arithmetic is analytic differentiation, not finite differences.
.zinb_jet_spec <- local({
  cache <- vector("list", 5L)
  function(order) {
    if (!is.null(cache[[order + 1L]])) {
      return(cache[[order + 1L]])
    }
    a <- as.matrix(expand.grid(0:order, 0:order, 0:order))
    a <- a[rowSums(a) <= order, , drop = FALSE]
    a <- a[order(rowSums(a), a[, 3], a[, 2]), , drop = FALSE]
    key <- apply(a, 1, paste, collapse = ",")
    pairs <- lapply(seq_len(nrow(a)), function(k) {
      i <- which(apply(sweep(a, 2, a[k, ], "<="), 1, all))
      j <- match(apply(matrix(a[k, ], length(i), 3, byrow = TRUE) - a[i, , drop = FALSE], 1, paste, collapse = ","), key)
      cbind(i, j)
    })
    cache[[order + 1L]] <<- list(a = a, pairs = pairs, fact = apply(factorial(a), 1, prod), order = order)
    cache[[order + 1L]]
  }
})
.zinb_jmul <- function(x, y, s) {
  out <- x * 0
  for (k in seq_len(ncol(x))) {
    ij <- s$pairs[[k]]
    out[, k] <- rowSums(x[, ij[, 1], drop = FALSE] * y[, ij[, 2], drop = FALSE])
  }
  out
}
.zinb_jcompose <- function(x, derivatives, s) {
  dx <- x
  dx[, 1] <- 0
  out <- x * 0
  out[, 1] <- derivatives[[1]]
  power <- x * 0
  power[, 1] <- 1
  if (s$order) {
    for (k in seq_len(s$order)) {
      power <- .zinb_jmul(power, dx, s)
      out <- out + power * (derivatives[[k + 1L]] / factorial(k))
    }
  }
  out
}
.zinb_jexp <- function(x, s) .zinb_jcompose(x, rep(list(exp(x[, 1])), s$order + 1L), s)
.zinb_jlog <- function(x, s) {
  .zinb_jcompose(x, c(
    list(log(x[, 1])),
    lapply(seq_len(s$order), function(k) (-1)^(k - 1) * factorial(k - 1) / x[, 1]^k)
  ), s)
}
.zinb_jsoftplus <- function(x, s) {
  p <- stats::plogis(x[, 1])
  q <- stats::plogis(-x[, 1])
  d <- list(pmax(x[, 1], 0) + log1p(exp(-abs(x[, 1]))), p, p * q, p * q * (q - p), p * q * (1 - 6 * p * q))
  .zinb_jcompose(x, d, s)
}
.zinb_jpoly <- function(x, coef, s) {
  out <- x * 0
  out[, 1] <- coef[length(coef)]
  if (length(coef) > 1) {
    for (k in (length(coef) - 1L):1L) {
      out <- .zinb_jmul(out, x, s)
      out[, 1] <- out[, 1] + coef[k]
    }
  }
  out
}
.zinb_jvar <- function(v, j, s) {
  out <- matrix(0, length(v), nrow(s$a))
  out[, 1] <- v
  if (s$order) out[, which(s$a[, j] == 1 & rowSums(s$a) == 1)] <- 1
  out
}

# NB link-scale jet, optionally with log(mu) subtracted before cancellation.
.zinb_nbjet <- function(y, u, v, s, strip = FALSE) {
  mu <- exp(u)
  theta <- exp(v)
  out <- matrix(0, length(y), nrow(s$a))
  out[, 1] <- stats::dnbinom(y, mu = mu, size = theta, log = TRUE) - if (strip) u else 0
  small <- y <= 100
  if (any(small)) {
    val <- (y[small] - as.integer(strip)) * u[small] - lgamma(y[small] + 1) -
      (theta[small] + y[small]) * log1p(mu[small] / theta[small])
    for (k in seq_len(max(y[small])) - 1L) val <- val + ifelse(y[small] > k, log1p(k / theta[small]), 0)
    out[small, 1] <- val
  }
  # The same convergent near-Poisson expansion as the NB derivative kernel.
  large <- !small & theta > 100 * pmax(1, y, mu)
  if (any(large)) {
    yy <- y[large]
    mm <- mu[large]
    tt <- theta[large]
    powers <- matrix(0, length(yy), 13L)
    powers[, 1] <- yy
    correction <- numeric(length(yy))
    for (r in 1:12) {
      sr <- yy * (yy / tt)^r
      for (j in 0:(r - 1L)) sr <- sr - choose(r + 1L, j) * powers[, j + 1L] / tt^(r - j)
      powers[, r + 1L] <- sr / (r + 1L)
      correction <- correction + (-1)^(r + 1L) *
        ((powers[, r + 1L] - yy * (mm / tt)^r) / r + mm * (mm / tt)^r / (r + 1L))
    }
    out[large, 1] <- stats::dpois(yy, mm, log = TRUE) + correction -
      if (strip) u[large] else 0
  }
  if (!s$order) {
    return(out)
  }
  dd <- .nbls_derivatives(y, mu, theta, max(2L, s$order))
  tri <- mgcv::trind.generator(2)
  de <- mgcv::gamlss.etamu(dd$l1, dd$l2, dd$l3, dd$l4, cbind(mu, theta),
    cbind(-1 / mu^2, -1 / theta^2), cbind(2 / mu^3, 2 / theta^3),
    cbind(-6 / mu^4, -6 / theta^4), tri$i2, tri$i3, tri$i4,
    deriv = max(0L, s$order - 1L)
  )
  for (k in 2:nrow(s$a)) {
    a <- s$a[k, ]
    r <- sum(a)
    if (a[2] == 0) out[, k] <- de[[paste0("l", r)]][, a[3] + 1L] / s$fact[k]
  }
  # Stable pure log-mean derivatives, including y=1 at tiny mu.
  p <- mu / (theta + mu)
  q <- theta / (theta + mu)
  pure <- list(
    (y - as.integer(strip)) - (theta + y) * p,
    -(theta + y) * p * q, -(theta + y) * p * q * (1 - 2 * p),
    -(theta + y) * p * q * (1 - 6 * p + 6 * p^2)
  )
  for (r in seq_len(s$order)) out[, which(s$a[, 1] == r & rowSums(s$a) == r)] <- pure[[r]] / factorial(r)
  out
}

.zinb_logjet <- function(y, eta, hurdle, order = 0L) {
  s <- .zinb_jet_spec(order)
  u <- .zinb_jvar(eta[, 1], 1, s)
  b <- .zinb_jvar(eta[, 2], 2, s)
  v <- .zinb_jvar(eta[, 3], 3, s)
  logp <- -.zinb_jsoftplus(-b, s)
  logabs <- -.zinb_jsoftplus(b, s)
  nb <- .zinb_nbjet(y, eta[, 1], eta[, 3], s, strip = hurdle)
  z <- y == 0
  if (!hurdle) {
    # log[(1-p)+p*q0] = log(1-p)+softplus(logit(p)+log(q0)).
    q <- .zinb_nbjet(rep(0, length(y)), eta[, 1], eta[, 3], s)
    out <- logp + nb
    mix <- logabs + .zinb_jsoftplus(b + q, s)
    out[z, ] <- mix[z, , drop = FALSE]
  } else {
    # log(1-q0) = log(mu) + log(log1p(r)/r) + h(a),
    # r=mu/theta, a=theta*log1p(r), h(a)=log((1-exp(-a))/a).
    # The two smooth correction functions admit stable power series at zero.
    r <- .zinb_jexp(u - v, s)
    lr <- .zinb_jsoftplus(u - v, s) # log1p(r)
    a <- .zinb_jmul(.zinb_jexp(v, s), lr, s)
    c1 <- r * 0
    lo <- r[, 1] < .01
    if (any(lo)) {
      # log(log1p(r)/r), evaluated by first forming its near-one argument.
      t <- .zinb_jpoly(r[lo, , drop = FALSE], (-1)^(0:12) / (1:13), s)
      c1[lo, ] <- .zinb_jlog(t, s)
    }
    if (any(!lo)) c1[!lo, ] <- .zinb_jlog(lr[!lo, , drop = FALSE], s) - (u - v)[!lo, , drop = FALSE]
    h <- a * 0
    lo <- a[, 1] < .1
    if (any(lo)) {
      h[lo, ] <- .zinb_jpoly(
        a[lo, , drop = FALSE],
        c(0, -.5, 1 / 24, 0, -1 / 2880, 0, 1 / 181440, 0, -1 / 9676800, 0, 1 / 479001600), s
      )
    }
    if (any(!lo)) {
      aa <- a[!lo, , drop = FALSE]
      t <- -.zinb_jexp(-aa, s)
      t[, 1] <- -expm1(-aa[, 1])
      h[!lo, ] <- .zinb_jlog(t, s) - .zinb_jlog(aa, s)
    }
    out <- logp + nb - c1 - h
    out[z, ] <- logabs[z, , drop = FALSE]
  }
  de <- list(l0 = out[, 1])
  if (order) {
    for (r in seq_len(order)) {
      # mgcv packs derivatives in lexicographic predictor-index order.
      ii <- which(rowSums(s$a) == r)
      ii <- ii[order(-s$a[ii, 1], -s$a[ii, 2])]
      de[[paste0("l", r)]] <- sweep(out[, ii, drop = FALSE], 2, s$fact[ii], "*")
    }
  }
  de
}

.zinb_ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0, d1b = 0, d2b = 0,
                     Hp = NULL, rank = 0, fh = NULL, D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
  if (is.list(X)) stop("zinb/zanb do not support discrete model matrices")
  n <- length(y)
  jj <- attr(X, "lpi")
  supplied_eta <- !is.null(eta)
  if (is.null(wt)) wt <- rep(1, n)
  if (is.null(eta)) {
    if (length(jj) != 3L) stop("zinb/zanb require three linear predictors")
    eta <- matrix(0, n, 3)
    for (j in 1:3) {
      eta[, j] <- drop(X[, jj[[j]], drop = FALSE] %*% coef[jj[[j]]])
      if (length(offset) >= j && !is.null(offset[[j]])) eta[, j] <- eta[, j] + offset[[j]]
    }
  }
  active <- wt > 0
  eta[!active, ] <- 0
  y[!active] <- 0
  if (any(!is.finite(eta)) || any(!is.finite(exp(eta[, c(1, 3)]))) ||
    any(exp(eta[, c(1, 3)]) == 0)) {
    return(list(
      l = -Inf, l0 = ifelse(active, -Inf, 0),
      lb = rep(0, ncol(X)), lbb = matrix(0, ncol(X), ncol(X))
    ))
  }
  order <- if (deriv == 0) 0L else if (deriv == 1) 2L else if (deriv < 4) 3L else 4L
  de <- .zinb_logjet(y, eta, family$hurdle, order)
  de <- lapply(de, function(x) x * wt)
  ret <- list()
  if (deriv && !is.null(jj) && !supplied_eta) {
    shared <- deriv == 4 && anyDuplicated(unlist(jj)) > 0
    tri <- family$tri
    ret <- mgcv::gamlss.gH(X, jj, de$l1, de$l2, tri$i2,
      l3 = de$l3, l4 = de$l4,
      i3 = tri$i3, i4 = tri$i4, d1b = d1b, d2b = d2b, deriv = if (shared) 2L else deriv - 1L,
      fh = fh, D = D, sandwich = sandwich
    )
    if (shared) ret$trHid2H <- .nbls_shared_trace(X, jj, de, tri, d1b, d2b, fh, D)
  }
  if (ncv) ret[c("l1", "l2", "l3")] <- de[c("l1", "l2", "l3")]
  ret$l0 <- de$l0
  ret$l <- sum(de$l0)
  ret
}

.zinb_initialize <- function(y, x, E, wt, offset, start, family) {
  if (is.null(wt)) wt <- rep(1, length(y))
  if (any(!is.finite(y) | y < 0 | y != floor(y))) stop("finite nonnegative integer counts required")
  if (any(!is.finite(wt) | wt < 0)) stop("finite nonnegative weights required")
  keep <- wt > 0
  if (!any(keep)) stop("positive-weight observations required")
  if (!any(y[keep] > 0)) stop("cannot fit an all-zero positive-weight response")
  if (all(y[keep] > 0)) warning("no observed zeros: probability may have a boundary estimate")
  if (all(y[keep & y > 0] == 1)) warning("positive counts are all one: mean and size may be unidentified")
  jj <- attr(x, "lpi")
  if (length(jj) != 3) stop("zinb/zanb require three linear predictors")
  if (!is.null(start)) {
    return(start)
  }
  off <- lapply(1:3, function(j) {
    if (length(offset) >= j && !is.null(offset[[j]])) {
      rep_len(offset[[j]], length(y))
    } else {
      numeric(length(y))
    }
  })
  # Positive-count regressions provide the count start; zeros receive a small
  # initialization weight only (never a changed likelihood weight).
  pstart <- (sum(wt[keep] * (y[keep] > 0)) + .5) / (sum(wt[keep]) + 1)
  if (!family$hurdle) pstart <- (1 + pstart) / 2
  target <- list(log(pmax(y, .5)), rep(stats::qlogis(pstart), length(y)), rep(log(2), length(y)))
  design <- rhs <- vector("list", 3)
  for (j in 1:3) {
    w <- sqrt(wt[keep] * if (j == 1) ifelse(y[keep] > 0, 1, .05) else 1)
    a <- matrix(0, sum(keep), ncol(x))
    a[, jj[[j]]] <- x[keep, jj[[j]], drop = FALSE]
    design[[j]] <- a * w
    rhs[[j]] <- (target[[j]][keep] - off[[j]][keep]) * w
  }
  if (length(E) && sum(E^2) > 0) E <- E * .01 * sqrt(sum(x[keep, , drop = FALSE]^2) / sum(E^2))
  ans <- qr.coef(qr(do.call(rbind, c(design, list(E)))), c(unlist(rhs), rep(0, nrow(E))))
  ans[!is.finite(ans)] <- 0
  ans
}

# Use mgcv's compiled deletion updates, but assemble the score locally:
# upstream QNCV assumes full-data prediction ordering and omits base offsets.
.zinb_ncv <- function(X, y, wt, nei, beta, family, llf, H = NULL, Hi = NULL, R = NULL,
                      offset = NULL, dH = NULL, db = NULL, deriv = FALSE, nt = 1) {
  if (is.null(wt)) wt <- rep(1, length(y))
  n <- length(y)
  jj <- attr(X, "lpi")
  m <- length(nei$d)
  off <- lapply(1:3, function(j) {
    if (length(offset) >= j && !is.null(offset[[j]])) {
      rep_len(offset[[j]], n)
    } else {
      numeric(n)
    }
  })
  f <- family
  f$qapprox <- FALSE
  lf <- llf
  lf$gamma <- 1
  native <- utils::getFromNamespace("gamlss.ncv", "mgcv")
  ret <- native(X, y, wt, nei, beta, f, lf,
    H = H, Hi = Hi, R = R, offset = off,
    dH = dH, db = db, deriv = deriv, nt = nt
  )
  ecv <- attr(ret$NCV, "eta.cv")
  dcv <- attr(ret$NCV, "deta.cv")
  eta <- sapply(1:3, function(j) drop(X[, jj[[j]], drop = FALSE] %*% beta[jj[[j]]]) + off[[j]])
  base <- family$ll(y, X, beta, wt, family, eta = eta, deriv = if (deriv > 0) 3 else 1, ncv = TRUE)
  ix <- nei$d
  gamma <- llf$gamma
  if (is.null(gamma)) gamma <- 1
  score0 <- -sum(base$l0[ix])
  delta <- ecv - eta[ix, , drop = FALSE]
  tri <- family$tri
  if (isTRUE(family$qapprox)) {
    Hd <- matrix(0, m, 3)
    for (j in 1:3) {
      for (k in 1:3) {
        Hd[, j] <- Hd[, j] + base$l2[ix, tri$i2[j, k]] * delta[, k]
      }
    }
    score <- score0 - gamma * sum(base$l1[ix, , drop = FALSE] * delta + .5 * Hd * delta)
  } else {
    cv <- family$ll(y[ix], X[ix, , drop = FALSE], beta, wt[ix], family, eta = ecv, deriv = 1, ncv = TRUE)
    score <- -gamma * cv$l + (1 - gamma) * score0
  }
  grad <- Vg <- NULL
  if (deriv > 0) {
    ns <- ncol(db)
    g <- matrix(0, m, ns)
    d0 <- lapply(1:3, function(j) X[ix, jj[[j]], drop = FALSE] %*% db[jj[[j]], , drop = FALSE])
    for (j in 1:3) {
      dc <- dcv[seq_len(m) + (j - 1L) * m, , drop = FALSE]
      g <- g - (1 - gamma) * base$l1[ix, j] * d0[[j]]
      if (isTRUE(family$qapprox)) {
        g <- g - gamma * (base$l1[ix, j] + Hd[, j]) * dc
        td <- numeric(m)
        for (k in 1:3) {
          for (l in 1:3) {
            td <- td + base$l3[ix, tri$i3[j, k, l]] * delta[, k] * delta[, l]
          }
        }
        g <- g - .5 * gamma * td * d0[[j]]
      } else {
        g <- g - gamma * cv$l1[, j] * dc
      }
    }
    grad <- colSums(g)
    Vg <- crossprod(g)
  }
  attr(score, "eta.cv") <- ecv
  if (deriv != 0) attr(score, "deta.cv") <- dcv
  list(NCV = score, NCV1 = grad, error = ret$error, Vg = Vg)
}

.zinb_logpositive <- function(mu, theta) log(-expm1(-theta * log1p(mu / theta)))
.zinb_quantile <- function(p, mu, hurdle) {
  n <- nrow(mu)
  p <- rep_len(p, n)
  prob <- mu[, 2]
  lp <- if (hurdle) .zinb_logpositive(mu[, 1], mu[, 3]) else rep(0, n)
  tail <- log1p(-p) - log(prob) + lp
  out <- numeric(n)
  use <- prob > 0 & !is.na(p) & p > 0
  # Tail probabilities >= P(N>0) correspond to the atom at zero.
  threshold <- .zinb_logpositive(mu[, 1], mu[, 3])
  use <- use & tail < threshold
  out[use] <- stats::qnbinom(tail[use],
    mu = mu[use, 1], size = mu[use, 3],
    lower.tail = FALSE, log.p = TRUE
  )
  out[p < 0 | p > 1 | is.na(p)] <- NaN
  out
}
.zinb_cdf <- function(q, mu, hurdle, logp = FALSE) {
  q <- rep_len(q, nrow(mu))
  ls <- log(mu[, 2]) + stats::pnbinom(q, mu = mu[, 1], size = mu[, 3], lower.tail = FALSE, log.p = TRUE)
  if (hurdle) ls <- ls - .zinb_logpositive(mu[, 1], mu[, 3])
  ls <- pmin(ls, 0)
  ans <- -expm1(ls)
  ans[q < 0] <- 0
  if (logp) log(ans) else ans
}
.zinb_moments <- function(mu, hurdle) {
  m <- mu[, 1]
  p <- mu[, 2]
  t <- mu[, 3]
  pos <- if (hurdle) exp(log(m) - .zinb_logpositive(m, t)) else m
  conditional <- pmax(0, pos * (1 + m * (1 + 1 / t) - pos))
  list(mean = p * pos, variance = p * conditional + p * (1 - p) * pos^2)
}
.zinb_saturated <- function(y, theta, hurdle) {
  out <- numeric(length(y))
  ii <- which(y > 0)
  if (!hurdle && length(ii)) {
    out[ii] <- .zinb_nbjet(
      y[ii], log(y[ii]),
      log(theta[ii]), .zinb_jet_spec(0)
    )[, 1]
  }
  if (hurdle) {
    for (i in which(y > 1)) {
      # Solve on log scale: tiny sizes can require mu far below machine epsilon.
      root <- stats::uniroot(
        function(u) {
          exp(u - .zinb_logpositive(exp(u), theta[i])) - y[i]
        },
        c(min(log(y[i]), log(theta[i])) - 40, log(y[i])),
        tol = 1e-10
      )$root
      out[i] <- .zinb_logjet(y[i], matrix(c(root, 0, log(theta[i])), 1), TRUE)$l0 + log(2)
    }
  }
  out
}
.zinb_residuals <- function(object, type = c("deviance", "pearson", "response")) {
  type <- match.arg(type)
  mu <- object$fitted.values
  w <- object$prior.weights
  mom <- .zinb_moments(mu, object$family$hurdle)
  r <- object$y - mom$mean
  if (type == "response") {
    return(r)
  }
  if (type == "pearson") {
    ans <- numeric(length(r))
    ii <- w > 0
    ans[ii] <- sqrt(w[ii]) * r[ii] / sqrt(mom$variance[ii])
    return(ans)
  }
  l <- .zinb_logjet(object$y, object$linear.predictors, object$family$hurdle)$l0
  sat <- .zinb_saturated(object$y, mu[, 3], object$family$hurdle)
  ans <- numeric(length(r))
  ii <- w > 0
  ans[ii] <- sign(r[ii]) * sqrt(pmax(0, 2 * w[ii] * (sat[ii] - l[ii])))
  ans
}
.zinb_null <- function(object, offset = NULL) {
  keep <- object$prior.weights > 0
  y <- object$y[keep]
  w <- object$prior.weights[keep]
  theta <- object$fitted.values[keep, 3]
  off <- lapply(1:2, function(j) {
    if (length(offset) >= j && !is.null(offset[[j]])) {
      rep_len(offset[[j]], length(keep))[keep]
    } else {
      numeric(sum(keep))
    }
  })
  initial <- c(
    stats::weighted.mean(object$linear.predictors[keep, 1] - off[[1]], w),
    stats::weighted.mean(object$linear.predictors[keep, 2] - off[[2]], w)
  )
  fn <- function(b) {
    eta <- cbind(b[1] + off[[1]], b[2] + off[[2]], log(theta))
    if (any(!is.finite(exp(eta[, 1]))) || any(exp(eta[, 1]) == 0)) {
      return(1e100)
    }
    val <- -sum(w * .zinb_logjet(y, eta, object$family$hurdle)$l0)
    if (is.finite(val)) val else 1e100
  }
  opt <- stats::optim(initial, fn, method = "BFGS", control = list(reltol = 1e-10, maxit = 200))
  if (opt$convergence != 0) warning("null-deviance optimization did not converge")
  2 * (sum(w * .zinb_saturated(y, theta, object$family$hurdle)) + opt$value)
}
