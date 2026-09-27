#' Beta-binomial mean and precision family
#'
#' A two-predictor family for grouped binomial data with extra-binomial
#' variation. The predictors model `logit(mu)` and `log(phi)`, where the latent
#' success probability has distribution `Beta(mu * phi, (1 - mu) * phi)`.
#'
#' @param link Two links, which must be `"logit"` and `"log"`, respectively.
#' @return A `general.family` object for [mgcv::gam()].
#' @md
#' @details Supply a list of two formulae, with the response in the first.
#'   Start with `~ 1` for the precision predictor. Both predictors support
#'   smooths and offsets. For trial total `m`, the count mean is `m * mu` and
#'   variance is `m * mu * (1 - mu) * (m + phi) / (1 + phi)`.
#'
#'   Supply either `cbind(successes, failures)`, with optional nonnegative
#'   likelihood weights, or a proportion with integer trial totals in
#'   `weights`. For proportions, weights are trial totals only: they do not
#'   multiply the likelihood again. Counts and implied successes must be
#'   integers (up to rounding tolerance). Zero-trial rows contribute nothing.
#'   Internally, `object$y` holds success counts, `object$family$trials` holds
#'   totals, and `object$prior.weights` holds likelihood weights.
#'
#'   Fitted values and response predictions have two columns, in the order
#'   `mu`, `phi`; link predictions contain their logit and log transforms.
#'   Multiply the first response column by trial totals for expected counts.
#'   Response standard errors use the delta method.
#'
#'   REML supports outer Newton and BFGS optimization with analytic derivatives
#'   through fourth order. NCV supports row-level neighbourhoods:
#'   each row is an entire binomial group. Correlation between rows requires
#'   appropriate neighbourhoods. The one-step NCV approximation can be poor
#'   when deleting a neighbourhood removes much of the precision information.
#'   ML, QNCV, `bam()`, and `gamm()` are not supported. Note that mgcv can silently
#'   change a request for ML into REML for a general family.
#'
#'   Bernoulli rows provide no precision information; at least one informative
#'   row with two or more trials is required. All-success or all-failure data
#'   are rejected. Large precision approaches the binomial limit; underdispersion
#'   cannot be represented. Precision is not artificially bounded. Separation,
#'   weak precision information, and boundary estimates can invalidate standard
#'   errors and the Laplace approximation used by REML. A warning is issued
#'   when the precision predictor has a standard error exceeding 5 on its log
#'   scale, or its information block is numerically singular. This is a
#'   diagnostic, not a guarantee that other fits are well identified.
#'
#'   Response and Pearson residuals are on the count scale. Deviance maximizes
#'   each row's likelihood over its mean with fitted precision held fixed;
#'   that mean is not generally the observed proportion. Null deviance fits a
#'   constant mean predictor with mean offsets and fitted precisions retained.
#'   These are conditional mean-fit diagnostics. Simulation, CDF, and quantile
#'   callbacks are not supplied. Mathematical details are in the source file
#'   `inst/maths/betabinomial.md`.
#' @export
#' @examples
#' set.seed(42)
#' dat <- data.frame(x = runif(300), z = runif(300), m = 20)
#' mu <- plogis(sin(2 * pi * dat$x))
#' phi <- exp(1 + dat$z)
#' dat$y <- rbinom(300, dat$m, rbeta(300, mu * phi, (1 - mu) * phi))
#' fit <- mgcv::gam(list(cbind(y, m - y) ~ s(x, k = 6), ~ s(z, k = 5)),
#'                  family = betabinomial(), method = "REML", data = dat)
#' head(predict(fit, type = "response"))
betabinomial <- function(link = list("logit", "log")) {
  if (length(link) != 2L ||
      !all(vapply(link, function(x) is.character(x) && length(x) == 1L &&
                    !is.na(x), logical(1))) ||
      link[[1L]] != "logit" || link[[2L]] != "log")
    stop('betabinomial requires links "logit" and "log"')
  linfo <- lapply(link, stats::make.link)
  # Do not use make.link's floors: likelihood derivatives must describe the
  # actual inverse links, including at extreme predictors.
  linfo[[1L]]$linkinv <- stats::plogis
  linfo[[1L]]$mu.eta <- function(eta) stats::plogis(eta) * stats::plogis(-eta)
  linfo[[2L]]$linkinv <- linfo[[2L]]$mu.eta <- function(eta) exp(eta)
  structure(list(family = "betabinomial", link = unlist(link), nlp = 2L,
    linfo = linfo, tri = mgcv::trind.generator(2), ll = .bb_ll,
    preinitialize = .bb_preinitialize,
    initialize = expression({
      n <- rep.int(1, nobs)
      start <- family$init(y, x, E, weights, offset, start, family)
    }),
    init = .bb_initialize,
    postproc = expression({
      colnames(object$fitted.values) <- c("mu", "phi")
      bbX <- getFromNamespace("Sl.initial.repara", "mgcv")(
        G$Sl, G$X, inverse = TRUE, cov = FALSE, both.sides = FALSE)
      object <- object$family$finish(object, bbX, G$offset)
    }),
    finish = .bb_finish, residuals = .bb_residuals, ncv = .bb_ncv,
    sandwich = function(y, X, coef, wt, family, offset = NULL) {
      family$ll(y, X, coef, wt, family, offset = offset,
                deriv = 1, sandwich = TRUE)$lbb
    },
    d2link = 1, d3link = 1, d4link = 1, ls = 1,
    available.derivs = 2L, discrete.ok = FALSE),
    class = c("general.family", "extended.family", "family"))
}

.bb_preinitialize <- function(G) {
  if (isTRUE(G$family$qapprox)) stop("betabinomial supports NCV but not QNCV")
  y <- G$y; wt <- G$w
  nr <- NROW(y)
  if (is.null(wt)) wt <- rep.int(1, nr)
  if (any(!is.finite(wt) | wt < 0))
    stop("betabinomial requires finite nonnegative weights")
  integerish <- function(z) all(is.finite(z) & z >= 0 &
    abs(z - round(z)) <= pmin(1e-7, 64 * .Machine$double.eps * pmax(1, abs(z))))
  if (is.matrix(y)) {
    if (ncol(y) != 2L || !integerish(y))
      stop("betabinomial requires two columns of nonnegative integer counts")
    y <- round(y)
    m <- rowSums(y); y <- y[, 1L]
  } else {
    if (any(!is.finite(y) | y < 0 | y > 1) || !integerish(wt) ||
        !integerish(y * wt))
      stop("betabinomial proportions require integer trial totals and implied successes")
    m <- round(wt); y <- round(y * wt); wt <- rep.int(1, nr)
  }
  if (any(!is.finite(m) | m > 2^53))
    stop("betabinomial trial totals must be finite and at most 2^53")
  keep <- wt > 0 & m > 0
  if (!any(wt > 0 & m >= 2))
    stop("betabinomial precision requires positive-weight rows with at least two trials")
  if (all(y[keep] == 0) || all(y[keep] == m[keep]))
    stop("betabinomial cannot fit all-failure or all-success informative responses")
  family <- G$family
  family$trials <- m
  # A scalar success-count response avoids mgcv's length(matrix) assumptions.
  # Keep totals in a per-fit family copy, never in shared mutable state.
  list(y = y, w = wt, family = family, offxset = G$offset)
}

.bb_softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

# F(t,n) = log Gamma(exp(t)+n) - log Gamma(exp(t)) - n*t.
# Columns contain F and its first four derivatives with respect to t.
.bb_rising <- function(t, n, order = 4L) {
  out <- matrix(0, length(n), order + 1L)
  ii <- which(n > 1 & n <= 64)
  if (length(ii)) {
    # Bound temporary matrix sizes while vectorizing the short recurrences.
    # In particular, saturated one-row likelihoods should not loop in R.
    for (first in seq.int(1L, length(ii), by = 256L)) {
      jj <- ii[seq.int(first, min(first + 255L, length(ii)))]
      k <- seq_len(max(n[jj]) - 1L)
      z <- outer(t[jj], log(k), function(a, b) b - a)
      z[outer(n[jj], k, `<=`)] <- -Inf
      r <- stats::plogis(z); s <- stats::plogis(-z)
      out[jj, 1L] <- rowSums(.bb_softplus(z))
      if (order > 0L) out[jj, 2L] <- -rowSums(r)
      if (order > 1L) out[jj, 3L] <- rowSums(r * s)
      if (order > 2L) out[jj, 4L] <- rowSums(r * s * (r - s))
      if (order > 3L) out[jj, 5L] <- rowSums(r * s * (1 - 6 * r * s))
    }
  }
  ii <- which(n > 64 & log(n) - t < log(0.05))
  if (length(ii)) {
    nn <- n[ii]; ratio <- exp(log(nn) - t[ii])
    # Normalized Faulhaber recurrence: powers[,r+1] = sum (k/n)^r.
    powers <- matrix(0, length(ii), 17L); powers[, 1L] <- nn
    for (r in 1:16) {
      sr <- nn
      for (j in 0:(r - 1L)) sr <- sr - choose(r + 1L, j) *
        powers[, j + 1L] / nn^(r - j)
      powers[, r + 1L] <- sr / (r + 1L)
      term <- (-1)^(r + 1L) * powers[, r + 1L] * ratio^r / r
      for (q in 0:order) out[ii, q + 1L] <- out[ii, q + 1L] + (-r)^q * term
    }
  }
  ii <- which(n > 64 & log(n) - t >= log(0.05))
  if (length(ii)) {
    a <- exp(t[ii]); nn <- n[ii]
    # Shifting away the k=0 term avoids singular polygamma values as a -> 0.
    out[ii, 1L] <- lgamma(a + nn) - lgamma(a + 1) - (nn - 1) * t[ii]
    large <- a >= 8
    if (any(large)) {
      aa <- a[large]; mm <- nn[large]
      correction <- function(x) {
        z <- 1 / x
        z / 12 - z^3 / 360 + z^5 / 1260 - z^7 / 1680 +
          z^9 / 1188 - 691 * z^11 / 360360 + z^13 / 156
      }
      out[ii[large], 1L] <- (aa + mm - 0.5) * log1p(mm / aa) - mm +
        correction(aa + mm) - correction(aa)
    }
    if (order > 0L) {
      # Stirling numbers of the second kind convert a-derivatives to t-derivatives.
      stirling <- list(1, c(1, 1), c(1, 3, 1), c(1, 7, 6, 1))
      dd <- lapply(seq_len(order), function(j)
        a^j * (psigamma(a + nn, deriv = j - 1L) - psigamma(a + 1, deriv = j - 1L)))
      for (q in seq_len(order)) {
        out[ii, q + 1L] <- Reduce(`+`, Map(`*`, dd[seq_len(q)], stirling[[q]]))
        if (q == 1L) out[ii, q + 1L] <- out[ii, q + 1L] - (nn - 1)
      }
    }
  }
  out
}

# Packed predictor derivatives: eta1^r, eta1^(r-1)*eta2, ..., eta2^r.
# Use partial Bell polynomials to compose F(log(mu)+eta2) and F(log(1-mu)+eta2).
.bb_derivatives <- function(y, m, eta, order = 4L) {
  u <- stats::plogis(eta[, 1L]); v <- stats::plogis(-eta[, 1L])
  lu <- -.bb_softplus(-eta[, 1L]); lv <- -.bb_softplus(eta[, 1L])
  fa <- .bb_rising(lu + eta[, 2L], y, order)
  fb <- .bb_rising(lv + eta[, 2L], m - y, order)
  fc <- .bb_rising(eta[, 2L], m, order)
  out <- list(l0 = lchoose(m, y) + y * lu + (m - y) * lv + fa[, 1L] + fb[, 1L] - fc[, 1L])
  if (!order) return(out)
  du <- list(v, -u * v, -u * v * (v - u), -u * v * (1 - 6 * u * v))
  dv <- du; dv[[1L]] <- -u
  bell <- function(d) list(list(1), list(d[[1L]]),
    list(d[[2L]], d[[1L]]^2),
    list(d[[3L]], 3 * d[[1L]] * d[[2L]], d[[1L]]^3),
    list(d[[4L]], 4 * d[[1L]] * d[[3L]] + 3 * d[[2L]]^2,
         6 * d[[1L]]^2 * d[[2L]], d[[1L]]^4))
  bu <- bell(du); bv <- bell(dv)
  for (r in seq_len(order)) {
    dd <- matrix(0, length(y), r + 1L)
    for (q in 0:r) {
      p <- r - q
      if (p == 0L) {
        dd[, q + 1L] <- fa[, q + 1L] + fb[, q + 1L] - fc[, q + 1L]
      } else {
        for (k in seq_len(p)) dd[, q + 1L] <- dd[, q + 1L] +
          bu[[p + 1L]][[k]] * fa[, q + k + 1L] + bv[[p + 1L]][[k]] * fb[, q + k + 1L]
        if (q == 0L) dd[, 1L] <- dd[, 1L] + y * du[[p]] + (m - y) * dv[[p]]
      }
    }
    out[[paste0("l", r)]] <- dd
  }
  out
}

.bb_ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                   d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                   D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
  if (is.list(X)) stop("betabinomial does not support discrete model matrices")
  jj <- attr(X, "lpi")
  if (length(jj) != 2L) stop("betabinomial requires two linear predictors")
  nr <- length(y); m <- family$trials
  if (length(m) != nr) stop("betabinomial trial totals are not aligned with the response")
  if (is.null(wt)) wt <- rep.int(1, nr)
  if (is.null(eta)) {
    eta <- matrix(0, nr, 2L)
    for (j in 1:2) {
      eta[, j] <- drop(X[, jj[[j]], drop = FALSE] %*% coef[jj[[j]]])
      if (length(offset) >= j && !is.null(offset[[j]])) eta[, j] <- eta[, j] + offset[[j]]
    }
  }
  active <- wt > 0 & m > 0
  eta[!active, ] <- 0
  if (any(!is.finite(eta)) || (any(eta[active, 2L] > log(.Machine$double.xmax)) ||
      any(exp(eta[active, 2L]) == 0))) {
    l0 <- numeric(nr); l0[active] <- -Inf
    return(list(l = -Inf, l0 = l0, lb = rep(0, ncol(X)),
                lbb = matrix(0, ncol(X), ncol(X))))
  }
  yy <- y; mm <- m; yy[!active] <- mm[!active] <- 0
  order <- if (deriv == 0L) 0L else if (deriv == 1L) 2L else if (deriv < 4L) 3L else 4L
  de <- .bb_derivatives(yy, mm, eta, order)
  for (k in seq_along(de)) de[[k]] <- de[[k]] * wt
  ret <- list()
  if (deriv > 0L) {
    tri <- family$tri
    shared4 <- deriv == 4L && anyDuplicated(unlist(jj)) > 0L
    ret <- mgcv::gamlss.gH(X, jj, de$l1, de$l2, tri$i2,
      l3 = de$l3, i3 = tri$i3, l4 = de$l4, i4 = tri$i4,
      d1b = d1b, d2b = d2b, deriv = if (shared4) 2L else deriv - 1L, fh = fh, D = D,
      sandwich = sandwich)
    if (shared4) ret$trHid2H <- .bb_shared_trace(X, jj, de, tri, d1b, d2b, fh, D)
    if (ncv) {
      ret$l1 <- de$l1; ret$l2 <- de$l2; ret$l3 <- de$l3
    }
  }
  ret$l <- sum(de$l0); ret$l0 <- de$l0
  ret
}

# gamlss.gH's fourth-order trace uses assignment into expanded design blocks.
# Those blocks overlap for shared coefficients, so explicitly sum all predictor
# pairs here. Other contractions still use mgcv's standard helper.
.bb_shared_trace <- function(X, jj, de, tri, d1b, d2b, fh, D) {
  if (is.list(fh)) {
    inv <- ifelse(fh$values > 0, 1 / fh$values, 0)
    Hi <- fh$vectors %*% (inv * t(fh$vectors))
  } else {
    ipiv <- order(attr(fh, "pivot"))
    Hi <- chol2inv(fh)[ipiv, ipiv, drop = FALSE]
  }
  Hi <- D * t(D * Hi)
  d1eta <- lapply(jj, function(j) X[, j, drop = FALSE] %*% d1b[j, , drop = FALSE])
  d2eta <- lapply(jj, function(j) X[, j, drop = FALSE] %*% d2b[j, , drop = FALSE])
  ans <- numeric(ncol(d2b))
  for (i in 1:2) for (j in i:2) {
    v <- rowSums((X[, jj[[i]], drop = FALSE] %*% Hi[jj[[i]], jj[[j]], drop = FALSE]) *
                   X[, jj[[j]], drop = FALSE])
    kk <- 0L
    for (k in seq_len(ncol(d1b))) for (l in k:ncol(d1b)) {
      kk <- kk + 1L
      value <- numeric(nrow(X))
      for (q in 1:2) {
        value <- value + de$l3[, tri$i3[i, j, q]] * d2eta[[q]][, kk]
        for (s in 1:2) value <- value + de$l4[, tri$i4[i, j, q, s]] *
          d1eta[[q]][, k] * d1eta[[s]][, l]
      }
      ans[kk] <- ans[kk] + (if (i == j) 1 else 2) * sum(v * value)
    }
  }
  ans
}

.bb_ncv <- function(X, y, wt, nei, beta, family, llf, H = NULL, Hi = NULL,
                    R = NULL, offset = NULL, dH = NULL, db = NULL,
                    deriv = FALSE, nt = 1) {
  response <- y; original <- family; lpi <- attr(X, "lpi")
  family$ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                        ..., eta = NULL) {
    attr(X, "lpi") <- lpi
    subfamily <- original
    subfamily$trials <- original$trials[y]
    # The helper's full-data score (deriv=0, used when gamma != 1) supplies
    # X beta without offsets. Reconstruct that predictor through the regular
    # likelihood path. Its held-out eta already includes offsets.
    if (deriv == 0L) eta <- NULL
    original$ll(response[y], X, coef, wt, subfamily, offset = offset,
                deriv = deriv, ..., eta = eta)
  }
  # gamlss.ncv indexes its response as a vector, including repeated nei$d.
  # Row identifiers let the callback subset successes AND totals together.
  getFromNamespace("gamlss.ncv", "mgcv")(X, seq_along(y), wt, nei, beta,
    family, llf, H = H, Hi = Hi, R = R, offset = offset, dH = dH,
    db = db, deriv = deriv, nt = nt)
}

.bb_initialize <- function(y, x, E, wt, offset, start, family) {
  if (!is.null(start)) return(start)
  m <- family$trials; keep <- wt > 0 & m > 0
  jj <- attr(x, "lpi")
  target <- list(stats::qlogis((y + 0.5) / (m + 1)), rep(log(10), length(y)))
  design <- lapply(jj, function(ind) {
    a <- matrix(0, sum(keep), ncol(x))
    a[, ind] <- x[keep, ind, drop = FALSE]
    a * sqrt(wt[keep])
  })
  ee <- E
  if (length(ee) && sum(ee^2) > 0) ee <- ee *
    (0.01 * sqrt(sum(x[keep, , drop = FALSE]^2) / sum(ee^2)))
  rhs <- unlist(lapply(1:2, function(j) {
    off <- if (length(offset) >= j && !is.null(offset[[j]]))
      rep_len(offset[[j]], length(y)) else numeric(length(y))
    (target[[j]][keep] - off[keep]) * sqrt(wt[keep])
  }))
  ans <- qr.coef(qr(do.call(rbind, c(design, list(ee)))), c(rhs, rep(0, nrow(ee))))
  ans[!is.finite(ans)] <- 0
  ans
}

.bb_saturated <- function(y, m, logphi) {
  vapply(seq_along(y), function(i) {
    if (m[i] == 0 || y[i] == 0 || y[i] == m[i]) return(0)
    objective <- function(eta) .bb_derivatives(y[i], m[i],
      matrix(c(eta, logphi[i]), 1L), 0L)$l0
    # For interior counts, the maximizing mean lies between y/m and 1/2.
    ends <- sort(c(stats::qlogis(y[i] / m[i]), 0))
    if (diff(ends) < 1e-10) return(objective(mean(ends)))
    stats::optimize(objective, ends, maximum = TRUE, tol = 1e-10)$objective
  }, numeric(1))
}

.bb_residuals <- function(object, type = c("deviance", "pearson", "response")) {
  type <- match.arg(type)
  y <- object$y; m <- object$family$trials
  mu <- object$fitted.values[, 1L]; phi <- object$fitted.values[, 2L]
  r <- y - m * mu
  if (type == "deviance") return(sign(r) * sqrt(object$bb.deviance))
  if (type == "pearson") {
    variance <- m * mu * (1 - mu) * (1 + (m - 1) / (1 + phi))
    keep <- object$prior.weights > 0 & m > 0
    r[!keep] <- 0
    r[keep] <- sqrt(object$prior.weights[keep]) * r[keep] / sqrt(variance[keep])
  }
  r
}

.bb_finish <- function(object, X, offset) {
  y <- object$y; m <- object$family$trials; wt <- object$prior.weights
  keep <- wt > 0 & m > 0
  eta <- object$linear.predictors
  sat <- numeric(length(y))
  sat[keep] <- .bb_saturated(y[keep], m[keep], eta[keep, 2L])
  loglik <- .bb_derivatives(y[keep], m[keep], eta[keep, , drop = FALSE], 0L)$l0
  object$bb.deviance <- numeric(length(y))
  object$bb.deviance[keep] <- pmax(0, 2 * wt[keep] * (sat[keep] - loglik))
  object$deviance <- sum(object$bb.deviance)
  off <- if (length(offset) && !is.null(offset[[1L]]))
    rep_len(offset[[1L]], length(y))[keep] else numeric(sum(keep))
  objective <- function(a) -sum(wt[keep] * .bb_derivatives(y[keep], m[keep],
    cbind(a + off, eta[keep, 2L]), 0L)$l0)
  initial <- stats::qlogis(sum(wt[keep] * y[keep]) / sum(wt[keep] * m[keep])) -
    stats::weighted.mean(off, wt[keep])
  null <- stats::optim(initial, objective, method = "BFGS", control = list(reltol = 1e-10))
  object$null.deviance <- max(0, 2 * (sum(wt[keep] * sat[keep]) + null$value))
  if (null$convergence != 0L) warning("betabinomial conditional null fit did not converge")
  ind <- attr(X, "lpi")[[2L]]
  info <- -object$family$ll(y, X, object$coefficients, wt, object$family,
    offset = offset, deriv = 1L)$lbb[ind, ind, drop = FALSE]
  eig <- eigen(info, symmetric = TRUE, only.values = TRUE)$values
  xp <- X[keep, attr(X, "lpi")[[2L]], drop = FALSE]
  vp <- object$Vp[attr(X, "lpi")[[2L]], attr(X, "lpi")[[2L]], drop = FALSE]
  se <- sqrt(pmax(0, rowSums((xp %*% vp) * xp)))
  if (any(!is.finite(se)) || any(se > 5) ||
      min(eig) <= max(1, max(abs(eig))) * 1e-10)
    warning("betabinomial precision is weakly identified; boundary estimates and standard errors may be unreliable")
  object
}
