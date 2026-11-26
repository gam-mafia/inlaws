#' Logistic-normal additive models for compositions
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' A general family for strictly positive compositions in [mgcv::gam()].
#' The additive log-ratios have a multivariate normal distribution with
#' separate mean and residual standard deviation predictors and a common,
#' unrestricted residual correlation matrix.
#'
#' @param D Number of response components, at least two.
#' @param refcat Reference component, given by its column index or unique name.
#'   The default is the first column.
#' @param lkj Shape of the LKJ-style correlation penalty, a finite number at
#'   least 1. The default of 1 leaves the correlations unpenalized. Larger
#'   values give stronger shrinkage toward zero residual correlations.
#' @param formula A two-sided formula specifying the matrix response and a
#'   common right-hand side for the separate log-ratio mean predictors.
#' @param family A family returned by `logistic_normal()`.
#' @param scale A one-sided formula used for each log-standard-deviation
#'   predictor, defaulting to an intercept.
#' @return `logistic_normal()` returns a general family for `gam()`.
#'   `logistic_normal_formula()` returns a list of formulas for that family.
#' @details
#' For `m = D - 1`, the first `m` formulas model `log(y_j / y_refcat)`
#' means, in non-reference response-column order. Give the matrix response on
#' the first formula. The next `m` formulas model log standard deviations in
#' the same order. The last `m * (m - 1) / 2` formulas are independent
#' intercepts. They correspond to the strictly lower triangle of a unit lower
#' triangular matrix in row order: `(2,1), (3,1), (3,2), ...`. Normalizing
#' each row gives the Cholesky factor of the correlation matrix. With standard
#' deviations in a diagonal matrix `S`, the residual covariance is
#' `S %*% correlation %*% S`. This is the same statistical parameterization
#' as `brms::logistic_normal()`, fitted here by penalized likelihood.
#'
#' With `lkj > 1`, the fit adds `(lkj - 1) * log(det(R))` to the fitting
#' objective once, where `R` is the residual correlation matrix. For two
#' parts there is only one log-ratio and the penalty is zero. NCV and QNCV
#' retain the penalty in each deletion update and score held-out rows using
#' the response likelihood. `logLik()` and the likelihood term in `AIC()`
#' use the response likelihood at the fitted coefficients. `gam()` computes
#' effective degrees of freedom and coefficient covariance using the
#' curvature of the penalized fit.
#'
#' The formula helper repeats the mean and scale right-hand sides and appends
#' the correlation intercepts. Each repetition has its own coefficients and
#' smoothing parameters. Supply an explicit list to model different terms in
#' different log-ratios. For three components that is, for example,
#' `list(y ~ s(x), ~ z, ~ 1, ~ s(x), ~ 1)`.
#'
#' Response rows must sum to one and all entries must be finite and strictly
#' between zero and one. Rounding-level discrepancies are normalized. There
#' must be at least `D` positive-weight rows. Prior weights multiply each
#' row's log likelihood. The likelihood includes the log-ratio Jacobian.
#'
#' `gam()` supports REML with outer BFGS, and NCV and QNCV with outer BFGS.
#' NCV deletes complete composition rows. Grouped neighbourhoods use
#' `jackknife = FALSE`; the jackknife covariance is available for the
#' default leave-one-out neighbourhoods.
#' On mgcv 1.9-3 and 1.9-4, unpenalized NCV or QNCV requires
#' `nei = list(jackknife = FALSE)`.
#' This family is intended for modest numbers of parts. Avoid large `D`:
#' the number of correlation parameters grows quadratically with `D` and
#' smoothing selection requires higher-order derivatives and dense matrix
#' calculations.
#'
#' Link predictions contain the mean, log-standard-deviation and correlation
#' predictors. The fitted correlation matrix is in
#' `fit$family$data$correlation`. `predict(..., type = "response")` returns
#' the inverse log-ratio transform of the fitted means, giving the
#' compositional center in response-column order. The arithmetic expected
#' proportions can be estimated by averaging jointly simulated compositions.
#' Response prediction standard errors use the delta method with the full
#' coefficient covariance. `fitted()` retains the link-scale predictors for
#' the family's joint simulation callback `rd(mu, wt, scale)`.
#'
#' Response residuals are differences from the compositional center. Pearson
#' residuals are unsigned Mahalanobis distances on the log-ratio scale,
#' multiplied by the square root of the prior weight. Deviance residuals
#' are undefined because the saturated likelihood with free scales is unbounded.
#' @template ncv-neighbourhoods
#' @seealso [dirichlet()], [mgcv::mvn()]
#' @export
#' @examples
#' set.seed(42)
#' dat <- data.frame(x = runif(120))
#' fam <- logistic_normal(D = 3, lkj = 2)
#' eta <- cbind(sin(2 * pi * dat$x), dat$x, -0.5, -0.8, 0.6)
#' dat$y <- fam$rd(eta, NULL, 1)
#' fit <- mgcv::gam(logistic_normal_formula(y ~ s(x, k = 6), fam),
#'                  family = fam, data = dat, method = "REML")
#' head(predict(fit, type = "response"))
logistic_normal <- function(D = 2, refcat = 1, lkj = 1) {
  if (length(D) != 1L || !is.numeric(D) || !is.finite(D) ||
      D < 2 || D != floor(D))
    stop("D must be an integer of at least two", call. = FALSE)
  named_reference <- is.character(refcat) && length(refcat) == 1L &&
    !is.na(refcat) && nzchar(refcat)
  indexed_reference <- is.numeric(refcat) && length(refcat) == 1L &&
    is.finite(refcat) && refcat == floor(refcat) && refcat >= 1 && refcat <= D
  if (!named_reference && !indexed_reference)
    stop("refcat must index a response column or name a response component", call. = FALSE)
  if (length(lkj) != 1L || !is.numeric(lkj) || !is.finite(lkj) || lkj < 1)
    stop("lkj must be a finite number at least 1", call. = FALSE)

  m <- D - 1L
  q <- 2L * m + choose(m, 2L)
  correlations <- seq_len(choose(m, 2L)) + 2L * m

  preinitialize <- function(G) {
    y <- G$y
    if (!is.matrix(y) || !is.numeric(y) || ncol(y) != D)
      stop("logistic_normal requires a numeric composition matrix with D columns", call. = FALSE)
    if (any(!is.finite(y) | y <= 0 | y >= 1))
      stop("logistic_normal response entries must be finite and strictly between zero and one", call. = FALSE)
    totals <- rowSums(y)
    if (any(abs(totals - 1) > sqrt(.Machine$double.eps)))
      stop("logistic_normal response rows must sum to one", call. = FALSE)
    if (any(!is.finite(G$w) | G$w < 0) || sum(G$w > 0) < D)
      stop("logistic_normal requires finite nonnegative weights and at least D positive-weight rows", call. = FALSE)

    ref <- if (named_reference) which(colnames(y) == refcat) else as.integer(refcat)
    if (length(ref) != 1L)
      stop("refcat must name a unique response column", call. = FALSE)

    jj <- attr(G$X, "lpi")
    for (j in correlations) {
      k <- jj[[j]]
      has_offset <- length(G$offset) >= j && any(G$offset[[j]] != 0)
      if (length(k) != 1L || sum(unlist(jj) == k) != 1L ||
          any(G$X[, k] != G$X[1L, k]) || G$X[1L, k] == 0 ||
          has_offset) {
        stop("logistic_normal correlation formulas must be separate intercepts without offsets", call. = FALSE)
      }
    }
    fam <- G$family
    fam$data <- list(refcat = ref, component.names = colnames(y))
    fam$rd <- .logistic_normal_rng(m, ref, colnames(y))
    fam$gratia_diagnostics <- .logistic_normal_diagnostics(fam$rd, m, ref)
    list(y = y / totals, family = fam)
  }

  ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                 d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                 D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
    supplied <- !is.null(eta)
    if (!supplied) eta <- .dirichlet_eta(X, coef, offset)
    logged <- log(y)
    ref <- family$data$refcat
    z <- logged[, -ref, drop = FALSE] - logged[, ref]
    order <- if (deriv == 0L) 0L else if (deriv == 1L) 2L else 3L
    tri <- family$tri
    out <- .logistic_normal_derivatives(z, -rowSums(logged), eta, wt, tri, order)
    result <- out[c("l", "l0")]
    if (!is.finite(out$l)) return(result)
    # The correlation penalty belongs to the fit, outside the row terms
    # used for NCV deletions, held-out scoring and the sandwich score.
    penalize <- family$lkj > 1 && length(correlations) > 0L && !supplied && !sandwich
    if (penalize) {
      penalty <- .logistic_normal_lkj(eta[1L, correlations], m, family$lkj, order)
      result$l <- result$l + penalty$l
    }
    if (!deriv || !is.finite(result$l)) return(result)
    if (!supplied) {
      jj <- attr(X, "lpi")
      gh <- .inlaws_gamlss_gh(
        X, jj, out$l1, out$l2, tri$i2, l3 = out$l3, i3 = tri$i3,
        d1b = d1b, deriv = deriv - 1L, fh = fh, sandwich = sandwich)
      if (penalize)
        gh <- .logistic_normal_penalty_gh(gh, penalty, X, jj[correlations], deriv, d1b, fh)
      result <- c(result, gh)
    }
    if (ncv) result[c("l1", "l2", "l3")] <- list(out$l1, out$l2, if (order == 3L) out$l3 else 0)
    result
  }

  initialize <- expression({
    nobs <- nrow(y)
    n <- rep(1, nobs)
    if (exists("scoreType", inherits = FALSE) && identical(scoreType, "NCV") &&
        !is.null(nei)) {
      if (family$grouped(nei) && isTRUE(nei$jackknife != 0))
        stop("logistic_normal grouped NCV neighbourhoods require jackknife = FALSE", call. = FALSE)
      if (nei$jackknife > 2 && !length(Sl) && utils::packageVersion("mgcv") < "1.9.5")
        stop("unpenalized logistic_normal NCV/QNCV on mgcv 1.9-3 and 1.9-4 requires nei with jackknife = FALSE", call. = FALSE)
    }
    if (is.null(start))
      start <- family$initial(y, x, weights, offset, E, family$data$refcat)
  })

  predict <- function(family, se = FALSE, eta = NULL, y = NULL, X = NULL,
                      beta = NULL, off = NULL, Vb = NULL) {
    if (is.null(eta)) eta <- .dirichlet_eta(X, beta, off) else se <- FALSE
    ref <- family$data$refcat
    keep <- seq_len(D)[-ref]
    center <- .logistic_normal_center(eta[, seq_len(m), drop = FALSE], ref)
    colnames(center) <- family$data$component.names
    if (!se) return(list(fit = center))

    jj <- attr(X, "lpi")
    se.fit <- center
    for (j in seq_len(D)) {
      gradient <- matrix(0, nrow(X), ncol(X))
      for (k in seq_len(m)) {
        d <- center[, j] * ((j == keep[k]) - center[, keep[k]])
        gradient[, jj[[k]]] <- gradient[, jj[[k]], drop = FALSE] +
          d * X[, jj[[k]], drop = FALSE]
      }
      se.fit[, j] <- sqrt(pmax(0, rowSums((gradient %*% Vb) * gradient)))
    }
    list(fit = center, se.fit = se.fit)
  }

  residuals <- function(object, type = c("deviance", "pearson", "response")) {
    .logistic_normal_residuals(object$y, object$fitted.values, object$prior.weights,
                              m, object$family$data$refcat, type)
  }

  postproc <- expression({
    object$deviance <- object$null.deviance <- NA_real_
    object$family$data <- local({
      data <- object$family$data
      m <- ncol(object$y) - 1L
      cor_cols <- seq.int(2L * m + 1L, length.out = choose(m, 2L))
      L <- object$family$chol(object$linear.predictors[1L, cor_cols], m)
      parts <- data$component.names[-data$refcat]
      data$correlation <- tcrossprod(L)
      dimnames(data$correlation) <- list(parts, parts)
      data
    })
    if (object$family$lkj > 1 && ncol(object$y) > 2L) {
      object$l <- object$family$ll(object$y, NULL, NULL, object$prior.weights,
        object$family, eta = object$linear.predictors)$l
      object$aic <- 2 * (sum(object$edf) - object$l)
    }
  })

  rd <- .logistic_normal_rng(m, refcat)
  family <- structure(list(
    family = "logistic_normal", nlp = q, link = rep("identity", q),
    linfo = replicate(q, stats::make.link("identity"), simplify = FALSE),
    ll = ll, tri = .dirmult_tri(q), ncv = .inlaws_gamlss_ncv,
    preinitialize = preinitialize, initialize = initialize,
    initial = .logistic_normal_initial, postproc = postproc, predict = predict,
    presetup = expression({
      if (family$grouped(nei)) {
        if (isTRUE(nei$jackknife > 0))
          stop("logistic_normal grouped NCV neighbourhoods require jackknife = FALSE", call. = FALSE)
        nei$jackknife <- FALSE
      }
    }),
    grouped = function(nei) !is.null(nei$ma) &&
      (length(nei$a) != length(nei$ma) ||
        (!is.null(nei$d) && !isTRUE(all.equal(nei$a, nei$d)))),
    rd = rd, residuals = residuals, chol = .logistic_normal_chol,
    gratia_diagnostics = .logistic_normal_diagnostics(rd, m, refcat),
    sandwich = function(y, X, coef, wt, family, offset = NULL)
      ll(y, X, coef, wt, family, offset, deriv = 1L, sandwich = TRUE)$lbb,
    data = list(refcat = refcat), component.count = D, lkj = lkj,
    available.derivs = 1L, d2link = 1, d3link = 1, d4link = 1, ls = 1,
    discrete.ok = FALSE, no.r.sq = TRUE),
    class = c("general.family", "extended.family", "family"))
  .inlaws_ncv_family(family)
}

.logistic_normal_penalty_gh <- function(gh, penalty, X, jj, deriv, d1b, fh) {
  keep <- lengths(jj) > 0L
  columns <- unlist(jj, use.names = FALSE)
  if (!length(columns)) return(gh)
  x <- X[1L, columns]
  xx <- tcrossprod(x)
  gh$lb[columns] <- gh$lb[columns] + x * penalty$gradient[keep]
  gh$lbb[columns, columns] <- gh$lbb[columns, columns, drop = FALSE] +
    xx * penalty$hessian[keep, keep, drop = FALSE]
  if (deriv > 1L) {
    third <- matrix(penalty$third[keep, keep, keep, drop = FALSE], nrow = length(xx))
    directions <- x * d1b[columns, , drop = FALSE]
    if (deriv == 2L) {
      # REML BFGS needs only the trace against the inverse penalized Hessian.
      trace_weight <- as.vector(xx * fh[columns, columns, drop = FALSE])
      gh$d1H <- gh$d1H + drop(crossprod(trace_weight, third) %*% directions)
    } else {
      # NCV needs each full coefficient Hessian derivative for its row downdates.
      blocks <- third %*% directions
      for (j in seq_along(gh$d1H)) {
        gh$d1H[[j]][columns, columns] <- gh$d1H[[j]][columns, columns, drop = FALSE] +
          xx * matrix(blocks[, j], nrow = length(columns))
      }
    }
  }
  gh
}

#' @rdname logistic_normal
#' @export
logistic_normal_formula <- function(formula, family, scale = ~ 1) {
  if (!inherits(formula, "formula") || length(formula) != 3L)
    stop("formula must be a two-sided formula with a matrix response", call. = FALSE)
  if (!inherits(scale, "formula") || length(scale) != 2L)
    stop("scale must be a one-sided formula", call. = FALSE)
  if (!inherits(family, "general.family") || !identical(family$family, "logistic_normal"))
    stop("family must be a logistic_normal family", call. = FALSE)

  m <- family$component.count - 1L
  mean_rhs <- stats::as.formula(as.call(list(as.name("~"), formula[[3L]])),
                                env = environment(formula))
  c(list(formula), rep(list(mean_rhs), m - 1L), rep(list(scale), m),
    rep(list(~ 1), choose(m, 2L)))
}

.logistic_normal_initial <- function(y, X, wt, offset, E, ref) {
  n <- nrow(y)
  m <- ncol(y) - 1L
  jj <- attr(X, "lpi")
  q <- length(jj)
  if (is.null(wt)) wt <- rep(1, n)
  off <- matrix(0, n, q)
  for (j in seq_len(q)) {
    if (length(offset) >= j && !is.null(offset[[j]])) off[, j] <- offset[[j]]
  }

  means <- seq_len(m)
  z <- log(y[, -ref, drop = FALSE]) - log(y[, ref]) - off[, means, drop = FALSE]
  mu <- colSums(z * wt) / sum(wt)
  r <- sweep(z, 2L, mu)
  empirical <- crossprod(r * sqrt(wt / sum(wt)))
  sigma <- sqrt(pmax(diag(empirical), .Machine$double.eps))
  rho <- empirical / outer(sigma, sigma)
  diag(rho) <- 1
  # An interior starting point is useful even for nearly singular residuals.
  L <- t(chol(0.9 * rho + 0.1 * diag(m)))
  raw <- L / diag(L)
  theta <- t(raw)[upper.tri(t(raw), diag = FALSE)]
  target <- c(mu, log(sigma), theta)
  off[, means] <- 0

  design <- matrix(0, n * q, ncol(X))
  rhs <- numeric(n * q)
  for (j in seq_len(q)) {
    rows <- (j - 1L) * n + seq_len(n)
    design[rows, jj[[j]]] <- X[, jj[[j]], drop = FALSE] * sqrt(wt)
    rhs[rows] <- (target[j] - off[, j]) * sqrt(wt)
  }
  start <- qr.coef(qr(rbind(design, E)), c(rhs, numeric(nrow(E))))
  start[!is.finite(start)] <- 0
  start
}

.logistic_normal_center <- function(mu, ref) {
  logits <- matrix(0, nrow(mu), ncol(mu) + 1L)
  logits[, -ref] <- mu
  p <- exp(logits - apply(logits, 1L, max))
  p / rowSums(p)
}

.logistic_normal_rng <- function(m, ref, component_names = NULL) {
  force(m); force(ref); force(component_names)
  function(mu, wt, scale) {
    if (is.character(ref))
      stop("a named refcat is resolved when logistic_normal is fitted", call. = FALSE)
    means <- seq_len(m)
    theta <- mu[, -seq_len(2L * m), drop = FALSE]
    noise <- matrix(stats::rnorm(nrow(mu) * m), ncol = m)
    for (i in seq_len(nrow(mu)))
      noise[i, ] <- drop(.logistic_normal_chol(theta[i, ], m) %*% noise[i, ])
    z <- mu[, means, drop = FALSE] + exp(mu[, m + means, drop = FALSE]) * noise
    y <- .logistic_normal_center(z, ref)
    colnames(y) <- component_names
    y
  }
}

.logistic_normal_residuals <- function(y, eta, wt, m, ref, type) {
  type <- match.arg(type, c("deviance", "pearson", "response"))
  if (type == "deviance")
    stop("logistic_normal deviance residuals are undefined; use type = 'pearson' or 'response'", call. = FALSE)
  means <- seq_len(m)
  mu <- eta[, means, drop = FALSE]
  if (type == "response") return(y - .logistic_normal_center(mu, ref))
  r <- (log(y[, -ref, drop = FALSE]) - log(y[, ref]) - mu) *
    exp(-eta[, m + means, drop = FALSE])
  theta <- eta[, -seq_len(2L * m), drop = FALSE]
  ans <- numeric(nrow(y))
  for (i in which(wt > 0)) {
    z <- forwardsolve(.logistic_normal_chol(theta[i, ], m), r[i, ])
    ans[i] <- sqrt(wt[i] * sum(z^2))
  }
  ans
}

.logistic_normal_diagnostics <- function(rd, m, ref) {
  force(rd); force(m); force(ref)
  function(mu, wt, scale, y, eta) {
    if (is.null(wt)) wt <- rep(1, nrow(mu))
    list(mean = matrix(NA_real_, nrow(mu), m + 1L),
         mean_available = FALSE, variance_available = FALSE,
         moment_reason = "Estimate response moments by averaging jointly simulated compositions.",
         independent = FALSE, continuous = TRUE,
         rd = function() rd(mu, wt, scale),
         diagnostic_residuals = function(y, type) {
           if (type != "response")
             stop("Use residuals(fit, type = 'pearson') for logistic-normal row-level distances", call. = FALSE)
           .logistic_normal_residuals(y, mu, wt, m, ref, type)
         })
  }
}
