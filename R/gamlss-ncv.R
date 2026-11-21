.inlaws_gamlss_ncv <- function(X, y, wt, nei, beta, family, llf, H = NULL,
                              Hi = NULL, R = NULL, offset = NULL, dH = NULL,
                              db = NULL, deriv = FALSE, nt = 1,
                              subset_family = NULL) {
  n <- nrow(X)
  jj <- attr(X, "lpi")
  predictors <- seq_along(jj)
  rows <- nei$d
  m <- length(rows)
  wt <- rep_len(if (is.null(wt)) 1 else wt, n)
  off <- lapply(predictors, function(j) {
    if (length(offset) < j || is.null(offset[[j]])) return(numeric(n))
    rep_len(offset[[j]], n)
  })
  # Use mgcv's compiled deletion updates. Compute the score and gradient
  # locally so the full and held-out predictors use the same offsets and rows.
  native_family <- family
  native_family$qapprox <- FALSE
  native_family$ll <- function(y, ..., eta) list(l = 0, l1 = eta * 0)
  native_ll <- llf
  native_ll$gamma <- 1
  updates <- utils::getFromNamespace("gamlss.ncv", "mgcv")(
    X, seq_len(n), wt, nei, beta, native_family, native_ll, H = H, Hi = Hi,
    R = R, offset = off, dH = dH, db = db, deriv = deriv, nt = nt)
  held_out <- attr(updates$NCV, "eta.cv")
  d_held_out <- attr(updates$NCV, "deta.cv")
  gamma <- if (is.null(llf$gamma)) 1 else llf$gamma
  quadratic <- isTRUE(family$qapprox)
  if (quadratic) {
    eta <- vapply(predictors, function(j)
      drop(X[rows, jj[[j]], drop = FALSE] %*% beta[jj[[j]]]) + off[[j]][rows], numeric(m))
    dim(eta) <- c(m, length(predictors))
    delta <- held_out - eta
    h_delta <- matrix(0, m, length(predictors))
    for (j in predictors) for (k in predictors) {
      h_delta[, j] <- h_delta[, j] + llf$l2[rows, family$tri$i2[j, k]] * delta[, k]
    }
    score <- -sum(llf$l0[rows]) -
      gamma * sum((llf$l1[rows, , drop = FALSE] + h_delta / 2) * delta)
  } else {
    response <- if (is.matrix(y)) y[rows, , drop = FALSE] else y[rows]
    bounds <- attr(y, "censor")
    if (!is.null(bounds)) attr(response, "censor") <- bounds[rows]
    held_family <- if (is.null(subset_family)) family else subset_family(family, rows)
    cv <- family$ll(response, X[rows, , drop = FALSE], beta, wt[rows], held_family,
      eta = held_out, deriv = as.integer(deriv > 0), ncv = TRUE)
    score <- -gamma * cv$l
    if (gamma != 1) score <- score - (1 - gamma) * sum(llf$l0[rows])
  }
  grad <- Vg <- NULL
  if (deriv > 0) {
    contributions <- matrix(0, m, ncol(db))
    # A deletion can put a held-out response outside parameter-dependent
    # support. The infinite score rejects the trial without an NA gradient.
    if (is.finite(score)) for (j in predictors) {
      dc <- d_held_out[seq_len(m) + (j - 1L) * m, , drop = FALSE]
      d0 <- X[rows, jj[[j]], drop = FALSE] %*% db[jj[[j]], , drop = FALSE]
      contributions <- contributions - (1 - gamma) * llf$l1[rows, j] * d0
      if (quadratic) {
        contributions <- contributions - gamma * (llf$l1[rows, j] + h_delta[, j]) * dc
        third <- numeric(m)
        for (k in predictors) for (l in predictors) {
          third <- third + llf$l3[rows, family$tri$i3[j, k, l]] * delta[, k] * delta[, l]
        }
        contributions <- contributions - (gamma / 2) * third * d0
      } else {
        contributions <- contributions - gamma * cv$l1[, j] * dc
      }
    }
    grad <- colSums(contributions)
    Vg <- crossprod(contributions)
  }
  attr(score, "eta.cv") <- held_out
  if (deriv != 0) attr(score, "deta.cv") <- d_held_out
  list(NCV = score, NCV1 = grad, error = updates$error, Vg = Vg)
}
