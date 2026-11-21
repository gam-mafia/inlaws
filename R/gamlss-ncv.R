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
  # Cholesky downdates can disagree with the full signed Hessian. Use the
  # native conjugate-gradient solver for predictions and smoothing derivatives.
  if (!is.null(R)) {
    H <- crossprod(R)
    Hi <- chol2inv(R)
  }
  updates <- if (deriv < 0) {
    .inlaws_gamlss_jackknife(X, nei, beta, llf, H, off)
  } else {
    utils::getFromNamespace("gamlss.ncv", "mgcv")(
      X, seq_len(n), wt, nei, beta, native_family, native_ll, H = H, Hi = Hi,
      R = NULL, offset = off, dH = dH, db = db, deriv = deriv, nt = nt)
  }
  cg_iterations <- updates$error
  if (any(!is.finite(cg_iterations) | cg_iterations >= 200L))
    stop(family$family, " NCV deletion solve failed to converge")
  held_out <- attr(updates$NCV, "eta.cv")
  d_held_out <- attr(updates$NCV, "deta.cv")
  if (any(!is.finite(held_out)) || (deriv != 0 && any(!is.finite(d_held_out))))
    stop(family$family, " NCV deletion solve returned non-finite predictions or derivatives")
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
  if (is.nan(score) || any(!is.finite(grad)))
    stop(family$family, " NCV criterion or gradient is non-finite")
  attr(score, "eta.cv") <- held_out
  if (deriv != 0) attr(score, "deta.cv") <- d_held_out
  # gam.fit5 expects a count of failed Cholesky downdates, not CG iterations.
  list(NCV = score, NCV1 = grad, error = 0L, Vg = Vg, cg.iterations = cg_iterations)
}

# The native CG helper has an uninitialized negative-deriv buffer in mgcv 1.9-4.
# Compute the coefficient deletion updates directly for jackknife output.
.inlaws_gamlss_jackknife <- function(X, nei, beta, llf, H, offset) {
  jj <- attr(X, "lpi")
  predictors <- seq_along(jj)
  eta <- matrix(0, length(nei$d), length(jj))
  delta <- matrix(0, length(nei$ma), length(beta))
  previous_a <- previous_d <- 0L
  for (i in seq_along(nei$ma)) {
    rows <- nei$a[previous_a + seq_len(nei$ma[i] - previous_a)]
    ix <- previous_d + seq_len(nei$md[i] - previous_d)
    pred <- nei$d[ix]
    held_hessian <- H
    score <- numeric(length(beta))
    column <- 0L
    for (j in predictors) {
      xj <- X[rows, jj[[j]], drop = FALSE]
      score[jj[[j]]] <- score[jj[[j]]] + drop(crossprod(xj, llf$l1[rows, j]))
      for (k in j:length(jj)) {
        column <- column + 1L
        block <- crossprod(xj, X[rows, jj[[k]], drop = FALSE] * llf$l2[rows, column])
        held_hessian[jj[[j]], jj[[k]]] <- held_hessian[jj[[j]], jj[[k]], drop = FALSE] + block
        if (j != k) held_hessian[jj[[k]], jj[[j]]] <-
          held_hessian[jj[[k]], jj[[j]], drop = FALSE] + t(block)
      }
    }
    delta[i, ] <- solve(held_hessian, score)
    for (j in predictors) {
      eta[ix, j] <- drop(X[pred, jj[[j]], drop = FALSE] %*%
        (beta[jj[[j]]] - delta[i, jj[[j]]])) + offset[[j]][pred]
    }
    previous_a <- nei$ma[i]
    previous_d <- nei$md[i]
  }
  list(NCV = structure(0, eta.cv = eta, deta.cv = delta), error = 0L)
}
