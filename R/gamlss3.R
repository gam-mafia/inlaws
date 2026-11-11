.inlaws_gamlss3_eta <- function(X, beta, offset, family) {
  if (!is.matrix(X)) stop(family, " requires ordinary model matrices")
  jj <- attr(X, "lpi")
  if (length(jj) != 3L) stop(family, " requires three formulae")
  if (anyDuplicated(unlist(jj))) stop(family, " does not support shared predictor coefficients")
  eta <- matrix(0, nrow(X), 3L)
  for (k in 1:3) {
    eta[, k] <- drop(X[, jj[[k]], drop = FALSE] %*% beta[jj[[k]]])
    if (length(offset) >= k && !is.null(offset[[k]])) eta[, k] <- eta[, k] + offset[[k]]
  }
  eta
}

# Assemble the mgcv likelihood from weighted multivariate Taylor coefficients.
.inlaws_gamlss_jet_ll <- function(ell, powers, X, jj, family, supplied, deriv,
                                 d1b, d2b, fh, D, ncv, sandwich) {
  ret <- list()
  if (deriv) {
    order <- max(rowSums(powers))
    # mgcv packs symmetric tensors lexicographically, not expand.grid order.
    packed <- lapply(seq_len(order), function(d) {
      ix <- which(rowSums(powers) == d)
      p <- powers[ix, , drop = FALSE]
      i <- base::order(-p[, 1L], -p[, 2L])
      factorials <- apply(factorial(p[i, , drop = FALSE]), 1L, prod)
      sweep(ell[, ix[i], drop = FALSE], 2L, factorials, `*`)
    })
    l3 <- if (order >= 3L) packed[[3]] else 0
    l4 <- if (order >= 4L) packed[[4]] else 0
    contracted <- deriv == 2L && any(lengths(jj) == 1L)
    # mgcv 1.9-4 drops matrix dimensions in its trace-only route
    # when a predictor has a single coefficient. Contract the full Hessian
    # derivative instead, preserving intercept-only scale/shape predictors.
    gh_deriv <- if (contracted) 2L else deriv - 1L
    if (!supplied && !is.null(jj)) ret <- mgcv::gamlss.gH(
      X, jj, packed[[1]], packed[[2]], family$tri$i2,
      l3 = l3, i3 = family$tri$i3, l4 = l4, i4 = family$tri$i4,
      d1b = d1b, d2b = d2b, deriv = gh_deriv, fh = fh, D = D, sandwich = sandwich)
    if (contracted && !supplied) ret$d1H <- vapply(ret$d1H, function(z) sum(z*fh), numeric(1))
    if (ncv) { ret$l1 <- packed[[1]]; ret$l2 <- packed[[2]]; ret$l3 <- l3
      if (order >= 4L) ret$l4 <- l4 }
  }
  ret$l0 <- ell[, 1]; ret$l <- sum(ret$l0); ret
}

# Three-predictor NCV/QNCV: compiled deletion updates, local score assembly
# for correct offsets, prediction subsets/reordering and gamma adjustment.
.inlaws_gamlss3_ncv <- function(X, y, wt, nei, beta, family, llf, H = NULL, Hi = NULL, R = NULL,
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
  # mgcv 1.9-4 Cholesky deletion updates can disagree with the actual
  # deleted Hessian for three predictors with negative cross-curvature.
  # Its conjugate-gradient route uses the full signed Hessian instead.
  if (!is.null(R)) { H <- crossprod(R); Hi <- chol2inv(R) }
  R <- NULL
  ret <- if (deriv < 0) .inlaws_gamlss3_jackknife(X, nei, beta, lf, H, off) else
    native(X, y, wt, nei, beta, f, lf, H = H, Hi = Hi, R = NULL,
           offset = off, dH = dH, db = db, deriv = deriv, nt = nt)
  cg_iterations <- ret$error
  if (any(!is.finite(cg_iterations) | cg_iterations >= 200L)) stop(family$family, " NCV deletion solve failed to converge")
  # The CG helper returns iteration counts, whereas gam.fit5 expects a
  # scalar count of failed Cholesky downdates. Keep these distinct.
  ret$error <- 0L
  ecv <- attr(ret$NCV, "eta.cv")
  dcv <- attr(ret$NCV, "deta.cv")
  if (any(!is.finite(ecv)) || (deriv != 0 && any(!is.finite(dcv))))
    stop(family$family, " NCV deletion solve returned non-finite predictions or derivatives")
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
  if (!is.finite(score) || any(!is.finite(grad)))
    stop(family$family, " NCV criterion or gradient is non-finite")
  attr(score, "eta.cv") <- ecv
  if (deriv != 0) attr(score, "deta.cv") <- dcv
  list(NCV = score, NCV1 = grad, error = ret$error, Vg = Vg, cg.iterations = cg_iterations)
}


# mgcv's CG helper in 1.9-4 does not initialize its negative-deriv output
# buffer. Assemble coefficient deletion updates directly for jackknife output.
.inlaws_gamlss3_jackknife <- function(X, nei, beta, llf, H, offset) {
  jj <- attr(X, "lpi"); p <- length(beta)
  ecv <- matrix(0, length(nei$d), 3L)
  delta <- matrix(0, length(nei$ma), p)
  previous_a <- previous_d <- 0L
  for (i in seq_along(nei$ma)) {
    a <- nei$a[previous_a + seq_len(nei$ma[i]-previous_a)]
    ix <- previous_d + seq_len(nei$md[i]-previous_d); pred <- nei$d[ix]
    hd <- H; score <- numeric(p); col <- 0L
    for (j in 1:3) {
      xj <- X[a, jj[[j]], drop = FALSE]
      score[jj[[j]]] <- drop(crossprod(xj, llf$l1[a, j]))
      for (k in j:3) {
        col <- col + 1L
        block <- crossprod(xj, X[a, jj[[k]], drop = FALSE]*llf$l2[a, col])
        hd[jj[[j]], jj[[k]]] <- hd[jj[[j]], jj[[k]], drop = FALSE] + block
        if (j != k) hd[jj[[k]], jj[[j]]] <- hd[jj[[k]], jj[[j]], drop = FALSE] + t(block)
      }
    }
    delta[i, ] <- solve(hd, score)
    for (j in 1:3) ecv[ix, j] <- drop(X[pred, jj[[j]], drop = FALSE] %*%
      (beta[jj[[j]]] - delta[i, jj[[j]]])) + offset[[j]][pred]
    previous_a <- nei$ma[i]; previous_d <- nei$md[i]
  }
  score <- 0
  attr(score, "eta.cv") <- ecv; attr(score, "deta.cv") <- delta
  list(NCV = score, error = 0L)
}
