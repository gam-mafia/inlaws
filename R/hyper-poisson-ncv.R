# Two-predictor signed-Hessian NCV/QNCV, following scatlss().
.hp_ncv <- function(X, y, wt, nei, beta, family, llf, H = NULL, Hi = NULL, R = NULL,
                      offset = NULL, dH = NULL, db = NULL, deriv = FALSE, nt = 1) {
  if (is.null(wt)) wt <- rep(1, length(y))
  n <- length(y)
  jj <- attr(X, "lpi")
  m <- length(nei$d)
  off <- lapply(1:2, function(j) {
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
  # mgcv Cholesky deletion updates can disagree with the actual
  # deleted Hessian when mixed curvature has either sign.
  # Its conjugate-gradient route uses the full signed Hessian instead.
  if (!is.null(R)) { H <- crossprod(R); Hi <- chol2inv(R) }
  R <- NULL
  ret <- if (deriv < 0) .hp_jackknife(X, nei, beta, lf, H, off) else
    native(X, y, wt, nei, beta, f, lf, H = H, Hi = Hi, R = NULL,
           offset = off, dH = dH, db = db, deriv = deriv, nt = nt)
  cg_iterations <- ret$error
  if (any(!is.finite(cg_iterations) | cg_iterations >= 200L)) stop("hyper_poissonls NCV deletion solve failed to converge")
  # The CG helper returns iteration counts, whereas gam.fit5 expects a
  # scalar count of failed Cholesky downdates. Keep these distinct.
  ret$error <- 0L
  ecv <- attr(ret$NCV, "eta.cv")
  dcv <- attr(ret$NCV, "deta.cv")
  if (any(!is.finite(ecv)) || (deriv != 0 && any(!is.finite(dcv))))
    stop("hyper_poissonls NCV deletion solve returned non-finite predictions or derivatives")
  eta <- sapply(1:2, function(j) drop(X[, jj[[j]], drop = FALSE] %*% beta[jj[[j]]]) + off[[j]])
  base <- family$ll(y, X, beta, wt, family, eta = eta, deriv = if (deriv > 0) 3 else 1, ncv = TRUE)
  ix <- nei$d
  gamma <- llf$gamma
  if (is.null(gamma)) gamma <- 1
  score0 <- -sum(base$l0[ix])
  delta <- ecv - eta[ix, , drop = FALSE]
  tri <- family$tri
  if (isTRUE(family$qapprox)) {
    Hd <- matrix(0, m, 2)
    for (j in 1:2) {
      for (k in 1:2) {
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
    d0 <- lapply(1:2, function(j) X[ix, jj[[j]], drop = FALSE] %*% db[jj[[j]], , drop = FALSE])
    for (j in 1:2) {
      dc <- dcv[seq_len(m) + (j - 1L) * m, , drop = FALSE]
      g <- g - (1 - gamma) * base$l1[ix, j] * d0[[j]]
      if (isTRUE(family$qapprox)) {
        g <- g - gamma * (base$l1[ix, j] + Hd[, j]) * dc
        td <- numeric(m)
        for (k in 1:2) {
          for (l in 1:2) {
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
    stop("hyper_poissonls NCV criterion or gradient is non-finite")
  attr(score, "eta.cv") <- ecv
  if (deriv != 0) attr(score, "deta.cv") <- dcv
  list(NCV = score, NCV1 = grad, error = ret$error, Vg = Vg, cg.iterations = cg_iterations)
}


# mgcv's CG helper in 1.9-4 does not initialize its negative-deriv output
# buffer. Assemble coefficient deletion updates directly for jackknife output.
.hp_jackknife <- function(X, nei, beta, llf, H, offset) {
  jj <- attr(X, "lpi"); p <- length(beta)
  ecv <- matrix(0, length(nei$d), 2L)
  delta <- matrix(0, length(nei$ma), p)
  previous_a <- previous_d <- 0L
  for (i in seq_along(nei$ma)) {
    a <- nei$a[previous_a + seq_len(nei$ma[i]-previous_a)]
    ix <- previous_d + seq_len(nei$md[i]-previous_d); pred <- nei$d[ix]
    hd <- H; score <- numeric(p); col <- 0L
    for (j in 1:2) {
      xj <- X[a, jj[[j]], drop = FALSE]
      score[jj[[j]]] <- score[jj[[j]]] + drop(crossprod(xj, llf$l1[a, j]))
      for (k in j:2) {
        col <- col + 1L
        block <- crossprod(xj, X[a, jj[[k]], drop = FALSE]*llf$l2[a, col])
        hd[jj[[j]], jj[[k]]] <- hd[jj[[j]], jj[[k]], drop = FALSE] + block
        if (j != k) hd[jj[[k]], jj[[j]]] <- hd[jj[[k]], jj[[j]], drop = FALSE] + t(block)
      }
    }
    delta[i, ] <- solve(hd, score)
    for (j in 1:2) ecv[ix, j] <- drop(X[pred, jj[[j]], drop = FALSE] %*%
      (beta[jj[[j]]] - delta[i, jj[[j]]])) + offset[[j]][pred]
    previous_a <- nei$ma[i]; previous_d <- nei$md[i]
  }
  score <- 0
  attr(score, "eta.cv") <- ecv; attr(score, "deta.cv") <- delta
  list(NCV = score, error = 0L)
}
