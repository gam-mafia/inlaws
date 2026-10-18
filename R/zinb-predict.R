# Native general.family prediction: eta is supplied for values alone; mgcv
# supplies X, its lpi blocks, offsets and the selected coefficient covariance
# for predictions with uncertainty. A shared coefficient must accumulate its
# contributions from every predictor before the covariance contraction.
.zinb_predict <- function(family, se = FALSE, eta = NULL, y = NULL, X = NULL,
                          beta = NULL, off = NULL, Vb = NULL) {
  if (is.null(eta)) {
    jj <- attr(X, "lpi")
    eta <- matrix(0, nrow(X), 3L)
    for (j in 1:3) {
      eta[, j] <- drop(X[, jj[[j]], drop = FALSE] %*% beta[jj[[j]]])
      if (length(off) >= j && !is.null(off[[j]]))
        eta[, j] <- eta[, j] + off[[j]]
    }
  } else se <- FALSE
  z <- .zinb_prediction(eta, family$hurdle)
  if (!se) return(list(fit = z$fit))
  ans <- z$fit * 0
  for (k in seq_len(ncol(ans))) {
    J <- matrix(0, nrow(X), ncol(X))
    for (j in 1:3) J[, jj[[j]]] <- J[, jj[[j]], drop = FALSE] +
      X[, jj[[j]], drop = FALSE] * z$gradient[, k, j]
    ans[, k] <- sqrt(pmax(0, rowSums((J %*% Vb) * J)))
  }
  list(fit = z$fit, se.fit = ans)
}

# Complete response moments and analytic derivatives on the predictor scales.
.zinb_prediction <- function(eta, hurdle) {
  m <- exp(eta[, 1]); p <- stats::plogis(eta[, 2])
  absent <- stats::plogis(-eta[, 2]); dp <- p * absent
  t <- exp(eta[, 3]); r <- m / t
  a <- t * log1p(r)
  # Preserve a ~ mu even if mu/theta underflows to zero.
  a[!is.na(r) & r == 0] <- m[!is.na(r) & r == 0]
  q0 <- exp(-a); positive <- -expm1(-a)
  au <- m / (1 + r)
  # D = log1p(r)/r - 1/(1+r). Its series avoids cancellation in
  # d a / d log(theta) = mu * D near the Poisson limit.
  D <- log1p(r) / r - 1 / (1 + r)
  small <- !is.na(r) & r < .01
  if (any(small)) {
    rr <- r[small]; dd <- rr * 0
    for (k in 12:1) dd <- rr * ((-1)^(k + 1L) * k / (k + 1L) + dd)
    D[small] <- dd
  }
  av <- m * D
  n <- nrow(eta)
  dm <- dv <- matrix(0, n, 3L)
  if (hurdle) {
    component <- m / positive
    excess <- component - 1
    dm[, 1] <- component * (1 - au * q0 / positive)
    dm[, 3] <- -component * av * q0 / positive
    # At tiny means the positive component tends to a point mass at one.
    # Evaluate its excess mean and derivatives without subtracting ones:
    # mu/(1-exp(-a)) = {1/L(r)} * {a/(1-exp(-a))}, L=log1p(r)/r.
    tiny <- small & !is.na(a) & a < .01
    if (any(tiny)) {
      rr <- r[tiny]; aa <- a[tiny]; lm1 <- rr * 0
      for (k in 12:1) lm1 <- rr * ((-1)^k / (k + 1L) + lm1)
      A <- -lm1 / (1 + lm1) # 1/L - 1
      Hm1 <- aa * (1/2 + aa * (1/12 + aa^2 * (-1/720 + aa^2/30240)))
      Hp <- 1/2 + aa * (1/6 + aa^2 * (-1/180 + aa^2/5040))
      excess[tiny] <- A + Hm1 + A * Hm1
      component[tiny] <- 1 + excess[tiny]
      dA <- D[tiny] / (1 + lm1)^2
      dm[tiny, 1] <- dA * (1 + Hm1) + (1 + A) * Hp * au[tiny]
      dm[tiny, 3] <- -dA * (1 + Hm1) + (1 + A) * Hp * av[tiny]
    }
    cv <- pmax(0, component * (m + r - excess))
    dv[, 1] <- dm[, 1] * (m + r - 1 - 2 * excess) + component * (m + r)
    dv[, 3] <- dm[, 3] * (m + r - 1 - 2 * excess) - component * r
    p0 <- absent
  } else {
    component <- m; cv <- m * (1 + r)
    dm[, 1] <- m
    dv[, 1] <- m * (1 + 2 * r); dv[, 3] <- -m * r
    p0 <- absent + p * q0
  }
  mn <- p * component
  vr <- p * cv + p * absent * component^2
  fit <- cbind(mu = m, p = p, theta = t, mean = mn, variance = vr, p0 = p0)
  g <- array(0, c(n, 6L, 3L))
  g[, 1, 1] <- m; g[, 2, 2] <- dp; g[, 3, 3] <- t
  for (j in 1:3) {
    dj <- if (j == 2L) dp else rep(0, n)
    g[, 4, j] <- p * dm[, j] + dj * component
    g[, 5, j] <- p * dv[, j] + dj * cv +
      dj * (absent - p) * component^2 + 2 * p * absent * component * dm[, j]
  }
  g[, 6, 2] <- -dp * if (hurdle) 1 else positive
  if (!hurdle) {
    g[, 6, 1] <- -p * q0 * au
    g[, 6, 3] <- -p * q0 * av
  }
  list(fit = fit, gradient = g)
}
