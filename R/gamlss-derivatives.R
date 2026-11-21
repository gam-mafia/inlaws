.inlaws_gamlss_gh <- function(X, jj, l1, l2, i2, l3 = 0, i3 = 0,
                             l4 = 0, i4 = 0, d1b = 0, d2b = 0, deriv = 0,
                             fh = NULL, D = NULL, sandwich = FALSE) {
  singleton <- deriv == 1L && any(lengths(jj) == 1L)
  shared <- deriv == 3L && anyDuplicated(unlist(jj)) > 0L
  # mgcv drops dimensions in the first trace for singleton predictor blocks.
  # Its fourth-order trace overwrites contributions from shared coefficients.
  native_deriv <- if (singleton) 0L else if (shared) 2L else deriv
  ret <- mgcv::gamlss.gH(X, jj, l1, l2, i2, l3 = l3, i3 = i3, l4 = l4, i4 = i4,
    d1b = d1b, d2b = d2b, deriv = native_deriv, fh = fh, D = D,
    sandwich = sandwich)
  if (singleton) ret$d1H <- .inlaws_first_trace(X, jj, l3, i3, d1b, fh)
  if (shared) ret$trHid2H <- .inlaws_shared_trace(X, jj, list(l3 = l3, l4 = l4),
    list(i3 = i3, i4 = i4), d1b, d2b, fh, D)
  ret
}

.inlaws_first_trace <- function(X, jj, l3, i3, d1b, V) {
  design <- lapply(jj, function(j) X[, j, drop = FALSE])
  d_eta <- lapply(seq_along(jj), function(j)
    design[[j]] %*% d1b[jj[[j]], , drop = FALSE])
  trace <- numeric(ncol(d1b))
  for (i in seq_along(jj)) for (j in i:length(jj)) {
    h <- rowSums((design[[i]] %*% V[jj[[i]], jj[[j]], drop = FALSE]) * design[[j]])
    mult <- if (i == j) 1 else 2
    for (k in seq_along(jj)) {
      trace <- trace + mult * colSums((h * l3[, i3[i, j, k]]) * d_eta[[k]])
    }
  }
  trace
}

.inlaws_shared_trace <- function(X, jj, de, tri, d1b, d2b, fh, D) {
  if (is.null(D)) D <- rep(1, ncol(X))
  if (is.list(fh)) {
    ev <- fh$values
    ev <- ifelse(ev > 0, 1 / ev, 0)
    V <- (fh$vectors %*% (ev * t(fh$vectors))) * outer(D, D)
  } else {
    pivot <- attr(fh, "pivot")
    if (is.null(pivot)) pivot <- seq_len(ncol(X))
    V <- chol2inv(fh)[order(pivot), order(pivot), drop = FALSE] * outer(D, D)
  }
  x <- lapply(jj, function(j) X[, j, drop = FALSE])
  a <- lapply(seq_along(jj), function(j) x[[j]] %*% d1b[jj[[j]], , drop = FALSE])
  b <- lapply(seq_along(jj), function(j) x[[j]] %*% d2b[jj[[j]], , drop = FALSE])
  out <- numeric(ncol(d2b))
  for (i in seq_along(jj)) for (j in i:length(jj)) {
    h <- rowSums((x[[i]] %*% V[jj[[i]], jj[[j]], drop = FALSE]) * x[[j]])
    kk <- 0L
    for (k in seq_len(ncol(d1b))) for (l in k:ncol(d1b)) {
      kk <- kk + 1L
      v <- numeric(nrow(X))
      for (q in seq_along(jj)) {
        v <- v + de$l3[, tri$i3[i, j, q]] * b[[q]][, kk]
        for (s in seq_along(jj)) v <- v + de$l4[, tri$i4[i, j, q, s]] *
          a[[q]][, k] * a[[s]][, l]
      }
      out[kk] <- out[kk] + (if (i == j) 1 else 2) * sum(h * v)
    }
  }
  out
}
