# Free entries in row order in a unit lower-triangular matrix. Dividing
# each row by its norm gives the Cholesky factor of the correlation matrix.
.logistic_normal_chol <- function(theta, m) {
  out <- diag(m)
  at <- 0L
  if (m > 1L) for (i in 2L:m) {
    ix <- seq_len(i - 1L)
    values <- c(theta[at + ix], 1)
    values <- values / max(abs(values))
    out[i, seq_len(i)] <- values / sqrt(sum(values^2))
    at <- at + i - 1L
  }
  out
}

.logistic_normal_row_norm <- function(theta) {
  scale <- do.call(pmax, c(list(1), unname(asplit(abs(theta), 2L))))
  scaled <- theta / scale
  squares <- rowSums(scaled^2)
  norm <- sqrt((1 / scale)^2 + squares)
  log_norm <- log(scale) + log(norm)
  unit_scale <- scale == 1
  log_norm[unit_scale] <- log1p(squares[unit_scale]) / 2
  list(log_norm = log_norm, normalized = scaled / norm)
}

.logistic_normal_lkj <- function(theta, m, lkj, order = 0L) {
  k <- length(theta)
  out <- list(l = 0)
  if (order >= 1L) out$gradient <- numeric(k)
  if (order >= 2L) out$hessian <- matrix(0, k, k)
  if (order >= 3L) out$third <- array(0, rep.int(k, 3L))
  if (k == 0L || lkj == 1) return(out)

  # Each diagonal of the correlation Cholesky factor is the reciprocal
  # of a row norm. This is a penalty on the correlation matrix itself.
  weight <- lkj - 1
  at <- 0L
  for (i in seq_len(m - 1L)) {
    ix <- at + seq_len(i)
    row <- .logistic_normal_row_norm(matrix(theta[ix], 1L))
    out$l <- out$l + weight * (-2 * row$log_norm)
    if (order >= 1L) {
      u <- drop(row$normalized)
      inv <- exp(-row$log_norm)
      weighted_inv <- weight * inv
      out$gradient[ix] <- weighted_inv * (-2 * u)
    }
    if (order >= 2L) {
      identity <- diag(i)
      product <- tcrossprod(u)
      weighted_inv2 <- weighted_inv * inv
      out$hessian[ix, ix] <- weighted_inv2 * (-2 * (identity - 2 * product))
    }
    if (order >= 3L) for (j in seq_len(i)) {
      slice <- u[j] * (8 * product - 2 * identity)
      slice[j, ] <- slice[j, ] - 2 * u
      slice[, j] <- slice[, j] - 2 * u
      out$third[ix, ix, ix[j]] <- (weighted_inv2 * inv) * (-2 * slice)
    }
    at <- at + i
  }
  out
}

.logistic_normal_derivatives <- function(z, jacobian, eta, wt, tri,
                                         order = 0L) {
  n <- nrow(z)
  m <- ncol(z)
  q <- ncol(eta)
  if (is.null(wt)) wt <- rep.int(1, n)
  active <- which(wt > 0)
  n2 <- q * (q + 1) / 2
  n3 <- n2 * (q + 2) / 3
  out <- list(l = 0, l0 = numeric(n))
  if (order >= 1L) out$l1 <- matrix(0, n, q)
  if (order >= 2L) out$l2 <- matrix(0, n, n2)
  if (order >= 3L) out$l3 <- matrix(0, n, n3)
  if (!length(active)) return(out)
  if (any(!is.finite(eta[active, , drop = FALSE]))) {
    out$l0[active] <- -Inf
    return(list(l = -Inf, l0 = out$l0))
  }

  lower_rows <- rep(seq_len(m), seq_len(m) - 1L)
  lower_cols <- sequence(seq_len(m) - 1L)
  predictors <- list(row = c(seq_len(m), seq_len(m), lower_rows),
    col = c(integer(2L * m), lower_cols),
    kind = rep.int(1:3, c(m, m, length(lower_rows))))
  # The result itself has packed derivatives per observation. Batch the
  # triangular solve for second derivatives to bound additional storage.
  batch_size <- if (order >= 3L) max(1L, floor(2^20 / (m * n2))) else
    length(active)
  for (first in seq.int(1L, length(active), by = batch_size)) {
    rows <- active[first:min(length(active), first + batch_size - 1L)]
    batch <- .logistic_normal_batch(z[rows, , drop = FALSE], jacobian[rows],
      eta[rows, , drop = FALSE], wt[rows], predictors, tri, order, n2, n3)
    if (is.null(batch)) {
      out$l0[rows] <- -Inf
      return(list(l = -Inf, l0 = out$l0))
    }
    out$l0[rows] <- batch$l0
    if (order >= 1L) out$l1[rows, ] <- batch$l1
    if (order >= 2L) out$l2[rows, ] <- batch$l2
    if (order >= 3L) out$l3[rows, ] <- batch$l3
  }
  out$l <- sum(out$l0)
  out
}

.logistic_normal_batch <- function(z, jacobian, eta, wt, predictors, tri,
                                   order, n2, n3) {
  n <- nrow(z)
  m <- ncol(z)
  q <- ncol(eta)
  means <- seq_len(m)
  log_sd <- eta[, m + means, drop = FALSE]
  cor_cols <- lapply(means, function(i) 2L * m + (i - 1L) * (i - 2L) / 2 +
    seq_len(i - 1L))
  log_norm <- matrix(0, n, m)
  normalized <- matrix(0, n, q)
  if (m > 1L) for (i in 2L:m) {
    ix <- cor_cols[[i]]
    row <- .logistic_normal_row_norm(eta[, ix, drop = FALSE])
    normalized[, ix] <- row$normalized
    log_norm[, i] <- row$log_norm
  }

  # C = diag(1 / s) T. For a standardized residual r, solve T x = s * r
  # with T unit lower triangular and s the vector of T's row norms.
  forward <- function(rhs) {
    if (m > 1L) for (i in 2L:m) {
      rhs[, i] <- rhs[, i] - rowSums(eta[, cor_cols[[i]], drop = FALSE] *
        rhs[, seq_len(i - 1L), drop = FALSE])
    }
    rhs
  }
  scale_residual <- exp(log_norm - log_sd)
  rhs <- (z - eta[, means, drop = FALSE]) * scale_residual
  white <- forward(rhs)
  out <- list(l0 = wt * (jacobian - m / 2 * log(2 * pi) +
    rowSums(log_norm - log_sd) - rowSums(white^2) / 2))
  if (any(!is.finite(out$l0))) return(NULL)
  if (!order) return(out)

  # The reverse solve gives x' T^-1 applied to any differentiated RHS.
  # This avoids solving for, or storing, third derivatives of x.
  adjoint <- white
  if (m > 1L) for (i in m:2L) {
    ix <- seq_len(i - 1L)
    adjoint[, ix] <- adjoint[, ix, drop = FALSE] -
      eta[, cor_cols[[i]], drop = FALSE] * adjoint[, i]
  }
  inv_norm <- exp(-log_norm)

  # For correlation variables a, b, c in one row, normalized contains
  # t/s. The returned derivatives are either of log(s) or of s divided
  # by s, so large norms need never be constructed on the original scale.
  norm_derivative <- function(ix, i, logarithm) {
    degree <- length(ix)
    a <- ix[1L]
    ha <- normalized[, a]
    inv <- inv_norm[, i]
    if (degree == 1L) return(ha * inv)
    b <- ix[2L]
    hb <- normalized[, b]
    if (degree == 2L) {
      factor <- if (logarithm) 2 else 1
      return(((a == b) - factor * ha * hb) * inv^2)
    }
    c <- ix[3L]
    hc <- normalized[, c]
    mixed <- (a == b) * hc + (a == c) * hb + (b == c) * ha
    if (logarithm) (8 * ha * hb * hc - 2 * mixed) * inv^3 else
      (3 * ha * hb * hc - mixed) * inv^3
  }
  rhs_derivative <- function(indices) {
    i <- predictors$row[indices[1L]]
    if (any(predictors$row[indices] != i)) return(0)
    kind <- predictors$kind[indices]
    nmean <- sum(kind == 1L)
    if (nmean > 1L) return(0)
    ans <- if (nmean == 1L) -scale_residual[, i] else rhs[, i]
    if (sum(kind == 2L) %% 2L) ans <- -ans
    ix <- indices[kind == 3L]
    if (length(ix)) ans <- ans * norm_derivative(ix, i, logarithm = FALSE)
    ans
  }
  logdet_derivative <- function(indices) {
    kind <- predictors$kind[indices]
    if (length(indices) == 1L && kind == 2L) return(-1)
    i <- predictors$row[indices[1L]]
    if (all(kind == 3L) && all(predictors$row[indices] == i))
      return(norm_derivative(indices, i, logarithm = TRUE))
    0
  }
  # T is linear in its free entries. Each T_a x lower-order derivative
  # contributes just one column of x to one row of the RHS.
  change_rhs <- function(indices, lower) {
    b <- matrix(0, n, m)
    i <- predictors$row[indices[1L]]
    b[, i] <- rhs_derivative(indices)
    for (j in seq_along(lower)) {
      a <- indices[j]
      k <- predictors$col[a]
      if (k > 0L) {
        i <- predictors$row[a]
        b[, i] <- b[, i] - lower[[j]][, k]
      }
    }
    b
  }
  dot <- function(a, b) rowSums(a * b)
  first <- vector("list", q)
  out$l1 <- matrix(0, n, q)
  for (a in seq_len(q)) {
    b <- change_rhs(a, list(white))
    first[[a]] <- forward(b)
    out$l1[, a] <- wt * (logdet_derivative(a) - dot(adjoint, b))
  }
  if (order >= 2L) {
    if (order >= 3L) second <- vector("list", n2)
    out$l2 <- matrix(0, n, n2)
    for (a in seq_len(q)) for (b in a:q) {
      ab <- tri$i2[a, b]
      indices <- c(a, b)
      v <- change_rhs(indices, list(first[[b]], first[[a]]))
      out$l2[, ab] <- wt * (logdet_derivative(indices) -
        dot(first[[a]], first[[b]]) - dot(adjoint, v))
      if (order >= 3L) second[[ab]] <- forward(v)
    }
  }
  if (order >= 3L) {
    out$l3 <- matrix(0, n, n3)
    for (a in seq_len(q)) for (b in a:q) for (c in b:q) {
      indices <- c(a, b, c)
      ab <- second[[tri$i2[a, b]]]
      ac <- second[[tri$i2[a, c]]]
      bc <- second[[tri$i2[b, c]]]
      v <- change_rhs(indices, list(bc, ac, ab))
      out$l3[, tri$i3[a, b, c]] <- wt * (logdet_derivative(indices) -
        dot(ab, first[[c]]) - dot(ac, first[[b]]) - dot(bc, first[[a]]) -
        dot(adjoint, v))
    }
  }
  if (!all(vapply(out, function(x) all(is.finite(x)), logical(1L)))) return(NULL)
  out
}
