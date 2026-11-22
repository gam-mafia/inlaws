.dirmult_prob <- function(eta) {
  z <- cbind(0, eta[, -ncol(eta), drop = FALSE])
  maximum <- cbind(seq_len(nrow(z)), max.col(z, ties.method = "first"))
  z <- z - z[maximum]
  mass <- exp(z)
  mass[maximum] <- 0
  lp <- z - log1p(rowSums(mass))
  list(p = exp(lp), lp = lp, logphi = eta[, ncol(eta)])
}

# EFS only needs i2. The fourth-order array has C^4 entries, so cache each
# higher-order index only when a fitting route requests it.
.dirmult_tri <- function(C) {
  tri <- new.env(parent = emptyenv())
  i2 <- matrix(0L, C, C)
  col <- 0L
  for (i in seq_len(C)) for (j in i:C) {
    col <- col + 1L
    i2[i, j] <- i2[j, i] <- col
  }
  tri$i2 <- i2
  delayedAssign("i3", {
    i3 <- array(0L, rep.int(C, 3L))
    col <- 0L
    for (i in seq_len(C)) for (j in i:C) for (k in j:C) {
      col <- col + 1L
      i3[i, j, k] <- i3[i, k, j] <- i3[j, i, k] <-
        i3[j, k, i] <- i3[k, i, j] <- i3[k, j, i] <- col
    }
    i3
  }, assign.env = tri)
  delayedAssign("i4", mgcv::trind.generator(C, reverse = FALSE)$i4, assign.env = tri)
  tri
}

# Derivatives of a common log shape of order 3 or 4 are the negatives of
# joint categorical cumulants. Grouping by index multiplicity preserves
# factors of 1 - p when probabilities approach one.
.dirmult_logshape <- function(index, p, complement) {
  if (any(index > ncol(p))) return(numeric(nrow(p)))
  distinct <- unique(index)
  multiplicity <- tabulate(match(index, distinct))
  distinct <- distinct[order(multiplicity, decreasing = TRUE)]
  a <- p[, distinct[1L]]
  q <- complement[, distinct[1L]]
  if (length(index) == 3L)
    return(switch(length(distinct),
      -a * q * (q - a),
      a * p[, distinct[2L]] * (q - a),
      -2 * a * p[, distinct[2L]] * p[, distinct[3L]]))
  switch(length(distinct),
    -a * q * (1 - 6 * a * q),
    {
      b <- p[, distinct[2L]]
      if (max(multiplicity) == 3L) a * b * (1 - 6 * a * q)
      else a * b * ((q - a) * (1 - 2 * b) + 2 * a * b)
    },
    -2 * a * p[, distinct[2L]] * p[, distinct[3L]] * (1 - 3 * a),
    6 * a * p[, distinct[2L]] * p[, distinct[3L]] * p[, distinct[4L]])
}

.dirmult_higher <- function(a, a0, sa, sd, d2j, logshape2, p, complement, order) {
  n <- nrow(p); C <- ncol(p) + 1L
  # J holds the first derivatives of all log shapes for each predictor.
  J <- lapply(seq_len(C - 1L), function(j) {
    z <- matrix(-p[, j], n, C)
    z[, j + 1L] <- complement[, j]
    z
  })
  J[[C]] <- matrix(1, n, C)
  d3 <- matrix(a$D3, n, C)
  l3 <- matrix(0, n, choose(C + 2L, 3L))
  col <- 0L
  for (i in seq_len(C)) for (j in i:C) for (k in j:C) {
    col <- col + 1L
    l3[, col] <- rowSums(d3 * J[[i]] * J[[j]] * J[[k]]) +
      logshape2[, i, j] * d2j[, k] + logshape2[, i, k] * d2j[, j] +
      logshape2[, j, k] * d2j[, i] + sa * .dirmult_logshape(c(i, j, k), p, complement)
  }
  l3[, col] <- l3[, col] - a0$D3
  if (order == 3L) return(list(l3 = l3))

  # The fourth-order chain rule partitions the four indices into two
  # pairs, a triple and a singleton, a pair and two singletons, or all
  # singletons. The last term differentiates a log shape four times.
  d3jj <- array(0, c(n, C, C))
  for (i in seq_len(C)) for (j in i:C)
    d3jj[, i, j] <- rowSums(d3 * J[[i]] * J[[j]])
  d4 <- matrix(a$D4, n, C)
  l4 <- matrix(0, n, choose(C + 3L, 4L))
  col <- 0L
  for (i in seq_len(C)) for (j in i:C) for (k in j:C) for (l in k:C) {
    col <- col + 1L
    hij <- logshape2[, i, j]; hik <- logshape2[, i, k]; hil <- logshape2[, i, l]
    hjk <- logshape2[, j, k]; hjl <- logshape2[, j, l]; hkl <- logshape2[, k, l]
    l4[, col] <- rowSums(d4 * J[[i]] * J[[j]] * J[[k]] * J[[l]]) +
      hij * d3jj[, k, l] + hik * d3jj[, j, l] + hil * d3jj[, j, k] +
      hjk * d3jj[, i, l] + hjl * d3jj[, i, k] + hkl * d3jj[, i, j] +
      (hij * hkl + hik * hjl + hil * hjk) * sd +
      .dirmult_logshape(c(i, j, k), p, complement) * d2j[, l] +
      .dirmult_logshape(c(i, j, l), p, complement) * d2j[, k] +
      .dirmult_logshape(c(i, k, l), p, complement) * d2j[, j] +
      .dirmult_logshape(c(j, k, l), p, complement) * d2j[, i] +
      sa * .dirmult_logshape(c(i, j, k, l), p, complement)
  }
  l4[, col] <- l4[, col] - a0$D4
  list(l3 = l3, l4 = l4)
}

.dirmult_log_mass <- function(y, pr, total, logshape, a, a0) {
  # Shift the first factor out at small shapes. This avoids cancelling
  # m * log(p) against a rising-factor correction of nearly the same size.
  correction <- matrix(a$r, nrow(y), ncol(y))
  denominator <- a0$r
  log_prob <- pr$lp
  log_prob[y == 0] <- 0
  power <- y
  small_total <- pr$logphi < 0 & total > 0
  small <- which(y > 0 & logshape < 0, arr.ind = TRUE)
  if (nrow(small)) {
    shape <- exp(logshape[small])
    m <- y[small]
    logphi <- pr$logphi[small[, 1L]]
    logphi[small_total[small[, 1L]]] <- 0
    correction[small] <- lgamma(shape + m) - lgamma(shape + 1) - (m - 1) * logphi
    power[small] <- 1
  }
  leading <- numeric(length(total))
  if (any(small_total)) {
    phi <- exp(pr$logphi[small_total])
    denominator[small_total] <- lgamma(phi + total[small_total]) - lgamma(phi + 1)
    leading[small_total] <- (rowSums(y[small_total, , drop = FALSE] > 0) - 1) *
      pr$logphi[small_total]
  }
  lgamma(total + 1) - rowSums(lgamma(y + 1)) +
    (rowSums(correction) - denominator) + rowSums(power * log_prob) + leading
}

# Observation derivatives are packed in the lexicographic order used by
# mgcv::trind.generator(). Only the likelihood contributions l0 are unweighted.
.dirmult_deriv <- function(y, eta, wt = NULL, deriv = TRUE,
                           order = if (deriv) 2L else 0L) {
  n <- nrow(y); C <- ncol(y)
  wt <- if (is.null(wt)) rep.int(1, n) else rep_len(wt, n)
  defined <- rowSums(!is.finite(eta)) == 0L
  eta[!defined, ] <- 0
  pr <- .dirmult_prob(eta)
  total <- rowSums(y)
  logshape <- pr$lp + pr$logphi
  a0 <- .dirmult_rising(pr$logphi, total, order)
  a <- .dirmult_rising(as.vector(logshape), as.vector(y), order)
  l0 <- .dirmult_log_mass(y, pr, total, logshape, a, a0)
  l0[!defined | is.nan(l0)] <- -Inf
  active <- wt > 0
  ret <- list(l = sum(wt[active] * l0[active]), l0 = l0)
  if (!order) return(ret)

  p <- pr$p[, -1L, drop = FALSE]
  complement <- -expm1(pr$lp[, -1L, drop = FALSE])
  A <- matrix(a$A, n, C)
  sa <- rowSums(A)
  l1 <- matrix(0, n, C)
  for (j in seq_len(C - 1L))
    l1[, j] <- A[, j + 1L] * complement[, j] - p[, j] * rowSums(A[, -(j + 1L), drop = FALSE])
  l1[, C] <- sa - a0$A
  # Near the multinomial limit use the small corrections to the total
  # count, since subtracting two almost equal scores loses precision.
  large <- a0$A > total / 2
  l1[large, C] <- (rowSums(matrix(a$dA, n, C)) - a0$dA)[large]
  derivatives <- list(l1 = l1)
  if (order > 1L) {
    d2 <- matrix(a$D2, n, C)
    sd <- rowSums(d2)
    d2j <- matrix(0, n, C)
    other_d2 <- matrix(0, n, C - 1L)
    for (j in seq_len(C - 1L)) {
      other_d2[, j] <- rowSums(d2[, -(j + 1L), drop = FALSE])
      d2j[, j] <- d2[, j + 1L] * complement[, j] - p[, j] * other_d2[, j]
    }
    d2j[, C] <- sd
    l2 <- matrix(0, n, choose(C + 1L, 2L))
    if (order > 2L) logshape2 <- array(0, c(n, C, C))
    col <- 0L
    for (j in seq_len(C)) for (k in j:C) {
      col <- col + 1L
      if (k == C) {
        l2[, col] <- if (j == C) sd - a0$D2 else d2j[, j]
      } else {
        h <- p[, j] * if (j == k) -complement[, j] else p[, k]
        direct <- if (j == k) d2[, j + 1L] * complement[, j]^2 + other_d2[, j] * p[, j]^2
          else -d2[, k + 1L] * p[, j] - p[, k] * d2j[, j]
        l2[, col] <- direct + h * sa
        if (order > 2L) logshape2[, j, k] <- h
      }
    }
    derivatives$l2 <- l2
    if (order > 2L) derivatives <- c(derivatives,
      .dirmult_higher(a, a0, sa, sd, d2j, logshape2, p, complement, order))
  }
  for (key in names(derivatives)) {
    value <- derivatives[[key]]
    value[!active | !is.finite(l0), ] <- 0
    ret[[key]] <- value * wt
  }
  ret
}


## Intercept-only reference with the same offsets and row weights.
.dirmult_null <- function(y,wt,offset,eta) {
  C <- ncol(y); n <- nrow(y)
  off <- matrix(0,n,C)
  for (j in seq_len(C)) if (length(offset)>=j && !is.null(offset[[j]]))
    off[,j] <- offset[[j]]
  p <- colSums(y*wt)+0.5; p <- p/sum(p)
  b <- c(log(p[-1]/p[1]),stats::weighted.mean(eta[,C],wt))
  fn <- function(b) -.dirmult_deriv(y,sweep(off,2,b,"+"),wt,FALSE)$l
  gr <- function(b) -colSums(.dirmult_deriv(y,sweep(off,2,b,"+"),wt)$l1)
  fit <- stats::optim(b,fn,gr,method="BFGS",control=list(maxit=200,reltol=1e-10))
  if (fit$convergence!=0) {
    warning("dirmult intercept-only reference did not converge")
    return(NA_real_)
  }
  2*fit$value
}
