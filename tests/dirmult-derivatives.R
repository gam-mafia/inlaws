library("inlaws")

dderiv <- getFromNamespace(".dirmult_deriv", "inlaws")
rising <- getFromNamespace(".dirmult_rising", "inlaws")
near <- function(x, y, tol = 1e-7, relative = FALSE) {
  error <- abs(x - y) / if (relative) abs(y) else (1 + abs(y))
  stopifnot(length(x) == length(y), all(is.finite(error)), all(error <= tol))
}

# The logarithm of the rising factorial is a finite sum. Each term's
# derivatives with respect to log concentration are Bernoulli cumulants.
rising_reference <- function(la, m) {
  if (m == 0) return(c(r = 0, A = 0, dA = 0, B = 0, D2 = 0, D3 = 0, D4 = 0))
  z <- log(seq_len(m - 1L)) - la
  q <- stats::plogis(-z)
  u <- stats::plogis(z)
  v <- q * u
  c(r = sum(pmax(z, 0) + log1p(exp(-abs(z)))),
    A = 1 + sum(q), dA = -sum(u), B = -1 - sum(q^2),
    D2 = sum(v), D3 = sum(v * (u - q)), D4 = sum(v * (1 - 6 * v)))
}
m <- c(0, 1, 64, 65, 65, 80, 80, 1000, 10000, 70, 70, 70)
la <- c(1, -900, -.4, -.4, 0, log(80 * 999), log(80 * 1001),
        log(.001), log(600), -900, -700, 700)
m <- c(m, 65, 65, 80, 80)
la <- c(la, log(c(31.99, 32.01, 80 / .0501, 80 / .0499)))
expected <- vapply(seq_along(m), function(i) rising_reference(la[i], m[i]), numeric(7))
actual <- rising(la, m, order = 4L)
for (name in rownames(expected)) near(actual[[name]], expected[name, ])
# Relative comparisons detect cancellation to zero in a representable tail.
for (name in c("D2", "D3", "D4"))
  near(actual[[name]][11:12], expected[name, 11:12], relative = TRUE)
near(actual$dA[12], expected["dA", 12], relative = TRUE)
large <- rising(c(-900, log(1e8), 700), c(1e8, 1e8, 1e8), order = 4L)
stopifnot(all(is.finite(unlist(large))))

# Independently enumerate nondecreasing tensor indices and use mgcv's index
# arrays to locate each mixed partial. Differentiate each preceding order.
tuples <- function(q, order) {
  x <- as.matrix(expand.grid(rep(list(seq_len(q)), order)))
  x[apply(x, 1L, function(i) all(diff(i) >= 0)), , drop = FALSE]
}
column <- function(tri, i) {
  if (length(i) == 1L) return(i)
  do.call(`[`, c(list(tri[[paste0("i", length(i))]]), as.list(i)))
}
mass <- function(y, eta) {
  logits <- cbind(0, eta[, -ncol(eta), drop = FALSE])
  alpha <- exp(eta[, ncol(eta)]) * exp(logits) / rowSums(exp(logits))
  n <- rowSums(y)
  phi <- rowSums(alpha)
  lgamma(n + 1) - rowSums(lgamma(y + 1)) +
    lgamma(phi) - lgamma(phi + n) + rowSums(lgamma(alpha + y) - lgamma(alpha))
}
set.seed(1806)
for (q in c(2L, 4L, 6L)) {
  y <- matrix(sample(0:35, 4L * q, replace = TRUE), 4L)
  y[1L, 1L] <- 0
  y[4L, ] <- 0
  eta <- matrix(stats::rnorm(4L * q, sd = .7), 4L)
  w <- c(.4, 2, 0, 1)
  tri <- mgcv::trind.generator(q)
  z <- dderiv(y, eta, w, order = 4L)
  target <- mass(y, eta)
  near(z$l0, target)
  near(z$l, sum(w * target))
  h <- 1e-4
  for (j in seq_len(q)) {
    ep <- em <- eta
    ep[, j] <- ep[, j] + h
    em[, j] <- em[, j] - h
    near(z$l1[, j], w * (mass(y, ep) - mass(y, em)) / (2 * h), 2e-6)
    pp <- dderiv(y, ep, w, order = 3L)
    mm <- dderiv(y, em, w, order = 3L)
    for (d in 1:3) {
      ids <- tuples(q, d)
      low <- paste0("l", d)
      high <- paste0("l", d + 1L)
      stopifnot(ncol(z[[high]]) == choose(q + d, d + 1L))
      for (k in seq_len(nrow(ids))) {
        i <- column(tri, ids[k, ])
        index <- column(tri, sort(c(ids[k, ], j)))
        near(z[[high]][, index], (pp[[low]][, i] - mm[[low]][, i]) / (2 * h), 2e-6)
      }
    }
  }
}

# A total of one is an ordinary categorical trial independently of dispersion.
# The exact logit derivatives are the signed Bernoulli cumulants.
y <- rbind(c(0, 1), c(1, 0), c(0, 0))
eta <- cbind(c(-.7, .4, 1000), c(-1000, 1000, -1000))
w <- c(.6, 2, 1)
p <- stats::plogis(eta[1:2, 1])
q <- stats::plogis(-eta[1:2, 1])
v <- p * q
z <- dderiv(y, eta, w, order = 4L)
near(z$l0, c(log(p[1]), log(q[2]), 0), 1e-12)
exact <- list(c(1, 0) - p, -v, -v * (q - p), -v * (1 - 6 * v))
for (d in 1:4) {
  ld <- z[[paste0("l", d)]]
  near(ld[, 1], c(w[1:2] * exact[[d]], 0), 1e-12)
  near(ld[, -1, drop = FALSE], rep(0, length(ld[, -1, drop = FALSE])), 1e-12)
}

# When a category shape vanishes, a positive count m contributes
# log(alpha) + lgamma(m). Subtract the exact total-concentration derivatives.
m <- 79L
log_phi <- .4
eta <- matrix(c(-900, log_phi), 1)
z <- dderiv(matrix(c(0, m), 1), eta, 1.3, order = 4L)
phi <- exp(log_phi)
near(z$l0, -900 + log_phi + lgamma(m) - sum(log(phi + 0:(m - 1L))), 1e-12)
tri <- mgcv::trind.generator(2)
a <- rising_reference(log_phi, m)
near(z$l1[1, ], 1.3 * c(1, 1 - a["A"]), 1e-12)
for (d in 2:4) {
  ld <- z[[paste0("l", d)]]
  pos <- column(tri, rep(2L, d))
  near(ld[1, pos], -1.3 * a[paste0("D", d)], 1e-11)
  near(ld[1, -pos], rep(0, ncol(ld) - 1L), 1e-12)
}

# Finite predictors may imply concentration or category shapes beyond the
# double-precision exponent range. Entirely unused rows may be non-finite.
y <- rbind(c(2, 0, 10000), c(10000, 4, 20), c(70, 5, 0), c(2, 30, 1))
eta <- rbind(c(-1000, 1000, -1000), c(1000, -1000, 1000),
             c(-.1, .5, 1000), c(Inf, -Inf, NaN))
z <- dderiv(y, eta, c(1, .3, 2, 0), order = 4L)
stopifnot(is.finite(z$l), all(is.finite(z$l0[1:3])),
          all(is.finite(unlist(z[paste0("l", 1:4)]))))
for (d in 1:4) stopifnot(all(z[[paste0("l", d)]][4, ] == 0))
p <- exp(c(0, eta[3, 1:2]))
near(z$l0[3], stats::dmultinom(y[3, ], prob = p, log = TRUE), 1e-12)

# Extreme dispersion leaves tiny but representable derivatives. In the
# one-category limit only the log probability of the observed type remains.
z <- dderiv(matrix(c(0, 20), 1), matrix(c(40, -900), 1), order = 4L)
q <- stats::plogis(-40)
near(z$l0, -log1p(exp(-40)), 1e-12, relative = TRUE)
for (d in 1:4)
  near(z[[paste0("l", d)]][1, 1], q * (-1)^(d - 1L), 1e-12, relative = TRUE)

# The leading multinomial correction is exp(-logphi) times the first
# normalized power sum, here 14*13 + 6*5 - 20*19/2 = 22.
z <- dderiv(matrix(c(14, 6), 1), matrix(c(0, 700), 1), order = 4L)
for (d in 1:4) {
  pos <- column(tri, rep(2L, d))
  near(z[[paste0("l", d)]][1, pos], 22 * exp(-700) * (-1)^d,
       1e-12, relative = TRUE)
}
