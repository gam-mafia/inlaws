library(inlaws)
quadrature <- getFromNamespace(".gauss_legendre", "inlaws")

# The n-point rule integrates polynomials through degree 2*n - 1 exactly.
# These checks provide independent references without a statmod dependency.
for (n in c(1L, 2L, 3L, 12L, 16L)) {
  rule <- quadrature(n)
  x <- rule$nodes
  w <- rule$weights
  stopifnot(length(x) == n, length(w) == n,
            all(is.finite(x)), all(is.finite(w)),
            all(abs(x) < 1), all(diff(x) > 0), all(w > 0),
            max(abs(x + rev(x))) < 1e-13,
            max(abs(w - rev(w))) < 1e-13)
  for (degree in 0:(2L * n - 1L)) {
    expected <- if (degree %% 2L) 0 else 2 / (degree + 1)
    stopifnot(abs(sum(w * x^degree) - expected) < 1e-13)
  }
  # Nodes are roots of the degree-n Legendre polynomial (three-term recurrence).
  previous <- rep(1, n)
  current <- x
  if (n > 1L) for (k in 2:n) {
    next_value <- ((2 * k - 1) * x * current - (k - 1) * previous) / k
    previous <- current
    current <- next_value
  }
  stopifnot(max(abs(current)) < 1e-12)

  # Cached rules are stable and cannot be corrupted by editing a returned list.
  stopifnot(identical(rule, quadrature(as.double(n))))
  rule$nodes[1] <- 100
  rule$weights[1] <- 100
  stopifnot(identical(quadrature(n), list(nodes = x, weights = w)))
}

# A known two-node rule and a smooth non-polynomial integral.
rule <- quadrature(2L)
stopifnot(max(abs(rule$nodes - c(-1, 1) / sqrt(3))) < 1e-14,
          max(abs(rule$weights - 1)) < 1e-14)
for (n in c(12L, 16L)) {
  rule <- quadrature(n)
  stopifnot(abs(sum(rule$weights * exp(rule$nodes)) - (exp(1) - exp(-1))) < 1e-13)
}
for (bad in list(0, -1, 1.5, NA_real_, Inf, numeric(), c(12, 16), "12", TRUE)) {
  stopifnot(inherits(try(quadrature(bad), silent = TRUE), "try-error"))
}
