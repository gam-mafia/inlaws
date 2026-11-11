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
