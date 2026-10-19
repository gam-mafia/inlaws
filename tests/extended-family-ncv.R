library("inlaws")

# Reject the affected upstream GAM route instead of maintaining a local solver.
if (packageVersion("mgcv") == "1.9.4") {
  set.seed(83)
  d <- data.frame(x = runif(24), y = rep(1:4, 6))
  grouped <- list(a = 1:24, ma = seq.int(3L, 24L, 3L),
    d = 1:24, md = seq.int(3L, 24L, 3L))
  theta <- c(cgamma = .4, clognormal = .5, cmp = 1.5, gpoisson = 2)
  for (name in names(theta)) for (fixed in c(FALSE, TRUE)) {
    for (method in c("NCV", "QNCV")) for (nei in list(NULL, grouped)) {
      family <- get(name)(if (fixed) theta[[name]] else NULL)
      initial <- family$getTheta(TRUE)
      err <- tryCatch(gam(y ~ s(x, k = 4), data = d, family = family,
        method = method, nei = nei), error = identity)
      stopifnot(inherits(err, "error"),
        grepl("mgcv 1.9-4", conditionMessage(err), fixed = TRUE),
        grepl("use REML or ML", conditionMessage(err), fixed = TRUE),
        identical(family$getTheta(TRUE), initial))
    }
    # A request for NCV diagnostics via nei must fail just as explicitly.
    err <- tryCatch(gam(y ~ s(x, k = 4), data = d,
      family = get(name)(if (fixed) theta[[name]] else NULL),
      method = "REML", nei = grouped), error = identity)
    stopifnot(inherits(err, "error"), grepl("mgcv 1.9-4", conditionMessage(err), fixed = TRUE))
  }
}
