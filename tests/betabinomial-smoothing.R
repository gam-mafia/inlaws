library(inlaws)

near <- function(x, y, tol = 5e-5) {
  stopifnot(length(x) == length(y), all(is.finite(x)), all(is.finite(y)),
            max(abs(x - y)) <= tol * max(1, abs(x), abs(y)))
}
set.seed(811)
d <- data.frame(x = runif(180), z = runif(180), m = sample(10:25, 180, TRUE),
                a = runif(180, -.3, .3), b = runif(180, -.2, .2))
mu <- plogis(sin(6 * d$x) + d$a); phi <- exp(1 + d$z + d$b)
d$y <- rbinom(180, d$m, rbeta(180, mu * phi, (1 - mu) * phi))
d$w <- seq(.5, 1.5, length.out = nrow(d)); d$w[7] <- 0
G <- gam(list(cbind(y, m - y) ~ s(x, k = 5) + offset(a),
              ~ s(z, k = 4) + offset(b)), data = d, weights = w, family = betabinomial(), fit = FALSE)
internal <- function(s) getFromNamespace(s, "mgcv")
G$Sl <- internal("Sl.setup")(G)
G$X <- internal("Sl.initial.repara")(G$Sl, G$X, both.sides = FALSE)
pini <- G$family$preinitialize(G); G[names(pini)] <- pini
G$family$qapprox <- FALSE
nei <- list(a = as.integer(c(1, 2, 3, 4, 5, 6)), ma = c(2L, 4L, 6L),
            d = c(1L, 2L, 3L, 3L, 6L), md = c(2L, 4L, 5L), jackknife = FALSE)
fit <- function(sp, method = "REML", deriv = 2, gamma = 1) internal("gam.fit5")(
  G$X, G$y, sp, G$Sl, weights = G$w, offset = G$offset, deriv = deriv,
  family = G$family, scoreType = method, gamma = gamma, control = gam.control(epsilon = 1e-10),
  nei = if (method == "NCV") nei else NULL)
# Test away from the optimum, where an incorrectly zero gradient cannot pass.
s <- c(-.3, .5); h <- 1e-3
for (method in c("REML", "NCV")) {
  f <- fit(s, method, if (method == "NCV") 1 else 2)
  for (j in 1:2) {
    ds <- c(0, 0); ds[j] <- h
    lo <- fit(s - ds, method, 1); hi <- fit(s + ds, method, 1)
    score <- if (method == "NCV") "NCV" else "REML"
    near((hi[[score]] - lo[[score]]) / (2 * h), f[[paste0(score, "1")]][j], 1e-5)
    if (method == "REML") near((hi$REML1 - lo$REML1) / (2 * h), f$REML2[, j])
  }
}
# NCV's gamma adjustment must differentiate the full-data term as well.
for (j in 1:2) {
  ds <- c(0, 0); ds[j] <- h
  f <- fit(s, "NCV", 1, gamma = 1.2)
  lo <- fit(s - ds, "NCV", 1, gamma = 1.2)
  hi <- fit(s + ds, "NCV", 1, gamma = 1.2)
  near((hi$NCV - lo$NCV) / (2 * h), f$NCV1[j], 1e-5)
}
