library("inlaws")

near <- function(actual, expected, tol = 1e-5) {
  stopifnot(isTRUE(all.equal(as.numeric(actual), as.numeric(expected), tolerance = tol)))
}

set.seed(29)
n <- 140L
d <- data.frame(x = runif(n), a = runif(n, -.1, .1), b = runif(n, -.1, .1))
mu <- exp(.6 + sin(6 * d$x) + d$a)
censor <- function(y) cbind(pmin(y, 3), ifelse(y > 3, Inf, y))
p <- plogis(-.3 + sin(6 * d$x) + d$a)
y <- rbeta(n, shape1 = 8 * p, shape2 = 8 * (1 - p))
total <- sample(8:20, n, replace = TRUE)
success <- rbinom(n, size = total, prob = y)
cases <- list(
  list(family = betabinomial, y = cbind(success, total - success)),
  list(family = cgammals, y = censor(rgamma(n, shape = 3, scale = mu / 3))),
  list(family = clognormalls, y = censor(rlnorm(n, log(mu) - .125, .5))),
  list(family = log_logisticls, y = log_logistic(.25)$rd(mu)),
  list(family = nbls, y = rnbinom(n, mu = mu, size = 2)),
  list(family = dirichlet, y = cbind(1 - y, y)),
  list(family = pigls, y = pig(.4)$rd(mu)),
  list(family = cmpls, y = cmp(.8)$rd(mu))
)
form <- list(y ~ s(x, k = 6) + offset(a), ~ offset(b))
nei <- list(
  a = c(4, 5, 6, 20, 21, 22, 4, 5, 6, 39, 40, 41), ma = c(3, 6, 9, 12),
  d = c(5, 21, 22, 5, 40), md = c(1, 3, 4, 5), jackknife = FALSE
)

for (case in cases) {
  fam <- case$family()
  d$y <- case$y
  fit <- gam(form, family = fam, data = d, method = "REML")
  if (fam$family == "log_logisticls") {
    bfgs <- gam(form, family = case$family(), data = d, method = "REML",
      optimizer = c("outer", "bfgs"))
    near(bfgs$gcv.ubre, fit$gcv.ubre, 1e-5)
  }
  for (method in c("NCV", "QNCV")) {
    cv_fit <- gam(form, family = case$family(), data = d, method = method, gamma = 1.3)
    stopifnot(cv_fit$method == method, is.finite(cv_fit$gcv.ubre),
      length(cv_fit$sp) == 1L, all(is.finite(coef(cv_fit))))
  }

  # Reordered, repeated predictions from a subset of the training rows.
  # Check the gradient with one smoothing direction and both predictor offsets.
  X <- predict(fit, type = "lpmatrix")
  beta <- coef(fit)
  offset <- list(d$a, d$b)
  wt <- seq(.7, 1.3, length.out = n)
  db <- matrix(seq(-.2, .2, length.out = length(beta)), ncol = 1L)
  fam <- fit$family
  score <- function(step, want_gradient = FALSE) {
    b <- beta + step * db[, 1L]
    ll <- fam$ll(fit$y, X, b, wt, fam, offset = offset,
      deriv = if (want_gradient) 3 else 1, d1b = db, ncv = TRUE)
    ll$gamma <- 1.3
    H <- -ll$lbb + diag(1, length(b))
    fam$ncv(X, fit$y, wt, nei, b, fam, ll, R = chol(H), offset = offset,
      dH = ll$d1H, db = db, deriv = want_gradient)
  }
  for (quadratic in c(FALSE, TRUE)) {
    fam$qapprox <- quadratic
    center <- score(0, TRUE)
    near(center$NCV1, (score(1e-4)$NCV - score(-1e-4)$NCV) / 2e-4)
    stopifnot(identical(dim(center$Vg), c(1L, 1L)))
  }
}
