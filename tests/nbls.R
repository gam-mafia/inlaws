library(inlaws)

near <- function(x, y, tol = 1e-6) {
  stopifnot(length(x) == length(y), all(is.finite(x)), all(is.finite(y)),
            max(abs(x - y)) <= tol * max(1, abs(x), abs(y)))
}
expect_error <- function(expr, pattern) {
  e <- tryCatch(expr, error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
}
set.seed(723)
n <- 700L
d <- data.frame(x = runif(n), z = runif(n), a = runif(n, -.4, .4),
                b = runif(n, -.2, .2))
d$mu <- exp(1 + .7 * sin(2 * pi * d$x) + d$a)
d$theta <- exp(.5 + 1.5 * d$z + d$b)
d$y <- rnbinom(n, mu = d$mu, size = d$theta)
form <- list(y ~ s(x, k = 7) + offset(a), ~ s(z, k = 5) + offset(b))
control <- gam.control(efs.tol = 1e-5)
fits <- lapply(list(c("outer", "newton"), c("outer", "bfgs"), "efs"), function(opt)
  gam(form, data = d, family = nbls(), method = "REML", optimizer = opt, control = control))
f <- fits[[1L]]
for (g in fits) {
  stopifnot(all(is.finite(coef(g))), all(is.finite(g$Vp)),
            all(g$fitted.values > 0), all(is.finite(g$sp)),
            identical(g$outer.info$conv, "full convergence"))
  near(g$gcv.ubre, f$gcv.ubre, tol = 1e-5)
  near(log(g$fitted.values), log(f$fitted.values), tol = .01)
  near(g$l, sum(dnbinom(d$y, mu = g$fitted.values[, 1],
                       size = g$fitted.values[, 2], log = TRUE)))
  near(predict(g, type = "response"), g$fitted.values)
}
stopifnot(cor(log(f$fitted.values[, 1]), log(d$mu)) > .9,
          cor(log(f$fitted.values[, 2]), log(d$theta)) > .8)

# Prediction standard errors agree with model-matrix covariance calculations.
nd <- d[1:11, ]
X <- predict(f, nd, type = "lpmatrix")
lp <- predict(f, nd, type = "link", se.fit = TRUE)
pr <- predict(f, nd, type = "response", se.fit = TRUE)
for (j in 1:2) {
  ind <- attr(X, "lpi")[[j]]
  Xi <- X[, ind, drop = FALSE]
  se <- sqrt(rowSums((Xi %*% f$Vp[ind, ind]) * Xi))
  near(lp$se.fit[, j], se)
  near(pr$fit[, j], exp(lp$fit[, j]))
  near(pr$se.fit[, j], pr$fit[, j] * se)
}
near(f$deviance, sum(residuals(f, type = "deviance")^2))
near(residuals(f, type = "response"), d$y - f$fitted.values[, 1])
near(residuals(f, type = "pearson"), (d$y - f$fitted.values[, 1]) /
       sqrt(f$fitted.values[, 1] + f$fitted.values[, 1]^2 / f$fitted.values[, 2]))
near(as.numeric(logLik(f)), f$l)
near(AIC(f), -2 * as.numeric(logLik(f)) + 2 * attr(logLik(f), "df"))
null <- optimize(function(a) -sum(dnbinom(d$y, mu = exp(a + d$a),
                  size = f$fitted.values[, 2], log = TRUE)), c(-10, 10), tol = 1e-9)
sat <- sum(dnbinom(d$y, mu = d$y, size = f$fitted.values[, 2], log = TRUE))
near(f$null.deviance, 2 * (sat + null$objective))

# Distribution callbacks at fitted and predicted mean/size parameters.
fam <- f$family
near(fam$cdf(nd$y, pr$fit), pnbinom(nd$y, mu = pr$fit[, 1], size = pr$fit[, 2]))
near(fam$cdf(nd$y, pr$fit, logp = TRUE),
     pnbinom(nd$y, mu = pr$fit[, 1], size = pr$fit[, 2], log.p = TRUE))
for (p in c(.2, .8))
  near(fam$qf(p, pr$fit), qnbinom(p, mu = pr$fit[, 1], size = pr$fit[, 2]))
set.seed(1)
a <- fam$rd(fitted(f))
set.seed(1)
near(a, rnbinom(n, mu = fitted(f)[, 1], size = fitted(f)[, 2]))
stopifnot(length(fam$rd(pr$fit[1, , drop = FALSE])) == 1L)

# Constant predictors give the joint MLE, independent of nb()'s REML treatment.
constant <- gam(list(y ~ 1, ~ 1), data = d, family = nbls(), method = "REML")
mle <- optim(coef(constant) + .1, function(b)
  -sum(dnbinom(d$y, mu = exp(b[1]), size = exp(b[2]), log = TRUE)),
  method = "BFGS", control = list(reltol = 1e-12))
near(coef(constant), mle$par, tol = 1e-5)

# Integer likelihood weights equal replication; zero-weight observations are
# omitted. These parametric models avoid changing data-dependent spline bases.
d$w <- sample(0:3, n, replace = TRUE)
weighted <- gam(list(y ~ x + offset(a), ~ z + offset(b)), data = d, weights = w,
                family = nbls(), method = "REML")
replicated <- gam(list(y ~ x + offset(a), ~ z + offset(b)),
                  data = d[rep(seq_len(n), d$w), ], family = nbls(), method = "REML")
near(coef(weighted), coef(replicated), tol = 1e-5)
omitted <- gam(list(y ~ x + offset(a), ~ z + offset(b)), data = d[d$w > 0, ],
               weights = w, family = nbls(), method = "REML")
near(coef(weighted), coef(omitted), tol = 1e-5)
near(residuals(weighted, "pearson"), sqrt(d$w) *
       (d$y - weighted$fitted.values[, 1]) / sqrt(weighted$fitted.values[, 1] +
       weighted$fitted.values[, 1]^2 / weighted$fitted.values[, 2]))

# Alternative links, user starts, and rank deficiency.
for (link in c("identity", "sqrt")) {
  alt <- gam(form, data = d, family = nbls(list(link, "log")), method = "REML")
  near(predict(alt, type = "response"), alt$fitted.values)
  near(alt$l, sum(dnbinom(d$y, mu = alt$fitted.values[, 1],
                        size = alt$fitted.values[, 2], log = TRUE)))
  alt2 <- gam(form, data = d, family = nbls(list(link, "log")), method = "REML",
               start = coef(alt), sp = alt$sp)
  near(coef(alt2), coef(alt), tol = 1e-4)
}
provided <- c(1, .2, .5, .3)
init <- getFromNamespace(".nbls_initialize", "inlaws")
X <- cbind(1, d$x, 1, d$z); attr(X, "lpi") <- list(1:2, 3:4)
stopifnot(identical(init(d$y, X, matrix(0, 0, 4), rep(1, n), NULL,
                         provided, nbls()), provided))
d$xcopy <- d$x
rankdef <- gam(list(y ~ x + xcopy, ~ z), data = d, family = nbls(), method = "REML")
fullrank <- gam(list(y ~ x, ~ z), data = d, family = nbls(), method = "REML")
near(rankdef$fitted.values, fullrank$fitted.values, tol = 1e-5)

# A real shared smooth exercises mgcv's setup, reparameterization, and the
# shared fourth-order trace path, rather than only duplicated design columns.
sharedform <- list(y ~ 1 + offset(a), ~ 1 + offset(b), 1 + 2 ~ s(x, k = 5) - 1)
shared <- gam(sharedform, data = d, family = nbls(), method = "REML")
sharedb <- gam(sharedform, data = d, family = nbls(), method = "REML",
               optimizer = c("outer", "bfgs"))
near(shared$gcv.ubre, sharedb$gcv.ubre, tol = 1e-5)
near(predict(shared, type = "response"), shared$fitted.values)
stopifnot(length(shared$family$rd(predict(shared, nd, type = "response"))) == 11L)

# Input validation and invalid line-search trials.
for (links in list("log", list("log", "identity"), list("inverse", "log"),
                  list(NA_character_, "log"))) expect_error(nbls(links), "two links")
expect_error(init(c(-1, 1), X[1:2, ], matrix(0, 0, 4), c(1, 1), NULL, NULL, nbls()), "counts")
expect_error(init(c(.5, 1), X[1:2, ], matrix(0, 0, 4), c(1, 1), NULL, NULL, nbls()), "counts")
expect_error(init(d$y, X, matrix(0, 0, 4), rep(-1, n), NULL, NULL, nbls()), "weights")
expect_error(init(d$y, X, matrix(0, 0, 4), rep(Inf, n), NULL, NULL, nbls()), "weights")
expect_error(init(d$y, X, matrix(0, 0, 4), rep(0, n), NULL, NULL, nbls()), "positive-weight")
expect_error(gam(list(y ~ 1, ~ 1), data = data.frame(y = rep(0, 20)), family = nbls()), "all-zero")
expect_error(gam(y ~ x, data = d, family = nbls()), "linear predictors")
for (link in c("log", "identity", "sqrt")) {
  fam <- nbls(list(link, "log"))
  eta <- cbind(rep(if (link == "log") 1000 else -1, n), rep(0, n))
  bad <- fam$ll(d$y, X, provided, rep(1, n), fam, eta = eta, deriv = 1)
  stopifnot(identical(bad$l, -Inf))
}
