library("inlaws")
near <- function(a, b, tol = 1e-5) stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = tol)))
fails <- function(expr) stopifnot(inherits(tryCatch(expr, error = identity), "error"))
evaluate <- getFromNamespace(".cgen_gamma_eval", "inlaws")
set.seed(31); n <- 200L
D <- data.frame(x = runif(n), a = runif(n, -.05, .05), b = runif(n, -.05, .05), c = runif(n, -.02, .02))
y <- gen_gammalss()$rd(cbind(exp(.5 + .5*sin(6*D$x) + D$a), exp(-.7 + D$b), exp(.15 + D$c)))
D$y <- cbind(y, y)
for (i in seq_len(n)) {
  if (i %% 4 == 0 && y[i] < 1) D$y[i, ] <- c(1, -Inf)
  if (i %% 4 == 1 && y[i] > 2) D$y[i, ] <- c(2, Inf)
  if (i %% 4 == 2) D$y[i, ] <- c(y[i]*.9, y[i]*1.1)
}
f <- cgen_gammalss(observe = censoring_process(1, 2))
form <- list(y ~ s(x, k = 5) + offset(a), ~ 1 + offset(b), ~ 1 + offset(c))
check <- function(fit, method) {
  stopifnot(fit$method == method, fit$outer.info$conv == "full convergence",
    all(is.finite(coef(fit))), is.finite(fit$gcv.ubre),
    identical(colnames(fitted(fit)), c("mu", "sigma", "r")), all(fitted(fit) > 0))
  p <- predict(fit, type = "response", se.fit = TRUE)
  l <- predict(fit, type = "link", se.fit = TRUE)
  near(p$fit, fitted(fit)); near(p$fit, exp(l$fit)); near(p$se.fit, p$fit*l$se.fit)
  near(logLik(fit), sum(evaluate(D$y, l$fit)$ell[, 1]), 1e-8)
  near(AIC(fit), -2*as.numeric(logLik(fit)) + 2*attr(logLik(fit), "df"))
  stopifnot(all(is.finite(residuals(fit))), fit$null.deviance > fit$deviance)
  exact <- D$y[, 1] == D$y[, 2]
  for (type in c("response", "pearson")) {
    r <- residuals(fit, type = type)
    stopifnot(all(is.na(r[!exact])), all(is.finite(r[exact])))
  }
}
a <- gam(form, data = D, family = f, method = "REML"); check(a, "REML")
b <- gam(form, data = D, family = f, method = "REML", optimizer = c("outer", "bfgs"))
check(b, "REML"); near(coef(a), coef(b), .001)
for (method in c("NCV", "QNCV")) check(gam(form, data = D, family = f, method = method), method)
nei <- list(a = seq_len(n), ma = seq.int(2L, n, 2L),
  d = as.vector(apply(matrix(seq_len(n), 2L), 2, rev)), md = seq.int(2L, n, 2L), jackknife = FALSE)
for (method in c("NCV", "QNCV")) check(gam(form, data = D, family = f, method = method, nei = nei), method)
# Exact observations reduce to the original family, for matrix or vector input.
exact <- D; exact$y <- y
u <- gam(form, data = exact, family = gen_gammalss(), method = "REML")
v <- gam(form, data = exact, family = cgen_gammalss(), method = "REML")
near(coef(u), coef(v), 1e-9); near(logLik(u), logLik(v), 1e-10)
near(residuals(u), residuals(v)); near(u$null.deviance, v$null.deviance)
# Replication weights agree with physical replication; zero weights are inert.
D$w <- rep(c(0, 1, 2, 1), length.out = n)
forms <- list(y ~ x, ~ 1, ~ 1)
u <- gam(forms, data = D, weights = w, family = f, method = "REML")
v <- gam(forms, data = D[rep(seq_len(n), D$w), ], family = f, method = "REML")
near(coef(u), coef(v)); near(logLik(u), logLik(v))
objective <- function(beta) {
  eta <- cbind(beta[1] + beta[2]*D$x, beta[3], beta[4])
  -sum(D$w*evaluate(D$y, eta)$ell[, 1])
}
opt <- optim(coef(u) + .01, objective, method = "BFGS", control = list(reltol = 1e-12))
near(coef(u), opt$par, 2e-4); near(logLik(u), -opt$value, 1e-7)
# Smooth nuisance predictors, missing rows, subsetting, and serialization.
smooth <- gam(list(y ~ s(x, k = 5), ~ s(x, k = 4), ~ s(x, k = 4)), data = D, family = f, method = "REML")
stopifnot(all(is.finite(coef(smooth))))
D$x[1] <- NA
missing <- gam(form, data = D, family = f, method = "REML", na.action = na.exclude)
stopifnot(nobs(missing) == n-1L, length(residuals(missing)) == n, is.na(residuals(missing)[1]))
sub <- gam(forms, data = D, subset = 2:150, family = f, method = "REML")
stopifnot(nobs(sub) == 149L, length(attr(sub$y, "censor")) == 149L)
copy <- unserialize(serialize(a, NULL))
near(predict(copy, D[2:8, ], type = "response"), predict(a, D[2:8, ], type = "response"))
# Distribution and diagnostics hooks describe the underlying uncensored law.
mu <- cbind(c(1, 2, 3), .7, c(.4, 1, 2))
g <- gen_gammalss()
near(f$cdf(c(.5, 1, 2), mu), g$cdf(c(.5, 1, 2), mu))
near(f$qf(c(.1, .5, .9), mu), g$qf(c(.1, .5, .9), mu))
set.seed(1); draw <- f$rd(mu); set.seed(1); near(draw, g$rd(mu))
dist <- copy$family$gratia_diagnostics(mu, 1, 1, c(1, 2, 3), NULL)
near(dist$mean, mu[, 1]); near(dist$variance, g$gratia_diagnostics(mu, 1, 1, NULL, NULL)$variance)
stopifnot(dist$continuous, !dist$variance_available, dist$reference_available)
observed <- dist$observe(c(.1, 1.5, 4))
near(observed, c(1, 1.5, 2)); near(attr(observed, "censor"), c(-Inf, 1.5, Inf))
near(dist$bounds(observed)$upper, c(1, 1.5, Inf))
stopifnot(!cgen_gammalss()$gratia_diagnostics(mu, 1, 1, NULL, NULL)$reference_available)
# All derivative orders vanish for irrelevant zero-weight and full-support rows.
X <- diag(3)[rep(1, 3), ]; X[,] <- 1; attr(X, "lpi") <- as.list(1:3)
yy <- cbind(c(.Machine$double.xmax, 0, 1), c(Inf, Inf, 1))
et <- matrix(0, 3, 3); et[1, ] <- 1000
ll <- f$ll(yy, X, numeric(3), c(0, 1, 1), f, eta = et, deriv = 4, ncv = TRUE)
stopifnot(all(ll$l0[1:2] == 0), all(ll$l1[1:2, ] == 0), all(ll$l4[1:2, ] == 0))
for (bad in list(cbind(0, 0), cbind(0, -Inf), cbind(2, 1), matrix(1, 2, 3)))
  fails(f$prepare.response(bad))
for (bad_weights in list(rep(-1, n), rep(0, n))) fails(gam(forms, data = D, weights = bad_weights, family = f))
fails(cgen_gammalss(link = list("identity", "log", "log")))
fails(cgen_gammalss(observe = 1))
for (discrete in c(FALSE, TRUE)) for (method in c("fREML", "NCV"))
  fails(suppressWarnings(bam(form, data = D, family = f, method = method, discrete = discrete)))
