library(inlaws)

near <- function(x, y, tol = 1e-6) {
  stopifnot(length(x) == length(y), all(is.finite(x)), all(is.finite(y)),
            max(abs(x - y)) <= tol * max(1, abs(x), abs(y)))
}
expect_error <- function(expr, pattern) {
  e <- tryCatch(expr, error = identity)
  stopifnot(inherits(e, "error"), grepl(pattern, conditionMessage(e)))
}
set.seed(724)
n <- 450L
d <- data.frame(x = runif(n), z = runif(n), m = sample(5:60, n, TRUE),
                a = runif(n, -.3, .3), b = runif(n, -.2, .2))
d$mu <- plogis(.2 + sin(2 * pi * d$x) + d$a)
d$phi <- exp(1 + d$z + d$b)
d$y <- rbinom(n, d$m, rbeta(n, d$mu * d$phi, (1 - d$mu) * d$phi))
form <- list(cbind(y, m - y) ~ s(x, k = 6) + offset(a), ~ s(z, k = 4) + offset(b))
family <- betabinomial()
f <- gam(form, data = d, family = family, method = "REML")
g <- gam(form, data = d, family = family, method = "REML", optimizer = c("outer", "bfgs"))
for (fit in list(f, g)) {
  stopifnot(fit$outer.info$conv == "full convergence", all(is.finite(fit$Vp)),
            fit$method == "REML", length(fit$y) == n)
}
near(f$gcv.ubre, g$gcv.ubre, 1e-5)
near(f$fitted.values, g$fitted.values, .002)
stopifnot(sqrt(mean((f$fitted.values[, 1] - d$mu)^2)) < .08,
          sqrt(mean((log(f$fitted.values[, 2]) - log(d$phi))^2)) < .5)
stopifnot(is.null(family$trials)) # fitting does not mutate the supplied family

propform <- form; propform[[1]] <- y / m ~ s(x, k = 6) + offset(a)
p <- gam(propform, data = d, weights = m, family = family, method = "REML")
near(coef(f), coef(p), 1e-7)
near(f$prior.weights, rep(1, n)); near(p$prior.weights, rep(1, n))
near(f$family$trials, d$m); near(p$family$trials, d$m)
near(p$deviance, f$deviance); near(p$null.deviance, f$null.deviance)

pr <- predict(f, newdata = d[1:8, ], type = "response", se.fit = TRUE)
pl <- predict(f, newdata = d[1:8, ], type = "link", se.fit = TRUE)
near(pr$fit[, 1], plogis(pl$fit[, 1])); near(pr$fit[, 2], exp(pl$fit[, 2]))
near(pr$se.fit[, 1], pl$se.fit[, 1] * plogis(pl$fit[, 1]) * plogis(-pl$fit[, 1]))
near(pr$se.fit[, 2], pl$se.fit[, 2] * pr$fit[, 2])
near(residuals(f, type = "response"), d$y - d$m * f$fitted.values[, 1])
near(sum(residuals(f, type = "deviance")^2), f$deviance)
stopifnot(all(is.finite(residuals(f, type = "pearson"))))
# Verify the conditional saturated optimum independently over a wider interval.
sat <- getFromNamespace(".bb_saturated", "inlaws")
for (i in 1:8) {
  if (d$y[i] %in% c(0, d$m[i])) next
  phi <- f$fitted.values[i, 2]
  objective <- function(p) lchoose(d$m[i], d$y[i]) +
    lbeta(d$y[i] + p * phi, d$m[i] - d$y[i] + (1 - p) * phi) -
    lbeta(p * phi, (1 - p) * phi)
  near(sat(d$y[i], d$m[i], log(phi)), optimize(objective, c(1e-9, 1 - 1e-9),
    maximum = TRUE, tol = 1e-10)$objective)
}

# Case weights equal replication for an unpenalized model.
d$w <- sample(0:3, n, TRUE)
f0 <- gam(list(cbind(y, m - y) ~ x, ~ z), data = d, weights = w,
          family = family, method = "REML")
dr <- d[rep(seq_len(n), d$w), ]
f1 <- gam(list(cbind(y, m - y) ~ x, ~ z), data = dr, family = family, method = "REML")
near(coef(f0), coef(f1)); near(f0$l, f1$l)
# Independent likelihood optimization for the unpenalized fit.
objective <- function(b) {
  mu <- plogis(b[1] + b[2] * d$x); phi <- exp(b[3] + b[4] * d$z)
  -sum(d$w * (lchoose(d$m, d$y) + lbeta(d$y + mu * phi, d$m - d$y + (1 - mu) * phi) -
                lbeta(mu * phi, (1 - mu) * phi)))
}
opt <- optim(c(0, 0, 1, 0), objective, method = "BFGS", control = list(reltol = 1e-11))
stopifnot(opt$convergence == 0L)
near(opt$par, coef(f0), 1e-4)

# Missing values, zero weights and zero trials preserve response alignment.
dz <- d; dz$m[1:3] <- dz$y[1:3] <- 0; dz$x[4] <- NA; dz$w[5] <- 0
fz <- gam(form, data = dz, weights = w, family = family, method = "REML")
near(fz$family$trials, dz$m[-4]); near(fz$prior.weights, dz$w[-4])
stopifnot(all(fz$bb.deviance[1:3] == 0), all(is.finite(residuals(fz, type = "pearson"))))
# Same data without zero-information rows, with fixed smoothing parameters.
fdrop <- gam(form, data = dz[dz$m > 0 & dz$w > 0 & !is.na(dz$x), ], weights = w,
             family = family, method = "REML", sp = fz$sp)
near(predict(fz, newdata = d, type = "response"), predict(fdrop, newdata = d, type = "response"), .01)

# NCV optimizes both predictors. Repeated prediction indices exercise the row
# adapter, and varying trial totals would expose accidental column indexing.
nei <- list(a = as.integer(c(1, 2, 3, 4, 5, 6)), ma = c(2L, 4L, 6L),
            d = c(1L, 2L, 3L, 3L, 6L), md = c(2L, 4L, 5L))
fn <- gam(form, data = d, family = family, method = "NCV")
stopifnot(fn$method == "NCV", fn$outer.info$conv == "full convergence",
          all(is.finite(fn$Vp)), is.finite(fn$gcv.ubre))
fp <- gam(propform, data = d, weights = m, family = family, method = "NCV")
near(coef(fn), coef(fp), 1e-6)
nc <- gam(form, data = d, family = family, method = "NCV", nei = nei, sp = f$sp)
exact <- 0
for (k in 1:3) {
  a <- nei$a[seq.int(c(0, nei$ma)[k] + 1L, nei$ma[k])]
  ix <- nei$d[seq.int(c(0, nei$md)[k] + 1L, nei$md[k])]
  fit <- gam(form, data = d[-a, ], family = family, method = "REML", sp = f$sp)
  pp <- predict(fit, newdata = d[ix, ], type = "response")
  aa <- pp[, 1] * pp[, 2]; bb <- (1 - pp[, 1]) * pp[, 2]
  exact <- exact - sum(lchoose(d$m[ix], d$y[ix]) +
    lbeta(d$y[ix] + aa, d$m[ix] - d$y[ix] + bb) - lbeta(aa, bb))
}
near(nc$gcv.ubre, exact, .002)
expect_error(gam(form, data = d, family = family, method = "QNCV"), "not QNCV")

expect_error(betabinomial(c("probit", "log")), "requires links")
expect_error(gam(list(cbind(y, 1 - y) ~ 1, ~ 1), data = data.frame(y = c(0, 1)),
                 family = family, method = "REML"), "at least two trials")
for (yy in c(0, 10)) expect_error(gam(list(cbind(y, 10 - y) ~ 1, ~ 1),
  data = data.frame(y = rep(yy, 10)), family = family, method = "REML"), "all-failure or all-success")
expect_error(gam(list(y ~ 1, ~ 1), data = data.frame(y = c(.2, .35)), weights = c(3, 4),
                 family = family, method = "REML"), "implied successes")
expect_error(gam(list(cbind(y, 10 - y) ~ 1, ~ 1), data = data.frame(y = c(1.5, 3)),
                 family = family, method = "REML"), "integer counts")
expect_error(gam(form, data = d, weights = rep(-1, n), family = family,
                 method = "REML"), "negative")
expect_error(bam(form, data = d, family = family), "general families")
# Precision near the binomial boundary must be diagnosed, not capped.
bound <- data.frame(y = rep(c(4, 5, 6, 5), 30))
warnings <- character()
fb <- withCallingHandlers(gam(list(cbind(y, 10 - y) ~ 1, ~ 1), data = bound,
  family = family, method = "REML"), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  })
stopifnot(any(grepl("precision is weakly identified", warnings)))

# A term shared by both predictors exercises the fourth-order overlap fix in
# the complete REML optimizer, rather than only in likelihood contractions.
set.seed(4)
ds <- data.frame(x = runif(300), z = runif(300), w = runif(300))
u <- plogis(.3 + sin(6 * ds$x) + .2 * sin(6 * ds$w))
p <- exp(1 + .4 * ds$z + .2 * sin(6 * ds$w))
ds$y <- rbinom(300, 20, rbeta(300, u * p, (1 - u) * p))
shared <- list(cbind(y, 20 - y) ~ s(x, k = 5), ~ s(z, k = 4), 1 + 2 ~ s(w, k = 4) - 1)
sa <- gam(shared, data = ds, family = family, method = "REML")
sb <- gam(shared, data = ds, family = family, method = "REML", optimizer = c("outer", "bfgs"))
sc <- gam(shared, data = ds, family = family, method = "NCV")
stopifnot(sa$outer.info$conv == "full convergence", sb$outer.info$conv == "full convergence",
          sc$outer.info$conv == "full convergence")
near(sa$gcv.ubre, sb$gcv.ubre, 1e-5)
sstart <- gam(shared, data = ds, family = family, method = "REML", sp = sa$sp, start = coef(sa))
near(coef(sstart), coef(sa), 1e-5)
