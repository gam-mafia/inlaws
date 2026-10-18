library("inlaws")
near <- function(x, y, tol = 1e-7) {
  stopifnot(identical(dim(x), dim(y)), length(x) == length(y),
            all(is.finite(x)), all(is.finite(y)),
            max(abs(x - y) / pmax(1, abs(x), abs(y))) < tol)
}
prediction <- getFromNamespace(".zinb_prediction", "inlaws")
columns <- c("mu", "p", "theta", "mean", "variance", "p0")

# The audit deliberately retains mgcv's default for these constructors.
defaults <- list(betabinomial, betals, cgamma, cgammals, clognormal,
                 clognormalls, cmp, gpoisson, log_logistic, log_logisticls,
                 nbls, pig)
stopifnot(all(vapply(defaults, function(f) is.null(f()$predict), logical(1))),
          is.function(cgen_gammalss()$predict))

# Independent PMF sums, not the package's residual/diagnostic moment formulas.
for (h in c(FALSE, TRUE)) {
  fam <- if (h) zanb() else zinb()
  for (par in list(c(.001, .3, .1), c(3, .7, 2), c(20, .9, 10))) {
    eta <- matrix(c(log(par[1]), qlogis(par[2]), log(par[3])), 1)
    z <- prediction(eta, h)
    yy <- 0:5000
    nb <- dnbinom(yy, mu = par[1], size = par[3])
    pmf <- par[2] * nb / if (h) -expm1(log(nb[1])) else 1
    pmf[1] <- 1 - par[2] + if (h) 0 else par[2] * nb[1]
    mn <- sum(yy * pmf)
    near(unname(z$fit[1, 4:6]), c(mn, sum((yy - mn)^2 * pmf), pmf[1]))
    stopifnot(identical(colnames(z$fit), columns))
    near(fam$predict(fam, eta = eta)$fit, z$fit)
  }
  # Independent numerical differentiation of each response column, including
  # both series boundaries, small sizes, large sizes and extreme mixture logits.
  eta <- as.matrix(expand.grid(u = c(-12, log(.00999), log(.01001), 0, 3),
                              b = c(-25, -.3, 25), v = c(-3, 0, 20)))
  z <- prediction(eta, h)
  for (j in 1:3) {
    ep <- em <- eta; ep[, j] <- ep[, j] + 1e-5; em[, j] <- em[, j] - 1e-5
    num <- (prediction(ep, h)$fit - prediction(em, h)$fit) / 2e-5
    near(z$gradient[, , j], num, 2e-6)
  }
  # Tiny-mean hurdle law tends to Bernoulli(p), including its uncertainty.
  tiny <- matrix(c(log(1e-16), qlogis(.3), log(2)), 1)
  tz <- prediction(tiny, h)
  if (h) {
    near(unname(tz$fit[1, 4:6]), c(.3, .21, .7), 1e-14)
    # Positive component variance is about mu*(1+1/theta)/2, not zero.
    tiny[1, 2] <- 40
    tt <- prediction(tiny, h)
    near(tt$fit[1, 5] / (7.5e-17 + plogis(-40)), 1, 1e-13)
    near(tt$gradient[1, 4, 1] / 7.5e-17, 1, 1e-13)
    near(tt$gradient[1, 4, 3] / -2.5e-17, 1, 1e-13)
    # The same limit when mu/theta underflows.
    tiny[1, c(1, 3)] <- c(log(1e-200), log(1e200))
    stopifnot(all(is.finite(prediction(tiny, h)$fit)))
  }
  # Near-Poisson NB component and rounded endpoint mixture probabilities.
  eta <- rbind(c(log(2), qlogis(.7), log(1e12)), c(0, -40, 0), c(0, 40, 0))
  zz <- prediction(eta, h)
  pos <- if (h) -expm1(-2) else 1
  mn <- .7 * 2 / pos
  near(unname(zz$fit[1, 4:6]), c(mn, .7 * (2 + 4) / pos - mn^2,
                                            .3 + if (h) 0 else .7 * exp(-2)), 1e-11)
  stopifnot(all(is.finite(zz$gradient)), zz$gradient[3, 2, 2] > 0)

  # Native contract with correlated coefficients and overlapping lpi blocks.
  X <- cbind(1, c(-.5, .3, .8), c(.2, -.4, .1), 1)
  attr(X, "lpi") <- list(1:2, 2:3, c(2L, 4L))
  beta <- c(.4, .2, -.3, .7)
  off <- list(c(.1, -.2, .05), .15, c(-.1, .1, .2))
  A <- matrix(c(1, .4, -.3, .2, .2, 1, .5, -.4, -.1, .3, 1, .2, .4, -.2, .3, 1), 4)
  V <- crossprod(A) / 20
  pr <- fam$predict(fam, TRUE, X = X, beta = beta, off = off, Vb = V)
  stopifnot(identical(colnames(pr$se.fit), columns))
  J <- array(0, c(3, 6, 4))
  for (j in seq_along(beta)) {
    bp <- bm <- beta; bp[j] <- bp[j] + 1e-5; bm[j] <- bm[j] - 1e-5
    J[, , j] <- (fam$predict(fam, X = X, beta = bp, off = off)$fit -
                  fam$predict(fam, X = X, beta = bm, off = off)$fit) / 2e-5
  }
  for (k in 1:6) near(pr$se.fit[, k], sqrt(rowSums((J[, k, ] %*% V) * J[, k, ])), 2e-7)
}

# Integration with predict.gam: offsets in every predictor, prediction blocks,
# omitted response, one row, missing rows, exclusion and unconditional Vb.
set.seed(180)
dat <- data.frame(x = runif(300, -1, 1), z = runif(300, -1, 1))
dat$o1 <- .1 * dat$x; dat$o2 <- -.2 * dat$z; dat$o3 <- .1 * dat$z
for (h in c(FALSE, TRUE)) {
  fam <- if (h) zanb() else zinb()
  dat$y <- fam$rd(cbind(exp(.7 + dat$x), plogis(.8 + .2 * dat$z), exp(1 + .2 * dat$z)))
  b <- gam(list(y ~ s(x, k = 5) + offset(o1), ~ z + offset(o2), ~ z + offset(o3)),
           data = dat, family = fam, method = "REML")
  nd <- dat[c(7, 2, 9, 1, 3), setdiff(names(dat), "y")]
  pr <- predict(b, nd, type = "response", se.fit = TRUE)
  eta <- predict(b, nd, type = "link")
  near(unname(pr$fit), unname(prediction(eta, h)$fit))
  near(pr$fit[, 1:3], unname(cbind(exp(eta[, 1]), plogis(eta[, 2]), exp(eta[, 3]))))
  for (bs in c(1, 2, 100)) {
    pp <- predict(b, nd, type = "response", se.fit = TRUE, block.size = bs)
    near(pp$fit, pr$fit); near(pp$se.fit, pr$se.fit)
  }
  one <- predict(b, nd[1, , drop = FALSE], type = "response", se.fit = TRUE)
  stopifnot(identical(dim(one$fit), c(1L, 6L)), identical(dim(one$se.fit), c(1L, 6L)))
  for (uc in c(FALSE, TRUE)) {
    pp <- predict(b, nd, type = "response", se.fit = TRUE, unconditional = uc)
    V <- vcov(b, unconditional = uc)
    J <- array(0, c(nrow(nd), 6, length(coef(b))))
    for (j in seq_along(coef(b))) {
      bp <- bm <- b; bp$coefficients[j] <- coef(b)[j] + 1e-5; bm$coefficients[j] <- coef(b)[j] - 1e-5
      J[, , j] <- (predict(bp, nd, type = "response") - predict(bm, nd, type = "response")) / 2e-5
    }
    for (k in 1:6) near(pp$se.fit[, k], sqrt(rowSums((J[, k, ] %*% V) * J[, k, ])), 2e-6)
  }
  # Ensure the new covariance calculation matters for these fitted moments.
  Jm <- J[, 4, ]
  stopifnot(max(abs(sqrt(rowSums((Jm %*% V) * Jm)) - sqrt(rowSums(Jm^2 *
                matrix(diag(V), nrow(nd), ncol(V), byrow = TRUE))))) > 1e-5)
  ex <- predict(b, nd, type = "response", exclude = "s(x)")
  le <- predict(b, nd, type = "link", exclude = "s(x)")
  near(unname(ex), unname(prediction(le, h)$fit))
  nd$x[2] <- NA
  missing <- predict(b, nd, type = "response", se.fit = TRUE, block.size = 1)
  stopifnot(all(is.na(missing$fit[2, ])), all(is.na(missing$se.fit[2, ])))
  near(missing$fit[-2, ], pr$fit[-2, ])
  dm <- dat; dm$y[2] <- NA; dm$x[3] <- NA
  bm <- gam(list(y ~ x + offset(o1), ~ z + offset(o2), ~ z + offset(o3)),
            data = dm, family = fam, method = "REML", na.action = na.exclude)
  pm <- predict(bm, type = "response", se.fit = TRUE)
  stopifnot(nrow(pm$fit) == nrow(dm), all(is.na(pm$fit[2:3, ])),
            all(is.na(pm$se.fit[2:3, ])), ncol(fitted(bm)) == 3L)
  # Real shared-coefficient fit, not just a synthetic callback design.
  shared <- gam(list(y ~ 1 + offset(o1), ~ z + offset(o2), ~ z + offset(o3), 1+3 ~ x-1),
                data = dat, family = fam, method = "REML")
  sp <- predict(shared, nd[-2, ], type = "response", se.fit = TRUE)
  stopifnot(all(is.finite(sp$se.fit)))
  for (k in 4:6) {
    J <- matrix(0, 4, length(coef(shared)))
    for (j in seq_along(coef(shared))) {
      bp <- bm <- shared; bp$coefficients[j] <- coef(shared)[j] + 1e-5; bm$coefficients[j] <- coef(shared)[j] - 1e-5
      J[, j] <- (predict(bp, nd[-2, ], type = "response")[, k] -
                   predict(bm, nd[-2, ], type = "response")[, k]) / 2e-5
    }
    near(sp$se.fit[, k], sqrt(rowSums((J %*% vcov(shared)) * J)), 2e-6)
  }
}
