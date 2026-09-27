library(inlaws)

near <- function(a, b, tol = 1e-6) {
  stopifnot(identical(dim(a), dim(b)),
            all(is.finite(a)), all(is.finite(b)),
            max(abs(a - b) / (1 + abs(b))) < tol)
}
fails <- function(expr, pattern) {
  z <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(z, "error"), grepl(pattern, conditionMessage(z)))
}
get_internal <- function(x) getFromNamespace(x, "inlaws")
derivatives <- get_internal(".dirichlet_derivatives")
jet_indices <- get_internal(".dirichlet_jets")
parameters <- get_internal(".dirichlet_parameters")

# Independently differentiate each preceding order. This covers every mixed
# partial, repeated index, concentration cross-derivative and packing position.
set.seed(305)
for (q in 2:5) {
  jets <- jet_indices(q)
  eta <- matrix(rnorm(q), 1L); eta[, q] <- 1.4
  y <- matrix(seq_len(q) / sum(seq_len(q)), 1L)
  z <- derivatives(y, eta, 2.3, jets, 4L)
  for (r in 1:4) {
    indices <- which(jets$degree == r)
    for (j in seq_along(indices)) {
      a <- jets$powers[indices[j], ]
      k <- which(a > 0)[1L]
      a[k] <- a[k] - 1L
      h <- 1e-4
      ep <- em <- eta; ep[, k] <- ep[, k] + h; em[, k] <- em[, k] - h
      low <- function(e) {
        zz <- derivatives(y, e, 2.3, jets, r - 1L)
        if (r == 1L) return(zz$l)
        idx <- which(jets$degree == r - 1L)
        pos <- which(apply(jets$powers[idx, , drop = FALSE], 1L,
                           function(b) all(a == b)))
        zz[[paste0("l", r - 1L)]][1L, pos]
      }
      near(as.numeric(z[[paste0("l", r)]][1L, j]),
           (low(ep) - low(em)) / (2 * h), 2e-6)
    }
  }
}

# Independent density and two-component beta likelihood, including constants.
y <- rbind(c(.2, .3, .5), c(.1, .6, .3))
e <- rbind(c(.5, -.7, 2), c(-.2, .3, 1))
a <- parameters(e)$mu * parameters(e)$phi
reference <- sum(c(2, .5) * (lgamma(rowSums(a)) - rowSums(lgamma(a)) +
                             rowSums((a - 1) * log(y))))
near(derivatives(y, e, c(2, .5), jet_indices(3), 0)$l, reference)
p <- parameters(matrix(c(.4, 2), 1))
near(derivatives(matrix(c(.3, .7), 1), matrix(c(.4, 2), 1), 1,
                 jet_indices(2), 0)$l,
     dbeta(.7, p$mu[2] * p$phi, p$mu[1] * p$phi, log = TRUE))

# Finite derivatives near the boundary and over a wide concentration range.
for (phi in c(.001, 1, 1e4)) {
  e <- matrix(c(-12, 2, log(phi)), 1)
  z <- derivatives(matrix(c(1e-8, .3, .7 - 1e-8), 1), e, 1, jet_indices(3), 4)
  stopifnot(all(is.finite(unlist(z))))
}
z <- derivatives(y, rbind(c(.5, -.7, 2), c(1000, -1000, 1000)),
                 c(1, 0), jet_indices(3), 4)
stopifnot(all(is.finite(unlist(z))), all(z$l1[2, ] == 0))
stopifnot(identical(derivatives(y[1,,drop=FALSE], matrix(c(0, 0, 1000), 1),
                              1, jet_indices(3), 0)$l, -Inf))

# Smooth composition and concentration, names, standard errors and recovery.
set.seed(512)
n <- 500
x <- runif(n); z <- runif(n)
eta <- cbind(.6 * sin(2 * pi * x), -.5 * x, 2.6 + .7 * z)
y <- dirichlet(2)$rd(eta, rep(1, n), 1)
colnames(y) <- c("a", "b", "c")
dat <- data.frame(x, z); dat$y <- y
b <- gam(list(y ~ s(x, k = 6), ~ s(x, k = 5), ~ s(z, k = 4)),
         data = dat, family = dirichlet(2), method = "REML")
stopifnot(identical(b$outer.info$conv, "full convergence"), is.finite(as.numeric(logLik(b))), is.na(b$deviance),
          is.na(b$null.deviance), all(is.finite(vcov(b))))
pr <- predict(b, type = "response", se.fit = TRUE)
stopifnot(identical(b$family$data$component.names, colnames(y)), all(pr$fit > 0),
          max(abs(rowSums(pr$fit) - 1)) < 1e-14,
          all(pr$se.fit > 0), all(is.finite(pr$se.fit)),
          sqrt(mean((pr$fit - parameters(eta)$mu)^2)) < .06,
          sqrt(mean((predict(b, type = "link")[, 3] - eta[, 3])^2)) < .4)

# Numerical prediction Jacobian against the entire coefficient covariance.
nd <- dat[1:3, ]
X <- predict(b, newdata = nd, type = "lpmatrix")
fam <- b$family; beta <- coef(b)
J <- matrix(0, 3 * 3, length(beta))
for (j in seq_along(beta)) {
  bp <- bm <- beta; bp[j] <- bp[j] + 1e-5; bm[j] <- bm[j] - 1e-5
  J[, j] <- as.vector((fam$predict(fam, X = X, beta = bp)$fit -
                       fam$predict(fam, X = X, beta = bm)$fit) / 2e-5)
}
se <- matrix(sqrt(rowSums((J %*% vcov(b)) * J)), 3, 3)
near(unname(predict(b, newdata = nd, type = "response", se.fit = TRUE)$se.fit), se)
near(unname(predict(b, newdata = nd[1,,drop=FALSE], type = "response")),
     unname(pr$fit[1,,drop=FALSE]))

# Coefficient MLE agrees with an independent optimization of the density.
u <- gam(list(y ~ x, ~ x, ~ z), data = dat, family = dirichlet(2), method = "REML")
Xm <- model.matrix(~ x, dat); Xp <- model.matrix(~ z, dat)
objective <- function(beta) {
  p <- parameters(cbind(Xm %*% beta[1:2], Xm %*% beta[3:4], Xp %*% beta[5:6]))
  a <- p$mu * p$phi
  -sum(lgamma(rowSums(a)) - rowSums(lgamma(a)) + rowSums((a - 1) * log(y)))
}
opt <- optim(coef(u) + .05, objective, method = "BFGS",
             control = list(reltol = 1e-11, maxit = 500))
stopifnot(opt$convergence == 0)
near(unname(coef(u)), unname(opt$par), 2e-4)
near(as.numeric(logLik(u)), -opt$value, 1e-7)

# Beta reduction with a single mean predictor and estimated constant precision.
bdat <- data.frame(x = x, p = y[, 2] / (y[, 1] + y[, 2]))
bdat$yy <- cbind(1 - bdat$p, bdat$p)
bb <- gam(p ~ x, data = bdat, family = betar(), method = "ML")
bd <- gam(list(yy ~ x, ~ 1), data = bdat, family = dirichlet(), method = "REML")
near(as.numeric(predict(bb, type = "response")),
     as.numeric(predict(bd, type = "response")[, 2]), 2e-4)
near(as.numeric(exp(coef(bd)[3])), bb$family$getTheta(TRUE), 2e-3)

# Integer likelihood weights equal row replication; zero weights equal omission.
dat$w <- rep(c(0, 1, 2), length.out = n)
w <- gam(list(y ~ x, ~ x, ~ 1), data = dat, weights = w,
         family = dirichlet(2), method = "REML")
dr <- dat[rep(seq_len(n), dat$w), ]
r <- gam(list(y ~ x, ~ x, ~ 1), data = dr, family = dirichlet(2), method = "REML")
near(unname(coef(w)), unname(coef(r)), 1e-6)
near(as.numeric(logLik(w)), as.numeric(logLik(r)), 1e-7)

# Offsets on every predictor; shared terms and their joint initialization.
dat$o1 <- .2 * x; dat$o2 <- -.3 * x; dat$o3 <- .4 * z
of <- gam(list(y ~ x + offset(o1), ~ x + offset(o2), ~ z + offset(o3)),
          data = dat, family = dirichlet(2), method = "REML")
near(unname(predict(of, type = "response")), unname(predict(u, type = "response")), 1e-5)
shared <- gam(list(y ~ 1, ~ 1, ~ 1, 1 + 2 ~ s(x, k = 5) - 1),
              data = dat, family = dirichlet(2), method = "REML")
stopifnot(all(is.finite(coef(shared))), all(is.finite(vcov(shared, sandwich = TRUE))))
ps <- predict(shared, newdata = nd, type = "response", se.fit = TRUE)
stopifnot(identical(dim(ps$fit), c(3L, 3L)), all(is.finite(ps$se.fit)))

# Missing response components remove whole rows; na.exclude restores rows.
dm <- dat; dm$y[2, 2] <- NA; dm$x[3] <- NA
bm <- gam(list(y ~ x, ~ x, ~ z), data = dm, family = dirichlet(2),
          method = "REML", na.action = na.exclude)
stopifnot(nrow(bm$y) == n - 2L, nrow(predict(bm, type = "response")) == n,
          all(is.na(predict(bm, type = "response")[2:3, ])))

# The rd callback consumes link-scale predictors and returns a simplex matrix.
set.seed(31)
sim <- b$family$rd(fitted(b), NULL, 1)
stopifnot(identical(dim(sim), c(as.integer(n), 3L)),
          max(abs(rowSums(sim) - 1)) < 1e-14)
set.seed(31)
near(sim, b$family$rd(predict(b, type = "link"), NULL, 1))
ss <- shared$family$rd(predict(shared, nd, type = "link"), NULL, 1)
stopifnot(identical(dim(ss), c(3L, 3L)), max(abs(rowSums(ss) - 1)) < 1e-14)
set.seed(81)
et <- matrix(rep(c(.3, -.5, log(10)), each = 15000), 15000)
draw <- dirichlet(2)$rd(et, NULL, 1)
pm <- parameters(et)$mu[1, ]
near(colMeans(draw), pm, .006)
near(apply(draw, 2, var), pm * (1 - pm) / 11, .002)

# Residuals and simulated QQ diagnostics use the same scalar discrepancy.
near(unname(residuals(b, type = "response")), unname(y - pr$fit))
stopifnot(length(residuals(b, type = "pearson")) == n)
pdf(file = tempfile(fileext = ".pdf"))
qq.gam(b, type = "pearson", rep = 5)
dev.off()
fails(residuals(b), "deviance residuals are undefined")

# Validation errors are informative and happen before fitting.
for (bad in list(0, -1, 1.5, NA_real_, Inf, c(1, 2), "2"))
  fails(dirichlet(bad), "positive integer")
fitbad <- function(yy, ww = rep(1, nrow(yy))) {
  dd <- data.frame(x = seq_len(nrow(yy)), w = ww); dd$y <- yy
  gam(list(y ~ 1, ~ 1, ~ 1), data = dd, weights = w,
      family = dirichlet(2), method = "REML")
}
yy <- y; yy[1, ] <- c(0, .5, .5); fails(fitbad(yy), "strictly between")
yy <- y; yy[1, 1] <- Inf; fails(fitbad(yy), "finite")
fails(fitbad(y * .9), "sum to one")
fails(fitbad(y[, 1:2]), "K \\+ 1 columns")
fails(fitbad(y, rep(0, n)), "positive weight")
fails(fitbad(y, rep(-1, n)), "negative|nonnegative")
fails(gam(list(y ~ 1, ~ 1), family = dirichlet(2), data = dat), "formula|predictor")
fails(bam(list(y ~ x, ~ x, ~ 1), family = dirichlet(2), data = dat),
      "general families not supported")
# Rounding-level discrepancies are normalized, without changing observations otherwise.
bn <- fitbad(y * (1 + 1e-10))
near(unname(bn$y), unname(y))
