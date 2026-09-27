library(inlaws)
near <- function(x, y, tol = 1e-5) stopifnot(length(x) == length(y), all(is.finite(x)), all(is.finite(y)), max(abs(x - y) / pmax(1, abs(x), abs(y))) < tol)
ns <- asNamespace("mgcv")
setup <- get("Sl.setup", ns)
repara <- get("Sl.initial.repara", ns)
fit5 <- get("gam.fit5", ns)
set.seed(908)
n <- 450
d <- data.frame(x = runif(n), z = runif(n), a = runif(n, -.2, .2), b = runif(n, -.2, .2), c = runif(n, -.2, .2))
mu <- cbind(exp(1 + .5 * sin(6 * d$x)), plogis(.3 + .4 * d$x), exp(.5 + .2 * d$z))
# Reordered and repeated prediction rows, with grouped deletions and m != n.
nei <- list(
  a = c(4, 5, 6, 10, 11, 12, 4, 5, 6, 19, 20, 21), ma = c(3, 6, 9, 12),
  d = c(5, 11, 12, 5, 20), md = c(1, 3, 4, 5), jackknife = FALSE
)
for (h in c(FALSE, TRUE)) {
  for (shared in c(FALSE, TRUE)) {
    fam <- if (h) zanb() else zinb()
    d$y <- fam$rd(mu)
    form <- if (shared) {
      list(y ~ 1 + offset(a), ~ x + offset(b), ~ z + offset(c), 1 + 3 ~ s(x, k = 5) - 1)
    } else {
      list(y ~ s(x, k = 5) + offset(a), ~ s(x, k = 4) + offset(b), ~ s(z, k = 4) + offset(c))
    }
    G <- gam(form, data = d, family = fam, fit = FALSE)
    Sl <- setup(G)
    X <- repara(Sl, G$X, both.sides = FALSE)
    rho <- rep(log(2), length(G$sp))
    start <- NULL
    call <- function(rho, method = "REML", deriv = 1, q = FALSE) {
      ff <- fam
      ff$qapprox <- q
      fit5(X, G$y, rho, Sl,
        weights = G$w, offset = G$offset, family = ff,
        deriv = deriv, scoreType = method, Mp = -1, nei = nei, gamma = 1.2, start = start, control = gam.control(epsilon = 1e-9)
      )
    }
    base <- call(rho, deriv = 2)
    start <- base$coefficients
    # A five-point stencil avoids amplifying inner-fit stopping noise when
    # differentiating the Laplace determinant of a mixture fit.
    for (j in seq_along(rho)) {
      step <- rep(0, length(rho)); step[j] <- .01
      fmm <- call(rho - 2*step); fm <- call(rho - step)
      fp <- call(rho + step); fpp <- call(rho + 2*step)
      near(base$REML1[j], (fmm$REML-8*fm$REML+8*fp$REML-fpp$REML)/.12, 5e-5)
      near(base$REML2[,j], (fmm$REML1-8*fm$REML1+8*fp$REML1-fpp$REML1)/.12, 2e-4)
    }
    for (q in c(FALSE, TRUE)) {
      f <- call(rho, "NCV", q = q)
      for (j in seq_along(rho)) {
        rp <- rm <- rho
        rp[j] <- rp[j] + 1e-4
        rm[j] <- rm[j] - 1e-4
        near(f$NCV1[j], (call(rp, "NCV", q = q)$NCV - call(rm, "NCV", q = q)$NCV) / 2e-4, 2e-4)
      }
      if (!q) {
        # Actual refits use the same full-data basis and fixed penalties.
        cv <- attr(f$NCV, "eta.cv")
        exact <- cv * 0
        prevA <- prevD <- 0
        for (k in seq_along(nei$ma)) {
          dropped <- nei$a[(prevA + 1):nei$ma[k]]
          pred <- nei$d[(prevD + 1):nei$md[k]]
          w <- G$w
          w[dropped] <- 0
          ff <- fit5(X, G$y, rho, Sl,
            weights = w, offset = G$offset, family = fam,
            deriv = 0, Mp = -1, gamma = 1.2, start = start, control = gam.control(epsilon = 1e-9)
          )
          exact[(prevD + 1):nei$md[k], ] <- ff$linear.predictors[pred, ]
          prevA <- nei$ma[k]
          prevD <- nei$md[k]
        }
        near(cv, exact, .015)
      }
    }
  }
}
