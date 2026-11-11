library("inlaws")

f <- bgevlss()
eval_density <- getFromNamespace(".bgev_eval", "inlaws")
eta <- c(.6, log(1.3), log(.27))
v <- function(p) -log(-log(p))
gev_q <- function(p, e) {
  xi <- exp(e[3L])
  e[1L] + exp(e[2L])*expm1(xi*(v(p) - v(.5))) /
    (exp(xi*(v(.25) - v(.5)))*expm1(xi*(v(.75) - v(.25))))
}
jet <- function(y, eta, order) eval_density(y, matrix(eta, 1), order, f$spec)
near <- function(a, b, tol) stopifnot(all(abs(a - b) <= tol*(1 + abs(b))))

# At ordinary left, blend and right points, check all mixed derivatives by
# differentiating the preceding Taylor order with a five-point stencil.
# Also differentiate the CDF in y to check the likelihood's density.
pars <- matrix(c(eta[1], exp(eta[2:3])), 1)
for (y in gev_q(c(.02, .09, .15, .4, .99), eta)) {
  h <- 1e-5
  dd <- (f$cdf(y + h, pars) - f$cdf(y - h, pars))/(2*h)
  near(exp(jet(y, eta, 0)$ell), dd, 2e-9)
  full <- jet(y, eta, 4L)
  for (order in 0:3) {
    powers <- jet(y, eta, order)$powers
    for (k in 1:3) {
      step <- numeric(3)
      step[k] <- 2e-5
      fd <- (jet(y, eta - 2*step, order)$ell - 8*jet(y, eta - step, order)$ell +
        8*jet(y, eta + step, order)$ell - jet(y, eta + 2*step, order)$ell)/(12*step[k])
      degree <- powers
      degree[, k] <- degree[, k] + 1L
      ix <- match(apply(degree, 1L, paste, collapse = ""), apply(full$powers, 1L, paste, collapse = ""))
      near(fd, full$ell[, ix]*degree[, k], 1e-7)
    }
  }
}

# The two joins move with all three predictors. At each boundary, approach
# from both branches along each predictor and compare the full fourth-order
# Taylor jets with the boundary value. Also difference the third-order jets
# across the boundary to verify the fourth-order terms at the join.
for (p in c(.05, .20)) {
  y <- gev_q(p, eta)
  near(f$cdf(y, pars), p, 1e-13)
  at <- jet(y, eta, 4L)
  for (k in 1:3) {
    step <- numeric(3)
    step[k] <- 1e-6
    minus <- jet(y, eta - step, 4L)$ell
    plus <- jet(y, eta + step, 4L)$ell
    near(2*jet(y, eta - step/2, 4L)$ell - minus, at$ell, 2e-6)
    near(2*jet(y, eta + step/2, 4L)$ell - plus, at$ell, 2e-6)
    lo <- jet(y, eta - step, 3L)
    hi <- jet(y, eta + step, 3L)
    degree <- lo$powers
    degree[, k] <- degree[, k] + 1L
    ix <- match(apply(degree, 1L, paste, collapse = ""), apply(at$powers, 1L, paste, collapse = ""))
    half <- jet(y, eta + step/2, 3L)$ell - jet(y, eta - step/2, 3L)$ell
    fd <- (2*half - (hi$ell - lo$ell)/2)/step[k]
    near(fd, at$ell[, ix]*degree[, k], 2e-6)
  }
}

# The log-shape derivatives vanish continuously as xi tends to zero;
# location/spread derivatives remain those of the Gumbel likelihood.
y <- gev_q(c(.015, .1, .4, .95), c(eta[1:2], log(1e-10)))
small <- eval_density(y, cbind(eta[1], eta[2], rep(log(1e-10), length(y))), 4L, f$spec)
gumbel <- eval_density(y, cbind(eta[1], eta[2], rep(-Inf, length(y))), 4L, f$spec)
sigma <- exp(eta[2])/(v(.75) - v(.25))
z <- (y - eta[1])/sigma + v(.5)
near(gumbel$ell[, 1L], -z - exp(-z) - log(sigma), 2e-14)
near(small$ell, gumbel$ell, 8e-8)
stopifnot(any(small$ell[, small$powers[, 3L] > 0] != 0))
