library("inlaws")
evaluate <- getFromNamespace(".cgammals_eval", "inlaws")
response <- getFromNamespace(".cgammals_response", "inlaws")
near <- function(x, y, tol = 1e-5) {
  err <- max(abs(x-y)/(1+abs(y)))
  if (!is.finite(err) || err > tol) stop("relative error: ", err)
}
take <- function(r, a, b) r$ell[, which(r$powers[,1] == a & r$powers[,2] == b)]*factorial(a)*factorial(b)
# Numerical derivatives are test oracles only. Include recurrence transitions,
# small and large shapes, both extreme tails, support endpoints and narrow bins.
for (shape in c(.03, .2, 1, 10, 100, 1000)) {
  bound <- c(.001, .1, .9, 1, 1.1, 2, 10, 100)
  y <- rbind(cbind(bound, -Inf), cbind(bound, Inf),
             cbind(bound, bound*(1+1e-9)), cbind(0, .2), cbind(0, Inf))
  u <- rep(0, nrow(y)); t <- rep(-log(shape), nrow(y))
  r <- evaluate(y, u, t, 4L)
  h <- 2e-5
  ap <- evaluate(y, u+h, t, 4L); am <- evaluate(y, u-h, t, 4L)
  tp <- evaluate(y, u, t+h, 4L); tm <- evaluate(y, u, t-h, 4L)
  for (degree in 1:4) for (b in 0:degree) {
    a <- degree-b
    if (a) near(take(r, a, b), (take(ap, a-1, b)-take(am, a-1, b))/(2*h), 2e-4)
    if (b) near(take(r, a, b), (take(tp, a, b-1)-take(tm, a, b-1))/(2*h), 2e-4)
  }
  near(r$ell[1:8,1], pgamma(bound, shape, rate=shape, log.p=TRUE), 1e-12)
  near(r$ell[9:16,1], pgamma(bound, shape, rate=shape, lower.tail=FALSE, log.p=TRUE), 1e-12)
  stopifnot(all(r$ell[nrow(y),] == 0))
}

# Independently integrate density derivatives for dispersion through order four.
# Here the derivative polynomials are assembled from the analytic gamma density,
# independently of the incomplete-gamma recurrence and Taylor arithmetic.
for (shape in c(.7, 2, 20)) for (bounds in list(c(0, .2), c(.3, 2), c(4, Inf))) {
  mu <- 1.3
  yy <- matrix(bounds, 1L)
  r <- evaluate(yy, log(mu), -log(shape), 4L)
  lp <- r$ell[1,1]
  density_derivative <- function(y, degree) {
    z <- y/mu
    A <- log(shape)+1-digamma(shape)+log(z)-z
    s1 <- -shape*A
    s2 <- -s1 + shape - shape^2*trigamma(shape)
    s3 <- -s2 - shape + 2*shape^2*trigamma(shape) + shape^3*psigamma(shape, 2)
    s4 <- -s3 + shape - 4*shape^2*trigamma(shape) -
      5*shape^3*psigamma(shape, 2) - shape^4*psigamma(shape, 3)
    bell <- switch(degree, s1, s2+s1^2, s3+3*s1*s2+s1^3,
                   s4+4*s1*s3+3*s2^2+6*s1^2*s2+s1^4)
    exp(dgamma(y, shape, scale=mu/shape, log=TRUE)-lp)*bell
  }
  d <- vapply(1:4, function(k) integrate(density_derivative, bounds[1], bounds[2],
                degree=k, rel.tol=1e-8, subdivisions=1000)$value, numeric(1))
  expected <- c(d[1], d[2]-d[1]^2, d[3]-3*d[1]*d[2]+2*d[1]^3,
                d[4]-4*d[1]*d[3]-3*d[2]^2+12*d[1]^2*d[2]-6*d[1]^4)
  near(vapply(1:4, function(k) take(r,0,k), numeric(1)), expected, 1e-6)
}

# Response contract and zero-weight safety, even with an extreme unused row.
fails <- function(expr) stopifnot(inherits(tryCatch(expr, error=identity), "error"))
for (y in list(0, -1, cbind(0,-Inf), cbind(2,1), cbind(NA,Inf), matrix(1,2,3)))
  fails(response(y))
fam <- cgammals()
y <- rbind(c(1,1), c(.Machine$double.xmax, Inf))
X <- matrix(1, 2, 2)
l <- fam$ll(y, X, c(0,0), c(1,0), fam, deriv=4, ncv=TRUE, eta=matrix(0,2,2))
stopifnot(l$l0[2] == 0, all(l$l1[2,] == 0), all(l$l4[2,] == 0))

# Near the mean, a large shape needs more than 10,000 series terms.
shape <- 2e6
lp <- stats::pgamma(shape, shape, log.p = TRUE)
expected <- -exp(log(shape) + stats::dgamma(shape, shape, log = TRUE) - lp)
large <- evaluate(cbind(1, -Inf), 0, -log(shape), 2L)
near(take(large, 1, 0), expected, 1e-7)
large <- cgamma()$Dd(cbind(1, -Inf), 1, -log(shape), 1, 2)
near(large$Dmu, -2 * expected, 1e-7)

# For integer shape, the upper tail is a finite Poisson sum. Its scaled
# terms are proportional to mu^j, so their index cumulants give the
# derivatives of the log sum without subtracting large raw boundary ratios.
x <- shape + 20 * sqrt(shape)
yy <- cbind(x / shape, Inf)
index <- 0:3000
term <- c(1, cumprod((shape - index[-1]) / x))
stopifnot(tail(term, 1) < 1e-18)
prob <- term / sum(term)
m <- sum(index * prob)
centred <- index - m
v <- sum(centred^2 * prob)
expected <- c(x - (shape - 1) + m, -x + v,
              x + sum(centred^3 * prob), -x + sum(centred^4 * prob) - 3 * v^2)
upper <- evaluate(yy, 0, -log(shape), 4L)
near(vapply(1:4, function(a) take(upper, a, 0), numeric(1)), expected, 1e-8)
d <- cgamma()$Dd(yy, 1, -log(shape), 1, 2)
partials <- unlist(d[c("Dmu", "Dmu2", "Dmu3", "Dmu4")])
change_coordinates <- rbind(c(1, 0, 0, 0), c(1, 1, 0, 0),
                            c(1, 3, 1, 0), c(1, 7, 6, 1))
near(drop(change_coordinates %*% partials) / -2, expected, 1e-8)

# A short interval can be nearly flat in density even when shape times its
# relative width is large. Integrate the differentiated density independently.
shape <- 1e7
yy <- cbind(1, 1 + .1 / sqrt(shape))
moments <- vapply(0:4, function(j) {
  integrate(function(z) {
    x <- shape + sqrt(shape) * z
    a <- x - shape; b <- -x; c <- x
    density_ratios <- list(1, a, b + a^2, c + 3*a*b + a^3,
                           -x + 4*a*c + 3*b^2 + 6*a^2*b + a^4)
    dgamma(x, shape) * density_ratios[[j + 1L]]
  }, 0, (yy[2] - 1) * sqrt(shape), rel.tol = 1e-10)$value
}, numeric(1))
ratios <- moments[-1] / moments[1]
a <- ratios[1]; b <- ratios[2]; c <- ratios[3]
expected <- c(a, b-a^2, c-3*a*b+2*a^3,
              ratios[4]-4*a*c-3*b^2+12*a^2*b-6*a^4)
interval <- evaluate(yy, 0, -log(shape), 4L)
near(vapply(1:4, function(a) take(interval, a, 0), numeric(1)), expected, 1e-8)
d <- cgamma()$Dd(yy, 1, -log(shape), 1, 2)
partials <- unlist(d[c("Dmu", "Dmu2", "Dmu3", "Dmu4")])
near(drop(change_coordinates %*% partials) / -2, expected, 1e-8)
