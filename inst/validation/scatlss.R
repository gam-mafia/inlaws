# Run with Rscript inst/validation/scatlss.R after installing inlaws.
library("inlaws")

# Reproduce the normalizer polynomial from Stirling's correction:
# .5*(log1p(h)/h - 1) + sum b_m*(2h)^m*((1+h)^(-m)-1).
a <- c(0, .5 * (-1)^(1:15) / (2:16))
b <- c(1/12, -1/360, 1/1260, -1/1680, 1/1188,
       -691/360360, 1/156, -3617/122400)
for (s in seq_along(b)) {
  m <- 2*s - 1
  if (m < 15) for (l in seq_len(15-m))
    a[m+l+1] <- a[m+l+1] + b[s]*2^m*(-1)^l*choose(m+l-1, l)
}
expected <- c(0, -1/4, 0, 1/24, 0, -1/20, 0, 17/112, 0, -31/36,
              0, 691/88, 0, -5461/52, 0, 929569/480)
stopifnot(max(abs(a-expected)) < 1e-10)
print(data.frame(power = 0:15, coefficient = expected))

# Independent density comparison and finite fourth-order derivative evaluation.
evaluate <- getFromNamespace(".scatlss_eval", "inlaws")
grid <- expand.grid(nu = c(3.01, 4, 40-1e-4, 40+1e-4, 1e3, 1e8, 1e100),
                    sigma = c(.03, 1, 10), z = c(0, .5, 3, 30))
eta <- cbind(.2, log(grid$sigma), log(grid$nu-3))
y <- .2 + grid$sigma*grid$z
r <- evaluate(y, eta, 4L)
reference <- dt(grid$z, grid$nu, log = TRUE) - log(grid$sigma)
error <- max(abs(r$ell[,1]-reference))
stopifnot(error < 1e-9, all(is.finite(r$ell)))
cat("Maximum density error:", format(error, digits = 5), "\n")
# Complete fitting, derivative and signed-deletion references are regression
# tests in tests/scatlss-{models,derivatives,distributions,ncv}.R.
