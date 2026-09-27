library(inlaws)
near <- function(x, y, tol = 1e-6) {
  stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol)))
}
error <- function(expr) stopifnot(inherits(tryCatch(expr, error = identity), "error"))
links <- c("logit", "probit", "cloglog", "loglog", "cauchit")

# All derivative components against independent central differences.
y <- rep(1:5, 3)
mu <- seq(-2, 2, length.out = length(y))
th <- log(c(0.7, 1.2, 0.9))
w <- seq(0.5, 2, length.out = length(y)); w[3] <- 0
h <- 1e-5
for (link in links) {
  f <- ocat_link(R = 5, link = link)
  dd <- function(mu = mu, theta = th) f$Dd(y, mu, theta, w, 2)
  a <- dd(mu, th)
  plus <- dd(mu+h, th); minus <- dd(mu-h, th)
  for (pair in list(c("D", "Dmu"), c("Dmu", "Dmu2"), c("Dmu2", "Dmu3"),
                    c("Dmu3", "Dmu4"), c("Dth", "Dmuth"),
                    c("Dmuth", "Dmu2th"), c("Dmu2th", "Dmu3th"),
                    c("Dth2", "Dmuth2"), c("Dmuth2", "Dmu2th2"))) {
    near((plus[[pair[1]]]-minus[[pair[1]]])/(2*h), a[[pair[2]]], 2e-5)
  }
  if (link == "cauchit") near((plus$EDmu2-minus$EDmu2)/(2*h), a$EDmu3, 2e-5)
  near(a$D, f$dev.resids(y, mu, w, th))
  for (j in seq_along(th)) {
    tp <- tm <- th; tp[j] <- tp[j]+h; tm[j] <- tm[j]-h
    dp <- dd(mu, tp); dm <- dd(mu, tm)
    near((dp$EDmu2-dm$EDmu2)/(2*h), a$EDmu2th[, j], 2e-5)
    for (pair in list(c("D", "Dth"), c("Dmu", "Dmuth"),
                      c("Dmu2", "Dmu2th"), c("Dmu3", "Dmu3th"))) {
      near((dp[[pair[1]]]-dm[[pair[1]]])/(2*h), a[[pair[2]]][, j], 2e-5)
    }
    col <- 0L
    for (k in seq_along(th)) for (l in k:length(th)) {
      col <- col+1L
      if (l == j) for (pair in list(c("Dth", "Dth2"), c("Dmuth", "Dmuth2"),
                                     c("Dmu2th", "Dmu2th2"))) {
        near((dp[[pair[1]]][, k]-dm[[pair[1]]][, k])/(2*h), a[[pair[2]]][, col], 2e-5)
      }
    }
  }
  # Distribution identities, integer quantiles, endpoints, log tails, recycling.
  f$putTheta(th)
  p <- seq(0.01, 0.99, length.out = length(mu))
  q <- f$qf(p, mu)
  stopifnot(all(q >= 1 & q <= 5), all(f$cdf(q, mu) >= p),
            all(f$cdf(q-1, mu) < p), all(f$qf(0, mu) == 1), all(f$qf(1, mu) == 5))
  near(f$qf(log(p), mu, log.p = TRUE), q)
  near(f$qf(1-p, mu, lower.tail = FALSE), q)
  near(f$qf(log1p(-p), mu, lower.tail = FALSE, log.p = TRUE), q)
  near(f$cdf(2.9, mu), f$cdf(2, mu))
  near(f$cdf(q, mu) + f$cdf(q, mu, lower.tail = FALSE), rep(1, length(mu)))
  near(exp(f$cdf(q, mu, log.p = TRUE)), f$cdf(q, mu))
  stopifnot(all(f$cdf(-Inf, mu) == 0), all(f$cdf(Inf, mu) == 1),
            is.na(f$cdf(NA_real_, 0)), is.na(f$qf(NA_real_, 0)),
            is.na(f$qf(0.5, NA_real_)), length(f$qf(numeric(), 0)) == 0)
  # Generalized inverse also holds exactly at a cumulative boundary.
  near(f$qf(f$cdf(2, mu), mu), rep(2, length(mu)))
  error(f$cdf(1:3, 1:2))
  error(f$qf(0.5, Inf))
  stopifnot(is.nan(suppressWarnings(f$qf(-0.1, 0))))
  set.seed(177)
  sim <- f$rd(rep(0, 25000))
  freq <- tabulate(sim, 5)/length(sim)
  near(freq, diff(f$cdf(0:5, 0)), 0.02)
  stopifnot(identical(mgcv::fix.family.rd(f)$rd, f$rd),
            identical(mgcv::fix.family.qf(f)$qf, f$qf))
}

# Logistic derivatives independently match mgcv's established family.
f <- ocat_link(R = 5)
a <- f$Dd(y, mu, th, w, 2)
b <- mgcv::ocat(R = 5)$Dd(y, mu, th, w, 2)
for (nm in intersect(names(a), names(b))) near(a[[nm]], b[[nm]], 2e-7)

set.seed(91)
n <- 500L
x <- runif(n, -2, 2)
d <- data.frame(x = x, y = as.integer(cut(0.6*x + rlogis(n),
                          c(-Inf, -1, 0.3, 1.4, Inf))))
d$ordered <- ordered(d$y)
for (link in links) {
  # Pure ML comparison: no smoothing penalties or REML integration differences.
  if (link != "cauchit") {
  f <- gam(y ~ x, data = d, family = ocat_link(R = 4, link = link), method = "ML")
  p <- MASS::polr(ordered ~ x, data = d,
                  method = if (link == "logit") "logistic" else link,
                  control = list(reltol = 1e-10))
  near(predict(f, type = "response"), predict(p, type = "probs"), 3e-5)
  near(f$family$getTheta(TRUE) - coef(f)[1], p$zeta, 3e-5)
  } else {
    # polr clips CDF arguments to +/-100, materially changing Cauchy tails.
    # Also avoid mgcv 1.9-4's native gam ML failure with negative curvature
    # and no penalized terms. Use bam and an independent exact likelihood.
    f <- bam(y ~ x, data=d, family=ocat_link(R=4,link=link), method="ML",
             control=gam.control(epsilon=1e-12,maxit=500))
    objective <- function(v) {
      cut <- c(-Inf, -1, -1+exp(v[3]), -1+sum(exp(v[3:4])), Inf)
      eta <- v[1]+v[2]*d$x
      -sum(log(pcauchy(cut[d$y+1L]-eta)-pcauchy(cut[d$y]-eta)))
    }
    ref <- optim(c(0,0.5,0,0),objective,method="BFGS",control=list(reltol=1e-12))
    near(c(coef(f),f$family$getTheta()),ref$par,2e-5)
  }
  # Every supported fitting route, including discrete prediction and SEs.
  for (engine in c("gam", "bam")) {
    methods <- if (engine == "gam") c("REML", "ML", "NCV") else c("REML", "ML", "fREML", "NCV")
    for (method in methods) {
      if (engine == "bam" && method == "NCV" && packageVersion("mgcv") < "1.9.4") next
      args <- list(formula = y ~ s(x, k = 6, bs = "cr"), data = d,
                   family = ocat_link(R = 4, link = link), method = method)
      if (engine == "bam") args$discrete <- method %in% c("fREML", "NCV")
      g <- do.call(get(engine), args)
      pr <- predict(g, newdata = d[1:12, ], type = "response", se.fit = TRUE)
      stopifnot(all(is.finite(pr$fit)), all(pr$fit >= 0), all(is.finite(pr$se.fit)),
                all(pr$se.fit >= 0), identical(g$method, method))
      near(rowSums(pr$fit), rep(1, 12), 1e-12)
      cp <- g$family$cdf(2, predict(g, newdata = d[1:12, ], type = "link"))
      near(cp, rowSums(pr$fit[, 1:2]))
    }
  }
  # Fixed thresholds and serialization retain the same distribution state.
  fixed <- gam(y ~ s(x, k = 6), data = d, method = "REML",
               family = ocat_link(theta = c(1.3, 1.1), link = link))
  near(fixed$family$getTheta(TRUE), c(-1, 0.3, 1.4))
  z <- unserialize(serialize(fixed, NULL))
  near(predict(z, type = "response"), predict(fixed, type = "response"))
  stopifnot(all(is.finite(residuals(fixed))), all(is.finite(residuals(fixed, type = "response"))))
}
# Additional bam fREML route without discretization.
g <- bam(y ~ s(x, k=6), data=d, family=ocat_link(R=4, link="probit"), method="fREML")
stopifnot(all(is.finite(coef(g))))
# Distribution callbacks accept fitted and predicted latent locations.
mu <- fitted(g)
stopifnot(length(g$family$rd(mu)) == n)
near(g$family$qf(.5, mu), g$family$qf(.5, predict(g, d, type = "link")))
stopifnot(length(g$family$rd(predict(g, d[1:7, ], type = "link"))) == 7L)

# Binary outcome uses the fixed first threshold and still needs higher derivatives.
d$binary <- as.integer(d$y > 2)+1L
for (link in c("logit", "probit", "cloglog", "cauchit")) {
  engine <- if (link == "cauchit") bam else gam
  f <- engine(binary ~ x, data=d, family=ocat_link(R=2,link=link), method="ML",
              control=gam.control(epsilon=1e-10,maxit=200))
  b <- glm(I(binary == 1) ~ x, data=d, family=binomial(link))
  near(predict(f,type="response")[,1], fitted(b), 2e-5)
}
# Fractional weights and zero-weight observations: compare with frequency expansion.
d$w <- rep(c(0,1,2,3),length.out=n)
f <- gam(y~x,data=d,weights=w,family=ocat_link(R=4,link="probit"),method="ML")
e <- d[rep(seq_len(n),d$w),]
b <- gam(y~x,data=e,family=ocat_link(R=4,link="probit"),method="ML")
near(predict(f,type="response"),predict(b,newdata=d,type="response"),2e-6)
for (bad in list(NULL, 1, 2.5, NA_real_, Inf, c(3,4))) error(ocat_link(R=bad))
error(ocat_link(theta=c(1,0)))
error(ocat_link(R=5,theta=c(1,2)))
error(ocat_link(R=4,link="identity"))
start_family <- ocat_link(theta=c(-0.4, 2), link="probit")
near(start_family$preinitialize(d$y, start_family)$Theta, log(c(.4,2)))
stopifnot(start_family$n.theta == 2L)
error(gam(y~x,data=transform(d,y=y+.1),family=ocat_link(R=4)))
error(gam(y~x,data=transform(d,y=ordered(y)),family=ocat_link(R=4)))

# Extreme-value tails: exact endpoint derivatives avoid catastrophic cancellation.
for (link in c("cloglog", "loglog")) {
  f <- ocat_link(R=4,link=link)
  m <- if(link == "cloglog") -100 else 100
  y <- if(link == "cloglog") 4L else 1L
  a <- f$Dd(y, m, c(0,0), 1, 2)
  z <- if(link == "cloglog") 1-m else 1+m
  for(k in 1:4) {
    name <- if(k == 1) "Dmu" else paste0("Dmu",k)
    sign <- if(link == "cloglog") (-1)^k else 1
    near(a[[name]]/(2*exp(z)),sign,1e-12)
  }
}
# Non-default leave-neighbourhood-out structures and one free threshold.
set.seed(70)
d3 <- data.frame(x=runif(240))
d3$y <- as.integer(cut(sin(5*d3$x)+rnorm(240),c(-Inf,-1,1,Inf)))
groups <- split(seq_len(nrow(d3)), rep(seq_len(80),each=3))
nei <- list(a=unlist(groups,use.names=FALSE),ma=seq(3,240,3),d=seq_len(240),md=seq(3,240,3))
for(engine in c("gam","bam")) {
  if(engine == "bam" && packageVersion("mgcv") < "1.9.4") next
  args <- list(formula=y~s(x,k=6),data=d3,family=ocat_link(R=3,link="probit"),method="NCV",nei=nei)
  if(engine == "bam") args$discrete <- TRUE
  g3 <- do.call(get(engine),args)
  stopifnot(all(is.finite(coef(g3))))
}
# Offset handling and response standard errors against manual delta method.
d3$o <- seq(-.2,.2,length.out=nrow(d3))
f <- gam(y~x+offset(o),data=d3,family=ocat_link(R=3,link="probit"),method="ML")
p <- MASS::polr(ordered(y)~x+offset(o),data=d3,method="probit",control=list(reltol=1e-10))
near(predict(f,type="response"),predict(p,type="probs"),3e-5)
pr <- predict(f,type="response",se.fit=TRUE)
lp <- predict(f,se.fit=TRUE)
a <- f$family$getTheta(TRUE)
near(pr$se.fit[,1],dnorm(a[1]-lp$fit)*lp$se.fit)
