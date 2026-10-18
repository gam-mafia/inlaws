# Integrate analytic log-response density jets over observed log bounds.
# Quadrature locations and scaling are fixed at the evaluation point: all
# parameter derivatives are taken under the integral, including at Q = 0.
.cgen_gamma_integral <- function(lo, hi, eta, order, span = NULL) {
  core <- .gen_gamma_core(eta, order); j <- core$j
  n <- nrow(eta); finite <- is.finite(lo) & is.finite(hi)
  boundary <- ifelse(is.finite(lo), lo, hi)
  direction <- ifelse(is.finite(lo), 1, -1)
  sigma <- core$sigma[, 1]; Q <- core$Q[, 1]; a <- core$location[, 1]
  z <- (boundary-a)/sigma
  slope <- -z
  nonzero <- Q != 0
  slope[nonzero] <- -expm1(Q[nonzero]*z[nonzero])/Q[nonzero]
  # In a super-exponential tail the integration width can be smaller than
  # one ulp of the log boundary. Laplace's boundary expansion then has
  # relative error O(Q^2 exp(-Q*z)); retain its full analytic Taylor jet.
  x <- Q*z
  asym <- !finite & direction*Q > 0 & x > 30 & Q^2*exp(-x) < 1e-12
  asym[is.na(asym)] <- FALSE
  if (any(asym)) {
    out <- j$constant(numeric(n)); ix <- which(asym); r <- core
    for (name in c("location", "sigma", "Q", "normalizer", "correction", "t"))
      r[[name]] <- core[[name]][ix, , drop = FALSE]
    zz <- j$multiply(j$constant(boundary[ix])-r$location, j$exp(-r$t))
    slope_jet <- j$multiply(j$expm1(j$multiply(r$Q, zz)),
      j$inverse(j$multiply(r$sigma, r$Q)))
    negative <- Q[ix] < 0
    slope_jet[negative, ] <- -slope_jet[negative, , drop = FALSE]
    out[ix, ] <- .gen_gamma_log_density(boundary[ix], r)-j$log(slope_jet)
    if (any(!asym)) out[!asym, ] <- .cgen_gamma_integral(lo[!asym], hi[!asym],
      eta[!asym, , drop = FALSE], order, if (is.null(span)) NULL else span[!asym])
    return(out)
  }
  width <- sigma/pmax(abs(slope), 1/(1+abs(Q)))
  width[finite] <- if (is.null(span)) hi[finite]-lo[finite] else span[finite]
  peak_at <- pmin(hi, pmax(lo, a))
  reference <- .gen_gamma_log_density(peak_at, core); peak <- reference[, 1]
  previous <- NULL
  for (size in c(16L, 32L, 64L, 128L, 256L)) {
    rule <- .gauss_legendre(size); u <- (rule$nodes+1)/2
    ix <- rep(seq_len(n), each = size); uu <- rep(u, n)
    at <- boundary[ix] + direction[ix]*width[ix]*uu/(1-uu)
    lj <- log(width[ix])-2*log1p(-uu)
    fi <- finite[ix]
    at[fi] <- lo[ix[fi]] + width[ix[fi]]*uu[fi]
    lj[fi] <- log(width[ix[fi]])
    r <- core
    for (name in c("location", "sigma", "Q", "normalizer", "correction", "t"))
      r[[name]] <- core[[name]][ix, , drop = FALSE]
    # Far-tail nodes can overflow exp(Q*z). Their contribution and all its
    # first four derivatives are below floating-point accuracy.
    zz <- (at-a[ix])/sigma[ix]; xx <- Q[ix]*zz
    penalty <- zz^2/2
    nz <- abs(xx) > 1e-5
    penalty[nz] <- (expm1(xx[nz])-xx[nz])/Q[ix[nz]]^2
    log_base <- core$normalizer[ix, 1]-eta[ix, 2]-penalty
    active <- is.finite(log_base) & log_base-peak[ix]+lj > -650
    at[!active] <- a[ix[!active]]
    log_f <- .gen_gamma_log_density(at, r)-reference[ix, , drop = FALSE]
    log_f[, 1] <- log_f[, 1]+lj
    log_f[!active, ] <- 0
    integrand <- j$exp(log_f)
    integrand[!active, ] <- 0
    integrand <- integrand * rep(rule$weights/2, n)
    total <- rowsum(integrand, ix, reorder = FALSE)
    answer <- j$log(total)+reference
    if (!is.null(previous) && all(is.finite(answer)) &&
        max(abs(answer-previous)/(1+abs(answer))) < 2e-10) return(answer)
    before <- previous; previous <- answer
  }
  if (any(!is.finite(answer))) return(answer * NaN)
  if (max(abs(answer-before)/(1+abs(answer))) > 1e-7)
    return(answer * NaN)
  answer
}

# Both log tails, evaluating the smaller tail directly and complementing its
# jet. This avoids subtracting probabilities close to one.
.cgen_gamma_tails <- function(bound, eta, order) {
  j <- .gen_gamma_jets(order); n <- length(bound)
  lp <- .gen_gamma_probability(bound, exp(eta), log.p = TRUE)
  up <- .gen_gamma_probability(bound, exp(eta), lower.tail = FALSE, log.p = TRUE)
  if (anyNA(lp) || anyNA(up)) return(list(
    lower = j$constant(rep(NaN, n)), upper = j$constant(rep(NaN, n))))
  if (!order) return(list(lower = matrix(lp, ncol = 1L),
    upper = matrix(up, ncol = 1L)))
  lower <- lp <= -log(2)
  small <- .cgen_gamma_integral(ifelse(lower, -Inf, log(bound)),
    ifelse(lower, log(bound), Inf), eta, order)
  complement <- j$constant(rep(1, n))-j$exp(small)
  complement[, 1] <- -expm1(small[, 1])
  complement <- j$log(complement)
  # If even the log of the small tail is beyond floating-point range,
  # its complement and every representable derivative equal zero.
  complement[ifelse(lower, lp, up) == -Inf, ] <- 0
  lo <- hi <- small
  lo[!lower, ] <- complement[!lower, , drop = FALSE]
  hi[lower, ] <- complement[lower, , drop = FALSE]
  list(lower = lo, upper = hi)
}

.cgen_gamma_eval <- function(y, eta, order = 0L) {
  response <- .cg_response(y); y <- response$y; j <- .gen_gamma_jets(order)
  ell <- j$constant(numeric(length(y)))
  # Reject numerical trial points before the public distribution callback
  # checks response-scale parameters (e.g. a null-model intercept step).
  if (any(!is.finite(eta)) || any(!is.finite(exp(eta)) | exp(eta) <= 0)) {
    ell[,] <- NaN; return(list(ell = ell, powers = j$powers))
  }
  exact <- response$exact
  if (any(exact)) ell[exact, ] <- .gen_gamma_eval(y[exact], eta[exact, , drop = FALSE], order)$ell
  # Full-support intervals have probability one and zero derivatives.
  left <- !exact & response$lo == 0 & is.finite(response$hi)
  right <- !exact & response$lo > 0 & is.infinite(response$hi)
  if (any(left)) ell[left, ] <- .cgen_gamma_tails(response$hi[left], eta[left, , drop = FALSE], order)$lower
  if (any(right)) ell[right, ] <- .cgen_gamma_tails(response$lo[right], eta[right, , drop = FALSE], order)$upper
  interval <- !exact & response$lo > 0 & is.finite(response$hi)
  # Short intervals use direct integration, retaining log probabilities even
  # when their CDF endpoints round to the same number.
  width <- log1p((response$hi-response$lo)/response$lo)
  short <- interval & width < exp(eta[, 2])
  if (any(short)) {
    lo <- log(response$lo[short])
    ell[short, ] <- .cgen_gamma_integral(lo, lo+width[short], eta[short, , drop = FALSE], order, span = width[short])
  }
  wide <- interval & !short
  if (any(wide)) {
    ix <- which(wide); ee <- eta[ix, , drop = FALSE]
    l <- .cgen_gamma_tails(response$lo[ix], ee, order)
    h <- .cgen_gamma_tails(response$hi[ix], ee, order)
    lower <- h$lower[, 1] < l$upper[, 1]; lower[is.na(lower)] <- TRUE
    big <- h$lower; small <- l$lower
    big[!lower, ] <- l$upper[!lower, , drop = FALSE]
    small[!lower, ] <- h$upper[!lower, , drop = FALSE]
    delta <- small-big
    diff <- j$constant(rep(1, length(ix)))-j$exp(delta)
    diff[, 1] <- -expm1(delta[, 1])
    ell[ix, ] <- big+j$log(diff)
  }
  list(ell = ell, powers = j$powers)
}

# Location maximizing a finite-interval probability solves f(log lo)=f(log hi).
# With d=Q*log(hi/lo)/sigma, a=log(lo)+(sigma/Q)*log(expm1(d)/d).
.cgen_gamma_saturated <- function(y, mu) {
  r <- .cg_response(y); y <- r$y; a <- .gen_gamma_distribution(mu)
  location <- log(r$y)
  interval <- !r$exact & r$lo > 0 & is.finite(r$hi)
  if (any(interval)) {
    ix <- which(interval); width <- log1p((r$hi[ix]-r$lo[ix])/r$lo[ix])
    d <- a$Q[ix]*width/a$sigma[ix]
    # log(expm1(d)/d)/d = 1/2 + d/24 - d^3/2880 + ...
    ratio <- .5+d/24-d^3/2880+d^5/181440
    reg <- abs(d) > .01
    x <- d[reg]
    log_ratio <- ifelse(x > 50, x+log1p(-exp(-x))-log(abs(x)),
      log(abs(expm1(x)))-log(abs(x)))
    ratio[reg] <- log_ratio/x
    location[ix] <- log(r$lo[ix])+width*ratio
  }
  eta <- log(mu); eta[, 1] <- location+a$correction
  out <- numeric(length(y)); exact <- r$exact
  if (any(exact)) out[exact] <- .gen_gamma_eval(y[exact], eta[exact, , drop = FALSE])$ell[, 1]
  if (any(interval)) {
    yy <- r$y[interval]; attr(yy, "censor") <- r$censor[interval]
    out[interval] <- .cgen_gamma_eval(yy, eta[interval, , drop = FALSE])$ell[, 1]
  }
  sign <- sign(location-a$location)
  sign[!exact & r$lo == 0] <- -1
  sign[!exact & is.infinite(r$hi)] <- 1
  sign[r$lo == 0 & is.infinite(r$hi)] <- 0
  list(loglik = out, sign = sign)
}
