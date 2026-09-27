#' Dirichlet-multinomial family for count compositions
#'
#' A multi-predictor family for [mgcv::gam()]. Each response row is a vector
#' of counts, conditional on its total. Category probabilities and concentration
#' can depend on different smooth predictors.
#'
#' @param K Positive integer: there are `K + 1` categories, `K` baseline-logit
#'   predictors, and one log-concentration predictor.
#' @return An object inheriting from `general.family`, `extended.family`, and
#'   `family`, suitable for `mgcv::gam()`.
#' @details
#' Supply a list of `K + 1` formulas. The first includes the matrix response
#' (for example `cbind(a, b, c)`). The first `K` predictors model
#' `log(p_j / p_1)` for categories 2 through `K + 1`. The last models
#' `log(concentration)`; use `~1` for constant concentration. Shared terms
#' can be specified as in [mgcv::formula.gam()]. Internally all links are
#' identity links because the likelihood is evaluated on the predictor scale.
#'
#' If `p` denotes the probabilities and `phi` the concentration, the Dirichlet
#' parameters are `phi * p`. The conditional mean is `N * p` and the marginal
#' variance in category `j` is
#' `N * p_j * (1 - p_j) * (1 + (N - 1) / (1 + phi))`.
#' Increasing concentration reduces overdispersion towards the multinomial
#' limit. Totals `N` are computed from the response, not from prior weights.
#' Prior weights multiply each row's complete log likelihood.
#'
#' Counts must be finite non-negative integers with positive row totals.
#' Individual zero counts are allowed, but each category must occur in at least
#' one positive-weight row. Concentration cannot be estimated if every
#' positive-weight row has total one. Missing data follow `gam()`'s `na.action`.
#' Nearly multinomial data can give weakly identified, very large concentration;
#' a warning is issued if fitted concentration exceeds one million. No upper
#' bound on concentration is imposed.
#'
#' Fitting uses REML with the extended Fellner-Schall optimizer (`optimizer =
#' "efs"`), which mgcv selects automatically for this family. Analytic first
#' and second likelihood derivatives are supplied. Higher-order outer Newton
#' optimization and neighbourhood cross-validation are not supported. Because
#' of mgcv's EFS smoothing-parameter mapping, smoothing parameters must be all
#' free and unlinked, or all fixed. Partly fixed parameters and shared smoothing
#' parameter IDs are rejected; shared coefficients across predictors are supported.
#' With mgcv 1.9-3, supply fixed smoothing parameters to a second call
#' `gam(G = G, sp = ...)`, where `G` was built with `gam(..., fit = FALSE)`,
#' to avoid that version's list-formula setup bug.
#' `mgcv::bam()` currently rejects general families, including this one;
#' a scalable BAM fitting engine is a separate future extension.
#'
#' `predict(..., type = "link")` returns the `K` logits followed by log
#' concentration. `type = "response"` returns all `K + 1` probabilities, in
#' response-column order, followed by concentration. mgcv may omit column names
#' on this matrix. When requesting a custom prediction `block.size`, supply
#' `newdata` explicitly (a matrix-response limitation in mgcv).
#' Multiply the probability columns by desired totals to get
#' expected counts. `se.fit = TRUE` uses the full coefficient covariance and
#' the multivariate delta method. These are parameter standard errors, not
#' predictive standard deviations. As for other mgcv general families,
#' `fitted()` contains the linear predictors; use `predict()` for response values.
#'
#' Response and Pearson residuals are matrices, with one column per category.
#' Pearson residuals use marginal variances and the square root of prior weights;
#' their columns are correlated. Default deviance residuals are unsigned
#' `sqrt(-2 * weight * log_probability)` for each row, relative to a
#' probability-one benchmark, not a saturated Dirichlet-multinomial fit.
#' The reported null deviance uses a fitted intercept-only model with the same
#' offsets and weights. The full count likelihood, including multinomial
#' coefficients, is used for `logLik()` and AIC.
#'
#' Simulate counts by drawing gamma variates with shapes `phi * p`, normalizing
#' them to probabilities, then calling `rmultinom()` with the desired total.
#' The example illustrates this at moderate concentration.
#' @seealso [mgcv::multinom()], [mgcv::gam()]
#' @export
#' @examples
#' set.seed(12)
#' n <- 200
#' dat <- data.frame(x = runif(n), z = runif(n))
#' eta <- cbind(0, sin(2 * pi * dat$x), 0.5 * dat$x)
#' p <- exp(eta) / rowSums(exp(eta))
#' phi <- exp(2 + dat$z)
#' totals <- sample(15:30, n, replace = TRUE)
#' dat$counts <- t(vapply(seq_len(n), function(i) {
#'   q <- rgamma(3, shape = phi[i] * p[i, ])
#'   as.vector(rmultinom(1, totals[i], prob = q))
#' }, numeric(3)))
#' fit <- mgcv::gam(list(counts ~ s(x, k = 6), ~s(x, k = 6), ~z),
#'                  family = dirmult(K = 2), data = dat,
#'                  method = "REML", optimizer = "efs")
#' pr <- predict(fit, type = "response")
#' head(pr)
#' expected_counts <- totals * pr[, 1:3]
dirmult <- function(K=1) {
  if (length(K)!=1L || !is.numeric(K) || !is.finite(K) || K<1 || K!=round(K))
    stop("K must be a positive integer")
  C <- K+1L
  linfo <- lapply(seq_len(C),function(i) {
    z <- mgcv::fix.family.link(structure(list(link="identity",canonical="none",
                linkfun=identity,mu.eta=function(eta) rep(1,length(eta))),class="family"))
    c(stats::make.link("identity"),z[c("d2link","d3link","d4link")])
  })
  validate <- function(y,wt) {
    if (!is.matrix(y) || !is.numeric(y) || ncol(y)!=C)
      stop("dirmult requires a numeric count matrix with K+1 columns")
    if (any(!is.finite(y)) || any(y<0) || any(y!=floor(y)))
      stop("dirmult counts must be finite non-negative integers")
    if (any(!is.finite(rowSums(y))) || any(rowSums(y)<=0)) stop("dirmult row totals must be positive")
    if (length(wt)!=nrow(y) || any(!is.finite(wt)) || any(wt<0) || !any(wt>0))
      stop("dirmult weights must be finite, non-negative and not all zero")
    if (any(colSums(y[wt>0,,drop=FALSE])==0))
      stop("dirmult categories must occur in positive-weight observations")
    if (all(rowSums(y)[wt>0]==1))
      stop("concentration is not identifiable when all informative totals equal one")
  }
  preinitialize <- function(G) {
    validate(G$y,G$w)
    ## mgcv's EFS driver does not apply a reduced smoothing-parameter map.
    ## Fully fixed smoothing parameters use its separate no-sps path.
    if (length(G$sp)>0 && (!is.null(G$L) &&
        !isTRUE(all.equal(unname(G$L),diag(length(G$sp)),check.attributes=FALSE)) ||
        any(G$lsp0!=0)))
      stop("dirmult EFS requires all smoothing parameters free and unlinked, or all fixed")
    if (length(attr(G$X,"lpi"))!=C)
      stop("dirmult requires K probability formulas and one concentration formula")
    fam <- G$family
    fam$training.offset <- G$offset
    fam$category.names <- colnames(G$y)
    if (is.null(fam$category.names)) fam$category.names <- paste0("p",seq_len(C))
    list(family=fam)
  }
  ll <- function(y,X,coef,wt,family,offset=NULL,deriv=0,d1b=0,d2b=0,
                 Hp=NULL,rank=0,fh=NULL,D=NULL,eta=NULL,ncv=FALSE,sandwich=FALSE) {
    if (deriv>1 || ncv) stop("dirmult supports EFS fitting only")
    jj <- attr(X,"lpi")
    if (is.null(eta)) {
      eta <- matrix(0,nrow(y),family$nlp)
      for (j in seq_len(family$nlp)) {
        eta[,j] <- X[,jj[[j]],drop=FALSE]%*%coef[jj[[j]]]
        if (!is.null(offset) && length(offset)>=j && !is.null(offset[[j]]))
          eta[,j] <- eta[,j]+offset[[j]]
      }
    }
    if (any(!is.finite(eta))) return(list(l=-Inf,lb=rep(0,length(coef)),
                                         lbb=matrix(0,length(coef),length(coef))))
    z <- .dirmult_deriv(y,eta,wt,deriv>0)
    if (deriv) {
      ret <- mgcv::gamlss.gH(X,jj,z$l1,z$l2,family$tri$i2,deriv=0,sandwich=sandwich)
    } else ret <- list()
    ret$l <- z$l; ret$l0 <- z$l0
    ret
  }
  initialize <- expression({
    nobs <- nrow(y); n <- rep(1,nobs)
    ## mgcv 1.9-3's EFS path omits the offset argument. Restore it here,
    ## before gam.fit5 constructs either the likelihood or fitted values.
    if (is.null(offset)) offset <- family$training.offset
    if (is.null(weights)) weights <- rep(1,nobs)
    family$validate(y,weights)
    if (is.null(start)) {
      jj <- attr(x,"lpi"); nc <- ncol(y)
      yt <- log((y[,-1,drop=FALSE]+0.5)/(y[,1]+0.5))
      ## Starting concentration only: the search interval is not a fit bound.
      p0 <- colSums(y*weights)+0.5; p0 <- p0/sum(p0)
      e0 <- matrix(log(p0[-1]/p0[1]),nobs,nc-1,byrow=TRUE)
      ph0 <- stats::optimize(function(z) -family$initial.ll(y,cbind(e0,z),weights),
                      c(-5,10))$minimum
      yt <- cbind(yt,ph0)
      ## Joint least squares handles coefficients shared by predictors.
      H0 <- matrix(0,ncol(x),ncol(x)); b0 <- numeric(ncol(x))
      for (j in seq_len(nc)) {
        ix <- jj[[j]]; xx <- x[,ix,drop=FALSE]
        yy <- yt[,j]
        if (!is.null(offset) && length(offset)>=j && !is.null(offset[[j]])) yy <- yy-offset[[j]]
        H0[ix,ix] <- H0[ix,ix]+crossprod(xx,weights*xx)
        b0[ix] <- b0[ix]+drop(crossprod(xx,weights*yy))
      }
      if (!is.null(E)) H0 <- H0+crossprod(E)
      eg <- eigen(H0,symmetric=TRUE)
      keep <- eg$values>max(eg$values)*.Machine$double.eps^.75
      start <- drop(eg$vectors[,keep,drop=FALSE]%*%
                    (crossprod(eg$vectors[,keep,drop=FALSE],b0)/eg$values[keep]))
    }
  })
  predict <- function(family,se=FALSE,eta=NULL,y=NULL,X=NULL,beta=NULL,off=NULL,Vb=NULL) {
    jj <- attr(X,"lpi")
    if (is.null(eta)) {
      eta <- matrix(0,nrow(X),family$nlp)
      for (j in seq_len(family$nlp)) {
        eta[,j] <- X[,jj[[j]],drop=FALSE]%*%beta[jj[[j]]]
        if (length(off)>=j && !is.null(off[[j]])) eta[,j] <- eta[,j]+off[[j]]
      }
    }
    pr <- .dirmult_prob(eta)
    fit <- cbind(pr$p,concentration=exp(pr$logphi))
    cn <- family$category.names
    if (is.null(cn)) cn <- paste0("p",seq_len(family$nlp))
    colnames(fit) <- make.unique(c(cn,"concentration"))
    ret <- list(fit=fit)
    if (se) {
      se.fit <- fit*0; nc <- family$nlp
      for (k in seq_len(nc+1L)) {
        J <- matrix(0,nrow(X),ncol(X))
        for (j in seq_len(nc)) {
          d <- if (k==nc+1L) {
            if (j==nc) fit[,k] else rep(0,nrow(X))
          } else if (j==nc) rep(0,nrow(X)) else
            pr$p[,k]*((k==j+1L)-pr$p[,j+1L])
          J[,jj[[j]]] <- J[,jj[[j]],drop=FALSE]+X[,jj[[j]],drop=FALSE]*d
        }
        se.fit[,k] <- sqrt(pmax(0,rowSums((J%*%Vb)*J)))
      }
      ret$se.fit <- se.fit
    }
    ret
  }
  residuals <- function(object,type=c("deviance","pearson","response")) {
    type <- match.arg(type)
    eta <- object$linear.predictors; y <- object$y
    if (type=="deviance") return(sqrt(pmax(0,-2*object$prior.weights*
                                .dirmult_deriv(y,eta,deriv=FALSE)$l0)))
    pr <- .dirmult_prob(eta); total <- rowSums(y)
    r <- y-total*pr$p
    if (type=="pearson") {
      rho <- stats::plogis(-pr$logphi)
      v <- total*pr$p*(1-pr$p)*(1+(total-1)*rho)
      r <- sqrt(object$prior.weights)*r/sqrt(v)
    }
    r
  }
  postproc <- expression({
    object$deviance <- -2*object$l
    object$null.deviance <- object$family$null.deviance(object$y,object$prior.weights,
                          object$family$training.offset,object$linear.predictors)
    eta <- object$linear.predictors
    colnames(eta) <- c(paste0("logit",seq_len(ncol(eta)-1L)),"log.concentration")
    object$linear.predictors <- object$fitted.values <- eta
    object$family$training.offset <- NULL
    if (any(eta[,ncol(eta)]>log(1e6)))
      warning("large dirmult concentration: the multinomial boundary may be weakly identified")
  })
  i2 <- matrix(0L,C,C); h <- 0L
  for (j in seq_len(C)) for (k in j:C) {
    h <- h+1L; i2[j,k] <- i2[k,j] <- h
  }
  structure(list(family="dirmult",nlp=C,link=rep("identity",C),linfo=linfo,
      ll=ll,initialize=initialize,preinitialize=preinitialize,validate=validate,
      postproc=postproc,predict=predict,residuals=residuals,
      initial.ll=function(y,eta,wt) .dirmult_deriv(y,eta,wt,FALSE)$l,
      null.deviance=.dirmult_null,
      tri=list(i2=i2),available.derivs=0,
      d2link=1,d3link=1,d4link=1,ls=1,discrete.ok=FALSE),
      class=c("general.family","extended.family","family"))
}
