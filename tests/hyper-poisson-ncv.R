library("inlaws")
near <- function(a,b,tol=1e-5) stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol)))
eta_fun<-getFromNamespace(".hp_eta","inlaws")
set.seed(407);n<-200L
d<-data.frame(x=runif(n),z=runif(n),a=runif(n,-.1,.1),b=runif(n,-.1,.1))
f<-hyper_poissonls()
d$y<-f$rd(cbind(exp(1+.5*sin(6*d$x)+d$a),exp(.8+.3*d$z+d$b)))
nei<-list(a=c(4,5,6,10,11,12,4,5,6,19,20,21),ma=c(3,6,9,12),
          d=c(5,11,12,5,20),md=c(1,3,4,5),jackknife=FALSE)
ns<-asNamespace("mgcv")
# Independent probability sums, separate from production likelihood numerics.
logmass<-function(y,eta) {
  k<-0:1500;s<-exp(eta[2]);mu<-exp(eta[1])
  prob<-function(a) {v<-k*a+lgamma(s)-lgamma(s+k);w<-exp(v-max(v));w/sum(w)}
  a<-uniroot(function(a) sum(k*prob(a))-mu,c(-40,log((mu+1)*(s+1))+4),tol=1e-13)$root
  log(prob(a)[y+1L])
}
for(shared in c(FALSE,TRUE)) {
  form<-if(shared) list(y~1+offset(a),~z+offset(b),1+2~s(x,k=5)-1) else
    list(y~s(x,k=5)+offset(a),~s(z,k=4)+offset(b))
  G<-gam(form,data=d,family=f,method="REML",fit=FALSE)
  Sl<-get("Sl.setup",ns)(G)
  X<-get("Sl.initial.repara",ns)(Sl,G$X,both.sides=FALSE)
  rho<-rep(log(2),length(G$sp));start<-NULL
  call<-function(rho,method="REML",deriv=1,q=FALSE,neighbours=nei) {
    ff<-f;ff$qapprox<-q
    get("gam.fit5",ns)(X,G$y,rho,Sl,weights=G$w,offset=G$offset,family=ff,
      deriv=deriv,scoreType=method,Mp=-1,nei=neighbours,gamma=1.2,start=start,
      control=gam.control(epsilon=1e-12))
  }
  base<-call(rho,deriv=2);start<-base$coefficients
  for(j in seq_along(rho)) {
    step<-numeric(length(rho));step[j]<-.01
    mm<-call(rho-2*step);m<-call(rho-step);p<-call(rho+step);pp<-call(rho+2*step)
    near(base$REML1[j],(mm$REML-8*m$REML+8*p$REML-pp$REML)/.12,1e-4)
    near(base$REML2[,j],(mm$REML1-8*m$REML1+8*p$REML1-pp$REML1)/.12,3e-4)
  }
  # Reference deletion updates build per-observation predictor design maps.
  # This accumulates overlapping coefficient blocks independently of the family.
  reference<-function(rho,q,neighbours=nei) {
    b<-call(rho);be<-b$coefficients;jj<-attr(X,"lpi");p<-ncol(X)
    eta<-eta_fun(X,be,G$offset)
    ll<-f$ll(G$y,X,be,G$w,f,G$offset,deriv=1,ncv=TRUE)
    H<-b$St-b$lbb;nd<-length(neighbours$d)
    cv<-matrix(0,nd,2);changes<-matrix(0,length(neighbours$ma),p)
    pa<-pd<-0L
    for(i in seq_along(neighbours$ma)) {
      a<-neighbours$a[seq.int(pa+1L,neighbours$ma[i])]
      id<-seq.int(pd+1L,neighbours$md[i]);pred<-neighbours$d[id]
      hd<-H;score<-numeric(p)
      for(row in a) {
        Z<-matrix(0,2,p)
        for(k in 1:2) Z[k,jj[[k]]]<-X[row,jj[[k]]]
        score<-score+drop(crossprod(Z,ll$l1[row,]))
        hd<-hd+crossprod(Z,matrix(ll$l2[row,f$tri$i2],2,2)%*%Z)
      }
      changes[i,]<-solve(hd,score)
      for(k in 1:2) cv[id,k]<-eta[pred,k]-drop(X[pred,jj[[k]],drop=FALSE]%*%changes[i,jj[[k]]])
      pa<-neighbours$ma[i];pd<-neighbours$md[i]
    }
    ix<-neighbours$d
    base<-G$w[ix]*vapply(seq_len(nd),function(i) logmass(G$y[ix[i]],eta[ix[i],]),numeric(1))
    if(!q) rows <- -1.2*G$w[ix]*vapply(seq_len(nd),function(i) logmass(G$y[ix[i]],cv[i,]),numeric(1))+.2*base else {
      delta<-cv-eta[ix,,drop=FALSE];rows<--base
      for(i in seq_len(nd)) rows[i]<-rows[i]-1.2*(sum(ll$l1[ix[i],]*delta[i,])+
        .5*drop(crossprod(delta[i,],matrix(ll$l2[ix[i],f$tri$i2],2,2)%*%delta[i,])))
    }
    list(score=sum(rows),rows=rows,eta=cv,changes=changes)
  }
  for(q in c(FALSE,TRUE)) {
    fit<-call(rho,"NCV",q=q);ref<-reference(rho,q)
    near(fit$NCV,ref$score,1e-8);near(attr(fit$NCV,"eta.cv"),ref$eta,1e-8)
    row_grad<-matrix(0,length(nei$d),length(rho))
    for(j in seq_along(rho)) {
      rp<-rm<-rho;rp[j]<-rp[j]+1e-4;rm[j]<-rm[j]-1e-4
      plus<-reference(rp,q);minus<-reference(rm,q)
      near(fit$NCV1[j],(plus$score-minus$score)/2e-4,2e-5)
      row_grad[,j]<-(plus$rows-minus$rows)/2e-4
    }
    near(fit$Vg,crossprod(row_grad),2e-5)
  }
  # Jackknife coefficient perturbations agree with independently signed solves.
  ll<-f$ll(G$y,X,start,G$w,f,G$offset,deriv=1,ncv=TRUE);ll$gamma<-1
  jk<-f$ncv(X,G$y,G$w,nei,start,f,ll,R=chol(base$St-base$lbb),offset=G$offset,deriv=-1)
  ref<-reference(rho,FALSE)
  near(attr(jk$NCV,"eta.cv"),ref$eta,1e-8)
  near(attr(jk$NCV,"deta.cv"),ref$changes,1e-8)
  pa<-pd<-0L
  for(i in seq_along(nei$ma)) {
    a<-nei$a[seq.int(pa+1L,nei$ma[i])];id<-seq.int(pd+1L,nei$md[i]);w<-G$w;w[a]<-0
    exact<-get("gam.fit5",ns)(X,G$y,rho,Sl,weights=w,offset=G$offset,family=f,
      deriv=0,Mp=-1,start=start,control=gam.control(epsilon=1e-12))
    stopifnot(max(abs(ref$eta[id,,drop=FALSE]-exact$linear.predictors[nei$d[id],,drop=FALSE]))<.08)
    pa<-nei$ma[i];pd<-nei$md[i]
  }
  # Fractional and zero weights affect every deletion contribution.
  G$w<-rep(c(0,.5,2,1),length.out=n)
  base<-call(rho,deriv=2);start<-base$coefficients
  for(q in c(FALSE,TRUE)) {
    fit<-call(rho,"NCV",q=q);ref<-reference(rho,q)
    near(fit$NCV,ref$score,1e-8)
    for(j in seq_along(rho)) {
      rp<-rm<-rho;rp[j]<-rp[j]+1e-4;rm[j]<-rm[j]-1e-4
      near(fit$NCV1[j],(reference(rp,q)$score-reference(rm,q)$score)/2e-4,3e-5)
    }
  }
}

# Complete-case and corrupted missing-row neighbourhoods use the shared guard.
form<-list(y~s(x,k=5)+offset(a),~z+offset(b))
dd<-d;dd$x[2]<-NA
bad<-tryCatch(gam(form,data=dd,family=f,method="NCV",nei=nei,fit=FALSE),error=identity)
if(packageVersion("mgcv")=="1.9.4") stopifnot(inherits(bad,"error"),grepl("neighbourhood preprocessing",conditionMessage(bad)))
whole<-list(a=1:n,ma=1:n,d=1:n,md=1:n)
safe<-gam(form,data=dd,family=f,method="NCV",nei=whole,fit=FALSE)
stopifnot(inherits(safe,"gam.prefit"))
