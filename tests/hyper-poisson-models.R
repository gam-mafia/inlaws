library("inlaws")
near <- function(a,b,tol=1e-5) stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol)))
error <- function(expr,pattern) {
  e<-tryCatch(expr,error=identity)
  stopifnot(inherits(e,"error"),grepl(pattern,conditionMessage(e)))
}
set.seed(407); n<-260L
d<-data.frame(x=runif(n),z=runif(n),a=runif(n,-.1,.1),b=runif(n,-.1,.1))
f<-hyper_poissonls()
d$y<-f$rd(cbind(exp(1+.5*sin(6*d$x)+d$a),exp(1+.5*d$z+d$b)))
form<-list(y~s(x,k=5)+offset(a),~s(z,k=4)+offset(b))
fits<-lapply(list(c("outer","newton"),c("outer","bfgs"),"efs"),function(opt)
  gam(form,data=d,family=f,method="REML",optimizer=opt,control=gam.control(efs.tol=1e-7)))
b<-fits[[1]]
for(g in fits) {
  stopifnot(g$outer.info$conv=="full convergence",all(is.finite(g$Vp)),all(g$sp>0))
  near(g$gcv.ubre,b$gcv.ubre,1e-5);near(g$fitted.values,b$fitted.values,.02)
}
for(method in c("NCV","QNCV")) {
  g<-gam(form,data=d,family=f,method=method)
  stopifnot(identical(g$method,method),g$outer.info$conv=="full convergence",
            all(is.finite(g$Vp)),is.finite(g$gcv.ubre))
}
nei<-list(a=c(4,5,6,10,11,12,4,5,6),ma=c(3,6,9),d=c(5,11,12,5),md=c(1,3,4),jackknife=FALSE)
for(method in c("NCV","QNCV")) {
  G<-gam(form,data=d,family=f,method="REML",fit=FALSE)
  g<-gam(G=G,sp=b$sp,method=method,nei=nei,gamma=1.2)
  stopifnot(identical(g$method,method),is.finite(g$gcv.ubre),all(is.finite(coef(g))))
}
near(b$deviance,sum(residuals(b,"deviance")^2))
near(residuals(b,"response"),d$y-fitted(b)[,1])
diag<-f$gratia_diagnostics(fitted(b),b$prior.weights,1,d$y,b$linear.predictors)
near(residuals(b,"pearson"),(d$y-fitted(b)[,1])/sqrt(diag$variance))
near(AIC(b),-2*as.numeric(logLik(b))+2*attr(logLik(b),"df"))
lp<-predict(b,newdata=d[1:7,],type="link",se.fit=TRUE)
pr<-predict(b,newdata=d[1:7,],type="response",se.fit=TRUE)
near(pr$fit,exp(lp$fit));near(pr$se.fit,pr$fit*lp$se.fit)
stopifnot(identical(colnames(fitted(b)),c("mu","sigma")),
          identical(dim(predict(b,newdata=d[1,,drop=FALSE],type="response")),c(1L,2L)))
clone<-unserialize(serialize(b,NULL))
near(predict(clone,type="response"),fitted(b))

# Constant dispersion, offsets, shared coefficients, and their prediction SEs.
constant<-gam(list(y~s(x,k=5)+offset(a),~1+offset(b)),data=d,family=f,method="REML",optimizer=c("outer","bfgs"))
stopifnot(all(is.finite(coef(constant))),constant$outer.info$conv=="full convergence")
shared<-list(y~1+offset(a),~z+offset(b),1+2~s(x,k=5)-1)
g<-gam(shared,data=d,family=f,method="REML")
stopifnot(g$outer.info$conv=="full convergence",all(is.finite(g$Vp)))
lp<-predict(g,newdata=d[1:7,],type="link",se.fit=TRUE)
pr<-predict(g,newdata=d[1:7,],type="response",se.fit=TRUE)
near(pr$fit,exp(lp$fit));near(pr$se.fit,pr$fit*lp$se.fit)

# Integer weights equal row replication in the same parametric basis.
d$w<-rep(c(0,1,2),length.out=n)
param<-list(y~x+offset(a),~z+offset(b))
fw<-gam(param,data=d,family=f,weights=w,method="REML")
fr<-gam(param,data=d[rep(seq_len(n),d$w),],family=f,method="REML")
near(coef(fw),coef(fr),2e-5)
dd<-d;dd$x[c(2,5)]<-NA
sub<-gam(param,data=dd,subset=seq_len(n)>10,family=f,method="REML")
explicit<-gam(param,data=dd[which(seq_len(n)>10 & complete.cases(dd)),],family=f,method="REML")
near(coef(sub),coef(explicit))

# All-fixed EFS is supported, reduced and partly fixed maps are rejected.
G<-gam(form,data=d,family=f,method="REML",fit=FALSE)
fixed<-gam(G=G,sp=b$sp,optimizer="efs",method="REML")
near(fitted(fixed),fitted(b),1e-4)
mixed<-G; mixed$L<-matrix(c(1,0),2,1);mixed$sp<-G$sp[1];mixed$lsp0<-c(0,log(2))
error(gam(G=mixed,optimizer="efs",method="REML"),"all smoothing parameters")
linked<-gam(list(y~s(x,k=5,id=1)+s(z,k=5,id=1),~s(z,k=4)),data=d,family=f,method="REML",fit=FALSE)
error(gam(G=linked,optimizer="efs",method="REML"),"all smoothing parameters")
# Reduced maps remain available with an outer optimizer.
outer<-gam(G=linked,optimizer=c("outer","bfgs"),method="REML")
stopifnot(all(is.finite(coef(outer))))

error(gam(param,data=d,family=f,method="ML"),"supports REML")
error(gam(list(y~1,~1),data=data.frame(y=c(0,.5)),family=f),"integer responses")
error(gam(list(y~1,~1),data=data.frame(y=c(0,0)),family=f),"all-zero")
error(gam(param,data=d,family=f,weights=rep(0,n)),"positive total")
error(gam(param,data=d,family=f,weights=rep(-1,n)),"nonnegative")
error(gam(y~x,data=d,family=f),"two linear predictors")
error(bam(form,data=d,family=f,method="REML"),"general|gam\\(\\)")

# Optional gratia integration uses its native distribution factory protocol.
if(requireNamespace("gratia",quietly=TRUE)) {
  # Older gratia releases do not yet consume native factories.
  has_factory<-any(vapply(ls(asNamespace("gratia"),all.names=TRUE),function(nm) {
    obj<-get(nm,asNamespace("gratia"))
    is.function(obj) && any(grepl("gratia_diagnostics",deparse(body(obj)),fixed=TRUE))
  },logical(1)))
  if(has_factory) {
    pit<-gratia::quantile_residuals(b,type="pit",seed=11)
    qr<-gratia::quantile_residuals(b,type="quantile",seed=11)
    stopifnot(length(pit)==n,all(is.finite(qr)),all(pit>=0 & pit<=1))
    stopifnot(inherits(gratia::qq_plot(b,method="simulate",type="deviance",n_simulate=3,seed=11),"ggplot"))
    stopifnot(inherits(gratia::worm_plot(b,method="simulate",type="pearson",n_simulate=3,seed=11),"ggplot"))
    stopifnot(inherits(gratia::appraise(b,method="simulate",n_simulate=3,seed=11),"patchwork"))
  }
}
