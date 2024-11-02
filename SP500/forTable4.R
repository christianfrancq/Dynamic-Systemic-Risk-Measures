#############################  ..."
##########################  ..."
install.packages("xtable", repos='https://ftp.belnet.be/mirror/CRAN/')
library("xtable")
install.packages("compiler", repos='https://ftp.belnet.be/mirror/CRAN/')
library("compiler")
install.packages("VGAM", repos='https://ftp.belnet.be/mirror/CRAN/')
library("VGAM")
install.packages("nleqslv", repos='https://ftp.belnet.be/mirror/CRAN/')
library("nleqslv")
#######################################################
#############
# True critical value of the Unconditional coverage test
###############
LRuc.crit<-function(n,p,alpha0=0.05){
LR.val<-rep(0,(n+1))
prob<-dbinom(c(0:n),n,p)
for(T1 in 0:n){
  T0<-n-T1
  num<-0
  if((1-p)>0)num <- num + (log(1-p))*T0 
  if(p>0)num <- num + log(p)*T1
  den<-0
  if((1-T1/n)>0)den<-den + (log(1-T1/n))*T0 
  if((T1/n)>0)den<-den + log(T1/n)*T1
  LR.val[T1+1]<- -2*(num - den)
}
ordre<-order(LR.val)
crit<-(LR.val[ordre])[min(which(cumsum(prob[ordre])>=(1-alpha0)))]
crit.plus<-(LR.val[ordre])[min(which(cumsum(prob[ordre])>=(1-alpha0)))+1]
prob.0<-(prob[ordre])[min(which(cumsum(prob[ordre])>=(1-alpha0)))]


crit.moins<-(LR.val[ordre])[min(which(cumsum(prob[ordre])>=(1-alpha0)))-1]
eps.plus<-(crit.plus-crit)/2
eps.moins<-(crit-crit.moins)/2
#1-pchisq(crit,df=1) cumsum(LR.val[ordre])
niv<-sum(prob[which(LR.val>crit)])
prob.crit<-(alpha0-niv)/prob.0
crit.approx<-qchisq(1-alpha0, df=1)
niv.crit.approx<-sum(prob[which(LR.val>crit.approx)])

#
# For the Kupiec LR Test of Unconditional Coverage, under the null, the probability of LR>crit (equivalently LR>crit+eps.plus) 
#   is equal  to niv (<= alpha0, but the closest to alpha0).
# To obtain a randomized test of level exactly alpha0, it suffices to also reject 
#    with probability prob.crit when LR=crit (equivalently crit-eps.moins<LR<crit+eps.plus)
# With the critical value obtained from the chisq approximation (crit.approx), the probability of rejection is niv.crit.approx instead of alpha0
# 
list(crit=crit,niv=niv,eps.plus=eps.plus, eps.moins=eps.moins, prob.crit=prob.crit, crit.approx=crit.approx, niv.crit.approx=niv.crit.approx)
}
######
# Unconditional coverage testing  
###############
LRuc.Viol<-function(Viol,p){
  n<-length(Viol)
  T0<-length(which(Viol==0))
  T1<-n-T0
  num<-0
  if((1-p)>0)num <- num + (log(1-p))*T0 
  if(p>0)num <- num + log(p)*T1
  den<-0
  if((1-T1/n)>0)den<-den + (log(1-T1/n))*T0 
  if((T1/n)>0)den<-den + log(T1/n)*T1
  LR<- -2*(num - den)
  pval<-1-pchisq(LR,df=1)
  list(LR=LR,pval=pval)
}
######
# Independence testing
######
LRind.Viol<-function(Viol){
  n<-length(Viol)
  T01<-0
  T00<-0
  T11<-0
  T10<-0
  for(t in 2:n){
    if(Viol[t-1]==0&Viol[t]==0)T00<-T00+1
    if(Viol[t-1]==0&Viol[t]==1)T01<-T01+1
    if(Viol[t-1]==1&Viol[t]==0)T10<-T10+1
    if(Viol[t-1]==1&Viol[t]==1)T11<-T11+1}
  T0<-T00+T10
  T1<-T01+T11
  pi01<-1
  if((T00+T01)>0)pi01<-T01/(T00+T01)
  pi11<-1
  if((T10+T11)>0)pi11<-T11/(T10+T11)
  num<-0
  if((1-T1/n)>0)num<-num+(log(1-T1/n))*(T0/n)
  if((T1/n)>0)num<-num+log(T1/n)*(T1/n)
  #num<-(log(1-T1/n))*(T0/n) + log(T1/n)*(T1/n)
  den<-0
  if((1-pi01)>0)den<-den+(log(1-pi01))*(T00/n)
  if(pi01>0)den<-den+log(pi01)*(T01/n)
  if((1-pi11)>0)den<-den+(log(1-pi11))*(T10/n) 
  if(pi11>0)den<-den+log(pi11)*(T11/n)
  #den<-(log(1-pi01))*(T00/n) +log(pi01)*(T01/n)+ (log(1-pi11))*(T10/n) +log(pi11)*(T11/n)
  LR<- -2*n*(num - den)
  pval<-1-pchisq(LR,df=1)
  list(LR=LR,pval=pval)
}
################

##########
u.Norm<-function(rho,alpha,alpha.prim,u0){

Norm2<-function(x,y,rho){
(1/(2*pi*sqrt(1-rho^2)))*exp(-(x^2+y^2-2*rho*x*y)/(2*(1-rho^2)))
}

xi<-qnorm(alpha.prim)

f.y<-function(y,rho,xi){
z=y
for(i in 1:length(y)){
z[i]<-integrate(Norm2,-Inf, xi, y=y[i], rho=rho)$value
}
z
}


obj.u<-function(u,alpha,alpha.prim,rho,xi){
abs(integrate(f.y,-Inf, u,  rho=rho,xi=xi)$value-alpha*alpha.prim)
}


u<-nlminb(u0,obj.u,alpha=alpha,alpha.prim=alpha.prim,rho=rho,xi=xi)$par

list(xi=xi,u=u)
}

#####################

u.Stud<-function(rho,nu, alpha,alpha.prim,xi0,u0){

St2<-function(x,y,rho,nu){
(1/(2*pi*sqrt(1-rho^2)))*(1+(x^2+y^2-2*rho*x*y)/(nu*(1-rho^2)))^(-(nu+2)/2)
}

densy<-function(y,rho,nu){
z=y
for(i in 1:length(y)){
z[i]<-integrate(St2,-Inf, Inf, y=y[i], rho=rho, nu=nu)$value
}
z
}


obj.xi<-function(y,alpha.prim,rho,nu){
abs(integrate(densy,-Inf, y,  rho=rho, nu=nu)$value-alpha.prim)
}

if(nu!=Inf){
xi<-nlminb(xi0,obj.xi,alpha.prim=alpha.prim,rho=rho,nu=nu)$par
}else{
xi<-qnorm(alpha.prim)
}

f.y<-function(y,rho,nu,xi){
z=y
for(i in 1:length(y)){
z[i]<-integrate(St2,-Inf, xi, y=y[i], nu=nu, rho=rho)$value
}
z
}





obj.u<-function(u,alpha,alpha.prim,rho,nu,xi){
abs(integrate(f.y,-Inf, u,  rho=rho, nu=nu, xi=xi)$value-alpha*alpha.prim)
}

if(nu!=Inf){
u<-nlminb(u0,obj.u,alpha=alpha,alpha.prim=alpha.prim,rho=rho,nu=nu, xi=xi)$par
cte<-sqrt((nu-2)/nu)
}else{
u<-u.Norm(rho,alpha,alpha.prim,u0)$u
cte<-1
}

list(xi=cte*xi,u=cte*u)
}

##### 




# racine carree d'une matrice symmetrique semi-definie positive
Sqrt <- function(Sigma) {
  n<- nrow(Sigma)
  sp<-eigen(Sigma)
  sp$vector%*%sqrt(diag(abs(sp$values)))%*%t(sp$vector)
}


# Inverse (g�n�ralis�e) d'une matrice symmetrique semi-definie positive
Inv <- function(Sigma,tol=sqrt(.Machine$double.eps)) {
  n<- nrow(Sigma)
  sp<-eigen(Sigma)
  sp$vector%*%diag(1/pmax(abs(sp$values),tol))%*%t(sp$vector)
}

######################################################################
# simule un GARCH bivarié avec bruit de Student (ou Gaussien) de correlation rho
######################################################################
#BivGarchSt.sim(n+P, omega,  mat.alpha, beta, rho, u.eta, nu)
BivGarchSt.sim<- function (n, omega,  mat.alpha, beta, rho, u.eta, nu, valinit=500) {
#
norm<-matrix(rnorm(2*(n+valinit)),nrow=(n+valinit),ncol=2)
norm[,2]<-rho*norm[,1]+sqrt(1-rho^2)*norm[,2]
#var(norm) plot(density(norm))
{if(nu==Inf|nu<2){
ct<-1}
else
{ct<-sqrt((nu-2)/nu)
}}
if(nu==Inf){
eta<-norm
}else{
eta<-ct*sqrt(nu/rchisq((n+valinit),df=nu))*norm}
#var(eta) plot(density(eta)) var(eta)

ht<-matrix(nrow=(n+valinit),ncol=2)
eps<-ht


ht[1,]<-omega
eps[1,]<-sqrt(ht[1,])*eta[1,]
for(t in 2:(n+valinit)){
ht[t,]<-omega+mat.alpha%*%(eps[t-1,]^2)+beta*ht[t-1,]
eps[t,]<-sqrt(ht[t,])*eta[t,]
}

#plot.ts(eps)

list(simu=eps[(valinit+1):(n+valinit),], vol=sqrt(ht[(valinit+1):(n+valinit),]))
}

BivGarchSt.sim<-cmpfun(BivGarchSt.sim)


##########################################################################
# estime un GARCH bivarié par QMLE et backtest 2 CoVaR
#########################################################################
#pour estimer le modèle de la premiere composante
objf.garch.biv1 <- function(vartheta, eps2, n, sig2init, r0=10){     
         omega <- vartheta[1]
         mat.alpha <- vartheta[2:3]
         beta <- vartheta[4]
sig2<-rep(0,n)
sig2[1]<-sig2init
for(t in 2:n){
sig2[t]<-omega+mat.alpha[1]*(eps2[t-1,1])+mat.alpha[2]*(eps2[t-1,2])+beta*sig2[t-1]
}
qml <- mean(eps2[(r0+1):n,1]/sig2[(r0+1):n]+log(sig2[(r0+1):n]))
qml }

objf.garch.biv1<-cmpfun(objf.garch.biv1)

#pour estimer le modèle de la seconde composante
objf.garch.biv2 <- function(vartheta, eps2, n, sig2init, r0=10){     
         omega <- vartheta[1]
         mat.alpha <- vartheta[2:3]
         beta <- vartheta[4]
sig2init<-mean(eps2[1:r0,2])
sig2<-rep(0,n)
sig2[1]<-sig2init
for(t in 2:n){
sig2[t]<-omega+mat.alpha[1]*(eps2[t-1,1])+mat.alpha[2]*(eps2[t-1,2])+beta*sig2[t-1]
}
qml <- mean(eps2[(r0+1):n,2]/sig2[(r0+1):n]+log(sig2[(r0+1):n]))
qml }

objf.garch.biv2<-cmpfun(objf.garch.biv2)

#res<-estimGARCHbiv(omega.init,mat.alpha.init,beta.init,sim)
# eps<-res.sim$simu; omega.init<-omega; mat.alpha.init<-mat.alpha; beta.init<-beta
##############################BacktestGARCHbiv.mod(omega,mat.alpha,beta,res.sim$simu,niv.alpha,niv.alpha.prim, P, crit.UC, crit.Ind)
 

BacktestGARCHbiv.mod<- function(omega.init,mat.alpha.init,beta.init,eps, niv.alpha, niv.alpha.prim, P, crit.UC, eps.plus, eps.moins, prob.crit, crit.Ind, eps.plus.ind, eps.moins.ind, prob.crit.ind,
             petit=sqrt(.Machine$double.eps), r0=10, alpha0=0.05)
     {
n <- length(eps[,1])-P
omega<-c(0,0)
beta<-omega
mat.alpha<-matrix(nrow=2,ncol=2)
valinit<-c(omega.init[1],mat.alpha.init[1,1:2],beta.init[1])
eps2<-eps^2
sig2init1<-mean(eps2[1:r0,1])

res <- nlminb(valinit,objf.garch.biv1,lower=c(petit^2,0,0,0),
   upper=c(Inf,Inf,Inf,1),eps2=eps2[1:n,],n=n,sig2init=sig2init1)
         omega[1] <- res$par[1]
         mat.alpha[1,1:2]<- res$par[2:3]
         beta[1] <- res$par[4]

valinit<-c(omega.init[2],mat.alpha.init[2,1:2],beta.init[2])
sig2init2<-mean(eps2[1:r0,2])

res <- nlminb(valinit,objf.garch.biv2,lower=c(petit^2,0,0,0),
   upper=c(Inf,Inf,Inf,1),eps2=eps2[1:n,],n=n,sig2init=sig2init2)
         omega[2] <- res$par[1]
         mat.alpha[2,1:2]<- res$par[2:3]
         beta[2] <- res$par[4]

#omega<-omega.init; mat.alpha<-mat.alpha.init; beta<-beta.init

ht<-matrix(nrow=(n+P),ncol=2)
eta<-eps

ht[1,]<-c(sig2init1,sig2init2)
eta[1,]<-eps[1,]/sqrt(pmax(ht[1,],petit)) 
for(t in 2:(n+P)){
ht[t,]<-omega+mat.alpha%*%eps2[t-1,]+beta*ht[t-1,]
eta[t,]<-eps[t,]/sqrt(pmax(ht[t,],petit))
}
################ semi-parametric method
#
xi2<-as.numeric(quantile(eta[1:n,2],probs=niv.alpha.prim,type=1))
select<-which(eta[1:n,2]<=xi2)
u.eta<-as.numeric(quantile(eta[select,1],probs=niv.alpha,type=1))
CoVaR<-rep(0,P)
VaR<-CoVaR
Viol<-CoVaR
for(i in 1:P){
VaR[i] <- -sqrt(ht[n+i,2])*xi2
CoVaR[i] <- -sqrt(ht[n+i,1])*u.eta
if((eps[n+i,1] < (-CoVaR[i]))&(eps[n+i,2] < (-VaR[i])))Viol[i]<-1 
}
#
# For the Kupiec LR Test of Unconditional Coverage, reject if UC.LR>crit.UC (equivalently UC.LR>crit.UC+eps.plus) 
#   and also reject  with probability prob.crit when UC.LR=crit.UC (equivalently crit.UC-eps.moins<UC.LR<crit.UC+eps.plus)
# 

UC.LR<-LRuc.Viol(Viol,niv.alpha*niv.alpha.prim)$LR
UC<-0
if(UC.LR>crit.UC+eps.plus)UC<-1
if((UC.LR>crit.UC-eps.moins)&(UC.LR<crit.UC+eps.plus))UC<-rbinom(1,1,prob.crit)

#
# For the Christoffersen LR Test of Independence, reject if Ind.LR>crit.Ind (equivalently Ind.LR>crit.Ind+eps.plus.ind) 
#   and also reject  with probability prob.crit.ind when Ind.LR=crit.Ind (equivalently crit.Ind-eps.moins.ind<Ind.LR<crit.Ind+eps.plus.ind)
# 

Ind.LR<-LRind.Viol(Viol)$LR
Ind<-0
if(Ind.LR>crit.Ind+eps.plus.ind)Ind<-1
if((Ind.LR>crit.Ind-eps.moins.ind)&(Ind.LR<crit.Ind+eps.plus.ind))Ind<-rbinom(1,1,prob.crit.ind)




################ Gaussian method
rho<-cor(eta)[1,2]
u.eta.norm<-u.Norm(rho,niv.alpha,niv.alpha.prim,-2.5)$u
xi.norm<-qnorm(niv.alpha.prim)
CoVaR.norm<-rep(0,P)
VaR.norm<-CoVaR.norm
Viol.norm<-CoVaR.norm
for(i in 1:P){
VaR.norm[i]<- -sqrt(ht[n+i,2])*xi.norm
CoVaR.norm[i]<- -sqrt(ht[n+i,1])*u.eta.norm
if((eps[n+i,1]< (-CoVaR.norm[i]))&(eps[n+i,2]< (-VaR.norm[i])))Viol.norm[i]<-1 
}



#
# For the Kupiec LR Test of Unconditional Coverage, reject if UCtest.norm>crit.UC (equivalently UCtest.norm>crit.UC+eps.plus) 
#   and also reject  with probability prob.crit when UCtest.norm=crit.UC (equivalently crit.UC-eps.moins<UCtest.norm<crit.UC+eps.plus)
# 


UCtest.norm<-LRuc.Viol(Viol.norm,niv.alpha*niv.alpha.prim)$LR
UC.norm<-0
if(UCtest.norm>crit.UC+eps.plus)UC.norm<-1
if((UCtest.norm>crit.UC-eps.moins)&(UCtest.norm<crit.UC+eps.plus))UC.norm<-rbinom(1,1,prob.crit)


#
# For the Christoffersen LR Test of Independence, reject if Indtest.norm>crit.Ind (equivalently Indtest.norm>crit.Ind+eps.plus.ind) 
#   and also reject  with probability prob.crit.ind when Indtest.norm=crit.Ind (equivalently crit.Ind-eps.moins.ind<Indtest.norm<crit.Ind+eps.plus.ind)
# 

Indtest.norm<-LRind.Viol(Viol.norm)$LR
Ind.norm<-0
if(Indtest.norm>crit.Ind+eps.plus.ind)Ind.norm<-1
if((Indtest.norm>crit.Ind-eps.moins.ind)&(Indtest.norm<crit.Ind+eps.plus.ind))Ind.norm<-rbinom(1,1,prob.crit.ind)



#####################
list(UC=UC,Ind=Ind,UC.norm=UC.norm,Ind.norm=Ind.norm,u.eta=u.eta,u.eta.norm=u.eta.norm)
}

BacktestGARCHbiv.mod<-cmpfun(BacktestGARCHbiv.mod)

##################simulation.mod(n, P, omega,  mat.alpha, beta, rho, u.eta, crit.UC, eps.plus, eps.moins, prob.crit, crit.Ind, eps.plus.ind, eps.moins.ind, prob.crit.ind, nu)

simulation.mod<- function(n, P, omega,  mat.alpha, beta, rho, u.eta, crit.UC, eps.plus, eps.moins, prob.crit, crit.Ind, eps.plus.ind, eps.moins.ind, prob.crit.ind, nu) {

res.sim<-BivGarchSt.sim(n+P, omega,  mat.alpha, beta, rho, u.eta, nu)
res<-BacktestGARCHbiv.mod(omega,mat.alpha,beta,res.sim$simu,niv.alpha,niv.alpha.prim, P, crit.UC, eps.plus, eps.moins, prob.crit, crit.Ind, eps.plus.ind, eps.moins.ind, prob.crit.ind)
c(res$UC,res$Ind,res$UC.norm,res$Ind.norm,res$u.eta,res$u.eta.norm)

}
simulation.mod<-cmpfun(simulation.mod)


##################
##################simul(n, P, niv.alpha, niv.alpha.prim ,alpha0, crit.UC)
simul<- function(n, P, niv.alpha, niv.alpha.prim,alpha0,crit.UC, eps.plus, eps.moins, prob.crit) {
Viol<-rep(0,P)
for(i in 1:P){
Viol[i]<-rbinom(1, 1, niv.alpha*niv.alpha.prim)
}
test<-LRuc.Viol(Viol,niv.alpha*niv.alpha.prim)
UCtest<-test$pval
UC<-0
if(UCtest<alpha0)UC<-1
UC.LR<-test$LR
UC.mod<-0
if(UC.LR>crit.UC+eps.plus/2)UC.mod<-1
if((UC.LR>crit.UC-eps.moins/2)&(UC.LR<crit.UC+eps.plus/2))UC.mod<-rbinom(1,1,prob.crit)
test<-LRind.Viol(Viol)
Indtest<-test$pval
Ind.LR<-test$LR
Ind<-0
if(Indtest<alpha0)Ind<-1

c(UC,Ind,UC.mod,UC.LR,Ind.LR)

}
simul<-cmpfun(simul)

##################simul(n, P, niv.alpha, niv.alpha.prim ,alpha0, crit.UC)
simula<- function(n, P, niv.alpha, niv.alpha.prim,alpha0,crit.UC, eps.plus, eps.moins, prob.crit, crit.Ind, eps.plus.ind, eps.moins.ind, prob.crit.ind) {
Viol<-rep(0,P)
for(i in 1:P){
Viol[i]<-rbinom(1, 1, niv.alpha*niv.alpha.prim)
}
######### test UC
test<-LRuc.Viol(Viol,niv.alpha*niv.alpha.prim)
UCtest<-test$pval
UC<-0
if(UCtest<alpha0)UC<-1
# modif
UC.LR<-test$LR
UC.mod<-0
if(UC.LR>crit.UC+eps.plus)UC.mod<-1
if((UC.LR>crit.UC-eps.moins)&(UC.LR<crit.UC+eps.plus))UC.mod<-rbinom(1,1,prob.crit)
############ test Ind
test<-LRind.Viol(Viol)
Indtest<-test$pval
Ind.LR<-test$LR
Ind<-0
if(Indtest<alpha0)Ind<-1
# modif
Ind.mod<-0
if(Ind.LR>crit.Ind+eps.plus.ind)Ind.mod<-1
if((Ind.LR>crit.Ind-eps.moins.ind)&(Ind.LR<crit.Ind+eps.plus.ind))Ind.mod<-rbinom(1,1,prob.crit.ind)

c(UC,Ind,UC.mod,Ind.mod,UC.LR,Ind.LR)

}
simula<-cmpfun(simula)

##################BacktestGARCHbiv(omega,mat.alpha,beta,res.sim$simu,niv.alpha,niv.alpha.prim)
#omega.init<-omega; mat.alpha.init<-mat.alpha; beta.init<-beta;eps<-res.sim$simu



BacktestGARCHbiv<- function(omega.init,mat.alpha.init,beta.init,eps, niv.alpha, niv.alpha.prim, P,
             petit=sqrt(.Machine$double.eps), r0=10, alpha0=0.05)
     {
n <- length(eps[,1])-P
omega<-c(0,0)
beta<-omega
mat.alpha<-matrix(nrow=2,ncol=2)
valinit<-c(omega.init[1],mat.alpha.init[1,1:2],beta.init[1])
eps2<-eps^2
sig2init1<-mean(eps2[1:r0,1])

res <- nlminb(valinit,objf.garch.biv1,lower=c(petit^2,0,0,0),
   upper=c(Inf,Inf,Inf,1),eps2=eps2[1:n,],n=n,sig2init=sig2init1)
         omega[1] <- res$par[1]
         mat.alpha[1,1:2]<- res$par[2:3]
         beta[1] <- res$par[4]

valinit<-c(omega.init[2],mat.alpha.init[2,1:2],beta.init[2])
sig2init2<-mean(eps2[1:r0,2])

res <- nlminb(valinit,objf.garch.biv2,lower=c(petit^2,0,0,0),
   upper=c(Inf,Inf,Inf,1),eps2=eps2[1:n,],n=n,sig2init=sig2init2)
         omega[2] <- res$par[1]
         mat.alpha[2,1:2]<- res$par[2:3]
         beta[2] <- res$par[4]

#omega<-omega.init; mat.alpha<-mat.alpha.init; beta<-beta.init

ht<-matrix(nrow=(n+P),ncol=2)
eta<-eps

ht[1,]<-c(sig2init1,sig2init2)
eta[1,]<-eps[1,]/sqrt(pmax(ht[1,],petit)) 
for(t in 2:(n+P)){
ht[t,]<-omega+mat.alpha%*%eps2[t-1,]+beta*ht[t-1,]
eta[t,]<-eps[t,]/sqrt(pmax(ht[t,],petit))
}
################ semi-parametric method
#
xi2<-as.numeric(quantile(eta[1:n,2],probs=niv.alpha.prim,type=1))
select<-which(eta[1:n,2]<=xi2)
u.eta<-as.numeric(quantile(eta[select,1],probs=niv.alpha,type=1))
CoVaR<-rep(0,P)
VaR<-CoVaR
Viol<-CoVaR
for(i in 1:P){
VaR[i] <- -sqrt(ht[n+i,2])*xi2
CoVaR[i] <- -sqrt(ht[n+i,1])*u.eta
if((eps[n+i,1] < (-CoVaR[i]))&(eps[n+i,2] < (-VaR[i])))Viol[i]<-1 
}
UCtest<-LRuc.Viol(Viol,niv.alpha*niv.alpha.prim)$pval
UC<-0
if(UCtest<alpha0)UC<-1
Indtest<-LRind.Viol(Viol)$pval
Ind<-0
if(Indtest<alpha0)Ind<-1


################ Gaussian method
rho<-cor(eta)[1,2]
u.eta.norm<-u.Norm(rho,niv.alpha,niv.alpha.prim,-2.5)$u
xi.norm<-qnorm(niv.alpha.prim)
CoVaR.norm<-rep(0,P)
VaR.norm<-CoVaR.norm
Viol.norm<-CoVaR.norm
for(i in 1:P){
VaR.norm[i]<- -sqrt(ht[n+i,2])*xi.norm
CoVaR.norm[i]<- -sqrt(ht[n+i,1])*u.eta.norm
if((eps[n+i,1]< (-CoVaR.norm[i]))&(eps[n+i,2]< (-VaR.norm[i])))Viol.norm[i]<-1 
}
UCtest.norm<-LRuc.Viol(Viol.norm,niv.alpha*niv.alpha.prim)$pval
UC.norm<-0
if(UCtest.norm<alpha0)UC.norm<-1
Indtest.norm<-LRind.Viol(Viol.norm)$pval
Ind.norm<-0
if(Indtest.norm<alpha0)Ind.norm<-1

#####################
list(UC=UC,Ind=Ind,UC.norm=UC.norm,Ind.norm=Ind.norm,u.eta=u.eta,u.eta.norm=u.eta.norm)
}

BacktestGARCHbiv<-cmpfun(BacktestGARCHbiv)

###########################
simulation<- function(n, P, omega,  mat.alpha, beta, rho, u.eta, nu) {

res.sim<-BivGarchSt.sim(n+P, omega,  mat.alpha, beta, rho, u.eta, nu)
#str(res.sim)
res<-BacktestGARCHbiv(omega,mat.alpha,beta,res.sim$simu,niv.alpha,niv.alpha.prim,P)
c(res$UC,res$Ind,res$UC.norm,res$Ind.norm,res$u.eta,res$u.eta.norm)

}
simulation<-cmpfun(simulation)

############################## 
####################################################

omega<-c(0.001,0.001); mat.alpha<-matrix(c(0.05,0.01,0.01,0.1),ncol=2,nrow=2); beta<-c(0.9,0.85) 
rho<-0.6; nu=6; niv.alpha<-0.05; niv.alpha.prim<-0.1
P<-250

u.eta<-u.Stud(rho,nu, niv.alpha,niv.alpha.prim,-1,-2)$u
u.eta.norm<-u.Norm(rho,niv.alpha,niv.alpha.prim,-2.5)$u
###########################
Nrep<-20000
p<-niv.alpha*niv.alpha.prim
testUC<-LRuc.crit(P,p)
crit.UC<-testUC$crit
eps.plus<-testUC$eps.plus 
eps.moins<-testUC$eps.moins 
prob.crit<-testUC$prob.crit
alpha0<-0.05
#qchisq(0.95, df=1)
set.seed(13)

res<-replicate(Nrep,simul(n, P, niv.alpha, niv.alpha.prim ,alpha0, crit.UC, eps.plus, eps.moins, prob.crit))

str(res)
unique(sort(res[4,]))


length(which(res[1,]==1))*100/Nrep
length(which(res[2,]==1))*100/Nrep
length(which(res[3,]==1))*100/Nrep # bon niveau
crit.Ind<-as.numeric(quantile(res[5,],1-alpha0,type=1))
eps.moins.ind<-(crit.Ind-unique(sort(res[5,]))[which(unique(sort(res[5,]))==crit.Ind)-1])/2
eps.plus.ind<-(-crit.Ind+unique(sort(res[5,]))[which(unique(sort(res[5,]))==crit.Ind)+1])/2
prob.crit.ind<-(alpha0-length(which(res[5,]>crit.Ind+eps.plus.ind))/Nrep)/(length(which((res[5,]<crit.Ind+eps.plus.ind)&(res[5,]>crit.Ind-eps.moins.ind)))/Nrep)

####################################
x<-read.table(file="complete_stocks")
x$x[1]
nb.stocks<-length(x$x)
res<-matrix(nrow=6,ncol=nb.stocks)
for(i.stock in 1:nb.stocks){
serie<-read.csv(file=paste0(x$x[i.stock], "_SP500.csv")) #str(serie)
series<-drop(as.matrix(serie[,2:3])) #str(series) series[1:5,]
n<-length(series[,1])
rend<-matrix(nrow=(n-1),ncol=2)
for(i in 1:(n-1)){
rend[i,1]<-series[i+1,1]/series[i,1]-1
#rend[i,1]<-100*log(series[i+1,1]/series[i,1])
rend[i,2]<-series[i+1,2]/series[i,2]-1
#rend[i,2]<-100*log(series[i+1,2]/series[i,2])
}
n<-length(rend[,1])-P
#print(n)
tests<-BacktestGARCHbiv.mod(omega,mat.alpha,beta,rend,niv.alpha,niv.alpha.prim, P, crit.UC, eps.plus, eps.moins, prob.crit, crit.Ind, eps.plus.ind, eps.moins.ind, prob.crit.ind)
res[,i.stock]<-c(tests$UC,tests$Ind,tests$UC.norm,tests$Ind.norm,tests$u.eta,tests$u.eta.norm)
}



a<-round(length(which(res[1,]==1))/nb.stocks,digits=3)
b<-round(length(which(res[2,]==1))/nb.stocks,digits=3)
c<-round(length(which(res[3,]==1))/nb.stocks,digits=3)
d<-round(length(which(res[4,]==1))/nb.stocks,digits=3)

for_tab<-data.frame(
			"Model"=c("Parametric Gaussian","Semi-Parametric"),
"Unconditional Coverage Test"=c(c,a),
"Independence Test"=c(d,b)
)
table_xtable <- xtable(
  for_tab,
)
table_xtable 

