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

##### 
#rho<-0.5; alpha<-0.1; alpha.prim<-0.1; u0<--2
#u.Norm(rho,alpha,alpha.prim,u0)
#n<-10000000
#norm<-matrix(rnorm(2*n),nrow=(n),ncol=2)
#norm[,2]<-rho*norm[,1]+sqrt(1-rho^2)*norm[,2]
#cov(norm)
#quantile(norm[,2],alpha.prim)
#indic<-which(norm[,2]<=qnorm(alpha.prim))
#quantile(norm[indic,1],alpha)

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
#rho<-0.5; alpha<-0.05; alpha.prim<-0.1; nu<-6; xi0<--1; u0<--2
#u.Stud(rho,nu, alpha,alpha.prim,xi0,u0)
#n<-10000000
#norm<-matrix(rnorm(2*n),nrow=n,ncol=2)
#norm[,2]<-rho*norm[,1]+sqrt(1-rho^2)*norm[,2]
#{if(nu==Inf|nu<2){
#ct<-1}
#else
#{ct<-sqrt((nu-2)/nu)
#}}
#if(nu==Inf){
#eta<-norm
#}else{
#eta<-ct*sqrt(nu/rchisq(n,df=nu))*norm}
##var(eta) 
#quantile(eta[,2],alpha.prim)
#xi<-quantile(eta[,2],alpha.prim)
#indic<-which(eta[,2]<=xi)
#quantile(eta[indic,1],alpha)





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

# Inverse (g�n�ralis�e) d'une racine carr�e d'une matrice symmetrique semi-definie positive
InvSqrt <- function(Sigma,tol=sqrt(.Machine$double.eps)) {
  n<- nrow(Sigma)
  sp<-eigen(Sigma)
  sp$vector%*%sqrt(diag(1/pmax(abs(sp$values),tol)))%*%t(sp$vector)
}

# # racine carree d'un matrice symmetrique semi-definie positive
# Sqrt <- function(Sigma) {
#   n<- nrow(Sigma)
#   sp<-eigen(Sigma)
#   sp$vector%*%sqrt(diag(sp$values))
# }

######################################################################
# simule un GARCH bivarié avec bruit de Student (ou Gaussien) de correlation rho
######################################################################
#sim0<-Garch.sim(n,omeg0,alph0,bet0) omega<-omeg0;mat.alpha<-alph0;beta<-bet0
BivGarchSt.sim<- function (n, omega,  mat.alpha, beta, rho, u.eta, nu, valinit=500) {
#omega et beta vecteurs de dimension 2, mat.alpha matrice 2x2
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
volat<-sqrt(omega+mat.alpha%*%eps[n+valinit,]^2+beta*ht[n+valinit,])
CoVaR.nplus1<-(-volat*u.eta)[1,1]

list(simu=eps[(valinit+1):(n+valinit),], vol=sqrt(ht[(valinit+1):(n+valinit),]), volat=volat, true.CoVaR=CoVaR.nplus1)
}

BivGarchSt.sim<-cmpfun(BivGarchSt.sim)


##########################################################################
# estime un GARCH bivarié par QMLE
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

#res<-estimGARCHbiv(omega,mat.alpha,beta,res.sim$simu,niv.alpha,niv.alpha.prim,res.sim$true.CoVaR)
# eps<-res.sim$simu; omega.init<-omega; mat.alpha.init<-mat.alpha; beta.init<-beta; true.CoVaR<-res.sim$true.CoVaR

#res<-estimGARCHbiv(omega,mat.alpha,beta,res.sim$simu,niv.alpha,niv.alpha.prim,res.sim$true.CoVaR)


estimGARCHbiv<- function(omega.init,mat.alpha.init,beta.init,eps, niv.alpha, niv.alpha.prim, true.CoVaR, 
             petit=sqrt(.Machine$double.eps), r0=10, alpha0=0.05, B=500, nb_try=5)
     {
n <- length(eps[,1])
omega<-c(0,0)
beta<-omega
mat.alpha<-matrix(nrow=2,ncol=2)
valinit<-c(omega.init[1],mat.alpha.init[1,1:2],beta.init[1])
eps2<-eps^2
sig2init1<-mean(eps2[1:r0,1])

res <- nlminb(valinit,objf.garch.biv1,lower=c(petit^2,0,0,0),
   upper=c(Inf,Inf,Inf,1),eps2=eps2,n=n,sig2init=sig2init1)
         omega[1] <- res$par[1]
         mat.alpha[1,1:2]<- res$par[2:3]
         beta[1] <- res$par[4]

valinit<-c(omega.init[2],mat.alpha.init[2,1:2],beta.init[2])
sig2init2<-mean(eps2[1:r0,2])

res <- nlminb(valinit,objf.garch.biv2,lower=c(petit^2,0,0,0),
   upper=c(Inf,Inf,Inf,1),eps2=eps2,n=n,sig2init=sig2init2)
         omega[2] <- res$par[1]
         mat.alpha[2,1:2]<- res$par[2:3]
         beta[2] <- res$par[4]

#omega<-omega.init; mat.alpha<-mat.alpha.init; beta<-beta.init

ht<-matrix(0,nrow=n,ncol=2)
eta<-eps
der.sig2.1<-matrix(0,nrow=n,ncol=4)
der.sig2.2<-der.sig2.1

ht[1,]<-c(sig2init1,sig2init2)
eta[1,]<-eps[1,]/sqrt(pmax(ht[1,],petit)) 
for(t in 2:n){
ht[t,]<-omega+mat.alpha%*%eps2[t-1,]+beta*ht[t-1,]
eta[t,]<-eps[t,]/sqrt(pmax(ht[t,],petit))
der.sig2.1[t,]<-c(1,eps2[t-1,1],eps2[t-1,2],ht[t-1,1])+beta[1]*der.sig2.1[t-1,]
der.sig2.2[t,]<-c(1,eps2[t-1,1],eps2[t-1,2],ht[t-1,2])+beta[2]*der.sig2.2[t-1,]
}

# var(eta)
s<-sqrt(diag(var(eta))) 
for(t in 1:n){
eta[t,1]<-(eta[t,1]-mean(eta[,1]))/s[1]
eta[t,2]<-(eta[t,2]-mean(eta[,2]))/s[2]
}
#var(eta) plot.ts(eta)


ht.nplus1<-omega+mat.alpha%*%eps2[n,]+beta*ht[n,]
################ semi-parametric method
#colMeans(eta^2)
D1<-der.sig2.1/(2*ht[,1])
J1<-t(D1)%*%D1/n
Omega1<-colMeans(D1)
#D2<-der.sig2.2/(2*ht[,2])
#J2<-t(D2)%*%D2/n
#kappa1<-mean(eta[,1]^4)
J1.inv<-Inv(J1)
#Sigma1<-((kappa1-1)/4)*J1.inv
#std<-sqrt(diag(Sigma1/n))
xi2<-as.numeric(quantile(eta[,2],probs=niv.alpha.prim,type=1))
select<-which(eta[,2]<=xi2)
u.eta<-as.numeric(quantile(eta[select,1],probs=niv.alpha,type=1))
#CoVaR<--sqrt(ht[,1])*u.eta
CoVaR.nplus1<--sqrt(ht.nplus1[1])*u.eta
delta<-rep(0,5)
delta[1:4]<-u.eta*(c(1,eps2[n,1],eps2[n,2],ht[n,1])+beta[1]*der.sig2.1[n,])/(sqrt(ht.nplus1[1])*2)
delta[5]<-sqrt(ht.nplus1[1])
Upsilon<-matrix(0,nrow=n,ncol=3)
select<-which(eta[,1] <= u.eta & eta[,2] <= xi2)
Upsilon[select,1]<-1
select<-which( eta[,2] <= xi2)
Upsilon[select,2]<-1
Upsilon[,3]<-eta[,1]^2
Sigma.Upsilon<-var(Upsilon)
#Upsilon[,1]<-Upsilon[,1]-niv.alpha*niv.alpha.prim
#Upsilon[,2]<-Upsilon[,2]-niv.alpha.prim
#Upsilon[,3]<-Upsilon[,3]-1

densite<-density(eta[,2])
g2.xi2<-densite$y[which.min(abs(densite$x-xi2))]
densite<-density(eta[select,1])
f1.u.xi2<-densite$y[which.min(abs(densite$x-u.eta))]
select<-which(eta[,1]<=u.eta)
densite<-density(eta[select,2])
f2.xi2.u<-densite$y[which.min(abs(densite$x-xi2))]
#G.1<-ecdf(eta[,1])
#G.1(u.eta)
G.1.u<-length(which(eta[,1]<=u.eta))/n
lambda<-c(-1/(niv.alpha.prim*f1.u.xi2),(G.1.u*f2.xi2.u)/(niv.alpha.prim*f1.u.xi2*g2.xi2),-u.eta/2)
#select<-which(eta[,1]<=u.eta)
#indic<-rep(0,n)
#indic[select]<-1
#select<-which(eta[,2]<=xi2)
#varrho<-mean((eta[,1]^2*indic)[select])-niv.alpha
#Nabla<-mean((eta[,1]^2)[select])-1
#
#Sigma.Up<-matrix(ncol=3,nrow=3)
#Sigma.Up[1,1]<-niv.alpha*niv.alpha.prim*(1-niv.alpha*niv.alpha.prim)
#Sigma.Up[1,2]<-niv.alpha*niv.alpha.prim*(1-niv.alpha.prim)
#Sigma.Up[2,1]<-Sigma.Up[1,2]
#Sigma.Up[1,3]<-niv.alpha.prim*varrho
#Sigma.Up[3,1]<-Sigma.Up[1,3]
#Sigma.Up[2,2]<-niv.alpha.prim*(1-niv.alpha.prim)
#Sigma.Up[2,3]<-niv.alpha.prim*Nabla
#Sigma.Up[3,2]<-Sigma.Up[2,3]
#Sigma.Up[3,3]<-kappa1-1
#
#Sigma.alpha.alpha.prim<-matrix(ncol=5,nrow=5)
#Sigma.alpha.alpha.prim[1:4,1:4]<-Sigma1
##eigen(Sigma.Upsilon) eigen(Sigma1)
#
#Sigma.alpha.alpha.prim[5,5]<-drop(t(lambda)%*%Sigma.Upsilon%*%lambda)
#Sigma.alpha.alpha.prim[1:4,5]<-J1.inv%*%Omega1%*%c(0,0,1)%*%Sigma.Up%*%lambda/2
#Sigma.alpha.alpha.prim[5,1:4]<-Sigma.alpha.alpha.prim[1:4,5]
#
bigUpsilon<-matrix(0,nrow=n,ncol=5)
#bigUpsilon[,1:4]<-(D1%*%J1.inv*eta[,1]^2)/2
bigUpsilon[,1:4]<-(D1%*%J1.inv*(eta[,1]^2-1))/2
bigUpsilon[,5]<-Upsilon%*%lambda

Sigma.alpha.alpha.prim<-var(bigUpsilon)

sd<-as.numeric(sqrt(delta%*%Sigma.alpha.alpha.prim%*%delta/n))
lower<-CoVaR.nplus1-qnorm(1-alpha0/2)*sd
upper<-CoVaR.nplus1+qnorm(1-alpha0/2)*sd

IW<-abs((upper-lower)/CoVaR.nplus1)
Cover<-0
if(true.CoVaR >= lower & true.CoVaR <= upper)Cover<-1 


################ Gaussian method
rho<-cor(eta)[1,2]
u.eta.norm<-u.Norm(rho,niv.alpha,niv.alpha.prim,-2.5)$u
xi.norm<-qnorm(niv.alpha.prim)
z.rho<-(xi.norm-rho*u.eta.norm)/sqrt(1-rho^2)
K.rho<-dnorm(z.rho)/pnorm(z.rho)
Sigma.N<-matrix(nrow=5,ncol=5)
Sigma.N[1:4,1:4]<-J1.inv/2
Sigma.N[1:4,5]<-(-rho*sqrt(1-rho^2)*K.rho/2)*J1.inv%*%Omega1
Sigma.N[5,1:4]<-Sigma.N[1:4,5]
Sigma.N[5,5]<-K.rho^2*(1-rho^2)
#eigen(Sigma.N)
delta[1:4]<-u.eta.norm*(c(1,eps2[n,1],eps2[n,2],ht[n,1])+beta[1]*der.sig2.1[n,])/(sqrt(ht.nplus1[1])*2)
sd.norm<-as.numeric(sqrt(delta%*%Sigma.N%*%delta/n))
CoVaR.nplus1.norm<--sqrt(ht.nplus1[1])*u.eta.norm
lower.norm<-CoVaR.nplus1.norm-qnorm(1-alpha0/2)*sd.norm
upper.norm<-CoVaR.nplus1.norm+qnorm(1-alpha0/2)*sd.norm
IW.norm<-abs((upper.norm-lower.norm)/CoVaR.nplus1.norm)
Cover.norm<-0
if(true.CoVaR >= lower.norm & true.CoVaR <= upper.norm)Cover.norm<-1 

########################## Bootstrap method
hat.alpha.prim<-length(which(eta[,2]<=xi2))/n
hat.alpha<-(length(which((eta[,2]<=xi2)&(eta[,1]<=u.eta)))/n)/hat.alpha.prim
hat.m2<-mean(eta[,1]^2)
CoVar.star<-rep(0,B)
sig2<-rep(0,n)
sig2[1]<-sig2init1
theta1<-c(omega[1],mat.alpha[1,1:2],beta[1])
for(b in 1:B){
bon<-0
try<-0
while(bon==0){
try<-try+1
if(try>nb_try)bon<-1
indic<-sample(c(1:n),replace=TRUE)
eta.star<-eta[indic,]
#var(eta.star)
Upsilon.star<-matrix(0,nrow=n,ncol=3)
Upsilon.star[which((eta.star[,2]<=xi2)&(eta.star[,1]<=u.eta)),1]<-1
Upsilon.star[,1]<-Upsilon.star[,1]-hat.alpha*hat.alpha.prim
Upsilon.star[which(eta.star[,2]<=xi2),2]<-1
Upsilon.star[,2]<-Upsilon.star[,2]-hat.alpha.prim
Upsilon.star[,3]<-eta.star[,1]^2-hat.m2
u.eta.star<-u.eta+drop(lambda%*%colMeans(Upsilon.star))
theta1.star<-theta1+(J1.inv/2)%*%colMeans(D1*(eta.star[,1]^2-hat.m2))
if((theta1.star[1]>petit)&(min(theta1.star[2:4])>=0)&(theta1.star[4]<0.99))bon<-1
}
omeg<-max(theta1.star[1],petit); alph1<-max(theta1.star[2],0); alph2<-max(theta1.star[3],0); bet<-max(theta1.star[4],0)
bet<-min(bet,0.99)

for(t in 2:n){
sig2[t]<-omeg+alph1*eps2[t-1,1]+alph2*eps2[t-1,2]+bet*sig2[t-1]
}
sig2.nplus1<-omeg+alph1*eps2[n,1]+alph2*eps2[n,2]+bet*sig2[n]

CoVar.star[b]<--sqrt(sig2.nplus1)*u.eta.star
}
lower.boot<-as.numeric(quantile(CoVar.star,alpha0/2))
upper.boot<-as.numeric(quantile(CoVar.star,1-alpha0/2))
IW.boot<-abs((upper.boot-lower.boot)/CoVaR.nplus1)
Cover.boot<-0
if(true.CoVaR >= lower.boot & true.CoVaR <= upper.boot)Cover.boot<-1 





#####################

#list(omega=omega,mat.alpha=mat.alpha,beta=beta,vol=sqrt(ht),residus=eta,Sigma1=Sigma1,std=std,xi2=xi2,u=u.eta,CoVaR=CoVaR,CoVaR.nplus1=CoVaR.nplus1)
#list(CoVaR.nplus1=CoVaR.nplus1,lower=lower,upper=upper,Cover=Cover,IW=IW,volat=sqrt(ht.nplus1),u.eta=u.eta)
list(Cover=Cover,IW=IW,Cover.norm=Cover.norm,IW.norm=IW.norm,Cover.boot=Cover.boot,IW.boot=IW.boot)
}

estimGARCHbiv<-cmpfun(estimGARCHbiv)

##################
simulation<- function(n, omega,  mat.alpha, beta, rho, u.eta.true, nu) {

res.sim<-BivGarchSt.sim(n, omega,  mat.alpha, beta, rho, u.eta.true, nu) #str(res.sim)
res<-estimGARCHbiv(omega,mat.alpha,beta,res.sim$simu,niv.alpha,niv.alpha.prim,res.sim$true.CoVaR)
c(res$Cover,res$IW,res$Cover.norm,res$IW.norm,res$Cover.boot,res$IW.boot)

}
simulation<-cmpfun(simulation)
#####################################################
#####################################################

omega<-c(0.001,0.001); mat.alpha<-matrix(c(0.05,0.01,0.01,0.1),ncol=2,nrow=2); beta<-c(0.9,0.85); rho<-0.6; nu=6; niv.alpha<-0.1; niv.alpha.prim<-0.2
nu<-6
n<-1000
#omega.init<-c(0.1,01); mat.alpha.init<-matrix(c(0.05,0.05,0.05,0.05),ncol=2,nrow=2); beta.init<-c(0.9,0.9)
Nrep<-2000
set.seed(13)


u.eta.true<-u.Stud(rho,nu, niv.alpha,niv.alpha.prim,-1,-3)$u



res<-replicate(Nrep,simulation(n, omega,  mat.alpha, beta, rho, u.eta.true, nu))



length(which(res[1,]==1))*100/Nrep
mean(res[2,])
length(which(res[3,]==1))*100/Nrep
mean(res[4,])
length(which(res[5,]==1))*100/Nrep
mean(res[6,])


set.seed(13)

n<-3000


res<-replicate(Nrep,simulation(n, omega,  mat.alpha, beta, rho, u.eta.true, nu))



length(which(res[1,]==1))*100/Nrep
mean(res[2,])
length(which(res[3,]==1))*100/Nrep
mean(res[4,])
length(which(res[5,]==1))*100/Nrep
mean(res[6,])



set.seed(13)

n<-5000


res<-replicate(Nrep,simulation(n, omega,  mat.alpha, beta, rho, u.eta.true, nu))



length(which(res[1,]==1))*100/Nrep
mean(res[2,])
length(which(res[3,]==1))*100/Nrep
mean(res[4,])
length(which(res[5,]==1))*100/Nrep
mean(res[6,])



