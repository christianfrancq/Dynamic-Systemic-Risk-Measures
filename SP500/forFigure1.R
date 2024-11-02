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
#################################
# eps<-rend

estimGARCHbiv<- function(P, omega.init,mat.alpha.init,beta.init,eps, niv.alpha, niv.alpha.prim, 
             petit=sqrt(.Machine$double.eps), r0=10, alpha0=0.10, B=500, nb_try=5)
     {
n <- length(eps[,1])
omega<-c(0,0)
beta<-omega
mat.alpha<-matrix(nrow=2,ncol=2)
valinit<-c(omega.init[1],mat.alpha.init[1,1:2],beta.init[1])
eps2<-eps^2
sig2init1<-mean(eps2[1:r0,1])

res <- nlminb(valinit,objf.garch.biv1,lower=c(petit^2,0,0,0),
   upper=c(Inf,Inf,Inf,1),eps2=eps2,n=n-P,sig2init=sig2init1)
         omega[1] <- res$par[1]
         mat.alpha[1,1:2]<- res$par[2:3]
         beta[1] <- res$par[4]

valinit<-c(omega.init[2],mat.alpha.init[2,1:2],beta.init[2])
sig2init2<-mean(eps2[1:r0,2])

res <- nlminb(valinit,objf.garch.biv2,lower=c(petit^2,0,0,0),
   upper=c(Inf,Inf,Inf,1),eps2=eps2,n=n-P,sig2init=sig2init2)
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

#var(eta[1:(n-P),])
s<-sqrt(diag(var(eta[1:(n-P),]))) 
for(t in 1:(n-P)){
eta[t,1]<-(eta[t,1]-mean(eta[1:(n-P),1]))/s[1]
eta[t,2]<-(eta[t,2]-mean(eta[1:(n-P),2]))/s[2]
}
#var(eta) plot.ts(eta)
################ semi-parametric method
#colMeans(eta[1:(n-P),]^2)
D1<-der.sig2.1[1:(n-P),]/(2*ht[1:(n-P),1])
J1<-t(D1)%*%D1/n
Omega1<-colMeans(D1)
J1.inv<-Inv(J1)
xi2<-as.numeric(quantile(eta[1:(n-P),2],probs=niv.alpha.prim,type=1))
select<-which(eta[1:(n-P),2]<=xi2)
u.eta<-as.numeric(quantile(eta[select,1],probs=niv.alpha,type=1))
CoVaR<--sqrt(ht[,1])*u.eta
Upsilon<-matrix(0,nrow=(n-P),ncol=3)
select<-which(eta[1:(n-P),1] <= u.eta & eta[1:(n-P),2] <= xi2)
Upsilon[select,1]<-1
select<-which( eta[1:(n-P),2] <= xi2)
Upsilon[select,2]<-1
Upsilon[1:(n-P),3]<-eta[1:(n-P),1]^2
Sigma.Upsilon<-var(Upsilon)

densite<-density(eta[1:(n-P),2])
g2.xi2<-densite$y[which.min(abs(densite$x-xi2))]
densite<-density(eta[select,1])
f1.u.xi2<-densite$y[which.min(abs(densite$x-u.eta))]
select<-which(eta[1:(n-P),1]<=u.eta)
densite<-density(eta[select,2])
f2.xi2.u<-densite$y[which.min(abs(densite$x-xi2))]
G.1.u<-length(which(eta[1:(n-P),1]<=u.eta))/n
lambda<-c(-1/(niv.alpha.prim*f1.u.xi2),(G.1.u*f2.xi2.u)/(niv.alpha.prim*f1.u.xi2*g2.xi2),-u.eta/2)
bigUpsilon<-matrix(0,nrow=(n-P),ncol=5)
bigUpsilon[1:(n-P),1:4]<-(D1%*%J1.inv*(eta[1:(n-P),1]^2-1))/2 
bigUpsilon[1:(n-P),5]<-Upsilon%*%lambda

Sigma.alpha.alpha.prim<-var(bigUpsilon)

delta<-rep(0,5)
lower.CoVaR<-CoVaR
upper.CoVaR<-CoVaR

for(t in 2:n){


delta[1:4]<-u.eta*(c(1,eps2[t-1,1],eps2[t-1,2],ht[t-1,1])+beta[1]*der.sig2.1[t-1,])/(sqrt(ht[t-1,1])*2)
delta[5]<-sqrt(ht[t-1,1])


sd<-as.numeric(sqrt(delta%*%Sigma.alpha.alpha.prim%*%delta/n))
lower.CoVaR[t]<-CoVaR[t]-qnorm(1-alpha0/2)*sd
upper.CoVaR[t]<-CoVaR[t]+qnorm(1-alpha0/2)*sd
}


################ Gaussian method
rho<-cor(eta[1:(n-P),])[1,2]
u.eta.norm<-u.Norm(rho,niv.alpha,niv.alpha.prim,-2.5)$u
CoVaR.norm<--sqrt(ht[,1])*u.eta.norm
xi.norm<-qnorm(niv.alpha.prim)
z.rho<-(xi.norm-rho*u.eta.norm)/sqrt(1-rho^2)
K.rho<-dnorm(z.rho)/pnorm(z.rho)
Sigma.N<-matrix(nrow=5,ncol=5)
Sigma.N[1:4,1:4]<-J1.inv/2
Sigma.N[1:4,5]<-(-rho*sqrt(1-rho^2)*K.rho/2)*J1.inv%*%Omega1
Sigma.N[5,1:4]<-Sigma.N[1:4,5]
Sigma.N[5,5]<-K.rho^2*(1-rho^2)
#eigen(Sigma.N)
lower.CoVaR.norm<-CoVaR.norm
upper.CoVaR.norm<-CoVaR.norm


for(t in 2:n){


delta[1:4]<-u.eta.norm*(c(1,eps2[t-1,1],eps2[t-1,2],ht[t-1,1])+beta[1]*der.sig2.1[t-1,])/(sqrt(ht[t-1,1])*2)
delta[5]<-sqrt(ht[t-1,1])


sd<-as.numeric(sqrt(delta%*%Sigma.N%*%delta/n))
lower.CoVaR.norm[t]<-CoVaR.norm[t]-qnorm(1-alpha0/2)*sd
upper.CoVaR.norm[t]<-CoVaR.norm[t]+qnorm(1-alpha0/2)*sd
}


########################## Bootstrap method
hat.alpha.prim<-length(which(eta[1:(n-P),2]<=xi2))/(n-P)
hat.alpha<-(length(which((eta[1:(n-P),2]<=xi2)&(eta[1:(n-P),1]<=u.eta)))/(n-P))/hat.alpha.prim

hat.m2<-mean(eta[1:(n-P),1]^2)
CoVaR.star<-matrix(nrow=n,ncol=B)
lower.boot<-CoVaR
upper.boot<-CoVaR
sig2<-rep(0,n)
sig2[1]<-sig2init1
theta1<-c(omega[1],mat.alpha[1,1:2],beta[1])
for(b in 1:B){
bon<-0
try<-0
while(bon==0){
try<-try+1
if(try>nb_try)bon<-1
indic<-sample(c(1:(n-P)),replace=TRUE)
eta.star<-eta[indic,]
#var(eta.star)
Upsilon.star<-matrix(0,nrow=(n-P),ncol=3)
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

CoVaR.star[,b]<--sqrt(sig2)*u.eta.star
}


for(t in 2:n){

lower.boot[t]<-as.numeric(quantile(CoVaR.star[t,],alpha0/2))
upper.boot[t]<-as.numeric(quantile(CoVaR.star[t,],1-alpha0/2))
}





#####################

list(lower.CoVaR=lower.CoVaR,CoVaR=CoVaR,upper.CoVaR=upper.CoVaR,lower.CoVaR.norm=lower.CoVaR.norm,CoVaR.norm=CoVaR.norm,upper.CoVaR.norm=upper.CoVaR.norm,lower.boot=lower.boot,upper.boot=upper.boot)
}

estimGARCHbiv<-cmpfun(estimGARCHbiv)


############################## 
####################################################

omega.init<-c(0.001,0.001); mat.alpha.init<-matrix(c(0.05,0.01,0.01,0.1),ncol=2,nrow=2); beta.init<-c(0.9,0.85) 

####################################
serie<-read.csv(file="AAPL_SP500.csv") #str(serie) 
last.for.est<-which(serie$date=="2020-12-31")
first<-which(serie$date=="2021-01-04")
last<-which(serie$date=="2021-12-31")
series<-drop(as.matrix(serie[1:last,2:3])) #str(series) series[1:5,]
n<-length(series[,1])
rend<-matrix(nrow=(n-1),ncol=2)
for(i in 1:(n-1)){
#rend[i,1]<-series[i+1,1]/series[i,1]-1
#rend[i,2]<-series[i+1,2]/series[i,2]-1
rend[i,1]<-log(series[i+1,1]/series[i,1])
rend[i,2]<-log(series[i+1,2]/series[i,2])
}

P<-last-first+1
niv.alpha<-0.05; niv.alpha.prim<-0.1
niv.alpha<-0.1; niv.alpha.prim<-0.2

est<-estimGARCHbiv(P, omega.init,mat.alpha.init,beta.init,rend, niv.alpha, niv.alpha.prim, 
             petit=sqrt(.Machine$double.eps), r0=10, alpha0=0.05, B=100)
#str(est)
# Required Libraries
library(ggplot2)

# Sample Data (Replace with actual data)
# Assuming you have the log-returns and VaR values in two separate vectors
time_series <- as.Date(serie$date[first:last]) # Time axis
log_returns <- rend[(first-1):(last-1),1]  # Replace with actual log-returns
# Define confidence intervals (replace with your calculated upper and lower bounds)
lower_bound <- est$lower.CoVaR.norm[(first-1):(last-1)]  # Example of a lower bound (replace with real data)
upper_bound <- est$upper.CoVaR.norm[(first-1):(last-1)]  # Example of an upper bound (replace with real data)
estimated_covar <- est$CoVaR.norm[(first-1):(last-1)] # Replace with actual CoVaR values

# Define confidence intervals (replace with your calculated upper and lower bounds)
lower_bound <- est$lower.CoVaR[(first-1):(last-1)]  # Example of a lower bound (replace with real data)
upper_bound <- est$upper.CoVaR[(first-1):(last-1)]  # Example of an upper bound (replace with real data)



# Define confidence intervals (replace with your calculated upper and lower bounds)
estimated_covar <- est$CoVaR[(first-1):(last-1)] # Replace with actual CoVaR values
lower_bound <- est$lower.boot[(first-1):(last-1)]  # Example of a lower bound (replace with real data)
upper_bound <- est$upper.boot[(first-1):(last-1)]  # Example of an upper bound (replace with real data)



# Create Data Frame for ggplot
data <- data.frame(
  Date = time_series,
  Log_Returns = log_returns,
  Estimated_CoVaR = -estimated_covar,
  Lower_Bound = -upper_bound,
  Upper_Bound = -lower_bound
)

# Plot
p<-ggplot(data, aes(x = Date)) +
  geom_line(aes(y = Log_Returns, color = "AAPL Log-Returns"), size = 0.5) + 
  geom_line(aes(y = Estimated_CoVaR, color = "-Estimated CoVaR"), size = 0.5) +
  geom_ribbon(aes(ymin = Lower_Bound, ymax = Upper_Bound, fill = "Confidence Interval"), alpha = 0.2) +
  labs(title = "",
       x = "Time (Backtest Period)",
       y = "Value",
       color = "Variables",
       fill = "") +
  scale_color_manual(values = c("AAPL Log-Returns" = "black",  "-Estimated CoVaR" = "blue")) +
  theme_bw() +
  theme(legend.text = element_text(size = 6))


ggsave("figure1.pdf", plot = p, width = 8, height = 4)  

