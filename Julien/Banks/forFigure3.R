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
install.packages("readxl", repos='https://ftp.belnet.be/mirror/CRAN/')
library("readxl")
install.packages("ggrepel", repos='https://ftp.belnet.be/mirror/CRAN/')
library("ggrepel")
install.packages("zoo", repos='https://ftp.belnet.be/mirror/CRAN/')
library("zoo")
install.packages("ggplot2", repos='https://ftp.belnet.be/mirror/CRAN/')
library(ggplot2)

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
#################################
# eps<-rend
# est<-estimGARCHbiv(omega.init,mat.alpha.init,beta.init,rend, niv.alpha, niv.alpha.prim, niv.alpha.second)


estimGARCHbiv<- function(omega.init,mat.alpha.init,beta.init,eps, niv.alpha, niv.alpha.prim, niv.alpha.second, 
             petit=sqrt(.Machine$double.eps), r0=10, alpha0=0.10)
     {
n <- length(eps[,1])#plot.ts(eps)
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
#var(eta)
################ 
xi2<-as.numeric(quantile(eta[,2],probs=niv.alpha.prim,type=1))
VaR<--sqrt(ht[,2])*xi2
select<-which(eta[,2]<=xi2)
u.eta<-as.numeric(quantile(eta[select,1],probs=niv.alpha,type=1))
med.sup<-as.numeric(quantile(eta[,2],probs=0.5+niv.alpha.second,type=1))
med.inf<-as.numeric(quantile(eta[,2],probs=0.5-niv.alpha.second,type=1))

select.median<-which((eta[,2]<=med.sup)&(eta[,2]>med.inf))
u.eta.med<-as.numeric(quantile(eta[select.median,1],probs=niv.alpha,type=1))

DeltaCoVaR<--sqrt(ht[,1])*(u.eta-u.eta.med)
#####################

list(omega=omega,mat.alpha=mat.alpha,beta=beta,VaR=mean(VaR),DeltaCoVaR=mean(DeltaCoVaR))
}

estimGARCHbiv<-cmpfun(estimGARCHbiv)


############################## 
####################################################
####################################
my_data <- read_excel("data_7.1.3.xlsx", sheet = "Prices (local currency)")
str(my_data)
time_series <- as.Date(my_data$Date) # Time axis
first<-which(time_series=="2004-01-02")
last<-which(time_series=="2015-12-01")
n<-last-first+1
names<-names(my_data)
data_frame<-as.data.frame(my_data[first:last,2:19]) 
str(data_frame)
series<-matrix(nrow=n,ncol=18)
for(i in 1:18){
series[,i]<-as.numeric(data_frame[,i])
}
for(i in 1:18){
if(is.na(series[1,i]))series[1,i]<-series[min(which(!is.na(series[,i]))),i]
if(is.na(series[n,i]))series[n,i]<-series[max(which(!is.na(series[,i]))),i]
}
#series[1:20,]
for(i in 1:18){
series[,i]<-na.approx(series[,i])
}
#series[1:20,]

cap <- read_excel("data_7.1.3.xlsx", sheet = "Market Cap (USD)")
str(cap)
time_series_caps <- as.Date(cap$Date) # Time axis
first_caps<-which(time_series=="2004-01-02")
last_caps<-which(time_series=="2015-12-01")
n_caps<-last_caps-first_caps+1

cap_frame<-as.data.frame(cap[first:last,2:19]) 
str(cap_frame)
caps<-matrix(nrow=n,ncol=18)
for(i in 1:18){
caps[,i]<-as.numeric(cap_frame[,i])
}
#caps[1:5,17:18]

for(i in 1:18){
if(is.na(caps[1,i]))caps[1,i] <-caps[min(which(!is.na(caps[ ,i]))),i]
if(is.na(caps[n,i]))caps[n,i] <-caps[max(which(!is.na(caps[ ,i]))),i]
}

for(i in 1:18){
caps[,i]<-na.approx(caps[,i])
}
#caps[1:20,] caps[n,] caps[1,]


returns<-matrix(0,nrow=n,ncol=18)

for(i in 1:18){
returns[2:n,i]<-100*log(series[2:n,i]/series[1:(n-1),i])
#returns[2:n,i]<-100*(series[2:n,i]-series[1:(n-1),i])/series[1:(n-1),i]
}

return_sys<-rep(0,n)
for(i in 1:n){
total_cap<-sum(caps[i,])
for(j in 1:18){
return_sys[i]<-return_sys[i]+returns[i,j]*caps[i,j]/total_cap
}
}

#plot.ts(cbind(return_sys,returns[,1:3]))

date.weights<-which(time_series=="2010-10-01")
total_cap<-sum(caps[date.weights-first+1,])
weights<-rep(0,18)
for(j in 1:18){
weights[j]<-caps[date.weights-first+1,j]/total_cap
}

names
100*weights

correl<-rep(0,18)
for(j in 1:18){
correl[j]<-cor(returns[2:n,j],return_sys[2:n])
}

names
order(correl)
names[order(correl)+1]

str(returns)


omega.init<-c(0.001,0.001); mat.alpha.init<-matrix(c(0.05,0.01,0.01,0.1),ncol=2,nrow=2); beta.init<-c(0.9,0.85) 
niv.alpha<-0.05; niv.alpha.prim<-0.1; niv.alpha.second<-0.3
niv.alpha<-0.05; niv.alpha.prim<-0.1; niv.alpha.second<-0.25
niv.alpha<-0.05; niv.alpha.prim<-0.1; niv.alpha.second<-0.2


VaR<-rep(0,18)
DeltaCoVaR<-rep(0,18)
Size<-VaR
for(i in 1:18){
rend<-cbind(return_sys[2:n],returns[2:n,i])# str(rend)
est<-estimGARCHbiv(omega.init,mat.alpha.init,beta.init,rend, niv.alpha, niv.alpha.prim, niv.alpha.second)
omega.init<-est$omega
mat.alpha.init<-est$mat.alpha
beta.init<-est$beta
VaR[i]<-est$VaR
DeltaCoVaR[i]<-est$DeltaCoVaR
Size[i]<-sum(caps[,i])
}
VaR
DeltaCoVaR
Size

# Sample Data (Replace with actual data)
# Names of the financial institutions
institutions <-names[2:19]

# Average Delta CoVaR (y-axis values)
avg_delta_covar <- DeltaCoVaR

# Institution size (point size in the plot)
sizes <- Size

# Create a DataFrame for plotting
data <- data.frame(
  Institution = institutions,
  VaR = VaR,
  DeltaCoVaR = avg_delta_covar,
  Size = sizes
)

# Create the plot
p<-ggplot(data, aes(x = VaR, y = DeltaCoVaR, size = Size, label = Institution)) +
  geom_point(color = "blue", alpha = 0.6) +  # Scatter plot with point size for institution size
 geom_text_repel(aes(label = Institution), vjust = -1.2, size = 2.5, max.overlaps = 10) +  # Add institution names as labels
  scale_size_continuous(range = c(2, 10)) +  # Adjust the size range of the points
  labs(title = "",
       x = "VaR",
       y = "Delta CoVaR",
       size = "Size") +  # Labels for axes and point sizes
  theme_bw() +  # Use a clean, minimal theme
  theme(plot.title = element_text(hjust = 0.1),
legend.position = "none" 
)  # Center the title


ggsave("figure3.pdf", plot = p, width = 8, height = 7)  

