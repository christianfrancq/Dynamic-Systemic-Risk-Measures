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
ratio<-rep(0,n)
for(i in 21:n) ratio[i]<-mean(DeltaCoVaR[(i-20):i]/VaR[(i-20):i])

list(omega=omega,mat.alpha=mat.alpha,beta=beta,VaR=mean(VaR),DeltaCoVaR=mean(DeltaCoVaR),ratio=ratio)
}

estimGARCHbiv<-cmpfun(estimGARCHbiv)


############################## 
####################################
my_data <- read_xlsx("data_7.2.xlsx")
str(my_data)
time_series <- as.Date(my_data$Date) # Time axis
time_series[1]
first<-1
last<-which(time_series=="2021-10-29")
n<-last-first+1
data_frame<-as.data.frame(my_data[first:last,]) 
str(data_frame)
names<-names(my_data)
series<-matrix(nrow=n,ncol=(length(names)-1))
#series[1:5,]
keep<-c(2:13)
names[keep]
data_frame<-as.data.frame(my_data[first:last,keep]) 
str(data_frame)
series<-matrix(nrow=n,ncol=length(keep))
for(i in 1:length(keep)){
series[,i]<-as.numeric(data_frame[,i])
}

for(i in 1:length(keep)){
if(is.na(series[1,i]))series[1,i]<-series[min(which(!is.na(series[,i]))),i]
if(is.na(series[n,i]))series[n,i]<-series[max(which(!is.na(series[,i]))),i]
}
#series[1:20,]
for(i in 1:length(keep)){
series[,i]<-na.approx(series[,i])
}
#series[1:20,]

data <- data.frame(
  Time = rep(time_series[first:last], times = length(keep)),  # Adjust `times` based on actual number of countries
  Value = as.vector(series),
  Country = rep(names[keep], each = n)
)

# Create the plot
ggplot(data, aes(x = Time, y = Value, color = Country)) +
  geom_line(size = 0.5) +  # Line plot for each country
  labs(title = "",
       x = "Time",
       y = "CDS Spread Value",
       color = "Country") +  # Labels for axes and legend
  scale_y_continuous(limits = c(0, 500)) +  # Set y-axis limits to match the range observed
  theme_minimal() +  # Minimal theme for a clean look
  theme(
    plot.title = element_text(hjust = 0.5),  # Center-align the plot title
    legend.text = element_text(size = 6)     # Adjust legend font size
  )


########### treatment of outliers
#str(series)
# Define the window size
window_size <- 5
series_filtered<-series
for(j in 1:ncol(series)){
moving_avg <- rollmedian(series[,j], k = window_size, fill = NA)
plot.ts(cbind(series[,j],moving_avg))
deviation <- abs(series[,j] - moving_avg)
threshold <- 2 * sd(deviation, na.rm = TRUE)
# Suppress outliers by replacing them with NA or the moving average
series_filtered[,j] <- ifelse(deviation > threshold, moving_avg, series[,j])
series_filtered[1:(floor(window_size/2)),j]<-  series[1:(floor(window_size/2)),j]
series_filtered[(nrow(series)-floor(window_size/2)):nrow(series),j]<-  series[(nrow(series)-floor(window_size/2)):nrow(series),j]
plot.ts(cbind(series[,j],series_filtered[,j]),plot.type="single",col=c("red","blue"))
}


data <- data.frame(
  Time = rep(time_series[first:last], times = length(keep)),  # Adjust `times` based on actual number of countries
  Value = as.vector(series_filtered),
  Country = rep(names[keep], each = n)
)




# Create the plot
p<-ggplot(data, aes(x = Time, y = Value, color = Country)) +
  geom_line(size = 0.5) +  # Line plot for each country
  labs(title = "",
       x = "Time",
       y = "CDS Spread Value",
       color = "Country") +  # Labels for axes and legend
  scale_y_continuous(limits = c(0, 550)) +  # Set y-axis limits to match the range observed
  theme_bw() +  # Minimal theme for a clean look
  theme(
    plot.title = element_text(hjust = 0.5),  # Center-align the plot title
    legend.text = element_text(size = 6)     # Adjust legend font size
  )
ggsave("figure5.pdf", plot = p, width = 9, height = 7)  


