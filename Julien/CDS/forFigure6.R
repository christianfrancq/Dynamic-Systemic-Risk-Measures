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
# estime un AR(1)-GARCH(1,1)  univarié par QMLE 
#########################################################################
#
objf.ar.garch <- function(vartheta, y, n, sig2init, r0=10){     
         mu <- vartheta[1]
         phi <- vartheta[2]
         omeg <- vartheta[3]
         alph <- vartheta[4]
         bet <- vartheta[5]
sig2<-rep(0,n)
eps<-sig2
for(t in 2:n){
eps[t]<-y[t]-mu-phi*y[t-1]
}
sig2[1]<-sig2init
for(t in 2:n){
sig2[t]<-omeg+alph*(eps[t-1]^2)+bet*sig2[t-1]
}
qml <- mean(eps[(r0+1):n]^2/sig2[(r0+1):n]+log(sig2[(r0+1):n]))
qml }

objf.ar.garch<-cmpfun(objf.ar.garch)

#################################
# eps<-rend
# est<-estimGARCHbiv(mu.init,phi.init,omeg.init,alph.init,bet.init,rend, niv.alpha, niv.alpha.prim, niv.alpha.second)


estimARGARCH<- function(mu.init,phi.init,omeg.init,alph.init,bet.init,CDS, niv.alpha, niv.alpha.prim, 
             petit=sqrt(.Machine$double.eps), r0=10)
     {

y<--CDS[,1]
n <- length(y)#plot.ts(y)
valinit<-c(mu.init,phi.init,omeg.init,alph.init,bet.init)
sig2init<-var(y[1:r0])

res <- nlminb(valinit,objf.ar.garch,lower=c(-Inf,-1+petit,petit^2,0,0),
   upper=c(Inf,1-petit,Inf,Inf,1-petit),y=y, n=n,sig2init=sig2init)
par1<-res$par
mu<-res$par[1];  phi<-res$par[2]; omeg <- res$par[3]; alph<- res$par[4]; bet<- res$par[5]

sig2<-rep(0,n)
eps<-sig2
for(t in 2:n){
eps[t]<-y[t]-mu-phi*y[t-1]
}
sig2[1]<-sig2init
for(t in 2:n){
sig2[t]<-omeg+alph*(eps[t-1]^2)+bet*sig2[t-1]
}

resid<-eps/sqrt(pmax(sig2,petit))
eta.1<-resid[(r0+1):n]  #plot.ts(eta.1)

####

y<--CDS[,2]
valinit<-c(mu,phi,omeg,alph,bet)  #valinit<-c(0.2,0.2,0.2,0.01,0.7)

sig2init<-var(y[1:r0])

res <- nlminb(valinit,objf.ar.garch,lower=c(-Inf,-1+petit,petit^2,0,0),
   upper=c(Inf,1-petit,Inf,Inf,1-petit),y=y,n=n,sig2init=sig2init)
par2<-res$par
mu<-res$par[1];  phi<-res$par[2]; omeg <- res$par[3]; alph<- res$par[4]; bet<- res$par[5]

sig2<-rep(0,n)
eps<-sig2
for(t in 2:n){
eps[t]<-y[t]-mu-phi*y[t-1]
}
sig2[1]<-sig2init
for(t in 2:n){
sig2[t]<-omeg+alph*(eps[t-1]^2)+bet*sig2[t-1]
}

resid<-eps/sqrt(pmax(sig2,petit))
eta.2<-resid[(r0+1):n]

#cov(eta.1,eta.2)

n<-length(eta.1)

xi2<-as.numeric(quantile(eta.2,probs=niv.alpha.prim,type=1))
select<-which(eta.2<=xi2)
u.eta<-as.numeric(quantile(eta.1[select],probs=niv.alpha,type=1))



densite<-density(eta.2)
g2.xi2<-densite$y[which.min(abs(densite$x-xi2))]
densite<-density(eta.1[select])
f1.u.xi2<-densite$y[which.min(abs(densite$x-u.eta))]
select<-which(eta.1<=u.eta)
densite<-density(eta.2[select])
f2.xi2.u<-densite$y[which.min(abs(densite$x-xi2))]

G.1.u<-length(which(eta.1<=u.eta))/n
lambda<-c(-1/(niv.alpha.prim*f1.u.xi2),(G.1.u*f2.xi2.u)/(niv.alpha.prim*f1.u.xi2*g2.xi2),-u.eta/2)

Upsilon<-matrix(0,nrow=length(eta.1),ncol=3)
select<-which(eta.1 <= u.eta & eta.2 <= xi2)
Upsilon[select,1]<-1
select<-which( eta.2 <= xi2)
Upsilon[select,2]<-1
Upsilon[,3]<-eta.1^2
Sigma.Upsilon<-var(Upsilon)

s2<-as.numeric(t(lambda)%*%Sigma.Upsilon%*%lambda/n)



list(par1=par1,par2=par2,u.eta=-u.eta,s2=s2)
}

estimARGARCH<-cmpfun(estimARGARCH)
#
########## simule un AR(1)-GARCH(1,1) avec bruit de Student (ou Gaussien)
sim.ar.garch <- function(mu,phi,omeg,alph,bet, n, nu, valinit=500){     
# sim du bruit
norm<-rnorm(n+valinit)
{if(nu==Inf|nu<2){
ct<-1}
else
{ct<-sqrt((nu-2)/nu)
}}
if(nu==Inf){
eta<-norm
}else{
eta<-ct*rt((n+valinit),df=nu)
}
# mean(eta^2)
#
y<-rep(0,(n+valinit))
sig2<-y
eps<-sig2
sig2[1]<-omeg
for(t in 2:(n+valinit)){
sig2[t]<-omeg+alph*(eps[t-1]^2)+bet*sig2[t-1]
eps[t]<-sqrt(sig2[t])*eta[t]
}
for(t in 2:(n+valinit)){
y[t]<-mu+phi*y[t-1]+eps[t]
}
y[(1+valinit):(n+valinit)]
}

sim.ar.garch <-cmpfun(sim.ar.garch )

#mu<-0.1;phi<-0.5;omeg<-0.01;alph<-0.05;bet<-0.9
# nu<-7;n<-50000
#y<-sim.ar.garch(mu,phi,omeg,alph,bet, n, nu)  plot.ts(y)   


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
  theme_minimal() +  # Minimal theme for a clean look
  theme(
    plot.title = element_text(hjust = 0.5),  # Center-align the plot title
    legend.text = element_text(size = 6)     # Adjust legend font size
  )
ggsave("figure5.pdf", plot = p, width = 9, height = 7)  

############### rendements
nr<-nrow(series_filtered)
nc<-ncol(series_filtered)
rend_series_filtered<-matrix(nrow=(nr-1),ncol=nc)
for(j in 1:nc){
rend_series_filtered[,j] <- log(series_filtered[2:nr,j]/series_filtered[1:(nr-1),j])
}
#plot.ts(rend_series_filtered[,1:6]) cor(rend_series_filtered)
#plot.ts(rend_series_filtered[,7:12])

#############################
rend_series_filtered<-as.data.frame(rend_series_filtered)
colnames(rend_series_filtered)<-names[keep]
rend_series_filtered<-rend_series_filtered[,c(1:10,12,11)]
str(rend_series_filtered)

niv.alpha<-0.5 
niv.alpha.prim<-0.5
mu.init<-0;phi.init<-0;omeg.init<-0.01;alph.init<-0.05;bet.init<-0.8
CDS<-matrix(nrow=nrow(rend_series_filtered),ncol=2)
CDS[,1]<-rend_series_filtered[,12]
u.chap<-rep(0,(ncol(rend_series_filtered)-1))
Mini<-u.chap
Maxi<-Mini
Q1<-Mini
Q3<-Q1
for(i in 1:(ncol(rend_series_filtered)-1)){
CDS[,2]<-rend_series_filtered[,i]
res<-estimARGARCH(mu.init,phi.init,omeg.init,alph.init,bet.init,CDS, niv.alpha, niv.alpha.prim)
mu.init<-res$par1[1];phi.init<-res$par1[2];omeg.init<-res$par1[3];alph.init<-res$par1[4];bet.init<-res$par1[5]
u.chap[i]<-res$u.eta
Q1[i]<-res$u.eta+sqrt(res$s2)*qnorm(0.25)
Q3[i]<-res$u.eta+sqrt(res$s2)*qnorm(0.75)
Mini[i]<-Q1[i]-1.5*(Q3[i]-Q1[i])
Maxi[i]<-Q3[i]+1.5*(Q3[i]-Q1[i])
}

############
# Sample data for boxplot values
data <- data.frame(
  Country = colnames(rend_series_filtered)[1:(nc-1)],
  Min = Mini,
  Q1 = Q1,
  Median = u.chap,
  Q3 = Q3,
  Max = Maxi
)

# Reorder Country factor by Median value
data$Country <- factor(data$Country, levels = data$Country[order(data$Median)])

# Create the boxplot with ggplot2
ggplot(data, aes(x = Country, ymin = Min, lower = Q1, middle = Median, upper = Q3, ymax = Max)) +
  geom_boxplot(stat = "identity", fill = "skyblue", color = "black") +
  coord_flip() +  # Horizontal boxplots
  labs(
    title = "CDS Spread Distributions by Country",
    x = "Country",
    y = "CDS Spread"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)  # Center the title
  )

p<-ggplot(data, aes(x = Country, ymin = Min, lower = Q1, middle = Median, upper = Q3, ymax = Max, fill = Median)) +
  geom_boxplot(stat = "identity", color = "black", width = .5) +
  coord_flip() +  # Horizontal boxplots
  scale_fill_gradient(low = "lightblue", high = "darkblue", guide = guide_legend(title = NULL)) +  # Gradient color based on Median
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)  # Center the title
  )
ggsave("figure6.pdf", plot = p, width = 8, height = 5)  


