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


estimGARCHbiv<- function(omega.init,mat.alpha.init,beta.init,eps, 
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

ht[1,]<-c(sig2init1,sig2init2)
eta[1,]<-eps[1,]/sqrt(pmax(ht[1,],petit)) 
for(t in 2:n){
ht[t,]<-omega+mat.alpha%*%eps2[t-1,]+beta*ht[t-1,]
eta[t,]<-eps[t,]/sqrt(pmax(ht[t,],petit))
}
#var(eta)
################ 
DeltaCoVaR<-rep(0,4)
CoVaR<-DeltaCoVaR
vec.alpha<-c(0.5,0.1,0.05,0.01)
vec.alpha.prim<-c(0.5,0.2,0.1,0.05)
vec.alpha.second<-c(0.25,0.25,0.25,0.25)
for(iter in 1:4){
niv.alpha<-vec.alpha[iter]
niv.alpha.prim<-vec.alpha.prim[iter]
niv.alpha.second<-vec.alpha.second[iter]
xi2<-as.numeric(quantile(eta[,2],probs=niv.alpha.prim,type=1))
select<-which(eta[,2]<=xi2)
u.eta<-as.numeric(quantile(eta[select,1],probs=niv.alpha,type=1))
med.sup<-as.numeric(quantile(eta[,2],probs=0.5+niv.alpha.second,type=1))
med.inf<-as.numeric(quantile(eta[,2],probs=0.5-niv.alpha.second,type=1))

select.median<-which((eta[,2]<=med.sup)&(eta[,2]>med.inf))
u.eta.med<-as.numeric(quantile(eta[select.median,1],probs=niv.alpha,type=1))

DeltaCoVaR[iter]<-mean(-sqrt(ht[,1])*(u.eta-u.eta.med))
CoVaR[iter]<-mean(-sqrt(ht[,1])*u.eta)
}
#####################

list(omega=omega,mat.alpha=mat.alpha,beta=beta,CoVaR=mean(CoVaR),DeltaCoVaR=mean(DeltaCoVaR))
}

estimGARCHbiv<-cmpfun(estimGARCHbiv)


############################## 
####################################################
####################################
my_data <- read.csv("CRSP_log_return.csv")
str(my_data)
time_series <- as.Date(my_data$date) # Time axis
first<-which(time_series=="2009-02-05")
last<-which(time_series=="2021-01-05")
n<-last-first+1

US_Banks<-c("S&P Global Inc", "Wells Fargo & Co", "Marsh & McLennan Inc", "JPMorgan Chase & Co", "American Express Co", 
"Bank of America Corp", "Progressive Corp", "Morgan Stanley", "Citigroup Inc","Charles Schwab Corp", "Chubb Ltd",
  "Goldman Sachs Group Inc",    "BlackRock Inc",  "Berkshire Hathaway Inc Class B")
tickers<-c("SPGI",  "WFC",   "MMC",   "JPM",   "AXP",   "BAC",   "PGR",   "MS",    "C",    
"SCHW",  "CB",    "GS",    "BLK",   "BRK.B")



my_ETFdata <- read_xlsx("IXG.xlsx")
str(my_ETFdata)
Sys.setlocale("LC_TIME", "C")
time_series2 <- rev(as.Date(my_ETFdata$Date, format = "%b %d, %Y")) # Time axis
ETF<-rev(my_ETFdata$ETF)
first_ETF<-which(time_series2=="2009-02-05")
last_ETF<-which(time_series2=="2021-01-05")
n_ETF<-last_ETF-first_ETF+1


# time_series2 = time_series
as.Date(setdiff(time_series2[first_ETF:last_ETF],time_series[first:last]))
as.Date(setdiff(time_series[first:last],time_series2[first_ETF:last_ETF]))

data_banks_US<-as.data.frame(my_data[first:last,]) 
str(data_banks_US)

data_banks_US$date<-time_series2[first_ETF:last_ETF]
ETF<-ETF[(first_ETF-1):last_ETF]
rend.ETF<-log(ETF[2:length(ETF)]/ETF[1:(length(ETF)-1)])
data_banks_US$ETF<-rend.ETF


CoVaR<-rep(0,14)
DeltaCoVaR<-CoVaR
rho<-CoVaR


for(i in 1:14){
omega.init<-c(0.01,0.01); mat.alpha.init<-matrix(c(0.05,0.01,0.01,0.1),ncol=2,nrow=2); beta.init<-c(0.9,0.85) 
rend<-100*cbind(data_banks_US$ETF,data_banks_US[,(i+1)])# str(rend)
est<-estimGARCHbiv(omega.init,mat.alpha.init,beta.init,rend)
omega.init<-est$omega
mat.alpha.init<-est$mat.alpha
beta.init<-est$beta
CoVaR[i]<-est$CoVaR
DeltaCoVaR[i]<-est$DeltaCoVaR
rho[i]<-cor(data_banks_US$ETF,data_banks_US[,(i+1)])
}
order_CoVaR<-order(CoVaR,decreasing = TRUE)
order_DeltaCoVaR<-order(DeltaCoVaR,decreasing = TRUE)
DeltaCoVaR[order_DeltaCoVaR]
tickers[order_DeltaCoVaR]
rho


ranked_data_banks_US<-data_banks_US
ranked_tickers<-tickers[order_DeltaCoVaR]

names.rank<-US_Banks
for(i in 1:14){
names.rank[i]<-US_Banks[order_DeltaCoVaR[i]]
ranked_data_banks_US[,(i+1)]<-data_banks_US[,(order_DeltaCoVaR[i]+1)]
}
colnames(ranked_data_banks_US)[2:15]<-ranked_tickers
str(data_banks_US)
str(ranked_data_banks_US)


for(i in 1:14){
rend<-100*cbind(ranked_data_banks_US$ETF,ranked_data_banks_US[,(i+1)])# str(rend)
omega.init<-c(0.01,0.01); mat.alpha.init<-matrix(c(0.05,0.01,0.01,0.1),ncol=2,nrow=2); beta.init<-c(0.9,0.85) 
est<-estimGARCHbiv(omega.init,mat.alpha.init,beta.init,rend)
omega.init<-est$omega
mat.alpha.init<-est$mat.alpha
beta.init<-est$beta
CoVaR[i]<-est$CoVaR
DeltaCoVaR[i]<-est$DeltaCoVaR
rho[i]<-cor(ranked_data_banks_US$ETF,ranked_data_banks_US[,(i+1)])
}
order_CoVaR<-order(CoVaR,decreasing = TRUE)
order_DeltaCoVaR<-order(DeltaCoVaR,decreasing = TRUE)
order(rho,decreasing = TRUE)

sort(DeltaCoVaR,decreasing = TRUE)
sort(rho)

### 
data <- data.frame(
  `Company Name` = names.rank,
  `Rank wrt. rho` = order(rho[order_DeltaCoVaR],decreasing = TRUE),
  `Rank wrt. CoVaR` = order(CoVaR[order_DeltaCoVaR],decreasing = TRUE),
  `Rank wrt. Delta CoVaR` =order(DeltaCoVaR[order_DeltaCoVaR],decreasing = TRUE)
)


# Create the xtable object
table_xtable <- xtable(
  data,
  caption = "Ranking of Financial Institutions by Systemic Impact (i.e., $\\overline{\\Delta CoVaR}$)"
)

# Customize the column names
colnames(table_xtable) <- c("Company Name", 
                            "Rank wrt. $\\rho$", "Rank wrt. $\\overline{CoVaR}$", "Rank wrt. $\\overline{\\Delta CoVaR}$")

# Print the LaTeX code for the table with formatting options
print(table_xtable,
      include.rownames = FALSE,          # Remove row numbers
      sanitize.text.function = identity, # Allows LaTeX math expressions (e.g., $\\rho$)
      tabular.environment = "tabular",   # Sets the table environment
      floating = TRUE,                   # Allows the table to float
      booktabs = TRUE,                   # Enables cleaner table lines
      hline.after = c(-1, 0, nrow(data)) # Adds lines at specific rows
)


####################################################################################################
############################## 
####################################################
####################################
my_data <- read.csv("19_log_return.csv")
str(my_data)
time_series <- as.Date(my_data$X) # Time axis
first<-which(time_series=="2009-02-05")
last<-which(time_series=="2021-01-05")
n<-last-first+1

data_banks_US_ETF<-data_banks_US
str(data_banks_US)
data_banks_US$ETF
# All the elements in data_banks$date are in time_series
as.Date(setdiff(data_banks_US_ETF$date,time_series[first:last]))
# reverse not true
as.Date(setdiff(time_series[first:last],data_banks_US_ETF$date))
tickers_sup<-c("RY.TO",  "ALV.DE",   "TD.TO",   "HSBA.L", "CBA.AX")

Bank_names<-c(US_Banks,"Royal Bank of Canada", "Allianz", "Toronto Dominion Bank",  
 "HSBC Holdings PLC", "Commonwealth Bank of Australia")



b_sup<-matrix(nrow=length(data_banks_US_ETF$date),ncol=length(tickers_sup))
for(i in 1:length(data_banks_US_ETF$date)){
for(j in 1:length(tickers_sup)){
b_sup[i,j]<-my_data[which(time_series[first:last]==data_banks_US_ETF$date[i]),tickers_sup[j]]
}}

str(b_sup)

for(i in 1:length(tickers_sup)){
if(is.na(b_sup[1,i]))b_sup[1,i]<-b_sup[min(which(!is.na(b_sup[,i]))),i]
if(is.na(b_sup[nrow(b_sup),i]))b_sup[nrow(b_sup),i]<-b_sup[max(which(!is.na(b_sup[,i]))),i]
}


for(i in 1:length(tickers_sup)){
b_sup[,i]<-na.approx(b_sup[,i])
}

str(b_sup)

b_sup<-as.data.frame(b_sup)

colnames(b_sup)<-tickers_sup
str(b_sup)

data_banks<-cbind(data_banks_US_ETF[,1:15],b_sup)
data_banks$ETF<-data_banks_US_ETF$ETF

str(data_banks)
str(data_banks_US)


CoVaR<-rep(0,19)
DeltaCoVaR<-CoVaR
rho<-CoVaR


for(i in 1:19){
rend<-100*cbind(data_banks$ETF,data_banks[,(i+1)])# str(rend)
omega.init<-c(0.01,0.01); mat.alpha.init<-matrix(c(0.05,0.01,0.01,0.1),ncol=2,nrow=2); beta.init<-c(0.9,0.85) 
est<-estimGARCHbiv(omega.init,mat.alpha.init,beta.init,rend)
omega.init<-est$omega
mat.alpha.init<-est$mat.alpha
beta.init<-est$beta
CoVaR[i]<-est$CoVaR
DeltaCoVaR[i]<-est$DeltaCoVaR
rho[i]<-cor(data_banks$ETF,data_banks[,(i+1)])
}

order_CoVaR<-order(CoVaR,decreasing = TRUE)
order_DeltaCoVaR<-order(DeltaCoVaR,decreasing = TRUE)
names.rank<-Bank_names[order_DeltaCoVaR]

sort(DeltaCoVaR,decreasing = TRUE)


### 
data <- data.frame(
  `Company Name` = names.rank,
  `Rank wrt. rho` = order(rho[order_DeltaCoVaR],decreasing = TRUE),
  `Rank wrt. CoVaR` = order(CoVaR[order_DeltaCoVaR],decreasing = TRUE),
  `Rank wrt. Delta CoVaR` =order(DeltaCoVaR[order_DeltaCoVaR],decreasing = TRUE)
)


# Create the xtable object
table_xtable <- xtable(
  data,
  caption = "Ranking of Financial Institutions by Systemic Impact (i.e., $\\overline{\\Delta CoVaR}$)"
)

# Customize the column names
colnames(table_xtable) <- c("Company Name", 
                            "Rank wrt. $\\rho$", "Rank wrt. $\\overline{CoVaR}$", "Rank wrt. $\\overline{\\Delta CoVaR}$")

# Print the LaTeX code for the table with formatting options
print(table_xtable,
      include.rownames = FALSE,          # Remove row numbers
      sanitize.text.function = identity, # Allows LaTeX math expressions (e.g., $\\rho$)
      tabular.environment = "tabular",   # Sets the table environment
      floating = TRUE,                   # Allows the table to float
      booktabs = TRUE,                   # Enables cleaner table lines
      hline.after = c(-1, 0, nrow(data)) # Adds lines at specific rows
)










