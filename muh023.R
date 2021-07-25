set.seed(123)
rm(list=ls())
library(mgcv)
library(Matrix)
library(here)
library(dplyr)
library(tidyverse)
library(quantreg)
library(splines)



# In this Part, I am classifying cholestrol into different level. Then I build data ready to use by gam funnction for smoothing
r1_IFG2 <- read.csv("C:/Users/Muhannad/Desktop/multvariat/R_codes/prediabetes/pre_sub/Data/ldl2.csv",header=TRUE, sep=",")
#r1_IFG2<- read.csv(here::here("data", "ldl2.csv"),header=TRUE, sep=",")
r1_IFG=r1_IFG2 
#remove duplicated rows and remove na data.
#r1_IFG=r1_IFG[!duplicated(r1_IFG[ , c("SEQN")]),]
#r1_IFG$HDL=r1_IFG$LBXTC-r1_IFG$LBDLDL

#r1_IFG[is.na(r1_IFG$HDL),] <- 0
#r1_IFG$LBDHDD[is.na(r1_IFG$LBDHDD)] <- 0
#
#r1_IFG$LBXTC[is.na(r1_IFG$LBXTC)] <- 0
#r1_IFG$LBDLDL[is.na(r1_IFG$LBDLDL)] <- 0
#r1_IFG$HDL=r1_IFG$LBDHDD+r1_IFG$HDL
#r1_IFG$LBXTC=r1_IFG$LBDLDL+r1_IFG$HDL
#
#r1_IFG=r1_IFG[r1_IFG$HDL!=0,]
#r1_IFG=r1_IFG[r1_IFG$LBXTC!=0,]
r1_IFG$LR=r1_IFG$LBXTC/r1_IFG$LBDHDD
summary(r1_IFG$LR)
###
#to remove rows that have zero values
#r1_IFG$[r1_IFG$ != 0, ]

# To consider BMWXT. below 18.5 ? you're in the underweight range
#between 18.5 and 24.9 ? you're in the healthy weight range
#between 25 and 29.9 ? you're in the overweight range
#between 30 and 39.9 ? you're in the obese range 
#r1_IFG=r1_IFG[r1_IFG$BMXBMI<29.9,]
#& r1_IFG$BMXBMI<29.9

#r1_IFG=t
#for medicated people
#r1_IFG=r1_IFG[(r1_IFG$BPQ090D==1 & r1_IFG$BPQ100D==1),]
# Not medicated people
#r1_IFG=r1_IFG[(r1_IFG$BPQ090D==2 | r1_IFG$BPQ100D==2),]

#r1_IFG=r1_IFG[,c(3,5:9)]

r1_IFG=r1_IFG[r1_IFG$RIDAGEYR>39 & r1_IFG$RIDAGEYR<60,]
#   Poeple who took medication TCRx=1, otherwise TCRx=0,36->MEC weight
r1_IFG=r1_IFG[r1_IFG$TCRx==0,]
# Statin names :atorvastatin (Lipitor),fluvastatin (Lescol, Lescol XL),lovastatin (Mevacor, Altoprev),
#pravastatin (Pravachol),rosuvastatin (Crestor),simvastatin (Zocor),pitavastatin (Livalo).
#Statin Users code :  LOVASTATIN   d00280,SIMVASTATIN d00746,ROSUVASTATIN d04851, ATORVASTATIN d04105,PRAVASTATIN d00348, EZIT\SIMVASTATIN d05348, AMILop\ATROV d05048,

#r1_IFG=r1_IFG[r1_IFG$RXDDRGID=='d00280'|r1_IFG$RXDDRGID=='d00746'|r1_IFG$RXDDRGID=='d04851'|r1_IFG$RXDDRGID=='d04105'|r1_IFG$RXDDRGID=='d00348'|r1_IFG$RXDDRGID=='d05348'|r1_IFG$RXDDRGID=='d05048',]  
########
# RIDRETH1 (race, 1=mexican, 2=otherHispanic, 3= Non-Hispanic White,4=Non-Hispanic Black ,5=Other Race - Including Multi-Racial

#r1_IFG=r1_IFG[r1_IFG$RIDRETH1==3,]
#for male,RIAGENDR=1, for female, RIAGENDR=2
r1_IFG=r1_IFG[r1_IFG$RIAGENDR==1,]

r1_IFG=r1_IFG[r1_IFG$LBXGLU>99.99 & r1_IFG$LBXGLU<126,]
#r1_IFG=r1_IFG[r1_IFG$LBXGLU<125.99,]
#32-> total chol,33->TCRx,  34->LDL

#To use gender
#r1_IFG=r1_IFG[,c(5,6,22,31,32,37)]
#r1_IFG=r1_IFG[,c(6,22,31,32,37,42)]
head(r1_IFG)
r1_IFG=r1_IFG%>% select (RIDAGEYR,LR)

r1_IFG=unique(r1_IFG)
r1_IFG=na.omit(r1_IFG)

#r1_IFG$WTMEC6YR=((r1_IFG$WTMEC6YR-min(r1_IFG$WTMEC6YR))/(max(r1_IFG$WTMEC6YR)-min(r1_IFG$WTMEC6YR)))*100
# r1_IFG <- r1_IFG[rep(row.names(r1_IFG), round(r1_IFG$WTMEC6YR)), 1:2]


#r1_IFG <- r1_IFG[rep(row.names(r1_IFG), round(r1_IFG$WTMEC6YR)), 1:2]

# cuts LDL: Good LDL<100, border LDL<140, High LDL>160 
#       TC: Good<100, Border TC<240, High TC > 240, Statin effects 0.19815 (KURT)
#1,3.047,4.74703, 6.0139
temp <- cut(r1_IFG$LR, breaks = c(1,3.84,5.5,8.7),
            #temp <- cut(r1_IFG$LR, breaks = c(1,4.552,6.589,10.42),
            # temp <- cut(r1_IFG$LR, breaks = c(1,3.8,5.3,9),
            labels = c(0,1,2),
            right = FALSE)
data_male=r1_IFG

temp=as.numeric(as.character(temp))
data_male=cbind(r1_IFG,temp)

data_male=na.omit(data_male)
K=3
head(data_male)


##############
### Creating factor for age
#Age <-seq(1,length = nrow(data_male))
#f<- function(x) { (x%/%10) }
# temp <- lapply(Age , function(x) f(data_male$RIDAGEYR[x]) )
# temp <- do.call(rbind, temp)
#Age <- as.factor(temp)

#To use weight uncome
#bmi_1=data_male[,c(1,3)]
bmi_1=data_male[,c(1,3)]
xxyy=table(bmi_1)
xxy=as.data.frame.matrix(xxyy)


##############

# Fitting MSSMs to cross-sectional data from the Netherlands: data

# Ardo, UCL 2016
library(msm)

# Prelim:
digits<-3

# (Extended) Data taken from van de Kasteele et al (SIM)
# (Only subset so that mortality is no problem):
shift <- -40
age=as.numeric(rownames(xxy))+ shift
#age   <- 30:70 

# age=40:80
n=length(age)
xxy=as.data.frame.matrix(table(bmi_1))
#
freq1 <-xxy$`0`  
freq2 <- xxy$`1`
freq3 <- xxy$`2`
size  <- freq1+freq2+freq3
dta   <- cbind(age=age,freq1=freq1,freq2=freq2,freq3=freq3,size=size)
n     <- length(age)
####################


# Print data:
cat("\nAge transformed by age + shift with shift = ",shift,"\n\n")
cat("Data set:\n")
print(dta)

# Plot data:
# Prelim:
lwd <- 4
pch <- 16
col <- c(1,2,3)
# Plot framwork:

plot(c(age[1]-shift,age[n]-shift),c(0,1),type="n",xlab="Age",
     ylab="Distribution 3 States")
# Plot lines:
for(i in 1:3){
  y <- dta[,i+1]/size
  points(age-shift,y,pch=16,col=col[i])
  lines(age-shift,y,lwd=lwd,col=col[i])
}
# Add legend:
legend(36, 1, c("State 1", "State 2", "State 3"), col = col,
       text.col = "black", lty = 1, cex = 1.5, lwd=lwd ,
       merge = TRUE, bg = "gray")

# Fitting MSSMs to cross-sectional data from the Netherlands: models

# Ardo, UCL 2016


# Prelim:
digits <- 2

# Model:
Parameters <- function(p){
  beta <- c(p[1],p[2],p[1],p[2])
  xi   <- c(p[3],0,p[3],0)
  list(beta=beta,xi=xi)
}
pnames <- c("q12","q21","q23","q32","xiF")
p0     <- c(-1,-2,0,-1,0)

# Q matrix:
Qmatrix <- function(beta,xi,t){
  q12 <- exp(beta[1]+xi[1]*t)
  q21 <- exp(beta[2]+xi[2]*t)
  q23 <- exp(beta[3]+xi[1]*t)
  q32 <- exp(beta[4]+xi[2]*t)
  matrix(c(-q12,q12,0, q21,-(q21+q23),q23, 0,q32,-q32),3,3,byrow=TRUE)
}

# Loglikelihood:      
loglikelihood <- function(p){
  # Extract parameters:
  param <- Parameters(p)
  beta  <- param$beta
  xi    <- param$xi
  # Contribution per year:
  loglik <- 0
  for(i in 2:n){
    # Prevalence:
    pie <- c(freq1[i-1],freq2[i-1],freq3[i-1])/size[i-1]
    x   <- c(freq1[i],freq2[i],freq3[i])
    # Qmatrix:
    Q <- Qmatrix(beta,xi,age[i-1])
    # One-step probs:
    P <- MatrixExp(Q,t=1)
    # Update:
    loglik <- loglik+log(dmultinom(x=x,prob=pie%*%P))
  }
  # Return:
  -loglik
}

# Minimise:
max <- optim(par=p0, fn=loglikelihood, method = "Nelder-Mead",
             control=list(maxit=5000,trace=FALSE),hessian=TRUE)
cat("\nStarting values =", p0,"\n")
cat("-2log(likelihood) =", round(2*max$value,digits),"\n")
conv <- max$convergence
cat("Convergence code =", conv,"\n")
p    <- max$par 
p.se <- sqrt(diag(solve(max$hessian)))
print(cbind(parameter=pnames,p=round(p,digits),se=round(p.se,digits)),
      quote=FALSE)

# Extract parameters:
param <- Parameters(p)
beta  <- param$beta
xi    <- param$xi
# Construct estimated one-year P matrices:
P <- array(NA,c(n,3,3))
for(i in 1:n){
  # Qmatrix:
  Q <- Qmatrix(beta,xi,age[i])
  # One-step probs:
  P[i,,] <- MatrixExp(Q,t=1)
}

# Prediction:
pie <- c(freq1[1],freq2[1],freq3[1])/size[1]
pie.pred     <-matrix(NA,n,3)
pie.pred[1,] <-pie
P.pred       <-array(NA,c(n,3,3))
P.pred[1,,] <-diag(3)
for(i in 2:n){
  # Qmatrix:
  Q <- Qmatrix(beta,xi,age[i-1])
  # One-step probs:
  P.pred[i,,] <- P.pred[i-1,,]%*%MatrixExp(Q,t=1)
  # Prediction:
  pie.pred[i,] <- pie%*%P.pred[i,,]
}

# Plot prediction?:
png(file = here::here("images", "prevmns01.png"),
    res = 400, height = 9, width = 16, units = "in")

cex.lab <- 1.2
plot(c(age[1]-shift,age[n]-shift),c(0,1),type="n",xlab="Age",
     ylab="Prevalence",cex.lab=cex.lab, main="")
pch <- c(15,17,19)
greys <- gray.colors(n=3, start = 0.2, end = 0.7, gamma = 2.2)
for(i in 1:3){
  lines(age-shift,pie.pred[,i],lwd=lwd,col=greys[i],lty=1)
  y <- dta[,i+1]/size
  points(age-shift,y,pch=pch[i],col=greys[i])
}
legend(60, 1, c("LOw  ","Middle  ","High"), col = greys,
       text.col = 1, lty = NULL, pch = pch, cex=1, bg = "white",
       text.width = strwidth("000000"))


dev.off()
age.i <-  age
age.pred=age-shift

n.pred <- length(age.pred)
# Example P-matrix:
# age.i <- c(25,30,35,40,45,50,55,60,65,70,75) + shift
for(i in 1:length(age.i)){
  cat("\nFor age", age.i[i] -shift,"one-year P-matrix is: \n")
  Q <- Qmatrix(beta,xi,age.i[i])
  print(round(MatrixExp(Q,t=1),digits))
}

####################

K=3

trans      <- list(E.trans = array(NA, dim = c(K, K, n.pred - 1)))

#for (i in 1:(n.pred - 1)) {

for(i in 1:length(age.i)){
  cat("\nFor age", age.i[i]-shift,"one-year P-matrix is: \n")
  Q <- Qmatrix(beta,xi,age.i[i])
  print(round(MatrixExp(Q,t=1),digits))
  trans$E.trans[, , i]<-round(MatrixExp(Q,t=1),digits)
}


##############


###########

png(file = here::here("images", "tranmns01.png"),
    res = 400, height = 8, width = 10, units = "in")
K=3

names.states <- c("LR1", "LR2", "LR3" )
# Plot transition probabilities in K x K matrix
#

# Set graphical parameters

par(
  mar = c(3, 3, 1.5, 0.1),
  mfrow = c(K, K),
  mgp = c(1.8, 0.7, 0))

# 
for (i in 1:K) {
  for (j in 1:3) {
    plot.new()
    with(trans, {
      # If simulations are not available
      # if (!exists.sim) {
      # Only plot the expected values
      plot.window(
        xlim = range(age.pred),
        ylim = range(E.trans[i, j, ]))
      lines(x = age.pred[-1], y = E.trans[i, j, ], lty = 1, col = 1)
      
      
      # } 
      #else {
      # Plot expected values and confidence intervals
      #  plot.window(
      #    xlim = range(age.pred),
      #   ylim = range(c(l.trans[i, j, ], u.trans[i, j, ])))
      # polygon(
      #   x = c(age.pred[-n.pred], rev(age.pred[-1])),
      #   y = c(l.trans[i, j, ], rev(u.trans[i, j, ])),
      #   col = grey(0.95), border = NA)
      # lines(x = age.pred[-1], y = l.trans[i, j, ], lty = 3, col = grey(0.5))
      # lines(x = age.pred[-1], y = u.trans[i, j, ], lty = 3, col = grey(0.5))
      # lines(x = age.pred[-1], y = E.trans[i, j, ], lty = 1, col = 1)
      # }
      # Axes and titles
      axis(1)
      axis(2)
      box( )
      title(
        main = paste(names.states[i], "to", names.states[j]),
        xlab = "Age",
        ylab = "Transition")
    })
  }
}

dev.off()
#######################


rm(list=ls())
# In this Part, I am classifying cholestrol into different level. Then I build data ready to use by gam funnction for smoothing
r1_IFG <- read.csv("C:/Users/Muhannad/Desktop/multvariat/R_codes/prediabetes/pre_sub/Data/ldl2.csv",header=TRUE, sep=",")
r1_IFG$SBP=(r1_IFG$BPXSY1+r1_IFG$BPXSY2+r1_IFG$BPXSY3)/3

r1_IFG=data.frame(r1_IFG$RIDAGEYR ,r1_IFG$SBP)
r1_IFG=na.omit(r1_IFG)
#loading the Splines Packages
require(splines)
#ISLR contains the Dataset
require(ISLR)
#attach(Wage) #attaching Wage dataset
agelims<-range(r1_IFG$r1_IFG.RIDAGEYR)
#c(min(as.vector(r1_IFG.RIDAGEYR)), as.vector(max(r1_IFG.RIDAGEYR)))

#Generating Test Data
age.grid<-seq(from=agelims[1], to = agelims[2])
#Fitting smooth cubic splin
#3 cutpoints at ages 25 ,50 ,60
fit<-lm(r1_IFG[,2] ~ bs(r1_IFG[,1],knots = c(25,40,60)),data = r1_IFG )
summary(fit)


#Plotting the Regression Line to the scatterplot   
plot(r1_IFG$r1_IFG.RIDAGEYR,r1_IFG$r1_IFG.SBP,col="grey",xlab="Age",ylab="Wages")
points(age.grid,predict(fit,newdata = list(age=age.grid)),col="darkgreen",lwd=2,type="l")
#adding cutpoints
abline(v=c(25,40,60),lty=2,col="darkgreen")

###############

#fitting smoothing splines using smooth.spline(X,Y,df=...)
fit1<-smooth.spline(r1_IFG[,1],r1_IFG[,2],df=16) #16 degrees of freedom
#Plotting both cubic and Smoothing Splines 
plot(r1_IFG[,1],r1_IFG[,2],col="grey",xlab="Age",ylab="Wages")
points(age.grid,predict(fit,newdata = list(age=age.grid)),col="darkgreen",lwd=2,type="l")
#adding cutpoints
abline(v=c(25,40,60),lty=2,col="darkgreen")
lines(fit1,col="red",lwd=2)
legend("topright",c("Smoothing Spline with 16 df","Cubic Spline"),col=c("red","darkgreen"),lwd=2)

fit2<-smooth.spline(r1_IFG[,1],r1_IFG[,2],cv = TRUE)
fit2

#It selects $\lambda=0.0279$ and df = 6.794596 as it is a Heuristic and can take various values for how rough the #function is
plot(age,wage,col="grey")
#Plotting Regression Line
lines(fit2,lwd=2,col="purple")
legend("topright",("Smoothing Splines with 6.78 df selected by CV"),col="purple",lwd=2)



#############


Farm_CHD=function(age,sbp,LR)
{
  mu_hat=15.5305+28.4441-(1.479+14.4588)*log(age)+1.8515*(log(age))^2-0.9119*log(sbp)-0.2767-0.7181*log(LR)
  
  sigma_hat=exp(0.9145-0.2784*mu_hat)
  u_hat=(log(1)-mu_hat)/sigma_hat
  prob=1-exp(-exp(u_hat))
  return(prob)
}
Farm_CHD(55,135,230/48)



########
#loading the Splines Packages
require(splines)
#ISLR contains the Dataset
require(ISLR)
attach(Wage) #attaching Wage dataset
#?Wage #for more details on the dataset
agelims<-range(age)
#Generating Test Data
age.grid<-seq(from=agelims[1], to = agelims[2])


#3 cutpoints at ages 25 ,50 ,60
fit<-lm(wage ~ bs(age,knots = c(25,40,60)),data = Wage )
summary(fit)

#Plotting the Regression Line to the scatterplot   
plot(age,wage,col="grey",xlab="Age",ylab="Wages")
points(age.grid,predict(fit,newdata = list(age=age.grid)),col="darkgreen",lwd=2,type="l")
#adding cutpoints
abline(v=c(25,40,60),lty=2,col="darkgreen")
