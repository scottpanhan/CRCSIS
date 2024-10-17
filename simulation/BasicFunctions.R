### Loading required packages ###

library("survival")
library("mvtnorm")
library("lars")
library("crrp")
library("MFSIS")
library("evd")
library("cmprsk")
library("pec")
library("riskRegression")

##########  Basic Functions #############
###########################################################################
## RiskSet
##
#       function getting the risk set
#
# Input: - data : dataset containing X, delta*epsilon, and Z
#        - i: index of covariate
#
# Output: - Risk set for i-th covariate
#
###########################################################################


RiskSet <- function(data,i){
  which(data$X>= data$X[i]|(data$X<= data$X[i]&data$DC!= 1))
}

###########################################################################
## RSets
##
#       function getting the risk sets
#
# Input: - data : dataset containing X, delta*epsilon, and Z
#
# Output: - Risk sets for every covariate
#
###########################################################################

RSets <- function(data){
  n <- dim(data)[1]
  sapply(1:n,function(i) RiskSet(data,i))
}

###########################################################################
## ESet
##
#       function getting the E_i set defined 
#       in Equation(3.8) of Li et. al.(2018)
#
# Input: - data : dataset containing X, delta*epsilon, and Z
#        - i: index of covariate
#
# Output: - E_i set for i-th covariate
#
###########################################################################

ESet <- function(data,i){
  which(data$DC[i]  == 2|(data$X<= data$X[i]&data$DC[i]!= 2))
}


###########################################################################
## ESets
##
#       function getting all E_i sets
#
# Input: - data : dataset containing X, delta*epsilon, and Z
#
# Output: - E_i sets for every covariate
#
###########################################################################
ESets <- function(data){
  n <- dim(data)[1]
  sapply(1:n,function(i) ESet(data,i))
}

###########################################################################
## Kappa
##
#       function getting all observed failures due to risk 1
#
# Input: - data : dataset containing X, delta*epsilon, and Z
#
# Output: - indicator of failure due to risk 1
#
###########################################################################
Kappa <- function(data){
  data$DC  == 1
}

###########################################################################
## WeightMatrix
##
#       function calculating weights w_j(X_k) defined 
#       in Equation(2.2) of Li et. al(2018)
#
# Input: - data : dataset containing X, delta*epsilon, and Z
#
# Output: - Weight matrix
#
###########################################################################
WeightMatrix <- function(data){
  n <- dim(data)[1]
  
  X <- data$X
  Del <- data$Del
  
  rEquiv <- function(i,t){
    (t<= X[i]|(Del[i]  == 1&X[i]<= t))
  }
  
  time <- data$X
  event <- 1-data$Del
  kmsurvival  <-  survfit(Surv(data$X,) ~ 1)
  KMresults <- summary(kmsurvival)
  
  KMestimate <- function(t){
    KMresults$surv[max(which((KMresults$time)<= max(t,min(KMresults$time))))]
  }
  
  WeightFunction <- function(i,t){
    w <- NULL
    if(KMestimate(min(X[i],t))>0)
      w  = rEquiv(i,t)*KMestimate(t)/KMestimate(min(X[i],t))
    else w  = rEquiv(i,t)
    w
  }
  
  sapply(1:n,function(k){
    sapply(1:n,function(i){WeightFunction(i,X[k])})})
}

###########################################################################
## eta.beta
##
#       function calculating eta defined 
#       in Equation(3.3) of Li et. al.(2018)
#
# Input: - data : dataset containing X, delta*epsilon, and Z
#        - beta : coefficent
#
# Output: - Weight matrix
#
###########################################################################
eta.beta <- function(data,beta){
  z <- as.matrix(data[,-(1:3)])
  z%*%beta
}