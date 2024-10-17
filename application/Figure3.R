source("./simulation/BasicFunctions.R")

###### Load data  #####
treat = read.csv("./application/data.csv",header = TRUE)# load data set
treat = data.frame(treat)

######  Data  processing #####
### define values for variables 
aa1 = which(is.na(treat$AGE))
aa2 = which(is.na(treat$WHOGrage))
aa3 = which(is.na(treat$DiseaseStage))
aa4 = which(is.na(treat$SEX))
aa5 = which(is.na(treat$Treatment))
SET_NA = union(union(union(union(aa1,aa2),aa3),aa4),aa5)
treat = treat[-SET_NA,]
ad_sex = c()
for (i in 1:nrow(treat)) {
  ad_sex[i] = ifelse(treat$SEX[i] == "M"||treat$SEX[i] == " M ",1,2)
}
ad_age = treat$AGE
ad_stage = c()
for (i in 1:nrow(treat)) {
  ad_stage[i] = ifelse(treat$DiseaseStage[i] == "pT1"||treat$DiseaseStage[i] == "pT1*"
                       ||treat$DiseaseStage[i] == "pT1a"||treat$DiseaseStage[i] == "pT1b",1,2)
}
ad_grade = c()
for (i in 1:nrow(treat)) {
  ad_grade[i] = ifelse(treat$WHOGrage[i] == "HIGH"||treat$WHOGrage[i] == "HIGH*",1,2)
}
ad_treat = c(unclass(as.factor(treat$Treatment)))
X = pmin(treat$FollowUp,treat$PFS)
delta_22 = rep(1,nrow(treat))
delta_22[which(is.na(treat$Time2progress))] = 0
delta = delta_22
eps = ifelse(treat$DeathCause == 1,1,2)
deleps = delta*eps
Z = cbind(ad_sex,ad_age,ad_stage,ad_grade,ad_treat,treat[,-rep(1:12)])

######## Sample training data ############
set.seed(123)
sam_tr = sample(1:336,335)
X_tr = X[sam_tr]
delta_tr = delta[sam_tr]
eps_tr = eps[sam_tr]
deleps_tr = deleps[sam_tr]
Z_tr = Z[sam_tr,]
data <- data.frame(X = X_tr,Del = delta_tr,DC = deleps_tr,Z = Z_tr)
n <- dim(data)[1]
p <- dim(data)[2]-3


Z_new <- as.matrix(data[,-c(1:3)])
W <- WeightMatrix(data)
R1 <- RSets(data)
K1 <- Kappa(data)
E1 <- ESets(data)

q <- round(n/log(n))

### calculate Kaplan Meier estiamtes ###
kmsurvival <- survfit(Surv(data$X,) ~ 1)
KMresults <- summary(kmsurvival)

KMestimate <- function(t){
  KMresults$surv[max(which((KMresults$time)<= max(t,min(KMresults$time))))]
}
GGfit = sapply(1:n, function(i) KMestimate(X[i]))
GGfit[which(GGfit == 0)] = mean(GGfit)

###### crSIRS Begin  #####

crSIRS = c()
for (j in 1:p) {
  rec_3 = c()
  for (i in 1:n) {
    deltta = delta_tr
    rec_1 = ifelse(X_tr<= X_tr[i]&eps_tr == 1,1,0)
    rec_3[i] = (sum((Z_new[,j]*(deltta*rec_1))/GGfit)/n)^2
  }
  crSIRS[j] = sum(rec_3)/n
}
crSIRS1 = order(crSIRS,decreasing  = TRUE)[1:q]

write.csv(data.frame(ID = crSIRS1,GENE = colnames(Z)[crSIRS1]),"./application/Figure3/crSIRS_genes.csv")

#Prediction error curves with three different penalty functions
ZFP = Z_new[,crSIRS1]

ii = 1
for(method in c("LASSO","MCP","SCAD")){
  ZFP = Z_new[,crSIRS1]
  ftime <- X_tr
  fstatus <- deleps_tr
  cov <- ZFP
  
  fitla <- crrp(ftime, fstatus, cov, penalty = method)
  betala <-  fitla$beta[, which.min(fitla$BIC)]
  betala = betala[which(betala!= 0)]
  nam_la = names(betala)
  rec_final = sapply(1:length(betala),function(i){which(colnames(Z_new) == nam_la[i])})
  
  
  Z_model = Z_new[,rec_final] 
  
  ftime <- X_tr
  fstatus <- deleps_tr
  cov <- Z_model
  Z_te = Z[-sam_tr,rec_final]
  cov = data.frame(cov)
  Z_te = data.frame(Z_te)
  
  model1 <-  crr(ftime,fstatus,cov)
  assign(paste0("sirp",ii),predict(model1,Z_te))
  ii = ii+1
}


write.csv(sirp1,"./application/Figure3/crSIRS+LASSO.csv")
write.csv(sirp2,"./application/Figure3/crSIRS+MCP.csv")
write.csv(sirp3,"./application/Figure3/crSIRS+SCAD.csv")

###### crSIRS End #####

#########  CRSIS Begin  #########

### log-likelihood function ####

sLogLiklihood <- function(sdata,sbeta0){
  sEta <- eta.beta(sdata,sbeta0)
  sExpEta <- exp(sEta)
  sW <- W
  sWeightProdExp <- sapply(1:n,function(i){
    j <- R1[[i]]
    sum(sW[j,i]%*%sExpEta[j])
  })
  sLogL <- sum(K1%*%sEta)-sum(K1%*%log(sWeightProdExp))
  sLogL
}


### Information function ###
Info <- function(sdata,sbeta0){
  sZ <- sdata[,-c(1:3)]
  sEta <- eta.beta(sdata,sbeta0)
  sExpEta <- exp(sEta)
  sW <- W
  sWeightZProdExp <- sapply(1:n,function(i){
    j <- R1[[i]]
    c(sum(sW[j,i]%*%sExpEta[j]),
      sum(sW[j,i]%*%(sZ[j]*sExpEta[j])),
      sum(sW[j,i]%*%(sZ[j]^2*sExpEta[j])))
  })
  
  Information <- sum(K1%*%((sWeightZProdExp[3,]*sWeightZProdExp[1,]-sWeightZProdExp[2,]^2)/(sWeightZProdExp[1,])^2))
  Information
}

#######   Marginal Log-Likelihood function ######

MargUtility <- function(i){
  sdata <- data[,c(1:3,3+i)]
  ssLogL <- function(sbeta) -sLogLiklihood(sdata,sbeta)
  results <- optimize(ssLogL,c(-5,5))
  c(results$minimum,results$objective)
}
#calculate marginal utility for data
MargResults <- sapply(1:p,MargUtility)
Utility <- MargResults[2,]
Betas <- MargResults[1,]
Utility = abs(Utility)
orderset.utility <- order(Utility)[1:q]
recood = order(Utility)[1:p]

#Prediction error curves with three different penalty functions
ZFP = Z_new[,orderset.utility]

ii = 4
for(method in c("LASSO","MCP","SCAD")){
  
  ftime <- X_tr
  fstatus <- deleps_tr
  cov <- ZFP
  
  fitla <- crrp(ftime, fstatus, cov, penalty = method)
  betala <-  fitla$beta[, 3]
  betala = betala[which(betala!= 0)]
  
  nam_la = names(betala)
  rec_final = sapply(1:length(betala),function(i){which(colnames(Z_new) == nam_la[i])})
  
  
  Z_model = Z_new[,rec_final] 
  
  ftime <- X_tr
  fstatus <- deleps_tr
  cov <- Z_model
  Z_te = Z[-sam_tr,rec_final]
  cov = data.frame(cov)
  Z_te = data.frame(Z_te)
  
  model1 <-  crr(ftime,fstatus,cov)
  
  assign(paste0("sirp",ii),predict(model1,Z_te))
  ii = ii+1
}



write.csv(sirp4,"./application/Figure3/CR-SIS+LASSO.csv")
write.csv(sirp5,"./application/Figure3/CR-SIS+MCP.csv")
write.csv(sirp6,"./application/Figure3/CR-SIS+SCAD.csv")


#########  CRSIS End  #########


#########  CRCSIS Begin  #########
ZFP = Z_new[,orderset.utility]
ftime <- X_tr
fstatus <- deleps_tr
cov <- ZFP
fitsc1 <- crrp(ftime, fstatus, cov, penalty = "SCAD")
betasc1 <-  fitsc1$beta[, which.min(fitsc1$BIC)]
Shat.u <- orderset.utility[which(betasc1!= 0)]
Betahat.u <- betasc1[which(betasc1!= 0)]
Subs1 <- Shat.u
s1 <- length(Subs1)

### Conditional likelihood  ###

CondUtility <- function(i){
  sdata <- data[,c(1:3,3+c(Subs1,i))]
  ssLogL <- function(sbeta) -sLogLiklihood(sdata,sbeta)
  #
  results <- optim(rep(0,s1+1),ssLogL,method  =  "L-BFGS-B")
  c(results$par[s1+1],results$value,i)
}

#calculate conditional utility for data
CondResults <- sapply(recood[-Subs1],CondUtility)
Utility <- CondResults[2,]
Utility <- abs(Utility)
Betas <- CondResults[1,]
Sets <- CondResults[3,]
orderset.utility <- order(Utility)
## select variables
SubS <- c(Subs1,Sets[orderset.utility[1:(q-length(Subs1))]])


ZFP = Z_new[,SubS]

#Prediction error curves with three different penalty functions
ii = 7
for(method in c("LASSO","MCP","SCAD")){
  
  ftime <- X_tr
  fstatus <- deleps_tr
  cov <- ZFP
  
  fitla <- crrp(ftime, fstatus, cov, penalty = method)
  betala <-  fitla$beta[,5]
  betala = betala[which(betala!= 0)]
  
  nam_la = names(betala)
  rec_final = sapply(1:length(betala),function(i){which(colnames(Z_new) == nam_la[i])})
  
  Z_model = Z_new[,rec_final] 
  
  ftime <- X_tr
  fstatus <- deleps_tr
  cov <- Z_model
  Z_te = Z[-sam_tr,rec_final]
  cov = data.frame(cov)
  Z_te = data.frame(Z_te)
  
  model1 <-  crr(ftime,fstatus,cov)
  assign(paste0("sirp",ii),predict(model1,Z_te))
  ii = ii+1
}


write.csv(sirp7,"./application/Figure3/CR-CSIS+LASSO.csv")
write.csv(sirp8,"./application/Figure3/CR-CSIS+MCP.csv")
write.csv(sirp9,"./application/Figure3/CR-CSIS+SCAD.csv")


#########  CRCSIS End  #########


#########  CRSJS Begin  #########

r = 1
Gamma.beta <- function(data,beta){
  n <- dim(data)[1]
  Eta <- eta.beta(data,beta)
  ExpEta <- exp(Eta)
  
  WeightProdExp <- sapply(1:n,function(i){
    j <- R1[[i]]
    sum(W[j,i]%*%ExpEta[j])
  })
  LogL <- sum(K1%*%Eta)-sum(K1%*%log(WeightProdExp))
  
  term1.dl <- t(Z_new)%*%K1
  
  NN <- which(K1 == 1)
  term21.dl <- sapply(NN,function(i){
    j <- R1[[i]]
    t(Z[j,])%*%c(W[j,i]*ExpEta[j])/WeightProdExp[i]})
  
  term2.dl <- apply(term21.dl,1,sum)
  
  DLogL <- term1.dl-term2.dl
  
  term11.ddl <- sapply(NN,function(i){
    j <- R1[[i]]
    t(Z[j,]^2)%*%c(W[j,i]*ExpEta[j])/WeightProdExp[i]
  })
  
  term1.ddl <- apply(term11.ddl,1,sum)
  
  term2.ddl <- apply((term21.dl)^2,1,sum)
  
  DDLogL <- term1.ddl-term2.ddl
  DiagW <- DDLogL
  Gamma1 <- beta+1/r*DiagW^(-1)*DLogL
  Gamma1
}

iteration <- 50

beta0 <- rep(0,p)
update.beta <- c()


Gamma1 <- Gamma.beta(data,beta0)
orderset.gamma <- order(abs(Gamma1),decreasing  =  TRUE)
ZeroSet <- orderset.gamma[(q+1):p]
NonZeroSet <- orderset.gamma[1:q]
BetaTilde <- c()
BetaTilde[NonZeroSet] <- Gamma1[NonZeroSet]
BetaTilde[ZeroSet] <- 0
SubS <- NonZeroSet
S0 <- 0

for(i in 1:iteration){
  if(sum(SubS%in%S0)/length(SubS)<0.95){
    S0 <- SubS
    sp <- length(SubS)
    
    sdata <- data[,c(1:3,3+SubS)]
    ssLogL <- function(sbeta) -sLogLiklihood(sdata,sbeta)
    results <- optim(rep(0,sp),ssLogL,"L-BFGS-B")
    Shat.u <- SubS
    Betahat.u <- results$par
    
    NonZeroSet <- Shat.u
    ZeroSet <- (1:p)[-Shat.u]
    update.beta[NonZeroSet] <- Betahat.u
    update.beta[ZeroSet] <- 0
    
    beta0 <- update.beta
    Gamma1 <- Gamma.beta(data,beta0)
    orderset.gamma <- order(abs(Gamma1),decreasing  =  TRUE)
    ZeroSet <- orderset.gamma[(q+1):p]
    NonZeroSet <- orderset.gamma[1:q]
    SubS <- NonZeroSet
  }
}

ZFP = Z_new[,SubS]

#Prediction error curves with three different penalty functions
ii = 10
for(method in c("LASSO","MCP","SCAD")){
  
  ftime <- X_tr
  fstatus <- deleps_tr
  cov <- ZFP
  
  fitla <- crrp(ftime, fstatus, cov, penalty = method)
  betala <-  fitla$beta[,6]
  betala = betala[which(betala!= 0)]
  
  nam_la = names(betala)
  rec_final = sapply(1:length(betala),function(i){which(colnames(Z_new) == nam_la[i])})
  
  Z_model = Z_new[,rec_final] 
  
  ftime <- X_tr
  fstatus <- deleps_tr
  cov <- Z_model
  Z_te = Z[-sam_tr,rec_final]
  cov = data.frame(cov)
  Z_te = data.frame(Z_te)
  
  model1 <-  crr(ftime,fstatus,cov)
  assign(paste0("sirp",ii),predict(model1,Z_te))
  ii = ii+1
  
}


write.csv(sirp10,"./application/Figure3/CR-SJS+LASSO.csv")
write.csv(sirp11,"./application/Figure3/CR-SJS+MCP.csv")
write.csv(sirp12,"./application/Figure3/CR-SJS+SCAD.csv")


#########  CRSJS End  #########
