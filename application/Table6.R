source("./simulation/BasicFunctions.R")

### Load data  ###
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
set.seed(71)
sam_tr = sample(rep(1:336),268)
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

q = 25
r = 1


### calculate Kaplan Meier estiamtes ###
kmsurvival  <- survfit(Surv(data$X,) ~ 1)
KMresults <- summary(kmsurvival)

KMestimate <- function(t){
  KMresults$surv[max(which((KMresults$time)<= max(t,min(KMresults$time))))]
}
GGfit = sapply(1:n, function(i) KMestimate(X[i]))
GGfit[which(GGfit == 0)] = mean(GGfit)

### crSIRS Begin ###
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



Z_model2 = Z_new[,crSIRS1] 
Z_model = Z[-sam_tr,crSIRS1]
ne_data = data.frame(time = X_tr,status = deleps_tr,Z_tr = Z_model2)
ne_data2 = data.frame(time = X[-sam_tr],status = deleps[-sam_tr],Z_tr = Z_model)
names(ne_data2) = names(ne_data)
formm = as.formula(paste("Hist(time,status)~",paste(names(ne_data2)[-c(1:2)],sep = "",collapse  = "+")))
csc  <- CSC(formm,data = ne_data)
fgr  <- FGR(formm,data = ne_data,cause  = 1)

### get IBS All ###

score <- Score(list("Cause-specific Cox" = csc,"Fine-Gray" = fgr),
               formula  = Hist(time,status)~1,
               data = ne_data2,times  = sort(ne_data2$time),
               plots  = "calibration",
               summary  = "ibs")

output_IBS <- score[["Brier"]][["score"]][,7][nrow(score[["Brier"]][["score"]])]

###get C-index IBS All###

score <- Score(list("Cause-specific Cox" = csc,
                    "Fine-Gray" = fgr),
               formula  = Hist(time,status)~1,
               data = ne_data2,times  = 100,
               plots  = "calibration",
               summary  = "risks")

output_Cindex <- unlist(score$AUC[1])[which(names(unlist(score$AUC[1])) == "score.AUC2")]

write.csv(data.frame(IBS = output_IBS,Cindex = output_Cindex),"./application/Table6/crSIRS_ALL.csv")

ftime  <- X_tr
fstatus  <- deleps_tr
cov  <- Z_model2

for(method in c("LASSO","MCP","SCAD")){
  fitla  <- crrp(ftime, fstatus, cov, penalty = method)
  betala <- fitla$beta[, which.min(fitla$BIC)]
  nam_la = names(ne_data2)[-c(1,2)][which(betala!= 0)]
  formm = as.formula(paste("Hist(time,status)~",paste(nam_la,sep = "",collapse  = "+")))
  csc  <- CSC(formm,data = ne_data)
  fgr  <- FGR(formm,data = ne_data,cause  = 1)
  
  ### get IBS  for each method ###
  
  score <- Score(list("Cause-specific Cox" = csc,"Fine-Gray" = fgr),
                 formula  = Hist(time,status)~1,
                 data = ne_data2,times  = sort(ne_data2$time),
                 plots  = "calibration",
                 summary  = "ibs")
  output_IBS <- score[["Brier"]][["score"]][,7][nrow(score[["Brier"]][["score"]])]
  
  ###get C-index  for each method ###
  
  score <- Score(list("Cause-specific Cox" = csc,
                      "Fine-Gray" = fgr),
                 formula  = Hist(time,status)~1,
                 data = ne_data2,times  = 100,
                 plots  = "calibration",
                 summary  = "risks")
  output_Cindex <- unlist(score$AUC[1])[which(names(unlist(score$AUC[1])) == "score.AUC2")]
  write.csv(data.frame(IBS = output_IBS,Cindex = output_Cindex),paste("./application/Table6/crSIRS_",method,".csv",sep = ''))
  
}
###  crSIRS End ###

###  CR-SJS Begin ###

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
orderset.gamma <- order(abs(Gamma1),decreasing  = TRUE)
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
    orderset.gamma <- order(abs(Gamma1),decreasing  = TRUE)
    ZeroSet <- orderset.gamma[(q+1):p]
    NonZeroSet <- orderset.gamma[1:q]
    SubS <- NonZeroSet
  }
}

Z_model2 = Z_new[,SubS] 
Z_model = Z[-sam_tr,SubS]
ne_data = data.frame(time = X_tr,status = deleps_tr,Z_tr = Z_model2)
ne_data2 = data.frame(time = X[-sam_tr],status = deleps[-sam_tr],Z_tr = Z_model)
names(ne_data2) = names(ne_data)
formm = as.formula(paste("Hist(time,status)~",paste(names(ne_data2)[-c(1:2)],sep = "",collapse  = "+")))
csc  <- CSC(formm,data = ne_data)
fgr  <- FGR(formm,data = ne_data,cause  = 1)

###get IBS All ###

score <- Score(list("Cause-specific Cox" = csc,"Fine-Gray" = fgr),
               formula  = Hist(time,status)~1,
               data = ne_data2,times  = sort(ne_data2$time),
               plots  = "calibration",
               summary  = "ibs")
output_IBS <- score[["Brier"]][["score"]][,7][nrow(score[["Brier"]][["score"]])]

###get C-index All ###

score <- Score(list("Cause-specific Cox" = csc,
                    "Fine-Gray" = fgr),
               formula  = Hist(time,status)~1,
               data = ne_data2,times  = 100,
               plots  = "calibration",
               summary  = "risks")
output_Cindex <- unlist(score$AUC[1])[which(names(unlist(score$AUC[1])) == "score.AUC2")]

write.csv(data.frame(IBS = output_IBS,Cindex = output_Cindex),"./application/Table6/CRSJS_ALL.csv")

ftime  <- X_tr
fstatus  <- deleps_tr
cov  <- Z_model2

for(method in c("LASSO","MCP","SCAD")){
  fitla  <- crrp(ftime, fstatus, cov, penalty = method)
  betala <- fitla$beta[, which.min(fitla$BIC)]
  nam_la = names(ne_data2)[-c(1,2)][which(betala!= 0)]
  if(length(nam_la) == 0){
    betala <- fitla$beta[, 6]
    nam_la = names(ne_data2)[-c(1,2)][which(betala!= 0)]
  }
  
  formm = as.formula(paste("Hist(time,status)~",paste(nam_la,sep = "",collapse  = "+")))
  csc  <- CSC(formm,data = ne_data)
  fgr  <- FGR(formm,data = ne_data,cause  = 1)
  
  ###get  IBS for each method ###
  
  score <- Score(list("Cause-specific Cox" = csc,"Fine-Gray" = fgr),
                 formula  = Hist(time,status)~1,
                 data = ne_data2,times  = sort(ne_data2$time),
                 plots  = "calibration",
                 summary  = "ibs")
  output_IBS <- score[["Brier"]][["score"]][,7][nrow(score[["Brier"]][["score"]])]
  
  ###get C-index for each method ###
  
  score <- Score(list("Cause-specific Cox" = csc,
                      "Fine-Gray" = fgr),
                 formula  = Hist(time,status)~1,
                 data = ne_data2,times  = 100,
                 plots  = "calibration",
                 summary  = "risks")
  output_Cindex <- unlist(score$AUC[1])[which(names(unlist(score$AUC[1])) == "score.AUC2")]
  
  write.csv(data.frame(IBS = output_IBS,Cindex = output_Cindex),paste("./application/Table6/CRSJS_",method,".csv",sep = ''))
  
}
###  CR-SJS End ###

### CRSIS Begin ###

## log-likelihood function ##

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

## Information function ##

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

#######   Marginal LogLikelihood ######
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


Z_model2 = Z_new[,orderset.utility] 
Z_model = Z[-sam_tr,orderset.utility]
ne_data = data.frame(time = X_tr,status = deleps_tr,Z_tr = Z_model2)
ne_data2 = data.frame(time = X[-sam_tr],status = deleps[-sam_tr],Z_tr = Z_model)
names(ne_data2) = names(ne_data)
formm = as.formula(paste("Hist(time,status)~",paste(names(ne_data2)[-c(1:2)],sep = "",collapse  = "+")))
csc  <- CSC(formm,data = ne_data)
fgr  <- FGR(formm,data = ne_data,cause  = 1)

### get IBS  All ###

score <- Score(list("Cause-specific Cox" = csc,"Fine-Gray" = fgr),
               formula  = Hist(time,status)~1,
               data = ne_data2,times  = sort(ne_data2$time),
               plots  = "calibration",
               summary  = "ibs")
output_IBS <- score[["Brier"]][["score"]][,7][nrow(score[["Brier"]][["score"]])]

###get C-index All ##
score <- Score(list("Cause-specific Cox" = csc,
                    "Fine-Gray" = fgr),
               formula  = Hist(time,status)~1,
               data = ne_data2,times  = 100,
               plots  = "calibration",
               summary  = "risks")
output_Cindex <- unlist(score$AUC[1])[which(names(unlist(score$AUC[1])) == "score.AUC2")]

write.csv(data.frame(IBS = output_IBS,Cindex = output_Cindex),"./application/Table6/CRSIS_ALL.csv")


ftime  <- X_tr
fstatus  <- deleps_tr
cov  <- Z_model2
for(method in c("LASSO","MCP","SCAD")){
  fitla  <- crrp(ftime, fstatus, cov, penalty = method)
  betala <- fitla$beta[, which.min(fitla$BIC)]
  nam_la = names(ne_data2)[-c(1,2)][which(betala!= 0)]
  if(length(nam_la) == 0){
    betala <- fitla$beta[, 6]
    nam_la = names(ne_data2)[-c(1,2)][which(betala!= 0)]
  }
  formm = as.formula(paste("Hist(time,status)~",paste(nam_la,sep = "",collapse  = "+")))
  csc  <- CSC(formm,data = ne_data)
  fgr  <- FGR(formm,data = ne_data,cause  = 1)
  
  ###get IBS for each method ###
  
  score <- Score(list("Cause-specific Cox" = csc,"Fine-Gray" = fgr),
                 formula  = Hist(time,status)~1,
                 data = ne_data2,times  = sort(ne_data2$time),
                 plots  = "calibration",
                 summary  = "ibs")
  output_IBS <- score[["Brier"]][["score"]][,7][nrow(score[["Brier"]][["score"]])]
  
  ##get C-index for each method ##
  
  score <- Score(list("Cause-specific Cox" = csc,
                      "Fine-Gray" = fgr),
                 formula  = Hist(time,status)~1,
                 data = ne_data2,times  = 100,
                 plots  = "calibration",
                 summary  = "risks")
  output_Cindex <- unlist(score$AUC[1])[which(names(unlist(score$AUC[1])) == "score.AUC2")]
  write.csv(data.frame(IBS = output_IBS,Cindex = output_Cindex),paste("./application/Table6/CRSIS_",method,".csv",sep = ''))
}
###  CRSIS End ###

### CRCSIS Begin###

ftime  <- X_tr
fstatus  <- deleps_tr
cov  <- Z_model2
fitla  <- crrp(ftime, fstatus, cov, penalty = "SCAD")
betala <- fitla$beta[, which.min(fitla$BIC)]
Subs1 <- orderset.utility[which(betala!= 0)]
s1 <- length(Subs1)

### Conditional Likelihood ###

CondUtility <- function(i){
  sdata <- data[,c(1:3,3+c(Subs1,i))]
  ssLogL <- function(sbeta) -sLogLiklihood(sdata,sbeta)
  results <- optim(rep(0,s1+1),ssLogL,method  = "L-BFGS-B")
  c(results$par[s1+1],results$value,i)
}

#calculate conditional utility for data
CondResults <- sapply((1:p)[-Subs1],CondUtility)
Utility <- CondResults[2,]
Utility <- abs(Utility)
Betas <- CondResults[1,]
Sets <- CondResults[3,]
orderset.utility <- order(Utility)
## select variables
SubS <- c(Subs1,Sets[orderset.utility[1:(q-length(Subs1))]])

Z_model2 = Z_new[,SubS] 
Z_model = Z[-sam_tr,SubS]
ne_data = data.frame(time = X_tr,status = deleps_tr,Z_tr = Z_model2)
ne_data2 = data.frame(time = X[-sam_tr],status = deleps[-sam_tr],Z_tr = Z_model)
names(ne_data2) = names(ne_data)
formm = as.formula(paste("Hist(time,status)~",paste(names(ne_data2)[-c(1:2)],sep = "",collapse  = "+")))
csc  <- CSC(formm,data = ne_data)
fgr  <- FGR(formm,data = ne_data,cause  = 1)

###get IBS All ###

score <- Score(list("Cause-specific Cox" = csc,"Fine-Gray" = fgr),
               formula  = Hist(time,status)~1,
               data = ne_data2,times  = sort(ne_data2$time),
               plots  = "calibration",
               summary  = "ibs")
output_IBS <- score[["Brier"]][["score"]][,7][nrow(score[["Brier"]][["score"]])]

##get C-index All ##

score <- Score(list("Cause-specific Cox" = csc,
                    "Fine-Gray" = fgr),
               formula  = Hist(time,status)~1,
               data = ne_data2,times  = 100,
               plots  = "calibration",
               summary  = "risks")

output_Cindex <- unlist(score$AUC[1])[which(names(unlist(score$AUC[1])) == "score.AUC2")]
write.csv(data.frame(IBS = output_IBS,Cindex = output_Cindex),"./application/Table6/CRCSIS_ALL.csv")


ftime  <- X_tr
fstatus  <- deleps_tr
cov  <- Z_model2
for(method in c("LASSO","MCP","SCAD")){
  fitla  <- crrp(ftime, fstatus, cov, penalty = method)
  betala <- fitla$beta[, which.min(fitla$BIC)]
  nam_la = names(ne_data2)[-c(1,2)][which(betala!= 0)]
  if(length(nam_la) == 0){
    betala <- fitla$beta[, 6]
    nam_la = names(ne_data2)[-c(1,2)][which(betala!= 0)]
  }
  formm = as.formula(paste("Hist(time,status)~",paste(nam_la,sep = "",collapse  = "+")))
  csc  <- CSC(formm,data = ne_data)
  fgr  <- FGR(formm,data = ne_data,cause  = 1)
  
  ###get  IBS  for each method  ###
  
  score <- Score(list("Cause-specific Cox" = csc,"Fine-Gray" = fgr),
                 formula  = Hist(time,status)~1,
                 data = ne_data2,times  = sort(ne_data2$time),
                 plots  = "calibration",
                 summary  = "ibs")
  output_IBS <- score[["Brier"]][["score"]][,7][nrow(score[["Brier"]][["score"]])]
  
  ##get C-index for each method ###
  score <- Score(list("Cause-specific Cox" = csc,
                      "Fine-Gray" = fgr),
                 formula  = Hist(time,status)~1,
                 data = ne_data2,times  = 100,
                 plots  = "calibration",
                 summary  = "risks")
  output_Cindex <- unlist(score$AUC[1])[which(names(unlist(score$AUC[1])) == "score.AUC2")]
  write.csv(data.frame(IBS = output_IBS,Cindex = output_Cindex),paste("./application/Table6/CRCSIS_",method,".csv",sep = ''))
}


