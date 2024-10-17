###########################################################################
## Screening1
##
#       function simulating cases in Table 1
#
# Input: - sim : the number of replications
#        - n: sample size 
#        - p: dimension of covariates
#        - Rhoset: coefficient of covariates
#        - CensRateset: setup for censoring rate
#        - Boundset: bounded censoring time or not
#
# Output: - P_j and P_all in page 11 in the manuscript.
#
###########################################################################

#source basic functions 
source("./simulation/BasicFunctions.R")

Screening1 <- function(sim = 200,n = 200,p = 5000,Rhoset = 0.2,CensRateset = 60,Boundset = 1){
  set.seed(999)
  sle_scr1 = matrix(NA,nrow = sim,ncol = 6)#CR-SIS result matrix 
  sle_scr2 = matrix(NA,nrow = sim,ncol = 6)#CR-CSIS result matrix
  sle_scr3 = matrix(NA,nrow = sim,ncol = 6)#CR-SJS result matrix
  sle_scr4 = matrix(NA,nrow = sim,ncol = 6)#crCRS result matrix
  sle_scr5 = matrix(NA,nrow = sim,ncol = 6)#crSIRS result matrix
  sle_scr6 = matrix(NA,nrow = sim,ncol = 6)#SIS result matrix
  
  ### Setting covariance matrix begin ###
  
  Sig.X <- matrix(NA,nrow  = p,ncol  = p)
  rhos = Rhoset^(seq(0,p-1,1))
  
  for (i in 1:p) {
    Sig.X[i,] = c(rep(0,i-1),rhos[1:(p-i+1)])
  }
  Sig.X = Sig.X+t(Sig.X)
  diag(Sig.X) <- 1
  
  #set parameters 
  s = 5
  S000 <- 1:s
  S111 <- (s+1):p
  q <- round(n/log(n))
  beta1 <- c(4,4,4,-3*sqrt(2),-4.5,rep(0,(p-s)))
  beta2 <- -beta1
  
  ### Setting covariance matrix end ###
  
  m = 1
  while(m<= sim){
    try({
      
      ###   Generate data begin ###
      ## generate covariate matrix ##
      Z <- rmvnorm(n,mean = rep(0,p),sigma  = Sig.X,method  = "chol")
      
      ## generate competing risks indicator ##
      
      pp = 0.3
      p1 <- 1-(1-pp)^(exp(Z%*%beta1))
      eps <- c()
      eps <- sapply(1:n,function(i) rbinom(1,1,p1[i]))
      kappa1 <- which(eps == 1)
      kappa2 <- which(eps == 0)
      eps[kappa2] <- 2
      
      ######## generate true failure time ############
      
      #CIF1I: function to generate F_1(t|Z)
      CIF1I <- function(u,z){-log(1-1/pp*(1-exp(exp(-z%*%beta1)*log(1-u*(1-(1-pp)^(exp(z%*%beta1)))))))}
      #CIF2I: function to generate F_2(t|Z)
      CIF2I <- function(u,z){-log(1-u)*exp(-z%*%beta2)}
      TT <- c()
      u <- runif(n,0,1)
      TT[kappa1] <- sapply(kappa1,function(i) CIF1I(u[i],Z[i,]))
      TT[kappa2] <- sapply(kappa2,function(i) CIF2I(u[i],Z[i,]))
      
      ######## generate censoring time ############
      if(Boundset){
        if(CensRateset == 20)  C = runif(n,0,quantile(TT,0.95)) else if(CensRateset == 40)
          C = runif(n,0,quantile(TT,0.82)) else if(CensRateset == 60)
            C = runif(n,0,quantile(TT,0.60))
      }else{
        if(CensRateset == 20)  C = rexp(n,rate  = (1/quantile(TT,0.92))) else if(CensRateset == 40) 
          C = rexp(n,rate  = (1/quantile(TT,0.72))) else if(CensRateset == 60)
            C = rexp(n,rate  = (1/quantile(TT,0.55)))
      }
      
      X = pmin(TT,C)
      CR <- sum(TT>C)/n
      delta = as.numeric(X == TT)
      deleps = delta*eps
      
      data <- data.frame(X = X,Del = delta,DC = deleps,Z = Z)
      ### Generate data end ###
      
      n <- dim(data)[1]
      p <- dim(data)[2]-3
      
      Z <- as.matrix(data[,-c(1:3)])
      W <- WeightMatrix(data)
      R1 <- RSets(data)
      K1 <- Kappa(data)
      E1 <- ESets(data)
      p00=1:100
      #### CRSIS Begin ###
      
      ## log-likelihood function ##
      # Input: - sdata : a subset of data
      #        - sbeta0 : coefficient parameter
      # Output: value of log-likelihood 
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
      
      ## Information function ###
      # Input: - sdata : a subset of data
      #        - sbeta0 : coefficient parameter
      # Output: value of information 
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
      
      #######   Marginal LogLikelihood function ######
      # Input: - i : index of covariate
      # Output: - results$minimum : value of minimum of log-likelihood 
      #         - results$objective : value of minimized log-likelihood (utility)
      MargUtility <- function(i){
        sdata <- data[,c(1:3,3+i)]
        ssLogL <- function(sbeta) -sLogLiklihood(sdata,sbeta)
        results <- optimize(ssLogL,c(-5,5))
        c(results$minimum,results$objective)
      }
      
      #calculate marginal utilities for data
      MargResults <- sapply(1:p,MargUtility)
      Utility <- MargResults[2,]
      Betas <- MargResults[1,]
      orderset.utility <- order(Utility)
      MMS.u <- max(which(orderset.utility%in%S000)) 
      P1.u <- (sum(S000%in%orderset.utility) == s)
      #select variable
      SubS <- orderset.utility[1:q]
      
      
      for (i in 1:s) {
        sle_scr1[m,i] = ifelse(i%in%SubS,1,0)
      }
      sle_scr1[m,6] = ifelse(sum(sle_scr1[m,1:s]) == s,1,0)
      filenames = paste("./simulation/Table1/",'p',p,'+','rho',Rhoset*10,'+',
                        'C',CensRateset,'+','B',Boundset,"+",'CRSIS','.csv',sep = '')
      write.csv(data.frame(sle_scr1),filenames)
      
      ### CRSIS End ###
      
      ### CRCSIS Begin ###
      
      ZFP = Z[,SubS]
      ftime  <-  X
      fstatus  <-  deleps
      cov  <-  ZFP
      try(fitsc1  <-  crrp(ftime, fstatus, cov, penalty = "SCAD"),silent = T)
      
      betasc1 <-  fitsc1$beta[, which.min(fitsc1$BIC)]
      Shat.u <- SubS[which(betasc1!= 0)]
      Betahat.u <- betasc1[which(betasc1!= 0)]
      # get the conditional subset
      Subs1 <- Shat.u
      s1 <- length(Subs1)
      
      ## Conditional likelihood ###
      # Input: - i : index of covariate
      # Output: - results$minimum : value of minimum of conditional log-likelihood 
      #         - results$objective : value of minimized conditional log-likelihood (utility)
      #         - i : index of covariate
      CondUtility <- function(i){
        sdata <- data[,c(1:3,3+c(Subs1,i))]
        ssLogL <- function(sbeta) -sLogLiklihood(sdata,sbeta)
        results <- optim(rep(0,s1+1),ssLogL,method  = ifelse(CensRateset == 60,"Nelder-Mead","L-BFGS-B"))
        c(results$par[s1+1],results$value,i)
      }
      #calculate conditional utilities
      CondResults <- sapply((p00)[-Subs1],CondUtility)
      Utility <- CondResults[2,]
      Utility <- abs(Utility)
      Betas <- CondResults[1,]
      Sets <- CondResults[3,]
      orderset.utility <- order(Utility)
      #get the selected variables
      SubS <- c(Subs1,Sets[orderset.utility[1:(q-length(Subs1))]])
      
      for (i in 1:s) {
        sle_scr2[m,i] = ifelse(i%in%SubS,1,0)
      }
      sle_scr2[m,6] = ifelse(sum(sle_scr2[m,1:s]) == s,1,0)
      
      filenames = paste("./simulation/Table1/",'p',p,'+','rho',Rhoset*10,'+',
                        'C',CensRateset,'+','B',Boundset,"+",'CRCSIS','.csv',sep = '')
      write.csv(data.frame(sle_scr2),filenames)
      
      ### CRCSIS  End ###
      
      ###### CRSJS Begin  #########
      r = 1
      Gamma.beta <- function(data,beta){
        n <- dim(data)[1]
        #Z <- as.matrix(data[,-c(1:3)])
        Eta <- eta.beta(data,beta)
        ExpEta <- exp(Eta)
        
        WeightProdExp <- sapply(1:n,function(i){
          j <- R1[[i]]
          sum(W[j,i]%*%ExpEta[j])
        })
        LogL <- sum(K1%*%Eta)-sum(K1%*%log(WeightProdExp))
        
        term1.dl <- t(Z)%*%K1
        
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
        Gamma1 <- beta+1*DiagW^(-1)*DLogL
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
          results <- optim(rep(0,sp),ssLogL,method  = "Nelder-Mead")
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
      for (i in 1:s) {
        sle_scr3[m,i] = ifelse(i%in%SubS,1,0)
      }
      sle_scr3[m,6] = ifelse(sum(sle_scr3[m,1:s]) == s,1,0)
      
      filenames = paste("./simulation/Table1/",'p',p,'+','rho',Rhoset*10,'+',
                        'C',CensRateset,'+','B',Boundset,"+",'CRSJS','.csv',sep = '')
      write.csv(data.frame(sle_scr3),filenames)
      
      ### CRSJS End ###
      
      ###### crCRS  Begin #########
      
      kmsurvival  <-  survfit(Surv(data$X,) ~ 1)
      KMresults <- summary(kmsurvival)
      
      KMestimate <- function(t){
        KMresults$surv[max(which((KMresults$time)<= max(t,min(KMresults$time))))]
      }
      GGfit = sapply(1:n, function(i) KMestimate(X[i]))
      
      crCRS = c()
      for (j in 1:p) {
        rec_3 = c()
        for (i in 1:n) {
          deltta = ifelse(TT<= C,1,0)
          rec_1 = ifelse(X<= X[i]&eps == 1,1,0)
          rec_2 = 1-GGfit
          rec_3[i] = sum((deltta*rec_1)/rec_2)/n
        }
        crCRS[j] = (sum(Z[,j]*rec_3)/n)^2
      }
      crCRS1 = order(crCRS,decreasing  = TRUE)[1:q]
      
      for (i in 1:s) {
        sle_scr4[m,i] = ifelse(i%in%crCRS1,1,0)
      }
      sle_scr4[m,6] = ifelse(sum(sle_scr4[m,1:s]) == s,1,0)
      
      filenames = paste("./simulation/Table1/",'p',p,'+','rho',Rhoset*10,'+',
                        'C',CensRateset,'+','B',Boundset,"+",'crCRS','.csv',sep = '')
      write.csv(data.frame(sle_scr4),filenames)
      
      ### crCRS End ###
      
      ###### crSIRS  Begin #########
      
      GGfit[which(GGfit == 0)] = mean(GGfit)
      crSIRS = c()
      for (j in 1:p) {
        rec_3 = c()
        for (i in 1:n) {
          deltta = ifelse(TT<= C,1,0)
          rec_1 = ifelse(X<= X[i]&eps == 1,1,0)
          rec_3[i] = (sum((Z[,j]*(deltta*rec_1))/GGfit)/n)^2
        }
        crSIRS[j] = sum(rec_3)/n
      }
      crSIRS1 = order(crSIRS,decreasing  = TRUE)[1:q]
      for (i in 1:s) {
        sle_scr5[m,i] = ifelse(i%in%crSIRS1,1,0)
      }
      sle_scr5[m,6] = ifelse(sum(sle_scr5[m,1:s]) == s,1,0)
      
      filenames = paste("./simulation/Table1/",'p',p,'+','rho',Rhoset*10,'+',
                        'C',CensRateset,'+','B',Boundset,"+",'crSIRS','.csv',sep = '')
      write.csv(data.frame(sle_scr5),filenames)
      
      ###### crSIRS  End #########
      
      ###### SIS Begin  #########
      
      SISY = X[which(TT<= C)]
      SISX = Z[which(TT<= C),]  
      crSIS = MFSIS(SISX,SISY,method = "SIS",nsis  = q) 
      
      for (i in 1:s) {
        sle_scr6[m,i] = ifelse(i%in%crSIS,1,0)
      }
      sle_scr6[m,6] = ifelse(sum(sle_scr6[m,1:s]) == s,1,0)
      
      filenames = paste("./simulation/Table1/",'p',p,'+','rho',Rhoset*10,'+',
                        'C',CensRateset,'+','B',Boundset,"+",'SIS','.csv',sep = '')
      write.csv(data.frame(sle_scr6),filenames)
      
      ###### SIS End #########
      
      m = m+1
    },silent = TRUE)
    
    ########### Compute and Save Average results  ########
    
    Average <- data.frame(
      CRSIS = colMeans(na.omit(sle_scr1)),
      CRCSIS = colMeans(na.omit(sle_scr2)),
      CRSJS = colMeans(na.omit(sle_scr3)),
      crCRS = colMeans(na.omit(sle_scr4)),
      crSIRS = colMeans(na.omit(sle_scr5)),
      SIS = colMeans(na.omit(sle_scr6)),
      row.names  = c('P1','P2','P3','P4','P5','Pa')
    )
    filenames = paste("./simulation/Table1/",'p',p,'+','rho',Rhoset*10,'+',
                      'C',CensRateset,'+','B',Boundset,"+",'Average','.csv',sep = '')
    write.csv(Average,filenames)
    
  }
  
  ########### Compute and Save Average results  ########
  Average <- data.frame(
    CRSIS = colMeans(na.omit(sle_scr1)),
    CRCSIS = colMeans(na.omit(sle_scr2)),
    CRSJS = colMeans(na.omit(sle_scr3)),
    crCRS = colMeans(na.omit(sle_scr4)),
    crSIRS = colMeans(na.omit(sle_scr5)),
    SIS = colMeans(na.omit(sle_scr6)),
    row.names  = c('P1','P2','P3','P4','P5','Pa')
  )
  filenames = paste("./simulation/Table1/",'p',p,'+','rho',Rhoset*10,'+',
                    'C',CensRateset,'+','B',Boundset,"+",'Average','.csv',sep = '')
  write.csv(Average,filenames)
}






