###########################################################################
## FNR
##
#       function simulating cases in Figure 2
#
# Input: - sim : the number of replications
#        - n: sample size 
#        - p: dimension of covariates
#        - sigrho: coefficient of covariates
#        - CensRateset: setup for censoring rate
#        - Boundset: bounded censoring time or not
#        - q_gam: parameter of r_n in the remark of section 2.1, page 5. 
#
# Output: - P_j and P_all in page 11 in the manuscript.
#
###########################################################################

#source basic functions 
source("./simulation/BasicFunctions.R")

FNR <- function(sim = 200,n = 200,p = 5000,sigrho = 0.2,CensRateset = 20,Boundset = 1,q_gam){
  set.seed(999)
  
  T_sett = rep(1:20)
  gaa_rec1 = matrix(NA,nrow = sim,ncol = length(q_gam))
  gaa_rec2 = matrix(NA,nrow = sim,ncol = length(q_gam))
  
  m = 1
  while(m<= sim){
    try({
      
      ###   Generate data begin ###
      ## generate covariate matrix ##
      
      Sig.X <- matrix(NA,nrow =  p,ncol =  p)
      rhos = sigrho^(seq(0,p-1,1))
      
      for (i in 1:p) {
        Sig.X[i,] = c(rep(0,i-1),rhos[1:(p-i+1)])
      }
      Sig.X = Sig.X+t(Sig.X)
      diag(Sig.X) <- 1
      
      p00=1:50
      p01=1:500
      #parameter setting
      s = 20
      beta1 <- c(rep(0.5,s),rep(0,(p-s)))
      beta2 <- -beta1
      S000 <- 1:s
      S111 <- (s+1):p
      q1 <- round(1/2*(n/log(n)))
      
      Z <- rmvnorm(n,mean = rep(0,p),Sig.X)
      
      ## generate competing risks indicator ##
      
      
      pp = 0.5
      
      p1 <- 1-(1-pp)^(exp(Z%*%beta1))
      eps <- c()
      eps <- sapply(1:n,function(i) rbinom(1,1,p1[i]))
      kappa1 <- which(eps == 1)
      kappa2 <- which(eps == 0)
      eps[kappa2] <- 2
      
      
      ######## generate true failure time  ############
      
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
          C = runif(n,0,quantile(TT,0.80)) else if(CensRateset == 60)
            C = runif(n,0,quantile(TT,0.60))
      }else{
        if(CensRateset == 20)  C = rexp(n,rate =  (1/quantile(TT,0.92))) else if(CensRateset == 40) 
          C = rexp(n,rate =  (1/quantile(TT,0.72))) else if(CensRateset == 60)
            C = rexp(n,rate =  (1/quantile(TT,0.55)))
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
        #
        results <- optimize(ssLogL,c(-5,5))
        c(results$minimum,results$objective)
      }
      
      #calculate marginal utilities for data
      MargResults <- sapply(p01,MargUtility)
      Utility <- MargResults[2,]
      Betas <- MargResults[1,]
      InfoResults <- sapply(1:p,function(i){
        sdata <- data[,c(1:3,3+i)]
        sbeta <- Betas[i]
        Info(sdata,sbeta)
      })
      
      ## FNR for CRSIS Begin ##
      recgam1 = c()
      for (i in 1:length(q_gam)) {
        gamman <- qnorm(1-(q_gam[i]/2))
        SubS_gam <- which(abs(InfoResults^0.5*Betas)>= gamman)
        dif_gam = setdiff(T_sett,SubS_gam)
        recgam1[i] = length(dif_gam)/length(SubS_gam)
      }
      
      gaa_rec1[m,] = recgam1
      write.csv(gaa_rec1,paste("./simulation/Figure2/FNR",'p',p,'+','rho',sigrho*10,'+',
                               'C',CensRateset,'+','B',Boundset,'+','CR-SIS','.csv',sep = ''))
      
      ## FNR for CRSIS End ##
      
      ## CRCSIS Begin ##
      
      Utility = abs(Utility)
      orderset.utility <- order(Utility)
      MMS.u <- max(which(orderset.utility%in%S000)) 
      P1.u <- (sum(S000%in%orderset.utility) == s)
      SubS <- orderset.utility[1:q1]
      ZFP = Z[,SubS]
      ftime <- X
      fstatus <- deleps
      cov <- ZFP
      try(fitsc1 <- crrp(ftime, fstatus, cov, penalty = "SCAD"),silent = T)
      betasc1 <-  fitsc1$beta[, which.min(fitsc1$BIC)]
      Shat.u <- SubS[which(betasc1!= 0)]
      Betahat.u <- betasc1[which(betasc1!= 0)]
      #sget conditional subset
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
        results <- optim(rep(0,s1+1),ssLogL,method =  "L-BFGS-B")
        c(results$par[s1+1],results$value,i)
      }
      
      #calculate conditional utilities
      CondResults <- sapply((p00)[-Subs1],CondUtility)
      Utility2 <- CondResults[2,]
      Utility2 <- abs(Utility2)
      Betas2 <- CondResults[1,]
      Betas3 = rep(0,p)
      Subbbb = intersect(Subs1,rep(1:p))
      Betas3[Subbbb] = Betas[Subbbb]
      Sets <- CondResults[3,]
      Betas3[Sets] = Betas2
      
      InfoResults2 <- sapply((1:50)[-Subs1],function(i){
        sdata <- data[,c(1:3,3+i)]
        sbeta <- Betas3[i]
        Info(sdata,sbeta)
      })
      
      ## FNR for CRCSIS Begin ##
      
      recgam2 = c()
      for (i in 1:length(q_gam)) {
        gamman2 <- qnorm(1-(q_gam[i]/2))
        SubS_gam2 <- c(Subbbb,Sets[which(abs(InfoResults2^0.5*Betas2)>= gamman2)])
        dif_gam2 = setdiff(T_sett,SubS_gam2)
        recgam2[i] = length(dif_gam2)/length(SubS_gam2)
      }
      
      ## FNR for CRCSIS End ##
      gaa_rec2[m,] = recgam2  
      
      
      write.csv(gaa_rec2,paste("./simulation/Figure2/FNR",'p',p,'+','rho',sigrho*10,'+',
                               'C',CensRateset,'+','B',Boundset,'+','CR-CSIS','.csv',sep = ''))
      
      
      ########### Compute and Save Average results  ########
      
      write.csv(colMeans(na.omit(gaa_rec1)),paste("./simulation/Figure2/FNR",'p',p,'+','rho',sigrho*10,'+',
                                                  'C',CensRateset,'+','B',Boundset,"+",'CR-SIS+Average','.csv',sep = ''))
      write.csv(colMeans(na.omit(gaa_rec2)),paste("./simulation/Figure2/FNR",'p',p,'+','rho',sigrho*10,'+',
                                                  'C',CensRateset,'+','B',Boundset,"+",'CR-CSIS+Average','.csv',sep = ''))
      data.frame(CRSIS = colMeans(na.omit(gaa_rec1)),CRCSIS = colMeans(na.omit(gaa_rec2)))
      
      m = m+1},silent = TRUE)}
}

