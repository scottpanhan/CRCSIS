###########################################################################
## PostCRSIS
##
#       function simulating post-screening result for CRSIS in Table 4
#
# Input: - sim : the number of replications
#        - n: sample size 
#        - p: dimension of covariates
#        - Caseset: indicator for setup for 4 cases in example 5.4
#
# Output: - C, IC, C-fit and MMSE in page 11 in the manuscript.
#
###########################################################################

#source basic functions 
source("./simulation/BasicFunctions.R")

PostCRSIS <- function(sim = 200,n = 200,p = 3000,Caseset = 1){
  set.seed(1234)
  sle_scr = matrix(NA,nrow = sim,ncol = 7)# screening result matrix
  las = matrix(NA,nrow = sim,ncol = 4)#LASSO result matrix
  sc = matrix(NA,nrow = sim,ncol = 4)#SCAD result matrix
  mc = matrix(NA,nrow = sim,ncol = 4)#MCP result matrix
  or = matrix(NA,nrow = sim,ncol = 4)#Oracle result matrix
  
  m = 1
  while(m<= sim){
    try({
      
      ###   Generate data begin ###
      ## generate covariate matrix ##
      #parameter setting
      
      if(Caseset == 1){
        Sig.X <- diag(rep(1,p))
        #set parameters 
        s = 6
        beta1 <- c(rep(2,s),rep(0,p-s))
        beta2 <- -beta1
        S000 <- 1:s
        S111 <- (s+1):p
      }else if(Caseset == 2){
        Sig.X <- matrix(rep(0.5,p^2),nr = p)
        diag(Sig.X) <- 1
        #set parameters 
        s = 6
        beta1 <- c(rep(2,s),rep(0,p-s))
        beta2 <- -beta1
        S000 <- 1:s
        S111 <- (s+1):p
      }else if(Caseset == 3){
        Sig.X <- diag(rep(1,p))
        #set parameters 
        s = 6
        a = 4*log(n)/n^(1/2)
        u0 = rbinom(s,1,0.5)
        u0 <- 2*u0-1
        z0 = rnorm(s)
        beta1 <- c(u0*(a+abs(z0)/4),rep(0,p-s))
        beta2 <- -beta1
        S000 <- 1:s
        S111 <- (s+1):p
      }else if(Caseset == 4){
        Sig.X <- matrix(rep(0.5,p^2),nr = p)
        diag(Sig.X) <- 1
        #set parameters 
        s = 6
        a = 4*log(n)/n^(1/2)
        u0 = rbinom(s,1,0.5)
        u0 <- 2*u0-1
        z0 = rnorm(s)
        beta1 <- c(u0*(a+abs(z0)/4),rep(0,p-s))
        beta2 <- -beta1
        S000 <- 1:s
        S111 <- (s+1):p
      }else if(Caseset == 5){
        Sig4 <- matrix(rep(0.5,p^2),nr = p)
        Sig4[4,] <- 1/sqrt(2)
        Sig4[,4] <- 1/sqrt(2)
        Sig4[5,] <- 0
        Sig4[,5] <- 0
        diag(Sig4) <- 1
        Sig.X = Sig4
        
        #set parameters 
        s = 5
        beta1 <- c(4,4,4,-3*sqrt(2),4,rep(0,p-s))
        beta2 <- -beta1
        S000 <- 1:s
        S111 <- (s+1):p
      }else if(Caseset == 6){
        Sig3 <- matrix(rep(0.5,p^2),nr = p)
        Sig3[4,] <- 1/sqrt(2)
        Sig3[,4] <- 1/sqrt(2)
        diag(Sig3) <- 1
        Sig.X = Sig3
        
        #set parameters 
        s = 4
        beta1 <- c(4,4,4,-4,rep(0,p-s))
        beta2 <- -beta1
        S000 <- 1:s
        S111 <- (s+1):p
      }
      
      
      q <- floor(n/log(n))
      Z <- rmvnorm(n,mean = rep(0,p),sigma  = Sig.X,method  = "chol")
      
      ## generate competing risks indicator ##
      pp = 0.3
      
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
      
      if(Caseset == 1){
        C = runif(n,0,quantile(TT,0.90))
      }else if(Caseset == 2){
        C = runif(n,0,quantile(TT,0.85))
      }else if(Caseset == 3){
        C = runif(n,0,quantile(TT,0.90))
      }else if(Caseset == 4){
        C = runif(n,0,quantile(TT,0.85))
      }else if(Caseset == 5){
        C = runif(n,0,quantile(TT,0.85))
      }else if(Caseset == 6){
        C = runif(n,0,quantile(TT,0.85))
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
        results <- optimize(ssLogL,c(-5,5))
        c(results$minimum,results$objective)
      }
      #calculate marginal utilities for data
      MargResults <- sapply(1:p,MargUtility)
      Utility <- MargResults[2,]
      Betas <- MargResults[1,]
      orderset.utility <- order(Utility)
      MMS.u <- max(which(orderset.utility%in%S000)) ## MMS
      P1.u <- (sum(S000%in%orderset.utility) == s)
      #select variable
      SubS <- orderset.utility[1:q]
      
      
      for (i in 1:s) {
        sle_scr[m,i] = ifelse(i%in%SubS,1,0)
      }
      sle_scr[m,7] = ifelse(sum(sle_scr[m,1:s]) == s,1,0)
      
      filenames = paste("./simulation/Table4/",'n',n,'+','p',p,'+','Case',Caseset,'CRSIS.csv',sep = '')
      write.csv(data.frame(sle_scr),filenames)
      
      ### CRSIS End ###
      
      ### post CRSIS Begin ###
      
      #generate subset of data
      ZFP = Z[,SubS]
      ftime <- X
      fstatus <- deleps
      cov <- ZFP
      
      T_sett = sapply(1:s,function(x) paste('Z.',x,sep = ''))
      
      ## LASSO Begin ##
      fitla <- crrp(ftime, fstatus, cov, penalty = "LASSO")
      betala <-  fitla$beta[, which.min(fitla$BIC)]
      betala = betala[which(betala!= 0)]
      nam_la = names(betala)
      # compute C
      las[m,1] = length(intersect(nam_la,T_sett))
      # compute IC
      las[m,2] = length(setdiff(nam_la,T_sett))
      # compute C-fit
      las[m,3] = ifelse(las[m,1] == s&las[m,2] == 0,1,0)
      las_beta = rep(0,p)
      names(las_beta) = paste("Z",1:p,sep  = ".")
      las_beta[sapply(1:length(nam_la),function(i){which(names(las_beta) == nam_la[i])})] = betala
      # compute MMSE
      las[m,4] = t(as.matrix(las_beta-beta1))%*%Sig.X%*%as.matrix(las_beta-beta1)
      
      filenames = paste("./simulation/Table4/",'n',n,'+','p',p,'+','Case',Caseset,'+CRSIS+LASSO.csv',sep = '')
      write.csv(data.frame(las),filenames)
      
      ## LASSO End ##
      
      ## SCAD Begin ##
      
      fitsc <- crrp(ftime, fstatus, cov, penalty = "SCAD")
      betasc <-  fitsc$beta[, which.min(fitsc$BIC)]
      betasc = betasc[which(betasc!= 0)]
      nam_sc = names(betasc)
      # compute C
      sc[m,1] = length(intersect(nam_sc,T_sett))
      # compute IC
      sc[m,2] = length(setdiff(nam_sc,T_sett))
      # compute C-fit
      sc[m,3] = ifelse(sc[m,1] == s&sc[m,2] == 0,1,0)
      sc_beta = rep(0,p)
      names(sc_beta) = paste("Z",1:p,sep  = ".")
      sc_beta[sapply(1:length(nam_sc),function(i){which(names(sc_beta) == nam_sc[i])})] = betasc
      # compute MMSE
      sc[m,4] = t(as.matrix(sc_beta-beta1))%*%Sig.X%*%as.matrix(sc_beta-beta1)
      
      filenames = paste("./simulation/Table4/",'n',n,'+','p',p,'+','Case',Caseset,'+CRSIS+SCAD.csv',sep = '')
      write.csv(data.frame(sc),filenames)
      
      ## SCAD End ##
      
      ## MCP Begin ##
      
      fitmc <- crrp(ftime, fstatus, cov, penalty = "MCP")
      betamc <-  fitmc$beta[, which.min(fitmc$BIC)]
      betamc = betamc[which(betamc!= 0)]
      nam_mc = names(betamc)
      # compute C
      mc[m,1] = length(intersect(nam_mc,T_sett))
      # compute IC
      mc[m,2] = length(setdiff(nam_mc,T_sett))
      # compute C-fit
      mc[m,3] = ifelse(mc[m,1] == s&mc[m,2] == 0,1,0)
      mc_beta = rep(0,p)
      names(mc_beta) = paste("Z",1:p,sep  = ".")
      mc_beta[sapply(1:length(nam_mc),function(i){which(names(mc_beta) == nam_mc[i])})] = betamc
      # compute MMSE
      mc[m,4] = t(as.matrix(mc_beta-beta1))%*%Sig.X%*%as.matrix(mc_beta-beta1)
      
      filenames = paste("./simulation/Table4/",'n',n,'+','p',p,'+','Case',Caseset,'+CRSIS+MCP.csv',sep = '')
      write.csv(data.frame(mc),filenames)
      
      ## MCP End ##
      
      ## Oracle Begin ##
      
      ftime <- X
      fstatus <- deleps
      cov1 = Z[,1:s]
      coeff = crr(ftime, fstatus, cov1)$coef
      or_beta = c(coeff,rep(0,(p-s)))
      # compute MMSE
      or[m,4] = t(as.matrix(or_beta-beta1))%*%Sig.X%*%as.matrix(or_beta-beta1)
      # get C, IC and C-fit
      or[m,1:3] <- c(s,0,1)
      
      filenames = paste("./simulation/Table4/",'n',n,'+','p',p,'+','Case',Caseset,'+CRSIS+Oracle.csv',sep = '')
      write.csv(data.frame(or),filenames)
      
      ## Oracle End ##
      
      ########### Compute and Save Average results  ########
      
      Average <- data.frame(
        LASSO = c(colMeans(na.omit(las[,1:3])),median(na.omit(las[,4]))),
        SCAD = c(colMeans(na.omit(sc[,1:3])),median(na.omit(sc[,4]))),
        MCP = c(colMeans(na.omit(mc[,1:3])),median(na.omit(mc[,4]))),
        Oracle = c(colMeans(na.omit(or[,1:3])),median(na.omit(or[,4]))),
        row.names = c("C","IC",'C-fit','MMSE'))
      filenames = paste("./simulation/Table4/",'n',n,'+','p',p,'+','Case',Caseset,'+CRSIS+Average.csv',sep = '')
      write.csv(Average,filenames)
      
      
      m = m+1},silent = TRUE)}
}

