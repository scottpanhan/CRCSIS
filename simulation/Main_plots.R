

### Plot Figure 1 Begin  #####

### Generate Example5.1 result for Figure1 #######

###################################################
# Parameters specified #
###################################################
Rhosets = c(0.2,0.8)
CensRatesets = c(20,40,60)
Boundsets = c(1,0)
n = 200
pset = 5000

Pa = data.frame()

# Load Example 5.1 results 
for(Boundset in Boundsets){
  for(Rhoset in Rhosets){
    for(CensRateset in CensRatesets){
      filenames = paste("./simulation/Table1/",'p',pset,'+','rho',Rhoset*10,'+',
                        'C',CensRateset,'+','B',Boundset,"+",'Average','.csv',sep = '')
      Average = read.csv(filenames,header = TRUE)
      r = which(Average$X == 'Pa')
      Pa = rbind(Pa,Average[r,])
    }
  }
}

write.csv(Pa,'./simulation/Figure1/Example5_1.csv')

### Generate Example5.1 result End #######

### Generate Example5.2 results for Figure 1 #######
###################################################
# Parameters specified #
###################################################
n = 200
pset = 5000

Rhosets = c(0.2,0.8)
CensRatesets = c(20,40,60)

# Load Example5.2 Results
Pa = data.frame()
for(Rhoset in Rhosets){
  for(CensRateset in CensRatesets){
    filenames = paste("./simulation/Table2/",'p',pset,'+','rho',Rhoset*10,'+',
                      'C',CensRateset,"+",'Average','.csv',sep = '')
    Average = read.csv(filenames,header = TRUE)
    r = which(Average$X == 'Pa')
    Pa = rbind(Pa,Average[r,])
  }
}

write.csv(Pa,'./simulation/Figure1/Example5_2.csv')

### Generate Example5.2 results  End #######

### Generate Example5.3 results  Begin #######
###################################################
# Parameters specified #
###################################################
n = 200
pset = 5000

CensRatesets = c(20,40,60)
Errorsets = c(1,0)

Pa = data.frame()

# Load  Example5.3 results
for(Error in Errorsets){
  for(CensRateset in CensRatesets){
    filenames = paste("./simulation/Table2(2)/",'p',pset,'+',
                      'C',CensRateset,'+','Norm',Error,"+",'Average','.csv',sep = '')
    Average = read.csv(filenames,header = TRUE)
    r = which(Average$X == 'Pa')
    Pa = rbind(Pa,Average[r,])
  }
}


write.csv(Pa,'./simulation/Figure1/Example5_3.csv')

### Generate Example5.3 results  End #######

### Generate Example5.4 results  Begin #######
###################################################
# Parameters specified #
###################################################
Cases = c(1,2,3,4,5,6)
nps = c(200,200,500,3000)

Pa = data.frame()
# Load Example5.4 results
for(i in 1:2){
  n = nps[i]
  p = nps[i+2]
  for(Caseset in Cases){
    filenames = paste("./simulation/Table3/",'n',n,'+','p',p,'+','Case',Caseset,"+",'Average','.csv',sep = '')
    Average = read.csv(filenames,header = TRUE)
    r = which(Average$X == 'Pa')
    Pa = rbind(Pa,Average[r,])
  }
}
write.csv(Pa,'./simulation/Figure1/Example5_4.csv')

### Generate Example5.4 results  End #######

### Plot Figure 1 Begin  #####

par(mgp = c(2,1,0))

# load results
a1 = read.csv('./simulation/Figure1/Example5_1.csv',header = TRUE)[,-(1:2)]
setEPS()
postscript('./simulation/Figure1/box_a.eps',width = 6, height = 4) 
boxplot(a1,names = c('CR-SIS','CR-CSIS','CR-SJS','crCRS','crSIRS','SIS'),outlty  =  3,outpch = 20,main = '(a)',boxwex = 0.7,whisklty  =  3)

dev.off()
# load results
a2 = read.csv('./simulation/Figure1/Example5_2.csv',head = TRUE)[,-(1:2)]

postscript('./simulation/Figure1/box_b.eps',width = 6, height = 4) 
boxplot(a2,names = c('CR-SIS','CR-CSIS','CR-SJS','crCRS','crSIRS','SIS'),outlty  =  3,outpch = 20,main = '(b)',boxwex = 0.7,whisklty  =  3)

dev.off()

# load results
a3 = read.csv('./simulation/Figure1/Example5_3.csv',head = TRUE)[,-(1:2)]
postscript('./simulation/Figure1/box_c.eps',width = 6, height = 4) 

boxplot(a3,names = c('CR-SIS','CR-CSIS','CR-SJS','crCRS','crSIRS','SIS'),outlty  =  3,outpch = 20,main = '(c)',boxwex = 0.7,whisklty  =  3)

dev.off()
# load results
a4 = read.csv('./simulation/Figure1/Example5_4.csv',head = TRUE)[,-(1:2)]

postscript('./simulation/Figure1/box_d.eps',width = 6, height = 4) 
boxplot(a4,names = c('CR-SIS','CR-CSIS','CR-SJS','crCRS','crSIRS','SIS'),outlty  =  3,outpch = 20,main = '(d)',boxwex = 0.7,whisklty  =  3)

dev.off()

### Plot Figure 1 End  #####







### For Figure 2, we firstly need to simulate Example 5.5. #####


## Simulate Example 5.5 Begin ###

###################################################
# Parameters specified #
###################################################
n = 200
sim = 200 
pset = 5000

Rhosets = c(0.2,0.8)
CensRatesets = c(20,40,60)
Boundsets = c(1,0)

q_gam = c(1/1000,seq(0.01,0.2,0.02),0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)


# Run simulations
Results <- data.frame(row.names  =  q_gam)
for(Boundset in Boundsets){
  for(Rhoset in Rhosets){
    for(CensRateset in CensRatesets){
      set.seed(999)
      Results1 <- FNR(sim,n,pset,Rhoset,CensRateset,Boundset,q_gam)
      Results = cbind(Results,Results1)
      write.csv(Results,"./simulation/Figure2/FNR.csv")
    }
  }
}
### Simulate Example 5.5 End ###

### Generate Example5.5 result for Figure 2#######
###################################################
# Parameters specified #
###################################################
n = 200
sim = 200 
pset = 5000

Rhosets = c(0.2,0.8)
CensRatesets = c(20,40,60)
Boundsets = c(1,0)

q_gam = c(1/1000,seq(0.01,0.2,0.02),0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)

# combine FNR results for example 5.5
Results <- data.frame(row.names  =  q_gam)
for(Boundset in Boundsets){
  for(Rhoset in Rhosets){
    for(CensRateset in CensRatesets){
      set.seed(999)
      R_SIS <- read.csv(paste("./simulation/Figure2/FNR",'p',pset,'+','rho',Rhoset*10,'+',
                              'C',CensRateset,'+','B',Boundset,"+",'CR-SIS+Average','.csv',sep = ''))
      R_CSIS <- read.csv(paste("./simulation/Figure2/FNR",'p',pset,'+','rho',Rhoset*10,'+',
                               'C',CensRateset,'+','B',Boundset,"+",'CR-CSIS+Average','.csv',sep = ''))
      Results = cbind(Results,data.frame(CRSIS = R_SIS[,2],CRCSIS = R_CSIS[,2]))
      write.csv(Results,"./simulation/Figure2/FNR.csv")
    }
  }
}
### Generate Example5.5 result end #######

## Plot Figure 2 Begin ###
# load results
FNR = read.csv("./simulation/Figure2/FNR.csv",header  =  TRUE)
par(pin=c(5,5),mgp = c(2,1,0))

setEPS()
postscript('./simulation/Figure2/b_01.eps') 
plot(FNR[,1],FNR[,2],ylim = c(0,1),xlim = c(0,1),ylab  = "False negative rate",xlab  =  expression(r[n]),main = expression(paste(paste(rho," = 0.2"),", 20% censoring, Bounded")), cex.main = 2, cex.axis = 2, cex.lab = 2, type  =  "l",lty = 2,lwd = 1.5)
lines(FNR[,1],FNR[,3],type = "l",lty = 1,lwd = 1.5)
legend("topright", legend  =  c("CR-SIS", "CR-CSIS"),
       lty  =  c(2,1))
dev.off()

postscript('./simulation/Figure2/b_02.eps') 
plot(FNR[,1],FNR[,4],ylim = c(0,1),xlim = c(0,1),ylab  = "False negative rate",xlab  =  expression(r[n]),main = expression(paste(paste(rho," = 0.2"),", 40% censoring, Bounded")), cex.main = 2, cex.axis = 2, cex.lab = 2, type  =  "l",lty = 2,lwd = 1.5)
lines(FNR[,1],FNR[,5],type = "l",lty = 1,lwd = 1.5)
dev.off()

postscript('./simulation/Figure2/b_03.eps') 
plot(FNR[,1],FNR[,6],ylim = c(0,1),xlim = c(0,1),ylab  = "False negative rate",xlab  =  expression(r[n]),main = expression(paste(paste(rho," = 0.2"),", 60% censoring, Bounded")), cex.main = 2, cex.axis = 2, cex.lab = 2, type  =  "l",lty = 2,lwd = 1.5)
lines(FNR[,1],FNR[,7],type = "l",lty = 1,lwd = 1.5)
dev.off()

postscript('./simulation/Figure2/b_04.eps') 
plot(FNR[,1],FNR[,8],ylim = c(0,1),xlim = c(0,1),ylab  = "False negative rate",xlab  =  expression(r[n]),main = expression(paste(paste(rho," = 0.8"),", 20% censoring, Bounded")), cex.main = 2, cex.axis = 2, cex.lab = 2, type  =  "l",lty = 2,lwd = 1.5)
lines(FNR[,1],FNR[,9],type = "l",lty = 1,lwd = 1.5)
dev.off()

postscript('./simulation/Figure2/b_05.eps') 
plot(FNR[,1],FNR[,10],ylim = c(0,1),xlim = c(0,1),ylab  = "False negative rate",xlab  =  expression(r[n]),main = expression(paste(paste(rho," = 0.8"),", 40% censoring, Bounded")), cex.main = 2, cex.axis = 2, cex.lab = 2, type  =  "l",lty = 2,lwd = 1.5)
lines(FNR[,1],FNR[,11],type = "l",lty = 1,lwd = 1.5)
dev.off()

postscript('./simulation/Figure2/b_06.eps') 
plot(FNR[,1],FNR[,12],ylim = c(0,1),xlim = c(0,1),ylab  = "False negative rate",xlab  =  expression(r[n]),main = expression(paste(paste(rho," = 0.8"),", 60% censoring, Bounded")), cex.main = 2, cex.axis = 2, cex.lab = 2, type  =  "l",lty = 2,lwd = 1.5)
lines(FNR[,1],FNR[,13],type = "l",lty = 1,lwd = 1.5)
dev.off()

postscript('./simulation/Figure2/ub_01.eps') 
plot(FNR[,1],FNR[,14],ylim = c(0,1),xlim = c(0,1),ylab  = "False negative rate",xlab  =  expression(r[n]),main = expression(paste(paste(rho," = 0.2"),", 20% censoring, Unbounded")), cex.main = 2, cex.axis = 2, cex.lab = 2, type  =  "l",lty = 2,lwd = 1.5)
lines(FNR[,1],FNR[,15],type = "l",lty = 1,lwd = 1.5)
dev.off()

postscript('./simulation/Figure2/ub_02.eps') 
plot(FNR[,1],FNR[,16],ylim = c(0,1),xlim = c(0,1),ylab  = "False negative rate",xlab  =  expression(r[n]),main = expression(paste(paste(rho," = 0.2"),", 40% censoring, Unbounded")), cex.main = 2, cex.axis = 2, cex.lab = 2, type  =  "l",lty = 2,lwd = 1.5)
lines(FNR[,1],FNR[,17],type = "l",lty = 1,lwd = 1.5)
dev.off()

postscript('./simulation/Figure2/ub_03.eps') 
plot(FNR[,1],FNR[,18],ylim = c(0,1),xlim = c(0,1),ylab  = "False negative rate",xlab  =  expression(r[n]),main = expression(paste(paste(rho," = 0.2"),", 60% censoring, Unbounded")), cex.main = 2, cex.axis = 2, cex.lab = 2, type  =  "l",lty = 2,lwd = 1.5)
lines(FNR[,1],FNR[,19],type = "l",lty = 1,lwd = 1.5)
dev.off()

postscript('./simulation/Figure2/ub_04.eps') 
plot(FNR[,1],FNR[,20],ylim = c(0,1),xlim = c(0,1),ylab  = "False negative rate",xlab  =  expression(r[n]),main = expression(paste(paste(rho," = 0.8"),", 20% censoring, Unbounded")), cex.main = 2, cex.axis = 2, cex.lab = 2, type  =  "l",lty = 2,lwd = 1.5)
lines(FNR[,1],FNR[,21],type = "l",lty = 1,lwd = 1.5)
dev.off()

postscript('./simulation/Figure2/ub_05.eps') 
plot(FNR[,1],FNR[,22],ylim = c(0,1),xlim = c(0,1),ylab  = "False negative rate",xlab  =  expression(r[n]),main = expression(paste(paste(rho," = 0.8"),", 40% censoring, Unbounded")), cex.main = 2, cex.axis = 2, cex.lab = 2, type  =  "l",lty = 2,lwd = 1.5)
lines(FNR[,1],FNR[,23],type = "l",lty = 1,lwd = 1.5)
dev.off()

postscript('./simulation/Figure2/ub_06.eps') 
plot(FNR[,1],FNR[,24],ylim = c(0,1),xlim = c(0,1),ylab  = "False negative rate",xlab  =  expression(r[n]),main = expression(paste(paste(rho," = 0.8"),", 60% censoring, Unbounded")), cex.main = 2, cex.axis = 2, cex.lab = 2, type  =  "l",lty = 2,lwd = 1.5)
lines(FNR[,1],FNR[,25],type = "l",lty = 1,lwd = 1.5)

dev.off()

### Plot Figure 2 End  #####