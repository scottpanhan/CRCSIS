## Required functions ###

source("./simulation/BasicFunctions.R")
source("./simulation/Table1.R")
source("./simulation/Table2.R")
source("./simulation/Table2(2).R")
source("./simulation/Table3.R")
source("./simulation/post_CRSIS.R")
source("./simulation/post_CRCSIS.R")
source("./simulation/FNR.R")


### Simulate Table 1 Begin ####

###################################################
# Parameters to specify #
###################################################
# sim : the number of replications
# n: sample size 
# pset: dimension of covariates
# Rhosets: coefficient of covariates
# CensRatesets: setup for censoring rate
# Boundsets: bounded censoring time or not

n = 200
sim = 200 
pset = 5000

Rhosets = c(0.2,0.8)
CensRatesets = c(20,40,60)
Boundsets = c(1,0)

for(Boundset in Boundsets){
  for(Rhoset in Rhosets){
    for(CensRateset in CensRatesets){
      set.seed(999)
      Screening1(sim,n,pset,Rhoset,CensRateset,Boundset)
    }
  }
}

## Load results and generate Table1.csv
Tab1 = data.frame()
for(Rhoset in Rhosets){
  for(Boundset in Boundsets){
    for(CensRateset in CensRatesets){
      filenames = paste("./simulation/Table1/",'p',pset,'+','rho',Rhoset*10,'+',
                        'C',CensRateset,'+','B',Boundset,"+",'Average','.csv',sep = '')
      res = cbind(Boundset,Rhoset,CensRateset,t(read.csv(filenames)[,-1]))
      colnames(res) = c("B","rho","cr","P1","P2","P3",'P4','P5','Pall')
      Tab1 = rbind(Tab1,res)
    }
  }
}

Table1 <- cbind(Tab1[1:36,],Tab1[37:72,])
write.csv(Table1,"./simulation/Table1.csv")

### Simulate Table 1 End ####




## Simulate Table 2 part 1 Begin ###

###################################################
# Parameters to specify #
###################################################
# sim : the number of replications
# n: sample size 
# pset: dimension of covariates
# Rhosets: coefficient of covariates
# CensRatesets: setup for censoring rate


n = 200
sim = 200 
pset = 5000

Rhosets = c(0.2,0.8)
CensRatesets = c(20,40,60)

for(Rhoset in Rhosets){
  for(CensRateset in CensRatesets){
    set.seed(999)
    Screening2(sim,n,pset,Rhoset,CensRateset)
  }
}


## Load results and generate Table2_1.csv
Tab2 <- data.frame()
for(Rhoset in Rhosets){
  for(CensRateset in CensRatesets){
    filenames = paste("./simulation/Table2/",'p',pset,'+','rho',Rhoset*10,'+',
                      'C',CensRateset,"+",'Average','.csv',sep = '')
    res = cbind(Rhoset,CensRateset,t(read.csv(filenames)[,-1]))
    colnames(res) = c("rho","cr","P1","P2","P3",'P4','P5','Pall')
    Tab2 = rbind(Tab2,res)
  }
}
Tab2 <- cbind(Tab2[1:18,],Tab2[19:36,])
write.csv(Tab2,"./simulation/Table2_1.csv")

## Simulate Table 2 part 1 End ###




## Simulate Table 2 part 2 Begin ###


###################################################
# Parameters to specify #
###################################################
# sim : the number of replications
# n: sample size 
# pset: dimension of covariates
# CensRatesets: setup for censoring rate
# Errorsets: Error follows from standard normal or standard Gumbel distributions
sim = 200 
n = 200
pset = 5000

CensRatesets = c(20,40,60)
Errorsets = c(1,0)

for(Error in Errorsets){
  for(CensRateset in CensRatesets){
    set.seed(999)
    Screening22(sim,n,pset,CensRateset,Error)
  }
}

## Load results and generate Table2_2.csv
Tab22 <- data.frame()
for(Error in Errorsets){
  for(CensRateset in CensRatesets){
    filenames = paste("./simulation/Table2(2)/",'p',pset,'+',
                      'C',CensRateset,'+','Norm',Error,"+",'Average','.csv',sep = '')
    res = cbind(Rhoset,CensRateset,t(read.csv(filenames)[,-1]))
    colnames(res) = c("Error","cr","P1","P2","P3",'P4','P5','Pall')
    Tab22 = rbind(Tab22,res)
  }
}
Tab22 <- cbind(Tab22[1:18,],Tab22[19:36,])

write.csv(Tab22,"./simulation/Table2_2.csv")
## Simulate Table 2 part 2 End ###




## Simulate Table 3 Begin ###
###################################################
# Parameters to specify #
###################################################
# sim : the number of replications
# n: sample size 
# pset: dimension of covariates
# Caseset: indicator for setup for 4 cases in example 5.4

sim = 200
Cases = c(1,2,3,4,5,6)

n = 200
pset = 500

for(Caseset in Cases){
  Screening3(sim,n,pset,Caseset)
}

n = 200
pset = 3000
for(Caseset in Cases){
  Screening3(sim,n,pset,Caseset)
}


## Simulate Table 3 End ###


## Load results and generate Table3.csv
Cases = c(1,2,3,4,5,6)
ss = c(6,6,6,6,5,4)
n=200
ps = c(500,3000)

Table3 = data.frame()
Rowname = c('P1','P2','P3','P4','P5','P6','Pall')

for(i in 1:2){
  for(j in 1:6){
    s = ss[j]
    Caseset = Cases[j]
    p = ps[i]
    filenames = paste("./simulation/Table3/",'n',n,'+','p',p,'+','Case',Caseset,"+",'Average','.csv',sep = '')
    Results = t(read.csv(filenames)[,-1])
    Results = round(Results,2)
    if(Caseset %in% 1:4){
      res = cbind(p,Caseset,s,Results)
      colnames(res) = c('p','case','s',Rowname)
      Table3 = rbind(Table3,res)
    }else if(Caseset == 5){
      Resultsnew = cbind(Results[,1:5],'-',Results[,6])
      res = cbind(p,Caseset,s,Resultsnew)
      colnames(res) = c('p','case','s',Rowname)
      Table3 = rbind(Table3,res)
    }else if(Caseset == 6){
      Resultsnew = cbind(Results[,1:4],'-','-',Results[,5])
      res = cbind(p,Caseset,s,Resultsnew)
      colnames(res) = c('p','case','s',Rowname)
      Table3 = rbind(Table3,res)
    }
  }
}
Table3 <- cbind(Table3[1:36,],Table3[37:72,])
write.csv(Table3,'./simulation/Table3.csv')


## Simulate Table 4 Begin ###
###################################################
# Parameters to specify #
###################################################
# sim : the number of replications
# n: sample size 
# pset: dimension of covariates
# Caseset: indicator for setup for 4 cases in example 5.4

n = 200
pset = 3000
sim = 200

Cases = c(1,2,3,4,5,6)


# generate post-screening resultS for CR-SIS
for(Caseset in Cases){
  PostCRSIS(sim,n,pset,Caseset)
}

# generate post-screening resultS for CR-cSIS
for(Caseset in Cases){
  PostCRCSIS(sim,n,pset,Caseset)
}



## Simulate Table 4 End ###
## Load results and generate Table4.csv
Cases = c(1,2,3,4,5,6)

n = 200

postSIS = data.frame()
postCSIS = data.frame()
for(p in c(500,3000)){
  for(Caseset in Cases){
    filename1 = paste("./simulation/Table4/",'n',n,'+','p',p,'+','Case',Caseset,'+CRSIS+Average.csv',sep = '')
    filename2 = paste("./simulation/Table4/",'n',n,'+','p',p,'+','Case',Caseset,'+CRCSIS+Average.csv',sep = '')
    res1 <- cbind(p,'CRSIS',Caseset,t(read.csv(filename1)[,-1]))
    res2 <- cbind(p,'CRCSIS',Caseset,t(read.csv(filename2)[,-1]))
    postSIS = rbind(postSIS,res1)
    postCSIS = rbind(postCSIS,res2)
  }
}
colnames(postSIS) = c('p',"method",'case',"C","IC",'C-fit','MMSE')
colnames(postCSIS) = c('p',"method",'case',"C","IC",'C-fit','MMSE')

Table4 <- cbind(postSIS,postCSIS)

write.csv(Table4,'./simulation/Table4.csv')



