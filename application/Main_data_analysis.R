source("./simulation/BasicFunctions.R")

######## Generate Table S1  Begin ##########

source("./application/Table5.R")

tables1 <- matrix(NA,nr = 58,nc = 6)
Method <- c("CRSIS","CRCSIS","CRSJS","crCRS","crSIRS","SIS")
for(i in 1:6){
  mm = Method[i]
  value <- read.csv(paste("./application/TableS1/",mm,"-genes.csv",sep = ""))[,-1]
  tables1[,i] <- value
}
tables1df <- data.frame(tables1)
colnames(tables1df) <- Method
write.csv(tables1df,"./application/TableS1/TableS1.csv")

######## Generate Table S1 End ##########

######## Generate Table 5 Begin ##########

source("./application/Table5.R")

Method <- c("CRSIS","CRCSIS","CRSJS","crSIRS")
Method2 <- c("LASSO","SCAD","MCP")

for(i in Method){
  t1 <- read.csv(paste("./application/Table5/",i,'_LASSOestimates.csv',sep = ''))
  colnames(t1) <- c('gene','LASSO')
  t2 <- read.csv(paste("./application/Table5/",i,'_SCADestimates.csv',sep = ''))
  colnames(t2) <- c('gene','SCAD')
  t3 <- read.csv(paste("./application/Table5/",i,'_MCPestimates.csv',sep = ''))
  colnames(t3) <- c('gene','MCP')
  table5 <- merge(merge(t1,t2,by = 'gene',all = TRUE),t3,by = 'gene',all = TRUE)
  write.csv(table5,paste("./application/Table5/Table5_",i,'.csv',sep = ''))
}

######## Generate Table 5 End ##########


######## Generate Table 6 Begin ##########

source("./application/Table6.R")

result <- matrix(NA,nc = 3,nr = 16)
R_screen <- c("CRSIS","CRCSIS","CRSJS","crSIRS")
R_method <- c("LASSO","SCAD","MCP","ALL")
for(i in 1:4){
  mm1 <- R_screen[i]
  for(j in 1:4){
    mm2 <- R_method[j]
    value <- read.csv(paste("./application/Table6/",mm1,"_",mm2,".csv",sep = ''))[,-1]
    result[4*(i-1)+j,] <- c(paste(mm1,'_',mm2,sep = ''),as.matrix(value))
  }
}
resultdf <- data.frame(result)
colnames(resultdf) <- c("Method","IBS","C-index")

write.csv(resultdf,"./application/Table6/Table6.csv")

######## Generate Table 6 End ##########

######## Plot Figure3 Begin ##########

source("./application/Figure3.R")

sirp1 <- read.csv("./application/Figure3/crSIRS+LASSO.csv")[,-1]
sirp2 <- read.csv("./application/Figure3/crSIRS+MCP.csv")[,-1]
sirp3 <- read.csv("./application/Figure3/crSIRS+SCAD.csv")[,-1]
sirp4 <- read.csv("./application/Figure3/CR-SIS+LASSO.csv")[,-1]
sirp5 <- read.csv("./application/Figure3/CR-SIS+MCP.csv")[,-1]
sirp6 <- read.csv("./application/Figure3/CR-SIS+SCAD.csv")[,-1]
sirp7 <- read.csv("./application/Figure3/CR-CSIS+LASSO.csv")[,-1]
sirp8 <- read.csv("./application/Figure3/CR-CSIS+MCP.csv")[,-1]
sirp9 <- read.csv("./application/Figure3/CR-CSIS+SCAD.csv")[,-1]
sirp10 <- read.csv("./application/Figure3/CR-SJS+LASSO.csv")[,-1]
sirp11 <- read.csv("./application/Figure3/CR-SJS+MCP.csv")[,-1]
sirp12 <- read.csv("./application/Figure3/CR-SJS+SCAD.csv")[,-1]

par(mgp = c(2,1,0))

setEPS()
postscript('./application/Figure3/LASSO.eps') 

plot(sirp1,main  =  "", xlab  =  "time (months)",ylab = "prediction CIF",xlim = c(0,80),ylim = c(0,0.15),lty = 2,type = 'l',col = "grey",lwd = 2)
lines(sirp4,lty = 2)
lines(sirp7,lty = 1)
lines(sirp10,lty = 1,col = "grey")
legend("topleft", legend = c("CR-SIS+LASSO", "CR-CSIS+LASSO","CR-SJS+LASSO","crSIRS+LASSO"),
       lty = c(2,1,1,2), col = c("black", "black","grey","grey"),cex = 0.7)
dev.off()

postscript('./application/Figure3/MCP.eps') 
plot(sirp2,main  =  "", xlab  =  "time (months)",ylab = "prediction CIF",xlim = c(0,80),ylim = c(0,0.15),type = 'l',lty = 2,col = "grey",lwd = 2)
lines(sirp5,lty = 2)
lines(sirp8,lty = 1)
lines(sirp11,lty = 1,col = "grey")
legend("topleft", legend = c("CR-SIS+MCP", "CR-CSIS+MCP","CR-SJS+MCP","crSIRS+MCP"),
       lty = c(2,1,1,2), col = c("black", "black","grey","grey"),cex = 0.7)
dev.off()

postscript('./application/Figure3/SCAD.eps') 
plot(sirp3,main  =  "", xlab  =  "time (months)",ylab = "prediction CIF",xlim = c(0,80),ylim = c(0,0.15),type = 'l',lty = 2,col = "grey",lwd = 2)
lines(sirp6,lty = 2)
lines(sirp9,lty = 1)
lines(sirp12,lty = 1,col = "grey")
legend("topleft", legend = c("CR-SIS+SCAD", "CR-CSIS+SCAD","CR-SJS+SCAD","crSIRS+SCAD"),
       lty = c(2,1,1,2), col = c("black", "black","grey","grey"),cex = 0.7)
dev.off()

######## Plot Figure3 End ##########
