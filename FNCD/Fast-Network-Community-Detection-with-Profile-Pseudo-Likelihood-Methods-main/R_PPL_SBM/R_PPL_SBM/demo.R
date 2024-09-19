#Enter "source('demo.R')" in the command line to run all the contents of the file

#Try to run in a higher version of R(4.02),but R-3.6 can also run

#need pkg Matrix and NMI
rm(list=ls())

#Example of SBM
#-------------------------------------------------------------------------------
#library(Matrix)
library(NMI)
source("SCP.R")
source("PPL_SBM.R")

#parameter setting
set.seed(123)
N_SBM <- 1000
Pi_SBM <- c(0.2,0.3,0.5)
K <- length(Pi_SBM)
P_in <- 0.012
P_bt <- 0.001
P_SBM <- (P_in-P_bt)*diag(K) + P_bt*matrix(1,K,K)


# data generation
Data_SBM <- Adj_Generating_SBM(N_SBM,Pi_SBM,P_SBM,sp=T,del_0d = T) 


# run the spectral method and PPL method
SCP_test <- SCP(Data_SBM$Adj,K)
PPL_test <- PPL_SBM(Data_SBM$Adj,K,SCP_test$C_Hat,MAX_ITER_Outer = 25,MAX_ITER_Inter = 100,DeltatMax = 0.01,sp = T)


# print the result of community detection
n <- length(Data_SBM$clusters)
X <- data.frame(1:n,Data_SBM$clusters)
Y_SCP <- data.frame(1:n,SCP_test$C_Hat)
Y_PPL <- data.frame(1:n,PPL_test$C_Hat)
print(c(NMI(X,Y_SCP)$value,NMI(X,Y_PPL)$value))
#result is 0.6992132 0.7629164
#-------------------------------------------------------------------------------

