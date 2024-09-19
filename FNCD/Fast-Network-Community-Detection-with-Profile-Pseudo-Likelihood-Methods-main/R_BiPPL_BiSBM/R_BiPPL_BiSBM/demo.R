#Enter "source('demo.R')" in the command line to run all the contents of the file

#Try to run in a higher version of R(4.02),but R-3.6 can also run

#need pkg Matrix and NMI
rm(list=ls())

# Example of SBM
#-------------------------------------------------------------------------------
#library(Matrix)
library(NMI)
source("SCP.R")
source("BiPPL_BiSBM.R")

# parameter setting
set.seed(123)
n_1 <- 1200
n_2 <- 1600
pi_1 <- c(0.5,0.5)
pi_2 <- c(0.5,0.5)
K_1 <- length(pi_1)
K_2 <- length(pi_2)
P <- matrix(c(0.16,0.2,0.2,0.16),2,2,byrow = T)

# data generation
Data_BiSBM <- Adj_Generating_BiSBM(n_1,n_2,pi_1,pi_2,P,sp = T)
Adj <- Data_BiSBM$Adj


# run the method of SCP and BiPPL
SCP_test_2 <- SCP(t(Adj)%*%Adj,K_2)
SCP_test_1 <- SCP(Adj%*%t(Adj),K_1)
BiPPL_test <- PPL_TwoSide_SBM(Adj,K_1,SCP_test_1$C_Hat,
                              K_2,SCP_test_2$C_Hat,
                              P_Init = NA,Pi_Init = NA,
                              MAX_ITER_Outer = 25,
                              MAX_ITER_Inter = 100,
                              DeltatMax = 0.01)

BiPPL_test_1 <- BiPPL_test$result_1
BiPPL_test_2 <- BiPPL_test$result_2


# print the result
n_1 <- length(Data_BiSBM$clusters_1)
n_2 <- length(Data_BiSBM$clusters_2)
X_1 <- data.frame(1:n_1,Data_BiSBM$clusters_1)
X_2 <- data.frame(1:n_2,Data_BiSBM$clusters_2)
Y_SCP_1 <- data.frame(1:n_1,SCP_test_1$C_Hat)
Y_SCP_2 <- data.frame(1:n_2,SCP_test_2$C_Hat)
Y_BiPPL_1 <- data.frame(1:n_1,BiPPL_test_1$C_Hat)
Y_BiPPL_2 <- data.frame(1:n_2,BiPPL_test_2$C_Hat)

print(c(NMI(X_1,Y_SCP_1)$value,NMI(X_1,Y_BiPPL_1)$value,
        NMI(X_2,Y_SCP_2)$value,NMI(X_2,Y_BiPPL_2)$value))
#result is 0.7584632 0.7935783 0.6978838 0.7491656
#-------------------------------------------------------------------------------


