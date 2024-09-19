#Enter "source('demo.R')" in the command line to run all the contents of the file

#Try to run in a higher version of R(4.02),but R-3.6 can also run

#need pkg Matrix and NMI
rm(list=ls())

# Example of SBM
#-------------------------------------------------------------------------------
#library(Matrix)
library(NMI)
source("SCP.R")
source("DCPPL_DCSBM.R")

# parameter setting
set.seed(123)
n <- 1200
pic <- c(0.2, 0.3, 0.5)
K <- length(pic)
Lambda <- (1*1e-2)*(matrix(1,K,K)+diag(c(2,3,4)))
m <- 4
Theta <- sample(c((m*2)/(m+1),2/(m+1)),size = n,replace = T,prob = c(0.5,0.5))

# data generation
Data_DCSBM <- Adj_Generating_DCSBM(n,pic,Lambda,Theta,sp = T,del_0d =T)

# run the method of SCP and DC_PPL
SCP_test <- SCP(Data_DCSBM$Adj,K)
DCPPL_test <- DCPPL(Data_DCSBM$Adj,K,SCP_test$C_Hat,MAX_ITER_Outer = 25,MAX_ITER_Inter = 50,sp = T) # this step takes about 16mins


# print the result
n <- length(Data_DCSBM$clusters)
X <- data.frame(1:n,Data_DCSBM$clusters)
Y_SCP <- data.frame(1:n,SCP_test$C_Hat)
Y_DCPPL <- data.frame(1:n,DCPPL_test$C_Hat)
print(c(NMI(X,Y_SCP)$value,NMI(X,Y_DCPPL)$value))
#result is 0.7137063 0.8039974
#-------------------------------------------------------------------------------

