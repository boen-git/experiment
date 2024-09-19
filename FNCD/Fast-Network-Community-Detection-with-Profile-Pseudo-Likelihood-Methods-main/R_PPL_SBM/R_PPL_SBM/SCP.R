
# the method of purmutated spectral clustering
SCP <- function(Adj,K){
  star_time<-Sys.time()
  n<-nrow(Adj)
  pho<-sum(Adj)/(n*(n-1))
  pAdj<-Adj+as.numeric(0.25*pho)*matrix(1,nrow=n,ncol=n) #purmutated Adj  
  #compute the Laplace Matrix
  D <- apply(pAdj,1,sum)^(-0.5)
  D[D==Inf] <- 0
  G <- diag(D)
  L <- G%*%pAdj%*%G
  
  #Take eigenvalues to reduce dimensionality
  U <- eigen(L,symmetric = T)
  U <- U$vectors
  U <- U[,1:K]
  
  C_Hat <- kmeans(U,K)$cluster
  
  #compute the Pi_Hat
  Pi_Hat <- numeric(K)
  for (i in 1:K){
    Pi_Hat[i] <- sum(C_Hat==i)/length(C_Hat)
  }
  
  #compute the P_Hat
  P_Hat <- diag(K)
  for (i in 1:K){
    for (j in i:K){
      row_ind <- C_Hat==i
      col_ind <- C_Hat==j
      sub_Adj <- Adj[row_ind,col_ind]
      P_Hat[i,j] <- sum(sub_Adj)/(sum(col_ind+0)*sum(row_ind+0))
      P_Hat[j,i] <- P_Hat[i,j]
    }
  }
  
  end_time<-Sys.time()
  timeused <- end_time-star_time
  
  #Putout
  re<-list(C_Hat,P_Hat,Pi_Hat,timeused)
  names(re) <- c('C_Hat','P_Hat','Pi_Hat','time')
  return(re)
}
