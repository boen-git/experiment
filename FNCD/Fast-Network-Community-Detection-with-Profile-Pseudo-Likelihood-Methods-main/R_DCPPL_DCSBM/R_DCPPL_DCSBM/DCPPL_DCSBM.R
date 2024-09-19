Adj_Generating_DCSBM <- function(N,Clusters_distribution,
                                 Lambda,Theta,
                                 sp = FALSE,del_0d = FALSE){
  num_clusters <- length(Clusters_distribution)
  label_clusters <- 1:num_clusters
  clusters <- sample(label_clusters,size=N,replace=T,
                     prob=Clusters_distribution)
  
  K_I <- diag(num_clusters)
  C_Mat <- K_I[clusters, ]
  
  Theta_Mat <- diag(Theta)
  
  P_A <- Theta_Mat%*%(C_Mat%*%Lambda%*%t(C_Mat))%*%Theta_Mat
  
  RM <- matrix(runif(N*N),N,N)
  RM[lower.tri(RM)] <- t(RM)[lower.tri(RM)] #Symmetrization
  
  Adj <- (RM < P_A)+0
  
  if(del_0d){
    Ind_0 <- apply(Adj,1,sum)!=0
    Adj <- Adj[Ind_0,Ind_0]
    clusters <- clusters[Ind_0] 
  }
  
  if(sp){
    library(Matrix)
    Adj <- Matrix(Adj,sparse=T)
  }
  
  re <- list(clusters,Adj)
  names(re) <- c('clusters','Adj')
  return(re)
}


ex_Ni<-function(Ni,NN,sp = FALSE){
  if(sp){Ni <- as.matrix(Ni)}
  Ni_ex<-matrix(rep(Ni,NN),NN,length(Ni),byrow = T)
  return(Ni_ex)
}

regl_m<-function(x){
  y <- x-max(x)
  return(y)
}

regl_s<-function(x){
  y <- x/(sum(x)+1e-30)
  return(y)
}


DCPPL<- function(
  Adj,K,Label_Init,
  MAX_ITER_Outer = 25,
  MAX_ITER_Inter = 50,
  #MAX_ITER_FPI =20,
  sp = F
){
  if(sp){
    library(Matrix)
  }
  
  star_time<-Sys.time()
  
  Ind_0 <- apply(Adj,1,sum)!=0
  Adj <- Adj[Ind_0,Ind_0]
  Label_Init <- Label_Init[Ind_0]
  
  Num_Nodes <- length(Adj[1,]) #compute the number of nodes
  
  dg <- apply(Adj,2,sum)
  Theta <- Num_Nodes*dg/sum(dg)
  
  K_I <- diag(K)
  if(sp){
    Q_2 <- Matrix(K_I[Label_Init , ],sparse = sp)
  }else{Q_2 <- K_I[Label_Init , ]}#initialize Q_2 = e_matrix_old
  Q_1 <- Q_2
  
  for(i in 1:MAX_ITER_Outer){
    B_Matrix <- Adj%*%Q_2 
    for(j in 1:MAX_ITER_Inter){
      Pi_Hat <- apply(Q_1,2,mean)
      
      a_title <- t(Q_1)%*%B_Matrix
      b_title <- t(Q_1)%*%Theta%*%(t(Theta)%*%Q_2)
      Lambda <- a_title/b_title
      
      #Adj_star <- Q_1%*%Lambda%*%t(Q_2)
      #for (k in 1:MAX_ITER_FPI){
      #Theta <- 2*dg/((Adj_star+t(Adj_star))%*%Theta)[,1]
      #Theta <- Num_Nodes*Theta/sum(Theta)
      #}
      for(k in 1:Num_Nodes){
        c_dk <- -dg[k]
        a_mf <- (2*(t(Q_1[k,])%*%Lambda%*%Q_2[k,]))[1,1]
        indx <- 1:Num_Nodes
        indx <- indx[-k]
        deltai <- (t(Theta[indx])%*%Q_2[indx,]%*%Lambda)[1,]
        b_key <- (t(Q_1[k,])%*%deltai)[1,1]
        
        Theta[k] <- (-b_key + sqrt((b_key)^2 - 4*a_mf*c_dk))/(2*a_mf)
      }
      Theta <- Num_Nodes*Theta/sum(Theta)
      
      T_1 <- ex_Ni(log(Pi_Hat+1e-30),Num_Nodes,sp)
      if(sp){
        T_2 <- Matrix(diag(Theta),sparse = sp)
      }else{T_2 <- diag(Theta)}
      T_3 <- ex_Ni(t(Theta)%*%Q_2%*%Lambda,Num_Nodes,sp)
      T_4 <- B_Matrix%*%log(Lambda+1e-30)
      Z <- T_1 - T_2%*%T_3 + T_4
      Z <- t(apply(Z,1,regl_m))
      U <- exp(Z)
      
      Q_1 <- t(apply(U,1,regl_s))
    }
    #update Q_2 
    T_5 <- ex_Ni(t(Theta)%*%Q_1%*%Lambda,Num_Nodes,sp)
    Q_2 <- -T_2%*%T_5 + t(Adj)%*%Q_1%*%log(Lambda)
    Q_2 <- apply(Q_2,1,which.max)
    Q_2 <- K_I[Q_2 , ]
  }
  #compute the value of vector C_Hat 
  #which is the types and label of nudes
  C_Hat <- apply(Q_2,1,which.max)
  
  Pi_Hat <- numeric(K)
  for (i in 1:K){
    Pi_Hat[i] <- sum(C_Hat==i)/length(C_Hat)
  }
  
  C_Hat_title <- numeric(length(Ind_0))
  C_Hat_title[Ind_0] <- C_Hat
  C_Hat <- C_Hat_title
  
  C_Hat[!Ind_0] <- sample(1:K,sum(((!Ind_0)+0)),replace = T,prob = Pi_Hat)
  
  end_time<-Sys.time()
  timeused <- end_time-star_time
  
  #putout the result
  re <- list(C_Hat,Pi_Hat,timeused)
  names(re) <- c('C_Hat','Pi_Hat','time')
  return(re)
}