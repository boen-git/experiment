Adj_Generating_SBM <- function(
  N,
  clusters_distribution,
  connection_matrix,
  sp = FALSE,del_0d = TRUE){
  
  num_clusters <- length(clusters_distribution)
  label_clusters <- 1:num_clusters
  clusters <- sample(label_clusters,size=N,replace=T,
                     prob=clusters_distribution)
  
  RM <- matrix(runif(N*N),N,N)
  RM[lower.tri(RM)] <- t(RM)[lower.tri(RM)] #Symmetrization
  detM <- connection_matrix[clusters,clusters]
  Adj <- (RM < detM) + 0
  
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

PPL_SBM <- function(
  Adj,K,Label_Init,
  P_Init = NA,Pi_Init = NA,
  MAX_ITER_Outer = 25,
  MAX_ITER_Inter = 100,
  DeltatMax = 0.01,
  sp = FALSE
){
  if(sp){
    library(Matrix)
  }
  
  star_time<-Sys.time()
  
  Ind_0 <- apply(Adj,1,sum)!=0
  Adj <- Adj[Ind_0,Ind_0]
  Label_Init <- Label_Init[Ind_0]
  
  Num_Nodes <- length(Adj[1,]) #compute the number of nodes
  
  if (any(is.na(P_Init)) | any(is.na(Pi_Init))){
    #compute the Pi_Hat
    Pi_Hat <- numeric(K)
    for (i in 1:K){
      Pi_Hat[i] <- sum(Label_Init==i)/length(Label_Init)
    }
    
    #compute the P_Hat
    P_Hat <- diag(K)
    for (i in 1:K){
      for (j in 1:K){
        row_ind <- Label_Init==i
        col_ind <- Label_Init==j
        sub_Adj <- Adj[row_ind,col_ind]
        P_Hat[i,j] <- sum(sub_Adj)/(sum(col_ind+0)*sum(row_ind+0))
      }
    }
  }else{
    P_Hat <- P_Init #initialize P_Hat
    Pi_Hat <- Pi_Init #initialize Pi_Hat
  }
  
  #Allocate space to record the number of iterations of 
  #the optimization lower bound J-function in article
  #and serve as the basis for stopping
  Iter <- numeric(MAX_ITER_Outer+1)
  Iter[1] <- MAX_ITER_Inter
  
  #Allocate space to record the value of likelihood
  #LH_value <- matrix(numeric(MAX_ITER_Inter*MAX_ITER_Outer),MAX_ITER_Outer,MAX_ITER_Inter)
  
  K_I <- diag(K)
  #Q_2 <- K_I[Label_Init , ] #initialize Q_2 = chat
  if(sp){
    Q_2 <- Matrix(K_I[Label_Init , ],sparse = sp)
  }else{Q_2 <- K_I[Label_Init , ]}#initialize Q_2 = e_matrix_old
  
  for(i in 1:MAX_ITER_Outer){
    #compute the matrix of the value of B-Function
    B_Matrix <- Adj%*%Q_2 
    
    #compute the vector of the value of N_2-Function
    N2_Matrix <- apply(Q_2,2,sum)
    
    Deltat <- Inf #record the new relative changes in likelihood
    Count_Inter <- 1 #the counter of the iterations times of 
    CONVERGED <- FALSE #the sign of convenience
    while((!CONVERGED)&&(Count_Inter<=MAX_ITER_Inter)){
      #compute the matrix of the value of F_1-Function
      Log_F1_Matrix <- B_Matrix%*%t(log(P_Hat+1e-30))+(ex_Ni(N2_Matrix,Num_Nodes,sp)-B_Matrix)%*%t(log(1-P_Hat+1e-30))+ex_Ni(log(Pi_Hat),Num_Nodes,sp)
      Log_F1_Matrix <- t(apply(Log_F1_Matrix,1,regl_m))
      F1_Matrix <- exp(Log_F1_Matrix)
      
      #update the matrix of the value of q_1-Function
      Q_1 <- t(apply(F1_Matrix,1,regl_s))
      
      #update the matrix of the value of N_1-Function
      N1_Matrix <- apply(Q_1,2,sum)
      
      #update the value of Pi_Hat by giving Q_1,Q_2
      Pi_Hat <- (N1_Matrix)/Num_Nodes
      
      #update the value of P_Hat by giving Q_1,Q_2
      P_Hat_up <- t(Q_1)%*%B_Matrix
      P_Hat_down <- N1_Matrix%*%t(N2_Matrix)
      P_Hat <- P_Hat_up/(P_Hat_down+1e-30)
      P_Hat[P_Hat>1] <- 1
      
      #compute the new value of likelihood
      LH_value_new <- sum((t(Q_1)%*%B_Matrix)*log(P_Hat+1e-30))+sum((t(Q_1)%*%(ex_Ni(N2_Matrix,Num_Nodes,sp)-B_Matrix))*log(1-P_Hat+1e-30))+t(N1_Matrix)%*%log(Pi_Hat)
      if (Count_Inter > 1){
        #compute the new relative changes in likelihood
        Deltat <- (LH_value_new - LH_value_old)/(LH_value_old+1e-30)
        #determine whether to converge
        CONVERGED <- Deltat < DeltatMax
      }
      #record the relative changes in likelihood
      #LH_value[i,Count_Inter]
      #update the old value of likelihood
      LH_value_old <- LH_value_new
      
      Count_Inter <- Count_Inter+1
    }
    
    #update Q_2
    B_Matrix_star <- t(Adj)%*%Q_1
    Q_2 <- B_Matrix_star%*%log(P_Hat+1e-30)+(ex_Ni(N1_Matrix,Num_Nodes,sp)-B_Matrix_star)%*%log(1-P_Hat+1e-30)
    Q_2 <- apply(Q_2,1,which.max)
    Q_2 <- K_I[Q_2 , ]
    
    #record the number of iterations
    Iter[i+1] <- Count_Inter
    if (Iter[i] <= 3) { break }#determine whether to end
  }
  #compute the value of vector C_Hat 
  #which is the types and label of nudes
  C_Hat <- apply(Q_2,1,which.max)
  
  Pi_hHat <- numeric(K)
  for (i in 1:K){
    Pi_hHat[i] <- sum(C_Hat==i)/length(C_Hat)
  }
  
  C_Hat_title <- numeric(length(Ind_0))
  C_Hat_title[Ind_0] <- C_Hat
  C_Hat <- C_Hat_title
  
  C_Hat[!Ind_0] <- sample(1:K,sum(((!Ind_0)+0)),replace = T,prob = Pi_hHat)
  
  end_time<-Sys.time()
  timeused <- end_time-star_time
  
  #putout the result
  re <- list(C_Hat,P_Hat,Pi_Hat,timeused)
  names(re) <- c('C_Hat','P_Hat','Pi_Hat','time')
  return(re)
}
