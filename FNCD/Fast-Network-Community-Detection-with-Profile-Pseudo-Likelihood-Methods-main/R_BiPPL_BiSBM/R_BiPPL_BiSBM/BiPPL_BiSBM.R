Adj_Generating_BiSBM <- function(
  N_1,N_2,
  clusters_distribution_1,
  clusters_distribution_2,
  connection_matrix,
  sp = FALSE){
  
  num_clusters_1 <- length(clusters_distribution_1)
  label_clusters_1 <- 1:num_clusters_1
  clusters_1 <- sample(label_clusters_1,size=N_1,replace=T,
                       prob=clusters_distribution_1)
  
  num_clusters_2 <- length(clusters_distribution_2)
  label_clusters_2 <- 1:num_clusters_2
  clusters_2 <- sample(label_clusters_2,size=N_2,replace=T,
                       prob=clusters_distribution_2)
  
  RM <- matrix(runif(N_1*N_2),N_1,N_2)
  detM <- connection_matrix[clusters_1,clusters_2]
  Adj <- (RM < detM) + 0
  
  Ind_1 <- apply(Adj,1,sum)!=0
  Ind_2 <- apply(Adj,2,sum)!=0
  Adj <- Adj[Ind_1,Ind_2]
  clusters_1 <- clusters_1[Ind_1] 
  clusters_2 <- clusters_2[Ind_2]
  
  if(sp){
    library(Matrix)
    Adj <- Matrix(Adj,sparse=T)
  }
  
  re <- list(clusters_1,clusters_2,Adj)
  names(re) <- c('clusters_1','clusters_2','Adj')
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


PPL_OneSide_SBM <- function(
  Adj,K_1,Label_Init_1,K_2,Label_Init_2,
  P_Init = NA,Pi_Init = NA,
  MAX_ITER_Outer = 25,
  MAX_ITER_Inter = 100,
  DeltatMax = 0.01
){
  star_time<-Sys.time()
  
  Ind_1 <- apply(Adj,1,sum)==0
  Ind_2 <- apply(Adj,1,sum)==0
  Ind_0 <- Ind_1 | Ind_2
  if(any(Ind_0)){
    print('Zero degree node here')
    return(NA)
  }
  
  Num_Nodes_1 <- length(Adj[,1]) #compute the number of nodes
  Num_Nodes_2 <- length(Adj[1,])
  
  if (any(is.na(P_Init)) | any(is.na(Pi_Init))){
    #compute the Pi_Hat
    Pi_Hat <- numeric(K_1)
    for (i in 1:K_1){
      Pi_Hat[i] <- sum(Label_Init_1==i)/length(Label_Init_1)
    }
    
    #compute the P_Hat
    P_Hat <- matrix(1:(K_1*K_2),K_1,K_2)
    for (i in 1:K_1){
      for (j in 1:K_2){
        row_ind <- Label_Init_1==i
        col_ind <- Label_Init_2==j
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
  
  K_I <- diag(K_2)
  Q_2 <- K_I[Label_Init_2 , ] #initialize Q_2 = chat
  
  for(i in 1:MAX_ITER_Outer){
    #compute the matrix of the value of B-Function in article
    B_Matrix <- Adj%*%Q_2 
    
    #compute the vector of the value of N_2-Function in article
    N2_Matrix <- apply(Q_2,2,sum)
    
    Deltat <- Inf #record the new relative changes in likelihood
    Count_Inter <- 1 #the counter of the iterations times of 
    CONVERGED <- FALSE #the sign of convenience
    while((!CONVERGED)&&(Count_Inter<=MAX_ITER_Inter)){
      #compute the matrix of the value of F_1-Function in article
      Log_F1_Matrix <- B_Matrix%*%t(log(P_Hat+1e-30))+(ex_Ni(N2_Matrix,Num_Nodes_1)-B_Matrix)%*%t(log(1-P_Hat+1e-30))+ex_Ni(log(Pi_Hat),Num_Nodes_1)
      Log_F1_Matrix <- t(apply(Log_F1_Matrix,1,regl_m))
      F1_Matrix <- exp(Log_F1_Matrix)
      
      #update the matrix of the value of q_1-Function in article
      Q_1 <- t(apply(F1_Matrix,1,regl_s))
      
      #update the matrix of the value of N_1-Function in article
      N1_Matrix <- apply(Q_1,2,sum)
      
      #update the value of Pi_Hat by giving Q_1,Q_2
      Pi_Hat <- (N1_Matrix)/Num_Nodes_1
      
      #update the value of P_Hat by giving Q_1,Q_2
      P_Hat_up <- t(Q_1)%*%B_Matrix
      P_Hat_down <- N1_Matrix%*%t(N2_Matrix)
      P_Hat <- P_Hat_up/(P_Hat_down+1e-30)
      P_Hat[P_Hat>1] <- 1
      
      #compute the new value of likelihood
      LH_value_new <- sum((t(Q_1)%*%B_Matrix)*log(P_Hat+1e-30))+sum((t(Q_1)%*%(ex_Ni(N2_Matrix,Num_Nodes_1)-B_Matrix))*log(1-P_Hat+1e-30))+t(N1_Matrix)%*%log(Pi_Hat)
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
    #update the matrix of the value of q_1-Function in article
    #Q_2 <- Q_1
    
    #update Q_2
    B_Matrix_star <- t(Adj)%*%Q_1
    Q_2 <- B_Matrix_star%*%log(P_Hat+1e-30)+(ex_Ni(N1_Matrix,Num_Nodes_2)-B_Matrix_star)%*%log(1-P_Hat+1e-30)
    Q_2 <- apply(Q_2,1,which.max)
    Q_2 <- K_I[Q_2 , ]
    
    #record the number of iterations
    Iter[i+1] <- Count_Inter
    if (Iter[i] <= 3) { break }#determine whether to end
  }
  #compute the value of vector C_Hat 
  #which is the types and label of nudes
  C_Hat <- apply(Q_2,1,which.max)
  
  end_time<-Sys.time()
  timeused <- end_time-star_time
  
  #putout the result
  re <- list(C_Hat,P_Hat,Pi_Hat,timeused)
  names(re) <- c('C_Hat','P_Hat','Pi_Hat','time')
  return(re)
}

PPL_TwoSide_SBM <- function(
  Adj,K_1,Label_Init_1,K_2,Label_Init_2,
  P_Init = NA,Pi_Init = NA,
  MAX_ITER_Outer = 25,
  MAX_ITER_Inter = 100,
  DeltatMax = 0.01){
  
  re_2 <- PPL_OneSide_SBM(Adj,K_1,Label_Init_1,K_2,Label_Init_2,
                          P_Init,Pi_Init,
                          MAX_ITER_Outer,MAX_ITER_Inter,DeltatMax)
  re_1 <- PPL_OneSide_SBM(t(Adj),K_2,Label_Init_2,K_1,Label_Init_1,
                          P_Init,Pi_Init,
                          MAX_ITER_Outer,MAX_ITER_Inter,DeltatMax)
  re <- list(re_1,re_2)
  names(re) <- c('result_1','result_2')
  return(re)
}