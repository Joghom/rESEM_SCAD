#### Title: Regularized ESEM method
#### Author: Jochem Breukers, Tra Le
#### Supervisor: Dr. Katrijn Van Deun 
#### Created: May 10, 2026
#### Last modified: June 22, 2026

#########################################################################################
#####                  rESEM-scad function for correlated factors                   #####
#########################################################################################
LOSS_SCAD <- function(DATA, SCORES, LOADINGS, LAMBDA, ASCAD, penalty){
  XHAT <- SCORES%*%t(LOADINGS)
  res <- sum(rowSums((XHAT-DATA)^2))
  if (penalty) {
    penalty_val <- numeric()
    for (j in 1:dim(LOADINGS)[1]){
      for (r in 1:dim(LOADINGS)[2]){
        penatlyLaod <- SCADPenalty(LOADINGS[j,r], LAMBDA, ASCAD)
        penalty_val <- c(penalty_val, penatlyLaod)
      }
    }
  }
  loss <- res+ sum(penalty_val)
  return(loss)
}

SCADPenalty <- function(LOADING, LAMBDA, ASCAD){
  if (abs(LOADING)<=LAMBDA) {
    penalty <- LAMBDA * abs(LOADING)
  }else if (abs(LOADING)>(ASCAD*LAMBDA)){
    penalty <- ((ASCAD+1) * LAMBDA^2)/2
  } else {
    penalty <- -((abs(LOADING)^2 - 2*ASCAD*LAMBDA*abs(LOADING) + LAMBDA^2)/
                   (2*(ASCAD-1)))
  }
  return(penalty)
}

#2. REGULARIZED ESEM WITH THE SCAD penalty
rESEM_SCAD <- function(DATA, R, lambda, ASCAD, RHO, MaxIter, eps, penalty){
  # initialize variables
  converged <- FALSE
  I <- dim(DATA)[1]
  J <- dim(DATA)[2]
  # 
  ssx <- sum(DATA^2)
  convAO <- 0
  iter <- 1
  Lossc <- 1
  Lossvec <- Lossc   #Used to compute convergence criterium
  svd1 <- svd(DATA, R, R)
  P2 <- svd1$v %*% diag(svd1$d[1:R])/sqrt(I)
  P1 <- matrix(rnorm(J*R), nrow = J, ncol = R) #initialize P 
  loadings <- P1*0.2+ P2*0.8
  scores <- svd1$u 
  diffT <- 0
  diffP <- 0
  diffQ <- 0
  
  #initialize matrices for ADMM
  Z <- t(loadings)
  U <- matrix(0, dim(Z)[1], dim(Z)[2])
  
  
  while (convAO == 0) {
    iter0 <- 1
    Losst <- 1
    Lossvec0 <- Losst
    convT0 <- 0
    Lossvec1 <- 1
    #1. Update component scores 
    while(convT0 == 0){
      Lossu1old <- LOSS_SCAD(DATA,scores,loadings,lambda, ASCAD, penalty)
      E <- DATA - scores%*%t(loadings) 
      for (r in 1:R){
        Er <- E + scores[,r]%*%t(loadings[,r])
        num <- Er%*%loadings[,r]
        scores[,r] <- sqrt(I)*num/sqrt(sum(num^2))
        Lossu1 <- LOSS_SCAD(DATA,scores,loadings,lambda, ASCAD, penalty)
        diffT <- c(diffT,Lossu1old-Lossu1)
        Lossu1old <- Lossu1
        Lossvec1 <- c(Lossvec1, Lossu1)
        
      }
      #t(scores)%*%scores
      #Calculate loss
      Lossu0 <- LOSS_SCAD(DATA,scores,loadings,lambda, ASCAD, penalty)/ssx
      Lossvec0 <- c(Lossvec0,Lossu0)
      
      # check convergence
      if (iter0 > MaxIter) {
        convT0 <- 1
      }
      if (abs(Losst-Lossu0) < eps){
        convT0 <- 1
      }
      iter0 <- iter0 + 1
      Losst <- Lossu0
      
    }
    Loss <- LOSS_SCAD(DATA, scores, loadings, lambda, ASCAD, penalty)/ssx
    
    
    #2. Update loadings
    iter1 <- 0
    Losst <- 1
    Lossvec0l <- Losst
    convT1 <- 0
    while(convT1 == 0){
      # apply ADMM method
      loadings <- t(2*t(scores)%*%DATA - RHO*(U-Z))%*%solve(2*t(scores)%*%scores+ diag(R)*RHO)
      
      for (j in 1:J){
        for (r in 1:R){
          x <- t(loadings)[r,j] + U[r,j]
          if (abs(x) <= (2*lambda)) {
            Z[r,j] <- sign(x)*max((abs(x)-lambda/RHO), 0)
          }else if (abs(x) > (ASCAD*lambda)) {
            Z[r,j] <- x
          }else {
            Z[r,j] <- ((sign(x)*ASCAD*(lambda/RHO)) + (ASCAD -1)*x)/(ASCAD-1-(1/RHO))
          }
        }
      }
      
      U = U + t(loadings) - Z
      
      
      
      Lossu0 <- LOSS_SCAD(DATA,scores,loadings,lambda, ASCAD, penalty)/ssx
      Lossvec0l <- c(Lossvec0,Lossu0)
      
      #calculate losses
      if (iter1 > MaxIter) {
        convT1 <- 1
      }
      if (abs(Losst-Lossu0) < eps){
        convT1 <- 1
      }
      iter1 <- iter1 + 1
      Losst <- Lossu0
      
    }
    
    Loss <- LOSS_SCAD(DATA, scores, loadings, lambda, ASCAD, penalty)/ssx
    
    #Calculate loss
    Lossu <- LOSS_SCAD(DATA,scores,loadings,lambda, ASCAD, penalty)/ssx
    Lossvec <- c(Lossvec,Lossu)
    if (iter > MaxIter) {
      convAO <- 1
    }
    #if (Lossc-Lossu < -1e-12) {
    #  warning('Increase in Loss')
    #  break
    #}
    if (abs(Lossc-Lossu) < eps) {
      convAO <- 1
    }
    
    iter <- iter + 1
    Lossc <- Lossu
  }
  if (iter < MaxIter) {
    converged <- TRUE
  }
  loadings <- t(Z)
  #print('loadings, Z')
  #print(t(loadings))
  #print(Z)
  
  return_rlslv <- list()
  return_rlslv$loadings <- loadings
  return_rlslv$scores <- scores
  return_rlslv$Loss <- Lossu
  return_rlslv$Lossvec <- Lossvec
  return_rlslv$converged <- converged
  
  return(return_rlslv)
}

#3. MULTISTART PROCEDURE

MULTISTART_rESEM_SCAD <- function(DATA, R, MaxIter, eps, nstarts, lambda, ASCAD, RHO, penalty){
  if(missing(nstarts)){
    nstarts <- 20
  } 
  if(missing(penalty)){
    penalty <- TRUE
  } 
  
  Pout3d <- list()
  Tout3d <- list()
  LOSS <- array()
  LOSSvec <- list()
  converged <- array()
  
  for (n in 1:nstarts){
    result <- rESEM_SCAD(DATA, R, lambda, ASCAD, RHO, MaxIter, eps, penalty)
    
    Pout3d[[n]] <- result$loadings
    Tout3d[[n]] <- result$scores
    LOSS[n] <- result$Loss
    LOSSvec[[n]] <- result$Lossvec
    converged[n] <- result$converged
  }
  
  # check how many times the minimum loss was achieved
  best <- min(LOSS)
  tol <- 1e-2
  n_best <- sum(abs(LOSS - best) < tol)
  n_distinct <- length(unique(round(LOSS, 2)))
  
  # choose solution with lowest loss value
  k <- which(LOSS == min(LOSS))
  if (length(k)>1){
    pos <- sample(1:length(k), 1)
    k <- k[pos]
  }
  
  return_varselect <- list()
  return_varselect$loadings <- Pout3d[[k]]
  return_varselect$scores <- Tout3d[[k]]
  return_varselect$Lossvec <- LOSSvec
  return_varselect$Loss <- LOSS[k]
  return_varselect$all_losses <- LOSS
  return_varselect$n_best <- n_best
  return_varselect$n_distinct <- n_distinct
  return_varselect$converged <- converged[k]
  
  
  return(return_varselect)
}

###function for recovery rate
num_correct <- function (TargetP, EstimatedP){
  total_vnumber <- dim(TargetP)[1] * dim(TargetP)[2]
  TargetP[which(TargetP != 0)] <- 1
  sum_select <- sum(TargetP)
  sum_zero <- total_vnumber - sum_select
  EstimatedP[which(EstimatedP != 0)] <- 1
  total_correct <- sum(TargetP == EstimatedP) # this is the total number of variables correctedly selected and zeros correctly retained
  prop_correct <- total_correct/total_vnumber
  return(prop_correct)
}

IS_rESEM_SCAD <- function(DATA, R, lambda, ASCAD, RHO, MaxIter, eps, nstarts){
  J <- dim(DATA)[2]
  
  VarSelect0 <- MULTISTART_rESEM_SCAD(DATA, R, lambda, ASCAD, RHO, MaxIter, eps, nstarts, penalty=FALSE)
  P_hat0 <- VarSelect0$loadings
  T_hat0 <- VarSelect0$scores
  
  V_oo <- sum(DATA^2)
  V_s <- sum((T_hat0%*%t(P_hat0))^2) 
  
  VarSelect <- MULTISTART_rESEM_SCAD(DATA, R, lambda, ASCAD, RHO, MaxIter, eps, nstarts)
  P_hat <- VarSelect$loadings
  T_hat <- VarSelect$scores
  
  card <- sum(P_hat != 0)
  
  V_a <- sum((T_hat %*% t(P_hat))^2)
  IS <- list()
  IS$value <- (V_a * V_s / V_oo^2) * (sum(round(P_hat,3) == 0) /(J*R))
  IS$vaf <- V_a/V_oo
  IS$propzero <- sum(round(P_hat,3) == 0)/(J*R)
  
  IS$smallestP <- ifelse(sum(rowSums(P_hat != 0)) < sum(card),
                         0, min(abs(P_hat[P_hat != 0]))
  )
  IS$maxsdP <- max(apply(P_hat, 2, sd))
  return(IS)
}