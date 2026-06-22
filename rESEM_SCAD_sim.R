#### Title: Regularized ESEM
#### Author: Tra Le
#### Created: may 19 2026
#### Last modified: 

#########################################################################################
######                         rESEM-l1 simulation analysis                         #####
#########################################################################################

# Clear workspace
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
rm(list = ls())

# load simulation setup
load("DATA/Info_simulation.RData")
Info_matrix = Infor_simulation$design_matrix_replication
Ndatasets = Infor_simulation$n_data_sets
MatrixInfo = Infor_simulation$design_matrix_replication
#ind  = missing_numbers

# register cores
library(doParallel)
no_cores <- detectCores()
c1 <- makePSOCKcluster(floor( no_cores*.4))
registerDoParallel(c1)

# run simulation
Simulation = foreach(i = 1:Ndatasets,
                     .combine = rbind)%dopar%{
                       
                       source("C:/Users/20233477/Documents/uni/jaar3/BEP/code/BEP_PRoject/Functions/rESEM-SCAD.R")
                       source("C:/Users/20233477/Documents/uni/jaar3/BEP/code/BEP_PRoject/Functions/IS-SCAD.R")
                       
                       # Load data
                       #i <- ind[id]
                       load(paste0("DATA/data", i, ".RData"))
                       
                       X <- out$X # not standardized
                       y <- out$y # not standardized
                       N <- dim(X)[1]
                       
                       # Population models
                       Ptrue = out$Ptrue
                       Htrue = out$Htrue
                       
                       nfactors <- dim(Ptrue)[2]
                       RHO <- 7
                       
                       lambda_max = find_lambda_max_SCAD(DATA = scale(X), R = nfactors, nstarts = 1, RHO = RHO)

                       for (j in 1:nrow(lambda_max$path)){
                         if(!lambda_max$path$feasible[j]){
                           lam = lambda_max$path$lambda[j-1]
                           break
                         }
                         if (lambda_max$path$feasible[j] & j==nrow(lambda_max$path)){
                           lam = lambda_max$path$lambda[j]
                         }
                       }
                       
                       # run model with the selected lambda
                       model <- MULTISTART_rESEM_SCAD(DATA = scale(X),
                                                      R = nfactors,
                                                      lambda = lam,
                                                      MaxIter = 5000,
                                                      eps = 10^-4,
                                                      nstarts = 50,
                                                      ASCAD = 3.7,
                                                      RHO = RHO)

                       Pmatrix <- model$loadings
                       Hmatrix <- model$scores
                       
                       
                       print(Pmatrix)
                       print(Ptrue)
                       
                       # zero/non-zero structure recovery rate + Tucker's congruence
                       perm <- gtools::permutations(nfactors, nfactors)
                       corrate <- vector(length = nrow(perm))
                       
                       for (p in 1:nrow(perm)) {
                         corrate[p] <-num_correct(Ptrue, round(Pmatrix[,perm[p,]],2))
                         
                       }
                       corrate_final <- max(corrate)
                       
                       
                       index <- which.max(corrate)
                       
                       Pmatrix <- Pmatrix[,perm[index,]]
                       corsign <- sign(diag(cor(Ptrue, Pmatrix)))
                       
                       Pmatrix <- Pmatrix%*%diag(corsign)
                       
                       # structural model
                       Hmatrix <- Hmatrix[,perm[index,]]%*%diag(corsign)
                       tucker <- sum(diag(abs(psych::factor.congruence(Htrue,Hmatrix))))/nfactors
                       beta_est <- lm(y ~ 0 + Hmatrix)
                       
                       # stability
                       n_best <- model$n_best
                       distinct <- model$n_distinct
                       
                       # convergence
                       convergence <- model$converged
                       
                       ## measure summary
                       measure_summary <- matrix(c(corrate_final, tucker, 
                                                   n_best, distinct, convergence), 
                                                 nrow = 1, byrow = T)
                       colnames(measure_summary) <- c("zero_rec", "tucker", 
                                                      "best_solution","distinct_loss", 
                                                      "convergence")
                       
                       # Saving performance measures
                       output_SCAD = list(measure_summary = measure_summary, 
                                        Pmatrix = Pmatrix, 
                                        Hmatrix = Hmatrix,
                                        
                                        beta_est = beta_est,
                                        #beta_known = beta_known,
                                        model = model)
                       
                       save(output_SCAD, file = paste0("SCAD_out1/scad", i, ".RData"))
                       
                     }
# stop Cluster
stopCluster(c1)