
# ============================================================
# SCAD lambda search via continuation path (robust version)
# ============================================================

# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------

enough_items <- function(P, min_items = 3, eps = 1e-6) {
  colSums(abs(P) > eps) >= min_items
}

degenerate_solution <- function(P, tol = 1e-8) {
  any(colSums(P^2) < tol)
}

# empirical SCAD scaling
lambda_upper_empirical_SCAD <- function(X) {
  G <- abs(crossprod(scale(X)))
  diag(G) <- 0
  quantile(G, 0.95)
}

# ------------------------------------------------------------
# Wrapper around your SCAD estimator
# ------------------------------------------------------------

fit_once_SCAD <- function(DATA,
                          R,
                          lambda,
                          RHO,
                          MaxIter = 250,
                          eps = 1e-5,
                          nstarts = 1)
{
  res <- MULTISTART_rESEM_SCAD(
    DATA     = DATA,
    R        = R,
    lambda   = lambda,
    ASCAD    = 3.7,
    RHO      = RHO,
    MaxIter  = MaxIter,
    eps      = eps,
    nstarts  = nstarts
  )
  
  res$loadings
}

# ------------------------------------------------------------
# SCAD feasibility checker
# ------------------------------------------------------------

is_feasible_SCAD <- function(P,
                             min_items = 3,
                             eps = 1e-6)
{
  if (inherits(P, "try-error"))
    return(FALSE)
  
  if (any(is.na(P)) || any(is.infinite(P)))
    return(FALSE)
  
  if (degenerate_solution(P))
    return(FALSE)
  
  all(enough_items(P, min_items, eps))
}

# ------------------------------------------------------------
# Robust SCAD lambda-max search
# ------------------------------------------------------------

find_lambda_max_SCAD <- function(DATA,
                                 R,
                                 RHO,
                                 MaxIter   = 250,
                                 eps       = 1e-5,
                                 nstarts   = 3,
                                 min_items = 3,
                                 grid_size = 40,
                                 verbose   = TRUE)
{
  # ----------------------------------------------------------
  # empirical upper bound
  # ----------------------------------------------------------
  
  lam_max <- lambda_upper_empirical_SCAD(DATA)
  
  # logarithmic continuation path
  lambda_grid <- exp(seq(log(1e-4),
                         log(lam_max),
                         length.out = grid_size))
  
  feasible_vec <- logical(grid_size)
  
  if (verbose) {
    cat("\n")
    cat("SCAD continuation search\n")
    cat("------------------------\n")
    cat("lambda upper empirical =", signif(lam_max, 4), "\n\n")
  }
  
  # ----------------------------------------------------------
  # continuation search
  # ----------------------------------------------------------
  
  for (k in seq_along(lambda_grid)) {
    
    lam <- lambda_grid[k]
    
    P <- try(
      fit_once_SCAD(
        DATA     = DATA,
        R        = R,
        lambda   = lam,
        RHO      = RHO,
        MaxIter  = MaxIter,
        eps      = eps,
        nstarts  = nstarts
      ),
      silent = TRUE
    )
    
    ok <- is_feasible_SCAD(
      P,
      min_items = min_items,
      eps       = eps
    )
    
    feasible_vec[k] <- ok
    
    if (verbose) {
      
      nz <- if (!inherits(P, "try-error"))
        paste(colSums(abs(P) > eps), collapse = " ")
      else
        "ERROR"
      
      cat(
        sprintf(
          "%2d) lambda = %-10.5f  %s   nonzeros: %s\n",
          k,
          lam,
          ifelse(ok, "FEASIBLE", "INFEASIBLE"),
          nz
        )
      )
    }
  }
  
  # ----------------------------------------------------------
  # largest feasible lambda
  # ----------------------------------------------------------
  
  feasible_idx <- which(feasible_vec)
  
  if (length(feasible_idx) == 0) {
    
    warning("No feasible SCAD solution found.")
    
    return(list(
      lambda_max = NA,
      path       = data.frame(
        lambda   = lambda_grid,
        feasible = feasible_vec
      )
    ))
  }
  
  best_lambda <- max(lambda_grid[feasible_idx])
  
  if (verbose) {
    cat("\n")
    cat("Largest feasible lambda =", signif(best_lambda, 5), "\n")
  }
  
  return(list(
    lambda_max = best_lambda,
    path = data.frame(
      lambda   = lambda_grid,
      feasible = feasible_vec
    )
  ))
}


