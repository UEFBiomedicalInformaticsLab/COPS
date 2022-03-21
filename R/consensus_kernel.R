#<HKH,C>+<HKH,D>+asumij(<Ci,HCjH>)+bsumij(<Ci,HDjH>)
# Solve the ecmc core optimization problem with Rmosek
ecmc_opt <- function(X, 
                     ret_sol = FALSE, 
                     ret_optime = FALSE, 
                     verbosity = 0,
                     parallel = 0) {
  n_mat <- length(X)
  #if (length(X) > 1) if(!all(sapply(X[-1], function(x) identical(X[[1]], x)))) {
  #  stop("Objectives are not of the same size")
  #}
  n_samples <- sapply(X, ncol) #dim(X[[1]])[1]
  # Specify the non-matrix variable part of the problem. 
  prob <- list(sense="max")
  prob$iparam <- list(NUM_THREADS = parallel)
  # Constraints
  prob$A <- Matrix::Matrix(nrow = n_mat, ncol = 0)
  prob$c <- numeric(0)
  # PSD variable constraint (trace == 1)
  prob$bc <- rbind(blc=rep(1, n_mat),
                   buc=rep(1, n_mat))
  prob$bx <- rbind(blx=numeric(0),
                   bux=numeric(0))
  # Specify semidefinite matrix variables
  prob$bardim <- n_samples
  # Block triplet format specifying the lower triangular part
  # of the symmetric coefficient matrix 'barc':
  # Here we make the sparse representation of X (given as input)
  
  n_ind <- (n_samples*n_samples+n_samples)/2 # sparse length
  
  barc_j <- list()
  barc_k <- list()
  barc_l <- list()
  barc_v <- list()
  
  barA_i <- list()
  barA_j <- list()
  barA_k <- list()
  barA_l <- list()
  barA_v <- list()
  for (i in 1:n_mat) {
    # column index
    l <- rep(0, n_ind[i])
    l[cumsum(n_samples[i]:2)+1] <- 1
    l <- cumsum(l) + 1
    # row index
    k <- Reduce("c", lapply(1:n_samples[i], function(x) x:n_samples[i]))
    
    barc_j[[i]] <- rep(i, n_ind[i])
    barc_k[[i]] <- k
    barc_l[[i]] <- l
    barc_v[[i]] <- X[[i]][cbind(k,l)]
    
    # Block triplet format specifying the lower triangular part
    # of the symmetric coefficient matrix 'barA':
    # In this case its the identity matrix
    barA_i[[i]]  <- rep(i, n_samples[i])
    barA_j[[i]]  <- rep(i, n_samples[i])
    barA_k[[i]]  <- 1:n_samples[i]
    barA_l[[i]]  <- 1:n_samples[i]
    barA_v[[i]]  <- rep(1, n_samples[i])
  }
  
  prob$barc$j <- Reduce('c', barc_j)
  prob$barc$k <- Reduce('c', barc_k)
  prob$barc$l <- Reduce('c', barc_l)
  prob$barc$v <- Reduce('c', barc_v)
  
  prob$barA$i <- Reduce('c', barA_i)
  prob$barA$j <- Reduce('c', barA_j)
  prob$barA$k <- Reduce('c', barA_k)
  prob$barA$l <- Reduce('c', barA_l)
  prob$barA$v <- Reduce('c', barA_v)
  
  # Solve the problem
  r <- Rmosek::mosek(prob, 
                     list(soldetail = as.integer(ret_sol), 
                          getinfo = ret_optime,
                          verbose = verbosity))
  
  # Convergence
  if (r$response$code != 0) {
    warning(paste("Received non-zero response from MOSEK:", r$response$msg))
  }
  if (r$sol$itr$solsta != "OPTIMAL") {
    stop(paste0("Non-optimal solution"))
  }
  
  # Return solution
  out <- list()
  out$solution <- list()
  for (i in 1:n_mat) {
    out$solution[[i]] <- new("dspMatrix", x = r$sol$itr$barx[[i]], 
                             uplo="L", Dim=c(n_samples[i],n_samples[i]))
  }
  if(ret_sol) out$val <- r$sol$itr$pobjval
  if(ret_optime) out$time <- r$dinfo$OPTIMIZER_TIME
  return(out)
}


#' Enhanced consensus multi-view clustering
#' 
#' Consensus kernel approach described by Cai & Li (2017).
#'
#' @param x list of kernel matrices
#' @param a consensus reward
#' @param b disagreement penalty
#' @param eps stopping threshold for mean difference in solution
#' @param solver only MOSEK for now
#' @param maxiter maximum number of iterations before terminating
#' @param parallel number of threads for optimization
#'
#' @return
#' @export
ECMC <- function(x, 
                 a, 
                 b, 
                 eps = 1e-6, 
                 solver = "MOSEK", 
                 maxiter = 10,
                 parallel = 0) {
  N_col <- sapply(x, ncol)
  N_row <- sapply(x, ncol)
  if (any(c(N_col, N_row) != N_col[1])) stop("Input dimensions do not match.")
  
  H <- diag(rep(1, times = N_col[1])) - 1 / N_col[1]
  
  C <- lapply(x, function(xi) xi - Matrix::Diagonal(N_col[1]) * 2)
  D <- lapply(x, function(xi) Matrix::Diagonal(N_col[1]) * 2)
  
  # Normalize initial matrices so that they are in the feasible set
  C <- lapply(x, function(xi) xi / Matrix::norm(xi, "F"))
  D <- lapply(x, function(xi) xi / Matrix::norm(xi, "F"))
  
  delta_m <- Inf
  solution_vals_c <- numeric()
  solution_vals_d <- numeric()
  
  iter <- 0
  while(delta_m > eps) {
    iter <- iter + 1
    if (iter > maxiter) {
      warning(paste0("Not converged in ", maxiter, ". Stopping ..."))
      return(list(consensus_score = NA, reconstruction_objective = NA))
    }
    delta_m <- 0
    
    D_sum <- Reduce("+", D)
    M <- list() 
    for (i in 1:length(x)) {
      M[[i]] <- H %*% (x[[i]] + 2 * a * Reduce("+", C[-i]) - b * D_sum) %*% H
    }
    mosek_sol <- ecmc_opt(M, ret_sol = TRUE, ret_optime = TRUE, parallel = parallel)
    Ct <- lapply(mosek_sol$solution, as.matrix)
    solution_vals_c <- c(solution_vals_c, mosek_sol$val)
    for (i in 1:length(x)) {
      delta_m <- delta_m + mean(abs(C[[i]] - Ct[[i]]))
    }
    C <- Ct
    
    C_sum <- Reduce("+", C)
    N <- list()
    for (i in 1:length(x)) {
      N[[i]] <- H %*% (x[[i]] - b * C_sum) %*% H
    }
    mosek_sol <- ecmc_opt(N, ret_sol = TRUE, ret_optime = TRUE)
    Dt <- lapply(mosek_sol$solution, as.matrix)
    solution_vals_d <- c(solution_vals_d, mosek_sol$val)
    
    for (i in 1:length(x)) {
      delta_m <- delta_m + mean(abs(C[[i]] - Ct[[i]]))
    }
    D <- Dt
    
    print(delta_m)
    flush.console()
  }
  
  c_score_1 <- sapply(1:length(x), function(xi) sum((H %*% x[[xi]] %*% H) * C[[xi]]))
  c_score_2 <- sapply(1:length(x), function(xi) sum((H %*% x[[xi]] %*% H) * (C[[xi]] + D[[xi]])))
  consensus_score <- c_score_1 / c_score_2
  
  reconstruction_objective <- 0
  for (i in 1:length(x)) {
    reconstruction_objective <- reconstruction_objective + sum(x[[i]] * (C[[i]] + D[[i]])) #
  }
  
  return(list(C_sum = C_sum, 
              D_sum = D_sum, 
              C = C, 
              D = D, 
              solution_vals_c = solution_vals_c, 
              solution_vals_d = solution_vals_d, 
              consensus_score = consensus_score, 
              reconstruction_objective = reconstruction_objective))
}