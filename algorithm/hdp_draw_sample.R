##---------------- MCMC Simulation --------------

## --------------- Simulation of mean and Sigma within each component -------
unique_params_update <- function(Y, J, Z, mu0, k0, omega0, Phi0){
  
  D <- length(Y)
  
  loop.result <- lapply(1:J,function(j) {
    
    N_j <- sum(unlist(Z)==j)
    
    if(N_j!=0){
      
      Y_j <- do.call(rbind,lapply(1:D,function(d) Y[[d]][Z[[d]]==j,]))
      if(N_j==1){
        Y_j <- matrix(Y_j, nrow=1)
      }
      
      k_n <- k0+N_j
      mu_n <- (k0*mu0+colSums(Y_j))/k_n
      Y_mean <- colMeans(Y_j)
      
      if(N_j==1){
        temp <- 0
      }else{
        temp <- lapply(1:N_j, function(i) {
          Matrix::tcrossprod(matrix(Y_j[i,]-Y_mean,ncol=1))
        })
        temp <- Reduce('+',temp)
      }
      
      Phi_n <- Phi0+temp+k0*N_j/k_n*Matrix::tcrossprod(matrix(Y_mean-mu0,ncol=1))
      omega_n <- N_j+omega0
      
      Sigma_j_new <- MCMCpack::riwish(omega_n, Phi_n)
      mu_j_new <- mniw::rmNorm(n=1, mu_n, Sigma_j_new/k_n)
      
    }else{
      # for empty clusters, draw from the prior
      Sigma_j_new <- MCMCpack::riwish(omega0, Phi0)
      mu_j_new <- mniw::rmNorm(1, mu0, Sigma_j_new/k0)
      
    }
    
    return(list(mu_j_new=mu_j_new,Sigma_j_new=Sigma_j_new))
    
  })
  
  
  mu_new <- do.call(rbind, lapply(loop.result, function(l) l$mu_j_new))
  Sigma_new <- lapply(loop.result, function(l) l$Sigma_j_new)
  
  return(list(mu_new=mu_new, Sigma_new=Sigma_new))
  
  
}


##--------------- Simulation of dataset-specific probabilities p_j_d ----------------

P_J_D_update <- function(Z, alpha, P){
  
  J <- length(P)
  D <- length(Z)
  
  
  loop.result <- lapply(1:D, function(d) {
    
    parameter <- as.vector(table(factor(Z[[d]], levels = c(1:J)))) + alpha*P
    P_J_d <- rdirichlet(n = 1, alpha = parameter)
    
    ## If any of the dataset-specific component probability equal to zero, then give it a small value
    P_J_d <- ifelse(P_J_d == 0, 0.001, P_J_d)
    P_J_d <- P_J_d/sum(P_J_d)
    
    return(P_J_d)
  })
  
  P_J_D <- matrix(unlist(loop.result), nrow=J, ncol=D)
  
  return(P_J_D)
}



##--------------- Simulation of allocations ----------------

allocation_variables_update <- function(Y, mu_star_1_J, Sigma_star_1_J, P_J_D){
  
  D <- length(Y)
  C <- unlist(lapply(Y, nrow))
  J <- nrow(mu_star_1_J)
  
  # Set up the list to save updated values
  Z <- NULL
  
  for(d in 1:D){
    
    # n * J
    P_mat <- matrix(rep(P_J_D[,d],C[d]),nrow=C[d],byrow = TRUE)
    
    # n * J
    LP1 <- log(P_mat)
    LP2 <- lapply(1:J,function(j) {
      temp <- mniw::dmNorm(as.matrix(Y[[d]]), mu=mu_star_1_J[j,], Sigma = Sigma_star_1_J[[j]], log = TRUE)
      return(temp)
    })
    LP2 <- do.call(cbind,LP2)
    LP <- LP1+LP2
    
    nc <- -apply(LP,1,max) # length=n
    # n * J
    LP_plus_nc <- LP+matrix(rep(nc,J),ncol=J,byrow = FALSE)
    
    P <- t(apply(LP_plus_nc, 1, function(x) {
      exp(x)/sum(exp(x))
    }))
    
    Z[[d]] <- apply(P, 1, function(x) rcat(1, prob=x))
  }
  return(Z)
}



##----------- Simulation of Component probabilities P----------------

component_log_prob <- function(P, P_J_D, alpha_0, alpha){
  
  J <- nrow(P_J_D)
  D <- ncol(P_J_D)
  
  lprod <- sum((alpha_0/J-1)*log(P))
  for(d in 1:D){
    lprod <- lprod + sum(alpha*P*log(P_J_D[,d])-lgamma(alpha*P))
  }
  
  # Returning outputs
  return(lprod)
}

# 2) Simulation
component_probabilities_update <- function(P, P_J_D, alpha_0, alpha, covariance,
                                           mean_x, tilde_s, iter_num, sd_P, MH.variance){
  
  J <- nrow(P_J_D)
  D <- ncol(P_J_D)
  
  
  # Define the inputs
  P_old <- P
  covariance_old <- covariance
  mean_x_old <- mean_x
  tilde_s_old <- tilde_s
  sd_P_old <- sd_P
  
  X_old <- log(P_old[1:(J-1)]/P_old[J]) # Length = J-1
  
  # Specify the iteration number
  n <- iter_num
  
  # Adaptive step
  if(n <= 100){
    X_new <- mvtnorm::rmvnorm(n = 1, mean = X_old, sigma = 0.01*diag(x=1,nrow = J-1, ncol = J-1))
  }else{
    X_new <- mvtnorm::rmvnorm(n = 1, mean = X_old, sigma = sd_P_old*(covariance_old + MH.variance*diag(1, nrow = J-1, ncol = J-1)))
  }
  
  # Compute P_new (Length = J) from X_new
  P_new <- c(exp(X_new)/(1+sum(exp(X_new))),1/(1+sum(exp(X_new))))
  
  if(any(P_new==0)){
    print(P_new)
  }
  
  # Compute acceptance probability
  log_acceptance <- component_log_prob(P_new, P_J_D, alpha_0, alpha) -
    component_log_prob(P_old, P_J_D, alpha_0, alpha) +
    sum(log(P_new)-log(P_old))
  acceptance_P <- exp(log_acceptance)
  acceptance_P <- min(1,acceptance_P)
  if(is.na(acceptance_P)) {
    print(acceptance_P)
  }
  
  sd_P_new <- exp(log(sd_P_old)+n^(-0.7)*(acceptance_P-0.234))
  if(sd_P_new>exp(50)) {sd_P_new <- exp(50)}
  if(sd_P_new<exp(-50)) {sd_P_new <- exp(-50)}
  # print(sd_P_new)
  
  outcome <- rbinom(n = 1, size = 1, prob=acceptance_P)
  if(is.na(outcome) == TRUE | outcome == 0){
    X_new <- X_old
    P_new <- P_old
    accept <- 0
  }else{
    accept <- 1
  }
  
  # Update covariance, mean_x and tilde_s
  tilde_s_new <- tilde_s_old + matrix(X_new, ncol = 1)%*%matrix(X_new, nrow = 1)
  mean_x_new <- mean_x_old*(1-1/n) + 1/n*matrix(X_new, nrow = 1)
  covariance_new <- 1/(n-1)*tilde_s_new - n/(n-1)*t(mean_x_new)%*%mean_x_new
  
  return(list(P_new=P_new, tilde_s_new=tilde_s_new, mean_x_new=mean_x_new, 
              covariance_new=covariance_new, accept=accept, accept_prob=acceptance_P, sd_P_new=sd_P_new))
}


## ----------------------- Simulation of Concentration parameter alpha -----------------------
alpha_log_prob <- function(P_J_D, P, alpha){
  
  D <- ncol(P_J_D)
  
  
  # Construct the log-probability
  lprod <- -alpha + D*lgamma(alpha) + sum(vapply(1:D, function(d) {
    sum(alpha*P*log(P_J_D[,d])-lgamma(alpha*P))
  }, FUN.VALUE = numeric(1)))
  
  
  # Return the log-probability
  return(lprod)
}

alpha_update <- function(P_J_D, P, alpha, X_mean, M_2, variance, iter_num, MH.variance){
  
  D <- ncol(P_J_D)
  
  # Defining the inputs
  alpha_old <- alpha; X_old <- log(alpha_old)
  X_mean_old <- X_mean
  M_2_old <- M_2
  variance_old <- variance
  
  # Defining the dimensions
  n <- iter_num
  
  # Apply AMH based on the iterative number of the current iteration
  # to simulated new value of X
  if(n <= 100){
    X_new <- rnorm(n = 1, mean = X_old, sd = 0.1)
  }else{
    X_new <- rnorm(n = 1, mean = X_old, sd = sqrt(2.4^2*variance_old + 2.4^2*MH.variance))
  }
  
  # Transform the new value of X back to new value of alpha, namely alpha_new
  alpha_new <- exp(X_new)
  
  # Compute log acceptance probability
  log_acceptance <- alpha_log_prob(P_J_D, P, alpha = alpha_new) -
    alpha_log_prob(P_J_D, P, alpha = alpha_old) +
    log(alpha_new) - log(alpha_old)
  acceptance_alpha <- exp(log_acceptance)
  acceptance_alpha <- min(1,acceptance_alpha)
  
  # Update X_alpha
  outcome <- rbinom(n = 1, size = 1, prob = acceptance_alpha)
  if(is.na(outcome) == TRUE | outcome == 0){
    X_new <- X_old
    alpha_new <- alpha_old
    accept <- 0
  }else{
    accept <- 1
  }
  
  X_mean_new <- (1-1/n)*X_mean_old + 1/n*X_new
  M_2_new <- M_2_old + (X_new-X_mean_old)*(X_new-X_mean_new)
  variance_new <- 1/(n-1)*M_2_new
  
  # Returning a list of outputs:
  return(list(alpha_new=alpha_new, X_mean_new=X_mean_new, M_2_new=M_2_new, 
              variance_new=variance_new, accept=accept, accept_prob=acceptance_alpha))
}


##---------------------- Simulation of concentration parameter alpha_0 -----------------
alpha_0_log_prob <- function(P,alpha_0){
  
  J <- length(P)
  
  lprob <- -alpha_0 + lgamma(alpha_0) - J*lgamma(alpha_0/J) + sum(alpha_0/J*log(P))
  
  # Return the log-probability
  return(lprob)
}

alpha_0_update <- function(P, alpha_0, X_mean, M_2, variance, iter_num, MH.variance){
  
  J <- length(P)
  
  alpha_0_old <- alpha_0; X_old <- log(alpha_0_old)
  variance_old <- variance
  M_2_old <- M_2
  X_mean_old <- X_mean
  
  # AMH
  n <- iter_num
  if(n <= 100){
    X_new <- rnorm(n = 1, mean = X_old, sd = 0.1)
  }else{
    X_new <- rnorm(n = 1, mean = X_old, sd = sqrt(2.4^2*variance_old + 2.4^2*MH.variance))
  }
  
  # Obtain the new simulated value for alpha_0
  alpha_0_new <- exp(X_new)
  
  # Compute acceptance probability
  log_acceptance <- alpha_0_log_prob(P,alpha_0 = alpha_0_new) -
    alpha_0_log_prob(P, alpha_0 = alpha_0_old) +
    log(alpha_0_new) - log(alpha_0_old)
  acceptance_alpha0 <- exp(log_acceptance)
  acceptance_alpha0 <- min(1,acceptance_alpha0)
  
  # Update X_alpha_0 and output alpha_0_new
  outcome <- rbinom(n = 1, size = 1, prob = acceptance_alpha0)
  if(is.na(outcome) == TRUE | outcome == 0){
    X_new <- X_old
    alpha_0_new <- alpha_0_old
    accept <- 0
  }else{
    accept <- 1
  }
  
  
  X_mean_new <- (1-1/n)*X_mean_old + 1/n*X_new
  M_2_new <- M_2_old + (X_new-X_mean_old)*(X_new-X_mean_new)
  variance_new <- 1/(n-1)*M_2_new
  
  return(list(alpha_0_new=alpha_0_new, X_mean_new=X_mean_new, M_2_new=M_2_new, 
              variance_new=variance_new, accept=accept, accept_prob=acceptance_alpha0))
}

