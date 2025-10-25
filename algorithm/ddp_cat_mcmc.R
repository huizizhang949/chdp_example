library(mvtnorm)
library(extraDistr) #for rcat
library(Matrix)

library(truncnorm) #draw from truncated normal
library(mclust)
library(mcclust.ext)
library(coda)

library(parallel) #mclapply

#----------------------------- Overall equation -------------------------
mvn_ddp1 <- function(Y, x, xc, niter, J, burn_in = 1000, thinning = 5, trunc_dirichlet_burn=50,
                     Z_fix = NULL, empirical_z = NULL,
                     mu0 = NULL, k0 = 0.01, omega0 = NULL, Phi0 = NULL,
                     mu_r = 0.5, sigma_r = 0.5, eta_1 = 5, eta_2 = 1, mu_h = -5,
                     sigma_h = 0.5, kappa_1 = 5, kappa_2 = 1, gamma = NULL, 
                     save.only.z = FALSE, MH.variance = 0.01, 
                     partial.save.name = NULL, save_frequency = 100, 
                     auto.save = FALSE, save_ind = NULL, verbose = TRUE){
  
  # Define Gaussian kernel (Radial basis function)
  rbf <- function(x, x_star, sigma_star_2) {exp(-(x-x_star)^2/2/sigma_star_2)}
  
  G <- ncol(Y)
  
  C <- nrow(Y)
  
  # number of unique values in categorical x
  R <- length(unique(xc))
  #------------------------ Step 0: Prepare for outputs -----------------
  Z_output <- NULL
  
  P_C_J_output <- NULL
  
  alpha_0_output <- c()
  mu_star_1_J_output <- NULL
  Sigma_star_1_J_output <- NULL
  
  Q_J_output <- NULL
  Xi_C_output <- NULL
  x_star_J_output <- NULL
  sigma_star_2_J_output <- NULL
  
  U_C_J_output <- NULL
  
  r_J_output <- NULL
  s_2_output <- c()
  h_J_output <- NULL
  m_2_output <- c()
  
  rho_output <- list()
  #----------------------- Step 1: Initial values in MCMC and prior setup -----------
  
  # ------ Initial Z --------
  if(is.null(Z_fix)){
    if(empirical_z==TRUE){
      
      Z_new <- kmeans(Y,centers=J)$cluster
      
    }else{
      Z_new <- sample(1:J, size = C,replace = TRUE,prob = rep(1/J,J))
    }
  }else{
    Z_new <- Z_fix
  }
  
  
  # ---------- Initial kernel and hyper-paramerters ---------
  # Compute the mean of x within each cluster to initialize x_star_J
  # There may be empty clusters
  # Fill the non-empty values into the initialization matrix, impute NA values with dataset-specific mean
  x_star_J_new <- rep(NA, J)
  # Mean of x within each cluster of each dataset
  val <- tapply(x, Z_new, mean)
  x_star_J_new[as.numeric(names(val))] <- val
  x_star_J_new[is.na(x_star_J_new)] <- mean(x)
  
  rm(val)
  
  # Initials for r_J
  r_J_new <- mean(x_star_J_new)
  
  # Initials for s^2: prior mean = 1/4
  # s_2_new <- eta_2/(eta_1-1)
  s_2_new <- var(x_star_J_new)
  
  # Initials for sigma_star_2_J
  sigma_star_2_J_new <- rep(NA, J)
  # Mean of x within each cluster of each dataset
  val <- tapply(x, Z_new, var)
  sigma_star_2_J_new[as.numeric(names(val))] <- val
  sigma_star_2_J_new[is.na(sigma_star_2_J_new)] <- var(x)
  
  rm(val)
  
  # Initials for h_J from log(sigma_star_2_J)
  h_J_new <- mean(log(sigma_star_2_J_new))
  
  # Initials for m^2: prior mean = 1/4
  # m_2_new <- kappa_2/(kappa_1-1)
  m_2_new <- var(log(sigma_star_2_J_new))
  
  
  # categorical x Dirichlet prior parameters
  if(is.null(gamma)){
    gamma <- rep(1,R)
  }
  # ----------- Initial Q ---------------
  # add 1 to avoid zero p if cluster is empty
  Q_J_new <- sapply(1:J, function(j) mean(Z_new==j))
  if(any(Q_J_new==0)){
    Q_J_new <- sapply(1:J, function(j) sum(Z_new==j)+1)/(C+J)
  }
  
  # ----------- Initial mu, Sigma and prior setup ----------
  
  if(is.null(mu0)) {
    mu0 <- colMeans(Y)
  }
  
  if(is.null(omega0)){
    omega0 <- G+2
  }
  
  if(is.null(Phi0)){
    Phi0 <- cov(Y)
  }

  # ----------- Initial alpha0 ---------
  
  alpha_0_new <- 1
  
  # ---------- Initial Xi and U --------
  
  rho_J_R_new <- rdirichlet(J, alpha=gamma)
  
  Xi_C_new <- Xi_C_update(Q_J = Q_J_new, C = C, x = x, xc = xc, x_star_J = x_star_J_new, 
                          sigma_star_2_J = sigma_star_2_J_new, rho=rho_J_R_new)
  
  U_C_J_new <- U_C_J_update(Xi_C = Xi_C_new, Q_J = Q_J_new, C = C, x = x, xc = xc,
                            x_star_J = x_star_J_new, sigma_star_2_J = sigma_star_2_J_new, rho = rho_J_R_new, rbf = rbf)
  
  #----------------------- Step 2:  Acceptance probability ----------------------------
  
  acceptance_count_avg <- data.frame(alpha_0_accept = rep(0,niter), sigma_accept = rep(0, niter))
  
  # Total count of acceptance
  alpha_0_count <- 0;  sigma_star_2_count <- 0
  
  #----------------------- Step 3: Prepare for the adaptive covariance update -----------------------------
  
  
  # 0) For sigma_star_2_J
  # X = -log(1/sigma_star - 1/upper)
  mean_X_sigma_star_2_new <- log(sigma_star_2_J_new)
  M_2_sigma_star_2_new <- rep(0, J)
  variance_sigma_star_2_new <- rep(0, J)
  
  for (j in 1:J) {
    ind <- c(1:C)[-log(U_C_J_new[,j])<Xi_C_new*Q_J_new[j]*rho_J_R_new[j,xc]]
    if(length(ind)!=0) {
      uppers <- -(x-x_star_J_new[j])^2/2/(log(-log(U_C_J_new[,j]))-log(Xi_C_new)-log(Q_J_new[j])-log(rho_J_R_new[j,xc]))
      upper <- min(uppers[ind])
      mean_X_sigma_star_2_new[j] <- -log(1/sigma_star_2_J_new[j] - 1/upper)
    }
  }
  
  
  # 1) For alpha_0
  mean_X_alpha_0_new <- log(alpha_0_new)
  M_2_alpha_0_new <- 0
  variance_alpha_0_new <- 0
  
  output_index <- 0
  
  start_time_mcmc <- Sys.time() 
  if(verbose) {
    print(paste('Start MCMC:', start_time_mcmc))
    cat('\n')
    
    pb <- txtProgressBar(min = 1, max = niter+1, style = 3)
  }
  #----------------------- Step 4: Updates -----------------------------
  # Iteration starts with iter_num = 2
  for(iter in 2:(niter+1)){
    
    # print(paste('Iteration', iter))
    if(verbose) {setTxtProgressBar(pb, iter)}
    # 0) Starting value of the output index = 1
    # If the current iteration is greater than the burn in and divisible by the thinning index
    
    # In the case of thinning, if we require to save some samples that are not devisible by thinning (for consesnsus clustering) or before burn-in
    if(is.null(save_ind)) {
      criterion1 <- FALSE
    }else{
      criterion1 <- any((iter-1)==save_ind)
    }
    
    criterion2 <- (iter-1 > burn_in & (iter-1-burn_in)%%thinning == 0)
    
    if(criterion1 | criterion2){
      output_index <- output_index + 1
      update <- TRUE
    }else{
      update <- FALSE
    }
    
    # 0) ------ Update mean and Sigma within each component -------------
    unique_params_output_sim <- unique_params_update(Y = Y, J = J, Z = Z_new, mu0 = mu0, k0 = k0, omega0 = omega0, Phi0 = Phi0)
    
    mu_star_1_J_new <- unique_params_output_sim$mu_new
    Sigma_star_1_J_new <- unique_params_output_sim$Sigma_new
    
    # 1) ------ Update latent/auxiliary variables Xi_C ---------
    Xi_C_new <- Xi_C_update(Q_J = Q_J_new, C = C, x = x, xc = xc, x_star_J = x_star_J_new, 
                            sigma_star_2_J = sigma_star_2_J_new, rho = rho_J_R_new)
    
    # 2) ------- Update dataset-specific vector q_j --------
    Q_J_new <- Q_J_update(Z = Z_new, alpha_0 = alpha_0_new, Xi_C = Xi_C_new,
                          x = x, xc = xc, x_star_J = x_star_J_new, sigma_star_2_J = sigma_star_2_J_new, rho = rho_J_R_new)
    
    # 3) ------- Update the allocation variable --------
    if(is.null(Z_fix)){
      Z_new <- allocation_variables_update(Y = Y, x = x, xc = xc, mu_star_1_J = mu_star_1_J_new, Sigma_star_1_J = Sigma_star_1_J_new,
                                           Q_J = Q_J_new, x_star_J = x_star_J_new,
                                           sigma_star_2_J = sigma_star_2_J_new, rho = rho_J_R_new)
    }
    
    # 4) ------- Update kernel parameters ----------
    # 4-1) ------ latent variable U_C_J ---------
    U_C_J_new <- U_C_J_update(Xi_C = Xi_C_new, Q_J = Q_J_new, C = C, x = x, xc = xc,
                              x_star_J = x_star_J_new, sigma_star_2_J = sigma_star_2_J_new, rho = rho_J_R_new, rbf = rbf)
    
    
    # 4-2) ------- x_star_J -------
    x_star_J_new <- x_star_J_update(r_J = r_J_new, s_2 = s_2_new, Z = Z_new, x = x, xc = xc, C = C,
                                    sigma_star_2_J = sigma_star_2_J_new, U_C_J = U_C_J_new,
                                    Xi_C = Xi_C_new, Q_J = Q_J_new, rho = rho_J_R_new)
    
    # update U again
    U_C_J_new <- U_C_J_update(Xi_C = Xi_C_new, Q_J = Q_J_new, C = C, x = x, xc = xc, 
                              x_star_J = x_star_J_new, sigma_star_2_J = sigma_star_2_J_new, rho = rho_J_R_new, rbf = rbf)
    
    
    
    # 4-3) ------- sigma_star_2_J -------
    sigma_star_output <- sigma_star_2_J_update(sigma_star_2_J = sigma_star_2_J_new,
                                               h_J = h_J_new, m_2 = m_2_new, Z = Z_new, x = x, xc = xc, C = C,
                                               x_star_J = x_star_J_new, U_C_J = U_C_J_new,
                                               Xi_C = Xi_C_new, Q_J = Q_J_new, rho = rho_J_R_new,
                                               X_mean = mean_X_sigma_star_2_new, M_2 = M_2_sigma_star_2_new,
                                               variance = variance_sigma_star_2_new, iter_num = iter,
                                               MH.variance = MH.variance)
    
    sigma_star_2_J_new <- sigma_star_output$sigma_star_2_J_new
    mean_X_sigma_star_2_new <- sigma_star_output$X_mean_new
    M_2_sigma_star_2_new <- sigma_star_output$M_2_new
    variance_sigma_star_2_new <- sigma_star_output$variance_new
    sigma_star_2_count <- sigma_star_2_count + sigma_star_output$accept
    acceptance_count_avg$sigma_accept[iter-1] <- sigma_star_2_count/((iter-1)*J)
    
    # 4-4) --------------- rho_J_R -----------
    rho_J_R_new <- rho_J_R_update(rho = rho_J_R_new, x = x, xc = xc, Z = Z_new, gamma = gamma, U_C_J = U_C_J_new,
                                  Xi_C = Xi_C_new, Q_J = Q_J_new, x_star_J = x_star_J_new, 
                                  sigma_star_2_J = sigma_star_2_J_new, rbf = rbf, trunc_dirichlet_burn = trunc_dirichlet_burn)
    
    # Not necessary, but update P_C_J
    P_C_J_new <- t(sapply(1:C,function(i) {
      xx <- x[i]
      xxc <- xc[i]
      # log-scale
      LP <- log(Q_J_new)-(xx-x_star_J_new)^2/2/sigma_star_2_J_new+log(rho_J_R_new[,xxc])
      nc <- -max(LP)
      P <- exp(LP+nc)/sum(exp(LP+nc))
      return(P)
    }))
    
    
    # 5) ------ Update hyper parameters in priors for x_star, sigma_star_2 ------
    # 5-1) r_J
    r_J_new <- r_J_update(x_star_J = x_star_J_new, mu_r = mu_r, sigma_r = sigma_r, s_2 = s_2_new)
    
    # 5-2) s_2
    s_2_new <- s_2_update(x_star_J = x_star_J_new, eta_1 = eta_1, eta_2 = eta_2, r_J = r_J_new)
    
    # 5-3) h_J
    h_J_new <- h_J_update(sigma_star_2_J = sigma_star_2_J_new, mu_h = mu_h, sigma_h = sigma_h, m_2 = m_2_new)
    
    # 5-4) m_2
    m_2_new <- m_2_update(sigma_star_2_J = sigma_star_2_J_new, kappa_1 = kappa_1, kappa_2 = kappa_2, h_J = h_J_new)
    
    
    # 6) ------- Update alpha_0 --------
    alpha_0_output_sim <- alpha_0_update(Q_J = Q_J_new, alpha_0 = alpha_0_new, 
                                         X_mean = mean_X_alpha_0_new,
                                         M_2 = M_2_alpha_0_new, 
                                         variance = variance_alpha_0_new, iter_num = iter,
                                         MH.variance = MH.variance)

    alpha_0_new <- alpha_0_output_sim$alpha_0_new
    mean_X_alpha_0_new <- alpha_0_output_sim$X_mean_new
    M_2_alpha_0_new <- alpha_0_output_sim$M_2_new
    variance_alpha_0_new <- alpha_0_output_sim$variance_new
    alpha_0_count <- alpha_0_count + alpha_0_output_sim$accept # cumulative count
    acceptance_count_avg$alpha_0_accept[iter-1] <- alpha_0_count/(iter-1) # update the acceptance probability

    
    #-------------------------- Step 5: Update simulated values ------------------------
    if(update == TRUE){
      if(is.null(Z_fix)){
        Z_output[[output_index]] <- Z_new
      }
      
      P_C_J_output[[output_index]] <- P_C_J_new
      
      alpha_0_output[output_index] <- alpha_0_new
      
      mu_star_1_J_output[[output_index]] <- mu_star_1_J_new
      Sigma_star_1_J_output[[output_index]] <- Sigma_star_1_J_new
      
      Q_J_output[[output_index]] <- Q_J_new

      x_star_J_output[[output_index]] <- x_star_J_new
      sigma_star_2_J_output[[output_index]] <- sigma_star_2_J_new
      
      r_J_output[[output_index]] <- as.vector(unname(r_J_new))
      s_2_output[output_index] <- s_2_new
      
      h_J_output[[output_index]] <- as.vector(unname(h_J_new)) 
      m_2_output[output_index] <- m_2_new
      
      rho_output[[output_index]] <- rho_J_R_new
      
    }
    
    if((iter-1) %% save_frequency == 0 && auto.save == TRUE){
      
      end_time <- Sys.time()
      diff_time <- difftime(end_time,start_time_mcmc)
      
      if(save.only.z){
        my_list <- list('Z_output' = Z_output, 'acceptance_count_avg' = acceptance_count_avg, 
                        'Z_new' = Z_new, 'P_C_J_new' = P_C_J_new,
                        'x_star_J_new' = x_star_J_new, 'r_J_new' = r_J_new, 's_2_new' = s_2_new, 'h_J_new' = h_J_new, 'm_2_new' = m_2_new,
                        'Q_J_new' = Q_J_new, 'rho_J_R_new' = rho_J_R_new,
                        'Xi_C_new' = Xi_C_new, 'U_C_J_new' = U_C_J_new,
                        'sigma_star_2_J_new' = sigma_star_2_J_new, 'mean_X_sigma_star_2_new' = mean_X_sigma_star_2_new,
                        'M_2_sigma_star_2_new' = M_2_sigma_star_2_new, 'variance_sigma_star_2_new' = variance_sigma_star_2_new,
                        'sigma_star_2_count' = sigma_star_2_count,
                        'alpha_0_new' = alpha_0_new, 'mean_X_alpha_0_new' = mean_X_alpha_0_new, 'M_2_alpha_0_new' = M_2_alpha_0_new, 
                        'variance_alpha_0_new' = variance_alpha_0_new, 'alpha_0_count' = alpha_0_count, 
                        'output_index' = output_index,
                        'mu0' = mu0, 'k0' = k0, 'omega0' = omega0, 'Phi0' = Phi0,
                        'time' = diff_time)
      }else{
        my_list <- list('Z_output' = Z_output, 'P_C_J_output' = P_C_J_output, 'alpha_0_output' = alpha_0_output, 
                        'mu_star_1_J_output' = mu_star_1_J_output, 'Sigma_star_1_J_output' = Sigma_star_1_J_output,
                        'Q_J_output' = Q_J_output, 'x_star_J_output' = x_star_J_output, 'sigma_star_2_J_output' = sigma_star_2_J_output, 
                        'r_J_output' = r_J_output, 's_2_output' = s_2_output, 'h_J_output' = h_J_output, 'm_2_output' = m_2_output, 
                        'rho_output' = rho_output,
                        'acceptance_count_avg' = acceptance_count_avg, 
                        'Z_new' = Z_new, 'P_C_J_new' = P_C_J_new, 'mu_star_1_J_new' = mu_star_1_J_new, 'Sigma_star_1_J_new' = Sigma_star_1_J_new,
                        'x_star_J_new' = x_star_J_new, 'r_J_new' = r_J_new, 's_2_new' = s_2_new, 'h_J_new' = h_J_new, 'm_2_new' = m_2_new,
                        'Q_J_new' = Q_J_new, 'rho_J_R_new' = rho_J_R_new,
                        'Xi_C_new' = Xi_C_new, 'U_C_J_new' = U_C_J_new,
                        'sigma_star_2_J_new' = sigma_star_2_J_new, 'mean_X_sigma_star_2_new' = mean_X_sigma_star_2_new,
                        'M_2_sigma_star_2_new' = M_2_sigma_star_2_new, 'variance_sigma_star_2_new' = variance_sigma_star_2_new,
                        'sigma_star_2_count' = sigma_star_2_count,
                        'alpha_0_new' = alpha_0_new, 'mean_X_alpha_0_new' = mean_X_alpha_0_new, 'M_2_alpha_0_new' = M_2_alpha_0_new, 
                        'variance_alpha_0_new' = variance_alpha_0_new, 'alpha_0_count' = alpha_0_count, 
                        'output_index' = output_index,
                        'mu0' = mu0, 'k0' = k0, 'omega0' = omega0, 'Phi0' = Phi0,
                        'time' = diff_time)
        
      }
      save(my_list, file=partial.save.name)
    }
    
  }
  
  end_time <- Sys.time()
  diff_time <- difftime(end_time,start_time_mcmc)
  
  if(verbose) {
    close(pb)
    
    cat('\n')
    print(paste('End:',end_time))
    cat('\n')
    print(paste('MCMC running time:', round(diff_time, digits = 3),units(diff_time)))
  }
  
  ## Return the list
  if(save.only.z){
    my_list <- list('Z_output' = Z_output, 'acceptance_count_avg' = acceptance_count_avg, 
                    'Z_new' = Z_new, 'P_C_J_new' = P_C_J_new,
                    'x_star_J_new' = x_star_J_new, 'r_J_new' = r_J_new, 's_2_new' = s_2_new, 'h_J_new' = h_J_new, 'm_2_new' = m_2_new,
                    'Q_J_new' = Q_J_new, 'rho_J_R_new' = rho_J_R_new,
                    'Xi_C_new' = Xi_C_new, 'U_C_J_new' = U_C_J_new,
                    'sigma_star_2_J_new' = sigma_star_2_J_new, 'mean_X_sigma_star_2_new' = mean_X_sigma_star_2_new,
                    'M_2_sigma_star_2_new' = M_2_sigma_star_2_new, 'variance_sigma_star_2_new' = variance_sigma_star_2_new,
                    'sigma_star_2_count' = sigma_star_2_count,
                    'alpha_0_new' = alpha_0_new, 'mean_X_alpha_0_new' = mean_X_alpha_0_new, 'M_2_alpha_0_new' = M_2_alpha_0_new, 
                    'variance_alpha_0_new' = variance_alpha_0_new, 'alpha_0_count' = alpha_0_count, 
                    'output_index' = output_index,
                    'mu0' = mu0, 'k0' = k0, 'omega0' = omega0, 'Phi0' = Phi0,
                    'time' = diff_time)
  }else{
    my_list <- list('Z_output' = Z_output, 'P_C_J_output' = P_C_J_output, 'alpha_0_output' = alpha_0_output, 
                    'mu_star_1_J_output' = mu_star_1_J_output, 'Sigma_star_1_J_output' = Sigma_star_1_J_output,
                    'Q_J_output' = Q_J_output, 'x_star_J_output' = x_star_J_output, 'sigma_star_2_J_output' = sigma_star_2_J_output, 
                    'r_J_output' = r_J_output, 's_2_output' = s_2_output, 'h_J_output' = h_J_output, 'm_2_output' = m_2_output, 
                    'rho_output' = rho_output,
                    'acceptance_count_avg' = acceptance_count_avg, 
                    'Z_new' = Z_new, 'P_C_J_new' = P_C_J_new, 'mu_star_1_J_new' = mu_star_1_J_new, 'Sigma_star_1_J_new' = Sigma_star_1_J_new,
                    'x_star_J_new' = x_star_J_new, 'r_J_new' = r_J_new, 's_2_new' = s_2_new, 'h_J_new' = h_J_new, 'm_2_new' = m_2_new,
                    'Q_J_new' = Q_J_new, 'rho_J_R_new' = rho_J_R_new,
                    'Xi_C_new' = Xi_C_new, 'U_C_J_new' = U_C_J_new,
                    'sigma_star_2_J_new' = sigma_star_2_J_new, 'mean_X_sigma_star_2_new' = mean_X_sigma_star_2_new,
                    'M_2_sigma_star_2_new' = M_2_sigma_star_2_new, 'variance_sigma_star_2_new' = variance_sigma_star_2_new,
                    'sigma_star_2_count' = sigma_star_2_count,
                    'alpha_0_new' = alpha_0_new, 'mean_X_alpha_0_new' = mean_X_alpha_0_new, 'M_2_alpha_0_new' = M_2_alpha_0_new, 
                    'variance_alpha_0_new' = variance_alpha_0_new, 'alpha_0_count' = alpha_0_count, 
                    'output_index' = output_index,
                    'mu0' = mu0, 'k0' = k0, 'omega0' = omega0, 'Phi0' = Phi0,
                    'time' = diff_time)
    
  }
  
  return(my_list)
  
}

