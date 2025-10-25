library(mvtnorm)
library(extraDistr) #for rcat
library(Matrix)

library(truncnorm) #draw from truncated normal
library(mclust)
library(mcclust.ext)
library(coda)

library(parallel) #mclapply

#----------------------------- Overall equation -------------------------
mvn_hdp <- function(Y, niter, J, burn_in = 1000, thinning = 5, 
                    Z_fix = NULL, empirical_z = NULL,
                    mu0 = NULL, k0 = 0.01, omega0 = NULL, Phi0 = NULL,
                    save.only.z = FALSE, MH.variance = 0.01, 
                    partial.save.name = NULL, save_frequency = 100, 
                    auto.save = FALSE, save_ind = NULL, verbose = TRUE){
  
  D <- length(Y); G <- ncol(Y[[1]])
  
  C <- sapply(Y, nrow)
  C_cum <- c(0,cumsum(C))
  
  #------------------------ Step 0: Prepare for outputs -----------------
  Z_output <- NULL
  
  
  P_output <- NULL
  alpha_output <- c()
  alpha_0_output <- c()
  mu_star_1_J_output <- NULL
  Sigma_star_1_J_output <- NULL
  
  P_J_D_output <- NULL
  
  sd_P_output <- c()
  #----------------------- Step 1: Initial values in MCMC and prior setup -----------
  
  # ------ Initial Z --------
  if(is.null(Z_fix)){
    if(empirical_z==TRUE){
      Y_all <- do.call(rbind, Y)
      km_cluster <- kmeans(Y_all,centers=J)$cluster
      Z_new <- lapply(1:D,
                      function(m) km_cluster[(C_cum[m]+1):C_cum[m+1]])
      
    }else{
      Z_new <- lapply(1:D, function(d) sample(1:J, size = C[[d]],replace = TRUE,prob = rep(1/J,J)))
    }
  }else{
    Z_new <- Z_fix
  }
  
  # ----------- Initial P and p_j_d ---------------
  # Component probabilities p_j, add 1 to avoid zero p if cluster is empty
  P_new <- sapply(1:J, function(j) mean(unlist(Z_new)==j))
  if(any(P_new==0)){
    # add 1 to avoid zero p if cluster is empty
    P_new <- sapply(1:J, function(j) sum(unlist(Z_new)==j)+1)/(sum(C)+J)
  }
  
  # Dataset-specific vector p_j_d
  P_J_D_new <- do.call(cbind, lapply(1:D, function(d) {
    val <- sapply(1:J,function(j) sum(Z_new[[d]]==j))+1
    return(val)
  }))
  
  # ----------- Initial mu, Sigma and prior setup ----------
  
  if(is.null(mu0)) {
    mu0 <- colMeans(do.call(rbind, Y))
  }
  
  if(is.null(omega0)){
    omega0 <- G+2
  }
  
  if(is.null(Phi0)){
    Phi0 <- cov(do.call(rbind, Y))
  }
  
  # ----------- Initial alpha, alpha0 ---------
  
  alpha_new <- 1
  alpha_0_new <- 1
  
  
  #----------------------- Step 2:  Acceptance probability ----------------------------
  
  acceptance_count_avg <- data.frame(P_accept = rep(0,niter), alpha_accept = rep(0,niter),
                                     alpha_0_accept = rep(0,niter))
  
  # Total count of acceptance
  P_count <- 0; alpha_count <- 0; alpha_0_count <- 0
  
  #----------------------- Step 3: Prepare for the adaptive covariance update -----------------------------
  
  # 1) For Component probabilities
  sd_P_new <- 0.001
  #X_n (n: indicate iterations)
  mean_X_component_new <- log(matrix(P_new[1:(J-1)]/P_new[J], nrow = 1)) #1x(J-1)
  tilde_s_component_new <- Matrix::crossprod(mean_X_component_new)
  covariance_component_new <- matrix(0, nrow = J-1, ncol = J-1)
  
  # 2) For alpha
  mean_X_alpha_new <- log(alpha_new)
  M_2_alpha_new <- 0
  variance_alpha_new <- 0
  
  # 3) For alpha_0
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
    
    
    # 1) ------- Update dataset-specific probabilities p_j_d --------
    P_J_D_new <- P_J_D_update(Z = Z_new, alpha = alpha_new, P = P_new)
    
    # 2) ------- Update the allocation variable --------
    if(is.null(Z_fix)){
      Z_new <- allocation_variables_update(Y = Y, mu_star_1_J = mu_star_1_J_new, Sigma_star_1_J = Sigma_star_1_J_new,
                                           P_J_D = P_J_D_new)
    }
    
    # 3) ------ Update the component probabilities P -------
    component_output <- component_probabilities_update(P = P_new, P_J_D = P_J_D_new, alpha_0 = alpha_0_new,
                                                       alpha = alpha_new, covariance = covariance_component_new,
                                                       mean_x = mean_X_component_new, 
                                                       tilde_s = tilde_s_component_new,
                                                       iter_num = iter, sd_P = sd_P_new, 
                                                       MH.variance = MH.variance)
    
    P_new <- component_output$P_new
    tilde_s_component_new <- component_output$tilde_s_new
    mean_X_component_new <- component_output$mean_x_new
    covariance_component_new <- component_output$covariance_new
    sd_P_new <- component_output$sd_P_new
    P_count <- P_count + component_output$accept # cumulative count
    acceptance_count_avg$P_accept[iter-1] <- P_count/(iter-1) # update the acceptance data.frame
    
    sd_P_output[iter-1] <- sd_P_new
    
    # 4) ------- Update alpha -------------
    alpha_output_sim <- alpha_update(P_J_D = P_J_D_new, P = P_new, alpha = alpha_new,
                                     X_mean = mean_X_alpha_new, M_2 = M_2_alpha_new, 
                                     variance = variance_alpha_new, iter_num = iter,
                                     MH.variance = MH.variance)
    
    alpha_new <- alpha_output_sim$alpha_new
    mean_X_alpha_new <- alpha_output_sim$X_mean_new
    M_2_alpha_new <- alpha_output_sim$M_2_new
    variance_alpha_new <- alpha_output_sim$variance_new
    alpha_count <- alpha_count + alpha_output_sim$accept # cumulative count
    acceptance_count_avg$alpha_accept[iter-1] <- alpha_count/(iter-1) # update the acceptance probability

    # 5) ------- Update alpha_0 --------
    alpha_0_output_sim <- alpha_0_update(P = P_new, alpha_0 = alpha_0_new, 
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
      
      P_output[[output_index]] <- P_new
      alpha_output[output_index] <- alpha_new
      alpha_0_output[output_index] <- alpha_0_new
      
      mu_star_1_J_output[[output_index]] <- mu_star_1_J_new
      Sigma_star_1_J_output[[output_index]] <- Sigma_star_1_J_new
      
      P_J_D_output[[output_index]] <- P_J_D_new
      
    }
    
    if((iter-1) %% save_frequency == 0 && auto.save == TRUE){
      
      end_time <- Sys.time()
      diff_time <- difftime(end_time,start_time_mcmc)
      
      if(save.only.z){
        my_list <- list('Z_output' = Z_output, 'acceptance_count_avg' = acceptance_count_avg, 
                        'Z_new' = Z_new, 'P_J_D_new' = P_J_D_new, 
                        'P_new' = P_new, 'tilde_s_component_new' = tilde_s_component_new, 'mean_X_component_new' = mean_X_component_new, 
                        'covariance_component_new' = covariance_component_new, 'P_count' = P_count,
                        'alpha_new' = alpha_new, 'mean_X_alpha_new' = mean_X_alpha_new, 'M_2_alpha_new' = M_2_alpha_new,
                        'variance_alpha_new' = variance_alpha_new, 'alpha_count' = alpha_count, 
                        'alpha_0_new' = alpha_0_new, 'mean_X_alpha_0_new' = mean_X_alpha_0_new, 'M_2_alpha_0_new' = M_2_alpha_0_new, 
                        'variance_alpha_0_new' = variance_alpha_0_new, 'alpha_0_count' = alpha_0_count, 
                        'output_index' = output_index,
                        'sd_P_new' = sd_P_new,
                        'mu0' = mu0, 'k0' = k0, 'omega0' = omega0, 'Phi0' = Phi0,
                        'time' = diff_time)
      }else{
        my_list <- list('Z_output' = Z_output, 'P_output' = P_output, 
                        'alpha_output' = alpha_output, 'alpha_0_output' = alpha_0_output, 
                        'mu_star_1_J_output' = mu_star_1_J_output, 'Sigma_star_1_J_output' = Sigma_star_1_J_output,
                        'P_J_D_output' = P_J_D_output, 
                        'acceptance_count_avg' = acceptance_count_avg, 
                        'Z_new' = Z_new, 'mu_star_1_J_new' = mu_star_1_J_new, 'Sigma_star_1_J_new' = Sigma_star_1_J_new,
                        'P_J_D_new' = P_J_D_new, 
                        'P_new' = P_new, 'tilde_s_component_new' = tilde_s_component_new, 'mean_X_component_new' = mean_X_component_new, 
                        'covariance_component_new' = covariance_component_new, 'P_count' = P_count, 
                        'alpha_new' = alpha_new, 'mean_X_alpha_new' = mean_X_alpha_new, 'M_2_alpha_new' = M_2_alpha_new,
                        'variance_alpha_new' = variance_alpha_new, 'alpha_count' = alpha_count, 
                        'alpha_0_new' = alpha_0_new, 'mean_X_alpha_0_new' = mean_X_alpha_0_new, 'M_2_alpha_0_new' = M_2_alpha_0_new, 
                        'variance_alpha_0_new' = variance_alpha_0_new, 'alpha_0_count' = alpha_0_count, 
                        'output_index' = output_index,
                        'sd_P_new' = sd_P_new,
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
                    'Z_new' = Z_new, 'P_J_D_new' = P_J_D_new, 
                    'P_new' = P_new, 'tilde_s_component_new' = tilde_s_component_new, 'mean_X_component_new' = mean_X_component_new, 
                    'covariance_component_new' = covariance_component_new, 'P_count' = P_count,
                    'alpha_new' = alpha_new, 'mean_X_alpha_new' = mean_X_alpha_new, 'M_2_alpha_new' = M_2_alpha_new,
                    'variance_alpha_new' = variance_alpha_new, 'alpha_count' = alpha_count, 
                    'alpha_0_new' = alpha_0_new, 'mean_X_alpha_0_new' = mean_X_alpha_0_new, 'M_2_alpha_0_new' = M_2_alpha_0_new, 
                    'variance_alpha_0_new' = variance_alpha_0_new, 'alpha_0_count' = alpha_0_count, 
                    'output_index' = output_index,
                    'sd_P_new' = sd_P_new,
                    'mu0' = mu0, 'k0' = k0, 'omega0' = omega0, 'Phi0' = Phi0,
                    'time' = diff_time)
  }else{
    my_list <- list('Z_output' = Z_output, 'P_output' = P_output, 
                    'alpha_output' = alpha_output, 'alpha_0_output' = alpha_0_output, 
                    'mu_star_1_J_output' = mu_star_1_J_output, 'Sigma_star_1_J_output' = Sigma_star_1_J_output,
                    'P_J_D_output' = P_J_D_output, 
                    'acceptance_count_avg' = acceptance_count_avg, 
                    'Z_new' = Z_new, 'mu_star_1_J_new' = mu_star_1_J_new, 'Sigma_star_1_J_new' = Sigma_star_1_J_new,
                    'P_J_D_new' = P_J_D_new, 
                    'P_new' = P_new, 'tilde_s_component_new' = tilde_s_component_new, 'mean_X_component_new' = mean_X_component_new, 
                    'covariance_component_new' = covariance_component_new, 'P_count' = P_count, 
                    'alpha_new' = alpha_new, 'mean_X_alpha_new' = mean_X_alpha_new, 'M_2_alpha_new' = M_2_alpha_new,
                    'variance_alpha_new' = variance_alpha_new, 'alpha_count' = alpha_count, 
                    'alpha_0_new' = alpha_0_new, 'mean_X_alpha_0_new' = mean_X_alpha_0_new, 'M_2_alpha_0_new' = M_2_alpha_0_new, 
                    'variance_alpha_0_new' = variance_alpha_0_new, 'alpha_0_count' = alpha_0_count, 
                    'output_index' = output_index,
                    'sd_P_new' = sd_P_new,
                    'mu0' = mu0, 'k0' = k0, 'omega0' = omega0, 'Phi0' = Phi0,
                    'time' = diff_time)
    
  }
  
  return(my_list)
  
}

