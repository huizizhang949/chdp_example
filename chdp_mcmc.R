library(mvtnorm)
library(extraDistr) #for rcat
library(Matrix)

library(truncnorm) #draw from truncated normal
library(mclust)
library(mcclust.ext)
library(coda)

library(parallel) #mclapply
library(pbmcapply) #pbmclapply for progress bar

#----------------------------- Overall equation -------------------------
mvn_chdp <- function(Y, x, niter, J, burn_in = 1000, thinning = 5, 
                     Z_fix = NULL, empirical_z = NULL,
                     mu0 = NULL, k0 = 0.01, omega0 = NULL, Phi0 = NULL,
                     mu_r = 0.5, sigma_r = 0.5, eta_1 = 5, eta_2 = 1, mu_h = -5,
                     sigma_h = 0.5, kappa_1 = 5, kappa_2 = 1, 
                     save.only.z = FALSE, MH.variance = 0.01, 
                     partial.save.name = NULL, save_frequency = 100, 
                     auto.save = FALSE, save_ind = NULL, verbose = TRUE){
  
  # Define Gaussian kernel (Radial basis function)
  rbf <- function(x, x_star, sigma_star_2) {exp(-(x-x_star)^2/2/sigma_star_2)}
  
  D <- length(Y); G <- ncol(Y[[1]])
  
  C <- sapply(Y, nrow)
  C_cum <- c(0,cumsum(C))
  
  #------------------------ Step 0: Prepare for outputs -----------------
  Z_output <- NULL
  
  P_C_J_D_output <- rep(list(NULL),D)
  
  P_output <- NULL
  alpha_output <- c()
  alpha_0_output <- c()
  mu_star_1_J_output <- NULL
  Sigma_star_1_J_output <- NULL
  
  Q_J_D_output <- NULL
  Xi_C_D_output <- NULL
  x_star_J_D_output <- NULL
  sigma_star_2_J_D_output <- NULL
  
  U_C_J_D_output <- rep(list(NULL),D)
  
  r_J_output <- NULL
  s_2_output <- c()
  h_J_output <- NULL
  m_2_output <- c()
  
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
  
  
  # ---------- Initial kernel and hyper-paramerters ---------
  # Compute the mean of x within each cluster in each dataset to initialize x_star_J_D
  # There may be empty clusters
  # Fill the non-empty values into the initialization matrix, impute NA values with dataset-specific mean
  x_star_J_D_list <- lapply(1:D, function(d) {
    
    val <- rep(NA, J)
    # Mean of x within each cluster of each dataset
    x_star_J <- tapply(x[[d]], Z_new[[d]], mean)
    val[as.numeric(names(x_star_J))] <- x_star_J
    val[is.na(val)] <- mean(x[[d]])
    
    return(val)
  })
  
  # matrix of [J, D]
  x_star_J_D_new <- do.call(cbind,x_star_J_D_list)
  rm(x_star_J_D_list)
  
  # Initials for r_J as an empirical mean from x_star_J_D_new
  r_J_new <- apply(x_star_J_D_new, 1, mean)
  
  # Initials for s^2
  s_2_new <- mean(apply(x_star_J_D_new, 1, var))
  
  # Initials for sigma_star_2_J_D
  sigma_star_2_J_D_list <- lapply(1:D, function(d) {
    
    val <- rep(NA, J)
    # Mean of x within each cluster of each dataset
    sigma_star_J <- tapply(x[[d]], Z_new[[d]], var)
    val[as.numeric(names(sigma_star_J))] <- sigma_star_J
    val[is.na(val)] <- var(x[[d]])
    
    return(val)
  })
  
  # matrix of [J, D]
  sigma_star_2_J_D_new <- do.call(cbind,sigma_star_2_J_D_list)
  rm(sigma_star_2_J_D_list)
  
  # Initials for h_J from log(sigma_star_2_J_D)
  h_J_new <- apply(log(sigma_star_2_J_D_new), 1, mean)
  
  # Initials for m^2
  m_2_new <- mean(apply(log(sigma_star_2_J_D_new), 1, var))
  
  # ----------- Initial P and q ---------------
  # Component probabilities p_j, add 1 to avoid zero p if cluster is empty
  P_new <- sapply(1:J, function(j) mean(unlist(Z_new)==j))
  if(any(P_new==0)){
    # add 1 to avoid zero p if cluster is empty
    P_new <- sapply(1:J, function(j) sum(unlist(Z_new)==j)+1)/(sum(C)+J)
  }
  
  # Dataset-specific vector q_j_d
  Q_J_D_new <- do.call(cbind, lapply(1:D, function(d) {
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
  #cluster-specific parameters in multivariate normal likelihood
  # mu_star_1_J_new <- do.call(rbind, lapply(1:J, function(j) {
  #   
  # }))
  # Sigma_star_1_J_new <- t(matrix(phi.estimate, nrow=G, ncol=J))
  
  # ----------- Initial alpha, alpha0 ---------
  
  alpha_new <- 1
  alpha_0_new <- 1
  
  # ---------- Initial Xi and U --------
  
  Xi_C_D_new <- Xi_C_D_update(Q_J_D = Q_J_D_new, C = C, x = x, x_star_J_D = x_star_J_D_new,
                              sigma_star_2_J_D = sigma_star_2_J_D_new)
  
  U_C_J_D_new <- U_C_J_D_update(Xi_C_D = Xi_C_D_new, Q_J_D = Q_J_D_new, C = C, x = x,
                                x_star_J_D = x_star_J_D_new, sigma_star_2_J_D = sigma_star_2_J_D_new,
                                rbf = rbf)
  
  #----------------------- Step 2:  Acceptance probability ----------------------------
  
  acceptance_count_avg <- data.frame(P_accept = rep(0,niter), alpha_accept = rep(0,niter),
                                     alpha_0_accept = rep(0,niter), sigma_accept = rep(0, niter))
  
  # Total count of acceptance
  P_count <- 0; alpha_count <- 0; alpha_0_count <- 0;  sigma_star_2_count <- 0
  
  #----------------------- Step 3: Prepare for the adaptive covariance update -----------------------------
  
  
  # 0) For sigma_star_2_J_D, matrix
  # X = -log(1/sigma_star - 1/upper)
  mean_X_sigma_star_2_new <- log(sigma_star_2_J_D_new)
  M_2_sigma_star_2_new <- matrix(0, nrow = J, ncol = D)
  variance_sigma_star_2_new <- matrix(0, nrow = J, ncol = D)
  
  for (d in 1:D) {
    for (j in 1:J) {
      ind <- c(1:C[d])[-log(U_C_J_D_new[[d]][,j])<Xi_C_D_new[[d]]*Q_J_D_new[j,d]]
      if(length(ind)!=0) {
        uppers <- -(x[[d]]-x_star_J_D_new[j,d])^2/2/(log(-log(U_C_J_D_new[[d]][,j]))-log(Xi_C_D_new[[d]])-log(Q_J_D_new[j,d]))
        upper <- min(uppers[ind])
        mean_X_sigma_star_2_new[j,d] <- -log(1/sigma_star_2_J_D_new[j,d] - 1/upper)
      }
    }
  }
  
  # 1) For Component probabilities
  sd_P_new <- 0.001
  #X_n (n: indicate iterations)
  mean_X_component_new <- log(matrix(P_new[1:(J-1)]/P_new[J], nrow = 1)) #1x(J-1)
  #tilde(S)_n = tilde(S)_(n-1)+t(X_n)%*%X_n (for n=1, tilde(S)_(n-1)=0)
  tilde_s_component_new <- Matrix::crossprod(mean_X_component_new)
  # t(mean_X_component_new)%*%mean_X_component_new
  # matrix(P_initial[1:J-1]/P_initial[J], ncol = 1)%*%
  # matrix(P_initial[1:J-1]/P_initial[J], nrow = 1)
  #At 1st iteration, the covariance based on the intial values are 0
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
    
    # 1) ------ Update latent/auxiliary variables Xi_C_D ---------
    Xi_C_D_new <- Xi_C_D_update(Q_J_D = Q_J_D_new, C = C, x = x, x_star_J_D = x_star_J_D_new,
                                sigma_star_2_J_D = sigma_star_2_J_D_new)
    
    # 2) ------- Update dataset-specific vector q_j_d --------
    Q_J_D_new <- Q_J_D_update(Z = Z_new, alpha = alpha_new, P = P_new, Xi_C_D = Xi_C_D_new,
                              x = x, x_star_J_D = x_star_J_D_new, 
                              sigma_star_2_J_D = sigma_star_2_J_D_new)
    
    # 3) ------- Update the allocation variable --------
    if(is.null(Z_fix)){
      Z_new <- allocation_variables_update(Y = Y, x = x, mu_star_1_J = mu_star_1_J_new, Sigma_star_1_J = Sigma_star_1_J_new,
                                           Q_J_D = Q_J_D_new, x_star_J_D = x_star_J_D_new,
                                           sigma_star_2_J_D = sigma_star_2_J_D_new)
    }
    
    # 4) ------- Update kernel parameters ----------
    # 4-1) ------ latent variable U_C_J_D ---------
    U_C_J_D_new <- U_C_J_D_update(Xi_C_D = Xi_C_D_new, Q_J_D = Q_J_D_new, C = C, x = x,
                                  x_star_J_D = x_star_J_D_new, sigma_star_2_J_D = sigma_star_2_J_D_new,
                                  rbf = rbf)
    
    
    # 4-2) ------- x_star_J_D -------
    x_star_J_D_new <- x_star_J_D_update(r_J = r_J_new, s_2 = s_2_new, Z = Z_new, x = x, C = C,
                                        sigma_star_2_J_D = sigma_star_2_J_D_new, U_C_J_D = U_C_J_D_new,
                                        Xi_C_D = Xi_C_D_new, Q_J_D = Q_J_D_new)
    
    # update U again
    U_C_J_D_new <- U_C_J_D_update(Xi_C_D = Xi_C_D_new, Q_J_D = Q_J_D_new, C = C, x = x,
                                  x_star_J_D = x_star_J_D_new, sigma_star_2_J_D = sigma_star_2_J_D_new,
                                  rbf = rbf)
    
    
    
    # 4-3) ------- sigma_star_2_J_D -------
    sigma_star_output <- sigma_star_2_J_D_update(sigma_star_2_J_D = sigma_star_2_J_D_new,
                                                 h_J = h_J_new, m_2 = m_2_new, Z = Z_new, x = x, C = C,
                                                 x_star_J_D = x_star_J_D_new, U_C_J_D = U_C_J_D_new,
                                                 Xi_C_D = Xi_C_D_new, Q_J_D = Q_J_D_new,
                                                 X_mean = mean_X_sigma_star_2_new, M_2 = M_2_sigma_star_2_new,
                                                 variance = variance_sigma_star_2_new, iter_num = iter,
                                                 MH.variance = MH.variance)
    
    # print('sigma_star_2_J_D updated')
    
    sigma_star_2_J_D_new <- sigma_star_output$sigma_star_2_J_D_new
    mean_X_sigma_star_2_new <- sigma_star_output$X_mean_new
    M_2_sigma_star_2_new <- sigma_star_output$M_2_new
    variance_sigma_star_2_new <- sigma_star_output$variance_new
    sigma_star_2_count <- sigma_star_2_count + sigma_star_output$accept
    acceptance_count_avg$sigma_accept[iter-1] <- sigma_star_2_count/((iter-1)*J*D)
    
    # Not necessary, but update P_C_J_D
    P_C_J_D_new <- lapply(1:D, function(d) {
      temp <- t(sapply(x[[d]],function(xx) {
        # log-scale
        LP <- log(Q_J_D_new[,d])-(xx-x_star_J_D_new[,d])^2/2/sigma_star_2_J_D_new[,d]
        # rbf(time,x_star_J_D_new[,d],sigma_star_2_J_D_new[,d])
        nc <- -max(LP)
        P <- exp(LP+nc)/sum(exp(LP+nc))
        return(P)
      }))
    })
    
    # 5) ------ Update hyper parameters in priors for x_star, sigma_star_2 ------
    # 5-1) r_J
    r_J_new <- r_J_update(x_star_J_D = x_star_J_D_new, mu_r = mu_r, sigma_r = sigma_r, s_2 = s_2_new)
    
    # 5-2) s_2
    s_2_new <- s_2_update(x_star_J_D = x_star_J_D_new, eta_1 = eta_1, eta_2 = eta_2, r_J = r_J_new)
    
    # 5-3) h_J
    h_J_new <- h_J_update(sigma_star_2_J_D = sigma_star_2_J_D_new, mu_h = mu_h, sigma_h = sigma_h, m_2 = m_2_new)
    
    # 5-4) m_2
    m_2_new <- m_2_update(sigma_star_2_J_D = sigma_star_2_J_D_new, kappa_1 = kappa_1, kappa_2 = kappa_2, h_J = h_J_new)
    
    
    # 6) ------ Update the component probabilities P -------
    component_output <- component_probabilities_update(P = P_new, Q_J_D = Q_J_D_new, alpha_0 = alpha_0_new,
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
    
    # 7) ------- Update alpha -------------
    alpha_output_sim <- alpha_update(Q_J_D = Q_J_D_new, P = P_new, alpha = alpha_new,
                                     X_mean = mean_X_alpha_new, M_2 = M_2_alpha_new, 
                                     variance = variance_alpha_new, iter_num = iter,
                                     MH.variance = MH.variance)
    
    # print('alpha updated')
    alpha_new <- alpha_output_sim$alpha_new
    mean_X_alpha_new <- alpha_output_sim$X_mean_new
    M_2_alpha_new <- alpha_output_sim$M_2_new
    variance_alpha_new <- alpha_output_sim$variance_new
    alpha_count <- alpha_count + alpha_output_sim$accept # cumulative count
    acceptance_count_avg$alpha_accept[iter-1] <- alpha_count/(iter-1) # update the acceptance probability
    # acceptance_prob$alpha_prob[iter-1] <- alpha_output_sim$accept_prob
    
    # 8) ------- Update alpha_0 --------
    alpha_0_output_sim <- alpha_0_update(P = P_new, alpha_0 = alpha_0_new, 
                                         X_mean = mean_X_alpha_0_new,
                                         M_2 = M_2_alpha_0_new, 
                                         variance = variance_alpha_0_new, iter_num = iter,
                                         MH.variance = MH.variance)
    # print('alpha_0 updated')
    
    alpha_0_new <- alpha_0_output_sim$alpha_0_new
    mean_X_alpha_0_new <- alpha_0_output_sim$X_mean_new
    M_2_alpha_0_new <- alpha_0_output_sim$M_2_new
    variance_alpha_0_new <- alpha_0_output_sim$variance_new
    alpha_0_count <- alpha_0_count + alpha_0_output_sim$accept # cumulative count
    acceptance_count_avg$alpha_0_accept[iter-1] <- alpha_0_count/(iter-1) # update the acceptance probability
    # acceptance_prob$alpha_0_prob[iter-1] <- alpha_0_output_sim$accept_prob
    
    
    #-------------------------- Step 5: Update simulated values ------------------------
    if(update == TRUE){
      if(is.null(Z_fix)){
        Z_output[[output_index]] <- Z_new
      }
      
      for(d in 1:D){
        P_C_J_D_output[[d]][[output_index]] <- P_C_J_D_new[[d]]
      }
      
      P_output[[output_index]] <- P_new
      alpha_output[output_index] <- alpha_new
      alpha_0_output[output_index] <- alpha_0_new
      
      mu_star_1_J_output[[output_index]] <- mu_star_1_J_new
      Sigma_star_1_J_output[[output_index]] <- Sigma_star_1_J_new
      
      Q_J_D_output[[output_index]] <- Q_J_D_new
      # Xi_C_D_output[[output_index]] <- Xi_C_D_new
      
      x_star_J_D_output[[output_index]] <- x_star_J_D_new
      sigma_star_2_J_D_output[[output_index]] <- sigma_star_2_J_D_new
      
      # for(d in 1:D){
      #   U_C_J_D_output[[d]][[output_index]] <- U_C_J_D_new[[d]]
      # }
      
      r_J_output[[output_index]] <- as.vector(unname(r_J_new))
      s_2_output[output_index] <- s_2_new
      
      h_J_output[[output_index]] <- as.vector(unname(h_J_new)) 
      m_2_output[output_index] <- m_2_new
      
    }
    
    if((iter-1) %% save_frequency == 0 && auto.save == TRUE){
      
      end_time <- Sys.time()
      diff_time <- difftime(end_time,start_time_mcmc)
      
      if(save.only.z){
        my_list <- list('Z_output' = Z_output, 'acceptance_count_avg' = acceptance_count_avg, 
                        'Z_new' = Z_new, 'P_C_J_D_new' = P_C_J_D_new,
                        'x_star_J_D_new' = x_star_J_D_new, 'r_J_new' = r_J_new, 's_2_new' = s_2_new, 'h_J_new' = h_J_new, 'm_2_new' = m_2_new,
                        'Q_J_D_new' = Q_J_D_new, 
                        'Xi_C_D_new' = Xi_C_D_new, 'U_C_J_D_new' = U_C_J_D_new,
                        'P_new' = P_new, 'tilde_s_component_new' = tilde_s_component_new, 'mean_X_component_new' = mean_X_component_new, 
                        'covariance_component_new' = covariance_component_new, 'P_count' = P_count, 
                        'sigma_star_2_J_D_new' = sigma_star_2_J_D_new, 'mean_X_sigma_star_2_new' = mean_X_sigma_star_2_new,
                        'M_2_sigma_star_2_new' = M_2_sigma_star_2_new, 'variance_sigma_star_2_new' = variance_sigma_star_2_new,
                        'sigma_star_2_count' = sigma_star_2_count,
                        'alpha_new' = alpha_new, 'mean_X_alpha_new' = mean_X_alpha_new, 'M_2_alpha_new' = M_2_alpha_new,
                        'variance_alpha_new' = variance_alpha_new, 'alpha_count' = alpha_count, 
                        'alpha_0_new' = alpha_0_new, 'mean_X_alpha_0_new' = mean_X_alpha_0_new, 'M_2_alpha_0_new' = M_2_alpha_0_new, 
                        'variance_alpha_0_new' = variance_alpha_0_new, 'alpha_0_count' = alpha_0_count, 
                        'output_index' = output_index,
                        'sd_P_new' = sd_P_new,
                        'mu0' = mu0, 'k0' = k0, 'omega0' = omega0, 'Phi0' = Phi0,
                        'time' = diff_time)
      }else{
        my_list <- list('Z_output' = Z_output, 'P_C_J_D_output' = P_C_J_D_output, 'P_output' = P_output, 
                        'alpha_output' = alpha_output, 'alpha_0_output' = alpha_0_output, 
                        'mu_star_1_J_output' = mu_star_1_J_output, 'Sigma_star_1_J_output' = Sigma_star_1_J_output,
                        'Q_J_D_output' = Q_J_D_output, 'x_star_J_D_output' = x_star_J_D_output, 'sigma_star_2_J_D_output' = sigma_star_2_J_D_output, 
                        'r_J_output' = r_J_output, 's_2_output' = s_2_output, 'h_J_output' = h_J_output, 'm_2_output' = m_2_output,
                        'acceptance_count_avg' = acceptance_count_avg, 
                        'Z_new' = Z_new, 'P_C_J_D_new' = P_C_J_D_new, 'mu_star_1_J_new' = mu_star_1_J_new, 'Sigma_star_1_J_new' = Sigma_star_1_J_new,
                        'x_star_J_D_new' = x_star_J_D_new, 'r_J_new' = r_J_new, 's_2_new' = s_2_new, 'h_J_new' = h_J_new, 'm_2_new' = m_2_new,
                        'Q_J_D_new' = Q_J_D_new, 
                        'Xi_C_D_new' = Xi_C_D_new, 'U_C_J_D_new' = U_C_J_D_new,
                        'P_new' = P_new, 'tilde_s_component_new' = tilde_s_component_new, 'mean_X_component_new' = mean_X_component_new, 
                        'covariance_component_new' = covariance_component_new, 'P_count' = P_count, 
                        'sigma_star_2_J_D_new' = sigma_star_2_J_D_new, 'mean_X_sigma_star_2_new' = mean_X_sigma_star_2_new,
                        'M_2_sigma_star_2_new' = M_2_sigma_star_2_new, 'variance_sigma_star_2_new' = variance_sigma_star_2_new,
                        'sigma_star_2_count' = sigma_star_2_count,
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
                    'Z_new' = Z_new, 'P_C_J_D_new' = P_C_J_D_new,
                    'x_star_J_D_new' = x_star_J_D_new, 'r_J_new' = r_J_new, 's_2_new' = s_2_new, 'h_J_new' = h_J_new, 'm_2_new' = m_2_new,
                    'Q_J_D_new' = Q_J_D_new, 
                    'Xi_C_D_new' = Xi_C_D_new, 'U_C_J_D_new' = U_C_J_D_new,
                    'P_new' = P_new, 'tilde_s_component_new' = tilde_s_component_new, 'mean_X_component_new' = mean_X_component_new, 
                    'covariance_component_new' = covariance_component_new, 'P_count' = P_count, 
                    'sigma_star_2_J_D_new' = sigma_star_2_J_D_new, 'mean_X_sigma_star_2_new' = mean_X_sigma_star_2_new,
                    'M_2_sigma_star_2_new' = M_2_sigma_star_2_new, 'variance_sigma_star_2_new' = variance_sigma_star_2_new,
                    'sigma_star_2_count' = sigma_star_2_count,
                    'alpha_new' = alpha_new, 'mean_X_alpha_new' = mean_X_alpha_new, 'M_2_alpha_new' = M_2_alpha_new,
                    'variance_alpha_new' = variance_alpha_new, 'alpha_count' = alpha_count, 
                    'alpha_0_new' = alpha_0_new, 'mean_X_alpha_0_new' = mean_X_alpha_0_new, 'M_2_alpha_0_new' = M_2_alpha_0_new, 
                    'variance_alpha_0_new' = variance_alpha_0_new, 'alpha_0_count' = alpha_0_count, 
                    'output_index' = output_index,
                    'sd_P_new' = sd_P_new,
                    'mu0' = mu0, 'k0' = k0, 'omega0' = omega0, 'Phi0' = Phi0,
                    'time' = diff_time)
  }else{
    my_list <- list('Z_output' = Z_output, 'P_C_J_D_output' = P_C_J_D_output, 'P_output' = P_output, 
                    'alpha_output' = alpha_output, 'alpha_0_output' = alpha_0_output, 
                    'mu_star_1_J_output' = mu_star_1_J_output, 'Sigma_star_1_J_output' = Sigma_star_1_J_output,
                    'Q_J_D_output' = Q_J_D_output, 'x_star_J_D_output' = x_star_J_D_output, 'sigma_star_2_J_D_output' = sigma_star_2_J_D_output, 
                    'r_J_output' = r_J_output, 's_2_output' = s_2_output, 'h_J_output' = h_J_output, 'm_2_output' = m_2_output,
                    'acceptance_count_avg' = acceptance_count_avg, 
                    'Z_new' = Z_new, 'P_C_J_D_new' = P_C_J_D_new, 'mu_star_1_J_new' = mu_star_1_J_new, 'Sigma_star_1_J_new' = Sigma_star_1_J_new,
                    'x_star_J_D_new' = x_star_J_D_new, 'r_J_new' = r_J_new, 's_2_new' = s_2_new, 'h_J_new' = h_J_new, 'm_2_new' = m_2_new,
                    'Q_J_D_new' = Q_J_D_new, 
                    'Xi_C_D_new' = Xi_C_D_new, 'U_C_J_D_new' = U_C_J_D_new,
                    'P_new' = P_new, 'tilde_s_component_new' = tilde_s_component_new, 'mean_X_component_new' = mean_X_component_new, 
                    'covariance_component_new' = covariance_component_new, 'P_count' = P_count, 
                    'sigma_star_2_J_D_new' = sigma_star_2_J_D_new, 'mean_X_sigma_star_2_new' = mean_X_sigma_star_2_new,
                    'M_2_sigma_star_2_new' = M_2_sigma_star_2_new, 'variance_sigma_star_2_new' = variance_sigma_star_2_new,
                    'sigma_star_2_count' = sigma_star_2_count,
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

