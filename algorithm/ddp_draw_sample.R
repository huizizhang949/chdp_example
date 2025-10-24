##---------------- MCMC Simulation --------------

## --------------- Simulation of mean and Sigma within each component -------
unique_params_update <- function(Y, J, Z, mu0, k0, omega0, Phi0){
  
  loop.result <- lapply(1:J,function(j) {
    
    N_j <- sum(Z==j)
    
    if(N_j!=0){
      
      Y_j <- Y[Z==j,]
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

##--------------- Simulation of latent Xi ----------------

Xi_C_update <- function(Q_J, C, x, x_star_J, sigma_star_2_J){
  J <- length(Q_J)
  
  #draw from Gamma (full conditional distribution)
  
  rates <- unlist(lapply(1:C, function(c) {
    log_rbf_value <- -(x[c]-x_star_J)^2/2/sigma_star_2_J
    log_temp <- log(Q_J)+log_rbf_value
    log_K <- max(log_temp)
    
    return(exp(log_K)*sum(exp(log_temp-log_K)))
  }))
  
  Xi <- rgamma(C, shape = 1, rate = rates)
  
  return(Xi)
}

##--------------- Simulation of dataset-specific vector q_j_d ----------------

Q_J_update <- function(Z, alpha_0, Xi_C, x, x_star_J, sigma_star_2_J){
  J <- length(x_star_J)
  
  shape_rate_param <- do.call(cbind,lapply(1:J, function(j) {
    shape <- sum(Z==j)+alpha_0/J
    
    log_rbf_value <- -(x-x_star_J[j])^2/2/sigma_star_2_J[j]
    log_temp <- log(Xi_C)+log_rbf_value
    log_K <- max(log_temp)
    rate <- 1+exp(log_K)*sum(exp(log_temp-log_K))

    return(c(shape, rate))
  }))
  
  Q <- rgamma(J, shape = shape_rate_param[1,], rate = shape_rate_param[2,])
  
  if(any(Q==0)) {
    Q[Q==0] <- shape_rate_param[1,Q==0]/shape_rate_param[2,Q==0]
  }
  
  
  return(Q)
}

##--------------- Simulation of allocations ----------------

allocation_variables_update <- function(Y, x, mu_star_1_J, Sigma_star_1_J, Q_J, 
                                        x_star_J, sigma_star_2_J){
  
  C <- nrow(Y)
  J <- nrow(mu_star_1_J)
  
  # n *J
  Q_mat <- matrix(rep(Q_J,C),nrow=C,byrow = TRUE)
  x_star_mat <- matrix(rep(x_star_J,C),nrow=C,byrow = TRUE)
  sigma_star_2_mat <- matrix(rep(sigma_star_2_J,C),nrow=C,byrow = TRUE)
  x_mat <- matrix(rep(x,J),nrow=C,byrow = FALSE)
  
  # n * J
  LP1 <- -(x_mat-x_star_mat)^2/2/sigma_star_2_mat+log(Q_mat)
  LP2 <- lapply(1:J,function(j) {
    temp <- mniw::dmNorm(as.matrix(Y), mu=mu_star_1_J[j,], Sigma = Sigma_star_1_J[[j]], log = TRUE)
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
  
  Z <- apply(P, 1, function(x) rcat(1, prob=x))
  
  return(Z)
}

##----------- Simulation of latent variables U_C_J for updating kernel parameters----------

U_C_J_update <- function(Xi_C, Q_J, C, x, x_star_J, sigma_star_2_J, rbf){
  J <- length(Q_J)
  
  U_C_J <- do.call(cbind, lapply(1:J, function(j) {
    K_C_J_D <- exp(-Xi_C*Q_J[j]*rbf(x,x_star_J[j],sigma_star_2_J[j]))
    val <- runif(C, min = 0, max = K_C_J_D) 
    
    return(val)
  }))
  
  return(U_C_J)
  
}


##----------- Simulation of kernel parameters x_star_J-------------

x_star_J_update <- function(r_J, s_2, Z, x, C, sigma_star_2_J, U_C_J, Xi_C, Q_J){
  J <- length(Q_J)
  
  # Set up the matrix to save updated values
  x_star <- rep(NA, J)
  
  # Find the truncation region
  for (j in 1:J){
    #find which cells to truncate regions
    ind <- c(1:C)[-log(U_C_J[,j])<Xi_C*Q_J[j]]
    if(length(ind)==0) {
      truncate <- FALSE
    }else{
      #if truncate, find truncated regions
      truncate <- TRUE
      i_all_complement <- lapply(ind, function(c) {
        temp <- -2*sigma_star_2_J[j]*(log(-log(U_C_J[c,j]))-log(Xi_C[c])-log(Q_J[j]))
        
        i_complement <- matrix(c(x[c]-sqrt(temp),x[c]+sqrt(temp)),nrow=1)
        return(i_complement)
      })
      i_all <- intervals::interval_complement(intervals::interval_union(intervals::Intervals(do.call(rbind,i_all_complement))))@.Data
    }
    
    
    #parameters in posterior distribution (Normal)
    ind1 <- (Z==j)
    sum_x_j <- sum(x[ind1])
    N_j <- sum(ind1)
    r_j_hat <- (r_J*sigma_star_2_J[j]+s_2*sum_x_j)/(sigma_star_2_J[j]+N_j*s_2)
    s_2_hat <- s_2*sigma_star_2_J[j]/(sigma_star_2_J[j]+N_j*s_2)
    
    if(!truncate){
      # No truncation, draw from the Normal posterior directly
      x_star[j] <- rnorm(1,mean = r_j_hat, sd = sqrt(s_2_hat))
    }else{
      N_interval <- nrow(i_all)
      # Compute the probability of lying in each interval
      log_p <- unlist(lapply(1:N_interval,function(l) {
        i <- i_all[l,]
        lower <- i[1]
        upper <- i[2]
        if(is.infinite(lower)){
          val <- pnorm(upper,mean = r_j_hat, sd = sqrt(s_2_hat),log.p = TRUE)
        }else if(is.infinite(upper)){
          val <- pnorm(lower,mean = r_j_hat, sd = sqrt(s_2_hat),log.p = TRUE, lower.tail = FALSE)
        }else{
          lp1 <- pnorm(upper,mean = r_j_hat, sd = sqrt(s_2_hat), log.p = TRUE)
          lp2 <- pnorm(lower,mean = r_j_hat, sd = sqrt(s_2_hat), log.p = TRUE)
          val <- lp1+log(1-exp(lp2-lp1))
        }
        return(val)
      }))
      log_K <- -max(log_p)
      # Select one truncated region (one interval)
      ind2 <- sample(1:N_interval, size = 1, prob = exp(log_p+log_K)) #breakpoint
      i_chosen <- i_all[ind2,]
      x_star[j] <- truncnorm::rtruncnorm(1,a=i_chosen[1],
                                         b=i_chosen[2],
                                         mean = r_j_hat, sd = sqrt(s_2_hat))
    }
  }
  
  
  return(x_star)
}

##----------- Simulation of kernel parameters sigma_star_2_J-------------

sigma_star_2_log_prob <- function(sigma_star_2_j, m_2, h_j, x_star_j, x_sub){
  lprod <- -log(sigma_star_2_j)-(log(sigma_star_2_j)-h_j)^2/2/m_2-sum((x_sub-x_star_j)^2)/2/sigma_star_2_j
  return(lprod)
}

sigma_star_2_J_update <- function(sigma_star_2_J, h_J, m_2, Z, x, C, x_star_J, U_C_J, Xi_C, Q_J,
                                  X_mean, M_2, variance, iter_num, MH.variance){
  
  J <- length(Q_J)
  
  sigma_star_2_J_old <- sigma_star_2_J
  X_mean_old <- X_mean
  M_2_old <- M_2
  variance_old <- variance
  
  # save updated covariance
  X_new <- rep(NA, J)
  X_mean_new <- rep(NA, J)
  M_2_new <- rep(NA, J)
  variance_new <- rep(NA, J)
  
  # save updated sigma_star_2_J
  sigma_star_2_J_new <- rep(NA, J)
  
  # Defining the dimensions
  n <- iter_num
  
  accept_count <- 0
  
  for (j in 1:J) {
    # compute truncation region if needed
    # find which cells to truncate regions
    ind <- c(1:C)[-log(U_C_J[,j])<Xi_C*Q_J[j]]
    if(length(ind)==0) {
      truncate <- FALSE
    }else{
      truncate <- TRUE
      uppers <- -(x-x_star_J[j])^2/2/(log(-log(U_C_J[,j]))-log(Xi_C)-log(Q_J[j]))
      upper <- min(uppers[ind]) #for sigma_star_2
      if(upper<=sigma_star_2_J_old[j]){
        print(upper)
      }
    }
    
    #for empty clusters, sample from the log-normal prior (may be truncated), always accept
    if(sum(Z==j)==0) {
      if(!truncate){ #no truncation
        sigma_star_2_J_new[j] <- rlnorm(1, meanlog = h_J, sdlog = sqrt(m_2))
        X_new[j] <- log(sigma_star_2_J_new[j])
        
        accept_count <- accept_count+1
      }else{ #truncate
        # X = -log(1/sigma_star - 1/upper)
        p.upper <- plnorm(upper, meanlog = h_J, sdlog = sqrt(m_2))
        u <- runif(1, 0, p.upper)
        sigma_star_2_J_new[j] <- qlnorm(u, meanlog = h_J, sdlog = sqrt(m_2))
        
        X_new[j] <- -log(1/sigma_star_2_J_new[j] - 1/upper)
        
        accept_count <- accept_count+1
      }
    }else{ #non-empty clusters, sample from adaptive MH
      if(!truncate){
        # No truncation, perform adaptive MH directly, propose sigma_star_j from log-normal, X=log(sigma_star_2)
        if(n <= 100){
          X_new[j] <- rnorm(n = 1, mean = log(sigma_star_2_J_old[j]), sd = 0.1)
        }else{
          X_new[j] <- rnorm(n = 1, mean = log(sigma_star_2_J_old[j]), sd = sqrt(2.4^2*variance_old[j]+2.4^2*MH.variance))
        }
        
        # Transform the new value of X back to new value of sigma_star_2_j_d
        sigma_star_2_J_new[j] <- exp(X_new[j])
        if(is.infinite(sigma_star_2_J_new[j])) {
          print(paste('No truncate:','cluster',j,'dataset',d,": sigma==Inf"))
          sigma_star_2_J_new[j] <- .Machine$double.xmax
        }
        if(sigma_star_2_J_new[j]==0) {
          print(paste('No truncate:','cluster',j,'dataset',d,": sigma==0"))
          sigma_star_2_J_new[j] <- .Machine$double.eps
        }
        
        # subset of x such that Z[c]=j
        x_sub <- x[Z==j]
        
        # Compute log acceptance probability
        log_acceptance <- sigma_star_2_log_prob(sigma_star_2_j = sigma_star_2_J_new[j], m_2 = m_2, 
                                                h_j = h_J, x_star_j = x_star_J[j], 
                                                x_sub = x_sub) - 
          sigma_star_2_log_prob(sigma_star_2_j = sigma_star_2_J_old[j], m_2, 
                                h_J, x_star_J[j], x_sub) -
          log(sigma_star_2_J_old[j]) + log(sigma_star_2_J_new[j]) 
        acceptance_sigma <- exp(log_acceptance)
        acceptance_sigma <- min(1,acceptance_sigma)
        
        if(is.na(acceptance_sigma)) {
          print(acceptance_sigma)
        }
        
        # Update sigma_star_2_J_new, X_new, X_mean_new, variance_new
        outcome <- rbinom(n = 1, size = 1, prob = acceptance_sigma)
        if(is.na(outcome) == TRUE | outcome == 0){
          X_new[j] <- log(sigma_star_2_J_old[j])
          sigma_star_2_J_new[j] <- sigma_star_2_J_old[j]
        }else{
          accept_count <- accept_count+1
        }
        
      }else{
        #truncate
        # X = -log(1/sigma_star - 1/upper)
        if(n <= 100){
          X_new[j] <- rnorm(n = 1, mean = -log(1/sigma_star_2_J_old[j] - 1/upper), sd = 0.1) 
        }else{
          X_new[j] <- rnorm(n = 1, mean = -log(1/sigma_star_2_J_old[j] - 1/upper), sd = sqrt(2.4^2*variance_old[j]+2.4^2*MH.variance))
        }
        
        # Transform the new value of X back to new value of sigma_star_2_j_d
        sigma_star_2_J_new[j] <- 1/(exp(-X_new[j])+1/upper)
        if(sigma_star_2_J_new[j]==0) { 
          print(paste('Truncate:','cluster',j,'dataset',d,": sigma==0"))
          sigma_star_2_J_new[j] <- .Machine$double.eps
        }
        
        # subset of x such that Z[c]=j
        x_sub <- x[Z==j]
        
        # Compute log acceptance probability
        log_acceptance <- sigma_star_2_log_prob(sigma_star_2_j = sigma_star_2_J_new[j], m_2 = m_2, 
                                                h_j = h_J, x_star_j = x_star_J[j], 
                                                x_sub = x_sub) - 
          sigma_star_2_log_prob(sigma_star_2_j = sigma_star_2_J_old[j], m_2, 
                                h_J, x_star_J[j], x_sub) - 
          log(sigma_star_2_J_old[j]) - log(upper-sigma_star_2_J_old[j]) + 
          log(sigma_star_2_J_new[j]) + log(upper-sigma_star_2_J_new[j])
        
        acceptance_sigma <- exp(log_acceptance)
        acceptance_sigma <- min(1,acceptance_sigma)
        
        if(is.na(acceptance_sigma)) {
          print(acceptance_sigma)
        }
        
        # Update sigma_star_2_J_new, X_new, X_mean_new, variance_new
        outcome <- rbinom(n = 1, size = 1, prob = acceptance_sigma)
        if(is.na(outcome) == TRUE | outcome == 0){
          X_new[j] <- -log(1/sigma_star_2_J_old[j] - 1/upper)
          sigma_star_2_J_new[j] <- sigma_star_2_J_old[j]
        }else{
          accept_count <- accept_count+1
        }
        
      }#end sampling in truncation case (non-empty clusters)
      
    }#end non-empty clusters
    
    # update variance in adaptive MH
    X_mean_new[j] <- (1-1/n)*X_mean_old[j] + 1/n*X_new[j]
    M_2_new[j] <- M_2_old[j] + (X_new[j]-X_mean_old[j])*(X_new[j]-X_mean_new[j])
    variance_new[j] <- 1/(n-1)*M_2_new[j]
    
  }#end for j in 1:J
  
  return(list(sigma_star_2_J_new=sigma_star_2_J_new, X_mean_new=X_mean_new, M_2_new=M_2_new, 
              variance_new=variance_new, accept=accept_count))
}


##----------- Simulation of hyper parameters r_J (in prior for x_star) ----------------
r_J_update <- function(x_star_J, mu_r, sigma_r, s_2){
  J <- length(x_star_J)
  
  # Update
  mu_r_hat <- (mu_r*s_2+sigma_r^2*sum(x_star_J))/(s_2+J*sigma_r^2)
  sigma_r_2_hat <- sigma_r^2*s_2/(s_2+J*sigma_r^2)
  
  return(rnorm(1, mean = mu_r_hat, sd = sqrt(sigma_r_2_hat)))
}

##----------- Simulation of hyper parameters s_2 (in prior for x_star) ----------------
s_2_update <- function(x_star_J, eta_1, eta_2, r_J){
  J <- length(x_star_J)
  
  # Update
  shape <- J/2+eta_1
  
  temp <- x_star_J-r_J
  
  rate <- eta_2+sum(temp^2)/2
  
  s_2 <- rinvgamma(1,alpha = shape,beta = rate)
  
  return(s_2)
}

##----------- Simulation of hyper parameters h_J (in prior for sigma_star_2) ----------------
h_J_update <- function(sigma_star_2_J, mu_h, sigma_h, m_2){
  J <- length(sigma_star_2_J)
  
  # Update
  mu_h_hat <- (mu_h*m_2+sigma_h^2*sum(log(sigma_star_2_J)))/(m_2+J*sigma_h^2)
  sigma_h_2_hat <- sigma_h^2*m_2/(m_2+J*sigma_h^2)
  
  return(rnorm(1, mean = mu_h_hat, sd = sqrt(sigma_h_2_hat)))
}

##----------- Simulation of hyper parameters m_2 (in prior for sigma_star_2) ----------------
m_2_update <- function(sigma_star_2_J, kappa_1, kappa_2, h_J){
  J <- length(sigma_star_2_J)
  
  # Update
  shape <- J/2+kappa_1
  
  temp <- log(sigma_star_2_J)-h_J
  
  rate <- kappa_2+sum(temp^2)/2
  
  m_2 <- rinvgamma(1,alpha = shape,beta = rate)
  
  return(m_2)
}


##---------------------- Simulation of concentration parameter alpha_0 -----------------
alpha_0_log_prob <- function(Q_J,alpha_0){
  
  J <- length(Q_J)
  
  lprob <- -alpha_0 - J*lgamma(alpha_0/J) + sum(alpha_0/J*log(Q_J))
  
  # Return the log-probability
  return(lprob)
}

alpha_0_update <- function(Q_J, alpha_0, X_mean, M_2, variance, iter_num, MH.variance){
  
  J <- length(Q_J)
  
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
  log_acceptance <- alpha_0_log_prob(Q_J, alpha_0 = alpha_0_new) -
    alpha_0_log_prob(Q_J, alpha_0 = alpha_0_old) +
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
