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

##--------------- Simulation of latent Xi ----------------

Xi_C_D_update <- function(Q_J_D, C, x, x_star_J_D, sigma_star_2_J_D){
  J <- nrow(Q_J_D)
  D <- ncol(Q_J_D)
  
  #set up the list to save updated values
  Xi <- NULL
  
  #draw from Gamma (full conditional distribution)
  for (d in 1:D) {
    
    rates <- unlist(lapply(1:C[d], function(c) {
      log_rbf_value <- -(x[[d]][c]-x_star_J_D[,d])^2/2/sigma_star_2_J_D[,d]
      log_temp <- log(Q_J_D[,d])+log_rbf_value
      log_K <- max(log_temp)
      
      return(exp(log_K)*sum(exp(log_temp-log_K)))
    }))
    
    Xi[[d]] <- rgamma(C[d], shape = 1, rate = rates)
  }
  
  return(Xi)
}

##--------------- Simulation of dataset-specific vector q_j_d ----------------

Q_J_D_update <- function(Z, alpha, P, Xi_C_D, x, x_star_J_D, sigma_star_2_J_D){
  J <- nrow(x_star_J_D)
  D <- ncol(x_star_J_D)
  
  #set up the matrix to save updated values
  Q <- matrix(NA, nrow = J, ncol = D)
  
  for (d in 1:D) {
    shape_rate_param <- do.call(cbind,lapply(1:J, function(j) {
      shape <- sum(Z[[d]]==j)+alpha*P[j]
      
      log_rbf_value <- -(x[[d]]-x_star_J_D[j,d])^2/2/sigma_star_2_J_D[j,d]
      log_temp <- log(Xi_C_D[[d]])+log_rbf_value
      log_K <- max(log_temp)
      rate <- 1+exp(log_K)*sum(exp(log_temp-log_K))
      # sum(Xi_C_D[[d]]*rbf(t[[d]], x_star_J_D[j,d], sigma_star_2_J_D[j,d]))
      
      return(c(shape, rate))
    }))
    
    Q[,d] <- rgamma(J, shape = shape_rate_param[1,], rate = shape_rate_param[2,])
    
    if(any(Q[,d]==0)) {
      Q[Q[,d]==0,d] <- shape_rate_param[1,Q[,d]==0]/shape_rate_param[2,Q[,d]==0]
    }
  }
  
  return(Q)
}

##--------------- Simulation of allocations ----------------

allocation_variables_update <- function(Y, x, mu_star_1_J, Sigma_star_1_J, Q_J_D, 
                                        x_star_J_D, sigma_star_2_J_D){
  
  D <- length(Y)
  C <- unlist(lapply(Y, nrow))
  J <- nrow(mu_star_1_J)
  
  # Set up the list to save updated values
  Z <- NULL
  
  for(d in 1:D){
    
    # n *J
    Q_mat <- matrix(rep(Q_J_D[,d],C[d]),nrow=C[d],byrow = TRUE)
    x_star_mat <- matrix(rep(x_star_J_D[,d],C[d]),nrow=C[d],byrow = TRUE)
    sigma_star_2_mat <- matrix(rep(sigma_star_2_J_D[,d],C[d]),nrow=C[d],byrow = TRUE)
    x_mat <- matrix(rep(x[[d]],J),nrow=C[d],byrow = FALSE)
    
    # n * J
    LP1 <- -(x_mat-x_star_mat)^2/2/sigma_star_2_mat+log(Q_mat)
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

##----------- Simulation of latent variables U_C_J_D for updating kernel parameters----------

U_C_J_D_update <- function(Xi_C_D, Q_J_D, C, x, x_star_J_D, sigma_star_2_J_D, rbf){
  J <- nrow(Q_J_D)
  D <- ncol(Q_J_D)
  
  # Set up the list to save updated values
  U_C_J_D <- NULL
  
  for (d in 1:D) {
    U_C_J_D[[d]] <- do.call(cbind, lapply(1:J, function(j) {
      K_C_J_D <- exp(-Xi_C_D[[d]]*Q_J_D[j,d]*rbf(x[[d]],x_star_J_D[j,d],sigma_star_2_J_D[j,d]))
      val <- runif(C[d], min = 0, max = K_C_J_D) 
      
      return(val)
    }))
  }
  
  return(U_C_J_D)
  
}


##----------- Simulation of kernel parameters x_star_J_D-------------

x_star_J_D_update <- function(r_J, s_2, Z, x, C, sigma_star_2_J_D, U_C_J_D, Xi_C_D, Q_J_D){
  J <- nrow(Q_J_D)
  D <- ncol(Q_J_D)
  
  # Set up the matrix to save updated values
  x_star <- matrix(NA, nrow = J, ncol = D)
  
  # Find the truncation region
  
  for (d in 1:D) {
    for (j in 1:J){
      #find which cells to truncate regions
      ind <- c(1:C[d])[-log(U_C_J_D[[d]][,j])<Xi_C_D[[d]]*Q_J_D[j,d]]
      if(length(ind)==0) {
        truncate <- FALSE
      }else{
        #if truncate, find truncated regions
        truncate <- TRUE
        i_all_complement <- lapply(ind, function(c) {
          temp <- -2*sigma_star_2_J_D[j,d]*(log(-log(U_C_J_D[[d]][c,j]))-log(Xi_C_D[[d]][c])-log(Q_J_D[j,d]))
          
          i_complement <- matrix(c(x[[d]][c]-sqrt(temp),x[[d]][c]+sqrt(temp)),nrow=1)
          return(i_complement)
        })
        i_all <- intervals::interval_complement(intervals::interval_union(intervals::Intervals(do.call(rbind,i_all_complement))))@.Data
      }
      
      
      #parameters in posterior distribution (Normal)
      ind1 <- (Z[[d]]==j)
      sum_x_j_d <- sum(x[[d]][ind1])
      N_j_d <- sum(ind1)
      r_j_hat <- (r_J[j]*sigma_star_2_J_D[j,d]+s_2*sum_x_j_d)/(sigma_star_2_J_D[j,d]+N_j_d*s_2)
      s_2_hat <- s_2*sigma_star_2_J_D[j,d]/(sigma_star_2_J_D[j,d]+N_j_d*s_2)
      
      if(!truncate){
        # No truncation, draw from the Normal posterior directly
        x_star[j,d] <- rnorm(1,mean = r_j_hat, sd = sqrt(s_2_hat))
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
        x_star[j,d] <- truncnorm::rtruncnorm(1,a=i_chosen[1],
                                             b=i_chosen[2],
                                             mean = r_j_hat, sd = sqrt(s_2_hat))
      }
    }
  }
  
  return(x_star)
}

##----------- Simulation of kernel parameters sigma_star_2_J_D-------------

sigma_star_2_log_prob <- function(sigma_star_2_j_d, m_2, h_j, x_star_j_d, x_sub){
  lprod <- -log(sigma_star_2_j_d)-(log(sigma_star_2_j_d)-h_j)^2/2/m_2-sum((x_sub-x_star_j_d)^2)/2/sigma_star_2_j_d
  return(lprod)
}

sigma_star_2_J_D_update <- function(sigma_star_2_J_D, h_J, m_2, Z, x, C, x_star_J_D, U_C_J_D, Xi_C_D, Q_J_D,
                                    X_mean, M_2, variance, iter_num, MH.variance){
  
  J <- nrow(Q_J_D)
  D <- ncol(Q_J_D)
  
  #all are JxD matrices
  sigma_star_2_J_D_old <- sigma_star_2_J_D
  X_mean_old <- X_mean
  M_2_old <- M_2
  variance_old <- variance
  
  # save updated covariance matrices
  X_new <- matrix(NA, nrow=J, ncol = D)
  X_mean_new <- matrix(NA, nrow=J, ncol = D)
  M_2_new <- matrix(NA, nrow=J, ncol = D)
  variance_new <- matrix(NA, nrow=J, ncol = D)
  
  # save updated sigma_star_2_J_D
  sigma_star_2_J_D_new <- matrix(NA, nrow=J, ncol = D)
  
  # Defining the dimensions
  n <- iter_num
  
  accept_count <- 0
  
  for (d in 1:D) {
    for (j in 1:J) {
      # compute truncation region if needed
      # find which cells to truncate regions
      ind <- c(1:C[d])[-log(U_C_J_D[[d]][,j])<Xi_C_D[[d]]*Q_J_D[j,d]]
      if(length(ind)==0) {
        truncate <- FALSE
      }else{
        truncate <- TRUE
        uppers <- -(x[[d]]-x_star_J_D[j,d])^2/2/(log(-log(U_C_J_D[[d]][,j]))-log(Xi_C_D[[d]])-log(Q_J_D[j,d]))
        upper <- min(uppers[ind]) #for sigma_star_2
        if(upper<=sigma_star_2_J_D_old[j,d]){
          print(upper)
        }
      }
      
      #for empty clusters in dataset d, sample from the log-normal prior (may be truncated), always accept
      if(sum(Z[[d]]==j)==0) {
        if(!truncate){ #no truncation
          sigma_star_2_J_D_new[j,d] <- rlnorm(1, meanlog = h_J[j], sdlog = sqrt(m_2))
          X_new[j,d] <- log(sigma_star_2_J_D_new[j,d])
          
          accept_count <- accept_count+1
        }else{ #truncate
          # X = -log(1/sigma_star - 1/upper)
          p.upper <- plnorm(upper, meanlog = h_J[j], sdlog = sqrt(m_2))
          u <- runif(1, 0, p.upper)
          sigma_star_2_J_D_new[j,d] <- qlnorm(u, meanlog = h_J[j], sdlog = sqrt(m_2))
          
          X_new[j,d] <- -log(1/sigma_star_2_J_D_new[j,d] - 1/upper)
          
          accept_count <- accept_count+1
        }
      }else{ #non-empty clusters, sample from adaptive MH
        if(!truncate){
          # No truncation, perform adaptive MH directly, propose sigma_star_j_d from log-normal, X=log(sigma_star_2)
          if(n <= 100){
            X_new[j,d] <- rnorm(n = 1, mean = log(sigma_star_2_J_D_old[j,d]), sd = 0.1)
          }else{
            X_new[j,d] <- rnorm(n = 1, mean = log(sigma_star_2_J_D_old[j,d]), sd = sqrt(2.4^2*variance_old[j,d]+2.4^2*MH.variance))
          }
          
          # Transform the new value of X back to new value of sigma_star_2_j_d
          sigma_star_2_J_D_new[j,d] <- exp(X_new[j,d])
          if(is.infinite(sigma_star_2_J_D_new[j,d])) {
            print(paste('No truncate:','cluster',j,'dataset',d,": sigma==Inf"))
            sigma_star_2_J_D_new[j,d] <- .Machine$double.xmax
          }
          if(sigma_star_2_J_D_new[j,d]==0) {
            print(paste('No truncate:','cluster',j,'dataset',d,": sigma==0"))
            sigma_star_2_J_D_new[j,d] <- .Machine$double.eps
          }
          
          # subset of x such that Z[[d]][c]=j
          x_sub <- x[[d]][Z[[d]]==j]
          
          # Compute log acceptance probability
          log_acceptance <- sigma_star_2_log_prob(sigma_star_2_j_d = sigma_star_2_J_D_new[j,d], m_2 = m_2, 
                                                  h_j = h_J[j], x_star_j_d = x_star_J_D[j,d], 
                                                  x_sub = x_sub) - 
            sigma_star_2_log_prob(sigma_star_2_j_d = sigma_star_2_J_D_old[j,d], m_2, 
                                  h_J[j], x_star_J_D[j,d], x_sub) -
            log(sigma_star_2_J_D_old[j,d]) + log(sigma_star_2_J_D_new[j,d]) 
          acceptance_sigma <- exp(log_acceptance)
          acceptance_sigma <- min(1,acceptance_sigma)
          
          if(is.na(acceptance_sigma)) {
            print(acceptance_sigma)
          }
          
          # Update sigma_star_2_J_D_new, X_new, X_mean_new, variance_new
          outcome <- rbinom(n = 1, size = 1, prob = acceptance_sigma)
          if(is.na(outcome) == TRUE | outcome == 0){
            X_new[j,d] <- log(sigma_star_2_J_D_old[j,d])
            sigma_star_2_J_D_new[j,d] <- sigma_star_2_J_D_old[j,d]
          }else{
            accept_count <- accept_count+1
          }
          
        }else{
          #truncate
          # X = -log(1/sigma_star - 1/upper)
          if(n <= 100){
            X_new[j,d] <- rnorm(n = 1, mean = -log(1/sigma_star_2_J_D_old[j,d] - 1/upper), sd = 0.1) 
          }else{
            X_new[j,d] <- rnorm(n = 1, mean = -log(1/sigma_star_2_J_D_old[j,d] - 1/upper), sd = sqrt(2.4^2*variance_old[j,d]+2.4^2*MH.variance))
          }
          
          # Transform the new value of X back to new value of sigma_star_2_j_d
          sigma_star_2_J_D_new[j,d] <- 1/(exp(-X_new[j,d])+1/upper)
          if(sigma_star_2_J_D_new[j,d]==0) { 
            print(paste('Truncate:','cluster',j,'dataset',d,": sigma==0"))
            sigma_star_2_J_D_new[j,d] <- .Machine$double.eps
          }
          
          # subset of x such that Z[[d]][c]=j
          x_sub <- x[[d]][Z[[d]]==j]
          
          # Compute log acceptance probability
          log_acceptance <- sigma_star_2_log_prob(sigma_star_2_j_d = sigma_star_2_J_D_new[j,d], m_2 = m_2, 
                                                  h_j = h_J[j], x_star_j_d = x_star_J_D[j,d], 
                                                  x_sub = x_sub) - 
            sigma_star_2_log_prob(sigma_star_2_j_d = sigma_star_2_J_D_old[j,d], m_2, 
                                  h_J[j], x_star_J_D[j,d], x_sub) - 
            log(sigma_star_2_J_D_old[j,d]) - log(upper-sigma_star_2_J_D_old[j,d]) + 
            log(sigma_star_2_J_D_new[j,d]) + log(upper-sigma_star_2_J_D_new[j,d])
          
          acceptance_sigma <- exp(log_acceptance)
          acceptance_sigma <- min(1,acceptance_sigma)
          
          if(is.na(acceptance_sigma)) {
            print(acceptance_sigma)
          }
          
          # Update sigma_star_2_J_D_new, X_new, X_mean_new, variance_new
          outcome <- rbinom(n = 1, size = 1, prob = acceptance_sigma)
          if(is.na(outcome) == TRUE | outcome == 0){
            X_new[j,d] <- -log(1/sigma_star_2_J_D_old[j,d] - 1/upper)
            sigma_star_2_J_D_new[j,d] <- sigma_star_2_J_D_old[j,d]
          }else{
            accept_count <- accept_count+1
          }
          
        }#end sampling in truncation case (non-empty clusters)
        
      }#end non-empty clusters
      
      # update variance in adaptive MH
      X_mean_new[j,d] <- (1-1/n)*X_mean_old[j,d] + 1/n*X_new[j,d]
      M_2_new[j,d] <- M_2_old[j,d] + (X_new[j,d]-X_mean_old[j,d])*(X_new[j,d]-X_mean_new[j,d])
      variance_new[j,d] <- 1/(n-1)*M_2_new[j,d]
      
    }#end for j in 1:J
    
    
  }
  
  return(list(sigma_star_2_J_D_new=sigma_star_2_J_D_new, X_mean_new=X_mean_new, M_2_new=M_2_new, 
              variance_new=variance_new, accept=accept_count))
}


##----------- Simulation of hyper parameters r_J (in prior for x_star) ----------------
r_J_update <- function(x_star_J_D, mu_r, sigma_r, s_2){
  J <- nrow(x_star_J_D)
  D <- ncol(x_star_J_D)
  
  # Update
  r_J <- sapply(1:J, function(j) {
    mu_r_hat <- (mu_r*s_2+sigma_r^2*sum(x_star_J_D[j,]))/(s_2+D*sigma_r^2)
    sigma_r_2_hat <- sigma_r^2*s_2/(s_2+D*sigma_r^2)
    return(rnorm(1, mean = mu_r_hat, sd = sqrt(sigma_r_2_hat)))
  })
}

##----------- Simulation of hyper parameters s_2 (in prior for x_star) ----------------
s_2_update <- function(x_star_J_D, eta_1, eta_2, r_J){
  J <- nrow(x_star_J_D)
  D <- ncol(x_star_J_D)
  
  # Update
  shape <- J*D/2+eta_1
  
  temp <- sapply(1:D, function(d) {
    return(x_star_J_D[,d]-r_J)
  })
  
  rate <- eta_2+sum(temp^2)/2
  
  s_2 <- rinvgamma(1,alpha = shape,beta = rate)
  
  return(s_2)
}

##----------- Simulation of hyper parameters h_J (in prior for sigma_star_2) ----------------
h_J_update <- function(sigma_star_2_J_D, mu_h, sigma_h, m_2){
  J <- nrow(sigma_star_2_J_D)
  D <- ncol(sigma_star_2_J_D)
  
  # Update
  h_J <- sapply(1:J, function(j) {
    mu_h_hat <- (mu_h*m_2+sigma_h^2*sum(log(sigma_star_2_J_D[j,])))/(m_2+D*sigma_h^2)
    sigma_h_2_hat <- sigma_h^2*m_2/(m_2+D*sigma_h^2)
    return(rnorm(1, mean = mu_h_hat, sd = sqrt(sigma_h_2_hat)))
  })
}

##----------- Simulation of hyper parameters m_2 (in prior for sigma_star_2) ----------------
m_2_update <- function(sigma_star_2_J_D, kappa_1, kappa_2, h_J){
  J <- nrow(sigma_star_2_J_D)
  D <- ncol(sigma_star_2_J_D)
  
  # Update
  shape <- J*D/2+kappa_1
  
  temp <- sapply(1:D, function(d) {
    return(log(sigma_star_2_J_D[,d])-h_J)
  })
  
  rate <- kappa_2+sum(temp^2)/2
  
  m_2 <- rinvgamma(1,alpha = shape,beta = rate)
  
  return(m_2)
}


##----------- Simulation of Component probabilities P----------------

component_log_prob <- function(P, Q_J_D, alpha_0, alpha){
  
  J <- nrow(Q_J_D)
  D <- ncol(Q_J_D)
  
  lprod <- sum((alpha_0/J-1)*log(P))
  for(d in 1:D){
    lprod <- lprod + sum(alpha*P*log(Q_J_D[,d])-lgamma(alpha*P))
  }
  
  # Returning outputs
  return(lprod)
}

# 2) Simulation
component_probabilities_update <- function(P, Q_J_D, alpha_0, alpha, covariance,
                                           mean_x, tilde_s, iter_num, sd_P, MH.variance){
  
  J <- nrow(Q_J_D)
  D <- ncol(Q_J_D)
  
  
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
  log_acceptance <- component_log_prob(P_new, Q_J_D, alpha_0, alpha) -
    component_log_prob(P_old, Q_J_D, alpha_0, alpha) +
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
alpha_log_prob <- function(Q_J_D, P, alpha){
  
  D <- ncol(Q_J_D)
  
  
  # Construct the log-probability
  lprod <- -alpha + sum(vapply(1:D, function(d) {
    sum(alpha*P*log(Q_J_D[,d])-lgamma(alpha*P))
  }, FUN.VALUE = numeric(1)))
  
  
  # Return the log-probability
  return(lprod)
}

alpha_update <- function(Q_J_D, P, alpha, X_mean, M_2, variance, iter_num, MH.variance){
  
  D <- ncol(Q_J_D)
  
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
  log_acceptance <- alpha_log_prob(Q_J_D, P, alpha = alpha_new) -
    alpha_log_prob(Q_J_D, P, alpha = alpha_old) +
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
  
  # Define the input value
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
