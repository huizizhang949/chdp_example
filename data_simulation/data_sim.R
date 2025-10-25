# 10 replicates, change rep_ind from 1 to 10
D <- 5; J <- 3

beta_J <- list(
  c(0, 20, 0), #cluster 1
  c(10, -20, 0), #cluster 2
  c(0, 0, 40) #cluster 3
)

set.seed(2345+rep_ind)
coefs <- lapply(1:D, function(d) {
  temp <- lapply(1:J, function(j) {
    
    mvtnorm::rmvnorm(n=1, mean=beta_J[[j]], sigma=5*diag(length(beta_J[[j]])))
    
  })
})

x <- rep(list(seq(0.01,1,length.out=300)),D)
px <- lapply(1:D,function(d) {
  mat=t(sapply(x[[d]],function(xx) {
    x_vec <- c(1,xx,xx^2)
    # log-scale
    LP <- sapply(1:J, function(j) {
      sum(coefs[[d]][[j]]*x_vec)
    })
    nc <- -max(LP)
    P <- exp(LP+nc)/sum(exp(LP+nc))
    
  }))
})

# ---------- simulation of data --------

set.seed(324+rep_ind)
z <- lapply(1:D, function(d) {
  
  apply(px[[d]], 1, function(x) extraDistr::rcat(1, prob=x))
  
  
})
table(unlist(z))


mu_J <- list(
  
  c(0,0),
  c(4,4),
  c(0,4)
  
)

Sigma_J <- list(
  
  diag(1,nrow = length(mu_J[[1]])),
  cbind(c(1,0.5),c(0.5,1)),
  cbind(c(1,-0.5),c(-0.5,1))
  
)

set.seed(434+rep_ind)
Y <- lapply(1:D, function(d) {
  
  temp <- do.call(rbind, lapply(1:length(z[[d]]), function(i) {
    
    mvtnorm::rmvnorm(n=1,mean=mu_J[[z[[d]][i]]],sigma = Sigma_J[[z[[d]][i]]])
    
  }))
  
})

save(Y,file=paste0('data_simulation/Y_',rep_ind,'.RData'))
save(z,file=paste0('data_simulation/z_',rep_ind,'.RData'))
save(px,file=paste0('data_simulation/px_',rep_ind,'.RData'))
