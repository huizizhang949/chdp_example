# Set up packages
pkgs_installed <- rownames(installed.packages())
if(!any(pkgs_installed == 'mcclust.ext')) devtools::install_github('https://github.com/sarawade/mcclust.ext')
if(!any(pkgs_installed == 'mniw')) devtools::install_github('https://github.com/mlysy/mniw.git')

library(coda)
library(mcclust.ext)
library(mclust)
library(scales) #for alpha()
library(RColorBrewer)
library(pbapply)
pkgs <- c('intervals', 'mvtnorm', 'extraDistr', 'Matrix', 'truncnorm', 'MCMCpack')
install.packages(setdiff(pkgs, pkgs_installed))  


setwd('~/chdp_example')

# -------------------------------------
# save posterior samples (both the full and the post-processing step), optimal clustering, ARI
if (!file.exists('result')){
  dir.create('result')
}
# save figures (convergence plot, covariate-dependent weights)
if (!file.exists('fig_result')){
  dir.create('fig_result')
}

# simulate data and save in the folder data_simulation
for (i in 1:10){
  rep_ind=i
  source('data_simulation/data_sim.R')
  rm(list=ls())
}
# ----- function to get optimal clustering ----
# return a vector
get_optimal_cl <- function(chain1, ind1, C){
  z_result1 <- chain1$Z_output[ind1]
  # niter*number of all obs in all datasets (dataset1 on the left, dataset 2 on the right)
  z_matrix1 <- t(sapply(z_result1,function(l) {
    unlist(l)
  }))
  rm(z_result1)
  psm1=comp.psm(z_matrix1)
  VI_test1=minVI(psm1,z_matrix1,method=('all'),include.greedy=FALSE)
  z.est=VI_test1$cl[1,]
  # make sure cluster labels are continuous
  uniq_cl <- unique(z.est)
  J.est <- length(uniq_cl)
  z.updated <- rep(0,sum(C))
  for(j in 1:J.est){
    z.updated[z.est==uniq_cl[j]]=j
  }
  
  return(z.updated)
}

# ----- param and data set up -------
D <- 5; J <- 3
x <- rep(list(seq(0.01,1,length.out=300)),D)
xc <- rep(1:D,each=300) # for ddp with groups, the categorical covariate

C <- sapply(x, length)
C_cum <- c(0,cumsum(C))

empirical_z = TRUE; auto.save = TRUE; save_frequency = 5000
JJ = 8; niter = 15000; burn_in = 12000; thinning = 3; verbose = FALSE

# iteration index for posterior inference
ind1 <- seq(1,1000,by=1)
# for plot
indd1 <- ind1

# ----- run parallel -------
cl <- parallel::makeCluster(5,type='FORK')

start <- Sys.time()
print(start)
ari.all <- pbapply::pblapply(1:5,cl=cl,function(l) {
  
  rep_ind <- l
  load(file=paste0(file='data_simulation/Y_',rep_ind,'.RData'))
  load(file=paste0(file='data_simulation/z_',rep_ind,'.RData'))
  load(file=paste0(file='data_simulation/px_',rep_ind,'.RData'))
  
  
  print(paste0('------- Start CHDP ',rep_ind,' ---------'))
  # need to set local = TRUE to read current local environment
  source('algorithm/chdp_run.R',local=TRUE)
  print(paste0('------- Start HDP ',rep_ind,' ---------'))
  source('algorithm/hdp_run.R',local=TRUE)
  print(paste0('------- Start DDP ',rep_ind,' ---------'))
  source('algorithm/ddp_run.R',local=TRUE)
  print(paste0('------- Start DDP (with group)',rep_ind,' ---------'))
  source('algorithm/ddp_cat_run.R',local=TRUE)
  
  ari.df <- data.frame(chdp=ari.chdp,hdp=ari.hdp,ddp=ari.ddp,ddp_cat=ari.ddp.cat)
  
  return(ari.df)
})
parallel::stopCluster(cl)
end <- Sys.time()
print(end)

df.ari <- do.call(rbind,ari.all)

