# ----- chdp --------
source("algorithm/chdp_mcmc.R")
source("algorithm/chdp_draw_sample.R")


set.seed(4+rep_ind)
chdp.run <- mvn_chdp(Y = Y, x = x, J = JJ, niter = niter, burn_in = burn_in, thinning = thinning, empirical_z = empirical_z,
                     partial.save.name = paste0('result/',rep_ind,'.full.chdp.RData'), save_frequency = save_frequency, 
                     auto.save = auto.save, verbose = verbose)
# get optimal clustering
chain1 <- chdp.run
z.est <- get_optimal_cl(chain1, ind1, C)
J.est <- length(unique(z.est))
## clean up
ari.chdp <- adjustedRandIndex(z.est,unlist(z))
z.chdp <- lapply(1:D, function(d) z.est[(C_cum[d]+1):C_cum[d+1]])
save(ari.chdp,file=paste0('result/',rep_ind,'.ari.chdp.RData'))
save(z.chdp,file=paste0('result/',rep_ind,'.z.chdp.RData'))

# traceplot
pdf(file=paste0('fig_result/trace_full_chdp_',rep_ind,'.pdf'),w=12,h=12)
par(mfrow=c(2,2))
plot(mcmc(chain1$alpha_output[ind1]),auto.layout = FALSE,main='alpha')
plot(mcmc(chain1$alpha_0_output[ind1]), auto.layout = FALSE,main='alpha0')
dev.off()




set.seed(13+rep_ind)
chdp.run.post <- mvn_chdp(Y = Y, x = x, J = J.est, Z_fix = z.chdp, niter = niter, burn_in = burn_in, thinning = thinning,
                          partial.save.name = paste0('result/',rep_ind,'.post.chdp.RData'), save_frequency = save_frequency, 
                          auto.save = auto.save, verbose = verbose)
chain1 <- chdp.run.post
# traceplot
pdf(file=paste0('fig_result/trace_post_chdp_',rep_ind,'.pdf'),w=12,h=12)
par(mfrow=c(2,2))
plot(mcmc(chain1$alpha_output[ind1]),auto.layout = FALSE,main='alpha')
plot(mcmc(chain1$alpha_0_output[ind1]), auto.layout = FALSE,main='alpha0')
dev.off()

# covariate-dependent probability
pdf(file=paste0('fig_result/',rep_ind,'_px_chdp.pdf'),w=12,h=12)
par(mfrow=c(D,J.est),mar=c(5,4,2,2))
for(d in 1:D){
  for (j in 1:J.est) {
    plot(x[[d]],x[[d]],ylim=c(0,1),type='n',
         main=paste('Dataset',d,': Cluster',j, '( n = ',sum(z.chdp[[d]]==j),')'),xlab='x',ylab='prob')
    for (i in indd1) {
      lines(x[[d]],chain1$P_C_J_D_output[[d]][[i]][,j],col=alpha('grey',0.2))
    }
    for(jj in 1:J){
      lines(x[[d]],px[[d]][,jj],col=hue_pal()(J)[jj])
    }
  }
  
}
dev.off()   