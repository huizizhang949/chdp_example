# ------ ddp -------
source("algorithm/ddp_cat_mcmc.R")
source("algorithm/ddp_cat_draw_sample.R")

set.seed(4+rep_ind)
ddp.cat.run <- mvn_ddp1(Y = do.call(rbind, Y), x = unlist(x), xc = xc, J = JJ, 
                    niter = niter, burn_in = burn_in, thinning = thinning, empirical_z = empirical_z,
                    gamma = rep(1,D), trunc_dirichlet_burn = 50,
                    partial.save.name = paste0('result/',rep_ind,'.full.ddp.cat.RData'), save_frequency = save_frequency, 
                    auto.save = auto.save, verbose = verbose)
# get optimal clustering
chain1 <- ddp.cat.run
z.est <- get_optimal_cl(chain1, ind1, C)
J.est <- length(unique(z.est))
## clean up
ari.ddp.cat <- adjustedRandIndex(z.est,unlist(z))
z.ddp.cat <- lapply(1:D, function(d) z.est[(C_cum[d]+1):C_cum[d+1]])
save(ari.ddp.cat,file=paste0('result/',rep_ind,'.ari.ddp.cat.RData'))
save(z.ddp.cat,file=paste0('result/',rep_ind,'.z.ddp.cat.RData'))

# traceplot
pdf(file=paste0('fig_result/trace_full_ddp_cat_',rep_ind,'.pdf'),w=10,h=5)
par(mfrow=c(1,2))
plot(mcmc(chain1$alpha_0_output[ind1]), auto.layout = FALSE,main='alpha0')
dev.off()


set.seed(13+rep_ind)
ddp.cat.run.post <- mvn_ddp1(Y = do.call(rbind, Y), x = unlist(x), xc = xc, J = J.est, Z_fix = unlist(z.ddp.cat), 
                         niter = niter, burn_in = burn_in, thinning = thinning,
                         gamma = rep(1,D), trunc_dirichlet_burn = 50,
                         partial.save.name = paste0('result/',rep_ind,'.post.ddp.cat.RData'), save_frequency = save_frequency, 
                         auto.save = auto.save, verbose = verbose)
chain1 <- ddp.cat.run.post
# traceplot
pdf(file=paste0('fig_result/trace_post_ddp_cat_',rep_ind,'.pdf'),w=10,h=5)
par(mfrow=c(1,2))
plot(mcmc(chain1$alpha_0_output[ind1]), auto.layout = FALSE,main='alpha0')
dev.off()

# covariate-dependent probability
pdf(file=paste0('fig_result/',rep_ind,'_px_ddp_cat.pdf'),w=12,h=12)
par(mfrow=c(D,J.est),mar=c(5,4,2,2))
for(d in 1:D){
  for (j in 1:J.est) {
    plot(x[[d]],x[[d]],ylim=c(0,1),type='n',
         main=paste('Dataset',d,': Cluster',j, '( n = ',sum(z.ddp.cat[[d]]==j),')'),xlab='x',ylab='prob')
    for (i in indd1) {
      lines(x[[d]],chain1$P_C_J_output[[i]][(C_cum[d]+1):C_cum[d+1],j],col=alpha('grey',0.2))
    }
    for(jj in 1:J){
      lines(x[[d]],px[[d]][,jj],col=hue_pal()(J)[jj])
    }
  }
  
}
dev.off()   

