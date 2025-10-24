# ----- hdp ---------
source("algorithm/hdp_mcmc.R")
source("algorithm/hdp_draw_sample.R")

set.seed(4+rep_ind)
hdp.run <- mvn_hdp(Y = Y, J = JJ, niter = niter, burn_in = burn_in, thinning = thinning, empirical_z = empirical_z,
                   partial.save.name = paste0('result/',rep_ind,'.full.hdp.RData'), save_frequency = save_frequency, 
                   auto.save = auto.save, verbose = verbose)

# get optimal clustering
chain1 <- hdp.run
z.est <- get_optimal_cl(chain1, ind1, C)
J.est <- length(unique(z.est))
## clean up
ari.hdp <- adjustedRandIndex(z.est,unlist(z))
z.hdp <- lapply(1:D, function(d) z.est[(C_cum[d]+1):C_cum[d+1]])
save(ari.hdp,file=paste0('result/',rep_ind,'.ari.hdp.RData'))
save(z.hdp,file=paste0('result/',rep_ind,'.z.hdp.RData'))

# traceplot
pdf(file=paste0('fig_result/trace_full_hdp_',rep_ind,'.pdf'),w=12,h=12)
par(mfrow=c(2,2))
plot(mcmc(chain1$alpha_output[ind1]),auto.layout = FALSE,main='alpha')
plot(mcmc(chain1$alpha_0_output[ind1]), auto.layout = FALSE,main='alpha0')
dev.off()

set.seed(13+rep_ind)
hdp.run.post <- mvn_hdp(Y = Y, J = J.est, Z_fix = z.hdp, niter = niter, burn_in = burn_in, thinning = thinning,
                        partial.save.name = paste0('result/',rep_ind,'.post.hdp.RData'), save_frequency = save_frequency, 
                        auto.save = auto.save, verbose = verbose)

chain1 <- hdp.run.post
# traceplot
pdf(file=paste0('fig_result/trace_post_hdp_',rep_ind,'.pdf'),w=12,h=12)
par(mfrow=c(2,2))
plot(mcmc(chain1$alpha_output[ind1]),auto.layout = FALSE,main='alpha')
plot(mcmc(chain1$alpha_0_output[ind1]), auto.layout = FALSE,main='alpha0')
dev.off()

# probability
pdf(file=paste0('fig_result/',rep_ind,'_px_hdp.pdf'),w=12,h=12)
par(mfrow=c(D,J.est),mar=c(5,4,2,2))
for(d in 1:D){
  for (j in 1:J.est) {
    plot(x[[d]],x[[d]],ylim=c(0,1),type='n',
         main=paste('Dataset',d,': Cluster',j, '( n = ',sum(z.hdp[[d]]==j),')'),xlab='x',ylab='prob')
    for (i in indd1) {
      lines(x[[d]],rep(chain1$P_J_D_output[[i]][j,d],C[d]),col=alpha('grey',0.2))
    }
    for(jj in 1:J){
      lines(x[[d]],px[[d]][,jj],col=hue_pal()(J)[jj])
    }
  }
  
}
dev.off()




