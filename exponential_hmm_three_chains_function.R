# Need to load log conditional posterior functions
# Inputs: 
# chain.length - number of iterations
# inits - list containing starting values for parameters estimated by MH:
# inits$mean0, mean1, mean2, disps0, disps1, disps2; 
# as well as initial values for parameters estimated by GS for input to first iteration 
# of MH for other parameters: mean.prior.rate, disp.prior.rate
# counts matrix of counts (samples in rows, genes in columns)
# groups - vector of 1s and 2s indicating which sample is in which group


exp_hmm_3_chains <- function(counts, groups, chain.length, thin=1, inits1, inits2, inits3, 
                             mean.proposal.scales0=rep(0.2, ncol(counts)), 
                             mean.proposal.scales1=rep(0.2, ncol(counts)), 
                             mean.proposal.scales2=rep(0.2, ncol(counts)), 
                             disp.proposal.scales0=rep(0.5, ncol(counts)), 
                             disp.proposal.scales1=rep(0.5, ncol(counts)), 
                             disp.proposal.scales2=rep(0.5, ncol(counts))) {
  
  genes <- ncol(counts)
  counts1 <- counts[groups==1,]
  counts2 <- counts[groups==2,]
  samples0 <- nrow(counts)
  samples1 <- nrow(counts1)
  samples2 <- nrow(counts2)
  sample.means0 <- colMeans(counts)
  sample.means1 <- colMeans(counts1)
  sample.means2 <- colMeans(counts2)
  sqrt.mean.proposal.scales0 <- sqrt(mean.proposal.scales0)
  sqrt.mean.proposal.scales1 <- sqrt(mean.proposal.scales1)
  sqrt.mean.proposal.scales2 <- sqrt(mean.proposal.scales2)
  sqrt.disp.proposal.scales0 <- sqrt(disp.proposal.scales0)
  sqrt.disp.proposal.scales1 <- sqrt(disp.proposal.scales1)
  sqrt.disp.proposal.scales2 <- sqrt(disp.proposal.scales2)
  
  # Create empty posterior sample matrices and vectors ####
  posterior.means0.1 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.means1.1 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.means2.1 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.means0.2 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.means1.2 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.means2.2 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.means0.3 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.means1.3 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.means2.3 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.disps0.1 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.disps1.1 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.disps2.1 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.disps0.2 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.disps1.2 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.disps2.2 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.disps0.3 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.disps1.3 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.disps2.3 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.mean.prior.rate.1 <- numeric(chain.length/thin)
  posterior.mean.prior.rate.2 <- numeric(chain.length/thin)
  posterior.mean.prior.rate.3 <- numeric(chain.length/thin)
  posterior.disp.prior.rate.1 <- numeric(chain.length/thin)
  posterior.disp.prior.rate.2 <- numeric(chain.length/thin)
  posterior.disp.prior.rate.3 <- numeric(chain.length/thin)
  posterior.indicators.1 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.indicators.2 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.indicators.3 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.proportion.1 <- numeric(chain.length/thin)
  posterior.proportion.2 <- numeric(chain.length/thin)
  posterior.proportion.3 <- numeric(chain.length/thin)
  
  # Initial values ####
  current.posterior.means0.1 <- inits1$means0
  current.posterior.means1.1 <- inits1$means1
  current.posterior.means2.1 <- inits1$means2
  current.posterior.means0.2 <- inits2$means0
  current.posterior.means1.2 <- inits2$means1
  current.posterior.means2.2 <- inits2$means2
  current.posterior.means0.3 <- inits3$means0
  current.posterior.means1.3 <- inits3$means1
  current.posterior.means2.3 <- inits3$means2
  current.posterior.disps0.1 <- inits1$disps0
  current.posterior.disps1.1 <- inits1$disps1
  current.posterior.disps2.1 <- inits1$disps2
  current.posterior.disps0.2 <- inits2$disps0
  current.posterior.disps1.2 <- inits2$disps1
  current.posterior.disps2.2 <- inits2$disps2
  current.posterior.disps0.3 <- inits3$disps0
  current.posterior.disps1.3 <- inits3$disps1
  current.posterior.disps2.3 <- inits3$disps2
  current.posterior.mean.prior.rate.1 <- inits1$mean.prior.rate
  current.posterior.mean.prior.rate.2 <- inits2$mean.prior.rate
  current.posterior.mean.prior.rate.3 <- inits3$mean.prior.rate
  current.posterior.disp.prior.rate.1 <- inits1$mean.prior.rate
  current.posterior.disp.prior.rate.2 <- inits2$mean.prior.rate
  current.posterior.disp.prior.rate.3 <- inits3$mean.prior.rate
  current.posterior.indicators.1 <- numeric(genes)
  current.posterior.indicators.2 <- numeric(genes)
  current.posterior.indicators.3 <- numeric(genes)
  current.posterior.proportion.1 = 0.5
  current.posterior.proportion.2 = 0.5
  current.posterior.proportion.3 = 0.5
  
  # Create acceptance rate vectors/variables ####
  accept.means0.1 <- numeric(genes)
  accept.means1.1 <- numeric(genes)
  accept.means2.1 <- numeric(genes)
  accept.means0.2 <- numeric(genes)
  accept.means1.2 <- numeric(genes)
  accept.means2.2 <- numeric(genes)
  accept.means0.3 <- numeric(genes)
  accept.means1.3 <- numeric(genes)
  accept.means2.3 <- numeric(genes)
  accept.disps0.1 <- numeric(genes)
  accept.disps1.1 <- numeric(genes)
  accept.disps2.1 <- numeric(genes)
  accept.disps0.2 <- numeric(genes)
  accept.disps1.2 <- numeric(genes)
  accept.disps2.2 <- numeric(genes)
  accept.disps0.3 <- numeric(genes)
  accept.disps1.3 <- numeric(genes)
  accept.disps2.3 <- numeric(genes)
  
  # Run MCMC ####
  for (iter in 1:chain.length) {
    # Metropolis updates for per-gene overall means ####
    # Chain 1
    proposed.posterior.means0.1 <- rlnorm(n=genes, 
                                          meanlog=log(current.posterior.means0.1), 
                                          sdlog=sqrt.mean.proposal.scales0)
    replace <- log(runif(genes)) <= 
      mean.conditional.log.posterior(n=samples0, 
                                     sample.means=sample.means0,
                                     means=proposed.posterior.means0.1, 
                                     disps=current.posterior.disps0.1, 
                                     prior.rate=current.posterior.mean.prior.rate.1) + 
      (log(proposed.posterior.means0.1) + (log(proposed.posterior.means0.1) - 
                                             log(current.posterior.means0.1))^2 / (2*mean.proposal.scales0)) - 
      mean.conditional.log.posterior(n=samples0, 
                                     sample.means=sample.means0, 
                                     means=current.posterior.means0.1, 
                                     disps=current.posterior.disps0.1, 
                                     prior.rate=current.posterior.mean.prior.rate.1) - 
      (log(current.posterior.means0.1) + (log(current.posterior.means0.1) - 
                                            log(proposed.posterior.means0.1))^2 / (2*mean.proposal.scales0))
    current.posterior.means0.1[replace] <- proposed.posterior.means0.1[replace]
    accept.means0.1[replace] <- accept.means0.1[replace] + 1/chain.length
    
    # Chain 2
    proposed.posterior.means0.2 <- rlnorm(n=genes, 
                                          meanlog=log(current.posterior.means0.2), 
                                          sdlog=sqrt.mean.proposal.scales0)
    replace <- log(runif(genes)) <= 
      mean.conditional.log.posterior(n=samples0, 
                                     sample.means=sample.means0,
                                     means=proposed.posterior.means0.2, 
                                     disps=current.posterior.disps0.2, 
                                     prior.rate=current.posterior.mean.prior.rate.2) + 
      (log(proposed.posterior.means0.2) + (log(proposed.posterior.means0.2) - 
                                             log(current.posterior.means0.2))^2 / (2*mean.proposal.scales0)) - 
      mean.conditional.log.posterior(n=samples0, 
                                     sample.means=sample.means0, 
                                     means=current.posterior.means0.2, 
                                     disps=current.posterior.disps0.2, 
                                     prior.rate=current.posterior.mean.prior.rate.2) - 
      (log(current.posterior.means0.2) + (log(current.posterior.means0.2) - 
                                            log(proposed.posterior.means0.2))^2 / (2*mean.proposal.scales0))
    current.posterior.means0.2[replace] <- proposed.posterior.means0.2[replace]
    accept.means0.2[replace] <- accept.means0.2[replace] + 1/chain.length
    
    # Chain 3
    proposed.posterior.means0.3 <- rlnorm(n=genes, 
                                          meanlog=log(current.posterior.means0.3), 
                                          sdlog=sqrt.mean.proposal.scales0)
    replace <- log(runif(genes)) <= 
      mean.conditional.log.posterior(n=samples0, 
                                     sample.means=sample.means0,
                                     means=proposed.posterior.means0.3, 
                                     disps=current.posterior.disps0.3, 
                                     prior.rate=current.posterior.mean.prior.rate.3) + 
      (log(proposed.posterior.means0.3) + (log(proposed.posterior.means0.3) - 
                                             log(current.posterior.means0.3))^2 / (2*mean.proposal.scales0)) - 
      mean.conditional.log.posterior(n=samples0, 
                                     sample.means=sample.means0, 
                                     means=current.posterior.means0.3, 
                                     disps=current.posterior.disps0.3, 
                                     prior.rate=current.posterior.mean.prior.rate.3) - 
      (log(current.posterior.means0.3) + (log(current.posterior.means0.3) - 
                                            log(proposed.posterior.means0.3))^2 / (2*mean.proposal.scales0))
    current.posterior.means0.3[replace] <- proposed.posterior.means0.3[replace]
    accept.means0.3[replace] <- accept.means0.3[replace] + 1/chain.length
    
    # Metropolis updates for per-gene group 1 means ####
    # Chain 1
    proposed.posterior.means1.1 <- rlnorm(n=genes, 
                                          meanlog=log(current.posterior.means1.1), 
                                          sdlog=sqrt.mean.proposal.scales1)
    replace <- log(runif(genes)) <= 
      mean.conditional.log.posterior(n=samples1, 
                                     sample.means=sample.means1,
                                     means=proposed.posterior.means1.1, 
                                     disps=current.posterior.disps1.1, 
                                     prior.rate=current.posterior.mean.prior.rate.1) + 
      (log(proposed.posterior.means1.1) + (log(proposed.posterior.means1.1) - 
                                             log(current.posterior.means1.1))^2 / (2*mean.proposal.scales1)) - 
      mean.conditional.log.posterior(n=samples1, 
                                     sample.means=sample.means1, 
                                     means=current.posterior.means1.1, 
                                     disps=current.posterior.disps1.1, 
                                     prior.rate=current.posterior.mean.prior.rate.1) - 
      (log(current.posterior.means1.1) + (log(current.posterior.means1.1) - 
                                            log(proposed.posterior.means1.1))^2 / (2*mean.proposal.scales1))
    current.posterior.means1.1[replace] <- proposed.posterior.means1.1[replace]
    accept.means1.1[replace] <- accept.means1.1[replace] + 1/chain.length
    
    # Chain 2
    proposed.posterior.means1.2 <- rlnorm(n=genes, 
                                          meanlog=log(current.posterior.means1.2), 
                                          sdlog=sqrt.mean.proposal.scales1)
    replace <- log(runif(genes)) <= 
      mean.conditional.log.posterior(n=samples1, 
                                     sample.means=sample.means1,
                                     means=proposed.posterior.means1.2, 
                                     disps=current.posterior.disps1.2, 
                                     prior.rate=current.posterior.mean.prior.rate.2) + 
      (log(proposed.posterior.means1.2) + (log(proposed.posterior.means1.2) - 
                                             log(current.posterior.means1.2))^2 / (2*mean.proposal.scales1)) - 
      mean.conditional.log.posterior(n=samples1, 
                                     sample.means=sample.means1, 
                                     means=current.posterior.means1.2, 
                                     disps=current.posterior.disps1.2, 
                                     prior.rate=current.posterior.mean.prior.rate.2) - 
      (log(current.posterior.means1.2) + (log(current.posterior.means1.2) - 
                                            log(proposed.posterior.means1.2))^2 / (2*mean.proposal.scales1))
    current.posterior.means1.2[replace] <- proposed.posterior.means1.2[replace]
    accept.means1.2[replace] <- accept.means1.2[replace] + 1/chain.length
    
    # Chain 3
    proposed.posterior.means1.3 <- rlnorm(n=genes, 
                                          meanlog=log(current.posterior.means1.3), 
                                          sdlog=sqrt.mean.proposal.scales1)
    replace <- log(runif(genes)) <= 
      mean.conditional.log.posterior(n=samples1, 
                                     sample.means=sample.means1,
                                     means=proposed.posterior.means1.3, 
                                     disps=current.posterior.disps1.3, 
                                     prior.rate=current.posterior.mean.prior.rate.3) + 
      (log(proposed.posterior.means1.3) + (log(proposed.posterior.means1.3) - 
                                             log(current.posterior.means1.3))^2 / (2*mean.proposal.scales1)) - 
      mean.conditional.log.posterior(n=samples1, 
                                     sample.means=sample.means1, 
                                     means=current.posterior.means1.3, 
                                     disps=current.posterior.disps1.3, 
                                     prior.rate=current.posterior.mean.prior.rate.3) - 
      (log(current.posterior.means1.3) + (log(current.posterior.means1.3) - 
                                            log(proposed.posterior.means1.3))^2 / (2*mean.proposal.scales1))
    current.posterior.means1.3[replace] <- proposed.posterior.means1.3[replace]
    accept.means1.3[replace] <- accept.means1.3[replace] + 1/chain.length
    
    # Metropolis updates for per-gene group 2 means ####
    # Chain 1
    proposed.posterior.means2.1 <- rlnorm(n=genes, 
                                          meanlog=log(current.posterior.means2.1), 
                                          sdlog=sqrt.mean.proposal.scales2)
    replace <- log(runif(genes)) <= 
      mean.conditional.log.posterior(n=samples2, 
                                     sample.means=sample.means2,
                                     means=proposed.posterior.means2.1, 
                                     disps=current.posterior.disps2.1, 
                                     prior.rate=current.posterior.mean.prior.rate.1) + 
      (log(proposed.posterior.means2.1) + (log(proposed.posterior.means2.1) - 
                                             log(current.posterior.means2.1))^2 / (2*mean.proposal.scales2)) - 
      mean.conditional.log.posterior(n=samples2, 
                                     sample.means=sample.means2, 
                                     means=current.posterior.means2.1, 
                                     disps=current.posterior.disps2.1, 
                                     prior.rate=current.posterior.mean.prior.rate.1) - 
      (log(current.posterior.means2.1) + (log(current.posterior.means2.1) - 
                                            log(proposed.posterior.means2.1))^2 / (2*mean.proposal.scales2))
    current.posterior.means2.1[replace] <- proposed.posterior.means2.1[replace]
    accept.means2.1[replace] <- accept.means2.1[replace] + 1/chain.length
    
    # Chain 2
    proposed.posterior.means2.2 <- rlnorm(n=genes, 
                                          meanlog=log(current.posterior.means2.2), 
                                          sdlog=sqrt.mean.proposal.scales2)
    replace <- log(runif(genes)) <= 
      mean.conditional.log.posterior(n=samples2, 
                                     sample.means=sample.means2,
                                     means=proposed.posterior.means2.2, 
                                     disps=current.posterior.disps2.2, 
                                     prior.rate=current.posterior.mean.prior.rate.2) + 
      (log(proposed.posterior.means2.2) + (log(proposed.posterior.means2.2) - 
                                             log(current.posterior.means2.2))^2 / (2*mean.proposal.scales2)) - 
      mean.conditional.log.posterior(n=samples2, 
                                     sample.means=sample.means2, 
                                     means=current.posterior.means2.2, 
                                     disps=current.posterior.disps2.2, 
                                     prior.rate=current.posterior.mean.prior.rate.2) - 
      (log(current.posterior.means2.2) + (log(current.posterior.means2.2) - 
                                            log(proposed.posterior.means2.2))^2 / (2*mean.proposal.scales2))
    current.posterior.means2.2[replace] <- proposed.posterior.means2.2[replace]
    accept.means2.2[replace] <- accept.means2.2[replace] + 1/chain.length
    
    # Chain 3
    proposed.posterior.means2.3 <- rlnorm(n=genes, 
                                          meanlog=log(current.posterior.means2.3), 
                                          sdlog=sqrt.mean.proposal.scales2)
    replace <- log(runif(genes)) <= 
      mean.conditional.log.posterior(n=samples2, 
                                     sample.means=sample.means2,
                                     means=proposed.posterior.means2.3, 
                                     disps=current.posterior.disps2.3, 
                                     prior.rate=current.posterior.mean.prior.rate.3) + 
      (log(proposed.posterior.means2.3) + (log(proposed.posterior.means2.3) - 
                                             log(current.posterior.means2.3))^2 / (2*mean.proposal.scales2)) - 
      mean.conditional.log.posterior(n=samples2, 
                                     sample.means=sample.means2, 
                                     means=current.posterior.means2.3, 
                                     disps=current.posterior.disps2.3, 
                                     prior.rate=current.posterior.mean.prior.rate.3) - 
      (log(current.posterior.means2.3) + (log(current.posterior.means2.3) - 
                                            log(proposed.posterior.means2.3))^2 / (2*mean.proposal.scales2))
    current.posterior.means2.3[replace] <- proposed.posterior.means2.3[replace]
    accept.means2.3[replace] <- accept.means2.3[replace] + 1/chain.length
    
    # Metropolis updates for per-gene overall dispersions ####
    # Chain 1
    proposed.posterior.disps0.1 <- rlnorm(n=genes, 
                                          meanlog=log(current.posterior.disps0.1), 
                                          sdlog=sqrt.disp.proposal.scales0)
    replace <- log(runif(genes)) <= 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts, 
                                     n=samples0, 
                                     sample.means=sample.means0, 
                                     means=current.posterior.means0.1, 
                                     disps=proposed.posterior.disps0.1, 
                                     prior.rate=current.posterior.disp.prior.rate.1) + 
      (log(proposed.posterior.disps0.1) + (log(proposed.posterior.disps0.1) - 
                                             log(current.posterior.disps0.1))^2 / (2*disp.proposal.scales0)) - 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts, 
                                     n=samples0, 
                                     sample.means=sample.means0, 
                                     means=current.posterior.means0.1, 
                                     disps=current.posterior.disps0.1, 
                                     prior.rate=current.posterior.disp.prior.rate.1) - 
      (log(current.posterior.disps0.1) + (log(current.posterior.disps0.1) - 
                                            log(proposed.posterior.disps0.1))^2 / (2*disp.proposal.scales0))
    current.posterior.disps0.1[replace] <- proposed.posterior.disps0.1[replace]
    accept.disps0.1[replace] <- accept.disps0.1[replace] + 1/chain.length
    
    # Chain 2
    proposed.posterior.disps0.2 <- rlnorm(n=genes, 
                                          meanlog=log(current.posterior.disps0.2), 
                                          sdlog=sqrt.disp.proposal.scales0)
    replace <- log(runif(genes)) <= 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts, 
                                     n=samples0, 
                                     sample.means=sample.means0, 
                                     means=current.posterior.means0.2, 
                                     disps=proposed.posterior.disps0.2, 
                                     prior.rate=current.posterior.disp.prior.rate.2) + 
      (log(proposed.posterior.disps0.2) + (log(proposed.posterior.disps0.2) - 
                                             log(current.posterior.disps0.2))^2 / (2*disp.proposal.scales0)) - 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts, 
                                     n=samples0, 
                                     sample.means=sample.means0, 
                                     means=current.posterior.means0.2, 
                                     disps=current.posterior.disps0.2, 
                                     prior.rate=current.posterior.disp.prior.rate.2) - 
      (log(current.posterior.disps0.2) + (log(current.posterior.disps0.2) - 
                                            log(proposed.posterior.disps0.2))^2 / (2*disp.proposal.scales0))
    current.posterior.disps0.2[replace] <- proposed.posterior.disps0.2[replace]
    accept.disps0.2[replace] <- accept.disps0.2[replace] + 1/chain.length
    
    # Chain 3
    proposed.posterior.disps0.3 <- rlnorm(n=genes, 
                                          meanlog=log(current.posterior.disps0.3), 
                                          sdlog=sqrt.disp.proposal.scales0)
    replace <- log(runif(genes)) <= 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts, 
                                     n=samples0, 
                                     sample.means=sample.means0, 
                                     means=current.posterior.means0.3, 
                                     disps=proposed.posterior.disps0.3, 
                                     prior.rate=current.posterior.disp.prior.rate.3) + 
      (log(proposed.posterior.disps0.3) + (log(proposed.posterior.disps0.3) - 
                                             log(current.posterior.disps0.3))^2 / (2*disp.proposal.scales0)) - 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts, 
                                     n=samples0, 
                                     sample.means=sample.means0, 
                                     means=current.posterior.means0.3, 
                                     disps=current.posterior.disps0.3, 
                                     prior.rate=current.posterior.disp.prior.rate.3) - 
      (log(current.posterior.disps0.3) + (log(current.posterior.disps0.3) - 
                                            log(proposed.posterior.disps0.3))^2 / (2*disp.proposal.scales0))
    current.posterior.disps0.3[replace] <- proposed.posterior.disps0.3[replace]
    accept.disps0.3[replace] <- accept.disps0.3[replace] + 1/chain.length
    
    # Metropolis updates for per-gene group 1 dispersions ####
    # Chain 1
    proposed.posterior.disps1.1 <- rlnorm(n=genes, 
                                          meanlog=log(current.posterior.disps1.1), 
                                          sdlog=sqrt.disp.proposal.scales1)
    replace <- log(runif(genes)) <= 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts1, 
                                     n=samples1, 
                                     sample.means=sample.means1, 
                                     means=current.posterior.means1.1, 
                                     disps=proposed.posterior.disps1.1, 
                                     prior.rate=current.posterior.disp.prior.rate.1) + 
      (log(proposed.posterior.disps1.1) + (log(proposed.posterior.disps1.1) - 
                                             log(current.posterior.disps1.1))^2 / (2*disp.proposal.scales1)) - 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts1, 
                                     n=samples1, 
                                     sample.means=sample.means1, 
                                     means=current.posterior.means1.1, 
                                     disps=current.posterior.disps1.1, 
                                     prior.rate=current.posterior.disp.prior.rate.1) - 
      (log(current.posterior.disps1.1) + (log(current.posterior.disps1.1) - 
                                            log(proposed.posterior.disps1.1))^2 / (2*disp.proposal.scales1))
    current.posterior.disps1.1[replace] <- proposed.posterior.disps1.1[replace]
    accept.disps1.1[replace] <- accept.disps1.1[replace] + 1/chain.length
    
    # Chain 2
    proposed.posterior.disps1.2 <- rlnorm(n=genes, 
                                          meanlog=log(current.posterior.disps1.2), 
                                          sdlog=sqrt.disp.proposal.scales1)
    replace <- log(runif(genes)) <= 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts1, 
                                     n=samples1, 
                                     sample.means=sample.means1, 
                                     means=current.posterior.means1.2, 
                                     disps=proposed.posterior.disps1.2, 
                                     prior.rate=current.posterior.disp.prior.rate.2) + 
      (log(proposed.posterior.disps1.2) + (log(proposed.posterior.disps1.2) - 
                                             log(current.posterior.disps1.2))^2 / (2*disp.proposal.scales1)) - 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts1, 
                                     n=samples1, 
                                     sample.means=sample.means1, 
                                     means=current.posterior.means1.2, 
                                     disps=current.posterior.disps1.2, 
                                     prior.rate=current.posterior.disp.prior.rate.2) - 
      (log(current.posterior.disps1.2) + (log(current.posterior.disps1.2) - 
                                            log(proposed.posterior.disps1.2))^2 / (2*disp.proposal.scales1))
    current.posterior.disps1.2[replace] <- proposed.posterior.disps1.2[replace]
    accept.disps1.2[replace] <- accept.disps1.2[replace] + 1/chain.length
    
    # Chain 3
    proposed.posterior.disps1.3 <- rlnorm(n=genes, 
                                          meanlog=log(current.posterior.disps1.3), 
                                          sdlog=sqrt.disp.proposal.scales1)
    replace <- log(runif(genes)) <= 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts1, 
                                     n=samples1, 
                                     sample.means=sample.means1, 
                                     means=current.posterior.means1.3, 
                                     disps=proposed.posterior.disps1.3, 
                                     prior.rate=current.posterior.disp.prior.rate.3) + 
      (log(proposed.posterior.disps1.3) + (log(proposed.posterior.disps1.3) - 
                                             log(current.posterior.disps1.3))^2 / (2*disp.proposal.scales1)) - 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts1, 
                                     n=samples1, 
                                     sample.means=sample.means1, 
                                     means=current.posterior.means1.3, 
                                     disps=current.posterior.disps1.3, 
                                     prior.rate=current.posterior.disp.prior.rate.3) - 
      (log(current.posterior.disps1.3) + (log(current.posterior.disps1.3) - 
                                            log(proposed.posterior.disps1.3))^2 / (2*disp.proposal.scales1))
    current.posterior.disps1.3[replace] <- proposed.posterior.disps1.3[replace]
    accept.disps1.3[replace] <- accept.disps1.3[replace] + 1/chain.length
    
    # Metropolis updates for per-gene group 2 dispersions ####
    # Chain 1
    proposed.posterior.disps2.1 <- rlnorm(n=genes, 
                                          meanlog=log(current.posterior.disps2.1), 
                                          sdlog=sqrt.disp.proposal.scales2)
    replace <- log(runif(genes)) <= 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts2, 
                                     n=samples2, 
                                     sample.means=sample.means2, 
                                     means=current.posterior.means2.1, 
                                     disps=proposed.posterior.disps2.1, 
                                     prior.rate=current.posterior.disp.prior.rate.1) + 
      (log(proposed.posterior.disps2.1) + (log(proposed.posterior.disps2.1) - 
                                             log(current.posterior.disps2.1))^2 / (2*disp.proposal.scales2)) - 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts2, 
                                     n=samples2, 
                                     sample.means=sample.means2, 
                                     means=current.posterior.means2.1, 
                                     disps=current.posterior.disps2.1, 
                                     prior.rate=current.posterior.disp.prior.rate.1) - 
      (log(current.posterior.disps2.1) + (log(current.posterior.disps2.1) - 
                                            log(proposed.posterior.disps2.1))^2 / (2*disp.proposal.scales2))
    current.posterior.disps2.1[replace] <- proposed.posterior.disps2.1[replace]
    accept.disps2.1[replace] <- accept.disps2.1[replace] + 1/chain.length
    
    # Chain 2
    proposed.posterior.disps2.2 <- rlnorm(n=genes, 
                                          meanlog=log(current.posterior.disps2.2), 
                                          sdlog=sqrt.disp.proposal.scales2)
    replace <- log(runif(genes)) <= 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts2, 
                                     n=samples2, 
                                     sample.means=sample.means2, 
                                     means=current.posterior.means2.2, 
                                     disps=proposed.posterior.disps2.2, 
                                     prior.rate=current.posterior.disp.prior.rate.2) + 
      (log(proposed.posterior.disps2.2) + (log(proposed.posterior.disps2.2) - 
                                             log(current.posterior.disps2.2))^2 / (2*disp.proposal.scales2)) - 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts2, 
                                     n=samples2, 
                                     sample.means=sample.means2, 
                                     means=current.posterior.means2.2, 
                                     disps=current.posterior.disps2.2, 
                                     prior.rate=current.posterior.disp.prior.rate.2) - 
      (log(current.posterior.disps2.2) + (log(current.posterior.disps2.2) - 
                                            log(proposed.posterior.disps2.2))^2 / (2*disp.proposal.scales2))
    current.posterior.disps2.2[replace] <- proposed.posterior.disps2.2[replace]
    accept.disps2.2[replace] <- accept.disps2.2[replace] + 1/chain.length
    
    # Chain 3
    proposed.posterior.disps2.3 <- rlnorm(n=genes, 
                                          meanlog=log(current.posterior.disps2.3), 
                                          sdlog=sqrt.disp.proposal.scales2)
    replace <- log(runif(genes)) <= 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts2, 
                                     n=samples2, 
                                     sample.means=sample.means2, 
                                     means=current.posterior.means2.3, 
                                     disps=proposed.posterior.disps2.3, 
                                     prior.rate=current.posterior.disp.prior.rate.3) + 
      (log(proposed.posterior.disps2.3) + (log(proposed.posterior.disps2.3) - 
                                             log(current.posterior.disps2.3))^2 / (2*disp.proposal.scales2)) - 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts2, 
                                     n=samples2, 
                                     sample.means=sample.means2, 
                                     means=current.posterior.means2.3, 
                                     disps=current.posterior.disps2.3, 
                                     prior.rate=current.posterior.disp.prior.rate.3) - 
      (log(current.posterior.disps2.3) + (log(current.posterior.disps2.3) - 
                                            log(proposed.posterior.disps2.3))^2 / (2*disp.proposal.scales2))
    current.posterior.disps2.3[replace] <- proposed.posterior.disps2.3[replace]
    accept.disps2.3[replace] <- accept.disps2.3[replace] + 1/chain.length
    
    # Gibbs update for prior rate parameter for mean ####
    posterior.mean.prior.rate.params <- prior.rate.posterior.params(g=genes, 
                                                                    z=current.posterior.indicators.1, 
                                                                    hyperprior.rate=2, 
                                                                    parameters0=current.posterior.means0.1, 
                                                                    parameters1=current.posterior.means1.1, 
                                                                    parameters2=current.posterior.means2.1)
    current.posterior.mean.prior.rate.1 <- 
      rgamma(n=1, 
             shape=posterior.mean.prior.rate.params[1], 
             rate=posterior.mean.prior.rate.params[2])
    
    posterior.mean.prior.rate.params <- prior.rate.posterior.params(g=genes, 
                                                                    z=current.posterior.indicators.2, 
                                                                    hyperprior.rate=2, 
                                                                    parameters0=current.posterior.means0.2, 
                                                                    parameters1=current.posterior.means1.2, 
                                                                    parameters2=current.posterior.means2.2)
    current.posterior.mean.prior.rate.2 <- 
      rgamma(n=1, 
             shape=posterior.mean.prior.rate.params[1], 
             rate=posterior.mean.prior.rate.params[2])
    
    posterior.mean.prior.rate.params <- prior.rate.posterior.params(g=genes, 
                                                                    z=current.posterior.indicators.3, 
                                                                    hyperprior.rate=2, 
                                                                    parameters0=current.posterior.means0.3, 
                                                                    parameters1=current.posterior.means1.3, 
                                                                    parameters2=current.posterior.means2.3)
    current.posterior.mean.prior.rate.3 <- 
      rgamma(n=1, 
             shape=posterior.mean.prior.rate.params[1], 
             rate=posterior.mean.prior.rate.params[2])
    
    # Gibbs update for prior rate parameter for dispersion ####
    posterior.disp.prior.rate.params <- prior.rate.posterior.params(g=genes, 
                                                                    z=current.posterior.indicators.1, 
                                                                    hyperprior.rate=0.2, 
                                                                    parameters0=current.posterior.disps0.1, 
                                                                    parameters1=current.posterior.disps1.1, 
                                                                    parameters2=current.posterior.disps2.1)
    current.posterior.disp.prior.rate.1 <- 
      rgamma(n=1, 
             shape=posterior.disp.prior.rate.params[1], 
             rate=posterior.disp.prior.rate.params[2])
    
    posterior.disp.prior.rate.params <- prior.rate.posterior.params(g=genes, 
                                                                    z=current.posterior.indicators.2, 
                                                                    hyperprior.rate=0.2, 
                                                                    parameters0=current.posterior.disps0.2, 
                                                                    parameters1=current.posterior.disps1.2, 
                                                                    parameters2=current.posterior.disps2.2)
    current.posterior.disp.prior.rate.2 <- 
      rgamma(n=1, 
             shape=posterior.disp.prior.rate.params[1], 
             rate=posterior.disp.prior.rate.params[2])
    
    posterior.disp.prior.rate.params <- prior.rate.posterior.params(g=genes, 
                                                                    z=current.posterior.indicators.3, 
                                                                    hyperprior.rate=0.2, 
                                                                    parameters0=current.posterior.disps0.3, 
                                                                    parameters1=current.posterior.disps1.3, 
                                                                    parameters2=current.posterior.disps2.3)
    current.posterior.disp.prior.rate.3 <- 
      rgamma(n=1, 
             shape=posterior.disp.prior.rate.params[1], 
             rate=posterior.disp.prior.rate.params[2])
    
    # Gibbs updates for per-gene mixture components ####
    current.posterior.indicators.1 <- 
      rbinom(n=genes, 
             size=1, 
             prob=posterior.indicator.probabilities(genes=genes, 
                                                    counts=counts, 
                                                    counts1=counts1, 
                                                    counts2=counts2, 
                                                    n=samples0, 
                                                    n1=samples1, 
                                                    n2=samples2, 
                                                    sample.means0=sample.means0, 
                                                    sample.means1=sample.means1, 
                                                    sample.means2=sample.means2,
                                                    means0=current.posterior.means0.1, 
                                                    means1=current.posterior.means1.1, 
                                                    means2=current.posterior.means2.1, 
                                                    disps0=current.posterior.disps0.1, 
                                                    disps1=current.posterior.disps1.1, 
                                                    disps2=current.posterior.disps2.1, 
                                                    mean.prior.rate=current.posterior.mean.prior.rate.1, 
                                                    disp.prior.rate=current.posterior.disp.prior.rate.1, 
                                                    lambda=current.posterior.proportion.1))
    
    current.posterior.indicators.2 <- 
      rbinom(n=genes, 
             size=1, 
             prob=posterior.indicator.probabilities(genes=genes, 
                                                    counts=counts, 
                                                    counts1=counts1, 
                                                    counts2=counts2, 
                                                    n=samples0, 
                                                    n1=samples1, 
                                                    n2=samples2, 
                                                    sample.means0=sample.means0, 
                                                    sample.means1=sample.means1, 
                                                    sample.means2=sample.means2,
                                                    means0=current.posterior.means0.2, 
                                                    means1=current.posterior.means1.2, 
                                                    means2=current.posterior.means2.2, 
                                                    disps0=current.posterior.disps0.2, 
                                                    disps1=current.posterior.disps1.2, 
                                                    disps2=current.posterior.disps2.2, 
                                                    mean.prior.rate=current.posterior.mean.prior.rate.2, 
                                                    disp.prior.rate=current.posterior.disp.prior.rate.2, 
                                                    lambda=current.posterior.proportion.2))
    
    current.posterior.indicators.3 <- 
      rbinom(n=genes, 
             size=1, 
             prob=posterior.indicator.probabilities(genes=genes, 
                                                    counts=counts, 
                                                    counts1=counts1, 
                                                    counts2=counts2, 
                                                    n=samples0, 
                                                    n1=samples1, 
                                                    n2=samples2, 
                                                    sample.means0=sample.means0, 
                                                    sample.means1=sample.means1, 
                                                    sample.means2=sample.means2,
                                                    means0=current.posterior.means0.3, 
                                                    means1=current.posterior.means1.3, 
                                                    means2=current.posterior.means2.3, 
                                                    disps0=current.posterior.disps0.3, 
                                                    disps1=current.posterior.disps1.3, 
                                                    disps2=current.posterior.disps2.3, 
                                                    mean.prior.rate=current.posterior.mean.prior.rate.3, 
                                                    disp.prior.rate=current.posterior.disp.prior.rate.3, 
                                                    lambda=current.posterior.proportion.3))
    
    # Gibbs update for mixture proportion ####
    current.posterior.proportion.1 <- rbeta(n=1, 
                                            shape1=1+sum(current.posterior.indicators.1), 
                                            shape2=1+genes-sum(current.posterior.indicators.1))
    
    current.posterior.proportion.2 <- rbeta(n=1, 
                                            shape1=1+sum(current.posterior.indicators.2), 
                                            shape2=1+genes-sum(current.posterior.indicators.2))
    
    current.posterior.proportion.3 <- rbeta(n=1, 
                                            shape1=1+sum(current.posterior.indicators.3), 
                                            shape2=1+genes-sum(current.posterior.indicators.3))
    
    # Update posterior samples ####
    if (iter/thin==round(iter/thin)) {
      posterior.means0.1[iter/thin,] <- current.posterior.means0.1
      posterior.means1.1[iter/thin,] <- current.posterior.means1.1
      posterior.means2.1[iter/thin,] <- current.posterior.means2.1
      posterior.means0.2[iter/thin,] <- current.posterior.means0.2
      posterior.means1.2[iter/thin,] <- current.posterior.means1.2
      posterior.means2.2[iter/thin,] <- current.posterior.means2.2
      posterior.means0.3[iter/thin,] <- current.posterior.means0.3
      posterior.means1.3[iter/thin,] <- current.posterior.means1.3
      posterior.means2.3[iter/thin,] <- current.posterior.means2.3
      posterior.disps0.1[iter/thin,] <- current.posterior.disps0.1
      posterior.disps1.1[iter/thin,] <- current.posterior.disps1.1
      posterior.disps2.1[iter/thin,] <- current.posterior.disps2.1
      posterior.disps0.2[iter/thin,] <- current.posterior.disps0.2
      posterior.disps1.2[iter/thin,] <- current.posterior.disps1.2
      posterior.disps2.2[iter/thin,] <- current.posterior.disps2.2
      posterior.disps0.3[iter/thin,] <- current.posterior.disps0.3
      posterior.disps1.3[iter/thin,] <- current.posterior.disps1.3
      posterior.disps2.3[iter/thin,] <- current.posterior.disps2.3
      posterior.mean.prior.rate.1[iter/thin] <- current.posterior.mean.prior.rate.1
      posterior.mean.prior.rate.2[iter/thin] <- current.posterior.mean.prior.rate.2
      posterior.mean.prior.rate.3[iter/thin] <- current.posterior.mean.prior.rate.3
      posterior.disp.prior.rate.1[iter/thin] <- current.posterior.disp.prior.rate.1
      posterior.disp.prior.rate.2[iter/thin] <- current.posterior.disp.prior.rate.2
      posterior.disp.prior.rate.3[iter/thin] <- current.posterior.disp.prior.rate.3
      posterior.indicators.1[iter/thin,] <- current.posterior.indicators.1
      posterior.indicators.2[iter/thin,] <- current.posterior.indicators.2
      posterior.indicators.3[iter/thin,] <- current.posterior.indicators.3
      posterior.proportion.1[iter/thin] <- current.posterior.proportion.1
      posterior.proportion.2[iter/thin] <- current.posterior.proportion.2
      posterior.proportion.3[iter/thin] <- current.posterior.proportion.3
    }
    
  }
  
  return(list("chain.length"=chain.length, 
              "thin"=thin, 
              "inits1"=inits1, 
              "inits2"=inits2, 
              "inits3"=inits3, 
              "mean.proposal.scales0"=mean.proposal.scales0, 
              "mean.proposal.scales1"=mean.proposal.scales1, 
              "mean.proposal.scales2"=mean.proposal.scales2, 
              "disp.proposal.scales0"=disp.proposal.scales0, 
              "disp.proposal.scales1"=disp.proposal.scales1, 
              "disp.proposal.scales2"=disp.proposal.scales2, 
              "accept.means0.1"=accept.means0.1, 
              "accept.means1.1"=accept.means1.1, 
              "accept.means2.1"=accept.means2.1, 
              "accept.means0.2"=accept.means0.2, 
              "accept.means1.2"=accept.means1.2, 
              "accept.means2.2"=accept.means2.2, 
              "accept.means0.3"=accept.means0.3, 
              "accept.means1.3"=accept.means1.3, 
              "accept.means2.3"=accept.means2.3, 
              "accept.disps0.1"=accept.disps0.1, 
              "accept.disps1.1"=accept.disps1.1, 
              "accept.disps2.1"=accept.disps2.1, 
              "accept.disps0.2"=accept.disps0.2, 
              "accept.disps1.2"=accept.disps1.2, 
              "accept.disps2.2"=accept.disps2.2, 
              "accept.disps0.3"=accept.disps0.3, 
              "accept.disps1.3"=accept.disps1.3, 
              "accept.disps2.3"=accept.disps2.3, 
              "posterior.means0.1"=posterior.means0.1, 
              "posterior.means1.1"=posterior.means1.1, 
              "posterior.means2.1"=posterior.means2.1, 
              "posterior.means0.2"=posterior.means0.2, 
              "posterior.means1.2"=posterior.means1.2, 
              "posterior.means2.2"=posterior.means2.2, 
              "posterior.means0.3"=posterior.means0.3, 
              "posterior.means1.3"=posterior.means1.3, 
              "posterior.means2.3"=posterior.means2.3, 
              "posterior.disps0.1"=posterior.disps0.1, 
              "posterior.disps1.1"=posterior.disps1.1, 
              "posterior.disps2.1"=posterior.disps2.1, 
              "posterior.disps0.2"=posterior.disps0.2, 
              "posterior.disps1.2"=posterior.disps1.2, 
              "posterior.disps2.2"=posterior.disps2.2, 
              "posterior.disps0.3"=posterior.disps0.3, 
              "posterior.disps1.3"=posterior.disps1.3, 
              "posterior.disps2.3"=posterior.disps2.3, 
              "posterior.mean.prior.rate.1"=posterior.mean.prior.rate.1, 
              "posterior.mean.prior.rate.2"=posterior.mean.prior.rate.2, 
              "posterior.mean.prior.rate.3"=posterior.mean.prior.rate.3, 
              "posterior.disp.prior.rate.1"=posterior.disp.prior.rate.1, 
              "posterior.disp.prior.rate.2"=posterior.disp.prior.rate.2, 
              "posterior.disp.prior.rate.3"=posterior.disp.prior.rate.3, 
              "posterior.indicators.1"=posterior.indicators.1, 
              "posterior.indicators.2"=posterior.indicators.2, 
              "posterior.indicators.3"=posterior.indicators.3, 
              "posterior.proportion.1"=posterior.proportion.1, 
              "posterior.proportion.2"=posterior.proportion.2, 
              "posterior.proportion.3"=posterior.proportion.3))
  
}
