# Need to load log conditional posterior functions
# Inputs: 
# chain.length - number of iterations
# inits - list containing starting values for parameters estimated by MH:
# inits$mean0, mean1, mean2, disps0, disps1, disps2; 
# as well as initial values for parameters estimated by GS for input to first iteration 
# of MH for other parameters: mean.prior.rate, disp.prior.rate
# counts matrix of counts (samples in rows, genes in columns)
# groups - vector of 1s and 2s indicating which sample is in which group


exp_hmm_1_chain <- function(counts, groups, chain.length, thin=1, inits, 
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
  posterior.means0 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.means1 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.means2 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.disps0 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.disps1 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.disps2 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.mean.prior.rate <- numeric(chain.length/thin)
  posterior.disp.prior.rate <- numeric(chain.length/thin)
  posterior.indicators <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.proportion <- numeric(chain.length/thin)
  
  # Initial values ####
  current.posterior.means0 <- inits$means0
  current.posterior.means1 <- inits$means1
  current.posterior.means2 <- inits$means2
  current.posterior.disps0 <- inits$disps0
  current.posterior.disps1 <- inits$disps1
  current.posterior.disps2 <- inits$disps2
  current.posterior.mean.prior.rate <- inits$mean.prior.rate
  current.posterior.disp.prior.rate <- inits$mean.prior.rate
  current.posterior.indicators <- numeric(genes)
  current.posterior.proportion = 0.5
  
  # Create acceptance rate vectors/variables ####
  accept.means0 <- numeric(genes)
  accept.means1 <- numeric(genes)
  accept.means2 <- numeric(genes)
  accept.disps0 <- numeric(genes)
  accept.disps1 <- numeric(genes)
  accept.disps2 <- numeric(genes)
  
  # Run MCMC ####
  for (iter in 1:chain.length) {
    # Metropolis updates for per-gene overall means ####
    proposed.posterior.means0 <- rlnorm(n=genes, 
                                        meanlog=log(current.posterior.means0), 
                                        sdlog=sqrt.mean.proposal.scales0)
    replace <- log(runif(genes)) <= 
      mean.conditional.log.posterior(n=samples0, 
                                     sample.means=sample.means0,
                                     means=proposed.posterior.means0, 
                                     disps=current.posterior.disps0, 
                                     prior.rate=current.posterior.mean.prior.rate) + 
      (log(proposed.posterior.means0) + (log(proposed.posterior.means0) - 
                                           log(current.posterior.means0))^2 / (2*mean.proposal.scales0)) - 
      mean.conditional.log.posterior(n=samples0, 
                                     sample.means=sample.means0, 
                                     means=current.posterior.means0, 
                                     disps=current.posterior.disps0, 
                                     prior.rate=current.posterior.mean.prior.rate) - 
      (log(current.posterior.means0) + (log(current.posterior.means0) - 
                                          log(proposed.posterior.means0))^2 / (2*mean.proposal.scales0))
    current.posterior.means0[replace] <- proposed.posterior.means0[replace]
    accept.means0[replace] <- accept.means0[replace] + 1/chain.length
    
    # Metropolis updates for per-gene group 1 means ####
    proposed.posterior.means1 <- rlnorm(n=genes, 
                                        meanlog=log(current.posterior.means1), 
                                        sdlog=sqrt.mean.proposal.scales1)
    replace <- log(runif(genes)) <= 
      mean.conditional.log.posterior(n=samples1, 
                                     sample.means=sample.means1,
                                     means=proposed.posterior.means1, 
                                     disps=current.posterior.disps1, 
                                     prior.rate=current.posterior.mean.prior.rate) + 
      (log(proposed.posterior.means1) + (log(proposed.posterior.means1) - 
                                           log(current.posterior.means1))^2 / (2*mean.proposal.scales1)) - 
      mean.conditional.log.posterior(n=samples1, 
                                     sample.means=sample.means1, 
                                     means=current.posterior.means1, 
                                     disps=current.posterior.disps1, 
                                     prior.rate=current.posterior.mean.prior.rate) - 
      (log(current.posterior.means1) + (log(current.posterior.means1) - 
                                          log(proposed.posterior.means1))^2 / (2*mean.proposal.scales1))
    current.posterior.means1[replace] <- proposed.posterior.means1[replace]
    accept.means1[replace] <- accept.means1[replace] + 1/chain.length
    
    # Metropolis updates for per-gene group 2 means ####
    proposed.posterior.means2 <- rlnorm(n=genes, 
                                        meanlog=log(current.posterior.means2), 
                                        sdlog=sqrt.mean.proposal.scales2)
    replace <- log(runif(genes)) <= 
      mean.conditional.log.posterior(n=samples2, 
                                     sample.means=sample.means2,
                                     means=proposed.posterior.means2, 
                                     disps=current.posterior.disps2, 
                                     prior.rate=current.posterior.mean.prior.rate) + 
      (log(proposed.posterior.means2) + (log(proposed.posterior.means2) - 
                                           log(current.posterior.means2))^2 / (2*mean.proposal.scales2)) - 
      mean.conditional.log.posterior(n=samples2, 
                                     sample.means=sample.means2, 
                                     means=current.posterior.means2, 
                                     disps=current.posterior.disps2, 
                                     prior.rate=current.posterior.mean.prior.rate) - 
      (log(current.posterior.means2) + (log(current.posterior.means2) - 
                                          log(proposed.posterior.means2))^2 / (2*mean.proposal.scales2))
    current.posterior.means2[replace] <- proposed.posterior.means2[replace]
    accept.means2[replace] <- accept.means2[replace] + 1/chain.length
    
    # Metropolis updates for per-gene overall dispersions ####
    proposed.posterior.disps0 <- rlnorm(n=genes, 
                                        meanlog=log(current.posterior.disps0), 
                                        sdlog=sqrt.disp.proposal.scales0)
    replace <- log(runif(genes)) <= 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts, 
                                     n=samples0, 
                                     sample.means=sample.means0, 
                                     means=current.posterior.means0, 
                                     disps=proposed.posterior.disps0, 
                                     prior.rate=current.posterior.disp.prior.rate) + 
      (log(proposed.posterior.disps0) + (log(proposed.posterior.disps0) - 
                                           log(current.posterior.disps0))^2 / (2*disp.proposal.scales0)) - 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts, 
                                     n=samples0, 
                                     sample.means=sample.means0, 
                                     means=current.posterior.means0, 
                                     disps=current.posterior.disps0, 
                                     prior.rate=current.posterior.disp.prior.rate) - 
      (log(current.posterior.disps0) + (log(current.posterior.disps0) - 
                                          log(proposed.posterior.disps0))^2 / (2*disp.proposal.scales0))
    current.posterior.disps0[replace] <- proposed.posterior.disps0[replace]
    accept.disps0[replace] <- accept.disps0[replace] + 1/chain.length
    
    # Metropolis updates for per-gene group 1 dispersions ####
    proposed.posterior.disps1 <- rlnorm(n=genes, 
                                        meanlog=log(current.posterior.disps1), 
                                        sdlog=sqrt.disp.proposal.scales1)
    replace <- log(runif(genes)) <= 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts1, 
                                     n=samples1, 
                                     sample.means=sample.means1, 
                                     means=current.posterior.means1, 
                                     disps=proposed.posterior.disps1, 
                                     prior.rate=current.posterior.disp.prior.rate) + 
      (log(proposed.posterior.disps1) + (log(proposed.posterior.disps1) - 
                                           log(current.posterior.disps1))^2 / (2*disp.proposal.scales1)) - 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts1, 
                                     n=samples1, 
                                     sample.means=sample.means1, 
                                     means=current.posterior.means1, 
                                     disps=current.posterior.disps1, 
                                     prior.rate=current.posterior.disp.prior.rate) - 
      (log(current.posterior.disps1) + (log(current.posterior.disps1) - 
                                          log(proposed.posterior.disps1))^2 / (2*disp.proposal.scales1))
    current.posterior.disps1[replace] <- proposed.posterior.disps1[replace]
    accept.disps1[replace] <- accept.disps1[replace] + 1/chain.length
    
    # Metropolis updates for per-gene group 2 dispersions ####
    proposed.posterior.disps2 <- rlnorm(n=genes, 
                                        meanlog=log(current.posterior.disps2), 
                                        sdlog=sqrt.disp.proposal.scales2)
    replace <- log(runif(genes)) <= 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts2, 
                                     n=samples2, 
                                     sample.means=sample.means2, 
                                     means=current.posterior.means2, 
                                     disps=proposed.posterior.disps2, 
                                     prior.rate=current.posterior.disp.prior.rate) + 
      (log(proposed.posterior.disps2) + (log(proposed.posterior.disps2) - 
                                           log(current.posterior.disps2))^2 / (2*disp.proposal.scales2)) - 
      disp.conditional.log.posterior(genes=genes, 
                                     counts=counts2, 
                                     n=samples2, 
                                     sample.means=sample.means2, 
                                     means=current.posterior.means2, 
                                     disps=current.posterior.disps2, 
                                     prior.rate=current.posterior.disp.prior.rate) - 
      (log(current.posterior.disps2) + (log(current.posterior.disps2) - 
                                          log(proposed.posterior.disps2))^2 / (2*disp.proposal.scales2))
    current.posterior.disps2[replace] <- proposed.posterior.disps2[replace]
    accept.disps2[replace] <- accept.disps2[replace] + 1/chain.length
    
    # Gibbs update for prior rate parameter for mean ####
    posterior.mean.prior.rate.params <- prior.rate.posterior.params(g=genes, 
                                                                    z=current.posterior.indicators, 
                                                                    hyperprior.rate=2, 
                                                                    parameters0=current.posterior.means0, 
                                                                    parameters1=current.posterior.means1, 
                                                                    parameters2=current.posterior.means2)
    current.posterior.mean.prior.rate <- 
      rgamma(n=1, 
             shape=posterior.mean.prior.rate.params[1], 
             rate=posterior.mean.prior.rate.params[2])
    
    # Gibbs update for prior rate parameter for dispersion ####
    posterior.disp.prior.rate.params <- prior.rate.posterior.params(g=genes, 
                                                                    z=current.posterior.indicators, 
                                                                    hyperprior.rate=0.2, 
                                                                    parameters0=current.posterior.disps0, 
                                                                    parameters1=current.posterior.disps1, 
                                                                    parameters2=current.posterior.disps2)
    current.posterior.disp.prior.rate <- 
      rgamma(n=1, 
             shape=posterior.disp.prior.rate.params[1], 
             rate=posterior.disp.prior.rate.params[2])
    
    # Gibbs updates for per-gene mixture components ####
    current.posterior.indicators <- 
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
                                                    means0=current.posterior.means0, 
                                                    means1=current.posterior.means1, 
                                                    means2=current.posterior.means2, 
                                                    disps0=current.posterior.disps0, 
                                                    disps1=current.posterior.disps1, 
                                                    disps2=current.posterior.disps2, 
                                                    mean.prior.rate=current.posterior.mean.prior.rate, 
                                                    disp.prior.rate=current.posterior.disp.prior.rate, 
                                                    lambda=current.posterior.proportion))
    
    # Gibbs update for mixture proportion ####
    current.posterior.proportion <- rbeta(n=1, 
                                          shape1=1+sum(current.posterior.indicators), 
                                          shape2=1+genes-sum(current.posterior.indicators))
    
    # Update posterior samples ####
    if (iter/thin==round(iter/thin)) {
      posterior.means0[iter/thin,] <- current.posterior.means0
      posterior.means1[iter/thin,] <- current.posterior.means1
      posterior.means2[iter/thin,] <- current.posterior.means2
      posterior.disps0[iter/thin,] <- current.posterior.disps0
      posterior.disps1[iter/thin,] <- current.posterior.disps1
      posterior.disps2[iter/thin,] <- current.posterior.disps2
      posterior.mean.prior.rate[iter/thin] <- current.posterior.mean.prior.rate
      posterior.disp.prior.rate[iter/thin] <- current.posterior.disp.prior.rate
      posterior.indicators[iter/thin,] <- current.posterior.indicators
      posterior.proportion[iter/thin] <- current.posterior.proportion
    }
    
  }
  
  return(list("chain.length"=chain.length, 
              "thin"=thin, 
              "inits"=inits, 
              "mean.proposal.scales0"=mean.proposal.scales0, 
              "mean.proposal.scales1"=mean.proposal.scales1, 
              "mean.proposal.scales2"=mean.proposal.scales2, 
              "disp.proposal.scales0"=disp.proposal.scales0, 
              "disp.proposal.scales1"=disp.proposal.scales1, 
              "disp.proposal.scales2"=disp.proposal.scales2, 
              "accept.means0"=accept.means0, 
              "accept.means1"=accept.means1, 
              "accept.means2"=accept.means2, 
              "accept.disps0"=accept.disps0, 
              "accept.disps1"=accept.disps1, 
              "accept.disps2"=accept.disps2, 
              "posterior.means0"=posterior.means0, 
              "posterior.means1"=posterior.means1, 
              "posterior.means2"=posterior.means2, 
              "posterior.disps0"=posterior.disps0, 
              "posterior.disps1"=posterior.disps1, 
              "posterior.disps2"=posterior.disps2, 
              "posterior.mean.prior.rate"=posterior.mean.prior.rate, 
              "posterior.disp.prior.rate"=posterior.disp.prior.rate, 
              "posterior.indicators"=posterior.indicators, 
              "posterior.proportion"=posterior.proportion))
  
}
