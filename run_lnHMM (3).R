# This is the top-level function that runs the full algorithm and processes the results.

lnHMM <- function(counts, 
                  groups, 
                  chain.length=2000, 
                  initial.chain.length=100, 
                  adapt.chain.length=100, 
                  seed=1, 
                  return.raw.results=F) {
  
  # Set seed
  set.seed(seed)
  
  # Run MCMC algorithm
  raw.results <- ln_hmm_adapt_3_chains(counts=counts, 
                                       groups=groups, 
                                       chain.length=chain.length, 
                                       initial.chain.length=initial.chain.length, 
                                       adapt.chain.length=adapt.chain.length)
  
  if (return.raw.results) {
    return(raw.results)
  } else {
    
    # Extract and process results
    posterior.probabilities <- colMeans(raw.results$indicators)
    posterior.proportion <- mean(raw.results$proportion)
    threshold <- sort(posterior.probabilities, decreasing=T)[round(ncol(counts) * posterior.proportion)]
    bfdr <- bfdr(posterior.probabilities)
    mean.difference <- raw.results$means1 - raw.results$means2
    p.mean <- apply(mean.difference,2,hpd.pval)
    rm(mean.difference)
    mean.difference.log <- log(raw.results$means1) - log(raw.results$means2)
    p.mean.log <- apply(mean.difference.log,2,hpd.pval)
    rm(mean.difference.log)
    q.mean <- p.adjust(p.mean, method='BH')
    q.mean.log <- p.adjust(p.mean.log, method='BH')
    disp.difference <- raw.results$disps1 - raw.results$disps2
    p.disp <- apply(disp.difference,2,hpd.pval)
    rm(disp.difference)
    disp.difference.log <- log(raw.results$disps1) - log(raw.results$disps2)
    p.disp.log <- apply(disp.difference.log,2,hpd.pval)
    rm(disp.difference.log)
    q.disp <- p.adjust(p.disp, method='BH')
    q.disp.log <- p.adjust(p.disp.log, method='BH')
    
    # Return list with final results
    return(list(posterior.probabilities = posterior.probabilities,
                posterior.proportion = posterior.proportion,
                threshold = threshold,
                bfdr = bfdr,
                p.mean = p.mean,
                p.mean.log = p.mean.log,
                p.disp = p.disp,
                p.disp.log = p.disp.log,
                q.mean = q.mean,
                q.mean.log = q.mean.log,
                q.disp = q.disp,
                q.disp.log = q.disp.log))
  }
}
