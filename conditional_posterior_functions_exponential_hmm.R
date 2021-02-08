# Per gene conditional log posterior mean function
# Uses vectors of per-gene sample means and mean and dispersion estimates
mean.conditional.log.posterior <- function(n, sample.means, means, disps, prior.rate) {
  x <- 1/disps
  out <- rep(-1e20, length(means))
  index <- which(means>0)
  out[index] <- -n*(sample.means[index] + x[index]) * log(1 + means[index]*disps[index]) + 
    n*sample.means[index]*log(means[index]) - 
    prior.rate*means[index]
  return(out)
}


# Per gene conditional log posterior dispersion function
# Uses vector of counts for each gene to calculate sum over log gamma for each gene, 
# and vectors of per-gene sample means and mean and dispersion estimates
disp.conditional.log.posterior <- function(genes, counts, n, sample.means, 
                                           means, disps, prior.rate) {
  x <- 1/disps
  out <- rep(-1e20, genes)
  index <- which(disps>0)
  lgammasum <- colSums(lgamma(counts + rep(x, each = nrow(counts))))
  out[index] <- lgammasum[index] - 
    n*lgamma(x[index]) - 
    n*(sample.means[index] + x[index]) * log(1 + means[index]*disps[index]) + 
    n*sample.means[index]*log(disps[index]) - 
    prior.rate*disps[index]
  return(out)
}


# Function for rate of gamma conditional posteriors for exponential rate parameters
prior.rate.posterior.params <- function(g, z, hyperprior.rate, 
                                        parameters0, parameters1, parameters2) {
  return(c(1 + g + sum(z), 
           hyperprior.rate + sum((1-z)*parameters0 + z*(parameters1 + parameters2))))
}


# Function for log posterior probability z=0
# Uses vector of counts for each gene to calculate sum over log gamma for each gene, 
# and vectors of per-gene sample means and mean and dispersion estimates
pz0 <- function(genes, counts, n, sample.means, means, disps, 
                mean.prior.rate, disp.prior.rate, lambda) {
  x <- 1/disps
  lgammasum <- colSums(lgamma(counts + rep(x, each = nrow(counts))))
  return(lgammasum - 
           n*lgamma(x) - 
           n*(sample.means + x) * log(1 + means*disps) + 
           n*sample.means * log(means*disps) - 
           mean.prior.rate*means - 
           disp.prior.rate*disps + 
           log(1-lambda))
}


# Function for log posterior probability z=1
# Uses vector of counts for each gene to calculate sums over log gamma for each gene 
# for each group, 
# and vectors of per-gene sample means and mean and dispersion estimates for each group
pz1 = function(genes, counts1, counts2, n1, n2, 
               sample.means1, sample.means2, means1, means2, disps1, disps2, 
               mean.prior.rate, disp.prior.rate, lambda) {
  x1 <- 1/disps1
  x2 <- 1/disps2
  lgammasum1 <- colSums(lgamma(counts1 + rep(x1, each = nrow(counts1))))
  lgammasum2 <- colSums(lgamma(counts2 + rep(x2, each = nrow(counts2))))
  return(lgammasum1 - 
           n1*lgamma(x1) - 
           n1*(sample.means1 + x1) * log(1 + means1*disps1) + 
           n1*sample.means1 * log(means1*disps1) +
           lgammasum2 - 
           n2*lgamma(x2) - 
           n2*(sample.means2 + x2) * log(1 + means2*disps2) + 
           n2*sample.means2 * log(means2*disps2) + 
           log(mean.prior.rate) + 
           log(disp.prior.rate) - 
           mean.prior.rate * (means1 + means2) - 
           disp.prior.rate * (disps1 + disps2) + 
           log(lambda))
}


# Function to calculate probabilities for posterior Bernoulli distributions for 
# mixture components (exponentialtes and normalises calculated probabilities)
posterior.indicator.probabilities <- function(genes, counts0, counts1, counts2, n, n1, n2, 
                                              sample.means0, sample.means1, sample.means2, 
                                              means0, means1, means2, disps0, disps1, disps2, 
                                              mean.prior.rate, disp.prior.rate, lambda) {
  return(1 / (1 + exp(pz0(genes=genes, 
                          counts=counts0, 
                          n=n, 
                          sample.means=sample.means0, 
                          means=means0, 
                          disps=disps0, 
                          mean.prior.rate=mean.prior.rate, 
                          disp.prior.rate=disp.prior.rate, 
                          lambda=lambda) - 
                        pz1(genes=genes, 
                            counts1=counts1, 
                            counts2=counts2, 
                            n1=n1, 
                            n2=n2, 
                            sample.means1=sample.means1, 
                            sample.means2=sample.means2, 
                            means1=means1, 
                            means2=means2, 
                            disps1=disps1, 
                            disps2=disps2, 
                            mean.prior.rate=mean.prior.rate, 
                            disp.prior.rate=disp.prior.rate, 
                            lambda=lambda)
  )
  )
  )
}
