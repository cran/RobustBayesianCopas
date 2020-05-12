###################################
## ROBUST BAYESIAN META-ANALYSIS ##
## WITH NO PUBLICATION BIAS      ##
###################################
# Implements the robust Bayesian Copas (RBC) selection model when there is
# no publication bias, i.e. rho=0.
# See Bai et al. (2020): https://arxiv.org/abs/2005.02930

# INPUTS:
# y = given vector of observed treatment effects
# s = given vector of within-study standard errors
# het.dist = distribution of heterogeneity. Allows either "Cauchy" or "normal."
#            The default is "Cauchy."
# init = optional initialization values of (theta, tau)
#        If the user does not provide, these are estimated from the data.
# burn = number of burn-in samples. Default is 10,000
# nmc = number of posterior draws to be saved. Default is 10,000

# OUTPUTS:
# theta.hat = point estimate (posterior mean) for theta
# theta.samples = MCMC samples for theta after burn-in
# tau.hat = point estimate (posterior mean) for tau
# tau.samples = MCMC samples for tau after burnin
BayesNonBiasCorrected <- function(y, s, init = NULL, het.dist=c("Cauchy","normal"),
                                  burn=10000, nmc=10000){
  
  # Check that y and s are of the same length
  if(length(y) != length(s))
    stop("Please ensure that 'y' and 's' are of the same length.")
    
  # Sample size
  n = length(y)
  
  # Initial values 
  if(is.null(init)){
  
    ### obtain initial estimators using the MLE for 
    SMA = StandardMetaAnalysis(y, s)
    
    # List for initial values of (theta, tau)
    init.vals = list(theta = SMA$theta.hat, 
                     inv.tausq = (1/SMA$tau.hat)^2)
  } else if(!is.null(init)){
    
    if(length(init)!=2)
      stop("ERROR: Please enter initial values for (theta, tau) in this exact order.")
    
    if(init[2]<=0)
      stop("ERROR: Please enter initial tau > 0.")
    
    # List for initial values of (theta, tau)
    init.vals = list(theta = init[1],
                     inv.tausq = (1/init[2])^2)
  }
  
  # Coercion
  het.dist <- match.arg(het.dist)
  
  # If normal distribution is specified for the heterogeneity
  if(het.dist == "normal"){
    #####################
    # Specify the model #
    #####################
    model_string <- "model{

    # Likelihood
    for(i in 1:n){
      y[i] ~ dnorm(mu[i], 1/pow(s[i],2))
    }
  
    # Prior for mu[i]'s
    for(i in 1:n){
      mu[i] ~ dnorm(theta,inv.tausq) # normal distribution on random effects
    }
  
    # Prior for theta
    theta ~ dnorm(0,0.0001)
  
    # Prior for tau2
    inv.tausq ~ dgamma(0.0001,0.0001)
    tau <- pow(1/inv.tausq,0.5)

    }"
    
    # Compile the model in JAGS
    model <- rjags::jags.model(textConnection(model_string), 
                               data = list(y=y,s=s,n=n),
                               inits = init.vals,
                               n.chains=1,
                               quiet=TRUE)
    
    # Draw samples
    # First use function 'update' to draw 10,000 warm-up (burnin) samples
    stats::update(model, n.iter=burn, progress.bar="none") # Burnin for 10,000 samples
    
    # Next use the function 'coda.samples' to produce the next 10,000
    # samples which will ultimately used to approximate the posterior.
    full.samples <- rjags::coda.samples(model, 
                                  variable.names=c("theta","tau"), 
                                  n.iter=nmc, progress.bar="none")
    
    # Extract samples. Can be used for plotting, obtaining summary statistics, 
    # and obtaining the D measure.
    full.samples.mat <- as.matrix(full.samples[[1]])
    theta.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="theta")]
    tau.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="tau")]
  
  } else { # Otherwise use the default Cauchy distribution for u's
    
    #####################
    # Specify the model #
    #####################
    model_string <- "model{

    # Likelihood
    for(i in 1:n){
      y[i] ~ dnorm(mu[i], 1/pow(s[i],2))
    }
  
    # Prior for mu[i]'s
    for(i in 1:n){
      mu[i] ~ dt(theta,inv.tausq,1) # Cauchy is a t-distribution with 1 df
    }
  
    # Prior for theta
    theta ~ dnorm(0,0.0001)
  
    # Prior for tau2
    inv.tausq ~ dgamma(0.0001,0.0001)
    tau <- pow(1/inv.tausq,0.5)

    }"
  
    # Compile the model in JAGS
    model <- rjags::jags.model(textConnection(model_string), 
                               data = list(y=y,s=s,n=n),
                               inits = init.vals,
                               n.chains=1,
                               quiet=TRUE)
  
    # Draw samples
    # First use function 'update' to draw 10,000 warm-up (burnin) samples
    stats::update(model, n.iter=burn, progress.bar="none") # Burnin for 10,000 samples
  
    # Next use the function 'coda.samples' to produce the next 10,000
    # samples which will ultimately used to approximate the posterior.
    full.samples <- rjags::coda.samples(model, 
                                 variable.names=c("theta","tau"), 
                                 n.iter=nmc, progress.bar="none")
  
    # Extract samples. Can be used for plotting, obtaining summary statistics, 
    # and obtaining the D measure.
    full.samples.mat <- as.matrix(full.samples[[1]])
    theta.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="theta")]
    tau.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="tau")]
  
  }  
  
  # Extract point estimates
  theta.hat <- mean(theta.samples)
  tau.hat <- mean(tau.samples)
  
  # Return list   
  return(list(theta.hat = theta.hat, 
              theta.samples = theta.samples,
              tau.hat = tau.hat,
              tau.samples = tau.samples))
}  