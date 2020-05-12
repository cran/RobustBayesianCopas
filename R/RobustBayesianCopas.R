#################################################
## ROBUST BAYESIAN COPAS (RBC) SELECTION MODEL ##
#################################################
# Implements the full robust Bayesian Copas (RBC) selection model 
# of Bai et al. (2020): https://arxiv.org/abs/2005.02930

# INPUTS:
# y = given vector of observed treatment effects
# s = given vector of within-study standard errors 
# het.dist = distribution of heterogeneity. Allows either "Cauchy" or "normal."
#            The default is "Cauchy."
# init = optional initialization values of (theta, tau, rho, gamma0, gamma1)
#        If the user does not provide, these are estimated from the data.
# burn = number of burn-in samples. Default is 10,000
# nmc = number of posterior draws to be saved. Default is 10,000

# OUTPUTS:
# theta.hat = point estimate for theta. Returns posterior mean of theta
# theta.samples = MCMC samples for theta after burn-in. Used for plotting the 
#                 posterior for theta and performing inference for theta
# tau.hat = point estimate for tau. Returns posterior mean of tau
# tau.samples = MCMC samples for tau after burn-in. Used for plotting the 
#               posterior for tau and performing inference for tau 
# rho.hat = point estimate for rho. Returns posterior median of rho
# rho.samples = MCMC samples for rho after burn-in. Used for plotting the 
#               posterior for rho and performing inference for rho
# gamma0.hat = point estimate for gamma0. Returns posterior median of gamma0
# gamma0.samples = MCMC samples for gamma0 after burn-in. Used for plotting the 
#                  posterior for gamma0 and performing inference for gamma0 
# gamma1.hat = point estimate for gamma1. Returns posterior median of gamma1
# gamma1.samples = MCMC samples for gamma1 after burn-in. Used for plotting the
#                  posterior for gamma1 and performing inference for gamma1 

RobustBayesianCopas <- function(y, s, init=NULL, het.dist=c("Cauchy","normal"),
                                burn=10000, nmc=10000){
  
  # Check that y and s are of the same length
  if(length(y) != length(s))
    stop("ERROR: Please ensure that 'y' and 's' are of the same length.")

  # Sample size
  n = length(y)
  data = data.frame(cbind(y,s))
  max.s = max(s)
  
  # Initial values 
  if(is.null(init)){
  
    # obtain initial estimates for (theta,inv.tausq,rho,gamma0,gamma1) 
    theta.init=mean(data[,1])
    # moment estimator for tau
    inv.tausq.init=1/(abs(sd(data[,1])^2-(sd(data[,2])^2+mean(data[,2])^2)))
    gamma0.init=-1
    gamma1.init=0
    rho.init=0

    # List for initial values of (theta,inv.tausq,rho,gamma0,gamma1,omega,p)
    init.vals = list(theta = theta.init, 
                     inv.tausq = inv.tausq.init,
                     rho = rho.init,
                     gamma0 = gamma0.init,
                     gamma1 = gamma1.init)
    
  } else if(!is.null(init)){
    
    if(length(init) != 5)
      stop("ERROR: Please enter initial values for (theta, tau, rho, gamma0, gamma1) in this exact order.")
    if(init[2] <= 0)
      stop("ERROR: Please enter initial tau > 0.")
    if((init[3]<=-1) || (init[3]>= 1))
      stop("ERROR: Please enter initial rho between -1 and 1.")
    
    # Initialize (theta,inv.tausq,rho,gamma0,gamma1)
    theta.init = init[1]
    inv.tausq.init = (1/init[2])^2
    rho.init = max(init[3],0)

    if((init[4] > -2) & (init[4] < 2)){
       gamma0.init = init[4]
    } else if(init[4] <= -2){
      gamma0.init = -1.9999
    } else if (init[4] >= 2){
      gamma0.init = 1.9999
    }
    if(init[5] < max.s){
      gamma1.init = init[5]
    } else {
      gamma1.init = max.s-0.0001
    }
    
    # List for initial values of (theta,inv.tausq,rho,gamma0,gamma1,omega,p)
    init.vals = list(theta = theta.init, 
                     inv.tausq = inv.tausq.init,
                      rho = rho.init,
                      gamma0 = gamma0.init,
                      gamma1 = gamma1.init)
  }
  
  # Coercion
  het.dist <- match.arg(het.dist)
  
  # If normal distribution is specified for the heterogeneity
  if (het.dist == "normal"){ 
    
    #####################
    # Specify the model #
    #####################
    model_string <- "model{

    # Likelihood
    for(i in 1:n){
      # Propensity scores drawn from truncated normal to nonnegative reals
      z[i] ~ dnorm(omega[i],1) T(0,)
    
      mean.y[i] <- mu[i] + rho*s[i]*(z[i]-omega[i])
      var.y[i]<-pow(s[i],2)*(1-pow(rho,2))
     
      y[i] ~ dnorm(mean.y[i], 1/var.y[i])
    }
  
    # Prior for mu's
    for(i in 1:n){
      mu[i] ~ dnorm(theta,inv.tausq) # normal random effects
    }
  
    # Prior for theta
    theta ~ dnorm(0,0.0001)
  
    # Prior for tau2
    inv.tausq ~ dgamma(0.0001,0.0001)
    tau <- pow(1/inv.tausq,0.5)
 
    # Prior for gamma0
    gamma0 ~ dunif(-2,2)
  
    # Prior for gamma1. Truncated normal to nonnegative reals
    gamma1 ~ dunif(0,max.s)

    # Update Omega_i = gamma_0+gamma_1/s_i
    for(i in 1:n){
      omega[i] <- gamma0+gamma1/s[i]
    }
  
    # Prior for rho
    rho ~ dunif(-1,1)
 
    }"
    
    # Compile the model in JAGS
    model <- rjags::jags.model(textConnection(model_string), 
                               data = list(y=y,s=s,n=n,max.s=max.s),
                               inits = init.vals,
                               n.chains=1,
                               quiet=TRUE)
    
    
    # Draw samples
    # First use function 'update' to draw 10,000 warm-up (burnin) samples
    stats::update(model, n.iter=burn, progress.bar="none") # Burnin for 10,000 samples
    
    # Next use the function 'coda.samples' to produce the next 10,000
    # samples which will ultimately used to approximate the posterior.
    full.samples <- rjags::coda.samples(model, 
                        variable.names=c("theta","tau","rho","gamma0","gamma1"), 
                        n.iter=nmc, progress.bar="none")
    
    # Extract samples. Can be used for plotting, obtaining summary statistics, 
    # and obtaining the D measure.
    full.samples.mat <- as.matrix(full.samples[[1]])
    theta.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="theta")]
    tau.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="tau")]
    rho.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="rho")]
    gamma0.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="gamma0")]
    gamma1.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="gamma1")]
    
  } else { # Otherwise use the default Cauchy distribution for u's
   
    #####################
    # Specify the model #
    #####################
    model_string <- "model{

    # Likelihood
    for(i in 1:n){
      # Propensity scores drawn from truncated normal to nonnegative reals
      z[i] ~ dnorm(omega[i],1) T(0,)
    
      mean.y[i] <- mu[i] + rho*s[i]*(z[i]-omega[i])
      var.y[i]<-pow(s[i],2)*(1-pow(rho,2))
     
      y[i] ~ dnorm(mean.y[i], 1/var.y[i])
    }
  
    # Prior for mu's
    for(i in 1:n){
      mu[i] ~ dt(theta,inv.tausq,1) # Cauchy is a t-distribution w/ 1 df
    }
  
    # Prior for theta
    theta ~ dnorm(0,0.0001)
  
    # Prior for tau2
    inv.tausq ~ dgamma(0.0001,0.0001)
    tau <- pow(1/inv.tausq,0.5)
 
    # Prior for gamma0
    gamma0 ~ dunif(-2,2)
  
    # Prior for gamma1. Truncated normal to nonnegative reals
    gamma1 ~ dunif(0,max.s)

    # Update Omega_i = gamma_0+gamma_1/s_i
    for(i in 1:n){
      omega[i] <- gamma0+gamma1/s[i]
    }
  
    # Prior for rho
    rho ~ dunif(-1,1)
 
    }"
  
    # Compile the model in JAGS
    model <- rjags::jags.model(textConnection(model_string), 
                               data = list(y=y,s=s,n=n,max.s=max.s),
                               inits = init.vals,
                               n.chains=1,
                               quiet=TRUE)
  
  
    # Draw samples
    # First use function 'update' to draw 10,000 warm-up (burnin) samples
    stats::update(model, n.iter=burn, progress.bar="none") # Burnin for 10,000 samples

    # Next use the function 'coda.samples' to produce the next 10,000
    # samples which will ultimately used to approximate the posterior.
    full.samples <- rjags::coda.samples(model, 
                       variable.names=c("theta","tau","rho","gamma0","gamma1"), 
                       n.iter=nmc, progress.bar="none")
  
    # Extract samples. Can be used for plotting, obtaining summary statistics, 
    # and obtaining the D measure.
    full.samples.mat <- as.matrix(full.samples[[1]])
    theta.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="theta")]
    tau.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="tau")]
    rho.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="rho")]
    gamma0.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="gamma0")]
    gamma1.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="gamma1")]
  
  }
  
  # Extract point estimates
  theta.hat <- mean(theta.samples)
  tau.hat <- mean(tau.samples)
  rho.hat <- median(rho.samples)
  gamma0.hat <- median(gamma0.samples)
  gamma1.hat <- median(gamma1.samples)
  
  # Return list
  return(list(theta.hat = theta.hat, 
              theta.samples = theta.samples,
              tau.hat = tau.hat,
              tau.samples = tau.samples,
              rho.hat = rho.hat,
              rho.samples = rho.samples,
              gamma0.hat = gamma0.hat,
              gamma0.samples = gamma0.samples,
              gamma1.hat = gamma1.hat,
              gamma1.samples = gamma1.samples))
}  