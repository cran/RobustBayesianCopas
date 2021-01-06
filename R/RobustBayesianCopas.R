#################################################
## ROBUST BAYESIAN COPAS (RBC) SELECTION MODEL ##
#################################################
# Implements the full robust Bayesian Copas (RBC) selection model 
# of Bai et al. (2020): https://arxiv.org/abs/2005.02930

# INPUTS:
# y = given vector of observed treatment effects
# s = given vector of within-study standard errors 
# re.dist = distribution of random effects. Allows either "normal", "t", "Laplace", or "slash"
# t.df = degrees of freedom for t-distribution. Only used if "t" is specified for re.dist. 
#        Default is t.df=4.
# slash.shape = shape parameter in the slash distribution. Only used if "slash" is specified
#               for re.dist. Default is slash.shape=1.  
# init = optional initialization values of (theta, tau, rho, gamma0, gamma1)
#        If the user does not provide, these are estimated from the data.
# seed = seed. Set this if you want to reproduce exact results
# burn = number of burn-in samples. Default is 10,000
# nmc = number of posterior draws to be saved. Default is 10,000

# OUTPUTS:
# DIC = deviance information criterion (used for model selection)
# theta.hat = point estimate for theta. Returns posterior mean of theta
# theta.samples = theta samples. Used for plotting the posterior for theta
#                 and performing inference for theta
# tau.hat = point estimate for tau. Returns posterior mean of tau
# tau.samples = tau samples. Used for plotting the posterior for tau
#               and performing inference for tau 
# rho.hat = point estimate for rho. Returns posterior median of rho
# rho.samples = rho samples. Used for plotting the posterior for rho
#               and performing inference for rho
# gamma0.hat = point estimate for gamma0. Returns posterior median of gamma0
# gamma0.samples = gamma0 samples. Used for plotting the posterior for gamma0
#                  and performing inference for gamma0 
# gamma1.hat = point estimate for gamma1. Returns posterior median of gamma1
# gamma1.samples = gamma1 samples. Used for plotting the posterior for gamma1
#                  and performing inference for gamma1 

RobustBayesianCopas <- function(y, s, re.dist=c("normal", "StudentsT", "Laplace", "slash"),
                                t.df = 4, slash.shape = 1, init=NULL, seed=NULL, 
                                burn=10000, nmc=10000){
  
  # Check that y and s are of the same length
  if(length(y) != length(s))
    stop("ERROR: Please ensure that 'y' and 's' are of the same length.")
  if(any(s <= 0))
    stop("ERROR: Please ensure that all standard errors are strictly positive.")
  
  # Sample size
  n = length(y)
  data = data.frame(cbind(y,s))
  max.s = max(s)
  
  # Coersion
  re.dist <- match.arg(re.dist)
  
  # If Student's t is used, make sure that degrees of freedom is greater than or equal to 1 
  if(re.dist == "StudentsT"){
    if( (t.df-floor(t.df) != 0) || (t.df < 1) )
      stop("ERROR: Please enter degrees of freedom as an integer greater than or equal to 1.")
  }
  
  # If slash is used, make sure the slash shape parameter is greater than 0
  if(re.dist == "slash"){
    if(slash.shape <= 0)
      stop("ERROR: Please enter shape parameter strictly greater than 0.")
  }
  
  # Initial values 
  if(is.null(init)){
  
    # obtain initial estimates for (theta,inv.tausq,rho,gamma0,gamma1) 
    theta.init=mean(data[,1])
    # moment estimator for tau
    tau.init=abs(sd(data[,1])^2-(sd(data[,2])^2+mean(data[,2])^2))
    gamma0.init=-1
    gamma1.init=0
    rho.init=0
  
  } else if(!is.null(init)){
    
    if(length(init) != 5)
      stop("ERROR: Please enter initial values for (theta, tau, rho, gamma0, gamma1) in that exact order.")
    if(init[2] < 0)
      stop("ERROR: Please enter initial tau > 0.")
    
    if(init[3]>=1){
      init[3] = 0.99
    } else if(init[3]<=-1){
      init[3] = -0.99
    }
    
    # Initialize (theta,inv.tausq,rho,gamma0,gamma1)
    theta.init = init[1]
    tau.init = init[2]
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
  }
  
  # List for initial values of (theta,inv.tausq,rho,gamma0,gamma1)
  
  if(is.null(seed)){  
  
    # if no seed is specified
    init.vals = list(theta = theta.init, 
                     tau = tau.init,
                     rho = rho.init,
                     gamma0 = gamma0.init,
                     gamma1 = gamma1.init)
 
  } else if(!is.null(seed)){
    
    # if a seed is specified
    seed <- as.integer(seed)
    
    init.vals = list(.RNG.name = "base::Wichmann-Hill",
                     .RNG.seed = seed,
                     theta = theta.init, 
                     tau = tau.init,
                     rho = rho.init,
                     gamma0 = gamma0.init,
                     gamma1 = gamma1.init)
  }
  
  # Fit JAGS models
  
  if(re.dist=="normal") { 
    
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
      mu[i] ~ dnorm(theta,1/pow(tau,2)) # normal random effects
    }
  
    # Prior for theta
    theta ~ dnorm(0,0.0001)
  
    # Prior for tau. Half-Cauchy is a truncated standard t w/ 1 degree of freedom
    tau ~ dt(0,1,1) T(0,)
    
    # Prior for gamma0
    gamma0 ~ dunif(-2,2)
  
    # Prior for gamma1. Truncated normal to nonnegative reals
    gamma1 ~ dunif(0,max.s)

    # Update omega_i = gamma_0+gamma_1/s_i
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
                                        variable.names=c("mu","theta","tau","rho","gamma0","gamma1"), 
                                        n.iter=nmc, progress.bar="none")
    
  } else if(re.dist == "Laplace") {
  
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
      mu[i] ~ ddexp(theta,1/tau) # Laplace-distributed random effects
    }
  
    # Prior for theta
    theta ~ dnorm(0,0.0001)
  
    # Prior for tau. Half-Cauchy is a truncated standard t w/ 1 degree of freedom
    tau ~ dt(0,1,1) T(0,)
 
    # Prior for gamma0
    gamma0 ~ dunif(-2,2)
  
    # Prior for gamma1. Truncated normal to nonnegative reals
    gamma1 ~ dunif(0,max.s)

    # Update omega_i = gamma_0+gamma_1/s_i
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
                                        variable.names=c("mu","theta","tau","rho","gamma0","gamma1"), 
                                        n.iter=nmc, progress.bar="none")
    
  } else if(re.dist=="slash") { 
    
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
      mu[i] ~ dnorm(theta, w[i]/pow(tau,2)) # slash random effects
      w[i] ~ dbeta(slash.shape,1)
    }
  
    # Prior for theta
    theta ~ dnorm(0,0.0001)
  
    # Prior for tau. Half-Cauchy is a truncated standard t w/ 1 degree of freedom
    tau ~ dt(0,1,1) T(0,)
    
    # Prior for gamma0
    gamma0 ~ dunif(-2,2)
  
    # Prior for gamma1. Truncated normal to nonnegative reals
    gamma1 ~ dunif(0,max.s)

    # Update omega_i = gamma_0+gamma_1/s_i
    for(i in 1:n){
      omega[i] <- gamma0+gamma1/s[i]
    }
  
    # Prior for rho
    rho ~ dunif(-1,1)
 
    }"
    
    # Compile the model in JAGS
    model <- rjags::jags.model(textConnection(model_string), 
                               data = list(y=y,s=s,n=n,max.s=max.s,slash.shape=slash.shape),
                               inits = init.vals,
                               n.chains=1,
                               quiet=TRUE)
    
    # Draw samples
    # First use function 'update' to draw 10,000 warm-up (burnin) samples
    stats::update(model, n.iter=burn, progress.bar="none") # Burnin for 10,000 samples
    
    # Next use the function 'coda.samples' to produce the next 10,000
    # samples which will ultimately used to approximate the posterior.
    full.samples <- rjags::coda.samples(model, 
                                        variable.names=c("mu","theta","tau","rho","gamma0","gamma1"), 
                                        n.iter=nmc, progress.bar="none")
    
  } else {   # Default is t-distribution 
    
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
      mu[i] ~ dt(theta,1/pow(tau,2),t.df) # t-distributed random effects
    }
  
    # Prior for theta
    theta ~ dnorm(0,0.0001)
  
    # Prior for tau. Half-Cauchy is a truncated standard t w/ 1 degree of freedom
    tau ~ dt(0,1,1) T(0,)
 
    # Prior for gamma0
    gamma0 ~ dunif(-2,2)
  
    # Prior for gamma1. Truncated normal to nonnegative reals
    gamma1 ~ dunif(0,max.s)

    # Update omega_i = gamma_0+gamma_1/s_i
    for(i in 1:n){
      omega[i] <- gamma0+gamma1/s[i]
    }
  
    # Prior for rho
    rho ~ dunif(-1,1)
    
    }"
    
    # Compile the model in JAGS
    model <- rjags::jags.model(textConnection(model_string), 
                               data = list(y=y,s=s,n=n,max.s=max.s,t.df=t.df),
                               inits = init.vals,
                               n.chains=1,
                               quiet=TRUE)
    
    # Draw samples
    # First use function 'update' to draw 10,000 warm-up (burnin) samples
    stats::update(model, n.iter=burn, progress.bar="none") # Burnin for 10,000 samples
    
    # Next use the function 'coda.samples' to produce the next 10,000
    # samples which will ultimately used to approximate the posterior.
    full.samples <- rjags::coda.samples(model, 
                                        variable.names=c("mu","theta","tau","rho","gamma0","gamma1"), 
                                        n.iter=nmc, progress.bar="none")
    
  }
    
  # Extract samples. Can be used for plotting, obtaining summary statistics, 
  # and obtaining the D measure.
  full.samples.mat <- as.matrix(full.samples[[1]])
  mu.samples <- full.samples.mat[,c(1:n)]
  theta.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="theta")]
  tau.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="tau")]
  rho.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="rho")]
  gamma0.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="gamma0")]
  gamma1.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="gamma1")]
  
  # Extract point estimates
  mu.hat <- colMeans(mu.samples)
  theta.hat <- mean(theta.samples)
  tau.hat <- mean(tau.samples)
  
  # For computing DIC, should use the mean
  rho.mean <- mean(rho.samples)
  gamma0.mean <- mean(gamma0.samples)
  gamma1.mean <- mean(gamma1.samples)
  
  # Point estimates are the median
  rho.hat <- median(rho.samples)
  gamma0.hat <- median(gamma0.samples)
  gamma1.hat <- median(gamma1.samples)
  
  # Calculate the DIC
  v.samples = matrix(0, nrow=nmc, ncol=n)
  dev.summand.samples = matrix(0, nrow=nmc, ncol=n)

  for(i in 1:n){
    v.samples[,i] = (gamma0.samples + (gamma1.samples+rho.samples*(y[i]-mu.samples[,i]))/s[i])/(sqrt(1-rho.samples))
    v.samples[which(v.samples[,i] < -37),i] = -37
    v.samples[which(v.samples[,i] > 37),i] = 37
    dev.summand.samples[,i] = (y[i]-mu.samples[,i])^2/s[i]^2+2*log(pnorm(gamma0.samples+gamma1.samples/s[i]))-2*log(pnorm(v.samples[,i]))
  }
  dev.samples = rowSums(dev.summand.samples)
  
  v.hat = (gamma0.mean+(gamma1.mean+rho.mean*(y-mu.hat))/s)/(sqrt(1-rho.mean^2))
  v.hat[v.hat < -37] = -37
  v.hat[v.hat > 37] = 37
  dev.barparam = sum((y-mu.hat)^2/s^2 + 2*log(pnorm(gamma0.mean+gamma1.mean/s)) -2*log(pnorm(v.hat)))   
  
  DIC = 2*mean(dev.samples)-dev.barparam 
  
  # Return list
  return(list(DIC = DIC,
              theta.hat = theta.hat, 
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