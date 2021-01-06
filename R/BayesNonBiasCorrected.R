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
# re.dist = distribution of random effects. Allows either "normal", "t", "Laplace", or "slash"
# t.df = degrees of freedom for t-distribution. Only used if "t" is specified for re.dist. 
#        Default is t=4.
# slash.shape = shape parameter in the slash distribution. Only used if "slash" is specified
#               for re.dist. Default is slash.shape=1.  
# init = optional initialization values of (theta, tau)
#        If the user does not provide, these are estimated from the data.
# burn = number of burn-in samples. Default is 10,000
# nmc = number of posterior draws to be saved. Default is 10,000

# OUTPUTS:
# theta.hat = point estimate (posterior mean) for theta
# theta.samples = full MCMC samples for theta
# tau.hat = point estimate (posterior mean) for tau
# tau.samples = full MCMC samples for tau

BayesNonBiasCorrected <- function(y = y, s = s, re.dist = c("normal", "StudentsT", "Laplace", "slash"),
                                  t.df = 4, slash.shape = 1, init=NULL, seed=NULL,
                                  burn=10000, nmc=10000){
  
  # Check that y and s are of the same length
  if(length(y) != length(s))
    stop("Please ensure that 'y' and 's' are of the same length.")
    
  # Sample size
  n = length(y)
  
  # Coersion
  re.dist <- match.arg(re.dist)
  
  # If Student's t is used, make sure that degrees of freedom is greater than or equal to 3 
  if(re.dist == "StudentsT"){
    if( (t.df-floor(t.df) != 0) || (t.df < 3) )
      stop("ERROR: Please enter degrees of freedom as an integer greater than or equal to 3.")
  }
  
  # If slash is used, make sure the slash shape parameter is greater than 0
  if(re.dist == "slash"){
    if(slash.shape <= 0)
      stop("ERROR: Please enter shape parameter strictly greater than 0.")
  }
  
  # Initial values 
  if(is.null(init)){
  
    ### obtain initial estimators using the MLE for 
    SMA = StandardMetaAnalysis(y, s)
    
    # List for initial values of (theta, tau)
    theta.init = SMA$theta.hat
    tau.init = SMA$tau.hat
    
  } else if(!is.null(init)){
    
    if(length(init)!=2)
      stop("If providing initial values, please ensure they are provided for (theta, tau) in this exact order.")
    
    # List for initial values of (theta, tau)
    theta.init = init[1]
    tau.init = init[2]
  }
  
  if(is.null(seed)){  
    
    # if no seed is specified
    init.vals = list(theta = theta.init, 
                     tau = tau.init)
    
  } else if(!is.null(seed)){
    
    # if a seed is specified
    seed <- as.integer(seed)
    
    init.vals = list(.RNG.name = "base::Wichmann-Hill",
                     .RNG.seed = seed,
                     theta = theta.init, 
                     tau = tau.init)
  }
  
  # Fit JAGS models
  
  if(re.dist=="normal") { 
    
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
      mu[i] ~ dnorm(theta,1/pow(tau,2)) # normal distribution on random effects
    }
  
    # Prior for theta
    theta ~ dnorm(0,0.0001)
  
    # Prior for tau. Half-Cauchy is a truncated standard t w/ 1 degree of freedom
    tau ~ dt(0,1,1) T(0,)

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
  
  } else if(re.dist=="Laplace"){
    
    model_string <- "model{

    # Likelihood
    for(i in 1:n){
      y[i] ~ dnorm(mu[i], 1/pow(s[i],2))
    }
  
    # Prior for mu[i]'s
    for(i in 1:n){
      mu[i] ~ ddexp(theta,1/tau) # Laplace-distributed random effects
    }
  
    # Prior for theta
    theta ~ dnorm(0,0.0001)
  
    # Prior for tau. Half-Cauchy is a truncated standard t w/ 1 degree of freedom
    tau ~ dt(0,1,1) T(0,)

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
  
  } else if(re.dist=="slash"){
    
    model_string <- "model{

    # Likelihood
    for(i in 1:n){
      y[i] ~ dnorm(mu[i], 1/pow(s[i],2))
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

    }"
    
    # Compile the model in JAGS
    model <- rjags::jags.model(textConnection(model_string), 
                               data = list(y=y,s=s,n=n,slash.shape=slash.shape),
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
  
  } else {  # Default is t-distribution 
      
      model_string <- "model{

    # Likelihood
    for(i in 1:n){
      y[i] ~ dnorm(mu[i], 1/pow(s[i],2))
    }
  
    # Prior for mu[i]'s
    for(i in 1:n){
      mu[i] ~ dt(theta,1/pow(tau,2),t.df) # t-distributed random effects
    }
  
    # Prior for theta
    theta ~ dnorm(0,0.0001)
  
    # Prior for tau. Half-Cauchy is a truncated standard t w/ 1 degree of freedom
    tau ~ dt(0,1,1) T(0,)

    }"
      
      # Compile the model in JAGS
      model <- rjags::jags.model(textConnection(model_string), 
                                 data = list(y=y,s=s,n=n,t.df=t.df),
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
  }
    
  # Extract samples. Can be used for plotting, obtaining summary statistics, 
  # and obtaining the D measure.
  full.samples.mat <- as.matrix(full.samples[[1]])
  theta.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="theta")]
  tau.samples <- full.samples.mat[,which(colnames(full.samples.mat)=="tau")]
 
  # Extract point estimates
  theta.hat <- mean(theta.samples)
  tau.hat <- mean(tau.samples)
  
  # Return list   
  return(list(theta.hat = theta.hat, 
              theta.samples = theta.samples,
              tau.hat = tau.hat,
              tau.samples = tau.samples))
}  