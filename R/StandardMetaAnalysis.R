############################
## STANDARD META-ANALYSIS ##
############################
# Obtains maximum likelihood estimate (MLE) for (theta, tau)

# INPUTS:
# y = given vector of observed treatment effects
# s = given vector of within-study standard errors
# init = optional initialization values for (theta, tau)
#        If the user does not provide, these are estimated from the data.
# maxit = maximum number of iterations in the optimization

# OUTPUTS:
# theta.hat = estimate for theta
# tau.hat = estimate for tau
# H = estimate of Hessian matrix for (theta.hat, tau.hat).
#     Square root of diagonal entries can be used to estimate standard errors.
# conv = convergence. "1"=converged, "0"=failed to converge

StandardMetaAnalysis <- function(y, s, init = NULL, tol=1e-10, maxit=1000){
  
  # Check that y and s are of the same length
  if(length(y) != length(s))
    stop("Please ensure that 'y' and 's' are of the same length.")
  # Check that all values in s are > 0
  if(any(s <= 0))
    stop("Please ensure that all observed standard errors 's' are greater than zero.")
  
  # Bind 'y' and 's'
  data = data.frame(cbind(y,s))
  
  # If no initial values were specified
  if(is.null(init)){
  
    # initial estimates of (1) theta (2) tau
    init[1]=mean(data[,1])
    # moment estimator for tau
    init[2]=sqrt(abs(sd(data[,1])^2-(sd(data[,2])^2+mean(data[,2])^2)))
  }
  
  # If initial values were specified but not of correct length
  if(!is.null(init)){
    if(length(init) != 2)
      stop("Please enter a vector of length 2 for initial values of (theta, tau).")
  }
  
  # Maximize the log-likelihood
  output = stats::optim(init, RE.loglik, data=data, method="L-BFGS-B", lower=c(-Inf,1e-3), 
                        upper=c(Inf,Inf), control = list(maxit = maxit, factr=tol), 
                        hessian=TRUE)
  
  # Estimates
  par.est = output$par
  conv.iters = output$convergence
  if(conv.iters < maxit)
    conv = 1 ## convergence reached
  else
    conv = 0
  
  # Point estimates
  theta.hat = par.est[1]
  tau.hat = par.est[2]
  
  # Hessian matrix for (theta.hat, tau.hat)
  H = solve(output$hessian)
  
  # Return list
  return(list(theta.hat = theta.hat, 
              tau.hat = tau.hat,
              H = H,
              conv=conv)) 
} 