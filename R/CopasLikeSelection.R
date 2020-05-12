################################
## COPAS-LIKE SELECTION MODEL ##
################################
# Obtains maximum likelihood estimates based on the Copas-like 
# selection model in Ning et al. (2017)

# INPUTS:
# y = given vector of observed treatment effects
# s = given vector of within-study standard errors
# init = optional initialization values for (theta, tau, rho, gamma0, gamma1) 
#        If the user does not provide, these are estimated from the data.
# maxit = maximum number of iterations in the optimization

# OUTPUTS:
# theta.hat = estimate for theta
# tau.hat = estimate for tau
# rho.hat = estimate for rho
# gamma0.hat = estimate for gamma0
# gamma1.hat = estimate for gamma1
# H = estimate of Hessian matrix for (theta.hat, tau.hat, rho.hat, gamma0.hat, gamma1.hat)
#     Square root of diagonal entries can be used to estimate standard errors
# conv = convergence. "1"=converged, "0"=failed to converge

CopasLikeSelection <- function(y, s, init = NULL, tol=1e-20, maxit=1000){
  
  # Check that y and s are of the same length
  if(length(y) != length(s))
      stop("Please ensure that 'y' and 's' are of the same length.")
  # Check that all values in s are > 0
  if(any(s <= 0))
      stop("Please ensure that all observed standard errors 's' are greater than zero.")
  
  # Bind 'y' and 's'
  data = data.frame(cbind(y,s))
  max.s = max(s)
  
  # If no initial values were specified
  if(is.null(init)){
  
    ### for initial estimates of (1) theta (2) tau (3) rho (4) gamma0 and (5) gamma1
    init = rep(0,5)
    max.s = max(s)
    
    # initial estimate of theta
    init[1]=mean(data[,1])
    # moment estimator for tau
    init[2]=sqrt(abs(sd(data[,1])^2-(sd(data[,2])^2+mean(data[,2])^2)))
    
    rho.test=seq(-0.95,0.95,by=0.05)
    
    test = vector(mode = "list", length = length(rho.test))
    for(k in 1:length(test)){
      test[[k]]=optimal.gammas(data,theta=init[1],
                               tau=init[2],rho=rho.test[k])
    }
    # Find the rho that maximizes the log-likelihood
    loglik.max.index = 1
    loglik.max = as.double(test[[1]][2])
    for(k in 2:length(test)){
      if(as.double(test[[k]][2]) > loglik.max){
        loglik.max = as.double(test[[k]][2])
        loglik.max.index = k
      }
    }  
    # initial value fo rho
    init[3] = rho.test[loglik.max.index]
    # initial value for gamma0
    init[4] = as.double(test[[loglik.max.index]][1])
    # initial value for gamma1
    init[5] = as.double(test[[loglik.max.index]][2])
  }
  
  # If initial values were specified but not of correct length
  if(!is.null(init)){
    if(length(init) != 5)
      stop("ERROR: Please enter initial values for (theta, tau, rho, gamma0, gamma1) in this exact order.")
    if(init[2] <= 0)
      stop("ERROR: Please enter initial tau > 0.")
    if((init[3]<=-1) || (init[3]>= 1))
      stop("ERROR: Please enter initial rho between -1 and 1.")
  }

  par.new = init # Initialize values
  counter = 0    # start counter
  
  while (counter <= maxit) {
  
    counter = counter + 1
    par.old = par.new
    
    # E-step for EM algorithm
    myM= E.m(par.old,data)  

    # M-step for EM algorithm
    output = stats::optim(par.old,E.loglik,data=data,m=myM,method="L-BFGS-B",lower=c(-Inf,1e-3,-0.99,-Inf,0), 
                   upper=c(Inf,2,0.99,Inf,Inf), control = list(maxit = 1000),hessian=TRUE)
    
    
    par.new = output$par
    
    #cat("iter = ", counter, "P = ", par.new, "\n")
    
    if (max(abs(par.new-par.old)) < tol){  ###if convergence achieved
      
      # cat("\nSuccessfully Converged\n")
      # cat("iter = ", counter, "par = ", par.new, "\n")
      
      theta.hat = par.new[1]
      tau.hat = par.new[2]
      rho.hat = par.new[3]
      gamma0.hat = par.new[4]
      gamma1.hat = par.new[5]

      # Extract standard errors for theta and tau
      H = solve(output$hessian)
      
      return(list(theta.hat = theta.hat, 
                  tau.hat = tau.hat,
                  rho.hat = rho.hat,
                  gamma0.hat = gamma0.hat,
                  gamma1.hat = gamma1.hat,
                  H = H,
                  conv=1))
    
    }
  }
  
  # print("Convergence Failed")
  theta.hat = par.new[1]
  tau.hat = par.new[2]
  rho.hat = par.new[3]
  gamma0.hat = par.new[4]
  gamma1.hat = par.new[5]
  H = solve(output$hessian)

  return(list(theta.hat = theta.hat, 
              tau.hat = tau.hat,
              rho.hat = rho.hat,
              gamma0.hat = gamma0.hat,
              gamma1.hat = gamma1.hat,
              H = H,
              conv=0)) 
}  