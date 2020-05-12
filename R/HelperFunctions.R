######################
## HELPER FUNCTIONS ##
######################

# for the E-step in the EM algorithm for the Copas-like selection model
#' @keywords internal
E.m=function(para,data){ # used for E step in EM 
  theta=para[1] ;   tau=para[2] ;   rho=para[3]
  gamma0=para[4];    gamma1=para[5]
  theta.hat=data[,1];s=data[,2]
  temp1=gamma0+gamma1/s
  temp2=rho*s*(theta.hat-theta)/(tau^2+s^2)
  denom = sqrt(1-rho^2*s^2/(tau^2+s^2))
  #calculating P[Z>0|theta] = \Phi[v_i]
  p1=pnorm((temp1+temp2)/denom)
  p2=pnorm(-(temp1+temp2)/denom)
  #calculating E[m]=(1-p)/p (expected # of unpublished studies)
  return(p2/p1)  
}

# the log likelihood function for the M-step in the EM algorithm
# for the Copas-like selection model
#' @keywords internal
E.loglik=function(para,data,m){ 
  theta=para[1] ;   tau=para[2] ;   rho=para[3]
  gamma0=para[4];    gamma1=para[5]
  theta.hat=data[,1]; s=data[,2]
  junk1=gamma0 + gamma1/s + rho*s*(theta.hat-theta)/(tau^2+s^2)
  junk2=sqrt(1-rho^2*s^2/(tau^2+s^2))
  v=junk1/junk2
  ###bounding V?
  v[v < -37] <- -37
  v[v > 37] <- 37
  p1=pnorm(v)
  p2=pnorm(-v)
  ###full log likelihood 
  result=log(p1)+m*log(p2)-(m+1)/2*log(tau^2+s^2)-(m+1)/2*(theta.hat-theta)^2/(tau^2+s^2)
  return(-sum(result)) ###often prefer to minimize neg. log likelihood instead of max
}


# the log likelihood function for the standard random effects meta-analysis
#' @keywords internal
RE.loglik=function(para,data){ 
  theta=para[1]    
  tau.sq=para[2]^2 
  y.obs = data[,1]
  s.obs.sq = data[,2]^2
  
  ## full log likelihood 
  result=-log(s.obs.sq+tau.sq)-(y.obs-theta)^2/(s.obs.sq+tau.sq) 
  return(-sum(result)) ###often prefer to minimize neg. log likelihood instead of max
}


## Estimate optimal gamma0, gamma1 given other parameters (theta,tau,rho)
#' @keywords internal
optimal.gammas=function(data,theta,tau,rho){
  
  # Initial (gamma0,gamma1)
  gamma.params = c(-1,0.1)
  # Maximize w.r.t. (gamma0,gamma1)
  output = stats::optim(gamma.params,loglik.Copas,data=data,theta=theta,tau=tau,
                 rho=rho,method="L-BFGS-B",lower=c(-2,0.01), 
                 upper=c(1.5,0.9), control = list(maxit = 1000),hessian=TRUE)
  
  return(list(gamma0.hat=output$par[1],
              gamma1.hat=output$par[2],
              loglik=(-1)*output$val))
}

## Log-likelihood of original Copas selection model
#' @keywords internal
loglik.Copas=function(gamma.params,theta,tau,rho,data){
  
  theta.hat=data[,1]; s=data[,2];
  gamma0=gamma.params[1]
  gamma1=gamma.params[2]
  
  ## Compute v
  junk1 = gamma0+(gamma1/s)+(rho*s*(theta.hat-theta))/(tau^2+s^2)
  junk2 = sqrt(1-(rho^2*s^2)/(tau^2+s^2))
  ## Bound v
  v=junk1/junk2
  v[v < -37] <- -37
  v[v > 37] <- 37
  
  term1 = -0.5*log(tau^2+s^2)
  term2 = -((theta.hat-theta)^2)/(2*(tau^2+s^2))
  term3 = -log(pnorm(gamma0+gamma1/s))
  term4 = log(pnorm(v))
  
  # Minimize the negative of the sum
  return(-sum(term1+term2+term3+term4))
}