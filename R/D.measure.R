################################
## D MEASURE OF DISSIMILARITY ##
################################
# Obtains a measure of the dissimilarity in posterior distributions
# due to selection bias

# INPUTS:
# samples.selectionmodel = MCMC samples for theta under RBC model
# samples.nobiasmodel = MCMC samples for theta under no-bias model (rho=0)

# OUTPUT:
# D.measure = measure of dissimilarity in estimates from the two models

D.measure <- function(samples.RBCmodel, samples.nobiasmodel){

  fx = statip::densityfun(samples.RBCmodel)
  fy = statip::densityfun(samples.nobiasmodel)
  
  g <- function(z) (fx(z)^0.5 - fy(z)^0.5)^2
  
  # Compute Hellinger distance
  H.sq = stats::integrate(g, lower=-Inf, upper=Inf, subdivisions=500, stop.on.error=FALSE)$value/2
  H.hat = sqrt(H.sq)
  
  if(H.hat > 1){ # sometimes returns values slightly larger than 1
    H.hat = 1
  }
  
  # Return the Hellinger distance
  return(H.hat)

}