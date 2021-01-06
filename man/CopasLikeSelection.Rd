\name{CopasLikeSelection}
\alias{CopasLikeSelection}

\title{Copas-like selection model}

\description{This function performs maximum likelihood estimation (MLE) of \eqn{(\theta, \tau, \rho, \gamma_0, \gamma_1)} using the EM algorithm of Ning et al. (2017) for the Copas selection model,

\deqn{y_i | (z_i>0) = \theta + \tau u_i + s_i \epsilon_i,}
\deqn{z_i = \gamma_0 + \gamma_1 / s_i + \delta_i,}
\deqn{corr(\epsilon_i, \delta_i) = \rho,}

where \eqn{y_i} is the reported treatment effect for the \eqn{i}th study, \eqn{s_i} is the reported standard error for the \eqn{i}th study, \eqn{\theta} is the population treatment effect of interest, \eqn{\tau > 0} is a heterogeneity parameter, and \eqn{u_i}, \eqn{\epsilon_i}, and \eqn{\delta_i} are marginally distributed as \eqn{N(0,1)}, and \eqn{u_i} and \eqn{\epsilon_i} are independent. 

In the Copas selection model, \eqn{y_i} is published (selected) if and only if the corresponding propensity score \eqn{z_i} (or the propensity to publish) is greater than zero. The propensity score \eqn{z_i} contains two parameters: \eqn{\gamma_0} controls the overall probability of publication, and \eqn{\gamma_1} controls how the chance of publication depends on study sample size. The reported treatment effects and propensity scores are correlated through \eqn{\rho}. If \eqn{\rho=0}, then there is no publication bias and the Copas selection model reduces to the standard random effects meta-analysis model.  

This is called the "Copas-like selection model" because to find the MLE, the EM algorithm utilizes a latent variable \eqn{m} that is treated as missing data. See Ning et al. (2017) for more details. An alternative funtion for implementing the Copas selection model using a grid search for \eqn{(\gamma_0, \gamma_1)} is available in the \code{R} package \code{metasens}.
}

\usage{
CopasLikeSelection(y, s, init = NULL, tol=1e-20, maxit=1000)
}
\arguments{
  \item{y}{An \eqn{n \times 1} vector of reported treatment effects.}
  \item{s}{An \eqn{n \times 1} vector of reported within-study standard errors.}
  \item{init}{Optional initialization values for \eqn{(\theta, \tau, \rho, \gamma_0, \gamma_1)}. If specified, they must be provided in this exact order. If they are not provided, the program estimates initial values from the data.}
  \item{tol}{Convergence criterion for the Copas-like EM algorithm for finding the MLE. Default is \code{tol=1e-20}.}
  \item{maxit}{Maximum number of iterations for the Copas-like EM algorithm for finding the MLE. Default is \code{maxit=1000}.}
}

\value{The function returns a list containing the following components:
\item{theta.hat}{MLE of \eqn{\theta}. }
\item{tau.hat}{MLE of \eqn{\tau}. }
\item{rho.hat}{MLE of \eqn{\rho}. }
\item{gamma0.hat}{MLE of \eqn{\gamma_0}.}
\item{gamma1.hat}{MLE of \eqn{\gamma_1}.}
\item{H}{\eqn{5 \times 5} Hessian matrix for the estimates of \eqn{(\theta, \tau, \rho, \gamma_0, \gamma_1)}. The square root of the diagonal entries of \eqn{H} can be used to estimate the standard errors for \eqn{(\theta, \tau, \rho, \gamma_0, \gamma_1)}. }
\item{conv}{"1" if the optimization algorithm converged, "0" if algorithm did not converge. If \code{conv=0}, then using \eqn{H} to estimate the standard errors may not be reliable.}
}

\references{
Ning, J., Chen, Y., and Piao, J. (2017). "Maximum likelihood estimation and EM algorithm of Copas-like selection model for publication bias correction." \emph{Biostatistics}, \bold{18}(3):495-504.
}

\examples{

####################################
# Example on the Hackshaw data set #
####################################
data(Hackshaw1997)
attach(Hackshaw1997)
# Extract the log OR
y.obs = Hackshaw1997[,2]
# Extract the observed standard error
s.obs = Hackshaw1997[,3]

##################################
# Fit Copas-like selection model #
##################################

# First fit RBC model with normal errors
RBC.mod = RobustBayesianCopas(y=y.obs, s=s.obs, re.dist="normal", seed=123, burn=500, nmc=500)

# Fit CLS model with initial values given from RBC model fit.
# Initialization is not necessary but the algorithm will converge faster with initialization.
CLS.mod = CopasLikeSelection(y=y.obs, s=s.obs, init=c(RBC.mod$theta.hat, RBC.mod$tau.hat,
                                                       RBC.mod$rho.hat, RBC.mod$gamma0.hat,
                                                    RBC.mod$gamma1.hat))

# Point estimate for theta 
CLS.theta.hat = CLS.mod$theta.hat  

# Use Hessian to estimate standard error for theta
CLS.Hessian = CLS.mod$H
# Standard error estimate for theta
CLS.theta.se = sqrt(CLS.Hessian[1,1]) # 

# 95 percent confidence interval 
CLS.interval = c(CLS.theta.hat-1.96*CLS.theta.se, CLS.theta.hat+1.96*CLS.theta.se)

# Display results
CLS.theta.hat  
CLS.theta.se  
CLS.interval   

# Other parameters controlling the publication bias
CLS.mod$rho.hat 
CLS.mod$gamma0.hat
CLS.mod$gamma1.hat

}