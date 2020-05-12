\name{StandardMetaAnalysis}
\alias{StandardMetaAnalysis}

\title{Standard meta-analysis}

\description{This function performs maximum likelihood estimation (MLE) of \eqn{(\theta, \tau)} for the standard random effects meta-analysis model,

\deqn{y_i = \theta + \tau u_i + s_i \epsilon_i,} 

where \eqn{y_i} is the reported treatment effect for the \eqn{i}th study, \eqn{s_i} is the reported standard error for the \eqn{i}th study, \eqn{\theta} is the population treatment effect of interest, \eqn{\tau > 0} is a heterogeneity parameter, and \eqn{u_i} and \eqn{\epsilon_i} are independent and distributed as \eqn{N(0,1)}.
}

\usage{
StandardMetaAnalysis(y, s, init = NULL, tol=1e-10, maxit=1000)
}

\arguments{
  \item{y}{an \eqn{n \times 1} vector of reported treatment effects.}
  \item{s}{an \eqn{n \times 1} vector of reported within-study standard errors.}
  \item{init}{optional initialization values for \eqn{(\theta, \tau)}. If specified, they must be provided in this order. If they are not provided, the program estimates initial values from the data.}
  \item{tol}{Convergence criterion for the optimization algorithm for finding the MLE. Default is \code{tol=1e-10}.}
  \item{maxit}{Maximum number of iterations for the optimization algorithm for finding the MLE. Default is \code{maxit=1000}.}
}

\value{The function returns a list containing the following components:
\item{theta.hat}{MLE of \eqn{\theta}. }
\item{tau.hat}{MLE of \eqn{\tau}. }
\item{H}{\eqn{2 \times 2} Hessian matrix for the estimates of \eqn{(\theta, \tau)}. The square root of the diagonal entries of \eqn{H} can be used to estimate the standard errors for \eqn{(\theta, \tau)}. }
\item{conv}{"1" if the optimization algorithm converged, "0" if algorithm did not converge. If \code{conv=0}, then using \eqn{H} to estimate the standard errors may not be reliable.}
}

\references{
Bai, R., Lin, L., Boland, M. R., and Chen, Y. (2020). "A robust Bayesian Copas selection model for quantifying and correcting publication bias." \emph{arXiv preprint arXiv:2005.02930.}

Ning, J., Chen, Y., and Piao, J. (2017). "Maximum likelihood estimation and EM algorithm of Copas-like selection model for publication bias correction." \emph{Biostatistics}, \bold{18}(3):495-504.
}

\examples{

############################################
# Example on the antidepressants data set. #
# This is from Section 7.2 of the paper by #
# Bai et al. (2020).                       #
############################################
# Load the full data
data(antidepressants)
attach(antidepressants)

# Extract the 50 published studies
published.data = antidepressants[which(antidepressants$Published==1),]
# Observed treatment effect
y.obs = published.data$Standardized_effect_size
# Observed standard error
s.obs = published.data$Standardized_SE

#################################
# Fit a standard meta-analysis  #
# that ignores publication bias #
#################################
SMA.mod = StandardMetaAnalysis(y=y.obs, s=s.obs)

# Point estimate for theta
SMA.theta.hat = SMA.mod$theta.hat
# Use Hessian to estimate standard error for theta
SMA.Hessian = SMA.mod$H
# Standard error estimate for theta
SMA.theta.se = sqrt(SMA.Hessian[1,1]) 
# 95 percent confidence interval 
SMA.interval = c(SMA.theta.hat-1.96*SMA.theta.se, SMA.theta.hat+1.96*SMA.theta.se)

# Display results
SMA.theta.hat
SMA.theta.se
SMA.interval

\donttest{
######################################
# Example on the Barlow2014 data set #
######################################
data(Barlow2014)
attach(Barlow2014)
# Observed treatment effect
y.obs = Barlow2014[,1]
# Observed standard error
s.obs = Barlow2014[,2]

#################################
# Fit a standard meta-analysis  #
# that ignores publication bias #
#################################
SMA.mod = StandardMetaAnalysis(y=y.obs, s=s.obs)

# Point estimate for theta
SMA.theta.hat = SMA.mod$theta.hat

# Use Hessian to estimate standard error for theta
SMA.Hessian = SMA.mod$H
# Standard error estimate for theta
SMA.theta.se = sqrt(SMA.Hessian[1,1]) 

# 95 percent confidence interval 
SMA.interval = c(SMA.theta.hat-1.96*SMA.theta.se, SMA.theta.hat+1.96*SMA.theta.se)

# Display results
SMA.theta.hat
SMA.theta.se
SMA.interval
}
}