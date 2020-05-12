\name{BayesNonBiasCorrected}
\alias{BayesNonBiasCorrected}

\title{Non-bias-corrected robust Bayesian meta-analysis model}

\description{This function implements the \emph{non}-bias-corrected Robust Bayesian Copas selection model of Bai et al. (2020) when there is no publication bias (i.e. \eqn{\rho=0}). If \eqn{\rho=0}, then there is no publication bias and the Copas selection model reduces to the standard random effects meta-analysis model:

\deqn{y_i = \theta + \tau u_i + s_i \epsilon_i,}

where \eqn{y_i} is the reported treatment effect for the \eqn{i}th study, \eqn{s_i} is the reported standard error for the \eqn{i}th study, \eqn{\theta} is the population treatment effect of interest, \eqn{\tau > 0} is a heterogeneity parameter, \eqn{\epsilon_i} is distributed as \eqn{N(0,1)}, and \eqn{u_i} and \eqn{\epsilon_i} are independent.

For the \emph{non}-bias-corrected model, we place noninformative priors on \eqn{(\theta, \tau^2)} (see Bai et al. (2020) for details). For the heterogeneity \eqn{u_i, i=1, \ldots, n}, the default RBC approach uses heavy-tailed standard Cauchy priors \eqn{C(0,1)}. However, we also give the option for using \eqn{u_i \sim N(0,1), i=1, \ldots, n}, for "conventional" meta-analysis.
}

\usage{
BayesNonBiasCorrected(y, s, init = NULL, het.dist=c("Cauchy","normal"),
                      burn=10000, nmc=10000)}    

\arguments{
  \item{y}{an \eqn{n \times 1} vector of reported treatment effects.}
  \item{s}{an \eqn{n \times 1} vector of reported within-study standard errors.}
  \item{init}{optional initialization values for \eqn{(\theta, \tau)}. If specified, they must be provided in this exact order. If they are not provided, the program estimates initial values from the data.}
  \item{het.dist}{Distribution for the heterogeneity \eqn{u_i, i=1, \ldots, n}. The user may specify either Cauchy or normal priors on the hetereogeneity. The default is \code{Cauchy}.}
  \item{burn}{Number of burn-in samples. Default is \code{burn=10000}.}
  \item{nmc}{Number of posterior samples to save. Default is \code{nmc=10000}.}
}

\value{The function returns a list containing the following components:
\item{theta.hat}{posterior mean for \eqn{\theta}. }
\item{theta.samples}{MCMC samples for \eqn{\theta} after burn-in. Used for plotting the posterior for \eqn{\theta} and performing inference of \eqn{\theta}.}
\item{tau.hat}{posterior mean for \eqn{\tau}. }
\item{tau.samples}{MCMC samples for \eqn{\tau} after burn-in. Used for plotting the posterior for \eqn{\tau} and performing inference of \eqn{\tau}.}
}

\references{
Bai, R., Lin, L., Boland, M. R., and Chen, Y. (2020). "A robust Bayesian Copas selection model for quantifying and correcting publication bias." \emph{arXiv preprint arXiv:2005.02930.}
}

\examples{

######################################
# Example on the Barlow2014 data set #
######################################
data(Barlow2014)
attach(Barlow2014)
# Observed treatment effect
y.obs = Barlow2014[,1]
# Observed standard error
s.obs = Barlow2014[,2]

####################################
# Fit the non-bias-corrected model #
####################################

# NOTE: Use default burn-in (burn=10000) and post-burn-in samples (nmc=10000)
RBCNoBias.mod = BayesNonBiasCorrected(y=y.obs, s=s.obs, burn=500, nmc=500)

# Point estimate for theta
theta.hat.RBCNoBias = RBCNoBias.mod$theta.hat 
# Standard error for theta
theta.se.RBCNoBias = sd(RBCNoBias.mod$theta.samples)
# 95% posterior credible interval for theta
theta.cred.int = quantile(RBCNoBias.mod$theta.samples, probs=c(0.025,0.975))

# Display results
theta.hat.RBCNoBias
theta.se.RBCNoBias
theta.cred.int

# Plot the posterior for theta
hist(RBCNoBias.mod$theta.samples)

}