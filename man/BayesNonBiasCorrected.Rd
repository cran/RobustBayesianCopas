\name{BayesNonBiasCorrected}
\alias{BayesNonBiasCorrected}

\title{Non-bias-corrected robust Bayesian meta-analysis model}

\description{This function implements the \emph{non}-bias-corrected Robust Bayesian Copas selection model of Bai et al. (2020) when there is no publication bias (i.e. \eqn{\rho=0}). In this case, the Copas selection model reduces to the standard random effects meta-analysis model:

\deqn{y_i = \theta + \tau u_i + s_i \epsilon_i,}

where \eqn{y_i} is the reported treatment effect for the \eqn{i}th study, \eqn{s_i} is the reported standard error for the \eqn{i}th study, \eqn{\theta} is the population treatment effect of interest, \eqn{\tau > 0} is a heterogeneity parameter, \eqn{\epsilon_i} is distributed as \eqn{N(0,1)}, and \eqn{u_i} and \eqn{\epsilon_i} are independent.

For the \emph{non}-bias-corrected model, we place noninformative priors on \eqn{(\theta, \tau^2)} (see Bai et al. (2020) for details). For the random effects \eqn{u_i, i=1, \ldots, n}, we give the option for using normal, Student's t, Laplace, or slash distributions for the random effects. If this function is being run in order to quantify publication bias with the robust Bayesian Copas selection model, then the practitioner should use the \emph{same} random effects distribution that they used for \code{RobustBayesianCopas}.  
}

\usage{
BayesNonBiasCorrected(y, s, re.dist=c("normal", "StudentsT", "Laplace", "slash"),
                      t.df = 4, slash.shape = 1, init=NULL, seed=NULL,
                      burn=10000, nmc=10000)}    

\arguments{
  \item{y}{An \eqn{n \times 1} vector of reported treatment effects.}
  \item{s}{An \eqn{n \times 1} vector of reported within-study standard errors.}
  \item{re.dist}{Distribution for the between-study random effects \eqn{u_i, i=1, \ldots, n}. The user may specify the normal, Student's t, Laplace, or slash distribution. The default is \code{StudentsT} with 4 degrees of freedom.}
    \item{t.df}{Degrees of freedom for t-distribution. Only used if \code{StudentsT} is specified for \code{re.dist}. Default is 4.}
    \item{slash.shape}{Shape parameter in the slash distribution. Only used if \code{slash} is specified for \code{re.dist}. Default is 1.}
    \item{init}{Optional initialization values for \eqn{(\theta, \tau)}. If specified, they must be provided in this exact order. If they are not provided, the program estimates initial values from the data.}
    \item{seed}{Optional seed. This needs to be specified if you want to reproduce the exact results of your analysis. }
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
# Fit the model with Laplace errors.
RBCNoBias.mod = BayesNonBiasCorrected(y=y.obs, s=s.obs, re.dist="Laplace", burn=500, nmc=500)

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