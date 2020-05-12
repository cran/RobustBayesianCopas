\name{RobustBayesianCopas}
\alias{RobustBayesianCopas}

\title{Robust Bayesian Copas selection model}

\description{This function implements the Robust Bayesian Copas selection model of Bai et al. (2020) for the Copas selection model,

\deqn{y_i | (z_i>0) = \theta + \tau u_i + s_i \epsilon_i,}
\deqn{z_i = \gamma_0 + \gamma_1 / s_i + \delta_i,}
\deqn{corr(\epsilon_i, \delta_i) = \rho,}

where \eqn{y_i} is the reported treatment effect for the \eqn{i}th study, \eqn{s_i} is the reported standard error for the \eqn{i}th study, \eqn{\theta} is the population treatment effect of interest, \eqn{\tau > 0} is a heterogeneity parameter, and \eqn{\epsilon_i}, and \eqn{\delta_i} are marginally distributed as \eqn{N(0,1)} and \eqn{u_i}, and \eqn{\epsilon_i} are independent. 

In the Copas selection model, \eqn{y_i} is published (selected) if and only if the corresponding propensity score \eqn{z_i} (or the propensity to publish) is greater than zero. The propensity score \eqn{z_i} contains two parameters: \eqn{\gamma_0} controls the overall probability of publication, and \eqn{\gamma_1} controls how the probability of publication depends on study sample size. The reported treatment effects and propensity scores are correlated through \eqn{\rho}. If \eqn{\rho=0}, then there is no publication bias and the Copas selection model reduces to the standard random effects meta-analysis model.  

The RBC model places noninformative priors on \eqn{(\theta, \tau^2, \rho, \gamma_0, \gamma_1)} (see Bai et al. (2020) for details). For the heterogeneity \eqn{u_i, i=1, \ldots, n}, the default RBC approach uses heavy-tailed standard Cauchy priors \eqn{C(0,1)}. However, we also give the option for using \eqn{u_i \sim N(0,1), i=1, \ldots, n}, for "conventional" meta-analysis.
}

\usage{
RobustBayesianCopas(y, s, init=NULL, het.dist=c("Cauchy","normal"), 
                    burn=10000, nmc=10000)}

\arguments{
  \item{y}{an \eqn{n \times 1} vector of reported treatment effects.}
  \item{s}{an \eqn{n \times 1} vector of reported within-study standard errors.}
  \item{init}{optional initialization values for \eqn{(\theta, \tau, \rho, \gamma_0, \gamma_1)}. If specified, they must be provided in this exact order. If they are not provided, the program estimates initial values from the data.}
  \item{het.dist}{Distribution for the heterogeneity \eqn{u_i, i=1, \ldots, n}. The user may specify either Cauchy or normal priors on the hetereogeneity. The default is \code{Cauchy}.}
  \item{burn}{Number of burn-in samples. Default is \code{burn=10000}.}
  \item{nmc}{Number of posterior samples to save. Default is \code{nmc=10000}.}
}

\value{The function returns a list containing the following components:
\item{theta.hat}{posterior mean for \eqn{\theta}. }
\item{theta.samples}{MCMC samples for \eqn{\theta} after burn-in. Used for plotting the posterior for \eqn{\theta} and performing inference of \eqn{\theta}.}
\item{tau.hat}{posterior mean for \eqn{\tau}. }
\item{tau.samples}{MCMC samples for \eqn{\tau} after burn-in. Used for plotting the posterior for \eqn{\tau} and performing inference of \eqn{\tau}.}
\item{rho.hat}{posterior median for \eqn{\rho}. }
\item{rho.samples}{MCMC samples for \eqn{\rho} after burn-in. Used for plotting the posterior for \eqn{\rho} and performing inference of \eqn{\rho}.}
\item{gamma0.hat}{posterior median for \eqn{\gamma_0}. }
\item{gamma0.samples}{MCMC samples for \eqn{\gamma_0} after burn-in. Used for plotting the posterior for \eqn{\gamma_0} and performing inference of \eqn{\gamma_0}.}
\item{gamma1.hat}{posterior median for \eqn{\gamma_1}. }
\item{gamma1.samples}{MCMC samples for \eqn{\gamma_1} after burn-in. Used for plotting the posterior for \eqn{\gamma_1} and performing inference of \eqn{\gamma_1}.}
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

#####################
# Fit the RBC model #
#####################
# NOTE: Use default burn-in (burn=10000) and post-burn-in samples (nmc=10000)
RBC.mod = RobustBayesianCopas(y=y.obs, s=s.obs, burn=500, nmc=500)

# Point estimate for rho
rho.hat.RBC = RBC.mod$rho.hat

# Plot posterior for rho
hist(RBC.mod$rho.samples)

# Point estimate for theta
theta.hat.RBC = RBC.mod$theta.hat 
# Standard error for theta
theta.se.RBC = sd(RBC.mod$theta.samples)
# 95% posterior credible interval for theta
theta.cred.int = quantile(RBC.mod$theta.samples, probs=c(0.025,0.975))

# Display results
theta.hat.RBC
theta.se.RBC
theta.cred.int

# Plot the posterior for theta
hist(RBC.mod$theta.samples)

\donttest{
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

#####################
# Fit the RBC model #
######################
RBC.mod = RobustBayesianCopas(y=y.obs, s=s.obs)

# Point estimate for rho
rho.hat.RBC = RBC.mod$rho.hat

# Plot posterior for rho
hist(RBC.mod$rho.samples)

# Point estimate for theta
theta.hat.RBC = RBC.mod$theta.hat 
# 95% posterior credible interval for theta
theta.cred.int = quantile(RBC.mod$theta.samples, probs=c(0.025,0.975))
# Plot the posterior for theta
hist(RBC.mod$theta.samples)

}
}