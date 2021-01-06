\name{RobustBayesianCopas}
\alias{RobustBayesianCopas}

\title{Robust Bayesian Copas selection model}

\description{This function implements the Robust Bayesian Copas selection model of Bai et al. (2020) for the Copas selection model,

\deqn{y_i | (z_i>0) = \theta + \tau u_i + s_i \epsilon_i,}
\deqn{z_i = \gamma_0 + \gamma_1 / s_i + \delta_i,}
\deqn{corr(\epsilon_i, \delta_i) = \rho,}

where \eqn{y_i} is the reported treatment effect for the \eqn{i}th study, \eqn{s_i} is the reported standard error for the \eqn{i}th study, \eqn{\theta} is the population treatment effect of interest, \eqn{\tau > 0} is a heterogeneity parameter, and \eqn{\epsilon_i}, and \eqn{\delta_i} are marginally distributed as \eqn{N(0,1)} and \eqn{u_i}, and \eqn{\epsilon_i} are independent. 

In the Copas selection model, \eqn{y_i} is published (selected) if and only if the corresponding propensity score \eqn{z_i} (or the propensity to publish) is greater than zero. The propensity score \eqn{z_i} contains two parameters: \eqn{\gamma_0} controls the overall probability of publication, and \eqn{\gamma_1} controls how the probability of publication depends on study sample size. The reported treatment effects and propensity scores are correlated through \eqn{\rho}. If \eqn{\rho=0}, then there is no publication bias and the Copas selection model reduces to the standard random effects meta-analysis model.  

The RBC model places noninformative priors on \eqn{(\theta, \tau^2, \rho, \gamma_0, \gamma_1)} (see Bai et al. (2020) for details). For the random effects \eqn{u_i, i=1, \ldots, n}, we give the option for using normal, Student's t, Laplace, or slash distributions for the random effects. The function returns the Deviance Information Criterion (DIC), which can be used to select the appropriate distribution to use for the final analysis. 
}

\usage{
RobustBayesianCopas(y, s, re.dist=c("normal", "StudentsT", "Laplace", "slash"),
                    t.df = 4, slash.shape = 1, init=NULL, seed=NULL, 
                    burn=10000, nmc=10000)}

\arguments{
  \item{y}{An \eqn{n \times 1} vector of reported treatment effects.}
  \item{s}{An \eqn{n \times 1} vector of reported within-study standard errors.}
    \item{re.dist}{Distribution for the between-study random effects \eqn{u_i, i=1, \ldots, n}. The user may specify the normal, Student's T, Laplace, or slash distribution. The default is \code{StudentsT} with 4 degrees of freedom.}
        \item{t.df}{Degrees of freedom for t-distribution. Only used if \code{StudentsT} is specified for \code{re.dist}. Default is 4.}
    \item{slash.shape}{Shape parameter in the slash distribution. Only used if \code{slash} is specified for \code{re.dist}. Default is 1.}
  \item{init}{Optional initialization values for \eqn{(\theta, \tau, \rho, \gamma_0, \gamma_1)}. If specified, they must be provided in this exact order. If they are not provided, the program estimates initial values from the data.}
    \item{seed}{Optional seed. This needs to be specified if you want to reproduce the exact results of your analysis. }
  \item{burn}{Number of burn-in samples. Default is \code{burn=10000}.}
  \item{nmc}{Number of posterior samples to save. Default is \code{nmc=10000}.}
}

\value{The function returns a list containing the following components:
\item{DIC}{Deviance Information Criterion (DIC), a measure of model fit. This can be used to compare the results for different random effects distributions. The model that gives the lowest DIC gives the best fit to the data.}
\item{theta.hat}{Posterior mean for \eqn{\theta}. }
\item{theta.samples}{MCMC samples for \eqn{\theta} after burn-in. Used for plotting the posterior for \eqn{\theta} and performing inference of \eqn{\theta}.}
\item{tau.hat}{Posterior mean for \eqn{\tau}. }
\item{tau.samples}{MCMC samples for \eqn{\tau} after burn-in. Used for plotting the posterior for \eqn{\tau} and performing inference of \eqn{\tau}.}
\item{rho.hat}{Posterior median for \eqn{\rho}. }
\item{rho.samples}{MCMC samples for \eqn{\rho} after burn-in. Used for plotting the posterior for \eqn{\rho} and performing inference of \eqn{\rho}.}
\item{gamma0.hat}{Posterior median for \eqn{\gamma_0}. }
\item{gamma0.samples}{MCMC samples for \eqn{\gamma_0} after burn-in. Used for plotting the posterior for \eqn{\gamma_0} and performing inference of \eqn{\gamma_0}.}
\item{gamma1.hat}{Posterior median for \eqn{\gamma_1}. }
\item{gamma1.samples}{MCMC samples for \eqn{\gamma_1} after burn-in. Used for plotting the posterior for \eqn{\gamma_1} and performing inference of \eqn{\gamma_1}.}
}

\references{
Bai, R., Lin, L., Boland, M. R., and Chen, Y. (2020). "A robust Bayesian Copas selection model for quantifying and correcting publication bias." \emph{arXiv preprint arXiv:2005.02930.}
}

\examples{
######################################
# Example on the Barlow2014 data set #
######################################
# Load data
data(Barlow2014)
attach(Barlow2014)

# Observed treatment effect
y.obs = Barlow2014[,1]
# Observed standard error
s.obs = Barlow2014[,2]

###############################################
# Fit the RBC model with slash random effects #
###############################################
# NOTE: Use default burn-in (burn=10000) and post-burn-in samples (nmc=10000)
# Fit model with slash errors
RBC.mod = RobustBayesianCopas(y=y.obs, s=s.obs, re.dist="slash", burn=500, nmc=500)

# Point estimate for rho
rho.hat.RBC = RBC.mod$rho.hat
rho.hat.RBC
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
# Example on second-hand smoking data set. #
# This is from Section 6.1 of the paper by #
# Bai et al. (2020).                       #
############################################

# Set seed, so we can reproduce the exact same result as in the paper.
seed = 1234
set.seed(seed)

# Load the full data
data(Hackshaw1997)
attach(Hackshaw1997)
# Extract the log OR
y.obs = Hackshaw1997[,2]
# Extract the observed standard error
s.obs = Hackshaw1997[,3]


###################################################
# Fit the RBC model with different random effects #
# distributions and compare them using the DIC.   #
###################################################

# Normal
RBC.mod.normal = RobustBayesianCopas(y=y.obs, s=s.obs, re.dist="normal", seed=seed) 
RBC.mod.normal$DIC # DIC=429.7854
# Student's t
RBC.mod.StudentsT = RobustBayesianCopas(y=y.obs, s=s.obs, re.dist="StudentsT", seed=seed) 
RBC.mod.StudentsT$DIC # DIC=399.1955
# Laplace
RBC.mod.Laplace = RobustBayesianCopas(y=y.obs, s=s.obs, re.dist="Laplace", seed=seed) 
RBC.mod.Laplace$DIC # DIC=410.9086
# Slash
RBC.mod.slash = RobustBayesianCopas(y=y.obs, s=s.obs, re.dist="slash", seed=seed)
RBC.mod.slash$DIC # DIC=407.431


#######################################################
# Use the model with t-distributed random errors for  #
# the final analysis since it gave the lowest DIC.    #
#######################################################

# Point estimate for rho
rho.hat.RBC = RBC.mod.StudentsT$rho.hat # rho.hat=0.459 (moderate publication bias)
# Plot posterior for rho
hist(RBC.mod.StudentsT$rho.samples)

# Point estimate for theta
theta.hat.RBC = RBC.mod.StudentsT$theta.hat # theta.hat=0.1672
# 95% posterior credible interval for theta
theta.cred.int = quantile(RBC.mod.StudentsT$theta.samples, probs=c(0.025,0.975))
# Plot the posterior for theta
hist(RBC.mod.StudentsT$theta.samples)

# Obtain odds ratio estimates
OR.samples.RBC = exp(RBC.mod.StudentsT$theta.samples) # Samples of exp(theta)
# Posterior mean OR
OR.RBC.hat = mean(OR.samples.RBC) # OR.hat=1.185
# 95% posterior credible interval for OR
OR.RBC.credint = quantile(OR.samples.RBC, probs=c(0.025,0.975)) # (1.018, 1.350)


##############################################
# Use D measure to quantify publication bias #
##############################################

# Make sure that we specify the random effects as Student's t, since that is
# the distribution that we used for our final analysis.
Bayes.nobias.mod = BayesNonBiasCorrected(y=y.obs, s=s.obs, re.dist="StudentsT", seed=seed)

# Compute D measure based on posterior samples of theta
D.RBC = D.measure(RBC.mod.StudentsT$theta.samples, Bayes.nobias.mod$theta.samples)
D.RBC # D=0.33
}
}