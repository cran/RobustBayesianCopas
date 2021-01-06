\name{D.measure}
\alias{D.measure}

\title{D Measure for Quantifying Publication Bias}

\description{This function computes Bai's \eqn{D} measure for quantifying publication bias based on the robust Bayesian Copas (RBC) selection model. Let \eqn{\pi_{rbc}(\theta | y)} be the posterior distribution for \eqn{\theta} under the full RBC (bias-corrected) model, and let \eqn{\pi_{\rho=0} (\theta | y)} be the posterior distribution for \eqn{\theta} under the \emph{non}-bias corrected model (when \eqn{\rho} is fixed at \eqn{\rho=0}). The \eqn{D} measure is the Hellinger distance \eqn{H} between the bias-corrected and \emph{non}-bias-corrected posteriors.

\deqn{D = H(\pi_{rbc}(\theta | y), \pi_{\rho=0} (\theta | y)).} 
\eqn{D} is always between 0 and 1, with \eqn{D \approx 0} indicating negligible publication bias and \eqn{D \approx 1} indicating a very high magnitude of publication bias.

The random effects distributions for the bias-corrected and \emph{non}-bias-corrected models should be the same in order for the \eqn{D} measure to be valid. The posterior densities for \eqn{\pi_{rbc}(\theta | y)} and \eqn{\pi_{\rho=0} (\theta | y)} are approximated using MCMC samples. Numerical integration is used to compute the Hellinger distance between them.
}

\usage{
D.measure(samples.RBCmodel, samples.nobiasmodel)
}

\arguments{
  \item{samples.RBCmodel}{A vector of the MCMC samples from the RBC model. These can be obtained from the output of the function \code{RobustBayesianCopas}.}
  \item{samples.nobiasmodel}{A vector of the MCMC samples from the non-bias-corrected model. These can be obtained from the output of the function \code{BayesNonBiasCorrected}. In order for the \emph{D} measure to be valid, the random effects distribution used in \code{BayesNonBiasCorrected} should be the \emph{same} as the random effects distribution used in \code{RobustBayesianCopas}.}
}

\value{The function returns Bai's \eqn{D} measure, a value between 0 and 1. \eqn{D \approx 0} means negligible publication bias, and \eqn{D \approx 1} means a very high magnitude of publication bias.}

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

#############################################
# Compute D measure using posterior samples #
#############################################
# Fit RBC (bias-corrected) model with t-distributed random effects.
# NOTE: Use default burn-in (burn=10000) and post-burn-in samples (nmc=10000)
RBC.mod = RobustBayesianCopas(y=y.obs, s=s.obs, re.dist="StudentsT", burn=500, nmc=500)

# Fit non-bias-corrected model with t-distributed random effects.
# NOTE: Use default burn-in (burn=10000) and post-burn-in samples (nmc=10000)
RBCNoBias.mod = BayesNonBiasCorrected(y=y.obs, s=s.obs, re.dist="StudentsT", burn=500, nmc=500)

# Compute the D measure based on posterior samples for theta
D = D.measure(RBC.mod$theta.samples, RBCNoBias.mod$theta.samples)

\donttest{


############################################
# Example on the antidepressants data set. #
# This is from Section 6.2 of the paper    #
# by Bai et al. (2020).                    #
############################################

# Set seed, so we can reproduce the exact same results as in the paper.
seed = 123
set.seed(seed)

# Load the full data
data(antidepressants)
attach(antidepressants)

# Extract the 50 published studies
published.data = antidepressants[which(antidepressants$Published==1),]

# Observed treatment effect
y.obs = published.data$Standardized_effect_size

# Observed standard error
s.obs = published.data$Standardized_SE


################################
# Compute the D measure for    #
# quantifying publication bias #
################################

# Fit RBC model with different random effects distributions and compare them using DIC.

# Normal
RBC.mod.normal = RobustBayesianCopas(y=y.obs, s=s.obs, re.dist="normal", seed=seed) 
RBC.mod.normal$DIC # DIC=369.3126
# Student's t
RBC.mod.StudentsT = RobustBayesianCopas(y=y.obs, s=s.obs, re.dist="StudentsT", seed=seed)
RBC.mod.StudentsT$DIC # DIC=515.5151
# Laplace
RBC.mod.laplace = RobustBayesianCopas(y=y.obs, s=s.obs, re.dist="Laplace", seed=seed)
RBC.mod.laplace$DIC # DIC=475.5845
# Slash
RBC.mod.slash = RobustBayesianCopas(y=y.obs, s=s.obs, re.dist="slash", seed=seed)
RBC.mod.slash$DIC # DIC=381.2705

# The model with normal random effects gave the lowest DIC, so we use the model
# with normal study-specific effects for our final analysis.

RBC.mod.normal$rho.hat # rho.hat=0.804
# Plot posterior for rho, which suggests strong publication bias.
hist(RBC.mod.normal$rho.samples) 

# Fit non-biased-corrected model. Make sure that we specify the random effects as normal.
RBCNoBias.mod = BayesNonBiasCorrected(y=y.obs, s=s.obs, re.dist="normal", seed=seed)

# Compute D measure using posterior samples for theta
D.RBC = D.measure(RBC.mod.normal$theta.samples, RBCNoBias.mod$theta.samples) # D=0.95
}
}