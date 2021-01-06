\name{Hackshaw1997}
\alias{Hackshaw1997}

\title{A Meta-Analysis on the Relationship Between Second-hand Tobacco Smoke and Lung Cancer}
\description{This data set contains 37 studies analyzed by Hackshaw et al. (1997). Hackshaw et al. (1997) evaluated the risk of developing lung cancer in women who were lifelong nonsmokers but whose husbands smoked, compared to women whose husbands had never smoked.
 }
\usage{data(Hackshaw1997)}
\format{
A dataframe with 37 studies with the following four variables within each study. 
\describe{
    \item{\code{Study}:}{Study identifier.}
    \item{\code{log_OR}:}{The reported log-odds ratio. Use this as the treatment effect in meta-analysis.}
    \item{\code{SE}:}{The reported standard error.}
    \item{\code{weight}:}{The percent weight that the study contributes to the pooled log-odds ratio.}
  }
}

\source{ 
Hackshaw, A. K., Law, M. R., and Wald, N. J. (1997). "The accumulated evidence on lung cancer and environmental tobacco smoke." \emph{BMJ}, \bold{315}(7114):980-988.
}