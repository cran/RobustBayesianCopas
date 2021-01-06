\name{antidepressants}
\alias{antidepressants}

\title{A Meta-Analysis on the Efficacy of Antidepressants}
\description{This data set contains 73 studies with results on the effectiveness of antidepressants that were reported to the FDA. However, only 50 of these studies were subsequently published. Since studies reported their outcomes on different scales, effect sizes were all expressed as standardized mean differences by means of Hedges' \emph{g} scores, accompanied by corresponding variances. This data set was originally analyzed by Turner et al. (2008).
 }
\usage{data(antidepressants)}
\format{
A dataframe with 73 studies with the following seven variables. 
\describe{
    \item{\code{Drug}:}{Antidepressant name.}
    \item{\code{Study}:}{Study identifier.}
    \item{\code{Published}:}{A binary variable indicating whether the study was published: "1"=published, "0"=not published.}
    \item{\code{Nonstandardized_effect_size}:}{Estimated mean improvement in depression symptoms (nonstandardized).}
    \item{\code{Nonstandardized_SE}:}{Estimated standard error (nonstandardized).}
    \item{\code{Standardized_effect_size}}{Estimated mean improvement in depression symptoms (standardized). Note that the standardized values are only available for the \emph{published} studies (\code{NA} if not published). The \emph{non}-missing data in this column should be used as \code{y} in the selection model.}
        \item{\code{Standardized_SE}:}{Estimated standard error (standardized). Note that the standardized values are only available for the \emph{published} studies (\code{NA} if not published). The \emph{non}-missing data in this column should be used as \code{s} in the selection model. }
  }
}

\source{ 
Turner, E. H., Matthews, A. M., Linardatos, E., Tell, R. A., and Rosenthal, R. (2008). "Selective publication of antidepressant trials and its influence on apparent efficacy." \emph{New England Journal of Medicine}, \bold{358}(3):252-260.
}