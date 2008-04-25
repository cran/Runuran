\name{uqhinv}
\docType{function}

\alias{uqhinv}

\title{Quantile function using UNU.RAN method HINV}

\description{
  Evaluates quantile of distribution approximately using
  UNU.RAN method HINV.
}

\usage{
uqhinv(unur, U)
}
\arguments{
  \item{unur}{a "unuran" object using method \dQuote{HINV}.}
  \item{U}{vector of probabilities.}
}

\details{
  This routine evaluates the quantiles (inverser CDF) for a given
  (vector of) probabilities approximately.
  It requires a "unuran" object that uses method \dQuote{HINV} and thus
  contains the information about the target distribution.
}

\seealso{
  \code{\linkS4class{unuran}}.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

\examples{
## Compute quantiles of normal distribution
unr <- unuran.new("normal()","hinv")
uqhinv(unr,0.975)
uqhinv(unr,c(0.025,0.975))

## Compute quantiles of user-defined distribution
cdf <- function (x) { 1.-exp(-x) }
pdf <- function (x) { exp(-x) }
dist <- new("unuran.cont", cdf=cdf, pdf=pdf, lb=0, islog=FALSE)
unr <- unuran.new(dist, "hinv; u_resolution=1.e-12")
uqhinv(unr,seq(0,1,0.05))

}

\keyword{distribution}
\keyword{datagen}

