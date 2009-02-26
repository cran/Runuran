\name{uq}
\docType{function}

\alias{uq}

\title{Quantile function for "unuran" object}

\description{
  Evaluates quantile of distribution approximately using
  a \code{unuran} object that implements an inversion method.
  
  [Universal] -- Quantile Function.
}

\usage{
uq(unr, U)
}
\arguments{
  \item{unr}{a \code{unuran} object that implements an inversion menthod.}
  \item{U}{vector of probabilities.}
}

\details{
  The routine evaluates the quantiles (inverse CDF) for a given
  (vector of) probabilities approximately.
  It requires a \code{unuran} object that implements an inversion method.
  Currently these are
  \itemize{
    \item \sQuote{HINV}
    \item \sQuote{NINV}
    \item \sQuote{PINV}
  }
}

\seealso{
  \code{\linkS4class{unuran}},\code{\link{unuran.new}}.
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
## Compute quantiles of normal distribution using method 'HINV'
unr <- unuran.new("normal()","hinv")
uq(unr,0.975)
uq(unr,c(0.025,0.975))

## Compute quantiles of user-defined distributio using method 'HINV'
cdf <- function (x) { 1.-exp(-x) }
pdf <- function (x) { exp(-x) }
dist <- new("unuran.cont", cdf=cdf, pdf=pdf, lb=0, ub=Inf)
unr <- unuran.new(dist, "hinv; u_resolution=1.e-12")
uq(unr,seq(0,1,0.05))

## Compute quantiles of user-defined distribution using method PINV
pdf <- function (x) { exp(-x) }
dist <- new("unuran.cont", pdf=pdf, lb=0, ub=Inf)
unr <- unuran.new(dist, "pinv; u_resolution=1.e-12")
uq(unr,seq(0,1,0.05))

}

\keyword{distribution}
\keyword{datagen}

