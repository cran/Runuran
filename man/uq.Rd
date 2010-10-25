\name{uq}
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
  for continuous distributions and
  \itemize{
    \item \sQuote{DGT}
  }  
  for discrete distributions.
  
  \code{uq} returns the left boundary of the domain of the distribution
  if argument \code{U} is less than or equal to \code{0} and
  the right boundary if \code{U} is greater than or equal to \code{1}.
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
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Compute quantiles of normal distribution using method 'PINV'
gen <- pinv.new(pdf=dnorm, lb=-Inf, ub=Inf)
uq(gen,seq(0,1,0.05))

## Compute quantiles of user-defined distribution using method 'PINV'
pdf <- function (x) { exp(-x) }
gen <- pinv.new(pdf=pdf, lb=0, ub=Inf, uresolution=1.e-12)
uq(gen,seq(0,1,0.05))

## Compute quantiles of binomial distribution using method 'DGT'
gen <- dgt.new(pv=dbinom(0:1000,1000,0.4), from=0)
uq(gen,seq(0,1,0.05))


## Compute quantiles of normal distribution using method 'HINV'
## (using 'advanced' interface)
gen <- unuran.new("normal()","hinv")
uq(gen,0.975)
uq(gen,c(0.025,0.975))

## Compute quantiles of user-defined distributio using method 'HINV'
## (using 'advanced' interface)
cdf <- function (x) { 1.-exp(-x) }
pdf <- function (x) { exp(-x) }
dist <- new("unuran.cont", cdf=cdf, pdf=pdf, lb=0, ub=Inf)
gen <- unuran.new(dist, "hinv; u_resolution=1.e-12")
uq(gen,seq(0,1,0.05))

}

\keyword{distribution}
\keyword{datagen}
