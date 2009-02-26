\name{dari.new}
\docType{function}

\alias{dari.new}

\title{UNU.RAN generator based on Discrete Automatic Rejection Inversion (DARI)}

\description{
  UNU.RAN random variate generator for continuous distributions with
  given probability mass function (PMF).
  It is based on Discrete Automatic Rejection Inversion (\sQuote{DARI}).

  [Universal] -- Rejection Method.
}

\usage{
dari.new(pmf, lb, ub, mode=NA, sum=1, ...)
}
\arguments{
  \item{pmf}{probability mass function. (\R function)}
  \item{lb}{lower bound of domain;
    use \code{-Inf} if unbounded from left. (numeric, integer)}
  \item{ub}{upper bound of domain;
    use \code{Inf} if unbounded from right. (numeric, integer)}
  \item{mode}{mode of distribution. (integer)}
  \item{sum}{sum over all \dQuote{probabilities}. (numeric)}
  \item{...}{(optional) arguments for \code{pmf}}
}

\details{
  This function creates an \code{unuran} object based on \sQuote{DARI}
  (Discrete Automatic Rejection Inversion). It can be used to draw
  samples of a discrete random variate with given probability mass function
  using \code{\link{ur}}.
  
  The probabilities must be provided by a probability mass function
  \code{pmf} which must return non-negative numbers and which need not
  be normalized (i.e., it can be any multiple of a PMF).
  Moreover, the given function must be \eqn{T_{-0.5}}{T_(-0.5)}-concave;
  this includes all log-concave distributions.

  In addition the algorithm requires the location of the \code{mode}.
  If omitted then it is computed by a slow numerical search.

  If the sum over all probabilities is different from 1 then a rough
  estimate of this sum is required.
}

\seealso{
  \code{\link{ur}}, \code{\link{unuran.new}}, \code{\linkS4class{unuran}}.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg.
  See Section 10.2 (Tranformed Probability Rejection).

}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

\examples{
## Create a sample of size 100 for a Binomial distribution
## with 1000 number if observations and probability 0.2
gen <- dari.new(pmf=dbinom, lb=0, ub=1000, size=1000, prob=0.2)
x <- ur(gen,100)

## Create a sample from a distribution with PMF
##  p(x) = 1/x^3, x >= 1  (Zipf distribution)
zipf <- function (x) { 1/x^3 }
gen <- dari.new(pmf=zipf, lb=1, ub=Inf)
x <- ur(gen,100)

}

\keyword{distribution}
\keyword{datagen}

