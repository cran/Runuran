\name{urpareto}
\alias{urpareto}

\title{UNU.RAN Pareto random variate generator}

\description{
  UNU.RAN random variate generator for the Pareto distribution (of first
  kind) with shape parameters \code{k} and \code{a}.
  It also allows sampling from the truncated distribution.
  
  [Special Generator] -- Sampling Function: Pareto (of first kind).  %% Continuous
}

\usage{
urpareto(n, k, a, lb=k, ub=Inf)
}

\arguments{
  \item{n}{size of required sample.}
  \item{k}{(strictly positive) shape and location parameter.}
  \item{a}{(strictly positive) shape parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Pareto distribution with parameters \code{k} and \code{a} 
  has density
  \deqn{
    f(x) = a k^a x^{-(a+1)}
  }{f(x) = a * k^a * x^(-a-1)}
  for \eqn{x \ge k}, \eqn{k> 0} and \eqn{a>0}.

  The generation algorithm uses fast numerical inversion. The parameters
  \code{lb} and \code{ub} can be used to generate variates from 
  the Pareto distribution truncated to the interval (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation and \code{\linkS4class{unuran}} for the UNU.RAN class.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg

  N.L. Johnson, S. Kotz, and N. Balakrishnan (1994):
  Continuous Univariate Distributions,
  Volume 1. 2nd edition, 
  John Wiley & Sons, Inc., New York.
  Chap.20, p.574.  
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\note{
  This function is wrapper for the UNU.RAN class in \R.
}

\examples{
## Create a sample of size 1000
x <- urpareto(n=1000,k=2,a=3)
}

\keyword{distribution}
