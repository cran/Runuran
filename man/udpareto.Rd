\name{udpareto}
\alias{udpareto}

\title{UNU.RAN object for Pareto distribution}

\description{
  Create UNU.RAN object for a Pareto distribution (of first
  kind) with shape parameters \code{k} and \code{a}.

  [Distribution] -- Pareto (of first kind).  %% Continuous
}

\usage{
udpareto(k, a, lb=k, ub=Inf)
}

\arguments{
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

  The domain of the distribution can be truncated to the 
  interval (\code{lb},\code{ub}).
}

\value{
  An object of class \code{"unuran.cont"}.
}

\seealso{
  \code{\linkS4class{unuran.cont}}.
}

\references{
  N.L. Johnson, S. Kotz, and N. Balakrishnan (1994):
  Continuous Univariate Distributions, Volume 1.
  2nd edition, John Wiley & Sons, Inc., New York.
  Chap. 20, p. 574.  
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for Pareto distribution
distr <- udpareto(k=3,a=2)
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
