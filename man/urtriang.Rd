\name{urtriang}
\docType{function}

\alias{urtriang}

\title{UNU.RAN Triangular random variate generator}

\description{
  UNU.RAN random variate generator for the Triangular distribution with
  shape parameters \code{a}, \code{m} and \code{b}.
  It also allows sampling from the truncated distribution.
}

\usage{
urtriang(n, a, m, b, lb=a, ub=b)
}
\arguments{
  \item{n}{size of required sample.}
  \item{a,b}{left and right boundary of domain}
  \item{m}{mode of distribution}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Triangular distribution with domain \eqn{(a,b)} and mode \eqn{m}
  has density proportional to 
  \deqn{
    f(x) \sim (x-a)/(m-a)}{%
    f(x) ~ (x-a)/(m-a)}
  for \eqn{a \le x \le m}{a <= x <= m}, and
  \deqn{
    f(x) \sim (b-x)/(b-m)}{%
    f(x) ~ (b-x)/(b-m)}
  for \eqn{m \le x \le b}{m <= x <= b}.
  
  The generation algorithm uses fast numerical inversion. The parameters
  \code{lb} and \code{ub} can be used to generate variates from 
  the Triangular distribution truncated to the interval (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation and \code{\linkS4class{unuran}} for the UNU.RAN class.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

\note{
  This function is wrapper for the UNU.RAN class in \R.
}

\examples{
## Create a sample of size 1000
x <- urtriang(n=1000,a=-1,m=0,b=2)
}

\keyword{distribution}
\keyword{datagen}

