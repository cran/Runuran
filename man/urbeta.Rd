\name{urbeta}
\alias{urbeta}

\title{UNU.RAN Beta random variate generator}

\description{
  UNU.RAN random variate generator for the Beta distribution with
  parameters \code{shape1} and \code{shape2}.
  It also allows sampling from the truncated distribution.
  
  [Special Generator] -- Sampling Function: Beta.  %% Continuous
}

\usage{
urbeta(n, shape1, shape2, lb = 0, ub = 1)
}

\arguments{
  \item{n}{size of required sample.}
  \item{shape1,shape2}{positive shape parameters of the Beta distribution.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Beta distribution with parameters \code{shape1} \eqn{= a} and
  \code{shape2} \eqn{= b} has density
  \deqn{
    f(x) = \frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)}{x}^{a} {(1-x)}^{b}
  }{f(x) = Gamma(a+b)/(Gamma(a)Gamma(b))x^(a-1)(1-x)^(b-1)}
  for \eqn{a > 0}, \eqn{b > 0} and \eqn{0 \le x \le 1}{0 <= x <= 1}.

  The generation algorithm uses fast numerical inversion. The parameters
  \code{lb} and \code{ub} can be used to generate variates from 
  the Beta distribution truncated to the interval (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation, \code{\linkS4class{unuran}} for the UNU.RAN class, and
  \code{\link{rbeta}} for the \R built-in generator.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\note{
  This function is wrapper for the UNU.RAN class in \R.
  Compared to \code{rbeta}, \code{urbeta} is faster, especially for
  larger sample sizes.
  However, in opposition to \code{rbeta} vector arguments are ignored,
  i.e. only the first entry is used.
}

\examples{
## Create a sample of size 1000
x <- urbeta(n=1000,shape1=2,shape2=5)
}

\keyword{distribution}
