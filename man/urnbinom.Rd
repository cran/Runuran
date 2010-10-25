\name{urnbinom}
\alias{urnbinom}

\title{UNU.RAN Negative Binomial random variate generator}

\description{
  UNU.RAN random variate generator for the Negative Binomial
  distribution with with parameters \code{size} and \code{prob}.
  It also allows sampling from the truncated distribution.
  
  [Special Generator] -- Sampling Function: Negative Binomial.  %% Discrete
}

\usage{
urnbinom(n, size, prob, lb = 0, ub = Inf)
}

\arguments{
  \item{n}{size of required sample.}
  \item{size}{target for number of successful trials, or dispersion
    parameter (the shape parameter of the gamma mixing distribution).
    Must be strictly positive.}
  \item{prob}{probability of success in each trial. \code{0 < prob <= 1}.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Negative Binomial distribution with \code{size} \eqn{= n} and
  \code{prob} \eqn{= p} has density
  \deqn{
    p(x) = \frac{\Gamma(x+n)}{\Gamma(n) x!} p^n (1-p)^x
  }{p(x) = Gamma(x+n)/(Gamma(n) x!) p^n (1-p)^x}
  for \eqn{x = 0, 1, 2, \ldots}, \eqn{n > 0} and \eqn{0 < p \le 1}.
  This represents the number of failures which occur in a sequence of
  Bernoulli trials before a target number of successes is reached.

  The generation algorithm uses guide table based inversion when the
  tails are not too heavy and method \sQuote{DARI} otherwise.
  The parameters \code{lb} and \code{ub} can be used to generate
  variates from the Negative Binomial distribution truncated to the
  interval (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation, \code{\linkS4class{unuran}} for the UNU.RAN class, and
  \code{\link{rnbinom}} for the \R built-in generator.
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
  Compared to \code{rnbinom}, \code{urnbinom} is faster, especially for
  larger sample sizes.
  However, in opposition to \code{rnbinom} vector arguments are ignored,
  i.e. only the first entry is used.
}

\examples{
## Create a sample of size 1000
x <- urnbinom(n=1000,size=10,prob=0.3)
}

\keyword{distribution}
