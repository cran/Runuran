\name{urbinom}
\docType{function}

\alias{urbinom}

\title{UNU.RAN Binomial random variate generator}

\description{
  UNU.RAN random variate generator for the Binomial distribution
  with parameters \code{size}and \code{prob}.
  It also allows sampling from the truncated distribution.
}

\usage{
urbinom(n, size, prob, lb = 0, ub = size)
}
\arguments{
  \item{n}{size of required sample.}
  \item{size}{number of trials (one or more).}
  \item{prob}{probability of success on each trial.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Binomial distribution with \code{size} \eqn{= n} and
  \code{prob} \eqn{= p} has density
  \deqn{p(x) = {n \choose x} {p}^{x} {(1-p)}^{n-x}}{
    p(x) = choose(n,x) p^x (1-p)^(n-x)}
  for \eqn{x = 0, \ldots, n}.

  The generation algorithm uses guide table based inversion. The parameters
  \code{lb} and \code{ub} can be used to generate variates from 
  the Binomial distribution truncated to the interval (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation, \code{\link{unuran-class}} for the UNU.RAN class, and
  \code{\link{rbinom}} for the \R built-in generator.
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
  Compared to \code{rbinom}, \code{urbinom} is faster, especially for
  larger sample sizes.
  However, in opposition to \code{rbinom} vector arguments are ignored,
  i.e. only the first entry is used.
}

\examples{
## Create a sample of size 1000
x <- urbinom(n=1000,size=10,prob=0.3)
}

\keyword{distribution}
\keyword{datagen}

