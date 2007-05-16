\name{urhyper}
\docType{function}

\alias{urhyper}

\title{UNU.RAN Hypergeometric random variate generator}

\description{
  UNU.RAN random variate generator for the Hypergeometric distribution.
  It also allows sampling from the truncated distribution.
}

\usage{
urhyper(nn, m, n, k, lb=max(0,k-n), ub=min(k,m))
}
\arguments{
  \item{nn}{number of observations.}
  \item{m}{the number of white balls in the urn.}
  \item{n}{the number of black balls in the urn.}
  \item{k}{the number of balls drawn from the urn.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Hypergeometric distribution is used for sampling \emph{without}
  replacement.  The density of this distribution with parameters
  \code{m}, \code{n} and \code{k} (named \eqn{Np}, \eqn{N-Np}, and
  \eqn{n}, respectively in the reference below) is given by
  \deqn{
    p(x) = \left. {m \choose x}{n \choose k-x} \right/ {m+n \choose k}%
  }{p(x) =      choose(m, x) choose(n, k-x) / choose(m+n, k)}
  for \eqn{x = 0, \ldots, k}{x = 0, ..., k}.

  The generation algorithm uses guide table based inversion. The parameters
  \code{lb} and \code{ub} can be used to generate variates from 
  the Hypergeometric distribution truncated to the interval (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation, \code{\link{unuran-class}} for the UNU.RAN class, and
  \code{\link{rhyper}} for the \R built-in generator.
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
  Compared to \code{rhyper}, \code{urhyper} is faster, especially for
  larger sample sizes.
  However, in opposition to \code{rhyper} vector arguments are ignored,
  i.e. only the first entry is used.
}

\examples{
## Create a sample of size 1000
x <- rhyper(nn=20,m=15,n=5,k=7)
}

\keyword{distribution}
\keyword{datagen}

