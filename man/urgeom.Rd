\name{urgeom}
\docType{function}

\alias{urgeom}

\title{UNU.RAN Geometric random variate generator}

\description{
  UNU.RAN random variate generator for the Geometric distribution
  with parameter \code{prob}.
  It also allows sampling from the truncated distribution.
  
  [Distribution] -- Sampling Function: Geometric.  %% Discrete
}

\usage{
urgeom(n, prob, lb = 0, ub = Inf)
}
\arguments{
  \item{n}{size of required sample.}
  \item{prob}{probability of success in each trial. \code{0 < prob <= 1}.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Geometric distribution with \code{prob} \eqn{= p} has density
  \deqn{p(x) = p {(1-p)}^{x}}{p(x) = p (1-p)^x}
  for \eqn{x = 0, 1, 2, \ldots}{x = 0, 1, 2, ...}, \eqn{0 < p \le 1}.

  The generation algorithm uses guide table based inversion for
  \eqn{p > 0.02} and method @code{"DARI"} otherwise.
  The parameters \code{lb} and \code{ub} can be used to generate
  variates from the Geometric distribution truncated to the interval
  (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation, \code{\linkS4class{unuran}} for the UNU.RAN class, and
  \code{\link{rgeom}} for the \R built-in generator.
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
  This function is a wrapper for the UNU.RAN class in \R.
  Compared to \code{rgeom}, \code{urgeom} is faster, especially for
  larger sample sizes.
  However, in opposition to \code{rgeom} vector arguments are ignored,
  i.e. only the first entry is used.
}

\examples{
## Create a sample of size 1000
x <- urgeom(n=1000,prob=0.2)
}

\keyword{distribution}
\keyword{datagen}

