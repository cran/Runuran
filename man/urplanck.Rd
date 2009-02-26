\name{urplanck}
\docType{function}

\alias{urplanck}

\title{UNU.RAN Planck random variate generator}

\description{
  UNU.RAN random variate generator for the Planck distribution with
  shape parameter \code{a}.
  It also allows sampling from the truncated distribution.
  
  [Distribution] -- Sampling Function: Planck.  %% Continuous
}

\usage{
urplanck(n, a, lb = 1.e-12, ub = Inf)
}
\arguments{
  \item{n}{size of required sample.}
  \item{a}{(strictly positive) shape parameter.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Planck distribution with parameter \code{a}
  has density proportional to 
  \deqn{
    f(x) \sim \frac{x^a}{\exp(x)-1}}{%
    f(x) ~ (x^a)/(exp(x)-1)}
  for \eqn{x \ge 0} and \eqn{a\ge 1}{a >= 1}.

  The generation algorithm uses transformed density rejection
  \sQuote{TDR}. The parameters \code{lb} and \code{ub} can be used to
  generate variates from the Planck distribution truncated to the
  interval (\code{lb},\code{ub}).
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
x <- urplanck(n=1000,a=2)
}

\keyword{distribution}
\keyword{datagen}

