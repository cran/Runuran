\name{urburr}
\docType{function}

\alias{urburr}

\title{UNU.RAN Burr random variate generator}

\description{
  UNU.RAN random variate generator for the Burr distribution with
  \code{shape1} and \code{shape2}.
  It also allows sampling from the truncated distribution.
  
  [Distribution] -- Sampling Function: Burr.  %% Continuous
}

\usage{
urburr(n, a, b, lb=0, ub=Inf)
}
\arguments{
  \item{n}{size of required sample.}
  \item{a,b}{positive shape parameters of the Burr distribution.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
 The Burr distribution with parameters \code{a} and \code{b} has density
  \deqn{
    f(x) = a (b-1) x^{a-1}/(1+x^a)^b}{%
    f(x) = a * (b-1) * x^(a-1)/(1+x^a)^b}
  for \eqn{x \ge 0}{x >= 0}, \eqn{a \ge 1}{a >= 1} and \eqn{b \ge 2}{b >= 2}.

  The generation algorithm uses transformed density rejection "TDR". The
  parameters \code{lb} and \code{ub} can be used to generate variates from 
  the Burr distribution truncated to the interval (\code{lb},\code{ub}).
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
  This function is a wrapper for the UNU.RAN class in \R.
}

\examples{
## Create a sample of size 1000
x <- urburr(n=1000,a=2,b=3)
}

\keyword{distribution}
\keyword{datagen}

