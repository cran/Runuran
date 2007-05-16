\name{urgamma}
\docType{function}

\alias{urgamma}

\title{UNU.RAN Gamma random variate generator}

\description{
  UNU.RAN random variate generator for the Gamma distribution with
  with parameters \code{shape} and \code{scale}.
  It also allows sampling from the truncated distribution.
}

\usage{
urgamma(n, shape, scale=1, lb=0, ub=Inf)
}
\arguments{
  \item{n}{size of required sample.}
  \item{shape}{(strictly positive) shape parameter.}
  \item{scale}{(strictly positive) scale parameter.}
  \item{lb}{lower bound of (truncated) distribution}
  \item{ub}{upper bound of (truncated) distribution}
}

\details{
  If \code{scale} is omitted, it assumes the default value of \code{1}.

  The Gamma distribution with parameters \code{shape} \eqn{=\alpha}{= a}
  and \code{scale} \eqn{=\sigma}{= s} has density
  \deqn{
    f(x)= \frac{1}{{\sigma}^{\alpha}\Gamma(\alpha)} {x}^{\alpha-1} e^{-x/\sigma}}{%
    f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s)}
  for \eqn{x \ge 0}, \eqn{\alpha > 0}{a > 0} and \eqn{\sigma > 0}{s > 0}.
  (Here \eqn{\Gamma(\alpha)}{Gamma(a)} is the function implemented by \R's
  \code{\link{gamma}()} and defined in its help.)

  The generation algorithm uses fast numerical inversion. The parameters
  \code{lb} and \code{ub} can be used to generate variates from 
  the Gamma distribution truncated to the interval (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation, \code{\link{unuran-class}} for the UNU.RAN class, and
  \code{\link{rgamma}} for the \R built-in generator.
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
  Compared to \code{rgamma}, \code{urgamma} is faster, especially for
  larger sample sizes.
  However, in opposition to \code{rgamma} vector arguments are ignored,
  i.e. only the first entry is used.
}

\examples{
## Create a sample of size 1000
x <- urgamma(n=1000,shape=2)
}

\keyword{distribution}
\keyword{datagen}

