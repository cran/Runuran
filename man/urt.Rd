\name{urt}
\docType{function}

\alias{urt}

\title{UNU.RAN Student t random variate generator}

\description{
  UNU.RAN random variate generator for the Student t distribution with
  with \code{df} degrees of freedom.
  It also allows sampling from the truncated distribution.
}

\usage{
urt(n, df, lb = -Inf, ub = Inf)
}
\arguments{
  \item{n}{size of required sample.}
  \item{df}{degrees of freedom (\eqn{> 0}, maybe non-integer).}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The \eqn{t} distribution with \code{df} \eqn{= \nu}{= n} degrees of
  freedom has density
  \deqn{
    f(x) = \frac{\Gamma ((\nu+1)/2)}{\sqrt{\pi \nu} \Gamma (\nu/2)}
    (1 + x^2/\nu)^{-(\nu+1)/2}%
  }{f(x) = Gamma((n+1)/2) / (sqrt(n pi) Gamma(n/2)) (1 + x^2/n)^-((n+1)/2)}
  for all real \eqn{x}.
  It has mean \eqn{0} (for \eqn{\nu > 1}{n > 1}) and
  variance \eqn{\frac{\nu}{\nu-2}}{n/(n-2)} (for \eqn{\nu > 2}{n > 2}).

  The generation algorithm uses fast numerical inversion. The parameters
  \code{lb} and \code{ub} can be used to generate variates from 
  the \eqn{t} distribution truncated to the interval (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation, \code{\linkS4class{unuran}} for the UNU.RAN class, and
  \code{\link{rt}} for the \R built-in generator.
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
  Compared to \code{rt}, \code{urt} is faster, especially for
  larger sample sizes.
  However, in opposition to \code{rt} vector arguments are ignored,
  i.e. only the first entry is used.
}

\examples{
## Create a sample of size 1000
x <- urt(n=1000,df=4)
}

\keyword{distribution}
\keyword{datagen}

