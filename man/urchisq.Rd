\name{urchisq}
\docType{function}

\alias{urchisq}

\title{UNU.RAN Chi-Squared random variate generator}

\description{
  UNU.RAN random variate generator for the Chi-Squared
  (\eqn{\chi^2}{chi^2}) distribution with \code{df} degrees of freedom.
  It also allows sampling from the truncated distribution.
}

\usage{
urchisq(n, df, lb=0, ub=Inf)
}
\arguments{
  \item{n}{size of required sample.}
  \item{df}{degrees of freedom (strictly positive, but can be non-integer).}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Chi-squared distribution with \code{df}\eqn{= n > 0} degrees of
  freedom has density
  \deqn{f_n(x) = \frac{1}{{2}^{n/2} \Gamma (n/2)} {x}^{n/2-1} {e}^{-x/2}}{%
    f_n(x) = 1 / (2^(n/2) Gamma(n/2))  x^(n/2-1) e^(-x/2)}
  for \eqn{x > 0}.

  The generation algorithm uses fast numerical inversion. The parameters
  \code{lb} and \code{ub} can be used to generate variates from 
  the Chi-squared distribution truncated to the interval (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation, \code{\linkS4class{unuran}} for the UNU.RAN class, and
  \code{\link{rchisq}} for the \R built-in generator.
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
  Compared to \code{rchisq}, \code{urchisq} is faster, especially for
  larger sample sizes.
  However, in opposition to \code{rchisq} vector arguments are ignored,
  i.e. only the first entry is used.
}

\examples{
## Create a sample of size 1000
x <- urchisq(n=1000,df=3)
}

\keyword{distribution}
\keyword{datagen}

