\name{udchi}
\docType{function}

\alias{udchi}

\title{UNU.RAN object for Chi distribution}

\description{
  Create UNU.RAN object for a Chi distribution
  with \code{df} degrees of freedom.
  
  [Distribution] -- Chi.  %% Continuous
}

\usage{
udchi(df, lb=0, ub=Inf)
}

\arguments{
  \item{df}{degrees of freedom (strictly positive).
    Non-integer values allowed.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Chi distribution with \code{df}\eqn{= n > 0} degrees of
  freedom has density
  \deqn{
    f(x) = x^{n-1} e^{-x^2/2}
  }{f(x) = x^(n-1) * exp(-x^2/2)}
  for \eqn{x > 0}.
  
  The domain of the distribution can be truncated to the 
  interval (\code{lb},\code{ub}).
}

\value{
  An object of class \code{"unuran.cont"}.
}

\seealso{
  \code{\linkS4class{unuran.cont}}.
}

\references{
  N.L. Johnson, S. Kotz, and N. Balakrishnan (1994):
  Continuous Univariate Distributions, Volume 1.
  2nd edition, John Wiley & Sons, Inc., New York.
  Chap. 18, p. 417.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for chi-squared distribution
distr <- udchi(df=5)
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
