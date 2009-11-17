\name{udt}
\docType{function}

\alias{udt}

\title{UNU.RAN object for Student t distribution}

\description{
  Create UNU.RAN object for a Student t distribution with
  with \code{df} degrees of freedom.
  
  [Distribution] -- t (Student).  %% Continuous
}

\usage{
udt(df, lb=-Inf, ub=Inf)
}

\arguments{
  \item{df}{degrees of freedom (strictly positive).
    Non-integer values allowed.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The \eqn{t} distribution with \code{df} \eqn{= \nu}{= n} degrees of
  freedom has density
  \deqn{
    f(x) = \frac{\Gamma ((\nu+1)/2)}{\sqrt{\pi \nu} \Gamma (\nu/2)}
    (1 + x^2/\nu)^{-(\nu+1)/2}
  }{f(x) = Gamma((n+1)/2) / (sqrt(n pi) Gamma(n/2)) (1 + x^2/n)^-((n+1)/2)}
  for all real \eqn{x}.
  It has mean \eqn{0} (for \eqn{\nu > 1}{n > 1}) and
  variance \eqn{\frac{\nu}{\nu-2}}{n/(n-2)} (for \eqn{\nu > 2}{n > 2}).

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
  N.L. Johnson, S. Kotz, and N. Balakrishnan (1995):
  Continuous Univariate Distributions, Volume 2.
  2nd edition, John Wiley & Sons, Inc., New York.
  Chap. 28, p. 362.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for t distribution
distr <- udt(df=4)
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
