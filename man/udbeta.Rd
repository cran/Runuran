\name{udbeta}
\docType{function}

\alias{udbeta}

\title{UNU.RAN object for Beta distribution}

\description{
  Create UNU.RAN object for a Beta distribution with
  with parameters \code{shape1} and \code{shape2}.
  
  [Distribution] -- Beta.  %% Continuous
}

\usage{
udbeta(shape1, shape2, lb=0, ub=1)
}

\arguments{
  \item{shape1,shape2}{positive shape parameters of the Beta distribution.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Beta distribution with parameters \code{shape1} \eqn{= a} and
  \code{shape2} \eqn{= b} has density
  \deqn{
    f(x) =
    \frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)}{x}^{a} {(1-x)}^{b}
  }{f(x) = Gamma(a+b)/(Gamma(a)Gamma(b))x^(a-1)(1-x)^(b-1)}
  for \eqn{a > 0}, \eqn{b > 0} and \eqn{0 \le x \le 1}{0 <= x <= 1}.

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
  Chap. 25, p. 210.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for beta distribution
distr <- udbeta(shape1=3,shape2=7)
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}