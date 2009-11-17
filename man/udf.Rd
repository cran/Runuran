\name{udf}
\docType{function}

\alias{udf}

\title{UNU.RAN object for F distribution}

\description{
  Create UNU.RAN object for an F distribution with mean
  with \code{df1} and \code{df2} degrees of freedom.
  
  [Distribution] -- F.  %% Continuous
}

\usage{
udf(df1, df2, lb=0, ub=Inf)
}

\arguments{
  \item{df1, df2}{(strictly positive) degrees of freedom.
    Non-integer values allowed.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The F distribution with \code{df1 =} \eqn{n_1}{n1} and \code{df2 =}
  \eqn{n_2}{n2} degrees of freedom has density
  \deqn{
    f(x) = \frac{\Gamma(n_1/2 + n_2/2)}{\Gamma(n_1/2)\Gamma(n_2/2)}
    \left(\frac{n_1}{n_2}\right)^{n_1/2} x^{n_1/2 -1}
    \left(1 + \frac{n_1 x}{n_2}\right)^{-(n_1 + n_2) / 2}
  }{f(x) = Gamma((n1 + n2)/2) / (Gamma(n1/2) Gamma(n2/2)) * (n1/n2)^(n1/2) x^(n1/2 - 1) * (1 + (n1/n2) x)^(-(n1 + n2)/2)}
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
  N.L. Johnson, S. Kotz, and N. Balakrishnan (1995):
  Continuous Univariate Distributions, Volume 2.
  2nd edition, John Wiley & Sons, Inc., New York.
  Chap. 27, p. 332
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for F distribution
distr <- udf(df1=3,df2=6)
## Generate generator object; use method PINV (inversion)
gen <- pinvd.new(distr)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
