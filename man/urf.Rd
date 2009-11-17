\name{urf}
\docType{function}

\alias{urf}

\title{UNU.RAN F random variate generator}

\description{
  UNU.RAN random variate generator for the F distribution with
  with \code{df1} and \code{df2} degrees of freedom.
  It also allows sampling from the truncated distribution.
  
  [Special Generator] -- Sampling Function: F.  %% Continuous
}

\usage{
urf(n, df1, df2, lb=0, ub=Inf)
}

\arguments{
  \item{n}{size of required sample.}
  \item{df1, df2}{(strictly positive) degrees of freedom. Non-integer
    values allowed.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The F distribution with \code{df1 =} \eqn{n_1}{n1} and \code{df2 =}
  \eqn{n_2}{n2} degrees of freedom has density
  \deqn{
    f(x) = \frac{\Gamma(n_1/2 + n_2/2)}{\Gamma(n_1/2)\Gamma(n_2/2)}
    \left(\frac{n_1}{n_2}\right)^{n_1/2} x^{n_1/2 -1}
    \left(1 + \frac{n_1 x}{n_2}\right)^{-(n_1 + n_2) / 2}%
  }{f(x) = Gamma((n1 + n2)/2) / (Gamma(n1/2) Gamma(n2/2)) * (n1/n2)^(n1/2) x^(n1/2 - 1) * (1 + (n1/n2) x)^-(n1 + n2)/2}
  for \eqn{x > 0}.

  The generation algorithm uses fast numerical inversion. The parameters
  \code{lb} and \code{ub} can be used to generate variates from 
  the F distribution truncated to the interval (\code{lb},\code{ub}).
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation, \code{\linkS4class{unuran}} for the UNU.RAN class, and
  \code{\link{rf}} for the \R built-in generator.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\note{
  This function is wrapper for the UNU.RAN class in \R.
  Compared to \code{rf}, \code{urf} is faster, especially for
  larger sample sizes.
  However, in opposition to \code{rf} vector arguments are ignored,
  i.e. only the first entry is used.
}

\examples{
## Create a sample of size 1000
x <- urf(n=1000,df1=3,df2=5)
}

\keyword{distribution}
