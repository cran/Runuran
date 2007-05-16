\name{urtdr}
\docType{function}

\alias{urtdr}

\title{UNU.RAN transformend density rejection generator}

\description{
  UNU.RAN random variate generator for continuous distributions with
  given probability density function (PDF).
  It is based on the Transformed Density Rejection method (TDR).
}

\usage{
urtdr(n, pdf, dpdf, lb=-Inf, ub=Inf, islog=TRUE, ...)
}
\arguments{
  \item{n}{size of required sample.}
  \item{pdf}{probability density function. (\R function)}
  \item{dpdf}{derivative of the \code{pdf}. (\R function)}
  \item{islog}{whether the \code{pdf} is given by tits logarithm (the
    \code{dpdf} is then the derivative of the logarithm). (boolean)}
  \item{lb, ub}{lower and upper bound of domain. (numeric)}
  \item{...}{(optional) arguments for \code{pdf}}
}

\details{
  This routine generates a sample of continuous random variates with
  given probability density function. This function must be provided by
  \code{pdf}, a function which must return non-negative numbers and
  which need not be normalized (i.e., it can be any multiple of a
  density function).
  Moreover, the given function must be \eqn{T_{-0.5}}{T_(-0.5)}-concave
  distributions (i.e., unimodal densities with tails not higher than
  \eqn{(1/x^2)}; this includes all log-concave distributions).
  The \code{dpdf} is optional. If omitted numerical derivation is used.
  Notice, however, that this might cause some round-off errors such that
  generation does not work. This is in particular the case when the
  density function is provided instead of its logarithm.
  
  It uses method \dQuote{TDR} (Transformed Density Rejection).
  The setup time of this method depends on the given PDF, whereas its
  marginal generation times are almost independent of the target
  distribution.
}

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation, \code{\linkS4class{unuran}} for the UNU.RAN class.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

\examples{
## Create a sample of size 10000 for a 
## Gaussian distribution (use logPDF)
logpdf <- function (x) { -0.5*x^2 }
dlogpdf <- function (x) { -x }
x <- urtdr(1000, pdf=logpdf, dpdf=dlogpdf)

## Create a sample of size 10000 for a 
## Gaussian distribution (use logPDF)
logpdf <- function (x) { -0.5*x^2 }
x <- urtdr(1000, pdf=logpdf, islog=TRUE)

## Create a sample of size 10000 for a 
## Gaussian(0,1) distribution (use dnorm)
x <- urtdr(1000, pdf=dnorm, mean=0, sd=1, islog=FALSE)

}

\keyword{distribution}
\keyword{datagen}

