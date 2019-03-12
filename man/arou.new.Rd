\name{arou.new}

\alias{arou.new}
\alias{aroud.new}

\title{UNU.RAN generator based on Automatic Ratio-Of-Uniforms method (AROU)}

\description{
  UNU.RAN random variate generator for continuous distributions with
  given probability density function (PDF).
  It is based on the Automatic Ratio-Of-Uniforms method (\sQuote{AROU}).

  [Universal] -- Rejection Method.
}

\usage{
arou.new(pdf, dpdf=NULL, lb, ub, islog=FALSE, \dots)
aroud.new(distr)
}
\arguments{
  \item{pdf}{probability density function. (\R function)}
  \item{dpdf}{derivative of \code{pdf}. (\R function)}
  \item{lb}{lower bound of domain;
    use \code{-Inf} if unbounded from left. (numeric)}
  \item{ub}{upper bound of domain;
    use \code{Inf} if unbounded from right. (numeric)}
  \item{islog}{whether \code{pdf} is given as log-density (the
    \code{dpdf} must then be the derivative of the log-density). (boolean)}
  \item{\dots}{(optional) arguments for \code{pdf}.}
  \item{distr}{distribution object. (S4 object of class \code{"unuran.cont"})}
}

\details{
  This function creates an \code{unuran} object based on \sQuote{AROU}
  (Automatic Ratio-Of-Uniforms method). It can be used to draw samples of a
  continuous random variate with given probability density function
  using \code{\link{ur}}.

  The density \code{pdf} must be positive but need not be normalized
  (i.e., it can be any multiple of a density function).
  The derivative \code{dpdf} of the (log-) density is optional. If
  omitted, numerical differentiation is used. Notice, however, that this
  might cause some round-off errors such that the algorithm fails.
  This is in particular the case when the density function is
  provided instead of the log-density.

  The given \code{pdf} must be \eqn{T_{-0.5}}{T_(-0.5)}-concave
  (with implies unimodal densities with tails not higher than
  \eqn{(1/x^2)}; this includes all log-concave distributions).
  
  It is recommended to use the log-density (instead of the density
  function) as this is numerically more stable.
  
  Alternatively, one can use function \code{aroud.new} where the object
  \code{distr} of class \code{"unuran.cont"} must contain all required
  information about the distribution.

  The setup time of this method depends on the given PDF, whereas its
  marginal generation times are almost independent of the target
  distribution.
}

\value{
  An object of class \code{"unuran"}.
}

\seealso{
  \code{\link{ur}},
  \code{\link{tdr.new}},
  \code{\linkS4class{unuran.cont}},
  \code{\link{unuran.new}},
  \code{\linkS4class{unuran}}.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg.
  See Section 4.8 (Automatic Ratio-Of-Uniforms).
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create a sample of size 100 for a Gaussian distribution
pdf <- function (x) { exp(-0.5*x^2) }
gen <- arou.new(pdf=pdf, lb=-Inf, ub=Inf)
x <- ur(gen,100)

## Create a sample of size 100 for a 
## Gaussian distribution (use logPDF)
logpdf <- function (x) { -0.5*x^2 }
gen <- arou.new(pdf=logpdf, islog=TRUE, lb=-Inf, ub=Inf)
x <- ur(gen,100)

## Same example but additionally provide derivative of log-density
## to prevent possible round-off errors
logpdf <- function (x) { -0.5*x^2 }
dlogpdf <- function (x) { -x }
gen <- arou.new(pdf=logpdf, dpdf=dlogpdf, islog=TRUE, lb=-Inf, ub=Inf)
x <- ur(gen,100)

## Draw sample from Gaussian distribution with mean 1 and
## standard deviation 2. Use 'dnorm'.
gen <- arou.new(pdf=dnorm, lb=-Inf, ub=Inf, mean=1, sd=2)
x <- ur(gen,100)

## Draw a sample from a truncated Gaussian distribution
## on domain [5,Inf)
logpdf <- function (x) { -0.5*x^2 }
gen <- arou.new(pdf=logpdf, lb=5, ub=Inf, islog=TRUE)
x <- ur(gen,100)

## Alternative approach
distr <- udnorm()
gen <- aroud.new(distr)
x <- ur(gen,100)

}

\keyword{datagen}
\keyword{distribution}
