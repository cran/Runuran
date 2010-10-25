\name{ars.new}

\alias{ars.new}
\alias{arsd.new}

\title{UNU.RAN generator based on Adaptive Rejection Sampling (ARS)}

\description{
  UNU.RAN random variate generator for continuous distributions with
  given probability density function (PDF).
  It is based on Adaptive Rejection Sampling (\sQuote{ARS}).

  [Universal] -- Rejection Method.
}

\usage{
ars.new(logpdf, dlogpdf=NULL, lb, ub, \dots)
arsd.new(distr)
}
\arguments{
  \item{logpdf}{log-density function. (\R function)}
  \item{dlogpdf}{derivative of \code{logpdf}. (\R function)}
  \item{lb}{lower bound of domain;
    use \code{-Inf} if unbounded from left. (numeric)}
  \item{ub}{upper bound of domain;
    use \code{Inf} if unbounded from right. (numeric)}
  \item{\dots}{(optional) arguments for \code{logpdf}.}
  \item{distr}{distribution object. (S4 object of class \code{"unuran.cont"})}
}

\details{
  This function creates a \code{unuran} object based on \sQuote{ARS}
  (Adaptive Rejection Sampling). It can be used to draw samples from
  continuous distributions with given probability density function
  using \code{\link{ur}}.

  Function \code{logpdf} is the logarithm the density function of the
  target distribution. It must be a concave function (i.e., the
  distribution must be log-concave).
  However, it need not be normalized (i.e., it can be a log-density plus
  some arbitrary constant).

  The derivative \code{dlogpdf} of the log-density is optional. If
  omitted, numerical differentiation is used. Notice, however, that this
  might cause some round-off errors such that the algorithm fails.

  Alternatively, one can use function \code{arsd.new} where the object
  \code{distr} of class \code{"unuran.cont"} must contain all required
  information about the distribution.

  The setup time of this method depends on the given PDF, whereas its
  marginal generation times are almost independent of the target
  distribution.
  
  \sQuote{ARS} is a special case of method \sQuote{TDR}
  (see \code{\link{tdr.new}}). It is a bit slower and less
  flexible but numerically more stable. In particular, it is useful if
  one wants to sample from truncated distributions with extreme
  truncation points; or when the integral of the given \dQuote{density}
  function is only known to be extremely large or small.
  However, this assumes that the log-density is computed
  analytically and not by just using \code{log(pdf(x))}.
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
  See Chapter 4 (Tranformed Density Rejection).
  
  W. R. Gilks and P. Wild (1992):
  Adaptive rejection sampling for Gibbs sampling.
  Applied Statistics 41(2), pp. 337--348.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create a sample of size 100 for a 
## Gaussian distribution (use logPDF)
lpdf <- function (x) { -0.5*x^2 }
gen <- ars.new(logpdf=lpdf, lb=-Inf, ub=Inf)
x <- ur(gen,100)

## Same example but additionally provide derivative of log-density
## to prevent possible round-off errors
lpdf <- function (x) { -0.5*x^2 }
dlpdf <- function (x) { -x }
gen <- ars.new(logpdf=lpdf, dlogpdf=dlpdf, lb=-Inf, ub=Inf)
x <- ur(gen,100)

## Draw a sample from a truncated Gaussian distribution
## on domain [100,Inf)
lpdf <- function (x) { -0.5*x^2 }
gen <- ars.new(logpdf=lpdf, lb=50, ub=Inf)
x <- ur(gen,100)

## Alternative approach
distr <- udnorm()
gen <- arsd.new(distr)
x <- ur(gen,100)

}

\keyword{datagen}
\keyword{distribution}
