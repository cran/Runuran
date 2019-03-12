\name{tabl.new}

\alias{tabl.new}
\alias{tabld.new}

\title{UNU.RAN generator based on TABLe based Rejection (TABL)}

\description{
  UNU.RAN random variate generator for continuous distributions with
  given probability density function (PDF).
  It is based on the TABLe based rejection method (\sQuote{TABL}).

  [Universal] -- Rejection Method.
}

\usage{
tabl.new(pdf, lb, ub, mode, islog=FALSE, \dots)
tabld.new(distr)
}
\arguments{
  \item{pdf}{probability density function. (\R function)}
  \item{lb}{lower bound of domain;
    use \code{-Inf} if unbounded from left. (numeric)}
  \item{ub}{upper bound of domain;
    use \code{Inf} if unbounded from right. (numeric)}
  \item{mode}{location of the mode. (numeric)}
  \item{islog}{whether \code{pdf} is given as log-density (the
    \code{dpdf} must then be the derivative of the log-density). (boolean)}
  \item{\dots}{(optional) arguments for \code{pdf}.}
  \item{distr}{distribution object. (S4 object of class \code{"unuran.cont"})}
}

\details{
  This function creates an \code{unuran} object based on \sQuote{TABL}
  (TABLe based rejection). It can be used to draw samples of a
  continuous random variate with given probability density function
  using \code{\link{ur}}.

  The density \code{pdf} must be positive but need not be normalized
  (i.e., it can be any multiple of a density function).

  The given \code{pdf} must be unimodal.
  
  Alternatively, one can use function \code{tabld.new} where the object
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
  See Section 5.1 (\dQuote{Ahrens Method}).
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create a sample of size 100 for a Gaussian distribution
pdf <- function (x) { exp(-0.5*x^2) }
gen <- tabl.new(pdf=pdf, lb=-Inf, ub=Inf, mode=0)
x <- ur(gen,100)

## Create a sample of size 100 for a 
## Gaussian distribution (use logPDF)
logpdf <- function (x) { -0.5*x^2 }
gen <- tabl.new(pdf=logpdf, islog=TRUE, lb=-Inf, ub=Inf, mode=0)
x <- ur(gen,100)

## Draw sample from Gaussian distribution with mean 1 and
## standard deviation 2. Use 'dnorm'.
gen <- tabl.new(pdf=dnorm, lb=-Inf, ub=Inf, mode=1, mean=1, sd=2)
x <- ur(gen,100)

## Draw a sample from a truncated Gaussian distribution
## on domain [5,Inf)
logpdf <- function (x) { -0.5*x^2 }
gen <- tabl.new(pdf=logpdf, lb=5, ub=Inf, mode=5, islog=TRUE)
x <- ur(gen,100)

## Alternative approach
distr <- udnorm()
gen <- tabld.new(distr)
x <- ur(gen,100)

}

\keyword{datagen}
\keyword{distribution}
