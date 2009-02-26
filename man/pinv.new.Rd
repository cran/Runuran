\name{pinv.new}
\docType{function}

\alias{pinv.new}

\title{UNU.RAN generator based on Polynomial interpolation of INVerse CDF (PINV)}

\description{
  UNU.RAN random variate generator for continuous distributions with
  given probability density function (PDF) or cumulative distribution
  function (CDF).
  It is based on the Polynomial interpolation of INVerse CDF (\sQuote{PINV}).

  [Universal] -- Inversion Method.
}

\usage{
pinv.new(pdf, cdf, lb, ub, islog=FALSE, center=0, uresolution=1.e-10, ...)
}
\arguments{
  \item{pdf}{probability density function. (\R function)}
  \item{cdf}{cumulative distribution function. (\R function)}
  \item{lb}{lower bound of domain;
    use \code{-Inf} if unbounded from left. (numeric)}
  \item{ub}{upper bound of domain;
    use \code{Inf} if unbounded from right. (numeric)}
  \item{islog}{whether \code{pdf} and \code{cdf} are given by their
    corresponding logarithms. (boolean)}
  \item{center}{\dQuote{typical} point of distribution. (numeric)}
  \item{uresolution}{maximal acceptable u-error. (numeric)}
  \item{...}{(optional) arguments for \code{pdf} and \code{cdf}.}
}

\details{
  This function creates an \code{unuran} object based on \sQuote{PINV}
  (Polynomial interpolation of INVerse CDF). It can be used to draw
  samples of a continuous random variate with given probability density
  function \code{pdf} or cumulative distribution function \code{cdf}
  by means of \code{\link{ur}}.
  It also allows to compute quantiles by means of \code{\link{uq}}.
  
  \code{pdf} is (a multiple of) a density, i.e., it must return
  non-negative numbers but need not integrate to 1.
  However, the set of points where the \code{pdf} is strictly positive
  must be connected.
  The algorithm then automatically computes the CDF using Gauss-Lobatto
  integration.
  
  If the \code{cdf} is given but not the \code{pdf} then the CDF is used
  instead of the PDF. However, we found in our experiments that using
  the PDF is numerically more stable.
  
  The \code{center} must be a point where the \code{pdf} is not too
  small, e.g., (a point near) the mode of the distribution.
  By default \code{0.} is assumed.
  
  The algorithms approximates the inverse of the CDF of the
  distribution. The distribution error is estimated by means of the the
  u-error, i.e., \eqn{|CDF(G(U)) - U|},
  where \eqn{G} denotes the approximation of the inverse CDF.
  The error can be controlled by means of argument \code{uresolution}.

  When sampling from truncated distributions with extreme
  truncation points, it is recommended to provide the log-density
  by setting \code{islog=TRUE}. Then the algorithm is numerically more
  stable.
  
  The setup time of this method depends on the given PDF, whereas its
  marginal generation times are independent of the target
  distribution.
}

\seealso{
  \code{\link{ur}}, \code{\link{uq}},
  \code{\link{unuran.new}}, \code{\linkS4class{unuran}}.
}

\references{
  G. Derflinger, W. H\"ormann, and J. Leydold (2009):
  in preparation.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

\examples{
## Create a sample of size 100 for a Gaussian distribution
pdf <- function (x) { exp(-0.5*x^2) }
gen <- pinv.new(pdf=pdf, lb=-Inf, ub=Inf)
x <- ur(gen,100)

## Create a sample of size 100 for a 
## Gaussian distribution (use logPDF)
logpdf <- function (x) { -0.5*x^2 }
gen <- pinv.new(pdf=logpdf, islog=TRUE, lb=-Inf, ub=Inf)
x <- ur(gen,100)

## Draw sample from Gaussian distribution with mean 1 and
## standard deviation 2. Use 'dnorm'.
gen <- pinv.new(pdf=dnorm, lb=-Inf, ub=Inf, mean=1, sd=2)
x <- ur(gen,100)

## Draw a sample from a truncated Gaussian distribution
## on domain [2,Inf)
gen <- pinv.new(pdf=dnorm, lb=2, ub=Inf)
x <- ur(gen,100)

## Improve the accuracy of the approximation
gen <- pinv.new(pdf=dnorm, lb=-Inf, ub=Inf, uresolution=1e-15)
x <- ur(gen,100)

## We have to provide a 'center' when PDF (almost) vanishes at 0.
gen <- pinv.new(pdf=dgamma, lb=0, ub=Inf, center=4, shape=5)
x <- ur(gen,100)
}

\keyword{distribution}
\keyword{datagen}

