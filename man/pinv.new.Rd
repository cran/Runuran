\name{pinv.new}

\alias{pinv.new}
\alias{pinvd.new}

\title{UNU.RAN generator based on Polynomial interpolation of INVerse CDF (PINV)}

\description{
  UNU.RAN random variate generator for continuous distributions with
  given probability density function (PDF) or cumulative distribution
  function (CDF).
  It is based on the Polynomial interpolation of INVerse CDF (\sQuote{PINV}).

  [Universal] -- Inversion Method.
}

\usage{
pinv.new(pdf, cdf, lb, ub, islog=FALSE, center=0,
         uresolution=1.e-10, smooth=FALSE, \dots)
pinvd.new(distr, uresolution=1.e-10, smooth=FALSE)
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
  \item{\dots}{(optional) arguments for \code{pdf} and \code{cdf}.}
  \item{distr}{distribution object. (S4 object of class \code{"unuran.cont"})}
  \item{uresolution}{maximal acceptable u-error. (numeric)}
  \item{smooth}{whether the inverse CDF is differentiable. (boolean)}
}

\details{
  This function creates an \code{unuran} object based on \sQuote{PINV}
  (Polynomial interpolation of INVerse CDF). It can be used to draw
  samples of a continuous random variate with given probability density
  function \code{pdf} or cumulative distribution function \code{cdf}
  by means of \code{\link{ur}}.
  It also allows to compute quantiles by means of \code{\link{uq}}.

  Function \code{pdf} must be positive but need not be normalized
  (i.e., it can be any multiple of a density function).
  The set of points where the \code{pdf} is strictly positive
  must be connected.
  The \code{center} is a point where the \code{pdf} is not too
  small, e.g., (a point near) the mode of the distribution.
  
  If the density \code{pdf} is given, then the algorithm automatically
  computes the CDF using Gauss-Lobatto integration.
  If the \code{cdf} is given but not the \code{pdf} then the CDF is used
  instead of the PDF. However, we found in our experiments that using
  the PDF is numerically more stable.
  
  Alternatively, one can use function \code{pinvd.new} where the object
  \code{distr} of class \code{"unuran.cont"} must contain all required
  information about the distribution.

  The algorithm approximates the inverse of the CDF of the
  distribution by means of Newton interpolation between carefully
  selected nodes. The approxiating functing is thus continuous.
  Argument \code{smooth} controls whether this function is also
  differentiable(\dQuote{smooth}) at the nodes.
  Using \code{smooth=TRUE} requires the \code{pdf} of the
  distribution. It results in a higher setup time and memory
  consumption. Thus using \code{smooth=TRUE} is not
  \emph{not recommended}, unless differentiability is important.

  The approximation error is estimated by means of the the
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

\section{Remark}{
  Using function \code{\link{up}} generator objects that implement
  method \sQuote{PINV} may also be used to approximate the cumulative
  distribution function of the given distribution when only the density
  is given. The approximation error is about one tenth of the requested
  \code{uresolution}.
}

\value{
  An object of class \code{"unuran"}.
}

\seealso{
  \code{\link{ur}}, \code{\link{uq}}, \code{\link{up}},
  \code{\linkS4class{unuran.cont}},
  \code{\link{unuran.new}}, 
  \code{\linkS4class{unuran}}.
}

\references{
  G. Derflinger, W. H\"ormann, and J. Leydold (2010):
  Random variate generation by numerical inversion when only the density
  is known. 
  ACM Trans. Model. Comput. Simul., 20:4, #18
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
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

## We also can force a smoother approximation
gen <- pinv.new(pdf=dnorm, lb=-Inf, ub=Inf, smooth=TRUE)
x <- ur(gen,100)

## Alternative approach
distr <- udnorm()
gen <- pinvd.new(distr)
x <- ur(gen,100)

}

\keyword{datagen}
\keyword{distribution}
