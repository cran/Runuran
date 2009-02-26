\name{tdr.new}
\docType{function}

\alias{tdr.new}

\title{UNU.RAN generator based on Transformed Density Rejection (TDR)}

\description{
  UNU.RAN random variate generator for continuous distributions with
  given probability density function (PDF).
  It is based on the Transformed Density Rejection method (\sQuote{TDR}).

  [Universal] -- Rejection Method.
}

\usage{
tdr.new(pdf, dpdf=NULL, lb, ub, islog=FALSE, ...)
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
  \item{...}{(optional) arguments for \code{pdf}}
}

\details{
  This function creates an \code{unuran} object based on \sQuote{TDR}
  (Transformed Density Rejection). It can be used to draw samples of a
  continuous random variate with given probability density function
  using \code{\link{ur}}.
  
  The density must be provided by a function \code{pdf} which must
  return non-negative numbers and which need not be normalized (i.e., it
  can be any multiple of a density function).
  Moreover, the given function must be \eqn{T_{-0.5}}{T_(-0.5)}-concave
  (i.e., unimodal densities with tails not higher than
  \eqn{(1/x^2)}; this includes all log-concave distributions).
  
  It is recommended to use the log-density (instead of the density
  function) as this is numerically more stable.
  
  The derivative \code{dpdf} is optional. If omitted, numerical
  differentiation is used. Notice, however, that this might cause some
  round-off errors such that generation does not work. This is in
  particular the case when the density function is provided instead of
  the log-density.
  
  The setup time of this method depends on the given PDF, whereas its
  marginal generation times are almost independent of the target
  distribution.

  There exists a variant of \sQuote{TDR} which is numerically more
  stable (albeit a bit slower and less flexible) which is avaible
  via the \code{\link{ars.new}} function.
}

\seealso{
  \code{\link{ur}}, \code{\link{ars.new}},
  \code{\link{unuran.new}}, \code{\linkS4class{unuran}}.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg.
  See Chapter 4 (Tranformed Density Rejection).
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

\examples{
## Create a sample of size 100 for a Gaussian distribution
pdf <- function (x) { exp(-0.5*x^2) }
gen <- tdr.new(pdf=pdf, lb=-Inf, ub=Inf)
x <- ur(gen,100)

## Create a sample of size 100 for a 
## Gaussian distribution (use logPDF)
logpdf <- function (x) { -0.5*x^2 }
gen <- tdr.new(pdf=logpdf, islog=TRUE, lb=-Inf, ub=Inf)
x <- ur(gen,100)

## Same example but additionally provide derivative of log-density
## to prevent possible round-off errors
logpdf <- function (x) { -0.5*x^2 }
dlogpdf <- function (x) { -x }
gen <- tdr.new(pdf=logpdf, dpdf=dlogpdf, islog=TRUE, lb=-Inf, ub=Inf)
x <- ur(gen,100)

## Draw sample from Gaussian distribution with mean 1 and
## standard deviation 2. Use 'dnorm'.
gen <- tdr.new(pdf=dnorm, lb=-Inf, ub=Inf, mean=1, sd=2)
x <- ur(gen,100)

## Draw a sample from a truncated Gaussian distribution
## on domain [5,Inf)
logpdf <- function (x) { -0.5*x^2 }
gen <- tdr.new(pdf=logpdf, lb=5, ub=Inf, islog=TRUE)
x <- ur(gen,100)

}

\keyword{distribution}
\keyword{datagen}

