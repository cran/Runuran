\name{itdr.new}
\docType{function}

\alias{itdr.new}

\title{UNU.RAN generator based on Inverse Transformed Density Rejection (ITDR)}

\description{
  UNU.RAN random variate generator for continuous distributions with
  given probability density function (PDF).
  It is based on the Inverse Transformed Density Rejection method (ITDR).

  [Universal] -- Rejection Method.
}

\usage{
itdr.new(pdf, dpdf, pole, lb=0, ub=Inf, islog=FALSE, ...)
}

\arguments{
  \item{pdf}{probability density function. (\R function)}
  \item{dpdf}{derivative of \code{pdf}. (\R function)}
  \item{pole}{pole of distribution. (numeric)}
  \item{lb, ub}{lower and upper bound of domain. (numeric)}
  \item{islog}{whether \code{pdf} is given as log-density (the
    \code{dpdf} must then be the derivative of the
    log-density). (boolean)}
  \item{...}{(optional) arguments for \code{pdf}}
}

\details{
  This function creates a \code{unuran} object based on \dQuote{ITDR}
  (Inverse Transformed Density Rejection). It can be used to draw samples of a
  continuous random variate with given probability density function
  using \code{\link{ur}}.
  
  The density must be provided by a function \code{pdf} which must
  return non-negative numbers and which need not be normalized (i.e., it
  can be any multiple of a density function).
  Moreover, it must be monotone on its domain.
  The algorithm is especially designed for distributions with unbounded
  densities. Thus the algorithm need the position of the \code{pole}.

  The derivative \code{dpdf} is essential. (Numerical derivation does
  not work as it results in serious round-off errors.)
  
  The setup time of this method depends on the given PDF, whereas its
  marginal generation times are almost independent of the target
  distribution.
}

\seealso{
  \code{\link{ur}}, \code{\link{unuran.new}}, \code{\linkS4class{unuran}}.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2007):
  Inverse transformed density rejection for unbounded monotone
  densities.
  ACM Trans. Model. Comput. Simul. 17(4), Article 18, 16 pages.
  DOI = 10.1145/1276927.1276931
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

\examples{
## Create a sample of size 100 for a Gamma(-0.5) distribution
pdf <- function (x) { x^(-0.5)*exp(-x) }
dpdf <- function (x) { (-x^(-0.5) - 0.5*x^(-1.5))*exp(-x) }
gen <- itdr.new(pdf=pdf, dpdf=dpdf, lb=0, pole=0)
x <- ur(gen,100)

}

\keyword{distribution}
\keyword{datagen}

