\name{srou.new}
\docType{function}

\alias{srou.new}
\alias{sroud.new}

\title{UNU.RAN generator based on Simple Ratio-Of-Uniforms Method (SROU)}

\description{
  UNU.RAN random variate generator for continuous distributions with
  given probability density function (PDF).
  It is based on the Simple Ratio-Of-Uniforms Method (\sQuote{SROU}).

  [Universal] -- Rejection Method.
}

\usage{
srou.new(pdf, lb, ub, mode, area, islog=FALSE, r=1, \dots)
sroud.new(distr, r=1)
}
\arguments{
  \item{pdf}{probability density function. (\R function)}
  \item{lb}{lower bound of domain;
    use \code{-Inf} if unbounded from left. (numeric)}
  \item{ub}{upper bound of domain;
    use \code{Inf} if unbounded from right. (numeric)}
  \item{mode}{location of the mode. (numeric)}
  \item{area}{area below \code{pdf}. (numeric)}
  \item{islog}{whether \code{pdf} is given as log-density (the
    \code{dpdf} must then be the derivative of the log-density). (boolean)}
  \item{\dots}{(optional) arguments for \code{pdf}.}
  \item{distr}{distribution object. (S4 object of class \code{"unuran.cont"})}
  \item{r}{adjust algorithm to heavy-tailed distribution. (numeric)}
}

\details{
  This function creates a \code{unuran} object based on \sQuote{SROU}
  (Simple Ratio-Of-Uniforms Method). It can be used to draw samples of a
  continuous random variate with given probability density function
  using \code{\link{ur}}.

  The density \code{pdf} must be positive but need not be normalized
  (i.e., it can be any multiple of a density function).
  It must be \eqn{T_c}{T_c}-concave for
  \eqn{c = -r/(r+1)}{c = -r/(r+1)}; this includes all log-concave
  distributions.

  The (exact) location of the \code{mode} and the area below
  the \code{pdf} are essential.

  Alternatively, one can use function \code{sroud.new} where the object
  \code{distr} of class \code{"unuran.cont"} must contain all required
  information about the distribution.
  
  The acceptance probability decreases with increasing parameter
  \code{r}. Thus it should be as small as possible. On the other hand it
  must be sufficiently large for heavy tailed distributions.
  If possible, use the default \code{r=1}.

  Compared to \code{\link{tdr.new}} it has much slower marginal
  generation times but has a faster setup and is numerically more
  robust. Moreover, It also works for unimodal distributions with tails
  that are heavier than those of the Cauchy distribution.
}

\value{
  An object of class \code{"unuran"}.
}

\seealso{
  \code{\link{ur}},
  \code{\linkS4class{unuran.cont}},
  \code{\link{unuran.new}},
  \code{\linkS4class{unuran}}.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg.
  Sections 6.3 and 6.4.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create a sample of size 100 for a Gaussian distribution.
pdf <- function (x) { exp(-0.5*x^2) }
gen <- srou.new(pdf=pdf, lb=-Inf, ub=Inf, mode=0, area=2.506628275)
x <- ur(gen,100)

## Create a sample of size 100 for a Gaussian distribution.
## Use 'dnorm'.
gen <- srou.new(pdf=dnorm, lb=-Inf, ub=Inf, mode=0, area=1)
x <- ur(gen,100)

## Alternative approach
distr <- udnorm()
gen <- sroud.new(distr)
x <- ur(gen,100)

}

\keyword{datagen}
\keyword{distribution}
