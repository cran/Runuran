\name{mixt.new}
\alias{mixt.new}

\title{UNU.RAN generator for finite mixture of distributions}

\description{
  UNU.RAN random variate generator for a finite mixture of continuous or
  discrete distributions.
  The components are given as \code{unuran} objects.

  [Universal] -- Composition Method.
}

\usage{
mixt.new(prob, comp, inversion=FALSE)
}
\arguments{
  \item{prob}{weights of mixture (\dQuote{probabilities}); these must be
    non-negative numbers but need not sum to 1. (numeric vector)}
  \item{comp}{components of mixture.
    (list of S4 object of class \code{"unuran"})}
  \item{inversion}{whether inversion method should be used. (boolean)}
}

\details{
  Given a set of probability density functions
  \eqn{p_1(x),\dots,p_n(x)} (called the mixture components) and weights
  \eqn{w_1,\dots,w_n} such that \eqn{w_i \ge 0}{w_i >= 0} and
  \eqn{\sum w_i = 1}{w_1+\dots+w_n=1},
  the sum 
  \deqn{q(x) = \sum_{i=1}^n \, w_i \, p_i(x)}{q(x) = w_1*p_1(x) + \dots + w_n*p_n(x)}
  is called the mixture density.

  Function \code{mixt.new} creates an \code{unuran} object for a finite
  mixture of continuous or discrete univariate distributions.
  It can be used to draw samples of a continuous random variate using
  \code{\link{ur}}.

  The weights \code{prob} must be a vector of non-negative numbers (not
  all equal to 0) but need not sum to 1.

  \code{comp} is a list of \code{"unuran"} generator objects. Each of
  which must sample from a continuous or discrete univariate
  distribution.

  If \code{inversion} is \code{TRUE}, then the inversion method is used
  for sampling from the mixture distribution.
  However, the following conditions must be satisfied:
  \itemize{
    \item
    Each component (\code{unuran} object) must use implement an
    inversion method (i.e., the quantile funtion \code{\link{uq}} must
    work).
    \item
    The domains of the components must not overlapping.
    \item
    The components must be order with respect to their domains.
  }
  If one of these conditions is violated, then initialization of the
  mixture object fails.
  
  The setup time is fast, whereas its marginal generation times strongly
  depend on the average generation times of its components.
}

\value{
  An object of class \code{"unuran"}.
}

\note{
  Each component in \code{comp} must correspond to a continuous or
  discrete univariate distribution. In particular this also includes
  mixtures of distributions. Thus mixtures can also be defined
  recursively.

  Moreover, none of these components must be packed
  (see \code{\link{unuran.packed}}).
}

\seealso{
  \code{\link{ur}}, \code{\link{uq}},
  \code{\link{unuran.new}},
  \code{\linkS4class{unuran}}.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg.
  See Section 2.3 (Composition).
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create a mixture of an Exponential and a Half-normal distribution
unr1 <- unuran.new(udnorm(lb=-Inf, ub=0))
unr2 <- unuran.new(udexp())
mix <- mixt.new( c(1,1), c(unr1, unr2) )
x <- ur(mix,100)

## Now use inversion method:
## It is important that
##  1. we use a inversion for each component
##  2. the domains to not overlap
##  3. the components are ordered with respect to their domains
unr1 <- pinvd.new(udnorm(lb=-Inf, ub=0))
unr2 <- pinvd.new(udexp())
mix <- mixt.new( c(1,1), c(unr1, unr2), inversion=TRUE )
x <- ur(mix,100)

## We also can compute the inverse distribution function
##x <- uq(mix,0.90)

## Create a mixture of Exponential and Geometric distrbutions
unr1 <- unuran.new(udexp())
unr2 <- unuran.new(udgeom(0.7))
mix <- mixt.new( c(0.6,0.4), c(unr1, unr2) )
x <- ur(mix,100)

}

\keyword{datagen}
\keyword{distribution}
