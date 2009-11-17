\name{dgt.new}
\docType{function}

\alias{dgt.new}
\alias{dgtd.new}

\title{UNU.RAN generator based on table guided discrete inversion (DGT)}

\description{
  UNU.RAN random variate generator for discrete distributions with given
  probability vector. 
  It applies the Guide-Table Method for discrete inversion (\sQuote{DGT}).

  [Universal] -- Inversion Method.
}

\usage{
dgt.new(pv, from=1)
dgtd.new(distr)
}
\arguments{
  \item{pv}{vector of non-negative numbers (need not sum to 1). (numeric vector)}
  \item{from}{index of first entry in vector. (integer)}
  \item{distr}{distribution object. (S4 object of class \code{"unuran.discr"})}
}

\details{
  This function creates an \code{unuran} object based on \sQuote{DGT}
  (Discrete Guide-Table method). It can be used to draw samples of a
  discrete random variate with given probability vector
  using \code{\link{ur}}.
  It also allows to compute quantiles by means of \code{\link{uq}}.

  Vector \code{pv} must be postive but need not be normalized
  (i.e., it can be any multiple of a probability vector).

  The method runs fast in constant time, i.e., marginal sampling
  times do not depend on the length of the given probability vector.
  Whereas their setup times grow linearly with this length.
  
  Notice that the range of random variates is
  \code{from:(from+length(pv)-1)}.

  Alternatively, one can use function \code{dgtd.new} where the object
  \code{distr} of class \code{"unuran.discr"} must contain all required
  information about the distribution.
}

\value{
  An object of class \code{"unuran"}.
}

\seealso{
  \code{\link{ur}}, \code{\link{uq}},
  \code{\linkS4class{unuran.discr}},
  \code{\link{unuran.new}},
  \code{\linkS4class{unuran}}.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg.
  See Section 3.1.2 (Indexed Search).
  
  H.C. Chen and Y. Asau (1974):
  On generating random variates from an empirical distribution.
  AIIE Trans. 6, pp.163--166.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create a sample of size 100 for a 
## binomial distribution with size=115, prob=0.5
gen <- dgt.new(pv=dbinom(0:115,115,0.5),from=0)
x <- ur(gen,100)

## Alternative approach
distr <- udbinom(size=100,prob=0.3)
gen <- dgtd.new(distr)
x <- ur(gen,100)

}

\keyword{datagen}
\keyword{distribution}
