\name{dau.new}
\docType{function}

\alias{dau.new}

\title{UNU.RAN generator based on the Alias method (DAU)}

\description{
  UNU.RAN random variate generator for discrete distributions with given
  probability vector.
  It applies the Alias-Urn method (\sQuote{DAU}).

  [Universal] -- Patchwork Method.
}

\usage{
dau.new(pv, from=1)
}
\arguments{
  \item{pv}{vector of non-negative numbers (need not sum to 1). (numeric vector)}
  \item{from}{index of first entry in vector. (integer)}
}

\details{
  This function creates a \code{unuran} object based on \sQuote{DAU}
  (Discrete Alias-Urn method). It can be used to draw samples of a
  discrete random variate with given probability vector
  using \code{\link{ur}}.
  
  The probabilities must be provided by a vector \code{pv} which must
  contain non-negative numbers and which need not necessarily sum up to
  1.

  The method runs fast in constant time, i.e., marginal sampling
  times do not depend on the length of the given probability vector.
  Whereas their setup times grow linearly with this length.
  
  Notice that the range of random variates is
  \code{from:(from+length(pv)-1)}.
}

\seealso{
  \code{\link{ur}}, \code{\link{unuran.new}}, \code{\linkS4class{unuran}}.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg.
  See Section 3.2 (The Alias Method).

  A.J. Walker (1977):
  An efficient method for generating discrete random variables with
  general distributions.
  ACM Trans. Model. Comput. Simul. 3, pp.253--256.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

\examples{
## Create a sample of size 100 for a 
## binomial distribution with size=115, prob=0.5
gen <- dau.new(pv=dbinom(0:115,115,0.5))
x <- ur(gen,100)
}

\keyword{distribution}
\keyword{datagen}

