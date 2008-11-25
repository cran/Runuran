\name{urdgt}
\docType{function}

\alias{urdgt}
\alias{urdau}

\title{UNU.RAN discrete random variate generator}

\description{
  UNU.RAN random variate generator for discrete distributions with given
  probability vector.
  It applies the Guide-Table Method (\code{urdgt}) for discrete
  inversion or the Alias-Urn Method (\code{urdau}).
  
  [ \strong{Deprecated}.
  Use \code{\link{dgt.new}} or \code{\link{dau.new}} instead. ]
}

\usage{
urdgt(n, probvector, from = 0, by = 1)
urdau(n, probvector, from = 0, by = 1)
}
\arguments{
  \item{n}{size of required sample.}
  \item{probvector}{vector of non-negative numbers (need not sum to 1).}
  \item{from}{number corresponding to the first probability in probvector.}
  \item{by}{"from + (k-1)*by" is the number corresponding to the
    \eqn{k}-th probability in probvector.}
}

\details{
  These routines generate a sample of discrete random variates with
  given probability vector. This vector must be provided by
  \code{probvector}, a vector of non-negative numbers which need not
  necessarily sum up to 1.

  Method \dQuote{DGT} uses a guide-table based inversion method.

  Method \dQuote{DAU} implements the \emph{Alias-Urn} method.

  Both methods run in constant time, i.e., the marginal sampling times
  do not depend on the length of the given probability vector.
  Whereas their setup times grow linearly with this length.
}

\seealso{
  \code{\link{dau.new}}, \code{dgt.new}, and \code{\link{ur}}
  for replacement.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg.

  H.C. Chen and Y. Asau (1974):
  On generating random variates from an empirical distribution.
  AIIE Trans. 6, pp.163--166.

  A.J. Walker (1977):
  An efficient method for generating discrete random variables with
  general distributions.
  ACM Trans. Model. Comput. Simul. 3, pp.253--256.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

\note{
urdgt() and urdau() are very fast for probvector not longer than about 1000.
}

%%\examples{}

\keyword{distribution}
\keyword{datagen}

