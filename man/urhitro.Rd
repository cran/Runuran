\name{urhitro}
\docType{function}

\alias{urhitro}

\title{UNU.RAN transformend density rejection generator}

\description{
  UNU.RAN random variate generator for continuous multivariate
  distributions with given probability density function (PDF).
  It is based on the Hit-and-Run algorithm in combinaton with the
  Ratio-of-Uniforms method (HITRO).
  
  [ \strong{Deprecated}.
  Use \code{\link{hitro.new}} instead. ]
}

\usage{
urhitro(n, dim=1, pdf, mode=NULL, center=NULL, ll=NULL, ur=NULL, thinning=1, burnin = 0, ...)
}
\arguments{
  \item{n}{size of required sample.}
  \item{dim}{number of dimensions of the distribution. (integer)}
  \item{pdf}{probability density function. (\R function)}
  \item{mode}{location of the mode -- optional. (numeric vector)}
  \item{center}{point in \dQuote{typical} region of distribution,
    e.g. the approximate location of the mode. If omitted the
    \code{mode} is implicitly used. If the \code{mode} is not given
    either, the origin is used. (numeric vector -- optional)}
  \item{ll,ur}{lower left and upper right vertex of a rectangular
    domain of the \code{pdf}. The domain is only set if both vertices
    are not \code{NULL}. Otherwise, the domain is unbounded by
    default. (numeric vectors)}
  \item{thinning}{thinning factor. (positive integer)}
  \item{burnin}{lenght of burin-in phase. (positive integer)}
  \item{...}{(optional) arguments for \code{pdf}}
}

\details{
  This method is \strong{deprecated}. Please use \code{\link{hitro.new}}
  together with \code{\link{ur}} instead.
  
  For further details see \code{\link{hitro.new}}.
}

\seealso{
  \code{\link{hitro.new}} and \code{\link{ur}} for replacement.
}

\references{
  R. Karawatzki, J. Leydold, and K. P\"otzelberger (2005):
  Automatic Markov Chain Monte Carlo Procedures for Sampling from
  Multivariate Distributions.
  Research Report Series / Department of Statistics and Mathematics,
  Nr. 27, December 2005 Department of Statistics and Mathematics,
  Wien, Wirtschaftsuniv., 2005.
  \url{http://epub.wu-wien.ac.at/dyn/virlib/wp/showentry?ID=epub-wu-01_8cb}
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

%% \examples{}

\keyword{distribution}
\keyword{datagen}

