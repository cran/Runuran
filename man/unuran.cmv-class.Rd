\name{unuran.cmv-class}
\docType{class}

\alias{unuran.cmv-class}
\alias{print,unuran.cmv-method}
\alias{show,unuran.cmv-method}
\alias{initialize,unuran.cmv-method}

\title{Experimental class "unuran.cmv"}

\description{%
  The class "unuran.cmv" provides am interface to UNU.RAN objects for
  continuous multivariate distributions.
  The interface might be changed in future releases.
  \strong{Do not use unnamed arguments!}
}

\details{
  Create a new instance of an "unuran.cmv" object using

  \code{new ("unuran.cmv", dim=1, pdf=NULL, mode=NULL, center=NULL, ll=NULL, ur=NULL)}.

  \describe{
    \item{dim}{number of dimensions of the distribution. (integer)}
    \item{pdf}{probability density function. (\R function)}
    \item{mode}{location of the mode. (numeric vector -- optional)}
    \item{center}{point in \dQuote{typical} region of distribution,
      e.g. the approximate location of the mode. It is used by several
      methods to locate the main region of the distribution.
      If omitted the \code{mode} is implicitly used. If the \code{mode}
      is not given either, the origin is used. (numeric vector --
      optional)}
    \item{ll,ur}{lower left and upper right vertex of a rectangular
      domain of the \code{pdf}. The domain is only set if both vertices
      are not \code{NULL}. Otherwise, the domain is unbounded by default.
      (numeric vectors)}
  }

  The user is responsible that the given informations are consistent.
  It depends on the chosen method which information must be given / are
  used. It is important, that the mode is contained in the (closure of
  the) domain.
} 

\examples{
## Create discrete distribution with given probability density function
mvpdf <- function (x) { exp(-sum(x^2)) }
mvdist <- new("unuran.cmv", dim=2, pdf=mvpdf)
## Make generator (using method HITRO)
mvunr <- unuran.new(mvdist, "hitro")
## Draw sample of size 100
x <- unuran.sample(mvunr, 100)
## make scatter plot of random sample
plot(x)

## now we draw a sample from the same density but restricted to
## the rectangle [0,1]x[0,1]
mvdist <- new("unuran.cmv", dim=2, pdf=mvpdf, mode=c(0,0), ll=c(0,0), ur=c(1,1))
mvunr <- unuran.new(mvdist, "hitro")
x <- unuran.sample(mvunr, 100)
plot(x)

## draw a sample from the same density but restricted to
## the the first quadrant
mvdist <- new("unuran.cmv", dim=2, pdf=mvpdf, mode=c(0,0), ll=c(0,0), ur=c(Inf,Inf))
mvunr <- unuran.new(mvdist, "hitro")
x <- unuran.sample(mvunr, 100)
plot(x)

## draw a sample from the same density but restricted to
## the rectangle [1,2]x[1,2]
## (Don't forget to provide a starting point for the HITRO method by
## setting 'center'
mvdist <- new("unuran.cmv", dim=2, pdf=mvpdf, center=c(1.1,1.1), ll=c(1,1), ur=c(2,2))
mvunr <- unuran.new(mvdist, "hitro")
x <- unuran.sample(mvunr, 100)
plot(x)

}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

\references{
  J. Leydold and W. H\"ormann (2000-2007):
  UNU.RAN User Manual, see
  \url{http://statmath.wu-wien.ac.at/unuran/}.
}

\seealso{
  \code{\linkS4class{unuran}}, \code{\link{unuran.new}}.
}

\keyword{classes}
\keyword{distribution}
\keyword{datagen}
