\name{unuran.discr-class}
\docType{class}

\alias{unuran.discr-class}
\alias{print,unuran.discr-method}
\alias{show,unuran.discr-method}
\alias{initialize,unuran.discr-method}

\title{Experimental class "unuran.discr"}

\description{%
  The class "unuran.discr" provides am interface to UNU.RAN objects for
  discrete distributions.
  The interface might be changed in future releases.
  \strong{Do not use unnamed arguments!}
}

\details{
  Create a new instance of an "unuran.discr" object using

  \code{new ("unuran.discr",pv=NULL, pmf=NULL, lb=0, ub=Inf)}.

  \describe{
    \item{pv}{probability vector. (numeric vector)}
    \item{pmf}{probability mass function. (\R function)}
    \item{lb, ub}{lower and upper bound of domain. (numeric, integers)}
  }

  The user is responsible that the given informations are consistent.
  It depends on the chosen method which information must be given / are
  used.
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

\examples{
## Create discrete distribution with given probability vector
pv <- c(1.,2.,1.5,0.,3.,1.2)     ## PV need not be normalized
dpv <- new("unuran.discr",pv=pv)
## Make generator (using method DGT)
unr <- unuran.new(dpv, "dgt")
## Draw sample of size 100
x <- unuran.sample(unr, 100)

## Create discrete distribution with given probability mass function
pmf <- function(x) dbinom(x,100,0.3)
dpmf <- new("unuran.discr",pmf=pmf,lb=0,ub=100)
## Make generator (using method DGT)
unr <- unuran.new(dpmf, "dgt")
## Draw sample of size 100
x <- unuran.sample(unr, 100)

}

\keyword{classes}
\keyword{distribution}
\keyword{datagen}
