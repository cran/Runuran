\name{unuran.details}
\alias{unuran.details}

\title{Information on a given "unuran" generator object}

\description{
  Prints type of \code{unuran} generator, data used from distribution, 
  parameter for algorithm, performance characteristic, and hints to
  adjust the performance of the generator.
  It also returns a list that contains some of these data.
  
  [Advanced] -- Print object.
}

\usage{
unuran.details(unr, show=TRUE, return.list=FALSE, debug=FALSE)
}

\arguments{
  \item{unr}{an \code{unuran} object.}
  \item{show}{whether the data are printed on the console. (boolean)}
  \item{return.list}{whether some of the data are returned in a
    list. (boolean)}
  \item{debug}{if TRUE, store additional data in returned list. This
    might be useful to examine a method. (boolean)}
}

\details{
  If \code{show} is \code{TRUE} then this routine prints data about the
  generator object to the console.

  If \code{return.list} is \code{TRUE} then a list that contains some of
  these data is returned. This an experimental feature and components of
  the list may be extended on request.

  The components of the returned list depend on the particular method.
  However, the following are common to all objects:
  
  \describe{
    \item{\code{method}}{string that contains the name of the generation
      method.}
    \item{\code{type}}{one of the following strings that describes the
    type of the generation method:
      \describe{
	\item{\code{"inv"}}{inversion method}
	\item{\code{"ar"}}{acceptance-rejection method}
	\item{\code{"iar"}}{acceptance-rejection whether inversion is
	  used for the proposal distribution}
	\item{\code{"mcmc"}}{Markov chain Monte Carlo sampler}
	\item{\code{"other"}}{none of the above methods}
      }
    }
  }

  In addition the following components may be available:
  \describe{
    \item{\code{area.pdf}}{area below density function of the
      distribution.}
    \item{\code{area.hat}}{area below hat function for an 
      acceptance-rejection method.}
    \item{\code{rejection.constant}}{rejection constant for an
      acceptance-rejection method.
      It given as the ratio \code{area.hat / area.pdf}.
    }
    \item{\code{area.squeeze}}{area below squeeze function for an 
      acceptance-rejection method.
      \code{area.hat / area.squeeze} can be used as upper bound for the
      rejection constant.}
    \item{\code{intervals}}{integer that contains the number of
      subintervals into which the domain of the target distribution is
      split for constructing a hat function / approximating function.}
    \item{\code{truncated.domain}}{vector of length 2 that contains
      upper and lower boundary of the \sQuote{computational domain} that
      is used for constructing an approximating function.}
  }
}

\seealso{%
  \code{\linkS4class{unuran}}.
}

%%\note{ }

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create a generator object
distr <- udnorm()
gen <- tdrd.new(distr)

## print data about object on console
unuran.details(gen)

## get list with some of these data
data <- unuran.details(gen,return.list=TRUE)

}

\keyword{distribution}
\keyword{datagen}

