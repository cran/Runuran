\name{unuran-class}
\docType{class}

\alias{unuran-class}
\alias{print,unuran-method}
\alias{show,unuran-method}
\alias{initialize,unuran-method}

\title{Class "unuran" -- Universal Non-Uniform RANdom variate generators}

\description{%
  The class "unuran" provides am interface to the UNU.RAN
  library for universal non-uniform random number generators.
  It uses the \R built-in uniform random number generator.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("unuran", distribution, method)}.

  \itemize{
    \item{distribution:}{A character string that describes the target
      distribution, see UNU.RAN User Manual.}
    \item{method:}{A character string that describes the chosen
      generation method, see UNU.RAN User Manual. If omitted method
      \code{"auto"} (automatic) is used.}
  }

  See \code{\link{unuran.new}} for short introduction and examples for
  this interface.
}

\section{Methods}{
  The class "unuran" provides the following methods for
  handling objects:
  
  \describe{
    \item{unuran.sample}{\code{signature(object = "unuran")}: 
      Get a random sample from the stream object.}
    \item{r}{\code{signature(object = "unuran")}:
      Same as \code{unuran.sample}.}
    \item{initialize}{\code{signature(.Object = "unuran")}:
      Initialize unuran object. (For Internal usage only).}
    \item{print}{\code{signature(x = "unuran")}:
      Print info about unuran object.}
    \item{show}{\code{signature(x = "unuran")}:
      Same as \code{print}.}
  }
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

\note{
  The interface has been changed compared to the DSC 2003 paper.
} 

\section{Warning }{
  "unuran" objects cannot be saved and restored in later \R sessions.
  They must be newly created!
}

\references{
  J. Leydold and W. H\"ormann (2000-2007):
  UNU.RAN User Manual, see
  \url{http://statmath.wu-wien.ac.at/unuran/}.
  
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg

  G. Tirler and J. Leydold (2003):
  Automatic Nonuniform Random Variate Generation in \R.
  In: K.~Hornik and F.~Leisch, 
  Proceedings of the 3rd International Workshop on
  Distributed Statistical Computing (DSC~2003),
  March 20--22, Vienna, Austria.
}

\seealso{
  \code{\link{unuran.new}} and 
  \code{\link{unuran.sample}}.
}

\examples{
## create a new unuran object
unr <- new("unuran","normal","auto")

## alternatively one can use
unr <- unuran.new("normal","auto")

## get random number
unuran.sample(unr)

## get a random sample of size 10
unuran.sample(unr,10)

## use command r
r(unr)

}

\keyword{classes}
\keyword{distribution}
\keyword{datagen}
