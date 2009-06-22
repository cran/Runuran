\name{unuran-class}
\docType{class}

\alias{unuran-class}
\alias{print,unuran-method}
\alias{show,unuran-method}
\alias{initialize,unuran-method}

\title{Class "unuran" -- Universal Non-Uniform RANdom variate generators}

\description{%
  The class \code{unuran} provides an interface to the UNU.RAN
  library for universal non-uniform random number generators.
  It uses the \R built-in uniform random number generator.
  
  [Advanced] -- UNU.RAN generator object.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("unuran", distribution, method)}.

  \itemize{
    \item{distribution:}{A character string that describes the target
      distribution (see UNU.RAN User Manual) or one of the S4 classes
      \code{\linkS4class{unuran.cont}},
      \code{\linkS4class{unuran.discr}}, or 
      \code{\linkS4class{unuran.cmv}}
      that holds information about the distribution.}
    \item{method:}{A character string that describes the chosen
      generation method, see UNU.RAN User Manual. If omitted method
      \code{"auto"} (automatic) is used.}
  }
  See \code{\link{unuran.new}} for short introduction and examples for
  this interface.
}

\section{Methods}{
  The class \code{\linkS4class{unuran}} provides the following methods
  for handling objects:
  
  \describe{
    \item{ur}{\code{signature(object = "unuran")}: 
      Get a random sample from the stream object.}
    \item{r}{\code{signature(object = "unuran")}:
      Same as \code{ur}.}
    \item{initialize}{\code{signature(.Object = "unuran")}:
      Initialize \code{unuran} object. (For Internal usage only).}
    \item{print}{\code{signature(x = "unuran")}:
      Print info about \code{unuran} object.}
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

\section{Warning}{
  \code{unuran} objects cannot be saved and restored in later \R
  sessions, nor is it possible to copy such objects to different nodes
  in a computer cluster.

  However, \code{unuran} objects for \emph{some} generation methods can be
  \dQuote{packed}, see \code{\link{unuran.packed}}.
  Then these objects can be handled like any other \R object
  (and thus saved and restored).

  All other objects \strong{must} be \strong{newly created} in a new \R
  session!
  (Using a restored object does not work as the \code{unuran} is then
  broken.)
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
  \code{\link{unuran.new}} and \code{\link{ur}} for faster creation and
  sampling routines, \code{\link{unuran.details}} for a more verbose
  version of \code{show}.

  \code{\link{unuran.packed}} allows to pack \emph{some} \code{unuran}
  objects.
  
  For distribution objects see \code{\linkS4class{unuran.cont}},
  \code{\linkS4class{unuran.discr}}, and \code{\linkS4class{unuran.cmv}}.
}

%% \examples{}

\keyword{classes}
\keyword{distribution}
\keyword{datagen}
