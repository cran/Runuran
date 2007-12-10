\name{Runuran-package}
\docType{package}

\alias{Runuran-package}
\alias{Runuran}
\alias{unuran}

\title{
  Runuran -- R interface to Universal Non-Uniform RANdom variate
  generators library}

\description{
  R interface to the UNU.RAN library for Universal Non-Uniform RANdom
  variate generators}

\details{
  \tabular{ll}{
    Package: \tab Runuran       \cr
    Type:    \tab Package       \cr
    Version: \tab 0.6           \cr
    Date:    \tab 2007-09-01    \cr
    License: \tab GPL 2 or later\cr
  }
  This package provides an interface to the UNU.RAN library for universal
  non-uniform random number generators.

  It provides the following two main calls:
  \describe{
    \item{\code{unuran.new(dist,method)}}{
      Create a UNU.RAN object. For details see
      \code{\link{unuran.new}}.
    }
    \item{\code{unuran.sample(unuran.object,sample.size)}}{
      Get a random sample from the UNU.RAN object.
      For details see \code{\link{unuran.sample}}.
    }
  }
  It uses the \R built-in uniform random number generator.

  The \pkg{Runuran} package makes use of the string API of UNU.RAN,
  i.e., distributions and methods are described as strings and passed
  directly to UNU.RAN. The \pkg{Runuran} uses S4 classes to implemented
  this interface to UNU.RAN.
  Alternatively distributions can also be passed by
  means of S4 classes that contain the required information.
  See \code{\link{unuran.new}} for a short introduction to the interface. 

  Additionally there are some shortcuts that allow provide an easy to
  use interface: \code{\link{urdgt}} and \code{\link{urdau}} for
  discrete distributions with given probability vectors,
  and \code{\link{urtdr}} for log-concave distributions with given
  density function.

  To facilitate the use of the package for standard distributions
  the package also includes direct generation functions for
  approximately 30 standard distributions.

  For continuous univariate distributions these calls are:
  \code{\link{urbeta}}, 
  \code{\link{urburr}}, 
  \code{\link{urcauchy}}, 
  \code{\link{urchi}}, 
  \code{\link{urchisq}}, 
  \code{\link{urexp}}, 
  \code{\link{urextremeI}}, 
  \code{\link{urextremeII}}, 
  \code{\link{urf}}, 
  \code{\link{urgamma}},
  \code{\link{urgig}}, 
  \code{\link{urhyperbolic}}, 
  \code{\link{urlaplace}}, 
  \code{\link{urlnorm}}, 
  \code{\link{urlogis}}, 
  \code{\link{urlomax}}, 
  \code{\link{urnorm}}, 
  \code{\link{urpareto}}, 
  \code{\link{urplanck}}, 
  \code{\link{urpowerexp}}, 
  \code{\link{urrayleigh}}, 
  \code{\link{urt}}, 
  \code{\link{urtriang}}, and 
  \code{\link{urweibull}}.

  For discrete distributions these are:
  \code{\link{urbinom}}, 
  \code{\link{urgeom}},  
  \code{\link{urhyper}}, 
  \code{\link{urlogarithmic}}, 
  \code{\link{urnbinom}}, and
  \code{\link{urpois}}.

  All these calls have similar syntax to the analogous \R built-in
  generating functions (if these exist) but have an optional domain
  arguments \code{lb} and \code{ub},
  i.e., these calls also allow to draw samples from truncated
  distributions:

  \code{ur...(n, distribution parameters, lb , ub)}

  Compared to the corresponding \R functions these \code{ur} functions
  have a slightly different behavior.
  \itemize{
    \item{}{\code{ur} functions are often much faster (e.g. a factor) of
      about 5 for the \eqn{t} distribution when the sample size is
      large. For small samples the \code{ur} functions are slow.}
    \item{}{All \code{ur} functions allow to sample from truncated
      versions of the original distributions. Therefore the arguments
      \code{lb} (lower border) and \code{ub} (upper border) are
      available for all \code{ur} functions.}
    \item{}{Allmost all \code{ur} functions are based on fast numerical
      inversion algorithms. This is important for example for generating
      order statistics or random vectors from copulas.}
    \item{}{All \code{ur} functions do \strong{not} allow vectors as
      arguments (to be more precise: they only use the first element of
      the vector).}
  }
}

\note{
  The interface has been changed compared to the DSC 2003 paper.
} 

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

\references{
  J. Leydold and W. H\"ormann (2000-2007):
  UNU.RAN User Manual, see
  \url{http://statmath.wu-wien.ac.at/unuran/}.
  
  W.~H\"ormann, J.~Leydold, and G.~Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg

  G.~Tirler and J.~Leydold (2003):
  Automatic Nonuniform Random Variate Generation in \R.
  In: K.~Hornik and F.~Leisch, 
  Proceedings of the 3rd International Workshop on
  Distributed Statistical Computing (DSC~2003),
  March 20--22, Vienna, Austria.
}

\seealso{
  \code{\linkS4class{unuran}}, \code{\link{unuran.new}},
  \code{\link{unuran.sample}}, \code{\linkS4class{unuran.distr}}.
}

\keyword{distribution}
\keyword{datagen}

