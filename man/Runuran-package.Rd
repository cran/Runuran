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
    Version: \tab 0.32          \cr
    Date:    \tab 2020-10-28    \cr
    License: \tab GPL 2 or later\cr
  }
  \pkg{Runuran} provides an interface to the UNU.RAN library for
  universal non-uniform random number generators.
  It provides a collection of so called automatic methods for
  non-uniform random variate generation. Thus it is possible to draw
  samples from uncommon distributions. Nevertheless, (some of) these
  algorithms are also well suited for standard distribution like the
  normal distribution. Moreover, sampling from distributions like the
  generalized hyperbolic distribution is very fast. Such distributions
  became recently popular in financial engineering.

  \pkg{Runuran} compiles four sets of functions of increasing power (and
  thus complexity):
  
  \describe{
    \item{}{\bold{[Special Generator]} --
      Generators for particular distributions. Their syntax is similar to
      the corresponding \R built-in functions.}
    \item{}{\bold{[Universal]} --
      Functions that offer an interface to a carefully selected
      collection of UNU.RAN methods with their most important
      parameters.}
    \item{}{\bold{[Distribution]} --
      Functions that create objects for important distributions.
      These objects can then be used in combination with one of the
      universal methods which is best suited for a particular
      application.}
    \item{}{\bold{[Advanced]} --
      Wrapper to the UNU.RAN string API. This gives access to all
      UNU.RAN methods and their variants.}
  }
  
  We have marked all functions in their corresponding help page by one
  these four tags.

  An introduction to \pkg{Runuran} with examples together with a very
  short survey on non-uniform random variate generation can be found in
  the package vignette (which can be displayed using
  \code{vignette("Runuran")}). 
}

\section{[Special Generator]}{
  These functions have similar syntax to the analogous \R built-in
  generating functions (if these exist) but have optional domain
  arguments \code{lb} and \code{ub}, i.e., these calls also allow to
  draw samples from truncated distributions:

  \code{ur\ldots(n, distribution parameters, lb , ub)}

  Compared to the corresponding \R functions these \code{ur\ldots}
  functions have a different behavior:
  \itemize{
    \item \code{ur\ldots} functions are often much faster for large
      samples (e.g., a factor of about 5 for the \eqn{t} distribution).
      For small samples they are slow.
    \item All \code{ur\ldots} functions allow to sample from truncated
      versions of the original distributions. Therefore the arguments
      \code{lb} (lower border) and \code{ub} (upper border) are
      available for all these functions.
    \item Almost all \code{ur\ldots} functions are based on fast numerical
      inversion algorithms. This is important for example for generating
      order statistics, quasi-Monte Carlo methods or random vectors from
      copulas.
    \item All \code{ur\ldots} functions do \strong{not} allow vectors as
      arguments (to be more precise: they only use the first element of
      the vector).
  }

  However, we recommend to use the more flexible approach described in
  the next sections below.

  A list of all available special generators can be found in
  \code{\link{Runuran.special.generators}}.
}

\section{[Universal]}{
  These functions allow access to a selected collection of UNU.RAN
  methods. They require some data about the target distribution as
  arguments and return an instance of a UNU.RAN generator object that is
  implemented as an S4 class \code{\linkS4class{unuran}}. These can then
  be used to draw samples from the desired distribution by means of 
  function \code{\link{ur}}.
  Methods that implement an \emph{inversion} method can also be
  used for quantile function \code{\link{uq}}.
  
  Currently the following methods are available by such functions.
  
  Continuous Univariate Distributions:
  
  \tabular{lcl}{ 
    \emph{Function}         \tab        \tab \emph{Method} \cr
    \code{\link{arou.new}}  \tab \ldots \tab Automatic Ratio-of-Uniforms method \cr
    \code{\link{ars.new}}   \tab \ldots \tab Adaptive Rejection Sampling \cr
    \code{\link{itdr.new}}  \tab \ldots \tab Inverse Transformed Density Rejection \cr
    \code{\link{pinv.new}}  \tab \ldots \tab Polynomial interpolation of INVerse CDF \cr
    \code{\link{srou.new}}  \tab \ldots \tab Simple Ratio-Of-Uniforms method \cr
    \code{\link{tabl.new}}  \tab \ldots \tab TABLe based rejection \cr
    \code{\link{tdr.new}}   \tab \ldots \tab Transformed Density Rejection \cr
  }
  
  Discrete Distributions:
  
    \tabular{lcl}{ 
    \emph{Function}         \tab        \tab \emph{Method} \cr
    \code{\link{dari.new}}  \tab \ldots \tab Discrete Automatic Rejection Inversion \cr
    \code{\link{dau.new}}   \tab \ldots \tab Alias-Urn Method \cr
    \code{\link{dgt.new}}   \tab \ldots \tab Guide-Table Method for discrete inversion \cr
  }
  
  Multivariate Distributions:
  
  \tabular{lcl}{ 
    \emph{Function}         \tab        \tab \emph{Method} \cr
    \code{\link{hitro.new}} \tab \ldots \tab Hit-and-Run with Ratio-of-Uniforms method \cr
    \code{\link{vnrou.new}} \tab \ldots \tab Multivariate Naive Ratio-Of-Uniforms method \cr
  }
}
  
\section{[Distribution]}{
  Coding the required functions for particular distributions can be
  tedious. Thus we have compiled a set of functions that create UNU.RAN
  distribution objects that can directly be used with the functions from
  section [Universal].

  A list of all available distributions can be found in
  \code{\link{Runuran.distributions}}.
}

\section{[Advanced]}{
  This interface provides the most flexible access to UNU.RAN.
  It requires three steps:
  
  \enumerate{
    \item
    Create a \code{\linkS4class{unuran.distr}} object that contains
    all required information about the \emph{target distribution}.
    We have three types of distributions:
  
    \tabular{lcl}{ 
      \emph{Function}                \tab        \tab \emph{Type of distribution} \cr
      \code{\link{unuran.cont.new}}  \tab \ldots \tab continuous distributions \cr
      \code{\link{unuran.discr.new}} \tab \ldots \tab discrete distributions \cr
      \code{\link{unuran.cmv.new}}   \tab \ldots \tab multivariate continuous distributions \cr
    }

    The functions from section [Distribution] creates such objects for
    particular distributions.
    
    \item
    Choose a \emph{generation method} and create a
    \code{\linkS4class{unuran}} object using function
    \code{\link{unuran.new}}.
    This function takes two argument: the distribution object created
    in Step 1, and a string that contains the chosen UNU.RAN method
    and (optional) some parameters to adjust this method to the given
    target distribution. We refer to the UNU.RAN for more details on
    this \dQuote{method string}.

    \item
    Use this object to \emph{draw samples} from the target
    distribution using \code{\link{ur}} or \code{\link{uq}}.
  
    \tabular{lcl}{ 
      \emph{Function}              \tab        \tab \cr
      \code{\link{ur}}             \tab \ldots \tab draw sample \cr
      \code{\link{uq}}             \tab \ldots \tab compute quantile (inverse CDF) \cr
      \code{\link{unuran.details}} \tab \ldots \tab show \code{unuran} object \cr
    }
  }
}

\section{Density and distribution function}{
  UNU.RAN distribution objects and generator objects may also be used to
  compute density and distribution function for a given distribution by
  means of \code{\link{ud}} and \code{\link{up}}.
}

\section{Uniform random numbers}{
  All UNU.RAN methods use the \R built-in random number generator as
  source of (pseudo-) random numbers. Thus the generated samples depend
  on the state \code{\link{.Random.seed}} and can be controlled by the
  \R functions \code{\link{RNGkind}} and \code{\link{set.seed}}.
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
  (Using a restored object does not work as the \code{"unuran"} object
  is then broken.)
}

\note{
  The interface has been changed compared to the DSC 2003 paper.
} 

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\references{
  J. Leydold and W. H\"ormann (2000-2008):
  UNU.RAN User Manual, see
  \url{http://statmath.wu.ac.at/unuran/}.
  
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg

  G. Tirler and J. Leydold (2003):
  Automatic Nonuniform Random Variate Generation in \R.
  In: K. Hornik and F. Leisch, 
  Proceedings of the 3rd International Workshop on
  Distributed Statistical Computing (DSC 2003),
  March 20--22, Vienna, Austria.
}

\seealso{
  All objects are implemented as respective S4 classes
  \code{\linkS4class{unuran}},
  \code{\linkS4class{unuran.distr}},
  \code{\linkS4class{unuran.cont}},
  \code{\linkS4class{unuran.discr}},
  \code{\linkS4class{unuran}}.

  See \code{\link{Runuran.special.generators}} for an overview of
  special generators and 
  \code{\link{Runuran.distributions}} for a list of ready-to-use
  distributions suitable for the automatic methods.
}

\keyword{distribution}
\keyword{datagen}
