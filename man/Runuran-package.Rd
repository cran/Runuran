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
    Version: \tab 0.10.0        \cr
    Date:    \tab 2009-06-22    \cr
    License: \tab GPL 2 or later\cr
  }
  \pkg{Runuran} provides an interface to the UNU.RAN library for
  universal non-uniform random number generators.
  It compiles three sets of functions of increasing power (and thus
  complexity):
  
  \describe{
    \item{}{\bold{[Distribution]} --
      Generators for paricular distributions. Their syntax is similar to
      the corresponding \R built-in functions.}
    \item{}{\bold{[Universal]} --
      Functions that offer an interface to a carefully selected
      collection of UNU.RAN methods with their most important
      parameters.}
    \item{}{\bold{[Advanced]} --
      Wrapper to the UNU.RAN string API. This gives access to all
      UNU.RAN methods and their variants.}
  }
  
  We have marked all functions in their corresponding help page by one
  these three tags.

  An introduction to \pkg{Runuran} with examples together with a very
  short survey on non-uniform random variate generation can be found in
  the package vignette (which can be displayed using
  \code{vignette("Runuran")}). 
}

\section{[Distribution]}{
  These functions have similar syntax to the analogous \R built-in
  generating functions (if these exist) but have an optional domain
  arguments \code{lb} and \code{ub}, i.e., these calls also allow to
  draw samples from truncated distributions:

  \code{ur...(n, distribution parameters, lb , ub)}

  Compared to the corresponding \R functions these \code{ur} functions
  have a different behavior:
  \itemize{
    \item{--}{\code{ur} functions are often much faster for large
      samples (e.g., a factor of about 5 for the \eqn{t} distribution).
      For small samples they are slow.}
    \item{--}{All \code{ur} functions allow to sample from truncated
      versions of the original distributions. Therefore the arguments
      \code{lb} (lower border) and \code{ub} (upper border) are
      available for all these functions.}
    \item{--}{Allmost all \code{ur} functions are based on fast numerical
      inversion algorithms. This is important for example for generating
      order statistics or random vectors from copulas.}
    \item{--}{All \code{ur} functions do \strong{not} allow vectors as
      arguments (to be more precise: they only use the first element of
      the vector).}
  }
  
  Currently generators for the following distributions are implemented.
  
  %% -- begin: list of distributions --
  Continuous Univariate Distributions (24):

  \tabular{lcl}{ 
    \emph{Function} \tab \tab \emph{Distribution} \cr
    \code{\link{urbeta}} \tab \ldots \tab Beta \cr
    \code{\link{urburr}} \tab \ldots \tab Burr \cr
    \code{\link{urcauchy}} \tab \ldots \tab Cauchy \cr
    \code{\link{urchi}} \tab \ldots \tab Chi \cr
    \code{\link{urchisq}} \tab \ldots \tab Chi-square \cr
    \code{\link{urexp}} \tab \ldots \tab Exponential \cr
    \code{\link{urextremeI}} \tab \ldots \tab Gumbel (extreme value type I) \cr
    \code{\link{urextremeII}} \tab \ldots \tab Frechet (extreme value type II) \cr
    \code{\link{urf}} \tab \ldots \tab F \cr
    \code{\link{urgamma}} \tab \ldots \tab Gamma \cr
    \code{\link{urgig}} \tab \ldots \tab GIG (generalized inverse Gaussian) \cr
    \code{\link{urhyperbolic}} \tab \ldots \tab Hyperbolic \cr
    \code{\link{urlaplace}} \tab \ldots \tab Laplace \cr
    \code{\link{urlnorm}} \tab \ldots \tab Log-Normal \cr
    \code{\link{urlogis}} \tab \ldots \tab Logistic \cr
    \code{\link{urlomax}} \tab \ldots \tab Lomax \cr
    \code{\link{urnorm}} \tab \ldots \tab Normal (Gaussian) \cr
    \code{\link{urpareto}} \tab \ldots \tab Pareto (of first kind) \cr
    \code{\link{urplanck}} \tab \ldots \tab Planck \cr
    \code{\link{urpowerexp}} \tab \ldots \tab Powerexponential (Subbotin) \cr
    \code{\link{urrayleigh}} \tab \ldots \tab Rayleigh \cr
    \code{\link{urt}} \tab \ldots \tab t (Student) \cr
    \code{\link{urtriang}} \tab \ldots \tab Triangular \cr
    \code{\link{urweibull}} \tab \ldots \tab Weibull \cr
  }

  Discrete Distributions (6):

  \tabular{lcl}{ 
    \emph{Function} \tab \tab \emph{Distribution} \cr
    \code{\link{urbinom}} \tab \ldots \tab Binomial \cr
    \code{\link{urgeom}} \tab \ldots \tab Geometric \cr
    \code{\link{urhyper}} \tab \ldots \tab Hypergeometric \cr
    \code{\link{urlogarithmic}} \tab \ldots \tab Logarithmic \cr
    \code{\link{urnbinom}} \tab \ldots \tab Negative Binomial \cr
    \code{\link{urpois}} \tab \ldots \tab Poisson \cr
  }

  %% -- end: list of distributions --
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
    \code{\link{ars.new}}   \tab \ldots \tab Adaptive Rejection Sampling \cr
    \code{\link{itdr.new}}  \tab \ldots \tab Inverse Transformed Density Rejection \cr
    \code{\link{pinv.new}}  \tab \ldots \tab Polynomial interpolation of INVerse CDF \cr
    \code{\link{srou.new}}  \tab \ldots \tab Simple Ratio-Of-Uniforms method \cr
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
  
\section{[Advanced]}{
  This interface provides the most flexible access to UNU.RAN.
  It requires three steps:
  
  \itemize{
    \item{1.}{
      Create a \code{\linkS4class{unuran.distr}} object that contains
      all required information about the \emph{target distribution}.
      We have three types of distribuions:}
  }
  \tabular{lcl}{ 
    \emph{Function}                \tab        \tab \emph{Type of distribution} \cr
    \code{\link{unuran.cont.new}}  \tab \ldots \tab continuous distributions \cr
    \code{\link{unuran.discr.new}} \tab \ldots \tab discrete distributions \cr
    \code{\link{unuran.cmv.new}}   \tab \ldots \tab multivariate continuous distributions \cr
  }
  \itemize{
    \item{2.}{
      Choose a \emph{generation method} and create a
      \code{\linkS4class{unuran}} object using function
      \code{\link{unuran.new}}.
      This function takes two argument: the distribution object created
      in Step 1, and a string that contains the chosen UNU.RAN method
      and (optional) some parameters to adjust this method to the given
      target distribution. We refer to the UNU.RAN for more details on
      this \dQuote{method string}.}
    \item{3.}{
      Use this object to \emph{draw samples} from the target
      distribution using \code{\link{ur}} or \code{\link{uq}}.}
  }
  \tabular{lcl}{ 
    \emph{Function}              \tab        \tab \cr
    \code{\link{ur}}             \tab \ldots \tab draw sample \cr
    \code{\link{uq}}             \tab \ldots \tab compute quantile (inverse CDF) \cr
    \code{\link{unuran.details}} \tab \ldots \tab show \code{unuran} object \cr
  }
}

\section{Uniform random numbers}{
  All UNU.RAN methods use the \R built-in random number generator as
  source of (pseudo-) random numbers. Thus the generated samples depend
  on the state \code{\link{.Random.seed}} and can be controlled by the
  \R\ functions \code{\link{RNGkind}} and \code{\link{set.seed}}.
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

\note{
  The interface has been changed compared to the DSC 2003 paper.
} 

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

\references{
  J. Leydold and W. H\"ormann (2000-2008):
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
  All objects are implemented as respective S4 classes
  \code{\linkS4class{unuran}},
  \code{\linkS4class{unuran.distr}},
  \code{\linkS4class{unuran.cont}},
  \code{\linkS4class{unuran.discr}},
  \code{\linkS4class{unuran}}.
}

\keyword{distribution}
\keyword{datagen}

