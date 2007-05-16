\name{unuran.new}
\docType{function}

\alias{unuran.new}

\title{Create a UNU.RAN object}

\description{
  Create a new UNU.RAN object in package \pkg{Runuran} that can be used
  for sampling from the specified distribution.
  The function \code{unuran.sample} can then be used to draw a random
  sample.
}

\usage{
unuran.new(distr,method="auto")
}
\arguments{
  \item{distr}{a string or an S4 class describing the distribution.}
  \item{method}{a string describing the random variate generation
    method.}
}

\details{
  The UNU.RAN package works in the following way:
  \describe{
    \item{Initialization Step:}{
      Create an UNU.RAN object where the distribution object and the
      method are given independendly. This step is performed by a
      \code{unuran.new} call. Thus it is possible to choose a method
      that corresponds best to ones application. However, not all methods
      work in combination with all distribution. Thus the result of
      \code{unuran.new} can be empty and must be checked.
      Both the distribution and the method are described by character
      strings, see below.}
    \item{Generation Step:}{Once a UNU.RAN object has been created it can
      be used to draw random samples from the target distribution using
      \code{\link{unuran.sample}}.
    }
  }
}

\note{
  \code{unuran.new(distr,method)} is an alias for 
  \code{new("unuran", distr, method)}.
}

\section{The Distribution String}{
  The argument \code{distr} describes the target distribution from which
  wants to draw one (or more) random sample(s).
  It is passed through to the UNU.RAN interface, see the UNU.RAN manual
  for a detailed description at
  \url{http://statmath.wu-wien.ac.at/unuran/doc/unuran.html\#StringDistr}.
  Notice that the key-value pairs have to be separated by semicolons
  \code{;}.
  
  \describe{
    \item{Standard distributions:}{
      These can be given by a string with the parameters as a
      comma-separated list in parenthesis, e.g.,
      \code{"normal(1,2)"} for the normal distribution with mean 1 and
      standard deviation 2. We remark, however, that we have prepared
      wrapper functions with function names starting with \code{ur} that
      combine the initialization and sampling steps for theses
      distributions. See the man pages and source of these functions for
      further reference. It is possible to truncate standard
      distributions with the \code{domain} keyword, e.g. to sample from
      the half-normal distribution one has to use
      \code{"normal(); domain=(0,Inf)"}.
    }

    \item{Continuous univariate distributions:}{
      The string has to start with \code{"cont"} to indicate this type
      of distribution. Then the probability density function
      (\code{pdf}) or cumulative distribution function (\code{cdf}) must
      be given. Moreover, the domain must be set as otherwise
      \code{(-Inf,Inf)} is assumed. 
      E.g., if the target distribution has a density proportional to
      \eqn{1/x^2} on \eqn{x\ge 1}{x >= 1} we use
      \code{"cont; pdf='1/x^2'; domain=(1,Inf)"}.
      Notice the single quotes \code{'} that encapsulates the function
      string. Note that the given PDF must be a non-negative function with
      finite integral. It need not integrate to 1, i.e., it can be any
      multiple of density function. The syntax of the function string
      (\code{"1/x^2"} in our example) is very similar to the \R syntax;
      see the UNU.RAN manual for detailed description at
      \url{http://statmath.wu-wien.ac.at/unuran/doc/unuran.html\#StringFunct}.

      Notice, that there also exist an indicator function for densities.
      For example, to get a triangular distribution on \eqn{(-1,1)} with
      mode at \eqn{0} we write
      \code{"cont; pdf='(x<0)*(1+x) + (x<=0)*(1-x)'; domain=(-1,1)"}.

      The CDF can be specified analogously:
      \code{"cont; cdf='exp(-x)'; domain=(0,Inf)"}.
    }

    \item{Discrete distribution:}{
      The string has to start with \code{"discr"}. Then one either can
      provide a (multiple of a) probability mass function
      (\dQuote{density}) using the keyword \code{pmf} or a (multiple of
      a) probability vector using \code{pv}.

      Here is an example for a discrete distribution on \eqn{(1,100)}
      with probability vector proportional to \eqn{1/x}:
      \code{"discr; pmf='1/x'; domain=(1,100)"}.

      A discrete distribution with a probability vector proportional to
      \eqn{(1,2,3,4,5,4,3,2,1)} is given by
      \code{"discr; pv=(1,2,3,4,5,4,3,2,1)"} (there is no necessity for
      providing a domain except it does not start at \eqn{0}).

      Notice that for the latter case there exist wrapper functions
      \code{\link{urdgt}} and \code{\link{urdau}} for drawing a random
      sample.
    }

    \item{Domain:}{
      Domains can be set using the \code{domain}
      keyword. E.g., when the domain should be set to
      \eqn{(-1,\infty)}{(-1,Inf)}
      we write \code{domain=(-1,Inf)}.
      If non-integer domains are given for discrete distribution, then
      the given numbers are simply truncated to integers.
    }

    \item{Mode:}{
      For some sampling methods the location of the mode of a
      distribution can be helpful. It can be given using the \code{mode}
      keyword: \code{"mode=1.5"}.
    }

    \item{Continuous multivariate distributions:}{
      The string interface for distributions currently does not support
      multivariate distributions. Please use the corresponding S4 class
      \code{"unuran.cmv"} instead.
    }
  }
} %% -- end "The Distribution String"

\section{The S4 Classes for Distributions}{
  The package \pkg{Runuran} provides an alternative approach to describe
  distributions. Thus the following S4 classes can be used.
  \describe{
    \item{class "unuran.cont"}{
      for univariate continuous distributions,
      see \code{\linkS4class{unuran.cont}}.
    }
    \item{class "unuran.discr"}{
      for discrete distributions,
      see \code{\linkS4class{unuran.discr}}.
    }
    \item{class "unuran.cmv"}{
      for multivariate continuous distributions,
      see \code{\linkS4class{unuran.cmv}}.
    }
  }
  This interface allows to use \R functions and object to define
  distributions.
} %% -- end "The S4 Classes for Distributions"

\section{The Method String}{
  The argument \code{method} describes generation method for drawing
  the random samples. UNU.RAN provides various such method for all
  types of distributions. This is indentional as not all methods work
  for all distributions. Morever, the choice of the generation method
  also might depend on the particular application (stochastic
  simulation problem) for which the sample has to be drawn.
  We give her short list of methods that are most useful.
  for a detailed description see the UNU.RAN manual at
  \url{http://statmath.wu-wien.ac.at/unuran/doc/unuran.html\#StringMethod}.
  
  \describe{
    \item{\code{"AUTO"}:}{The automatic method is obviously the
      simplest one (but is still experimental). It can be used for any
      type of distribution and tries to find an appropriate method among
      all that are available.
      This is the default method when no \code{method} string is
      provided by the user.
    }
  }

  All other methods only work for continuous univariate, discrete, or
  continuous multivariate distributions.
  (We remark, that UNU.RAN also provides methods for other types of
  distributions which are currently only partially supported by the
  \pkg{Runuran} interface.)
  All methods have optional parameters that allow to adjust the method
  for the target distribution. However, most of these are only useful
  in rare situations with distributions with extreme properties. See
  the UNU.RAN manual if you need a the full list of paramters.
  
  \emph{For continuous univariate distributions} the most important
  methods are:
  
  \describe{
    \item{\code{"HINV"}:}{
      Fast numerical inversion using a piecewise polynomial
      interpolation of the inverse CDF (Hermite interpolation). It
      requires the CDF of the distritbution.
      
      The accuracy of the method can be controlled by the parameter
      \code{u_resolution} which gives the 
      maximal tolerated \eqn{u}-error (i.e., \eqn{|U-CDF(X)|})
      of the numeric inversion (Default is \code{1.e-10}).
      Use \code{"hinv; u_resolution=1.e-12"} to use this method with
      increased accuracy.
    }
    
    \item{\code{"TDR"}:}{
      Transformed density rejection. It requires the density of the
      distribution and currently works only for
      \eqn{T_{-0.5}}{T_(-0.5)}-concave distributions (i.e., unimodal
      densities with tails not higher than \eqn{(1/x^2)}). It 
      decomposes the density into subintervals and uses a seperate hat
      proportional to \eqn{1/(c+x^2)} for each subinterval; the number
      of subintervals is selected such that the acceptance probability
      is about 99\%.
      
      It is possible to choose a piecewise exponential hat by means of
      parameter \code{c}: \code{"TDR; c=0"}. It then only works for
      log-concave distributions.
      
      The acceptance probability can be changed using parameter
      \code{"max_sqhratio"}. For example if an acceptance probability
      of 90\% is acceptable (which results in a faster setup) use
      \code{"TDR; max_sqhratio=0.95"}.
      
      The marginal generation speed can be increased using
      parameter \code{"variant_ia"} which cases TDR to use a techique
      called \dQuote{immediate acceptance}:
      \code{"TDR; variant_ia"}.
    }
  }

  \emph{For discrete univariate distributions} the most important
  methods are:
  
  \describe{
    \item{\code{"DGT"}:}{
      Guide-Table method for discrete inversion. Works for
      distributions with given probability vectors and for
      distributions with given PMF when the number of non-negligible
      probabilities is not too large.
    }
    \item{\code{"DAU"}:}{
      Alias-Urn method. Works for the same class of distributions as
      method \code{"DGT"}.
    }
    \item{\code{"DARI"}:}{
      Automatic Rejection Inversion. It is based on discrete
      rejection and requires the PMF of an  
      unimodal discrete distributions with tails not higher than
      \eqn{1/x^2}. Optionally the mode could be provided.
    }
  }

  \emph{For continuous multivariate distributions} the most important
  methods are:
  \describe{
    \item{\code{"HITRO"}:}{
      An easy to use Markov chain sampling algorithms based on
      the Hit-and-Run sample in combination with the ratio-of-Uniforms
      methods. It requires the density of the distribution.

      Notice that MCMC samplers do not produce a sequence of independent
      variates. The drawn sample follows the target distribution only
      appropriately.

      \strong{Beware: MCMC sampling can be dangerous!}
      
      Method HITRO provides two parameters that increase the quality the
      generated sample: \code{"thinning"} and \code{"burnin"}.
      By the first parameter only a fraction of the generated sequence
      is returned while the second one discards the first part of the
      sequence (\dQuote{burn-in phase}):
      \code{"HITRO; thinning=3; burnin=1000"}.

      There are some more parameter to control the behavior of the
      algorithm. 
      See the UNU.RAN manual and \code{\link{urhitro}} for more
      informations.
    }
    \item{\code{"VNROU"}:}{
      A naiv implementation of the multivariate Ratio-of-Uniforms
      method. It uses rejection from the (minimal) bounded rectangle.
      It requires the density function (or any multiple of it) which
      should be log-concave. Otherwise the computation of the bounding
      rectangle might fail and should be provided by the user.
      This method only works for a small number of dimensions (up to
      10). The acceptance rate (and the thus the marginal generation
      times) can be prohibitive small when the number of dimensions are
      too high or the shape the density function is not well-suited.
    }
  }
} %% -- end "The Method String"

\seealso{
  \code{\link{runif}} and \code{\link{.Random.seed}} about random number
  generation, \code{\linkS4class{unuran}} for the UNU.RAN class,
  \code{\link{unuran.sample}}, \code{\link{urtdr}}, \code{\link{urdgt}}
  and \code{\link{urdau}}.

  For distribution objects see \code{\linkS4class{unuran.cont}},
  \code{\linkS4class{unuran.discr}}, and
  \code{\linkS4class{unuran.cmv}}.
}

\references{
  W. H\"ormann, J. Leydold, and G. Derflinger (2004):
  Automatic Nonuniform Random Variate Generation.
  Springer-Verlag, Berlin Heidelberg.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu-wien.ac.at}.
}

\examples{
## Make histograms for samples of size 10000 for some non-standard distributions

unr <- new("unuran", "cont; pdf='log(x)*(x<1.3098)+(x>1.3098)*exp(-x)'; domain = (1, Inf)", "TDR");
x <- unuran.sample(unr,10000);
hist(x)

unr <- new("unuran", "distr=cont; pdf='1/x^3'; domain = (0.5, Inf)", "TDR");
x <- unuran.sample(unr,10000);
hist(x)

unr <- new("unuran", "distr=discr; pmf='1/(2+x)^3'; domain = (1, Inf);pmfsum=0.05","DARI");
x <- unuran.sample(unr,10000);
hist(x)
}

\keyword{distribution}
\keyword{datagen}

