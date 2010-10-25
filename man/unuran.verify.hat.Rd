\name{unuran.verify.hat}
\alias{unuran.verify.hat}

\title{Verify hat and squeezes in a "unuran" generator object}

\description{
  Verify hat function and squeezes in a given \code{unuran} generator
  that implements a rejection method.
  
  [Advanced] -- Verify rejection method.
}

\usage{
unuran.verify.hat(unr, n=1e5, show=TRUE)
}

\arguments{
  \item{unr}{an \code{unuran} object.}
  \item{n}{sample size. (integer)}
  \item{show}{whether the result is printed on the console. (boolean)}
}

\details{
  UNU.RAN is a library of so called black-box algorithms. For algorithms
  based on the rejection method this means that hat and squeezes are
  created automatically during the setup. Obviously not all algorithms
  work for all distribution. Then usually the setup fails (which is
  good, since then one does not silently obtain a random sample from a
  distribution other then the requested.) 

  Although we have tested these algorithms with a lot of distributions
  (including those with extreme properties) there is still some (minor)
  chance that hat and squeezes are computed without any warnings, but
  are incorrect, i.e., the inequalities
  \deqn{squeeze(x) \le density(x) \le hat(x)}
  are not satisfied for all \eqn{x}.
  This might happen due to serious round-off errors for densities with
  extreme properties (e.g., sharp and narrow peaks). But it also might
  be caused by some incorrect additional information about the
  distribution given by the user which has not been detected by various
  checks during the setup. If one is unsure about his or her chosen
  generation method one can check these inequalities.

  Routine \code{unuran.verify.hat} allows to run generator \code{unr}
  and check whether the two inequalities are violated. This is done for 
  every point \eqn{x} that is sampled from the hat distribution. This
  includes also those points that are rejected.
  The function counts the occurrences of such evaluations and 
  returns the ratio of this number and the sample size \code{n}.
  (It is thus a little bit too high since the total number of generated
  but rejected points is not known.)
  Yet, it does not provide any information about the magnitude of
  violation of the inequality.

  If \code{show} is \code{TRUE} then this routine prints this ratio and
  some diagnostcs to the console.

  Routine \code{unuran.verify.hat} does not work for algorithms that do
  not implement a rejection method.
}

\value{
  Ratio of number occurrences where the hat and squeezes violate the
  inequality and the sample size.
}

\note{
  Due to round-off errors there might exist a few points where the ratio
  \eqn{density(x) / hat(x)} is slightly larger than 1.
  In our experiments we observed a few cases where this ratio was as
  large as \eqn{1+10^{-8}}{1+1.e-8} for some points although we could
  proof (using real numbers instead of floating point numbers) that hat
  and squeeze are computed correctly.

  On the other hand, there are cases where, due to the limitation of
  floating point arithmetic, it is not possible to sample from the
  target distribution at all. 
  The Gamma distribution with extremely small shape parameter, say
  0.0001, is such an example. Then the continuous Gamma distribution
  degenerates to a point distribution with only a few points with
  significant mass.
}

\seealso{%
  \code{\linkS4class{unuran}}.
}


\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create a generator object that implements a rejection method
unr <- tdrd.new(udnorm())

## Verify hat and squeeze
unuran.verify.hat(unr)

}

\keyword{distribution}
\keyword{datagen}

