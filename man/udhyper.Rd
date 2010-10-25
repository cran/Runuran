\name{udhyper}
\alias{udhyper}

\title{UNU.RAN object for Hypergeometric distribution}

\description{
  Create UNU.RAN object for a Hypergeometric distribution with
  parameters \code{m}, \code{n}, and \code{k}.
  
  [Distribution] -- Hypergeometric.  %% Discrete
}

\usage{
udhyper(m, n, k, lb=max(0,k-n), ub=min(k,m))
}

\arguments{
  \item{m}{the number of white balls in the urn.}
  \item{n}{the number of black balls in the urn.}
  \item{k}{the number of balls drawn from the urn.}
  \item{lb}{lower bound of (truncated) distribution.}
  \item{ub}{upper bound of (truncated) distribution.}
}

\details{
  The Hypergeometric distribution is used for sampling \emph{without}
  replacement.  The density of this distribution with parameters
  \code{m}, \code{n} and \code{k} (named \eqn{Np}, \eqn{N-Np}, and
  \eqn{n}, respectively in the reference below) is given by
  \deqn{
    p(x) = \left. {m \choose x}{n \choose k-x} \right/ {m+n \choose k}
  }{p(x) =      choose(m, x) choose(n, k-x) / choose(m+n, k)}
  for \eqn{x = 0, \ldots, k}{x = 0, \ldots, k}.

  The domain of the distribution can be truncated to the 
  interval (\code{lb},\code{ub}).
}

\value{
  An object of class \code{"unuran.discr"}.
}

\references{
  N.L. Johnson, S. Kotz, and A.W. Kemp (1992):
  Univariate Discrete Distributions.
  2nd edition, John Wiley & Sons, Inc., New York.
  Chap. 6, p. 237.
}

\author{
  Josef Leydold and Wolfgang H\"ormann
  \email{unuran@statmath.wu.ac.at}.
}

\examples{
## Create distribution object for Hypergeometric distribution
dist <- udhyper(m=15,n=5,k=7)
## Generate generator object; use method DGT (inversion)
gen <- dgtd.new(dist)
## Draw a sample of size 100
x <- ur(gen,100)

}

\keyword{distribution}
