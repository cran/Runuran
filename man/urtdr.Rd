\name{urtdr}
\docType{function}

\alias{urtdr}

\title{UNU.RAN transformend density rejection generator}

\description{
  UNU.RAN random variate generator for continuous distributions with
  given probability density function (PDF).
  It is based on the Transformed Density Rejection method (TDR).

  [ \strong{Deprecated}.
  Use \code{\link{tdr.new}} instead. ]
}

\usage{
urtdr(n, pdf, dpdf=NULL, lb=-Inf, ub=Inf, islog=TRUE, ...)
}
\arguments{
  \item{n}{size of required sample.}
  \item{pdf}{probability density function. (\R function)}
  \item{dpdf}{derivative of the \code{pdf}. (\R function)}
  \item{islog}{whether the \code{pdf} is given by tits logarithm (the
    \code{dpdf} is then the derivative of the logarithm). (boolean)}
  \item{lb, ub}{lower and upper bound of domain. (numeric)}
  \item{...}{(optional) arguments for \code{pdf}}
}

\details{
  This method is \strong{deprecated}. Please use \code{\link{tdr.new}}
  together with \code{\link{ur}} instead.
  
  For further details see \code{\link{tdr.new}}.
}

\seealso{
  \code{\link{tdr.new}} and \code{\link{ur}} for replacement.
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

%% \examples{}

\keyword{distribution}
\keyword{datagen}

