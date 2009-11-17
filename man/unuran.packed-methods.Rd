\name{unuran.packed-method}
\docType{methods}

\alias{unuran.packed}
\alias{unuran.packed-method}
\alias{unuran.packed,unuran-method}
\alias{unuran.packed<-}
\alias{unuran.packed<--methods}
\alias{unuran.packed<-,unuran-method}

\title{Pack "unuran" object}

\description{
  Pack \code{unuran} object in package \pkg{Runuran} and report its
  status (packed/unpacked).

  Packed \code{unuran} objects can be saved and loaded or sent to other
  nodes in a computer cluster (which is not possible for unpacked
  object).

  FIXME
}

\usage{
\S4method{unuran.packed}{unuran}(unr)
unuran.packed(unr) <- value
}

\arguments{
  \item{unr}{a \code{unuran} object.}
  \item{value}{\code{TRUE} to pack the object.}
}

\details{
  A \code{unuran} object contains a pointer to an external object in
  library UNU.RAN. Thus it cannot be saved and restored in later \R
  sessions, nor is it possible to copy such an object to different nodes
  in a computer cluster.

  By \dQuote{packing} an \code{unuran} object all required data are
  copied from the external object into an \R list and stored in the
  \code{unuran} object while the external UNU.RAN object is
  destroyed. Thus the object can be handled like any other \R
  object. Moreover, it can be still used as argument for
  \code{\link{ur}} and \code{\link{uq}} (which may have even faster
  execution times then).
  Packed \code{unuran} objects cannot be unpacked any more.

  Notice that currently only objects that implement method \sQuote{PINV}
  can be packed.
}

\section{Methods}{
  Currently only objects that implement method \sQuote{PINV} can be
  packed.
}

\note{
  Note that due to limitations of floating point arithmetic the output
  of a \code{\link{uq}} call with the same input value for \code{u} may
  slightly differ for the packed and unpacked version.
}
  
\seealso{%
  \code{\linkS4class{unuran}}, \code{\link{pinv.new}}.
}

\examples{
## create a unuran object for half-normal distribution using methed 'PINV'
gen <- pinv.new(dnorm,lb=0,ub=Inf)

## status of object 
unuran.packed(gen)

## draw a random sample of size 10
x <- ur(gen,10)

## pack unuran object
unuran.packed(gen) <- TRUE
unuran.packed(gen)

## draw a random sample of size 10
x <- ur(gen,10)

\dontrun{
## unpacking is not supported
unuran.packed(gen) <- FALSE    ## results in error 
}
}

\keyword{methods}
\keyword{distribution}
\keyword{datagen}
