
#############################################################################
##                                                                         ##
##   Tests for class 'unuran'                                              ##
##                                                                         ##
#############################################################################

## --- Test Parameters ------------------------------------------------------

## size of sample for test
samplesize <- 1.e7

## --- Load library ---------------------------------------------------------

library(Runuran)

## --------------------------------------------------------------------------

gu <- pinv.new(dnorm,lb=-Inf,ub=Inf)
gp <- pinv.new(dnorm,lb=-Inf,ub=Inf)
unuran.packed(gp) <- TRUE

u <- (-1:(samplesize+1))/samplesize
xu <- uq(gu,u)
xp <- uq(gp,u)

if (! identical(xu,xp))
  stop("packed and unpacked version of PINV differ !")

rm(gu)
rm(gp)

## --------------------------------------------------------------------------

gu <- pinv.new(dnorm,lb=0,ub=Inf)
gp <- pinv.new(dnorm,lb=0,ub=Inf)
unuran.packed(gp) <- TRUE

u <- (-1:(samplesize+1))/samplesize
xu <- uq(gu,u)
xp <- uq(gp,u)

if (! identical(xu,xp))
  stop("packed and unpacked version of PINV differ !")

rm(gu)
rm(gp)

## --------------------------------------------------------------------------

gu <- pinv.new(dnorm,lb=0,ub=1)
gp <- pinv.new(dnorm,lb=0,ub=1)
unuran.packed(gp) <- TRUE

u <- (-1:(samplesize+1))/samplesize
xu <- uq(gu,u)
xp <- uq(gp,u)

if (! identical(xu,xp))
  stop("packed and unpacked version of PINV differ !")

rm(gu)
rm(gp)

## --- End ------------------------------------------------------------------

silent <- gc()
detach("package:Runuran",unload = TRUE)

## --------------------------------------------------------------------------
