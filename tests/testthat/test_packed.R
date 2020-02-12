## --------------------------------------------------------------------------
##
## Check function Runuran.packed()
##
## --------------------------------------------------------------------------

## --- Test Parameters ------------------------------------------------------

## size of sample for test
samplesize <- 1.e6

## --------------------------------------------------------------------------

context("[packed] - pack and unpack 'Runuran' objects")

## --------------------------------------------------------------------------

test_that("[packed-01] compare packed and unpacked object (unbounded domain)", {
    ## create 'Runuran' objects
    gu <- pinv.new(dnorm,lb=-Inf,ub=Inf)
    gp <- pinv.new(dnorm,lb=-Inf,ub=Inf)
    ## pack object 'gp'
    unuran.packed(gp) <- TRUE

    ## create samples
    u <- (-1:(samplesize+1))/samplesize
    msg <- "\\[UNU\\.RAN - warning\\] argument out of domain: U not in \\[0,1\\]"
    expect_output(xu <- uq(gu,u), msg)
    expect_warning(xp <- uq(gp,u), msg)

    ## compare
    expect(isTRUE(all.equal(xu,xp)), "packed and unpacked version of PINV differ !")
})

## --------------------------------------------------------------------------

test_that("[packed-02] compare packed and unpacked object (bounded domain)", {
    ## create 'Runuran' objects
    gu <- pinv.new(dnorm,lb=0,ub=1)
    gp <- pinv.new(dnorm,lb=0,ub=1)
    ## pack object 'gp'
    unuran.packed(gp) <- TRUE

    ## create samples
    u <- (-1:(samplesize+1))/samplesize
    msg <- "\\[UNU\\.RAN - warning\\] argument out of domain: U not in \\[0,1\\]"
    expect_output(xu <- uq(gu,u), msg)
    expect_warning(xp <- uq(gp,u), msg)

    ## compare
    expect(isTRUE(all.equal(xu,xp)), "packed and unpacked version of PINV differ !")
})

## --- End ------------------------------------------------------------------
