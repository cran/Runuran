## --------------------------------------------------------------------------
##
## Check functions for auxiliary random number generator:
##   set.aux.seed(), use.aux.urng() 
##
## --------------------------------------------------------------------------

## --- Test Parameters ------------------------------------------------------

## size of sample for test
samplesize <- 1.e5

SEED <- 123456
AUX.SEED <- 78910

## --------------------------------------------------------------------------

context("[aux-gen] - use auxiliary generator")

## --------------------------------------------------------------------------

test_that("[aux-gen-01] calling use.aux.urng(): aux.gen not available", {
    gen <- unuran.new(udnorm(), "srou");
    msg <- "use.aux.urng() must return NA"
    expect_true(is.na(use.aux.urng(gen)), msg)
})
    
## --------------------------------------------------------------------------

test_that("[aux-gen-02] calling use.aux.urng(): aux.gen available", {
    gen1 <- unuran.new(udnorm(), "tdr; cpoints=2; max_sqhratio=0.5; usedars=on")
    gen2 <- unuran.new(udexp(), "tdr; cpoints=2; max_sqhratio=0.5; usedars=on")

    ## there should be no correlation
    set.seed(SEED)
    x1 <- ur(gen1,samplesize)
    set.seed(SEED)
    x2 <- ur(gen2,samplesize)
    expect_lt(abs(cor(x1,x2)), 0.5, label="Correlation between streams without auxgen")

    ## aux.gen available but not used
    expect_false(use.aux.urng(gen1))

    ## use aux.gen
    use.aux.urng(gen1) <- TRUE
    use.aux.urng(gen2) <- TRUE

    ## aux.gen is used now
    expect_true(use.aux.urng(gen1))

    ## there should be high correlation now
    set.seed(SEED); set.aux.seed(AUX.SEED)
    x1 <- ur(gen1,samplesize)
    set.seed(SEED); set.aux.seed(AUX.SEED)
    x2 <- ur(gen2,samplesize)
    expect_gt(abs(cor(x1,x2)), 0.5, label="Correlation between streams *with* auxgen")

    ## switch off aux.gen
    use.aux.urng(gen1) <- FALSE
    use.aux.urng(gen2) <- FALSE

    ## aux.gen is not used any more
    expect_false(use.aux.urng(gen1))
    set.seed(SEED)
    x1 <- ur(gen1,samplesize)
    set.seed(SEED)
    x2 <- ur(gen2,samplesize)
    expect_lt(abs(cor(x1,x2)), 0.5, label="Correlation between streams without auxgen")
})
    
## --------------------------------------------------------------------------

test_that("[aux-gen-01] calling use.aux.urng(): reseed aux.gen", {
    gen <- unuran.new(udnorm(), "tdr; cpoints=2; max_sqhratio=0.5; usedars=on")
    use.aux.urng(gen) <- TRUE

    ## aux.gen not reseeded
    set.seed(SEED)
    x1 <- ur(gen,samplesize)
    set.seed(SEED)
    x2 <- ur(gen,samplesize)
    expect(!isTRUE(all.equal(x1,x2)), "aux.gen *not* reseeded: streams are equal but should differ")

    ## aux.gen is reseeded
    set.seed(SEED); set.aux.seed(AUX.SEED)
    x1 <- ur(gen,samplesize)
    set.seed(SEED); set.aux.seed(AUX.SEED)
    x2 <- ur(gen,samplesize)
    expect(isTRUE(all.equal(x1,x2)), "aux.gen reseeded: streams differ but should be equal")
})
    
## --------------------------------------------------------------------------

context("[aux-gen] - Invalid arguments")

## --------------------------------------------------------------------------

test_that("[aux-gen-i01] calling set.aux.seed() invalid arguments", {

    msg <- "argument \"seed\" is missing, with no default"
    expect_error(set.aux.seed(), msg, label="set.aux.seed() with missing seed")

    msg <- "seed must be positive integer"
    expect_error(set.aux.seed(0), msg, label="set.aux.seed(0) with invvalid seed=0")
    expect_error(set.aux.seed(-1), msg, label="set.aux.seed(0) with invvalid seed<0")
})

## --- End ------------------------------------------------------------------
