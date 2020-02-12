#############################################################################
##                                                                         ##
##   Tests for unuran.details                                              ##
##   (all skipped on CRAN)                                                 ##
##                                                                         ##
#############################################################################

## --- Test Parameters ------------------------------------------------------

## whether Rout files are updated
##update.Rout <- TRUE
update.Rout <- FALSE

SEED <- 123456

mkmsg.e <- function(...) { makemsg.e("unuran\\.details",...) }

## --------------------------------------------------------------------------

context("[details] - print 'Runuran' objects")

## --- Auxiliary functions --------------------------------------------------

test_unuran.details <- function(distr, method, name=toupper(method), skip.on.cran=TRUE) {
    ## run unuran.details() with various arguments
    ##
    ## if skip.on.cran == TRUE then the test is only run if
    ## environment variable NOT_CRAN is set to "true"
    ## (e.g., by Sys.setenv(NOT_CRAN="true") in a calling R script)
    ##
    ## Using skip.on.cran=TRUE is recommended if the output of unuran.details()
    ## shows too many digits of accuracy, e.g., for method PINV,
    ## or numbers like INF which may look weird (1.#INF) on Windows.
    ## ----------------------------------------------------------------------
    
    set.seed(SEED)
    unr <- unuran.new(distr, method)
    test_that(paste0("[details-",name,"]"), {
        if (isTRUE(skip.on.cran)) { skip_on_cran() }
        expect_known_output( {
            print(unuran.details(unr,show=TRUE, return.list=FALSE))
            print(unuran.details(unr,show=FALSE,return.list=TRUE))
            print(unuran.details(unr,show=FALSE,debug=TRUE))
        },
        file=file.path("saves", paste0(name,".Rout")),
        update=update.Rout)
    })
}

## --- Discrete distributions -----------------------------------------------

## DARI
test_unuran.details("binomial(20,0.5)", "dari")

## DAU
test_unuran.details("binomial(20,0.5)", "dau")

## DGT
test_unuran.details("binomial(20,0.5)", "dgt")

## DSROU
test_unuran.details("binomial(20,0.5)", "dsrou")

## DSS
test_unuran.details("binomial(20,0.5)", "dss")

## DSTD
test_unuran.details("binomial(20,0.5)", "dstd")

## --- Continuous distributions ---------------------------------------------

## AROU
test_unuran.details("normal(1,2)", "arou")

## ARS
test_unuran.details("normal(1,2)", "ars")

## CSTD
test_unuran.details("normal(1,2)", "cstd")

## HINV
test_unuran.details("normal(1,2)", "hinv", skip.on.cran=TRUE)

## HRB
## test_unuran.details("normal(1,2)", "hrb")

## HRD
## test_unuran.details("normal(1,2)", "hrd")

## HRI
## test_unuran.details("normal(1,2)", "hri")

## ITDR
test_unuran.details("gamma(0.5)", "itdr")

## NINV
test_unuran.details("normal(1,2)", "ninv")

## NROU
test_unuran.details("normal(1,2)", "nrou")

## PINV
test_unuran.details("normal(1,2)", "pinv", skip.on.cran=TRUE)

## SROU
test_unuran.details("normal(1,2)", "srou")

## SSR
test_unuran.details("normal(1,2)", "ssr")

## TABL
test_unuran.details("normal(1,2)", "tabl")

## TDR
test_unuran.details("normal(1,2)", "tdr")

## UTDR
test_unuran.details("normal(1,2)", "utdr")

## MIXT

## --- Continuous Multivariate distributions --------------------------------

## --------------------------------------------------------------------------

context("[details] - Invalid arguments")

## --------------------------------------------------------------------------

test_that("[details-i01] calling unuran.details with invalid arguments", {

    msg <- mkmsg.e("Argument 'unr' must be of class 'unuran'")
    expect_error( unuran.details(1),  msg)

})

## --- End ------------------------------------------------------------------

