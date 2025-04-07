## Load package 'testthat' 
if(! require("testthat", quietly=TRUE)) {
  message("\nCannot run unit tests -- package 'testthat' is not available!\n")
}

##Sys.setenv(NOT_CRAN="true")

## run testthat tests
test_check("Runuran", reporter=c("summary","check","fail"))

