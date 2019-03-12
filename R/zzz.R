
.onAttach <- function(library, pkg)
{
    ## we can't do this in .onLoad
    unlockBinding(".Runuran.Options", asNamespace("Runuran"))
    unlockBinding(".unuran.error.level", asNamespace("Runuran"))

    ## print startup message
    desc <- utils::packageDescription("Runuran", fields = c("Package", "Version"))
    msg <- paste0("Package ", sQuote(desc$Package), ", ",
                  "version ",  desc$Version, ": ",
                  "type 'help(",desc$Package,")' for summary information")
    if (interactive()) packageStartupMessage(msg)
    invisible()
}

.onUnload <- function(libpath) {
    ## call garbage collector to run the finalizer
    ## for all Runuran objects.
    gc(verbose=FALSE)
    ## unload DLL
    library.dynam.unload("Runuran", libpath)
}
