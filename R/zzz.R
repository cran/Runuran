
.onLoad <- function(lib, pkg) {
    library.dynam("Runuran", pkg, lib)
}

.onUnload <- function(libpath) {
    ## call garbage collector to run the finalizer
    ## for all Runuran objects.
    gc(verbose=FALSE)
    ## unload DLL
    library.dynam.unload("Runuran", libpath)
}
