#############################################################################
##                                                                         ##
##   Tests for special routines                                            ##
##                                                                         ##
#############################################################################

## --- Load library ---------------------------------------------------------

library(Runuran)

## --- HINV -----------------------------------------------------------------

unr <- unuran.new("normal()","hinv; u_resolution=1.e-13")

# single double as argument
Tmax <- 0
for (U in (0:20)/20) {
        T <- pnorm( uqhinv(unr,U) ) - U
        print (T)
        Tmax <- max(abs(T),Tmax)
}
cat("Max. error =",Tmax,"\n")
if (Tmax > 1.e-12) stop ("Max. error exceeds limit")

# special arguments
for (U in c(-1,-0.001,0,0.5,1.,1.001,NA,NaN) ) {
        cat ("U =",U,"\tX =",uqhinv(unr,U),"\n")
}

# vector argument
U <- (0:20)/20
T <- pnorm( uqhinv(unr,U) ) - U
print (T)
Tmax <- max(abs(T))
cat("Max. error =",Tmax,"\n")
if (Tmax > 1.e-12) stop ("Max. error exceeds limit")

# special arguments
U <- c(-1,-0.001,0,0.5,1.,1.001,NA,NaN)
T <- uqhinv(unr,U)
rbind(U,T)

uqhinv(unr,numeric())

## -- End -------------------------------------------------------------------
