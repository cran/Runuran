#! /bin/sh

# ---------------------------------------------------------------------------
# Run this script in the top-level Runuran directory.
# The script removes all generated files and all UNU.RAN files.
# ---------------------------------------------------------------------------

# Check directory
PACKAGE=`grep -s "Package: Runuran" DESCRIPTION`
test -f DESCRIPTION -a "$PACKAGE" || { \
    echo "You must run this script in the parent directory of the top-level Runuran directory"
    exit 1
}

# Remove autotools files
rm -rf configure config.log config.status autom4te.cache
(cd ./src && rm -rf Makevars config.h* )

# Remove compiled files
(cd ./src && rm -rf Runuran.so *.o )
(cd ./inst/doc && \
    rm -f *.aux *.bbl *.blg *.log *.out *.toc && \
    rm -f Runuran.R Runuran.pdf Runuran.tex )

# Remove exported header files files
(cd ./inst/include && rm -f Runuran_ext.h unuran.h )

# Remove UNU.RAN files
(cd ./src/unuran-src && rm -rf * )

# End
exit 0
