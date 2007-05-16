# Makefile for Unices

# --- Package name -----------------------------------------------------------

pkg = Runuran

# --- Location of UNU.RAN ----------------------------------------------------

UNURAN_DIR := unuran-src
UNURAN_SRC := $(wildcard ${UNURAN_DIR}/*/*.c)
UNURAN_CPPFLAGS := -I. -I${UNURAN_DIR} -DHAVE_CONFIG_H

# --- File lists -------------------------------------------------------------

R_UNURAN_SRC := Runuran.c Runuran_distr.c ${UNURAN_SRC}
R_UNURAN_OBJ := $(R_UNURAN_SRC:.c=.o)
R_UNURAN_SHLIB := $(pkg)$(SHLIB_EXT)

# --- Compiler and linker flags ----------------------------------------------

PKG_CPPFLAGS = ${UNURAN_CPPFLAGS}

# --- Create DLL -------------------------------------------------------------

all: $(R_UNURAN_SHLIB)

$(R_UNURAN_SHLIB): $(R_UNURAN_OBJ)
	$(SHLIB_LINK) -o $@ $(R_UNURAN_OBJ) $(ALL_LIBS)

# --- Clear working space ----------------------------------------------------

mostlyclean: clean

clean:
	@echo "Clear working space"
	@find ./ -type f -name '*.o' -o -name '*.so' -exec rm -v {} ';'

# --- End --------------------------------------------------------------------
