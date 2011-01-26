/* Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

UNUR_GEN *unur_str2gen( const char *string );
UNUR_DISTR *unur_str2distr( const char *string );
UNUR_GEN *unur_makegen_ssu( const char *distrstr, const char *methodstr, UNUR_URNG *urng );
UNUR_GEN *unur_makegen_dsu( const UNUR_DISTR *distribution, const char *methodstr, UNUR_URNG *urng );
UNUR_PAR *_unur_str2par( const UNUR_DISTR *distribution, const char *method, struct unur_slist **mlist );
