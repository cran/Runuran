/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#define UNUR_STDGEN_DEFAULT   0        
#define UNUR_STDGEN_INVERSION (~0u)    
#define UNUR_STDGEN_FAST      (0)      
UNUR_PAR *unur_cstd_new( const UNUR_DISTR *distribution );
int unur_cstd_set_variant( UNUR_PAR *parameters, unsigned variant );
int unur_cstd_chg_truncated( UNUR_GEN *generator, double left, double right );
double unur_cstd_eval_invcdf( const UNUR_GEN *generator, double u );
