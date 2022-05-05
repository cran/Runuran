/* Copyright (c) 2000-2022 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef UNUR_DEPRECATED_DISTR_H_SEEN
#define UNUR_DEPRECATED_DISTR_H_SEEN
int unur_distr_cvec_set_stdmarginals( UNUR_DISTR *distribution, UNUR_DISTR *marginal );
int unur_distr_cvec_set_stdmarginal_array( UNUR_DISTR *distribution, UNUR_DISTR **marginals );
int unur_distr_cvec_set_stdmarginal_list( UNUR_DISTR *distribution, ... );
const UNUR_DISTR *unur_distr_cvec_get_stdmarginal( const UNUR_DISTR *distribution, int n );
#endif  
