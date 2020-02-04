/* Copyright (c) 2000-2020 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/distr.h>
#include <distr/cvec.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
static const char distr_name[] = "copula";
#define DISTR distr->data.cvec
struct unur_distr *
unur_distr_copula( int dim, const double *rankcorr )
{
  struct unur_distr *distr;
  struct unur_distr *marginal;
  distr = unur_distr_cvec_new(dim);
  if (distr == NULL) {
    return NULL;
  }
  distr->id = UNUR_DISTR_COPULA;
  distr->name = distr_name;
  DISTR.init = NULL;
  if ( unur_distr_cvec_set_rankcorr(distr,rankcorr)!=UNUR_SUCCESS ) {
    unur_distr_free( distr );
    return NULL;
  }
  marginal = unur_distr_uniform(NULL,0);
  unur_distr_cvec_set_marginals(distr,marginal);
  unur_distr_free(marginal);
  distr->set |= ( 0 );
  return distr;
} 
#undef DISTR
