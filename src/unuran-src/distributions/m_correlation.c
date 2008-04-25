/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <distr/distr_source.h>
#include <distr/distr.h>
#include <distr/matr.h>
#include <specfunct/unur_specfunct_source.h>
#include "unur_distributions.h"
#include "unur_distributions_source.h"
#include "unur_stddistr.h"
static const char distr_name[] = "correlation";
#define DISTR distr->data.matr
struct unur_distr *
unur_distr_correlation( int n )
{
  struct unur_distr *distr;
  distr = unur_distr_matr_new(n,n);
  if (distr == NULL) {
    return NULL;
  }
  distr->id = UNUR_DISTR_MCORRELATION;
  distr->name = distr_name;
  DISTR.init = NULL;
  return distr;
} 
#undef DISTR
