/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include "urng.h"
double
unur_urng_sample (UNUR_URNG *urng)
{
  if (urng == NULL) 
    urng = unur_get_default_urng();
  return _unur_call_urng(urng);
}  
double
unur_sample_urng (UNUR_GEN *gen)
{
  struct unur_urng *urng = (gen) ? gen->urng : unur_get_default_urng(); 
  return _unur_call_urng(urng);
} 
int
unur_urng_sample_array (UNUR_URNG *urng, double *X, int dim)
{
  if (urng == NULL) 
    urng = unur_get_default_urng();
  if (urng->samplearray) {
    return (urng->samplearray(urng->state,X,dim));
  }
  else {
    int i;
    for (i=0; i<dim; i++) 
      X[i] = _unur_call_urng(urng);
    return dim;
  }
}  
int
unur_urng_reset (UNUR_URNG *urng)
{
  if (urng == NULL) 
    urng = unur_get_default_urng();
#ifdef UNUR_URNG_UNURAN
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);
  if (urng->reset != NULL) {
    urng->reset (urng->state);
    return UNUR_SUCCESS;
  }
  if (urng->setseed != NULL && urng->seed != ULONG_MAX) {
    unur_urng_seed(urng,urng->seed);
    return UNUR_SUCCESS;
  }
  _unur_error("URNG",UNUR_ERR_URNG_MISS,"reset");
  return UNUR_ERR_URNG_MISS;
#else
  return _unur_call_reset(urng);
#endif
}  
