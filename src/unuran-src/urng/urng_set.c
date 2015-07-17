/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include "urng.h"
int
unur_set_urng( struct unur_par *par, UNUR_URNG *urng )
{
  _unur_check_NULL( NULL, par, UNUR_ERR_NULL );
  _unur_check_NULL("URNG", urng, UNUR_ERR_NULL);
  par->urng = urng;
  if (par->urng_aux) par->urng_aux = urng;
  return UNUR_SUCCESS;
} 
UNUR_URNG *
unur_get_urng( struct unur_gen *gen )
{
  CHECK_NULL(gen,NULL);
  return gen->urng;
} 
UNUR_URNG *
unur_chg_urng( struct unur_gen *gen, UNUR_URNG *urng )
{
  UNUR_URNG *urng_old;
  CHECK_NULL(gen,NULL);
  CHECK_NULL(urng,NULL);
  urng_old = gen->urng;
  gen->urng = urng;
  if (gen->gen_aux)
    unur_chg_urng(gen->gen_aux,urng);
  if (gen->gen_aux_list && gen->n_gen_aux_list) {
    int i;
    for (i=0; i<gen->n_gen_aux_list; i++) {
      if (gen->gen_aux_list[i])
	unur_chg_urng(gen->gen_aux_list[i],urng);
    }
  }
  if (gen->urng_aux) gen->urng_aux = urng;
  return urng_old;
} 
int
unur_set_urng_aux( struct unur_par *par, UNUR_URNG *urng_aux )
{
  _unur_check_NULL( NULL, par, UNUR_ERR_NULL );
  _unur_check_NULL("URNGaux", urng_aux, UNUR_ERR_NULL);
  if (par->urng_aux == NULL)
    return UNUR_ERR_GENERIC;
  par->urng_aux = urng_aux;
  return UNUR_SUCCESS;
} 
UNUR_URNG *
unur_get_urng_aux( struct unur_gen *gen )
{
  CHECK_NULL(gen,NULL);
  return gen->urng_aux;
} 
UNUR_URNG *
unur_chg_urng_aux( struct unur_gen *gen, UNUR_URNG *urng_aux )
{
  UNUR_URNG *urng_aux_old;
  CHECK_NULL(gen,NULL);
  CHECK_NULL(urng_aux,NULL);
  if (gen->urng_aux == NULL) 
    return NULL;
  urng_aux_old = gen->urng_aux;
  gen->urng_aux = urng_aux;
  if (gen->gen_aux)
    unur_chg_urng_aux(gen->gen_aux,urng_aux);
  if (gen->gen_aux_list && gen->n_gen_aux_list) {
    int i;
    for (i=0; i<gen->n_gen_aux_list; i++) {
      if (gen->gen_aux_list[i])
	unur_chg_urng_aux(gen->gen_aux_list[i],urng_aux);
    }
  }
  return urng_aux_old;
} 
int
unur_use_urng_aux_default( UNUR_PAR *par )
{
  if (par->urng_aux == NULL)
    return UNUR_ERR_GENERIC;
  par->urng_aux = unur_get_default_urng_aux();
  return UNUR_SUCCESS;
} 
int
unur_chgto_urng_aux_default( UNUR_GEN *gen )
{
  if (gen->urng_aux == NULL)
    return UNUR_ERR_GENERIC;
  gen->urng_aux = unur_get_default_urng_aux();
  return UNUR_SUCCESS;
} 
