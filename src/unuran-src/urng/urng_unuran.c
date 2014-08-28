/* Copyright (c) 2000-2014 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include "urng.h"
#ifdef UNUR_URNG_UNURAN
UNUR_URNG *
unur_urng_new( double (*sampleunif)(void *state), void *state )
{
  UNUR_URNG *urng = NULL;         
  _unur_check_NULL( "URNG", sampleunif, NULL );
  urng = _unur_xmalloc( sizeof(struct unur_urng) );
  urng->sampleunif = sampleunif;
  urng->state      = state;
  urng->samplearray = NULL;
  urng->sync     = NULL;
  urng->seed     = ULONG_MAX;
  urng->setseed  = NULL;
  urng->delete   = NULL;
  urng->reset    = NULL;
  urng->nextsub  = NULL;
  urng->resetsub = NULL;
  urng->anti     = NULL;
  COOKIE_SET(urng,CK_URNG);
  return urng;
} 
int
unur_urng_set_sample_array( UNUR_URNG *urng, 
			   unsigned int (*samplearray)(void *state, double *X, int dim) )
{
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);
  urng->samplearray = samplearray;
  return UNUR_SUCCESS;
} 
int
unur_urng_set_sync( UNUR_URNG *urng, void (*sync)(void *state) )
{
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);
  urng->sync = sync;
  return UNUR_SUCCESS;
} 
int
unur_urng_set_seed( UNUR_URNG *urng, void (*setseed)(void *state, unsigned long seed) )
{
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);
  urng->setseed = setseed;
  return UNUR_SUCCESS;
} 
int
unur_urng_set_anti( UNUR_URNG *urng, void (*setanti)(void *state, int anti) )
{
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);
  urng->anti = setanti;
  return UNUR_SUCCESS;
} 
int
unur_urng_set_reset( UNUR_URNG *urng, void (*reset)(void *state) )
{
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);
  urng->reset = reset;
  return UNUR_SUCCESS;
} 
int
unur_urng_set_nextsub( UNUR_URNG *urng, void (*nextsub)(void *state) )
{
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);
  urng->nextsub = nextsub;
  return UNUR_SUCCESS;
} 
int
unur_urng_set_resetsub( UNUR_URNG *urng, void (*resetsub)(void *state) )
{
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);
  urng->resetsub = resetsub;
  return UNUR_SUCCESS;
} 
int
unur_urng_set_delete( UNUR_URNG *urng, void (*delete)(void *state) )
{
  _unur_check_NULL( "URNG", urng, UNUR_ERR_NULL );
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);
  urng->delete  = delete;
  return UNUR_SUCCESS;
} 
int
unur_urng_sync (UNUR_URNG *urng)
{
  if (urng == NULL) 
    urng = unur_get_default_urng();
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);
  if (urng->sync == NULL) {
    _unur_error("URNG",UNUR_ERR_URNG_MISS,"sync");
    return UNUR_ERR_URNG_MISS;
  }
  urng->sync (urng->state);
  return UNUR_SUCCESS;
}  
int
unur_urng_seed (UNUR_URNG *urng, unsigned long seed)
{
  if (urng == NULL) 
    urng = unur_get_default_urng();
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);
  if (urng->setseed == NULL) {
    _unur_error("URNG",UNUR_ERR_URNG_MISS,"seeding function");
    return UNUR_ERR_URNG_MISS;
  }
  urng->setseed (urng->state,seed);
  urng->seed = seed;
  return UNUR_SUCCESS;
} 
int
unur_urng_anti (UNUR_URNG *urng, int anti)
{
  if (urng == NULL) 
    urng = unur_get_default_urng();
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);
  if (urng->anti == NULL) {
    _unur_error("URNG",UNUR_ERR_URNG_MISS,"antithetic flag");
    return UNUR_ERR_URNG_MISS;
  }
  urng->anti (urng->state,anti);
  return UNUR_SUCCESS;
}  
int
unur_urng_nextsub (UNUR_URNG *urng)
{
  if (urng == NULL) 
    urng = unur_get_default_urng();
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);
  if (urng->nextsub == NULL) {
    _unur_error("URNG",UNUR_ERR_URNG_MISS,"next substream");
    return UNUR_ERR_URNG_MISS;
  }
  urng->nextsub (urng->state);
  return UNUR_SUCCESS;
}  
int
unur_urng_resetsub (UNUR_URNG *urng)
{
  if (urng == NULL) 
    urng = unur_get_default_urng();
  COOKIE_CHECK(urng,CK_URNG,UNUR_ERR_COOKIE);
  if (urng->resetsub == NULL) {
    _unur_error("URNG",UNUR_ERR_URNG_MISS,"reset substream");
    return UNUR_ERR_URNG_MISS;
  }
  urng->resetsub (urng->state);
  return UNUR_SUCCESS;
}  
void
unur_urng_free (UNUR_URNG *urng)
{
  if (urng == NULL) return;  
  COOKIE_CHECK(urng,CK_URNG,RETURN_VOID);
  if (urng->delete != NULL) urng->delete (urng->state);
  free (urng);
  urng = NULL;
  return;
}  
int
unur_gen_sync (UNUR_GEN *gen)
{
  _unur_check_NULL( "URNG", gen, UNUR_ERR_NULL );
  return unur_urng_sync(gen->urng);
}  
int
unur_gen_seed (UNUR_GEN *gen, unsigned long seed)
{
  _unur_check_NULL( "URNG", gen, UNUR_ERR_NULL );
  return unur_urng_seed(gen->urng, seed);
}  
int
unur_gen_anti (UNUR_GEN *gen, int anti)
{
  _unur_check_NULL( "URNG", gen, UNUR_ERR_NULL );
  return unur_urng_anti(gen->urng, anti);
}  
int
unur_gen_reset (UNUR_GEN *gen)
{
  _unur_check_NULL( "URNG", gen, UNUR_ERR_NULL );
  return unur_urng_reset(gen->urng);
}  
int
unur_gen_nextsub (UNUR_GEN *gen)
{
  _unur_check_NULL( "URNG", gen, UNUR_ERR_NULL );
  return unur_urng_nextsub(gen->urng);
}  
int
unur_gen_resetsub (UNUR_GEN *gen)
{
  _unur_check_NULL( "URNG", gen, UNUR_ERR_NULL );
  return unur_urng_resetsub(gen->urng);
}  
#endif    
