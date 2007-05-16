/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
unsigned _unur_default_debugflag = UNUR_DEBUGFLAG_DEFAULT;
int
unur_set_debug( struct unur_par *par, unsigned debug )
{
  _unur_check_NULL( NULL,par,UNUR_ERR_NULL );
#ifdef UNUR_ENABLE_LOGGING
  par->debug = debug;
  return UNUR_SUCCESS;
#else
  _unur_warning("DEBUG",UNUR_ERR_COMPILE,"debugging not enabled");
  return UNUR_ERR_COMPILE;
#endif
} 
int
unur_chg_debug( struct unur_gen *gen, unsigned debug )
{
  CHECK_NULL( gen, UNUR_ERR_NULL );
#ifdef UNUR_ENABLE_LOGGING
  gen->debug = debug;
  return UNUR_SUCCESS;
#else
  _unur_warning("DEBUG",UNUR_ERR_COMPILE,"debugging not enabled");
  return UNUR_ERR_COMPILE;
#endif
} 
int
unur_set_default_debug( unsigned debug )
{
  _unur_default_debugflag = debug;
  return UNUR_SUCCESS;
} 
char * 
_unur_make_genid( const char *gentype )
{
  static int count = 0;   
  char *genid;
  genid = _unur_xmalloc(sizeof(char)*(strlen(gentype) + 6));
  ++count; count %= 1000;      
  sprintf(genid,"%s.%03d",gentype,count);
  return genid;
} 
