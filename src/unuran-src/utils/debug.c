/* Copyright (c) 2000-2022 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
unsigned _unur_default_debugflag = UNUR_DEBUGFLAG_DEFAULT;
int
unur_set_debug( struct unur_par *par ATTRIBUTE__UNUSED,
		unsigned debug ATTRIBUTE__UNUSED )
{
#ifdef UNUR_ENABLE_LOGGING
  _unur_check_NULL( NULL,par,UNUR_ERR_NULL );
  par->debug = debug;
  return UNUR_SUCCESS;
#else
  _unur_warning("DEBUG",UNUR_ERR_COMPILE,"debugging not enabled");
  return UNUR_ERR_COMPILE;
#endif
} 
int
unur_chg_debug( struct unur_gen *gen ATTRIBUTE__UNUSED,
		unsigned debug ATTRIBUTE__UNUSED )
{
#ifdef UNUR_ENABLE_LOGGING
  CHECK_NULL( gen, UNUR_ERR_NULL );
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
  size_t len;
  len = strlen(gentype);
  genid = _unur_xmalloc(sizeof(char)*(len+5));
  ++count; count %= 1000;      
#if HAVE_DECL_SNPRINTF
  snprintf(genid, len+5, "%s.%03d", gentype, count);
#else
  sprintf(genid, "%s.%03d", gentype, count);
#endif
  return genid;
} 
