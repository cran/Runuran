/* Copyright (c) 2000-2015 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

char *_unur_make_genid( const char *gentype );
#define _unur_set_genid(gentype) _unur_make_genid(gentype)
#define _unur_free_genid(gen)    do {if (gen->genid) free((gen)->genid);} while(0)
extern unsigned _unur_default_debugflag;     
#define _unur_print_if_default(par,flag)   if(!((par)->set & (flag))) fprintf(LOG,"  [default]")
#ifdef UNUR_ENABLE_CHECKNULL
#define CHECK_NULL(ptr,rval)             \
  if (!(ptr)) {                          \
    _unur_error(NULL,UNUR_ERR_NULL,"");  \
    return rval;                         \
  }
#else               
#define CHECK_NULL(ptr,rval)  do {} while(0)
#endif
#define _unur_check_NULL(gid,ptr,rval)    \
  if (!(ptr)) {                           \
    _unur_error((gid),UNUR_ERR_NULL,"");  \
    return rval;                          \
  }
#define RETURN_VOID ;
