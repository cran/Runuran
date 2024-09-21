/* Copyright (c) 2000-2024 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <stdarg.h>
#include <ctype.h>
#include <unur_source.h>
#include <methods/unur_methods_source.h>
#include <methods/x_gen.h>
#include "parser_source.h"
#include "parser.h"
#include <urng/urng.h>
#include <distr/distr_source.h>
#include <distr/distr.h>
#include <distr/cemp.h>
#include <distr/cont.h>
#include <distr/corder.h>
#include <distr/cvemp.h>
#include <distr/discr.h>
#include <distributions/unur_distributions.h>
#include <methods/arou.h>
#include <methods/ars.h>
#include <methods/auto.h>
#include <methods/cstd.h>
#include <methods/dari.h>
#include <methods/dau.h>
#include <methods/dgt.h>
#include <methods/dsrou.h>
#include <methods/dss.h>
#include <methods/dstd.h>
#include <methods/empk.h>
#include <methods/empl.h>
#include <methods/gibbs.h>
#include <methods/hinv.h>
#include <methods/hist.h>
#include <methods/hitro.h>
#include <methods/hrb.h>
#include <methods/hrd.h>
#include <methods/hri.h>
#include <methods/itdr.h>
#include <methods/mcorr.h>
#include <methods/mvstd.h>
#include <methods/mvtdr.h>
#include <methods/ninv.h>
#include <methods/norta.h>
#include <methods/nrou.h>
#include <methods/pinv.h>
#include <methods/srou.h>
#include <methods/ssr.h>
#include <methods/tabl.h>
#include <methods/tdr.h>
#include <methods/unif.h>
#include <methods/utdr.h>
#include <methods/vempk.h>
#include <methods/vnrou.h>
#if defined(UNUR_URNG_UNURAN) && defined(UNURAN_HAS_PRNG)
#include <uniform/urng_prng.h>
#endif
#define GENTYPE "STRING"       
static struct unur_distr *_unur_str_distr( char *str_distr );
static struct unur_distr *_unur_str_distr_new( char *distribution );
static struct unur_distr *_unur_str_distr_make_os( UNUR_DISTR *distr, 
						   const char *key,
						   char *type_args, char **args );
static int _unur_str_distr_set( UNUR_DISTR **ptr_distr, const char *key, char *value );
typedef int distr_set_i( UNUR_DISTR *distr, int i );
typedef int distr_set_ii( UNUR_DISTR *distr, int i1, int i2 );
typedef int distr_set_d( UNUR_DISTR *distr, double d );
typedef int distr_set_dd( UNUR_DISTR *distr, double d1, double d2 );
typedef int distr_set_Di( UNUR_DISTR *distr, const double *array, int size );
typedef int distr_set_C( UNUR_DISTR *distr, const char *string );
static int _unur_str_distr_set_i( UNUR_DISTR *distr, const char *key, char *type_args, char **args,
				  distr_set_i set );
static int _unur_str_distr_set_ii( UNUR_DISTR *distr, const char *key, char *type_args, char **args,
				   distr_set_ii set );
static int _unur_str_distr_set_d( UNUR_DISTR *distr, const char *key, char *type_args, char **args,
				  distr_set_d set );
static int _unur_str_distr_set_dd( UNUR_DISTR *distr, const char *key, char *type_args, char **args,
				   distr_set_dd set );
static int _unur_str_distr_set_Di( UNUR_DISTR *distr, const char *key, char *type_args, char **args,
				   distr_set_Di set );
static int _unur_str_distr_set_C( UNUR_DISTR *distr, const char *key, char *type_args, char **args,
				  distr_set_C set );
static struct unur_par *_unur_str_par( char *str_method, const UNUR_DISTR *distr,
				       struct unur_slist *mlist );
static struct unur_par *_unur_str_par_new( const char *method, const UNUR_DISTR *distr );
static int _unur_str_par_set( UNUR_PAR *par, const char *key, char *value,
			      struct unur_slist *mlist );
typedef int par_set_void( UNUR_PAR *par );
typedef int par_set_i( UNUR_PAR *par, int i );
typedef int par_set_ii( UNUR_PAR *par, int i1, int i2 );
typedef int par_set_u( UNUR_PAR *par, unsigned u );
typedef int par_set_d( UNUR_PAR *par, double d );
typedef int par_set_dd( UNUR_PAR *par, double d1, double d2 );
typedef int par_set_iD( UNUR_PAR *par, int size, const double *array );
typedef int par_set_Di( UNUR_PAR *par, const double *array, int size );
static int _unur_str_par_set_void( UNUR_PAR *par, const char *key, char *type_args, char **args,
				   par_set_void set );
static int _unur_str_par_set_i( UNUR_PAR *par, const char *key, char *type_args, char **args,
				par_set_i set );
static int _unur_str_par_set_ii( UNUR_PAR *par, const char *key, char *type_args, char **args,
				   par_set_ii set );
static int _unur_str_par_set_u( UNUR_PAR *par, const char *key, char *type_args, char **args,
				par_set_u set );
static int _unur_str_par_set_d( UNUR_PAR *par, const char *key, char *type_args, char **args,
				par_set_d set );
static int _unur_str_par_set_dd( UNUR_PAR *par, const char *key, char *type_args, char **args,
				 par_set_dd set );
static int _unur_str_par_set_iD( UNUR_PAR *par, const char *key, char *type_args, char **args,
				 par_set_iD set, struct unur_slist *mlist );
static int _unur_str_par_set_Di( UNUR_PAR *par, const char *key, char *type_args, char **args,
				 par_set_Di set, struct unur_slist *mlist );
static UNUR_URNG *_unur_str2urng( char *str_urng );
static int _unur_str_set_args( char *value, char *type_args, char **args, int max_args );
static int _unur_parse_ilist( char *liststr, int **iarray );
static int _unur_parse_dlist( char *liststr, double **darray );
static int _unur_atoi ( const char *str );
static unsigned _unur_atou ( const char *str );
static double _unur_atod ( const char *str );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_str_debug_string( int level, const char *key, const char *value );
static void _unur_str_debug_distr( int level, const char *name, double *params, int n_params );
static void _unur_str_debug_set( int level, const char *key, const char *type, ... );
#endif
static void _unur_str_error_unknown( const char *file, int line, const char *key, const char *type );
static void _unur_str_error_invalid( const char *file, int line, const char *key, const char *type );
static void _unur_str_error_args( const char *file, int line, const char *key );
#define _unur_error_unknown(key,what) \
   do { \
     _unur_str_error_unknown( __FILE__,__LINE__, (key), (what) ); \
   } while (0)
#define _unur_error_invalid(key,what) \
   do { \
     _unur_str_error_invalid( __FILE__,__LINE__, (key), (what) ); \
   } while (0)
#define _unur_error_args(key) \
   do { \
     _unur_str_error_args( __FILE__,__LINE__, (key) ); \
   } while (0)
#define MAX_SET_ARGS  (10)  
#include "stringparser_lists.ch"
struct unur_gen *
unur_str2gen (const char *string)
{
  UNUR_DISTR *distr = NULL;       
  UNUR_PAR   *par   = NULL;       
  UNUR_GEN   *gen   = NULL;       
  UNUR_URNG  *urng  = NULL;       
  char *str_distr   = NULL;       
  char *str_method  = NULL;       
  char *str_urng    = NULL;       
  char *str = NULL;               
  char *token;
  struct unur_slist *mlist;       
  _unur_check_NULL( GENTYPE,string,NULL );
  mlist = _unur_slist_new();
  str = _unur_parser_prepare_string( string );
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag)
    _unur_str_debug_string(0,"[input]",str);
#endif
  str_distr = strtok(str, "&");
  for ( token  = strtok(NULL, "&"); 
        token != NULL;     
        token  = strtok(NULL, "&") ) {
    if (!strncmp(token,"method=",(size_t)7)) {
      str_method = token;
    }
    else if (!strncmp(token,"urng=",(size_t)5)) {
      str_urng = token;
    }
    else {
      _unur_error_unknown(token,"category");
      _unur_slist_free( mlist );
      if (str) free(str);
      return NULL;
    }
  }
  distr = _unur_str_distr(str_distr);
  if ( distr == NULL ) {
    _unur_slist_free( mlist );
    if (str) free(str);
    return NULL;
  }
  if ( str_method != NULL )
    par = _unur_str_par(str_method, distr, mlist);
  else
    par = unur_auto_new(distr);
  gen = unur_init(par);
  unur_distr_free(distr);
  if ( str_urng != NULL )
    if (gen != NULL) {
      if ((urng = _unur_str2urng(str_urng)) != NULL )
	unur_chg_urng(gen, urng);
    }
  _unur_slist_free(mlist);
  if (str) free(str);
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag)
    _unur_str_debug_string(0,"",NULL);
#endif
  return gen;
} 
struct unur_par *
_unur_str2par (const struct unur_distr *distr, const char *string, struct unur_slist **mlist )
{
  UNUR_PAR *par = NULL;           
  char *str = NULL;               
  _unur_check_NULL( GENTYPE,distr,NULL );
  _unur_check_NULL( GENTYPE,string,NULL );
  *mlist = _unur_slist_new();
  str = _unur_parser_prepare_string( string );
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag)
    _unur_str_debug_string(0,"[input]",str);
#endif
  par = _unur_str_par(str, distr, *mlist);
  if (str) free(str);
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag)
    if ( par != NULL )
      _unur_str_debug_string(0,"",NULL);
#endif
  return par;
} 
struct unur_distr *
unur_str2distr (const char *string)
{
  UNUR_DISTR *distr = NULL;       
  char *str = NULL;               
  _unur_check_NULL( GENTYPE,string,NULL );
  str = _unur_parser_prepare_string( string );
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag)
    _unur_str_debug_string(0,"[input]",str);
#endif
  distr = _unur_str_distr(str);
  if (str) free(str);
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag)
    if ( distr != NULL )
      _unur_str_debug_string(0,"",NULL);
#endif
  return distr;
} 
struct unur_gen *
unur_makegen_ssu( const char *distrstr, const char *methodstr, UNUR_URNG *urng )
{
  UNUR_DISTR *distr = NULL;       
  UNUR_PAR   *par   = NULL;       
  UNUR_GEN   *gen   = NULL;       
  char *str_distr   = NULL;       
  char *str_method  = NULL;       
  struct unur_slist *mlist;       
  _unur_check_NULL( GENTYPE, distrstr, NULL );
  mlist = _unur_slist_new();
  str_distr = _unur_parser_prepare_string( distrstr );
  str_method = (methodstr) 
    ? _unur_parser_prepare_string( methodstr )
    : NULL;
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag) {
    _unur_str_debug_string(0,"[input-distr]",str_distr);
    _unur_str_debug_string(0,"[input-method]",
			   (str_method) ? str_method : "(NULL)" );
  }
#endif
  do {
    distr = _unur_str_distr(str_distr);
    if (distr == NULL) break;
    if ( str_method != NULL && strlen(str_method)>0 )
      par = _unur_str_par(str_method, distr, mlist);
    else
      par = unur_auto_new(distr);
    if (par == NULL) break;
    gen = unur_init(par);
    if (gen == NULL) break;
    if (urng != NULL) 
      unur_chg_urng(gen, urng);
  } while (0);
  unur_distr_free(distr);
  _unur_slist_free(mlist);
  if (str_distr) free(str_distr);
  if (str_method) free(str_method);
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag)
    _unur_str_debug_string(0,"",NULL);
#endif
  return gen;
} 
struct unur_gen *
unur_makegen_dsu( const struct unur_distr *distr, const char *methodstr, UNUR_URNG *urng )
{
  UNUR_PAR   *par   = NULL;       
  UNUR_GEN   *gen   = NULL;       
  char *str_method  = NULL;       
  struct unur_slist *mlist;       
  _unur_check_NULL( GENTYPE,distr,NULL );
  mlist = _unur_slist_new();
  str_method = (methodstr)
    ? _unur_parser_prepare_string( methodstr )
    : NULL;
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag) {
    _unur_str_debug_string(0,"[input-distr]","(distribution object)");
    _unur_str_debug_string(0,"[input-method]",
			   (str_method) ? str_method : "(NULL)" );
  }
#endif
  do {
    if ( str_method != NULL && strlen(str_method)>0 )
      par = _unur_str_par(str_method, distr, mlist);
    else
      par = unur_auto_new(distr);
    if (par == NULL) break;
    gen = unur_init(par);
    if (gen == NULL) break;
    if (urng != NULL)
      unur_chg_urng(gen, urng);
  } while (0);
  _unur_slist_free(mlist);
  if (str_method) free(str_method);
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag)
    _unur_str_debug_string(0,"",NULL);
#endif
  return gen;
} 
struct unur_distr *
_unur_str_distr( char *str_distr )
{
  struct unur_distr *distr = NULL;
  char *token;             
  char *next;              
  char *key, *value;       
  for ( token = next = str_distr;
	next != NULL && *token != '\0';
	token = next ) {
    next = strchr(token,';');
    if (next != NULL) {
      *next = '\0';     
      next++;           
    }
    key = token;
    value = strchr(key, '=');
    if (value != NULL) {
      *value = '\0';    
      value++;          
    }
    if (key == str_distr) {
      if (value == NULL) {
	value = key;
      }
      else {
	if (*key != 'd') {
	  _unur_error(GENTYPE,UNUR_ERR_STR_SYNTAX,"key for distribution does not start with 'd'"); 
	  _unur_distr_free(distr);   
	  return NULL;
	}
      }
      if (distr != NULL) {
	_unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,""); 
	_unur_distr_free(distr);
      }
      distr = _unur_str_distr_new(value);
      if (distr == NULL) {
	return NULL;
      }
    }
    else {
      if (_unur_str_distr_set(&distr, key, value)!=UNUR_SUCCESS ) {
	_unur_distr_free(distr);
	return NULL;
      }
    }
  }
  return distr;
} 
struct unur_distr *
_unur_str_distr_make_os (UNUR_DISTR *distr, const char *key, char *type_args, char **args)
{
  int *iarray = NULL;                 
  struct unur_distr *os = NULL;       
  if ( !strcmp(type_args, "tt") ) {
    iarray = _unur_xmalloc( 2*sizeof(double) );
    iarray[0] = _unur_atoi( args[0] );
    iarray[1] = _unur_atoi( args[1] );
  }
  else if ( !strcmp(type_args, "L") ) {
    if ( _unur_parse_ilist( args[0], &iarray ) < 2 ) {
      free (iarray);
      iarray = NULL;
    }
  }
  if (iarray == NULL ) {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return NULL;
  }
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_set(2,key,"ii",iarray[0],iarray[1] );
#endif
  os =  unur_distr_corder_new( distr, iarray[0], iarray[1] );
  _unur_distr_free(distr);
  free (iarray);
  return os;
} 
int
_unur_str_distr_set_i (UNUR_DISTR *distr, const char *key, char *type_args, char **args, distr_set_i set)
{
  int iarg;
  if ( !strcmp(type_args, "t") ) {
    iarg = _unur_atoi( args[0] );
  }
  else if ( strlen(type_args) == 0 ) {
    iarg = 1;
  }
  else {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return UNUR_ERR_STR_INVALID;
  }
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_set(2,key,"i",iarg);
#endif
  return set(distr,iarg);
} 
int
_unur_str_distr_set_ii (UNUR_DISTR *distr, const char *key, char *type_args, char **args, distr_set_ii set)
{
  int *iarray = NULL;
  int iarg[2];
  int result;
  if ( !strcmp(type_args, "tt") ) {
    iarg[0] = _unur_atoi( args[0] );
    iarg[1] = _unur_atoi( args[1] );
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"ii",iarg[0],iarg[1]);
#endif
    return set( distr,iarg[0],iarg[1] );
  }
  else if ( !strcmp(type_args, "L") ) {
    if ( _unur_parse_ilist( args[0], &iarray ) < 2 ) {
      _unur_error_args(key);
      free (iarray);
#ifdef UNUR_ENABLE_LOGGING
      if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
	_unur_str_debug_set(2,key,"!");
#endif
      return UNUR_ERR_STR_INVALID;
    }
    result = set( distr,iarray[0],iarray[1] );
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"ii",iarray[0],iarray[1] );
#endif
    free (iarray);
    return result;
  }
  else {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return UNUR_ERR_STR_INVALID;
  }
} 
int
_unur_str_distr_set_d (UNUR_DISTR *distr, const char *key, char *type_args, char **args, distr_set_d set)
{
  double darg;
  if ( strcmp(type_args, "t") ) {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return UNUR_ERR_STR_INVALID;
  }
  darg = _unur_atod( args[0] );
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_set(2,key,"d",darg);
#endif
  return set(distr,darg);
} 
int
_unur_str_distr_set_dd (UNUR_DISTR *distr, const char *key, char *type_args, char **args, distr_set_dd set)
{
  double *darray = NULL;
  double darg[2];
  int result;
  if ( !strcmp(type_args, "tt") ) {
    darg[0] = _unur_atod( args[0] );
    darg[1] = _unur_atod( args[1] );
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"dd",darg[0],darg[1]);
#endif
    return set( distr,darg[0],darg[1] );
  }
  else if ( !strcmp(type_args, "L") ) {
    if ( _unur_parse_dlist( args[0], &darray ) < 2 ) {
      _unur_error_args(key);
      free (darray);
#ifdef UNUR_ENABLE_LOGGING
      if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
	_unur_str_debug_set(2,key,"!");
#endif
      return UNUR_ERR_STR_INVALID;
    }
    result = set( distr,darray[0],darray[1] );
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"dd",darray[0],darray[1] );
#endif
    free (darray);
    return result;
  }
  else {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return UNUR_ERR_STR_INVALID;
  }
} 
int
_unur_str_distr_set_Di (UNUR_DISTR *distr, const char *key, char *type_args, char **args, distr_set_Di set)
{
  int result;
  int t_size;
  int size = -1;
  double *darray = NULL;
  if ( !strcmp(type_args, "Lt") ) {
    t_size = _unur_atoi( args[1] );
    size = _unur_parse_dlist( args[0], &darray );
    if (size > t_size)
      size = t_size;
  }
  else if ( !strcmp(type_args, "L") ) {
    size = _unur_parse_dlist( args[0], &darray );
  }
  if ( !(size>0) ) {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    result = UNUR_ERR_STR_INVALID;
  }
  else {
    result = set( distr,darray,size );
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"Di",darray,size,size);
#endif
  }
  if (darray) free (darray);
  return result;
} 
int
_unur_str_distr_set_C (UNUR_DISTR *distr, const char *key, char *type_args, char **args, distr_set_C set)
{
  char *string;
  if ( strcmp(type_args, "s") ) {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return UNUR_ERR_STR_INVALID;
  }
  string = args[0];
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_set(2,key,"C",string);
#endif
  return set(distr,string);
} 
struct unur_par *
_unur_str_par( char *str_method, const UNUR_DISTR *distr, struct unur_slist *mlist )
{
  struct unur_par *par = NULL;     
  char *token;             
  char *next;              
  char *key, *value;       
  for ( token = next = str_method;
	next != NULL && *token != '\0';
	token = next ) {
    next = strchr(token,';');
    if (next != NULL) {
      *next = '\0';     
      next++;           
    }
    key = token;
    value = strchr(key, '=');
    if (value != NULL) {
      *value = '\0';    
      value++;          
    }
    if (key == str_method) {
      if (value == NULL) {
	value = key;
      }
      else {
	if (*key != 'm') {
	  _unur_error(GENTYPE,UNUR_ERR_STR_SYNTAX,"key for method does not start with 'm'"); 
	  return NULL;
	}
      }
      par = _unur_str_par_new(value,distr);
      if (par == NULL) {
	_unur_error(GENTYPE,UNUR_ERR_STR,"setting method failed"); 
	return NULL;
      }
    }
    else {
      if ( !_unur_str_par_set(par, key, value, mlist) ) {
	; 
      }
    }
  }
  return par;
} 
int
_unur_str_par_set_void (UNUR_PAR *par, const char *key, 
			char *type_args, char **args ATTRIBUTE__UNUSED, par_set_void set)
{
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_set(2,key,"v");
#endif
  if (*type_args != '\0')
    _unur_error_args(key);
  return set(par);
} 
int
_unur_str_par_set_i (UNUR_PAR *par, const char *key, char *type_args, char **args, par_set_i set)
{
  int iarg;
  if ( !strcmp(type_args, "t") ) {
    iarg = _unur_atoi( args[0] );
  }
  else if ( strlen(type_args) == 0 ) {
    iarg = 1;
  }
  else {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return UNUR_ERR_STR_INVALID;
  }
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_set(2,key,"i",iarg);
#endif
  return set(par,iarg);
} 
int
_unur_str_par_set_ii (UNUR_PAR *par, const char *key, char *type_args, char **args, par_set_ii set)
{
  int *iarray = NULL;
  int iarg[2];
  int result;
  if ( !strcmp(type_args, "tt") ) {
    iarg[0] = _unur_atoi( args[0] );
    iarg[1] = _unur_atoi( args[1] );
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"ii",iarg[0],iarg[1]);
#endif
    return set( par,iarg[0],iarg[1] );
  }
  else if ( !strcmp(type_args, "L") ) {
    if ( _unur_parse_ilist( args[0], &iarray ) < 2 ) {
      _unur_error_args(key);
      free (iarray);
#ifdef UNUR_ENABLE_LOGGING
      if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
	_unur_str_debug_set(2,key,"!");
#endif
      return UNUR_ERR_STR_INVALID;
    }
    result = set( par,iarray[0],iarray[1] );
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"ii",iarray[0],iarray[1] );
#endif
    free (iarray);
    return result;
  }
  else {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return UNUR_ERR_STR_INVALID;
  }
} 
int
_unur_str_par_set_u (UNUR_PAR *par, const char *key, char *type_args, char **args, par_set_u set)
{
  unsigned uarg;
  if ( strcmp(type_args, "t") ) {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return UNUR_ERR_STR_INVALID;
  }
  uarg = _unur_atou( args[0] );
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_set(2,key,"u",uarg);
#endif
  return set(par,uarg);
} 
int
_unur_str_par_set_d (UNUR_PAR *par, const char *key, char *type_args, char **args, par_set_d set)
{
  double darg;
  if ( strcmp(type_args, "t") ) {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return UNUR_ERR_STR_INVALID;
  }
  darg = _unur_atod( args[0] );
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_set(2,key,"d",darg);
#endif
  return set(par,darg);
} 
int
_unur_str_par_set_dd (UNUR_PAR *par, const char *key, char *type_args, char **args, par_set_dd set)
{
  double *darray = NULL;
  double darg[2];
  int result;
  if ( !strcmp(type_args, "tt") ) {
    darg[0] = _unur_atod( args[0] );
    darg[1] = _unur_atod( args[1] );
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"dd",darg[0],darg[1]);
#endif
    return set( par,darg[0],darg[1] );
  }
  else if ( !strcmp(type_args, "L") ) {
    if ( _unur_parse_dlist( args[0], &darray ) < 2 ) {
      _unur_error_args(key);
      free (darray);
#ifdef UNUR_ENABLE_LOGGING
      if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
	_unur_str_debug_set(2,key,"!");
#endif
      return UNUR_ERR_STR_INVALID;
    }
    result = set( par,darray[0],darray[1] );
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"dd",darray[0],darray[1] );
#endif
    free (darray);
    return result;
  }
  else {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    return UNUR_ERR_STR_INVALID;
  }
} 
int
_unur_str_par_set_iD (UNUR_PAR *par, const char *key, char *type_args, char **args,
		      par_set_iD set, struct unur_slist *mlist )
{
  int result;
  int t_size;
  int size = -1;
  double *darray = NULL;
  if ( !strcmp(type_args, "tL") ) {
    t_size = _unur_atoi( args[0] );
    size = _unur_parse_dlist( args[1], &darray );
    if (size > 0) {
      if (size > t_size)  size = t_size;
    }
    else {
      size = t_size;
      if (darray) free (darray);
      darray = NULL;
    }
  }
  else if ( !strcmp(type_args, "t") ) {
    size = _unur_atoi( args[0] );
    darray = NULL;
  }
  else if ( !strcmp(type_args, "L") ) {
    size = _unur_parse_dlist( args[0], &darray );
  }
  if ( !(size>0) ) {
    _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"!");
#endif
    result = UNUR_ERR_STR_INVALID;
  }
  else {
    result = set( par,size,darray );
#ifdef UNUR_ENABLE_LOGGING
    if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
      _unur_str_debug_set(2,key,"iD",size,darray,size);
#endif
  }
  if (darray) 
    _unur_slist_append( mlist, darray );
  return result;
} 
int
_unur_str_par_set_Di (UNUR_PAR *par, const char *key, char *type_args, char **args,
		      par_set_Di set, struct unur_slist *mlist )
{
  int result;
  int t_size;
  int size;
  double *darray = NULL;
  if ( !strcmp(type_args, "Lt") ) {
    t_size = _unur_atoi( args[1] );
    size = _unur_parse_dlist( args[0], &darray );
    if (size > 0) {
      result = set( par,darray,t_size );
#ifdef UNUR_ENABLE_LOGGING
      if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
	_unur_str_debug_set(2,key,"Di",darray,size,size);
#endif
      if (darray) 
	_unur_slist_append( mlist, darray );
      return result;
    }
  }
  _unur_error_args(key);
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_set(2,key,"!");
#endif
  return UNUR_ERR_STR_INVALID;
} 
UNUR_URNG *
_unur_str2urng( char *str_urng ATTRIBUTE__UNUSED)
{
#if defined(UNUR_URNG_UNURAN) && defined(UNURAN_HAS_PRNG)
  UNUR_URNG *urng = NULL;
  char *token;             
  char *next;              
  char *key, *value;       
  for ( token = next = str_urng;
	next != NULL && *token != '\0';
	token = next ) {
    next = strchr(token,';');
    if (next != NULL) {
      *next = '\0';     
      next++;           
    }
    key = token;
    value = strchr(key, '=');
    if (value != NULL) {
      *value = '\0';    
      value++;          
    }
    if ( !strcmp( key, "urng") ) {
      if (key == str_urng) {
	urng = unur_urng_prng_new(value);
	if (urng == NULL) {
	  _unur_error(GENTYPE,UNUR_ERR_STR,"setting URNG failed"); 
	  return NULL;
	}
	else {
#ifdef UNUR_ENABLE_LOGGING
	  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
	    _unur_str_debug_string(1,"urng",value);
#endif
	  ;
	}
      }
      else {
	if (urng) unur_urng_free (urng);
	_unur_error(GENTYPE,UNUR_ERR_STR_SYNTAX,"urng must be first key"); 
	return NULL;
      }
    }
    else {
      if (urng) unur_urng_free (urng);
      _unur_error_unknown(key,"parameter for uniform RNG");
      return NULL;
    }
  }
  return urng;
#else
  _unur_error(GENTYPE,UNUR_ERR_STR,"setting URNG requires PRNG library enabled"); 
  return NULL;
#endif   
} 
int
_unur_str_set_args( char *value, char *type_args, char **args, int max_args )
{
  char *token;             
  char *next;              
  int n_args;              
  n_args = 0;
  *type_args = '\0';
  *args = NULL;
  if (value && *value != '\0') {
    for ( token = next = value;
	  next != NULL && *token != '\0' && n_args < max_args;
	  token = next, n_args++ ) {
      if ( *token == '(' ) {
	type_args[n_args] = 'L';
	token++;
	args[n_args] = token;
	next = strchr(token,')');
	if (next != NULL) {
	  *next = '\0';      
	  token = ++next;    
	  if ( !(*token == ',' || *token == '\0') ) {
	    _unur_error(GENTYPE,UNUR_ERR_STR_SYNTAX,") not followed by comma"); 
	    return -1;
	  }
	}
	else {
	  token = NULL;
	}
      }
      else if ( *token == '"' ) {
	type_args[n_args] = 's';
	token++;
	args[n_args] = token;
	next = strchr(token,'"');
	if (next != NULL) {
	  *next = '\0';      
	  token = ++next;    
	  if ( !(*token == ',' || *token == '\0') ) {
	    _unur_error(GENTYPE,UNUR_ERR_STR_SYNTAX,"closing '\"' not followed by comma"); 
	    return -1;
	  }
	}
	else {
	  token = NULL;
	}
      }
      else {
	type_args[n_args] = 't';
	args[n_args] = token;
      }
      if (token) {
	next = strchr(token,',');
	if (next != NULL) {
	  *next = '\0';     
	  next++;           
	}
      }
      else {
	next = NULL;
      }
    }
    type_args[n_args] = '\0';
    if (n_args >= max_args) { 
      _unur_error( GENTYPE, UNUR_ERR_STR_SYNTAX, "too many arguments");
      return -1;
    }
  }
  return n_args;
} 
int 
_unur_parse_ilist( char *liststr, int **iarray )
{
  int *iarr = NULL;     
  int n_iarray = 0;     
  int n_alloc = 0;      
  char *token;
  if (liststr == NULL) {
    *iarray = iarr;
    return 0;
  }
  while ( *liststr == ',' || *liststr == '(' ){
    liststr++; 
  }
  for ( token  = strtok(liststr, ",)"); 
        token != NULL;     
        token  = strtok(NULL, ",)") ) {
    if (n_iarray >= n_alloc) {
      n_alloc += 100;
      iarr = _unur_xrealloc( iarr, n_alloc * sizeof(int) );
    }
    iarr[n_iarray++] = _unur_atoi(token);
  } 
  *iarray = iarr;
  return n_iarray;
} 
int 
_unur_parse_dlist( char *liststr, double **darray )
{
  double *darr = NULL;  
  int n_darray = 0;     
  int n_alloc = 0;      
  char *token;          
  char *next;           
  if (liststr == NULL) {
    *darray = NULL;
    return 0;
  }
  token = liststr;
  while (*token != '\0' && *token == '(')
    ++token;
  for ( next = token;
	next != NULL && *token != '\0' &&*token != ')';
	token = next ) {
    next = strchr(token,',');
    if (next != NULL) {
      *next = '\0';     
      next++;           
    }
    if (n_darray >= n_alloc) {
      n_alloc += 100;
      darr = _unur_xrealloc( darr, n_alloc * sizeof(double) );
    }
    darr[n_darray] = _unur_atod(token);
    n_darray++;
  }
  *darray = darr;
  return n_darray;
} 
int 
_unur_atoi ( const char *str )
{
  if ( !strcmp(str,"true") || !strcmp(str,"on") )
    return 1;
  else if ( !strcmp(str,"false") || !strcmp(str,"off") )
    return 0;
  else if ( !strncmp(str,"inf",(size_t)3) )
    return INT_MAX;
  else if ( !strncmp(str,"-inf",(size_t)4) )
    return INT_MIN;
  else
    return atoi(str);
} 
unsigned
_unur_atou ( const char *str )
{
  char *tail;
  if ( !strcmp(str,"true") || !strcmp(str,"on") )
    return 1u;
  else if ( !strcmp(str,"false") || !strcmp(str,"off") )
    return 0u;
  else
    return ((unsigned) strtoul(str, &tail, 16));
} 
double
_unur_atod ( const char *str )
{
  if ( !strncmp(str,"inf",(size_t)3) )
    return UNUR_INFINITY;
  else if ( !strncmp(str,"-inf",(size_t)4) )
    return -UNUR_INFINITY;
  else
    return atof(str);
} 
#ifdef UNUR_ENABLE_LOGGING
void 
_unur_str_debug_string( int level, const char *key, const char *value  )
{
  FILE *LOG;
  LOG = unur_get_stream();
  fprintf(LOG,"%s: ",GENTYPE);
  for (; level>0; level--) 
    fprintf(LOG,"\t");
  fprintf(LOG,"%s",key);
  if (value)
    fprintf(LOG,": %s",value);
  fprintf(LOG,"\n");
} 
void 
_unur_str_debug_distr( int level, const char *name, double *params, int n_params )
{
  FILE *LOG;
  int i;
  LOG = unur_get_stream();
  fprintf(LOG,"%s: ",GENTYPE);
  for (; level>0; level--) 
    fprintf(LOG,"\t");
  fprintf(LOG,"distribution = %s (",name);
  if (n_params >0) {
    fprintf(LOG,"%g",params[0]);
    for (i=1; i<n_params; i++)
      fprintf(LOG,", %g",params[i]);
  }
  fprintf(LOG,")\n");
} 
void
_unur_str_debug_set( int level, const char *key, const char *type, ... )
{
  va_list ap;
  FILE *LOG;
  va_start(ap, type);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: ",GENTYPE);
  for (; level>0; level--) 
    fprintf(LOG,"\t");
  fprintf(LOG,"%s: ",key);
  while (1) {
    switch (*type) {
    case 'v':
      fprintf(LOG," -none-");
      break;
    case 'd':
      fprintf(LOG," %g",va_arg(ap,double));
      break;
    case 'i':
      fprintf(LOG," %d",va_arg(ap,int));
      break;
    case 'u':
      fprintf(LOG," %x",va_arg(ap,unsigned int));
      break;
    case 'C':
      fprintf(LOG," %s",va_arg(ap,char *));
      break;
    case 'D': {
      int i,size;
      double *darray;
      darray = va_arg(ap,double *);
      size = va_arg(ap,int);
      if (size > 0 && darray != NULL) {
	fprintf(LOG," (%g",darray[0]);
	for (i=1; i<size; i++)
	  fprintf(LOG,",%g",darray[i]);
	fprintf(LOG,")");
      }
      else
	fprintf(LOG," (empty)");
      break;
    }
    case '!':
    default:
      fprintf(LOG," syntax error");
      break;
    }
    if ( *(++type) != '\0' )
      fprintf(LOG,",");
    else
      break;
  }
  fprintf(LOG,"\n");
  fflush(LOG);   
  va_end(ap);
} 
#endif   
void
_unur_str_error_unknown( const char *file, int line, const char *key, const char *type )
{
  struct unur_string *reason = _unur_string_new();
  _unur_string_append( reason, "unknown %s: '%s'", type, key );
  _unur_error_x( GENTYPE, file, line, "error", UNUR_ERR_STR_UNKNOWN,reason->text);
  _unur_string_free( reason );
} 
void
_unur_str_error_invalid( const char *file, int line, const char *key, const char *type )
{
  struct unur_string *reason = _unur_string_new();
  _unur_string_append( reason, "invalid data for %s '%s'", type, key );
  _unur_error_x( GENTYPE, file, line, "error", UNUR_ERR_STR_INVALID,reason->text);
  _unur_string_free( reason );
} 
void
_unur_str_error_args( const char *file, int line, const char *key )
{
  struct unur_string *reason = _unur_string_new();
  _unur_string_append( reason, "invalid argument string for '%s'", key );
  _unur_error_x( GENTYPE, file, line, "error", UNUR_ERR_STR_INVALID,reason->text);
  _unur_string_free( reason );
} 
