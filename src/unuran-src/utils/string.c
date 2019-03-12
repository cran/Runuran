/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
#include <time.h>
#include <stdarg.h>
#include <ctype.h>
#define MEMBLOCKSIZE     128    
#define MAXSTRINGSIZE   1024    
struct unur_string *
_unur_string_new ( void )
{
  struct unur_string *string;
  string = _unur_xmalloc(sizeof(struct unur_string));
  string->text = NULL;
  string->length = 0;
  string->allocated = 0;
  return string;
} 
int
_unur_string_append ( struct unur_string *string, const char *format, ... )
{
  size_t len;
  va_list ap;
  va_start(ap, format);
  while (string->length + MAXSTRINGSIZE + 1 > string->allocated) {
    string->allocated += MEMBLOCKSIZE;
    string->text = _unur_xrealloc( string->text, (size_t)string->allocated );
  }
#if HAVE_DECL_VSNPRINTF
  len = vsnprintf (string->text+string->length, (size_t)MAXSTRINGSIZE, format, ap);
#else
  len = vsprintf (string->text+string->length, format, ap);
  if (len >= MAXSTRINGSIZE) {
    _unur_error("UTIL",UNUR_ERR_SHOULD_NOT_HAPPEN,"string too long");
    exit (-1);   
  }
#endif
  string->length += len;
  va_end(ap);
  return UNUR_SUCCESS;
} 
int 
_unur_string_appendtext ( struct unur_string *string, const char *text )
{
  int len;
  len = strlen(text);
  while (string->length + len + 1 > string->allocated) {
    string->allocated += MEMBLOCKSIZE;
    string->text = _unur_xrealloc( string->text, (size_t)string->allocated );
  }
  strncpy( string->text+string->length, text, len+1 );
  string->length += len;
  return UNUR_SUCCESS;
} 
void
_unur_string_free ( struct unur_string *string )
{
  if (string) {
    if (string->text)  free (string->text);
    free (string);
    string = NULL;
  }
} 
void
_unur_string_clear ( struct unur_string *string )
{
  if (string) {
    string->length = 0;
    *(string->text) = '\0';
  }
} 
