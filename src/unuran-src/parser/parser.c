/* Copyright (c) 2000-2021 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <ctype.h>
#include <unur_source.h>
#include "parser.h"
#include "parser_source.h"
char *
_unur_parser_prepare_string( const char *str )
{
  char *tmp, *ptr;
  char *new;       
  size_t len;      
  len = strlen(str)+1;
  new = _unur_xmalloc( len * sizeof(char) );
  ptr = memcpy(new,str,len);
  for (tmp = ptr; *tmp != '\0'; tmp++)
    if ( !isspace(*tmp) ) {
      *ptr = tolower(*tmp);
      if (*ptr == '\'') *ptr = '"';
      ptr++;
    }
  *ptr = '\0';
  return new;
} 
