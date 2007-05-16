/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
void*
_unur_xmalloc(size_t size)
{
  register void *ptr;
  ptr = malloc( size );
  if (ptr == NULL) {
    _unur_error(NULL,UNUR_ERR_MALLOC,"");
    exit (EXIT_FAILURE);
  }
  return ptr;
} 
void*
_unur_xrealloc(void *ptr, size_t size)
{
  register void *new_ptr;
  new_ptr = realloc( ptr, size );
  if (new_ptr == NULL) {
    _unur_error(NULL,UNUR_ERR_MALLOC,"");
    exit (EXIT_FAILURE);
  }
  return new_ptr;
} 
