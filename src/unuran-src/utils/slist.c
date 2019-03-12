/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <unur_source.h>
struct unur_slist *
_unur_slist_new( void )
{
  struct unur_slist *slist;
  slist = _unur_xmalloc(sizeof(struct unur_slist));
  COOKIE_SET(slist,CK_SLIST);
  slist->ptr   = NULL;
  slist->n_ptr = 0;
  return slist;
} 
int
_unur_slist_length( const struct unur_slist *slist )
{
  CHECK_NULL(slist,0);
  COOKIE_CHECK(slist,CK_SLIST,0);
  if (slist->ptr==NULL)
    return 0;
  return (slist->n_ptr);
} 
void *
_unur_slist_get( const struct unur_slist *slist, int n )
{
  CHECK_NULL(slist,NULL);
  COOKIE_CHECK(slist,CK_SLIST,NULL);
  if (slist->ptr==NULL || n >= slist->n_ptr || n < 0) {
    _unur_warning("list",UNUR_ERR_GENERIC,"element does not exist");
    return NULL;
  }
  return (slist->ptr[n]);
} 
int
_unur_slist_append( struct unur_slist *slist, void *element )
{
  CHECK_NULL(slist,UNUR_ERR_NULL);
  COOKIE_CHECK(slist,CK_SLIST,UNUR_ERR_COOKIE);
  slist->ptr = _unur_xrealloc(slist->ptr,(slist->n_ptr+1)*sizeof(void *));
  slist->ptr[slist->n_ptr] = element;
  ++(slist->n_ptr);
  return UNUR_SUCCESS;
} 
void *
_unur_slist_replace( struct unur_slist *slist, int n, void *element )
{
  void *old_element; 
  CHECK_NULL(slist,NULL);
  COOKIE_CHECK(slist,CK_SLIST,NULL);
  if (slist->ptr==NULL || n >= slist->n_ptr || n < 0) {
    _unur_warning("list",UNUR_ERR_GENERIC,"element does not exist");
    return NULL;
  }
  old_element = slist->ptr[n];
  slist->ptr[n] = element;
  return old_element;
} 
void 
_unur_slist_free( struct unur_slist *slist )
{
  int i;
  if (slist == NULL) return;  
  COOKIE_CHECK(slist,CK_SLIST,RETURN_VOID);
  if ( slist->ptr != NULL ) {
    for (i=0; i < slist->n_ptr; i++)
      if (slist->ptr[i]) free(slist->ptr[i]); 
    free(slist->ptr);
    slist->ptr = NULL;
  }
  free (slist);
} 
