/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef SLIST_H_SEEN
#define SLIST_H_SEEN
struct unur_slist *_unur_slist_new( void );
int _unur_slist_append( struct unur_slist *slist, void *element );
int _unur_slist_length( const struct unur_slist *slist );
void *_unur_slist_get( const struct unur_slist *slist, int n );
void *_unur_slist_replace( struct unur_slist *slist, int n, void *element );
void _unur_slist_free( struct unur_slist *slist );
#endif  
