/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef FUNCTPARSER_SOURCE_H_SEEN
#define FUNCTPARSER_SOURCE_H_SEEN
struct ftreenode *_unur_fstr2tree ( const char *functstring );
struct ftreenode *_unur_fstr2tree_DefFunct ( const char *functstring );
double _unur_fstr_eval_tree ( const struct ftreenode *functtree_root, double x );
struct ftreenode *_unur_fstr_dup_tree (const struct ftreenode *functtree_root);
void _unur_fstr_free ( struct ftreenode *functtree_root );
char *_unur_fstr_tree2string ( const struct ftreenode *functtree_root,
			       const char *variable, const char *function, int spaces );
int _unur_fstr_tree2C ( FILE *out, const struct ftreenode *root,
			const char *variable, const char *function );
int _unur_fstr_tree2FORTRAN ( FILE *out, const struct ftreenode *root,
			      const char *variable, const char *function );
int _unur_fstr_tree2JAVA ( FILE *out, const struct ftreenode *root,
			   const char *variable, const char *function );
struct ftreenode *_unur_fstr_make_derivative ( const struct ftreenode *functtree_root );
#endif   
