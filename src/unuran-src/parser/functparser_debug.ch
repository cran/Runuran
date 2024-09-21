/* Copyright (c) 2000-2024 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifdef UNUR_ENABLE_LOGGING
void
_unur_fstr_debug_input ( const char *fstr )
{
  FILE *LOG = unur_get_stream();
  fprintf(LOG,"%s: Input string:\n",GENTYPE);
  fprintf(LOG,"%s:   %s\n",GENTYPE,fstr);
  fprintf(LOG,"%s:\n",GENTYPE);
} 
void
_unur_fstr_debug_token ( const struct parser_data *pdata )
{
  FILE *LOG = unur_get_stream();
  int i;
  CHECK_NULL(pdata,RETURN_VOID);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,RETURN_VOID);
  fprintf(LOG,"%s: Tokenized string (token separated by blanks):\n",GENTYPE);
  fprintf(LOG,"%s:   ",GENTYPE);
  for (i=0; i<pdata->n_tokens; i++)
    fprintf(LOG,"%s ",pdata->tpos[i]);
  fprintf(LOG,"\n%s:\n",GENTYPE);
} 
void
_unur_fstr_debug_tree( const struct parser_data *pdata,
		       const struct ftreenode *root )
{
  FILE *LOG = unur_get_stream();
  CHECK_NULL(pdata,RETURN_VOID);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,RETURN_VOID);
  CHECK_NULL(root,RETURN_VOID);   COOKIE_CHECK(root,CK_FSTR_TNODE,RETURN_VOID);
  fprintf(LOG,"%s: parse tree:  (left nodes above right nodes)\n",GENTYPE); 
  fprintf(LOG,"%s:\n",GENTYPE);
  _unur_fstr_debug_show_tree(pdata,root,0,0);
  fprintf(LOG,"%s:\n",GENTYPE);
} 
void
_unur_fstr_debug_show_tree(const struct parser_data *pdata,
			   const struct ftreenode *node,
			   int level, int location)
{ 
  FILE *LOG = unur_get_stream();
  const char *name;
  int i, mask; 
  CHECK_NULL(pdata,RETURN_VOID);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,RETURN_VOID);
  fprintf(LOG,"%s: ",GENTYPE); 
  for (i = 0, mask = 1; i < level; i++, mask <<= 1) 
    if (mask & location) 
      fprintf(LOG,"|   "); 
    else 
      fprintf(LOG,"    "); 
  if( node != NULL ) {
    COOKIE_CHECK(node,CK_FSTR_TNODE,RETURN_VOID);
    (mask & location) ? fprintf(LOG,"+--") : fprintf(LOG,"\\__");
    switch (node->type) {
    case S_SCONST:
      fprintf(LOG,"'%s'\t(const=%g)", node->symbol,node->val);  break;
    case S_UCONST:
      fprintf(LOG,"'%g'\t(const)", node->val);  break;
    case S_UIDENT:
      name = (pdata) ? pdata->variable_name : "x";
      fprintf(LOG,"'%s'\t(variable)", name);  break;
    case S_UFUNCT:
      name = (pdata) ? pdata->function_name : "f";
      fprintf(LOG,"'%s'\t(user function)", name);  break;
    default:
      fprintf(LOG,"'%s'", node->symbol);
    }
    fprintf(LOG,"\n");
    if ( node->left || node->right) {
      _unur_fstr_debug_show_tree(pdata,node->left, level+1,location|(mask<<1)); 
      _unur_fstr_debug_show_tree(pdata,node->right,level+1,location); 
    }
  }
  else {  
    (mask & location) ? fprintf(LOG,"+--") : fprintf(LOG,"\\__");
    fprintf(LOG,"(void)\n"); 
  } 
} 
void
_unur_fstr_debug_deriv (const struct ftreenode *funct, const struct ftreenode *deriv)
{
  FILE *LOG = unur_get_stream();
  char *str;
  CHECK_NULL(funct,RETURN_VOID);  COOKIE_CHECK(funct,CK_FSTR_TNODE,RETURN_VOID);
  CHECK_NULL(deriv,RETURN_VOID);  COOKIE_CHECK(deriv,CK_FSTR_TNODE,RETURN_VOID);
  fprintf(LOG,"%s: Derivative df/dx of \n",GENTYPE);
  str = _unur_fstr_tree2string(funct,"x","f",TRUE);
  fprintf(LOG,"%s:  f(x) = %s\n",GENTYPE,str);
  free (str);
  if (deriv) {
    str = _unur_fstr_tree2string(deriv,"x","df",TRUE);
    fprintf(LOG,"%s:  f'(x) = %s\n",GENTYPE,str);
    free (str);
  }
  else {
    fprintf(LOG,"%s:  f'(x) = (unknown)\n",GENTYPE);
  }
  fprintf(LOG,"%s:\n",GENTYPE);
} 
#endif   
