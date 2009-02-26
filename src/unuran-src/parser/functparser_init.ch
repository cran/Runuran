/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct ftreenode *
_unur_fstr2tree (const char *functstr)
{
  return _unur_fstr_2_tree( functstr, FALSE );
} 
struct ftreenode *
_unur_fstr2tree_DefFunct (const char *functstr)
{
  return _unur_fstr_2_tree( functstr, TRUE );
} 
struct ftreenode *
_unur_fstr_dup_tree (const struct ftreenode *root)
{
  struct ftreenode *dup;
  if (root==NULL) return NULL;
  COOKIE_CHECK(root,CK_FSTR_TNODE,NULL);
  dup = _unur_xmalloc(sizeof(struct ftreenode));
  memcpy(dup,root,sizeof(struct ftreenode));
  if (root->left)  dup->left  = _unur_fstr_dup_tree(root->left);
  if (root->right) dup->right = _unur_fstr_dup_tree(root->right);
  return dup;
} 
void
_unur_fstr_free (struct ftreenode *root)  
{ 
  if( root != NULL ) {
    COOKIE_CHECK(root,CK_FSTR_TNODE,RETURN_VOID);
    if (root->left)  _unur_fstr_free(root->left);
    if (root->right) _unur_fstr_free(root->right);
    free(root); 
  } 
} 
struct ftreenode *
_unur_fstr_2_tree (const char *functstr, int withDefFunct)
{ 
  struct parser_data *pdata;
  struct ftreenode *root;
  _unur_check_NULL( GENTYPE,functstr,NULL );
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag)
    _unur_fstr_debug_input(functstr);
#endif
  pdata = _unur_fstr_parser_init(functstr);
  if (pdata == NULL)
    return NULL;
  _unur_fstr_tokenize(pdata);
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag)
    _unur_fstr_debug_token(pdata);
#endif
  if (pdata->perrno) {
    _unur_fstr_parser_free(pdata);
    return NULL;
  }
  if (withDefFunct) {
    struct ftreenode *tmp = _unur_FunctDefinition(pdata);
    root = tmp->right;
    _unur_fstr_free(tmp->left);
    free(tmp);
  }
  else {
    root = _unur_Expression(pdata);
  }
#ifdef UNUR_ENABLE_LOGGING
  if ((_unur_default_debugflag & UNUR_DEBUG_SETUP) && root)
    _unur_fstr_debug_tree(pdata,root);
#endif
  if (pdata->tno < pdata->n_tokens && !pdata->perrno)
    _unur_fstr_error_parse(pdata,ERR_UNFINISHED,__LINE__); 
  if (pdata->perrno) {
    _unur_fstr_parser_free(pdata);
    _unur_fstr_free(root);
    return NULL;
  }
  _unur_fstr_parser_free(pdata);
  return root; 
} 
