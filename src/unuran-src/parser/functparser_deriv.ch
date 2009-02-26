/* Copyright (c) 2000-2009 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct ftreenode *
_unur_fstr_make_derivative ( const struct ftreenode *root )
{
  struct ftreenode *deriv = NULL;  
  int error = 0;                   
  _unur_check_NULL( GENTYPE,root,NULL );
  COOKIE_CHECK(root,CK_FSTR_TNODE,NULL);
  deriv = (root) ? (*symbol[root->token].dcalc)(root,&error) : NULL;
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag)
    _unur_fstr_debug_deriv(root,deriv);
#endif
  if (error == TRUE) {
    unur_errno = UNUR_ERR_FSTR_DERIV;
    if (deriv) _unur_fstr_free(deriv);
    return NULL;
  }
  return deriv;
} 
struct ftreenode *
d_error (const struct ftreenode *node, int *error)
{
  CHECK_NULL(node,NULL);  COOKIE_CHECK(node,CK_FSTR_TNODE,NULL);
  _unur_fstr_error_deriv(node,__LINE__);
  *error = TRUE;
  return NULL;
} 
struct ftreenode *
d_const (const struct ftreenode *node ATTRIBUTE__UNUSED, int *error ATTRIBUTE__UNUSED)
{
  CHECK_NULL(node,NULL);  COOKIE_CHECK(node,CK_FSTR_TNODE,NULL);
  return _unur_fstr_create_node(NULL,0.,s_uconst,NULL,NULL);
} 
struct ftreenode *
d_var (const struct ftreenode *node ATTRIBUTE__UNUSED, int *error ATTRIBUTE__UNUSED)
{
  CHECK_NULL(node,NULL);  COOKIE_CHECK(node,CK_FSTR_TNODE,NULL);
  return _unur_fstr_create_node(NULL,1.,s_uconst,NULL,NULL);
} 
struct ftreenode *
d_add (const struct ftreenode *node, int *error)
{
  struct ftreenode *left, *right;
  struct ftreenode *d_left, *d_right;
  int op;
  CHECK_NULL(node,NULL);  COOKIE_CHECK(node,CK_FSTR_TNODE,NULL);
  left = node->left;
  right = node->right;
  op = node->token;
  d_left  = (left)  ? (*symbol[left->token].dcalc) (left,error)  : NULL;
  d_right = (right) ? (*symbol[right->token].dcalc)(right,error) : NULL;
  return _unur_fstr_create_node(node->symbol,0.,op,d_left,d_right);
} 
struct ftreenode *
d_mul (const struct ftreenode *node, int *error)
{
  struct ftreenode *left, *right;
  struct ftreenode *d_left, *d_right;
  struct ftreenode *br_left, *br_right;
  CHECK_NULL(node,NULL);  COOKIE_CHECK(node,CK_FSTR_TNODE,NULL);
  left  = _unur_fstr_dup_tree(node->left);
  right = _unur_fstr_dup_tree(node->right);
  d_left  = (left)  ? (*symbol[left->token].dcalc) (left,error)  : NULL;
  d_right = (right) ? (*symbol[right->token].dcalc)(right,error) : NULL;
  br_left  = _unur_fstr_create_node("*",0.,s_mul,d_left,right);
  br_right = _unur_fstr_create_node("*",0.,s_mul,left,d_right);
  return _unur_fstr_create_node("+",0.,s_plus,br_left,br_right);
} 
struct ftreenode *
d_div (const struct ftreenode *node, int *error)
{
  struct ftreenode *left, *right;
  struct ftreenode *d_left, *d_right;
  struct ftreenode *br_left, *br_right, *two;
  struct ftreenode *numerator, *denominator; 
  CHECK_NULL(node,NULL);  COOKIE_CHECK(node,CK_FSTR_TNODE,NULL);
  left  = _unur_fstr_dup_tree(node->left);
  right = _unur_fstr_dup_tree(node->right);
  d_left  = (left)  ? (*symbol[left->token].dcalc) (left,error)  : NULL;
  d_right = (right) ? (*symbol[right->token].dcalc)(right,error) : NULL;
  two = _unur_fstr_create_node(NULL,2.,s_uconst,NULL,NULL);   
  denominator = _unur_fstr_create_node("^",0.,s_power,right,two);
  right = _unur_fstr_dup_tree(node->right);    
  br_left  = _unur_fstr_create_node("*",0.,s_mul,d_left,right);
  br_right = _unur_fstr_create_node("*",0.,s_mul,left,d_right);
  numerator= _unur_fstr_create_node("-",0.,s_minus,br_left,br_right);
  return _unur_fstr_create_node("/",0.,s_div,numerator,denominator);
} 
struct ftreenode *
d_power (const struct ftreenode *node, int *error)
{
  struct ftreenode *left, *right;
  struct ftreenode *d_left, *d_right;
  struct ftreenode *br_right;
  struct ftreenode *dup_node, *tmp1, *tmp2;
  CHECK_NULL(node,NULL);  COOKIE_CHECK(node,CK_FSTR_TNODE,NULL);
  left = node->left;
  right = node->right;
  if (right && (right->type == S_UCONST || right->type == S_SCONST) ) {
    d_left  = (left)  ? (*symbol[left->token].dcalc) (left,error)  : NULL;
    left  = _unur_fstr_dup_tree(node->left);
    right = _unur_fstr_dup_tree(node->right);
    tmp1     = _unur_fstr_create_node(NULL,right->val-1,s_uconst,NULL,NULL);
    tmp2     = _unur_fstr_create_node("^",0.,s_power,left,tmp1);
    br_right = _unur_fstr_create_node("*",0.,s_mul,right,tmp2);
    return _unur_fstr_create_node("*",0.,s_mul,d_left,br_right);
  }
  else if (left && (left->type == S_UCONST || left->type == S_SCONST) ) {
    int s_log = _unur_fstr_find_symbol("log",_ans_start,_ans_end);
    d_right = (right) ? (*symbol[right->token].dcalc) (right,error)  : NULL;
    left = _unur_fstr_dup_tree(node->left);
    dup_node = _unur_fstr_dup_tree(node);
    tmp1     = _unur_fstr_create_node("log",0.,s_log,NULL,left);
    br_right = _unur_fstr_create_node("*",0.,s_mul,tmp1,dup_node);
    return _unur_fstr_create_node("*",0.,s_mul,d_right,br_right);
  }
  else {
    _unur_fstr_error_deriv(node,__LINE__);
    *error = TRUE;
    return NULL;
  }
} 
struct ftreenode *
d_exp (const struct ftreenode *node, int *error)
{
  struct ftreenode *right;
  struct ftreenode *d_right;
  struct ftreenode *br_right;
  CHECK_NULL(node,NULL);  COOKIE_CHECK(node,CK_FSTR_TNODE,NULL);
  right = node->right;
  d_right = (right) ? (*symbol[right->token].dcalc)(right,error) : NULL;
  br_right = _unur_fstr_dup_tree(node);
  return _unur_fstr_create_node("*",0.,s_mul,d_right,br_right);
} 
struct ftreenode *
d_log (const struct ftreenode *node, int *error)
{
  struct ftreenode *right;
  struct ftreenode *d_right;
  CHECK_NULL(node,NULL);  COOKIE_CHECK(node,CK_FSTR_TNODE,NULL);
  right = _unur_fstr_dup_tree(node->right);
  d_right = (right) ? (*symbol[right->token].dcalc)(right,error) : NULL;
  return _unur_fstr_create_node("/",0.,s_div,d_right,right);
} 
struct ftreenode *
d_sin (const struct ftreenode *node, int *error)
{
  struct ftreenode *right;
  struct ftreenode *d_right;
  struct ftreenode *br_right;
  int s_cos = _unur_fstr_find_symbol("cos",_ans_start,_ans_end);
  CHECK_NULL(node,NULL);  COOKIE_CHECK(node,CK_FSTR_TNODE,NULL);
  right = _unur_fstr_dup_tree(node->right);
  d_right = (right) ? (*symbol[right->token].dcalc)(right,error) : NULL;
  br_right = _unur_fstr_create_node("cos",0.,s_cos,NULL,right);
  return _unur_fstr_create_node(NULL,0.,s_mul,d_right,br_right);
} 
struct ftreenode *
d_cos (const struct ftreenode *node, int *error)
{
  struct ftreenode *right;
  struct ftreenode *d_right;
  struct ftreenode *br_left, *br_right;
  struct ftreenode *zero;
  int s_sin = _unur_fstr_find_symbol("sin",_ans_start,_ans_end);
  CHECK_NULL(node,NULL);  COOKIE_CHECK(node,CK_FSTR_TNODE,NULL);
  right = _unur_fstr_dup_tree(node->right);
  d_right = (right) ? (*symbol[right->token].dcalc)(right,error) : NULL;
  br_right = _unur_fstr_create_node("sin",0.,s_sin,NULL,right);
  zero = _unur_fstr_create_node(NULL,0.,s_uconst,NULL,NULL);
  br_left  = _unur_fstr_create_node("-",0.,s_minus,zero,d_right);
  return _unur_fstr_create_node("*",0.,s_mul,br_left,br_right);
} 
struct ftreenode *
d_tan (const struct ftreenode *node, int *error)
{
  struct ftreenode *right;
  struct ftreenode *d_right;
  struct ftreenode *br_right, *sub_right;
  struct ftreenode *two;
  int s_sec = _unur_fstr_find_symbol("sec",_ans_start,_ans_end);
  CHECK_NULL(node,NULL);  COOKIE_CHECK(node,CK_FSTR_TNODE,NULL);
  right = _unur_fstr_dup_tree(node->right);
  d_right = (right) ? (*symbol[right->token].dcalc)(right,error) : NULL;
  two = _unur_fstr_create_node(NULL,2.,s_uconst,NULL,NULL);   
  sub_right = _unur_fstr_create_node("sec",0.,s_sec,NULL,right);
  br_right = _unur_fstr_create_node("^",0.,s_power,sub_right,two);
  return _unur_fstr_create_node("*",0.,s_mul,d_right,br_right);
} 
struct ftreenode *
d_sec (const struct ftreenode *node, int *error)
{
  struct ftreenode *right;
  struct ftreenode *d_right;
  struct ftreenode *br_right, *sub_right, *dup_node;
  int s_tan = _unur_fstr_find_symbol("tan",_ans_start,_ans_end);
  CHECK_NULL(node,NULL);  COOKIE_CHECK(node,CK_FSTR_TNODE,NULL);
  right = _unur_fstr_dup_tree(node->right);
  d_right = (right) ? (*symbol[right->token].dcalc)(right,error) : NULL;
  dup_node = _unur_fstr_dup_tree(node);
  sub_right = _unur_fstr_create_node("tan",0.,s_tan,NULL,right);
  br_right = _unur_fstr_create_node("*",0.,s_mul,sub_right,dup_node);
  return _unur_fstr_create_node("*",0.,s_mul,d_right,br_right);
} 
struct ftreenode *
d_sqrt (const struct ftreenode *node, int *error)
{
  struct ftreenode *right;
  struct ftreenode *d_right;
  struct ftreenode *br_right;
  struct ftreenode *two, *dup_tree;
  CHECK_NULL(node,NULL);  COOKIE_CHECK(node,CK_FSTR_TNODE,NULL);
  right = node->right;
  d_right = (right) ? (*symbol[right->token].dcalc)(right,error) : NULL;
  dup_tree = _unur_fstr_dup_tree(node);
  two = _unur_fstr_create_node(NULL,2.,s_uconst,NULL,NULL);   
  br_right = _unur_fstr_create_node("*",0.,s_mul,two,dup_tree);
  return _unur_fstr_create_node("/",0.,s_div,d_right,br_right);
} 
struct ftreenode *
d_abs (const struct ftreenode *node, int *error)
{
  struct ftreenode *right;
  struct ftreenode *d_right;
  struct ftreenode *br_right;
  int s_sgn = _unur_fstr_find_symbol("sgn",_ans_start,_ans_end);
  CHECK_NULL(node,NULL);  COOKIE_CHECK(node,CK_FSTR_TNODE,NULL);
  right = _unur_fstr_dup_tree(node->right);
  d_right = (right) ? (*symbol[right->token].dcalc)(right,error) : NULL;
  br_right = _unur_fstr_create_node("sgn",0.,s_sgn,NULL,right);
  return _unur_fstr_create_node("*",0.,s_mul,d_right,br_right);
} 
void
_unur_fstr_error_deriv (const struct ftreenode *node, int line)
{
  struct unur_string *reason;
  CHECK_NULL(node,RETURN_VOID);  COOKIE_CHECK(node,CK_FSTR_TNODE,RETURN_VOID);
  reason = _unur_string_new();
  _unur_string_append( reason, "cannot derivate subtree at '%s'", node->symbol);
  _unur_error_x( GENTYPE, __FILE__, line, "error", UNUR_ERR_FSTR_DERIV,reason->text);
  _unur_string_free( reason );
#ifdef UNUR_ENABLE_LOGGING
  _unur_fstr_debug_tree(NULL,node);
  _unur_log_printf_simple ( "%s:\n",GENTYPE );
#endif  
} 
