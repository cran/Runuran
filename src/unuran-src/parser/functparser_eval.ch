/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

double
_unur_fstr_eval_tree (const struct ftreenode *root, double x)
{  
  CHECK_NULL(root,INFINITY);   COOKIE_CHECK(root,CK_FSTR_TNODE,INFINITY);
  return _unur_fstr_eval_node( root, x );
} 
#define CHECK_INF(x)    if(_unur_FP_is_infinity((x))) return INFINITY;
#define CHECK_INFS(l,r) do { CHECK_INF((l)); CHECK_INF((r)); } while(0)
#define NODE_ARGS  double l ATTRIBUTE__UNUSED, double r ATTRIBUTE__UNUSED
double v_dummy  (NODE_ARGS) { return 0.; }
double v_const  (NODE_ARGS) { return 0.; }  
double v_less   (NODE_ARGS) { return (double)(l <  r); }
double v_equal  (NODE_ARGS) { return (double)(_unur_FP_same(l,r)); }
double v_greater(NODE_ARGS) { return (double)(l >  r); }
double v_less_or(NODE_ARGS) { return (double)(l <= r); }
double v_unequal(NODE_ARGS) { return (double)(!_unur_FP_same(l,r)); }
double v_grtr_or(NODE_ARGS) { return (double)(l >= r); }
double v_plus   (NODE_ARGS) { return (l + r); }
double v_minus  (NODE_ARGS) { return (l - r); }
double v_mul    (NODE_ARGS) { return (l * r); }
double v_div    (NODE_ARGS) { return (l / r); }
double v_power  (NODE_ARGS) { return pow(l,r); }
double v_mod    (NODE_ARGS) { return (double)((int)l % (int)r); }
double v_exp    (NODE_ARGS) { return exp(r); }
double v_log    (NODE_ARGS) { return (r<=0.) ? INFINITY : log(r); }
double v_sin    (NODE_ARGS) { CHECK_INF(r); return sin(r); }
double v_cos    (NODE_ARGS) { CHECK_INF(r); return cos(r); }
double v_tan    (NODE_ARGS) { CHECK_INF(r); return tan(r); }
double v_sec    (NODE_ARGS) { double cosr; CHECK_INF(r); cosr=cos(r); 
                                               return _unur_iszero(cosr) ? INFINITY : 1./cosr; }
double v_sqrt   (NODE_ARGS) { return (r<0.) ? INFINITY : sqrt(r); }
double v_abs    (NODE_ARGS) { return fabs(r); }
double v_sgn    (NODE_ARGS) { return ((r<0.) ? -1. : ((r>0.) ? 1. : 0.)); }
#undef CHECK_INF
#undef CHECK_INFS
#undef NODE_ARGS
double
_unur_fstr_eval_node (const struct ftreenode *node, double x)
{
  double val_l, val_r;
  CHECK_NULL(node,INFINITY);   COOKIE_CHECK(node,CK_FSTR_TNODE,INFINITY);
  switch (node->type) {
  case S_UCONST:
  case S_SCONST:
    return node->val;
  case S_UIDENT:
    return x;
  default:
    val_l = (node->left)  ? _unur_fstr_eval_node(node->left, x) : 0. ;
    val_r = (node->right) ? _unur_fstr_eval_node(node->right,x) : 0. ;
    return (*symbol[node->token].vcalc)(val_l,val_r);
  }
} 
