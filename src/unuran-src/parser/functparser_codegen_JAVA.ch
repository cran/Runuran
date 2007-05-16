/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

enum {
  J_FUNCT_ERROR = 0x10000000u,      
  J_FUNCT_SGN   = 0x00000001u,      
  J_FUNCT_SEC   = 0x00000002u,      
  J_FUNCT_LT    = 0x00000010u,      
  J_FUNCT_LE    = 0x00000020u,      
  J_FUNCT_GT    = 0x00000040u,      
  J_FUNCT_GE    = 0x00000080u,      
  J_FUNCT_EQ    = 0x00000100u,      
  J_FUNCT_NE    = 0x00000200u       
};
#define _unur_fstr_print_J(output,symb,number)   _unur_fstr_print_C((output),(symb),(number)) 
int 
_unur_fstr_tree2JAVA ( FILE *out, const struct ftreenode *root,
		       const char *variable, const char *funct_name )
{
  struct unur_string output = {NULL, 0, 0};
  unsigned rcode;
  _unur_check_NULL( GENTYPE, root, UNUR_ERR_NULL );
  _unur_check_NULL( GENTYPE, symbol[root->token].node2J, UNUR_ERR_NULL );
  rcode = symbol[root->token].node2J (&output,root,variable);
  if (rcode & J_FUNCT_ERROR) { 
    if (output.text) free(output.text);
    return UNUR_ERR_GEN_DATA;
  }
  _unur_fstr_J_specfunct (out,rcode);
  fprintf (out,"\tstatic double %s (double %s)\n",funct_name,variable );
  fprintf (out,"\t{\n\t\treturn (%s);\n\t}\n",output.text);
  free(output.text);
  return UNUR_SUCCESS;
} 
unsigned
J_error ( struct unur_string *output ATTRIBUTE__UNUSED,
	  const struct ftreenode *node ATTRIBUTE__UNUSED,
	  const char *variable ATTRIBUTE__UNUSED )
{
  _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
  return J_FUNCT_ERROR;
} 
unsigned
J_const ( struct unur_string *output, 
	  const struct ftreenode *node, const char *variable ATTRIBUTE__UNUSED )
{
  _unur_fstr_print_J( output, NULL, node->val );
  return 0u;
} 
unsigned
J_var ( struct unur_string *output, 
	const struct ftreenode *node ATTRIBUTE__UNUSED, const char *variable )
{
  _unur_fstr_print_J( output, variable, 0. );
  return 0u;
} 
unsigned
J_prefix_generic ( struct unur_string *output, const struct ftreenode *node,
		   const char *variable, const char *symb )
{
  unsigned rcode = 0u;
  struct ftreenode *left  = node->left;    
  struct ftreenode *right = node->right;   
  _unur_fstr_print_J( output, symb, 0. );
  _unur_fstr_print_J( output, "(", 0. );
  if (left) {
    rcode |= symbol[left->token].node2J (output,left,variable);
    _unur_fstr_print_J( output, ",", 0. );
  }
  if (right) {
    rcode |= symbol[right->token].node2J (output,right,variable);
  }
  _unur_fstr_print_J( output, ")", 0. );
  return rcode;
} 
unsigned
J_prefix ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  _unur_fstr_print_J( output, "Math.", 0. );
  return J_prefix_generic(output,node,variable,symbol[node->token].name);
} 
unsigned
J_lt ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return (J_FUNCT_LT | J_prefix_generic(output,node,variable,"RelLT"));
} 
unsigned
J_le ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return (J_FUNCT_LE | J_prefix_generic(output,node,variable,"RelLE"));
} 
unsigned
J_gt ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return (J_FUNCT_GT | J_prefix_generic(output,node,variable,"RelGT"));
} 
unsigned
J_ge ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return (J_FUNCT_GE | J_prefix_generic(output,node,variable,"RelGE"));
} 
unsigned
J_eq ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return (J_FUNCT_EQ | J_prefix_generic(output,node,variable,"RelEQ"));
} 
unsigned
J_ne ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return (J_FUNCT_NE | J_prefix_generic(output,node,variable,"RelNE"));
} 
unsigned
J_power ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return J_prefix_generic(output,node,variable,"Math.pow");
} 
unsigned
J_sec ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return (J_FUNCT_SEC | J_prefix_generic(output,node,variable,"sec"));
} 
unsigned
J_sgn ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return (J_FUNCT_SGN | J_prefix_generic(output,node,variable,"sgn"));
} 
unsigned
J_infix_generic ( struct unur_string *output, const struct ftreenode *node,
		  const char *variable, const char *symb )
{
  unsigned rcode = 0u;
  struct ftreenode *left  = node->left;    
  struct ftreenode *right = node->right;   
  if (left==NULL || right==NULL)
    return J_FUNCT_ERROR;
  _unur_fstr_print_J( output, "(", 0. );
  rcode |= symbol[left->token].node2J (output,left,variable);
  _unur_fstr_print_J( output, symb, 0. );
  rcode |= symbol[right->token].node2J (output,right,variable);
  _unur_fstr_print_J( output, ")", 0. );
  return rcode;
} 
unsigned
J_infix ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return J_infix_generic(output,node,variable,symbol[node->token].name);
} 
unsigned
J_minus ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  unsigned rcode = 0u;
  struct ftreenode *left  = node->left;    
  struct ftreenode *right = node->right;   
  if (left==NULL || right==NULL)
    return J_FUNCT_ERROR;
  _unur_fstr_print_J( output, "(", 0. );
  if (!(left->type == S_UCONST && _unur_iszero(left->val)))
    rcode |= symbol[left->token].node2J (output,left,variable);
  _unur_fstr_print_J( output, "-", 0. );
  rcode |= symbol[right->token].node2J (output,right,variable);
  _unur_fstr_print_J( output, ")", 0. );
  return rcode;
} 
unsigned
J_mod ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return J_infix_generic(output,node,variable,"%");
} 
int
_unur_fstr_J_specfunct ( FILE *out, unsigned flags )
{
  if (flags & J_FUNCT_SGN) {
    _unur_fstr_J_sgn(out);
  }
  if (flags & J_FUNCT_SEC) {
    fprintf(out,"\tstatic double sec (double x) { ");
    fprintf(out,"return (1./Math.cos(x)); }\n\n");
  }
  if (flags & J_FUNCT_LE) {
    fprintf(out,"\tstatic double RelLE (double x, double y) { ");
    fprintf(out,"return ((x<=y) ? 1. : 0.); }\n\n");
  }
  if (flags & J_FUNCT_GE) {
    fprintf(out,"\tstatic double RelGE (double x, double y) { ");
    fprintf(out,"return ((x>=y) ? 1. : 0.); }\n\n");
  }
  if (flags & J_FUNCT_LT) {
    fprintf(out,"\tstatic double RelLT (double x, double y) { ");
    fprintf(out,"return ((x<y) ? 1. : 0.); }\n\n");
  }
  if (flags & J_FUNCT_GT) {
    fprintf(out,"\tstatic double RelGT (double x, double y) { ");
    fprintf(out,"return ((x>y) ? 1. : 0.); }\n\n");
  }
  if (flags & J_FUNCT_EQ) {
    fprintf(out,"\tstatic double RelEQ (double x, double y) { ");
    fprintf(out,"return ((x==y) ? 1. : 0.); }\n\n");
  }
  if (flags & J_FUNCT_NE) {
    fprintf(out,"\tstatic double RelNE (double x, double y) { ");
    fprintf(out,"return ((x!=y) ? 1. : 0.); }\n\n");
  }
  return UNUR_SUCCESS;
} 
int
_unur_fstr_J_sgn ( FILE *out )
{
  fprintf(out,"\tstatic double sgn (double x)\n\t{\n");
  fprintf(out,"\t\tif (x<0.)  return -1.;\n");
  fprintf(out,"\t\tif (x>0.)  return  1.;\n");
  fprintf(out,"\t\t return  0.;\n");
  fprintf(out,"\t}\n\n");
  return UNUR_SUCCESS;
} 
