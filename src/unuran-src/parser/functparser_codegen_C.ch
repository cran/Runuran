/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

enum {
  C_FUNCT_ERROR = 0x10000000u,      
  C_FUNCT_SGN   = 0x00000001u,      
  C_FUNCT_SEC   = 0x00000002u       
};
int 
_unur_fstr_tree2C ( FILE *out, const struct ftreenode *root,
		    const char *variable, const char *funct_name )
{
  struct unur_string output = {NULL, 0, 0};
  unsigned rcode;
  _unur_check_NULL( GENTYPE, root, UNUR_ERR_NULL );
  _unur_check_NULL( GENTYPE, symbol[root->token].node2C, UNUR_ERR_NULL );
  rcode = symbol[root->token].node2C (&output,root,variable);
  if (rcode & C_FUNCT_ERROR) { 
    if (output.text) free(output.text);
    return UNUR_ERR_GEN_DATA;
  }
  _unur_fstr_C_specfunct (out,rcode);
  fprintf (out,"static double %s (double %s)\n",funct_name,variable );
  fprintf (out,"{\n\treturn (%s);\n}\n",output.text);
  free(output.text);
  return UNUR_SUCCESS;
} 
unsigned
C_error ( struct unur_string *output ATTRIBUTE__UNUSED, 
	  const struct ftreenode *node  ATTRIBUTE__UNUSED,
	  const char *variable ATTRIBUTE__UNUSED )
{
  _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
  return C_FUNCT_ERROR;
} 
unsigned
C_const ( struct unur_string *output, const struct ftreenode *node, 
	  const char *variable ATTRIBUTE__UNUSED )
{
  _unur_fstr_print_C( output, NULL, node->val );
  return 0u;
} 
unsigned
C_var ( struct unur_string *output,
	const struct ftreenode *node ATTRIBUTE__UNUSED, const char *variable )
{
  _unur_fstr_print_C( output, variable, 0. );
  return 0u;
} 
unsigned
C_prefix_generic ( struct unur_string *output, const struct ftreenode *node,
		   const char *variable, const char *symb )
{
  unsigned rcode = 0u;
  struct ftreenode *left  = node->left;    
  struct ftreenode *right = node->right;   
  _unur_fstr_print_C( output, symb, 0. );
  _unur_fstr_print_C( output, "(", 0. );
  if (left) {
    rcode |= symbol[left->token].node2C (output,left,variable);
    _unur_fstr_print_C( output, ",", 0. );
  }
  if (right) {
    rcode |= symbol[right->token].node2C (output,right,variable);
  }
  _unur_fstr_print_C( output, ")", 0. );
  return rcode;
} 
unsigned
C_prefix ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return C_prefix_generic(output,node,variable,symbol[node->token].name);
} 
unsigned
C_power ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return C_prefix_generic(output,node,variable,"pow");
} 
unsigned
C_sec ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return (C_FUNCT_SEC | C_prefix_generic(output,node,variable,"_acg_sec"));
} 
unsigned
C_abs ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return C_prefix_generic(output,node,variable,"fabs");
} 
unsigned
C_sgn ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return (C_FUNCT_SGN | C_prefix_generic(output,node,variable,"_acg_sgn"));
} 
unsigned
C_infix_generic ( struct unur_string *output, const struct ftreenode *node,
		  const char *variable, const char *symb )
{
  unsigned rcode = 0u;
  struct ftreenode *left  = node->left;    
  struct ftreenode *right = node->right;   
  if (left==NULL || right==NULL)
    return C_FUNCT_ERROR;
  _unur_fstr_print_C( output, "(", 0. );
  rcode |= symbol[left->token].node2C (output,left,variable);
  _unur_fstr_print_C( output, symb, 0. );
  rcode |= symbol[right->token].node2C (output,right,variable);
  _unur_fstr_print_C( output, ")", 0. );
  return rcode;
} 
unsigned
C_infix ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return C_infix_generic(output,node,variable,symbol[node->token].name);
} 
unsigned
C_equal ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return C_infix_generic(output,node,variable,"==");
} 
unsigned
C_unequal ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return C_infix_generic(output,node,variable,"!=");
} 
unsigned
C_minus ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  unsigned rcode = 0u;
  struct ftreenode *left  = node->left;    
  struct ftreenode *right = node->right;   
  if (left==NULL || right==NULL)
    return C_FUNCT_ERROR;
  _unur_fstr_print_C( output, "(", 0. );
  if (!(left->type == S_UCONST && _unur_iszero(left->val)))
    rcode |= symbol[left->token].node2C (output,left,variable);
  _unur_fstr_print_C( output, "-", 0. );
  rcode |= symbol[right->token].node2C (output,right,variable);
  _unur_fstr_print_C( output, ")", 0. );
  return rcode;
} 
unsigned
C_mod ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return C_infix_generic(output,node,variable,"%");
} 
int
_unur_fstr_C_specfunct ( FILE *out, unsigned flags )
{
  if (flags & C_FUNCT_SGN) {
    _unur_fstr_C_sgn(out);
  }
  if (flags & C_FUNCT_SEC) {
    _unur_fstr_C_sec(out);
  }
  return UNUR_SUCCESS;
} 
int
_unur_fstr_C_sgn ( FILE *out )
{
  fprintf(out,"#ifndef _ACG_FUNCT_SGN\n");
  fprintf(out,"#define _ACG_FUNCT_SGN\n");
  fprintf(out,"static double _acg_sgn(double x)\n{\n");
  fprintf(out,"\treturn ((x<0.) ? -1. : ((x>0.) ? 1. : 0.));\n");
  fprintf(out,"}\n");
  fprintf(out,"#endif \n\n");
  return UNUR_SUCCESS;
} 
int
_unur_fstr_C_sec ( FILE *out )
{
  fprintf(out,"#ifndef _ACG_FUNCT_SEC\n");
  fprintf(out,"#define _ACG_FUNCT_SEC\n");
  fprintf(out,"static double _acg_sec(double x)\n{\n");
  fprintf(out,"\tdouble cosx = cos(x);\n");
  fprintf(out,"\treturn ((cosx == 0.) ? HUGE_VAL : 1./cosx) ;\n");
  fprintf(out,"}\n");
  fprintf(out,"#endif \n\n");
  return UNUR_SUCCESS;
} 
int
_unur_fstr_print_C ( struct unur_string *output, const char *symb, double number )
{
  if (symb)
    _unur_string_appendtext( output, symb );
  else
    _unur_string_append( output, "%.20e", number);
  return UNUR_SUCCESS;
} 
