/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

enum {
  F_FUNCT_ERROR = 0x10000000u,      
  F_FUNCT_SGN   = 0x00000001u,      
  F_FUNCT_SEC   = 0x00000002u,      
  F_FUNCT_LT    = 0x00000010u,      
  F_FUNCT_LE    = 0x00000020u,      
  F_FUNCT_GT    = 0x00000040u,      
  F_FUNCT_GE    = 0x00000080u,      
  F_FUNCT_EQ    = 0x00000100u,      
  F_FUNCT_NE    = 0x00000200u       
};
int 
_unur_fstr_tree2FORTRAN ( FILE *out, const struct ftreenode *root,
			  const char *variable, const char *funct_name )
{
  struct unur_string output = {NULL, 0, 0};
  unsigned rcode;
  int line;
  _unur_check_NULL( GENTYPE, root, UNUR_ERR_NULL );
  _unur_check_NULL( GENTYPE, symbol[root->token].node2F, UNUR_ERR_NULL );
  rcode = symbol[root->token].node2F (&output,root,variable);
  if (rcode & F_FUNCT_ERROR) { 
    if (output.text) free(output.text);
    return UNUR_ERR_GEN_DATA;
  }
  fprintf (out,"      DOUBLE PRECISION FUNCTION %.6s(x)\n\n", funct_name);
  fprintf (out,"      IMPLICIT DOUBLE PRECISION (A-Z)\n");
  _unur_fstr_F_specfunct (out,rcode);
  fprintf (out,"C\n");
  fprintf (out,"C     compute PDF\n");
  fprintf (out,"C\n");
  fprintf (out,"      %.6s = \n", funct_name);
  for (line = 0; line < (output.length-1)/60 + 1; line++) {
    fprintf (out,"     $   %.60s\n", output.text+60*line);
  }
  fprintf (out,"      RETURN\n");
  fprintf (out,"\n");
  fprintf (out,"      END\n");
  fprintf (out,"\n");
  free(output.text);
  return UNUR_SUCCESS;
} 
unsigned
F_error ( struct unur_string *output ATTRIBUTE__UNUSED,
	  const struct ftreenode *node ATTRIBUTE__UNUSED,
	  const char *variable ATTRIBUTE__UNUSED )
{
  _unur_error(GENTYPE,UNUR_ERR_SHOULD_NOT_HAPPEN,"");
  return F_FUNCT_ERROR;
} 
unsigned
F_const ( struct unur_string *output, const struct ftreenode *node,
	  const char *variable ATTRIBUTE__UNUSED )
{
  _unur_fstr_print_F( output, NULL, node->val );
  return 0u;
} 
unsigned
F_var ( struct unur_string *output,
	const struct ftreenode *node ATTRIBUTE__UNUSED, const char *variable )
{
  _unur_fstr_print_F( output, variable, 0. );
  return 0u;
} 
unsigned
F_prefix_generic ( struct unur_string *output, const struct ftreenode *node,
		   const char *variable, const char *symb )
{
  unsigned rcode = 0u;
  struct ftreenode *left  = node->left;    
  struct ftreenode *right = node->right;   
  _unur_fstr_print_F( output, symb, 0. );
  _unur_fstr_print_F( output, "(", 0. );
  if (left) {
    rcode |= symbol[left->token].node2F (output,left,variable);
    _unur_fstr_print_F( output, ",", 0. );
  }
  if (right) {
    rcode |= symbol[right->token].node2F (output,right,variable);
  }
  _unur_fstr_print_F( output, ")", 0. );
  return rcode;
} 
unsigned
F_prefix ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return F_prefix_generic(output,node,variable,symbol[node->token].name);
} 
unsigned
F_lt ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return (F_FUNCT_LT | F_FUNCT_GE | F_prefix_generic(output,node,variable,"RelLT"));
} 
unsigned
F_le ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return (F_FUNCT_LE | F_prefix_generic(output,node,variable,"RelLE"));
} 
unsigned
F_gt ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return (F_FUNCT_GT | F_FUNCT_LE | F_prefix_generic(output,node,variable,"RelGT"));
} 
unsigned
F_ge ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return (F_FUNCT_GE | F_prefix_generic(output,node,variable,"RelGE"));
} 
unsigned
F_eq ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return (F_FUNCT_EQ | F_FUNCT_LE | F_FUNCT_GE | F_prefix_generic(output,node,variable,"RelEQ"));
} 
unsigned
F_ne ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return (F_FUNCT_NE | F_FUNCT_LE | F_FUNCT_GE | F_prefix_generic(output,node,variable,"RelNE"));
} 
unsigned
F_sec ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return (F_FUNCT_SEC | F_prefix_generic(output,node,variable,"sec"));
} 
unsigned
F_sgn ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return (F_FUNCT_SGN | F_prefix_generic(output,node,variable,"sgn"));
} 
unsigned
F_infix_generic ( struct unur_string *output, const struct ftreenode *node,
		  const char *variable, const char *symb )
{
  unsigned rcode = 0u;
  struct ftreenode *left  = node->left;    
  struct ftreenode *right = node->right;   
  if (left==NULL || right==NULL)
    return F_FUNCT_ERROR;
  _unur_fstr_print_F( output, "(", 0. );
  rcode |= symbol[left->token].node2F (output,left,variable);
  _unur_fstr_print_F( output, symb, 0. );
  rcode |= symbol[right->token].node2F (output,right,variable);
  _unur_fstr_print_F( output, ")", 0. );
  return rcode;
} 
unsigned
F_infix ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return F_infix_generic(output,node,variable,symbol[node->token].name);
} 
unsigned
F_minus ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  unsigned rcode = 0u;
  struct ftreenode *left  = node->left;    
  struct ftreenode *right = node->right;   
  if (left==NULL || right==NULL)
    return F_FUNCT_ERROR;
  _unur_fstr_print_F( output, "(", 0. );
  if (!(left->type == S_UCONST && _unur_iszero(left->val)))
    rcode |= symbol[left->token].node2F (output,left,variable);
  _unur_fstr_print_F( output, "-", 0. );
  rcode |= symbol[right->token].node2F (output,right,variable);
  _unur_fstr_print_F( output, ")", 0. );
  return rcode;
} 
unsigned
F_power ( struct unur_string *output, const struct ftreenode *node, const char *variable )
{
  return F_infix_generic(output,node,variable,"**");
} 
int
_unur_fstr_F_specfunct ( FILE *out, unsigned flags )
{
  if (flags & F_FUNCT_SEC) {
    fprintf (out,"      sec(a)=1.d0/cos(a)\n");
  }
  if (flags & F_FUNCT_SGN) {
    fprintf (out,"      sgn(a)=sign(1.d0,a)\n");
  }
  if (flags & F_FUNCT_LE) {
    fprintf (out,"      RelLE(a,b)=sign(0.5d0,b-a)+0.5d0\n");
  }
  if (flags & F_FUNCT_GE) {
    fprintf (out,"      RelGE(a,b)=sign(0.5d0,a-b)+0.5d0\n");
  }
  if (flags & F_FUNCT_LT) {
    fprintf (out,"      RelLT(a,b)=1.d0-RelGE(a,b)\n");
  }
  if (flags & F_FUNCT_GT) {
    fprintf (out,"      RelGT(a,b)=1.d0-RelLE(a,b)\n");
  }
  if (flags & F_FUNCT_EQ) {
    fprintf (out,"      RelEQ(a,b)=RelGE(a,b)*RelLE(a,b)\n");
  }
  if (flags & F_FUNCT_NE) {
    fprintf (out,"      RelNE(a,b)=1.d0-RelGE(a,b)*RelLE(a,b)\n");
  }
  return UNUR_SUCCESS;
} 
int
_unur_fstr_print_F ( struct unur_string *output, const char *symb, double number )
{
  char buf[128];
  char *here_is_e;
  if (symb) {
    _unur_string_appendtext( output, symb );
  }
  else {
    sprintf(buf,"%.20e",number);
    here_is_e = strchr(buf, 'e');
    if (here_is_e) *here_is_e = 'D';
    _unur_string_appendtext( output, buf );
  }
  return UNUR_SUCCESS;
} 
