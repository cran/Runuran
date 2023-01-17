/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

char *
_unur_fstr_tree2string ( const struct ftreenode *root,
			 const char *variable, const char *function, int spaces )
{
  struct unur_string output = {NULL, 0, 0};
  _unur_check_NULL( GENTYPE,root,NULL );
  _unur_fstr_node2string(&output,root,variable,function,spaces);
  return output.text;
} 
int
_unur_fstr_node2string ( struct unur_string *output, const struct ftreenode *node,
			 const char *variable, const char *function, int spaces )
{
  struct ftreenode *left  = node->left;    
  struct ftreenode *right = node->right;   
  const char *symb;                           
  int type = node->type;                   
  int priority = symbol[node->token].info; 
  int operator, parenthesis;               
  switch (type) {
  case S_UIDENT:    
    symb = variable;  break;
  case S_UFUNCT:    
    symb = function;  break;
  case S_UCONST:    
    symb = NULL;      break;  
  case S_SCONST:
  default:
    symb = node->symbol;
  }
  if (type == S_SFUNCT || type == S_UFUNCT) {
    _unur_fstr_print( output, symb, 0. );
    _unur_fstr_print( output, "(", 0. );
    if (left) {
      _unur_fstr_node2string(output,left,variable,function,spaces);
      _unur_fstr_print( output, ",", 0. );
    }
    if (right) {
      _unur_fstr_node2string(output,right,variable,function,spaces);
    }
    _unur_fstr_print( output, ")", 0. );
  }
  else if (symb && symb[0] == ',') {
    _unur_fstr_print( output, ",", 0. );
    if (left) {
      _unur_fstr_node2string(output,left,variable,function,spaces);
      _unur_fstr_print( output, ",", 0. );
    }
    if (right) {
      _unur_fstr_node2string(output,right,variable,function,spaces);
    }
  }    
  else {
    operator = (type==S_REL_OP || type==S_ADD_OP || type==S_MUL_OP);
    if (left) {
      parenthesis = 1;
      if (left->type == S_SCONST || left->type == S_UCONST || 
	  left->type == S_SFUNCT || left->type == S_UFUNCT || 
	  ( left->type == S_UIDENT && left->val >= 0. ) ||
	  ( priority < symbol[left->token].info && !isalpha(node->symbol[0]) ) ||
	  ( priority == symbol[left->token].info && (type == S_ADD_OP ) ) )
	parenthesis = 0;
      if (parenthesis) _unur_fstr_print( output, "(", 0. );
      if (left->type == S_UCONST && _unur_iszero(left->val) && node->symbol[0] == '-')
	 ;
      else
	_unur_fstr_node2string(output,left,variable,function,spaces);
      if (parenthesis) _unur_fstr_print( output, ")", 0. );
    }
    if (operator && spaces) _unur_fstr_print( output, " ", 0. );
    _unur_fstr_print( output, symb, node->val );
    if (operator && spaces) _unur_fstr_print( output, " ", 0. );
    if (right) {
      parenthesis = 1;
      if (right->type == S_SCONST || right->type == S_UCONST ||
	  right->type == S_SFUNCT || right->type == S_UFUNCT || 
	  ( right->type == S_UIDENT && right->val >= 0. ) ||
	  ( priority < symbol[right->token].info && !isalpha(node->symbol[0]) ) )
	parenthesis = 0;
      if (parenthesis) _unur_fstr_print( output, "(", 0. );
      _unur_fstr_node2string(output,right,variable,function,spaces);
      if (parenthesis) _unur_fstr_print( output, ")", 0. );
    }
  }
  return UNUR_SUCCESS;
} 
int
_unur_fstr_print ( struct unur_string *output, const char *symb, double number )
{
  if (symb)
    _unur_string_appendtext( output, symb );
  else
    _unur_string_append( output, "%.16g", number);
  return UNUR_SUCCESS;
} 
