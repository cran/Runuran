/* Copyright (c) 2000-2014 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct ftreenode *
_unur_FunctDefinition (struct parser_data *pdata)
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb;
  int              token; 
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);
  left = _unur_DefFunctDesignator(pdata);
  if (pdata->perrno) {
    _unur_fstr_free(left); return NULL;
  }
  if ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS ||
       strcmp(symb,"=") != 0 ) {
    _unur_fstr_free(left);
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_EQUAL,__LINE__); 
  }
  right = _unur_Expression(pdata);
  if (pdata->perrno) {
    _unur_fstr_free(left); _unur_fstr_free(right);
    return NULL;
  }
  node = _unur_fstr_create_node(symb,0.,token,left,right); 
  return node; 
} 
struct ftreenode *
_unur_DefFunctDesignator (struct parser_data *pdata) 
{ 
  struct ftreenode *node, *params; 
  char             *fsymb, *symb;
  int              n_params; 
  int              funct, token;
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);
  if ( _unur_fstr_next_token(pdata,&funct,&fsymb) != UNUR_SUCCESS ||
       symbol[funct].type != S_UFUNCT )
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_FUNCT,__LINE__); 
  if ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS ||
       symb[0] != '(' )
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_OPEN_P,__LINE__); 
  params = _unur_DefParameterlist(pdata,&n_params);
  if (pdata->perrno) {
    _unur_fstr_free(params);
    return NULL;
  }
  if ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS ||
       symb[0] != ')' ) {
    _unur_fstr_free(params);
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_CLOSE_P,__LINE__); 
  }
  node = _unur_fstr_create_node(fsymb,0.,funct,NULL,params); 
  return node;
} 
struct ftreenode *
_unur_DefParameterlist(struct parser_data *pdata, int *n_params) 
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb;
  int              token;
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);
  if ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS ||
       symbol[token].type != S_UIDENT )
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_VAR,__LINE__);
  node = _unur_fstr_create_node(symb,0.,token,NULL,NULL); 
  *n_params = 1;
  while ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS ||
	  symb[0] == ',' ) {
    if ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS ||
	 symbol[token].type != S_UIDENT ) {
      _unur_fstr_free(node);
      return _unur_fstr_error_parse(pdata,ERR_EXPECT_VAR,__LINE__);
    }
    left = node; 
    right = _unur_fstr_create_node(symb,0.,token,NULL,NULL); 
    (*n_params)++; 
    node = _unur_fstr_create_node(",",0.,s_comma,left,right); 
  }
  --(pdata->tno);
  return node; 
}  
struct ftreenode *
_unur_Expression (struct parser_data *pdata)
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb;
  int              token;
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);
  left = _unur_SimpleExpression(pdata);
  if (pdata->perrno) {
    _unur_fstr_free(left); return NULL; 
  }
  if ( _unur_fstr_next_token(pdata,&token,&symb) == UNUR_SUCCESS &&
       symbol[token].type == S_REL_OP ) {
    right = _unur_SimpleExpression(pdata);
    if (pdata->perrno) {
      _unur_fstr_free(left); _unur_fstr_free(right);
      return NULL;
    } 
    node = _unur_fstr_create_node(symb,0.,token,left,right); 
  }
  else {
    --(pdata->tno);
    node = left; 
  } 
  return node; 
} 
struct ftreenode *
_unur_SimpleExpression (struct parser_data *pdata)
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb;
  int              token;
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);
  node = _unur_STerm(pdata);
  if (pdata->perrno) {
    _unur_fstr_free(node); return NULL;
  }
  while ( _unur_fstr_next_token(pdata,&token,&symb) == UNUR_SUCCESS &&
	  symbol[token].type == S_ADD_OP) {
    left = node; 
    right = _unur_Term(pdata);
    if (pdata->perrno) {
      _unur_fstr_free(left); _unur_fstr_free(right);
      return NULL;
    }
    node = _unur_fstr_create_node(symb,0.,token,left,right); 
  }
  --(pdata->tno);
  return node; 
} 
struct ftreenode *
_unur_STerm (struct parser_data *pdata)
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb; 
  int              token;
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);
  if ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS)
    return NULL;
  if ( symb[0] == '-' ) {
    left = _unur_fstr_create_node(NULL,0.,s_uconst,NULL,NULL); 
    right = _unur_Term(pdata);
    if (pdata->perrno) {
      _unur_fstr_free(left); _unur_fstr_free(right);
      return NULL; 
    }
    node = _unur_fstr_create_node(symb,0.,token,left,right); 
  }
  else {
    if( symb[0] != '+' ) {
      --(pdata->tno);
    }
    node = _unur_Term(pdata);
    if (pdata->perrno) {
      _unur_fstr_free(node); return NULL; 
    }
  }
  return node; 
}  
struct ftreenode *
_unur_Term (struct parser_data *pdata)
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb;
  int              token;
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);
  node = _unur_Factor(pdata);
  if (pdata->perrno) {
    _unur_fstr_free(node); return NULL;
  }
  while ( _unur_fstr_next_token(pdata,&token,&symb) == UNUR_SUCCESS &&
	  symbol[token].type == S_MUL_OP ) {
     left = node; 
     right = _unur_Factor(pdata);
     if (pdata->perrno) {
       _unur_fstr_free(left); _unur_fstr_free(right);
       return NULL;
     }
     node = _unur_fstr_create_node(symb,0.,token,left,right); 
  }
  --(pdata->tno);
  return node; 
} 
struct ftreenode *
_unur_Factor (struct parser_data *pdata)
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb;
  int              token;
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);
  left = _unur_Bas_Exp(pdata);
  if (pdata->perrno) {
    _unur_fstr_free(left); return NULL;
  }
  if ( _unur_fstr_next_token(pdata,&token,&symb) == UNUR_SUCCESS &&
       symb[0] == '^' ) {
    right = _unur_Bas_Exp(pdata);
    if (pdata->perrno) {
      _unur_fstr_free(left); _unur_fstr_free(right);
      return NULL;
    }
    node = _unur_fstr_create_node(symb,0.,token,left,right); 
  }
  else {
    --(pdata->tno);
    node = left; 
  } 
  return node; 
} 
struct ftreenode *
_unur_Bas_Exp (struct parser_data *pdata)
{ 
  struct ftreenode *node; 
  char             *symb;
  int              token;
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);
  if ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS)
    return _unur_fstr_error_parse(pdata,7,__LINE__); 
  if( symbol[token].type==S_UCONST ||
      symbol[token].type==S_UIDENT ||
      symbol[token].type==S_SCONST ) {
    node = _unur_fstr_create_node(symb,0.,token,NULL,NULL); 
  }
  else if( symbol[token].type == S_SFUNCT ) {
    --(pdata->tno);
    node = _unur_FunctDesignator(pdata);
    if (pdata->perrno) {
      _unur_fstr_free(node); return NULL;
    }
  }
  else if( symb[0] == '(' ) {
    node = _unur_Expression(pdata); 
    if (pdata->perrno) {
      _unur_fstr_free(node); return NULL;
    }
    if ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS ||
	 symb[0] != ')' )
      return _unur_fstr_error_parse(pdata,ERR_EXPECT_CLOSE_P,__LINE__);
  }
  else {
    --(pdata->tno);
    return _unur_fstr_error_parse(pdata,ERR_UNKNOWN_SYMBOL,__LINE__);
  } 
  return node; 
}  
struct ftreenode *
_unur_FunctDesignator (struct parser_data *pdata)
{ 
  struct ftreenode *node, *params; 
  char             *fsymb, *symb;
  int              funct, token;
  int              n_params; 
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);
  if ( _unur_fstr_next_token(pdata,&funct,&fsymb) != UNUR_SUCCESS ||
       symbol[funct].type != S_SFUNCT )
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_FUNCT,__LINE__);
  n_params = symbol[funct].info;
  if ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS ||
       symb[0] != '(' )
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_OPEN_P,__LINE__);
  params = _unur_ActualParameterlist(pdata,n_params);
  if (pdata->perrno) {
    _unur_fstr_free(params); return NULL;
  }
  if ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS ||
       symb[0] != ')' ) {
    _unur_fstr_free(params);
    return _unur_fstr_error_parse(pdata,ERR_EXPECT_CLOSE_P,__LINE__);
  }
  node = _unur_fstr_create_node(fsymb,0.,funct,NULL,params); 
  return node; 
} 
struct ftreenode *
_unur_ActualParameterlist (struct parser_data *pdata, int n_params) 
{ 
  struct ftreenode *node, *left, *right; 
  char             *symb;
  int              token;
  int              c_params;   
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);
  node = _unur_Expression(pdata);
  if (pdata->perrno) {
    _unur_fstr_free(node); return NULL;
  }
  c_params = 1; 
  while ( _unur_fstr_next_token(pdata,&token,&symb) != UNUR_SUCCESS ||
	  symb[0] == ',' ) {
    c_params++; 
    if (c_params > n_params) {
      _unur_fstr_free(node);
      return _unur_fstr_error_parse(pdata,ERR_INVALID_N_PARAMS,__LINE__);
    }
    left = node; 
    right = _unur_Expression(pdata);
    if (pdata->perrno) {
      _unur_fstr_free(left); _unur_fstr_free(right);
      return NULL;
    }
    node = _unur_fstr_create_node(",",0.,s_comma,left,right); 
  }
  --(pdata->tno);
  if (c_params < n_params) {
    _unur_fstr_free(node);
    return _unur_fstr_error_parse(pdata,ERR_INVALID_N_PARAMS,__LINE__);
  }
  return node; 
} 
struct ftreenode *
_unur_fstr_simplification (const char *symb, int token,
			   struct ftreenode *left, struct ftreenode *right) 
{
  int l_const = left  && (left->type  == S_SCONST || left->type  == S_UCONST); 
  int r_const = right && (right->type == S_SCONST || right->type == S_UCONST); 
  int l_0     = (l_const && _unur_iszero(left->val));
  int l_1     = (l_const && _unur_isone(left->val));
  int r_0     = (r_const && _unur_iszero(right->val));
  int r_1     = (r_const && _unur_isone(right->val));
  int and;
  char s = symb[0];
  if ( left == NULL && right && right->symbol[0] == ',' ) {
    COOKIE_CHECK(right,CK_FSTR_TNODE,NULL);
    right->token  = token;
    right->symbol = symbol[token].name; 
    right->type   = symbol[token].type;
    return right; 
   }
  if ( (l_const || left==NULL) && r_const && s!=',') { 
    CHECK_NULL(right,NULL); COOKIE_CHECK(right,CK_FSTR_TNODE,NULL);
    right->val   = ( (left) 
		     ? (*symbol[token].vcalc)(left->val,right->val)
		     : (*symbol[token].vcalc)(0.,right->val) );
    right->token = s_uconst;
    right->type  = S_UCONST;
    right->left  = NULL; 
    right->right = NULL;
    _unur_fstr_free(left);
    return right; 
  } 
  if ( (l_0 && s=='+' ) || (l_1 && s=='*') ) { 
    CHECK_NULL(left,NULL);  COOKIE_CHECK(left,CK_FSTR_TNODE,NULL);
    CHECK_NULL(right,NULL); COOKIE_CHECK(right,CK_FSTR_TNODE,NULL);
    _unur_fstr_free(left);
    return right;
  } 
  if ( (r_0 && (s=='+' || s=='-')) ||
       (r_1 && (s=='*' || s=='/' || s=='^')) ) {
    CHECK_NULL(left,NULL);  COOKIE_CHECK(left,CK_FSTR_TNODE,NULL);
    CHECK_NULL(right,NULL); COOKIE_CHECK(right,CK_FSTR_TNODE,NULL);
    _unur_fstr_free(right);
    return left;
  }
  and = (strcmp(symb,"and")==0);
  if ( l_0 && (s=='*' || s=='/' || s=='^' || and) ) {
    CHECK_NULL(left,NULL);  COOKIE_CHECK(left,CK_FSTR_TNODE,NULL);
    CHECK_NULL(right,NULL); COOKIE_CHECK(right,CK_FSTR_TNODE,NULL);
    _unur_fstr_free(right);
    return left;
  }
  if (r_0 && (s=='*' || and ) ) {
    CHECK_NULL(left,NULL);  COOKIE_CHECK(left,CK_FSTR_TNODE,NULL);
    CHECK_NULL(right,NULL); COOKIE_CHECK(right,CK_FSTR_TNODE,NULL);
    _unur_fstr_free(left);
    return right;
  }
  if (r_0 && s=='^') {
    CHECK_NULL(left,NULL);  COOKIE_CHECK(left,CK_FSTR_TNODE,NULL);
    CHECK_NULL(right,NULL); COOKIE_CHECK(right,CK_FSTR_TNODE,NULL);
    _unur_fstr_free(left);
    right->val = 1.;
    return right;
  }
  if (l_1 && s=='^') {
    CHECK_NULL(left,NULL);  COOKIE_CHECK(left,CK_FSTR_TNODE,NULL);
    CHECK_NULL(right,NULL); COOKIE_CHECK(right,CK_FSTR_TNODE,NULL);
    _unur_fstr_free(right);
    return left;
  }
  if ( ( symb[0] == '/' &&
	 left  && left->left==NULL  && left->right==NULL  && 
	 right && right->left==NULL && right->right==NULL &&
	 strcmp(left->symbol,right->symbol)== 0 ) ) {
    CHECK_NULL(left,NULL);  COOKIE_CHECK(left,CK_FSTR_TNODE,NULL);
    CHECK_NULL(right,NULL); COOKIE_CHECK(right,CK_FSTR_TNODE,NULL);
    _unur_fstr_free(left);
    right->token = s_uconst;
    right->symbol= symbol[s_uconst].name; 
    right->val   = 1.;
    right->type  = S_UCONST;
    right->left  = NULL; 
    right->right = NULL;
    return right; 
  }
  return NULL; 
} 
int
_unur_fstr_reorganize (struct ftreenode *node) 
{
  struct ftreenode *left, *right, *tmp;
  char symb;
  int l_const, r_const;
  int rl_0, ll_0;
  CHECK_NULL(node,0);  COOKIE_CHECK(node,CK_FSTR_TNODE,0);
  left  = node->left;
  right = node->right;
  symb = node->symbol[0];
  l_const = left  && (left->type  == S_SCONST || left->type  == S_UCONST); 
  r_const = right && (right->type == S_SCONST || right->type == S_UCONST); 
  rl_0 = (right && right->left && right->left->type == S_UCONST && _unur_iszero(right->left->val));
  ll_0 = (left  && left->left  && left->left->type  == S_UCONST && _unur_iszero(left->left->val));
  if ( (l_const || left==NULL) && r_const && symb!=',') { 
    node->val   = ( (left) 
		    ? (*symbol[node->token].vcalc)(left->val,right->val)
		    : (*symbol[node->token].vcalc)(0.,right->val) );
    node->token = s_uconst;
    node->type  = S_UCONST;
    node->left  = NULL; 
    node->right = NULL;
    if (left)  free(left);
    if (right) free(right);
    return 1; 
  } 
  if ( rl_0 && symb=='+' && right->symbol[0]=='-' ) {
    node->symbol = symbol[s_minus].name; 
    node->token  = s_minus; 
    node->type   = symbol[s_minus].type; 
    node->right  = right->right; 
    free(right->left);
    free(right);
    return 1;
  }
  if ( rl_0 && symb=='-' && right->symbol[0]=='-' ) {
    node->symbol = symbol[s_plus].name; 
    node->token  = s_plus; 
    node->type   = symbol[s_plus].type; 
    node->right  = right->right; 
    free(right->left);
    free(right);
    return 1;
  }
  if ( ll_0 && symb=='+' && left->symbol[0]=='-' ) {
    node->symbol = symbol[s_minus].name; 
    node->token  = s_minus; 
    node->type   = symbol[s_minus].type; 
    tmp = node->right;
    node->right  = left->right; 
    node->left   = tmp;
    free(left->left);
    free(left);
    return 1;
  }
  if ( rl_0 && symb=='*' && right->symbol[0]=='-' ) {
    node->symbol = symbol[s_minus].name; 
    node->token  = s_minus; 
    node->type   = symbol[s_minus].type; 
    right->symbol= symbol[s_mul].name; 
    right->token = s_mul; 
    right->type  = symbol[s_mul].type; 
    tmp = left;
    node->left   = node->right->left;
    node->right  = tmp;
    return 1;
  }
  return 0;
} 
int
_unur_fstr_next_token (struct parser_data *pdata, int *token, char **symb)
{
  CHECK_NULL(pdata,UNUR_ERR_NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,UNUR_ERR_COOKIE);
  if (pdata->tno < pdata->n_tokens) {
    *token = pdata->token[pdata->tno];
    *symb = pdata->tpos[pdata->tno];
    ++(pdata->tno);
    return UNUR_SUCCESS;
  }
  else {
    ++(pdata->tno);
    return UNUR_ERR_SILENT;
  }
} 
struct ftreenode *
_unur_fstr_create_node (const char *symb, double val, int token, 
			struct ftreenode *left, struct ftreenode *right) 
{ 
  struct ftreenode *node; 
  if ( symb && (node = _unur_fstr_simplification(symb,token,left,right)) ) {
  }
  else {
    node = _unur_xmalloc(sizeof(struct ftreenode)); 
    COOKIE_SET(node,CK_FSTR_TNODE);
    node->symbol = symbol[token].name; 
    node->token  = token; 
    node->type   = symbol[token].type; 
    node->left   = left; 
    node->right  = right; 
    switch (symbol[token].type) {
    case S_UCONST:      
      node->val = (symb) ? atof(symb) : val;  break;
    case S_SCONST:      
      node->val = symbol[token].val;  break;
    default:
      node->val = 0.;
    }
  } 
  _unur_fstr_reorganize(node); 
  return node; 
} 
struct ftreenode *
_unur_fstr_error_parse ( struct parser_data *pdata, int perrno, int line )
{ 
  int i;
  struct unur_string *reason;
  CHECK_NULL(pdata,NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,NULL);
  if (!pdata->perrno) pdata->perrno = perrno;
  reason = _unur_string_new();
  _unur_string_append( reason, "%s: ", _unur_fstr_error_code(perrno) );
  for (i=0; i<pdata->tno-1; i++)
    _unur_string_append( reason, "%s ", pdata->tpos[i]);
  if (i<pdata->n_tokens)
    _unur_string_append( reason, " -->%s<--  ", pdata->tpos[i]);
  else
    _unur_string_append( reason, " <--  ");
  for (i++; i<pdata->n_tokens; i++)
    _unur_string_append( reason, "%s ",pdata->tpos[i]);
  _unur_error_x( GENTYPE, __FILE__, line, "error", UNUR_ERR_FSTR_SYNTAX,reason->text);
  _unur_string_free( reason );
  return NULL; 
}  
const char *
_unur_fstr_error_code ( int perrno )
{
  switch (perrno) {
  case ERR_UNFINISHED:
    return "incomplete. not all tokens parsed";
  case ERR_UNKNOWN_SYMBOL:
    return "unknown symbol in function string";
  case ERR_EXPECT_EQUAL:
    return "expected symbol: '='";
  case ERR_EXPECT_OPEN_P:
    return "expected symbol: '('";
  case ERR_EXPECT_CLOSE_P:
    return "expected symbol: ')'";
  case ERR_INVALID_N_PARAMS:
    return "invalid number of parameters for function";
  case ERR_EXPECT_FUNCT:
    return "function (name) expected";
  case ERR_EXPECT_VAR:
    return "user identifier (variable name) expected";
  case ERR_MISSING:
    return "more tokens expected";
  default:
    return "";
  }
} 
