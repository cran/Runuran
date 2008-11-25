/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct parser_data *
_unur_fstr_parser_init ( const char *fstr )
{
  static int symbols_initialized = FALSE;
  struct parser_data *pdata;
  if (symbols_initialized == FALSE) {
    _unur_fstr_symbols_init();
    symbols_initialized = TRUE;
  }
  pdata = _unur_xmalloc(sizeof(struct parser_data));
  COOKIE_SET(pdata,CK_FSTR_PDATA);
  pdata->fstr = _unur_parser_prepare_string(fstr);
  pdata->len_fstr = strlen(pdata->fstr);
  if (pdata->len_fstr <= 0) {
    _unur_error(GENTYPE,UNUR_ERR_STR,"empty string"); 
    if (pdata->fstr) free (pdata->fstr);
    free(pdata);
    return NULL;
  }
  pdata->token = _unur_xmalloc( (pdata->len_fstr+1) * sizeof(int) );
  pdata->tpos  = _unur_xmalloc( (pdata->len_fstr+1) * sizeof(char *) );
  pdata->tstr  = _unur_xmalloc( (2*pdata->len_fstr+1) * sizeof(char) );
  pdata->n_tokens = 0;
  memset(pdata->token,0  ,(size_t)pdata->len_fstr);
  memset(pdata->tpos, 0  ,(size_t)pdata->len_fstr);
  memset(pdata->tstr,'\0',(size_t)pdata->len_fstr);
  pdata->scanpos = 0;     
  pdata->lastpos = -1;
  pdata->perrno = 0;
  pdata->variable_name = NULL;
  pdata->function_name = NULL;
  pdata->tno = 0;
  return pdata;
} 
void
_unur_fstr_symbols_init (void)
{
  int i;
  char *s;
  _ros_start = 0;
  _nas_start = 0;
  _ans_start = 0;
  for (i=0; !_end; i++) {
    s = symbol[i].name;
    if (!_ros_start) {
      if ( strcmp(s,"_ROS") == 0) _ros_start = i;
      continue;
    }
    if (!_nas_start) {
      if ( strcmp(s,"_NAS") == 0) _nas_start = i;
      continue;
    }
    if (!_ans_start) {
      if ( strcmp(s,"_ANS") == 0) _ans_start = i;
      continue;
    }
    if (strcmp(s,"_END") == 0) _end = i;
  }
  _ros_end = _nas_start;
  _nas_end = _ans_start;
  _ans_end = _end;
  s_comma = _unur_fstr_find_symbol(",",_nas_start,_nas_end);
  s_minus = _unur_fstr_find_symbol("-",_nas_start,_nas_end);
  s_plus  = _unur_fstr_find_symbol("+",_nas_start,_nas_end);
  s_mul   = _unur_fstr_find_symbol("*",_nas_start,_nas_end);
  s_div   = _unur_fstr_find_symbol("/",_nas_start,_nas_end);
  s_power = _unur_fstr_find_symbol("^",_nas_start,_nas_end);
} 
void
_unur_fstr_parser_free ( struct parser_data *pdata )
{
  if (pdata) {
    COOKIE_CHECK(pdata,CK_FSTR_PDATA,RETURN_VOID);
    free(pdata->fstr);
    free(pdata->token);
    free(pdata->tpos);
    free(pdata->tstr);
    free(pdata);
  }
} 
int
_unur_fstr_tokenize (struct parser_data *pdata)
{
  int token;
  int n_token = 0;                    
  char *symb = pdata->tstr;           
  CHECK_NULL(pdata,UNUR_ERR_COOKIE);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,UNUR_ERR_COOKIE);
  while ((token = _unur_fstr_next_symbol(pdata,symb)) != S_NOSYMBOL) {
    pdata->token[n_token] = token;    
    pdata->tpos[n_token] = symb;  
    n_token++;                        
    symb += pdata->scanpos - pdata->lastpos + 1;   
  }
  pdata->n_tokens = n_token;
  pdata->tno = 0;
  return pdata->perrno;
} 
int
_unur_fstr_next_symbol (struct parser_data *pdata, char *symb)
{
  int token;
  int errcode = 0;
  char c;
  CHECK_NULL(pdata,S_NOSYMBOL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,S_NOSYMBOL);
  pdata->lastpos = pdata->scanpos;
  if (pdata->scanpos >= pdata->len_fstr)
    return S_NOSYMBOL;
  c = pdata->fstr[pdata->scanpos];
  if ( (c >= '0' && c <= '9') || c == '.') {
    _unur_fstr_UnsignedConstant(pdata,symb);
    token = s_uconst;
  }
  else if (c >=  'a' && c <= 'z') {
    _unur_fstr_Identifier(pdata,symb);
    if ( ( (token = _unur_fstr_find_symbol(symb,_ans_start,_ans_end)) == 0 ) && 
	 ( (token = _unur_fstr_find_user_defined(pdata,symb,pdata->fstr[pdata->scanpos])) <= 0 ) )
      errcode = ERR_UNKNOWN_SYMBOL;
  }
  else if ( c == '<' || c == '>' || c == '=' || c == '!' ) {
    _unur_fstr_RelationOperator(pdata,symb);
    if ((token = _unur_fstr_find_symbol(symb,_ros_start,_ros_end)) <= 0 )
      errcode = ERR_UNKNOWN_SYMBOL;
  }
  else {
    symb[0] = c; symb[1] = '\0';           
    (pdata->scanpos)++;
    if ((token = _unur_fstr_find_symbol(symb,_nas_start,_nas_end)) <= 0 )
      errcode = ERR_UNKNOWN_SYMBOL;
  }
  pdata->perrno = errcode;
  if (errcode) {
    _unur_fstr_error_scan (pdata,symb,__LINE__);
  }
  return (errcode) ? S_NOSYMBOL : token;
} 
int
_unur_fstr_find_symbol (const char *symb, int start, int end)
{
  int i;
  for (i = start + 1; i < end; i++)
    if (strcmp(symb,symbol[i].name) == 0) break;
  return ((i < end ) ? i : 0);
} 
int
_unur_fstr_find_user_defined (struct parser_data *pdata, char *symb, int next_char)
{
  CHECK_NULL(pdata,0);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,0);
  if (next_char == '(') {
    if (pdata->function_name == NULL) {
      pdata->function_name = symb;
      return s_ufunct;
    }
    else
      return (strcmp(pdata->function_name,symb) == 0) ? s_ufunct : 0;
  }
  else {
    if (pdata->variable_name == NULL) {
      pdata->variable_name = symb;
      return s_uident;
    }
    else
      return (strcmp(pdata->variable_name,symb) == 0) ? s_uident : 0;
  }
} 
int 
_unur_fstr_UnsignedConstant (struct parser_data *pdata, char *uc)
{
  int startpos = pdata->scanpos;
  CHECK_NULL(pdata,UNUR_ERR_NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,UNUR_ERR_COOKIE);
  _unur_fstr_DigitalSequence(pdata,uc);
  if( pdata->fstr[pdata->scanpos] == '.' ) {
    *(uc + pdata->scanpos - startpos) = '.';
    (pdata->scanpos)++;
    _unur_fstr_DigitalSequence(pdata, uc + pdata->scanpos - startpos);
  }
  if( pdata->fstr[pdata->scanpos] == 'e' ) {
    *(uc + pdata->scanpos - startpos) = 'e';
    (pdata->scanpos)++;
    _unur_fstr_ScaleFactor(pdata, uc + pdata->scanpos - startpos);
  }
  return UNUR_SUCCESS;
} 
int
_unur_fstr_DigitalSequence (struct parser_data *pdata, char *ds)
{
  CHECK_NULL(pdata,UNUR_ERR_NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,UNUR_ERR_COOKIE);
  while ( (*ds = pdata->fstr[pdata->scanpos]) >= '0' && *ds <= '9' ) {
     ds++;
     (pdata->scanpos)++;
  }
  *ds = '\0';
  return UNUR_SUCCESS;
} 
int 
_unur_fstr_ScaleFactor (struct parser_data *pdata, char *sf)
{
  CHECK_NULL(pdata,UNUR_ERR_NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,UNUR_ERR_COOKIE);
  if ( (sf[0] = pdata->fstr[pdata->scanpos]) == '+' || sf[0] == '-' ) {
     sf++;
     (pdata->scanpos)++;
  }
  _unur_fstr_DigitalSequence(pdata,sf);
  return UNUR_SUCCESS;
} 
int
_unur_fstr_Identifier (struct parser_data *pdata, char *id)
{
  CHECK_NULL(pdata,UNUR_ERR_NULL);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,UNUR_ERR_COOKIE);
  while ( ((*id = pdata->fstr[pdata->scanpos]) >= 'a' && *id <= 'z')
	  || *id == '_' 
	  || ( *id >= '0' && *id <= '9')) {
    id++;
    (pdata->scanpos)++;
  }
  *id = '\0';
  return UNUR_SUCCESS;
} 
int
_unur_fstr_RelationOperator (struct parser_data *pdata, char *ro)
{
  CHECK_NULL(pdata,1);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,1);
  while ((*ro = pdata->fstr[pdata->scanpos]) == '<' || *ro == '>' || *ro == '=' || *ro == '!' ) {
    ro++;
    (pdata->scanpos)++;
  }
  *ro = '\0';
  return UNUR_SUCCESS;
} 
void
_unur_fstr_error_scan (const struct parser_data *pdata, const char *symb, int line)
{
  struct unur_string *reason;
  char *c;
  CHECK_NULL(pdata,RETURN_VOID);  COOKIE_CHECK(pdata,CK_FSTR_PDATA,RETURN_VOID);
  reason = _unur_string_new();
  _unur_string_append( reason, "unknown symbol '%s': ", symb );
  for (c=pdata->fstr; c < pdata->fstr+pdata->lastpos; c++) 
    _unur_string_append( reason, "%c", *c );
  _unur_string_append( reason, "    %s", pdata->fstr + pdata->lastpos);
  _unur_error_x( GENTYPE, __FILE__, line, "error", UNUR_ERR_FSTR_SYNTAX,reason->text);
  _unur_string_free( reason );
} 
