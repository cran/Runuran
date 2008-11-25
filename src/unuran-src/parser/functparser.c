/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#include <ctype.h>
#include <unur_source.h>
#include "parser.h"
#include "parser_source.h"
#include "functparser_source.h"
#define GENTYPE "FSTRING"      
struct parser_data {
  char  *fstr;          
  int   *token;         
  char  *tstr;          
  char **tpos;          
  int    tno;           
  int    n_tokens;      
  char  *variable_name; 
  char  *function_name; 
  int    scanpos;       
  int    lastpos;       
  int    len_fstr;      
  int    perrno;        
#ifdef UNUR_COOKIES
  unsigned cookie;      
#endif
};
enum {
  SUCCESS = 0,          
  ERR_UNFINISHED,       
  ERR_UNKNOWN_SYMBOL,   
  ERR_EXPECT_EQUAL,     
  ERR_EXPECT_OPEN_P,    
  ERR_EXPECT_CLOSE_P,   
  ERR_INVALID_N_PARAMS, 
  ERR_EXPECT_FUNCT,      
  ERR_EXPECT_VAR,       
  ERR_MISSING            
};
#define PARSER
#include "functparser_symbols.h"
#undef PARSER
static struct ftreenode *_unur_fstr_2_tree (const char *functstr, int withDefFunct);
static struct parser_data *_unur_fstr_parser_init (const char *fstr);
static void _unur_fstr_symbols_init (void);
static void _unur_fstr_parser_free (struct parser_data *pdata);
static int _unur_fstr_tokenize (struct parser_data *pdata);
static int _unur_fstr_next_symbol (struct parser_data *pdata, char *symb);
static int _unur_fstr_find_symbol (const char *symb, int start, int end);
static int _unur_fstr_find_user_defined (struct parser_data *pdata, char *symb, int next_char);
static int _unur_fstr_UnsignedConstant (struct parser_data *pdata, char *uc);
static int _unur_fstr_DigitalSequence (struct parser_data *pdata, char *ds);
static int _unur_fstr_ScaleFactor (struct parser_data *pdata, char *sf);
static int _unur_fstr_Identifier (struct parser_data *pdata, char *id);
static int _unur_fstr_RelationOperator (struct parser_data *pdata, char *ro);
static void _unur_fstr_error_scan (const struct parser_data *pdata, const char *symb, int line);
static struct ftreenode *_unur_FunctDefinition (struct parser_data *pdata);
static struct ftreenode *_unur_DefFunctDesignator (struct parser_data *pdata);
static struct ftreenode *_unur_DefParameterlist (struct parser_data *pdata, int *n_params);
static struct ftreenode *_unur_Expression (struct parser_data *pdata);
static struct ftreenode *_unur_SimpleExpression (struct parser_data *pdata);
static struct ftreenode *_unur_STerm (struct parser_data *pdata);
static struct ftreenode *_unur_Term (struct parser_data *pdata);
static struct ftreenode *_unur_Factor (struct parser_data *pdata);
static struct ftreenode *_unur_Bas_Exp (struct parser_data *pdata);
static struct ftreenode *_unur_FunctDesignator (struct parser_data *pdata);
static struct ftreenode *_unur_ActualParameterlist (struct parser_data *pdata, int n_params);
static struct ftreenode *_unur_fstr_simplification (const char *symb, int token,
						    struct ftreenode *left,
						    struct ftreenode *right);
static int _unur_fstr_reorganize (struct ftreenode *node);
static int _unur_fstr_next_token (struct parser_data *pdata, int *token, char **symbol);
static struct ftreenode *_unur_fstr_create_node (const char *symb, double val, int token,
						 struct ftreenode *left,
						 struct ftreenode *right);
static struct ftreenode *_unur_fstr_error_parse ( struct parser_data *pdata, int perrno, int line );
static const char *_unur_fstr_error_code ( int perrno );
static double _unur_fstr_eval_node (const struct ftreenode *node, double x);
static void _unur_fstr_error_deriv (const struct ftreenode *node, int line);
static int _unur_fstr_node2string ( struct unur_string *output, const struct ftreenode *node,
				    const char *variable, const char *function, int spaces );
static int _unur_fstr_print ( struct unur_string *output, const char *symb, double number );
#ifdef UNUR_ENABLE_LOGGING
static void _unur_fstr_debug_input ( const char *fstr );
static void _unur_fstr_debug_token ( const struct parser_data *pdata );
static void _unur_fstr_debug_tree ( const struct parser_data *pdata,
				    const struct ftreenode *root );
static void _unur_fstr_debug_show_tree (const struct parser_data *pdata,
					const struct ftreenode *node,
					int level, int location);
static void _unur_fstr_debug_deriv (const struct ftreenode *funct,
				    const struct ftreenode *deriv);
#endif
#include "functparser_init.ch"
#include "functparser_scanner.ch"
#include "functparser_parser.ch"
#include "functparser_eval.ch"
#include "functparser_deriv.ch"
#include "functparser_stringgen.ch"
#include "functparser_debug.ch"
