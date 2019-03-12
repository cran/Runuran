/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#if defined(PARSER) && defined(CODEGEN)
#  error \
   You have to #define either 'PARSER' or 'CODEGEN' before including this header file!
#endif
#if !defined(PARSER) && !defined(CODEGEN)
#  error \
   You have to #define either 'PARSER' or 'CODEGEN' before including this header file!
#endif
#ifdef PARSER
static double v_dummy  (double l, double r);
static double v_less   (double l, double r);
static double v_equal  (double l, double r);
static double v_greater(double l, double r);
static double v_less_or(double l, double r);
static double v_unequal(double l, double r);
static double v_grtr_or(double l, double r);
static double v_plus   (double l, double r); 
static double v_minus  (double l, double r); 
static double v_mul    (double l, double r);
static double v_div    (double l, double r);
static double v_mod    (double l, double r);
static double v_power  (double l, double r);
static double v_const  (double l, double r);
static double v_exp    (double l, double r);
static double v_log    (double l, double r);
static double v_sin    (double l, double r); 
static double v_cos    (double l, double r);
static double v_tan    (double l, double r);
static double v_sec    (double l, double r); 
static double v_sqrt   (double l, double r);
static double v_abs    (double l, double r);
static double v_sgn    (double l, double r);
#endif
#ifdef PARSER
static struct ftreenode *d_error ( const struct ftreenode *node, int *error );
static struct ftreenode *d_const ( const struct ftreenode *node, int *error );
static struct ftreenode *d_var   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_add   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_mul   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_div   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_power ( const struct ftreenode *node, int *error );
static struct ftreenode *d_exp   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_log   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_sin   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_cos   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_tan   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_sec   ( const struct ftreenode *node, int *error );
static struct ftreenode *d_sqrt  ( const struct ftreenode *node, int *error );
static struct ftreenode *d_abs   ( const struct ftreenode *node, int *error );
#endif
#ifdef CODEGEN
static unsigned C_error  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned C_const  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned C_var    ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned C_prefix_generic ( struct unur_string *output, const struct ftreenode *node, 
				   const char *variable, const char *symbol );
static unsigned C_prefix ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned C_power  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned C_sec    ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned C_abs    ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned C_sgn    ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned C_infix_generic  ( struct unur_string *output, const struct ftreenode *node,
				   const char *variable, const char *symbol );
static unsigned C_infix  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned C_equal  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned C_unequal( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned C_minus  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned C_mod    ( struct unur_string *output, const struct ftreenode *node, const char *variable );
#endif
#ifdef CODEGEN
static unsigned F_error  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_const  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_var    ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_prefix_generic ( struct unur_string *output, const struct ftreenode *node, 
				   const char *variable, const char *symbol );
static unsigned F_prefix ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_lt     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_le     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_gt     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_ge     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_eq     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_ne     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_sec    ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_sgn    ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_infix_generic  ( struct unur_string *output, const struct ftreenode *node,
				   const char *variable, const char *symbol );
static unsigned F_infix  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_minus  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned F_power  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
#endif
#ifdef CODEGEN
static unsigned J_error  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_const  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_var    ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_prefix_generic ( struct unur_string *output, const struct ftreenode *node, 
				   const char *variable, const char *symbol );
static unsigned J_prefix ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_lt     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_le     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_gt     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_ge     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_eq     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_ne     ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_power  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_sec    ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_sgn    ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_infix_generic  ( struct unur_string *output, const struct ftreenode *node,
				   const char *variable, const char *symbol );
static unsigned J_infix  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_minus  ( struct unur_string *output, const struct ftreenode *node, const char *variable );
static unsigned J_mod    ( struct unur_string *output, const struct ftreenode *node, const char *variable );
#endif
enum {
  S_NOSYMBOL = 0,    
  S_SFUNCT,          
  S_SCONST,          
  S_UIDENT,          
  S_UFUNCT,          
  S_UCONST,          
  S_REL_OP,          
  S_ADD_OP,          
  S_MUL_OP,          
  S_HPR_OP,          
  S_OTHERS           
};
#define SYMBLENGTH 10           
struct symbols { 
  char   name[SYMBLENGTH];        
#ifdef PARSER
  int    type;                   
  int    info;                   
  double val;                    
  double (*vcalc)(double l, double r);
  struct ftreenode *(*dcalc)(const struct ftreenode *node, int *error); 
#endif
#ifdef CODEGEN
  unsigned (*node2C)(struct unur_string *output, const struct ftreenode *node, const char *variable );
  unsigned (*node2F)(struct unur_string *output, const struct ftreenode *node, const char *variable );
  unsigned (*node2J)(struct unur_string *output, const struct ftreenode *node, const char *variable );
#endif
};
#ifdef PARSER
#  define S(name,type,priority,value,funct,deriv) \
   name,type,priority,value,funct,deriv
#else
#  define S(name,type,priority,value,funct,deriv) \
   name,
#endif
#ifdef CODEGEN
#  define CG(cc,fortran,java)  cc,fortran,java
#else
#  define CG(cc,fortran,java)
#endif
static struct symbols symbol[] = {   
  {S (""    , S_NOSYMBOL, 0, 0.0 , v_dummy  , d_error) CG (C_error  , F_error , J_error  )},
  {S ("UCONST",S_UCONST , 9, 0.0 , v_const  , d_const) CG (C_const  , F_const , J_const  )},
  {S ("UFUNCT",S_UFUNCT , 0, 0.0 , v_dummy  , d_error) CG (C_error  , F_error , J_error  )},
  {S ("VAR" , S_UIDENT  , 9, 0.0 , v_dummy  , d_var  ) CG (C_var    , F_var   , J_var    )},
  {S ("_ROS", S_NOSYMBOL, 0, 0.0 , v_dummy  , d_error) CG (C_error  , F_error , J_error  )},
  {S ("<"   , S_REL_OP  , 1, 0.0 , v_less   , d_const) CG (C_infix  , F_lt    , J_lt     )},
  {S ("="   , S_REL_OP  , 1, 0.0 , v_equal  , d_const) CG (C_equal  , F_eq    , J_eq     )},
  {S ("=="  , S_REL_OP  , 1, 0.0 , v_equal  , d_const) CG (C_equal  , F_eq    , J_eq     )},
  {S (">"   , S_REL_OP  , 1, 0.0 , v_greater, d_const) CG (C_infix  , F_gt    , J_gt     )},
  {S ("<="  , S_REL_OP  , 1, 0.0 , v_less_or, d_const) CG (C_infix  , F_le    , J_le     )},
  {S ("<>"  , S_REL_OP  , 1, 0.0 , v_unequal, d_const) CG (C_unequal, F_ne    , J_ne     )},
  {S ("!="  , S_REL_OP  , 1, 0.0 , v_unequal, d_const) CG (C_unequal, F_ne    , J_ne     )},
  {S (">="  , S_REL_OP  , 1, 0.0 , v_grtr_or, d_const) CG (C_infix  , F_ge    , J_ge     )},
  {S ("_NAS", S_NOSYMBOL, 0, 0.0 , v_dummy  , d_error) CG (C_error  , F_error , J_error  )},
  {S ("("   , S_OTHERS  , 0, 0.0 , v_dummy  , d_error) CG (C_error  , F_error , J_error  )},
  {S (")"   , S_OTHERS  , 0, 0.0 , v_dummy  , d_error) CG (C_error  , F_error , J_error  )},
  {S (","   , S_OTHERS  , 0, 0.0 , v_dummy  , d_error) CG (C_error  , F_error , J_error  )},
  {S ("+"   , S_ADD_OP  , 2, 0.0 , v_plus   , d_add  ) CG (C_infix  , F_infix , J_infix  )},
  {S ("-"   , S_ADD_OP  , 2, 0.0 , v_minus  , d_add  ) CG (C_minus  , F_minus , J_minus  )},
  {S ("*"   , S_MUL_OP  , 4, 0.0 , v_mul    , d_mul  ) CG (C_infix  , F_infix , J_infix  )},
  {S ("/"   , S_MUL_OP  , 4, 0.0 , v_div    , d_div  ) CG (C_infix  , F_infix , J_infix  )},
  {S ("^"   , S_HPR_OP  , 5, 0.0 , v_power  , d_power) CG (C_power  , F_power , J_power  )},
  {S ("_ANS", S_NOSYMBOL, 0, 0.0 , v_dummy  , d_error) CG (C_error  , F_error , J_error  )},
  {S ("pi"  , S_SCONST  , 9, M_PI, v_const  , d_const) CG (C_const  , F_const , J_const  )},
  {S ("e"   , S_SCONST  , 9, M_E , v_const  , d_const) CG (C_const  , F_const , J_const  )},
  {S ("mod" , S_SFUNCT  , 2, 0.0 , v_mod    , d_const) CG (C_mod    , F_prefix, J_mod    )},
  {S ("exp" , S_SFUNCT  , 1, 0.0 , v_exp    , d_exp  ) CG (C_prefix , F_prefix, J_prefix )},
  {S ("log" , S_SFUNCT  , 1, 0.0 , v_log    , d_log  ) CG (C_prefix , F_prefix, J_prefix )},
  {S ("sin" , S_SFUNCT  , 1, 0.0 , v_sin    , d_sin  ) CG (C_prefix , F_prefix, J_prefix )},
  {S ("cos" , S_SFUNCT  , 1, 0.0 , v_cos    , d_cos  ) CG (C_prefix , F_prefix, J_prefix )},
  {S ("tan" , S_SFUNCT  , 1, 0.0 , v_tan    , d_tan  ) CG (C_prefix , F_prefix, J_prefix )},
  {S ("sec" , S_SFUNCT  , 1, 0.0 , v_sec    , d_sec  ) CG (C_sec    , F_sec   , J_sec    )},
  {S ("sqrt", S_SFUNCT  , 1, 0.0 , v_sqrt   , d_sqrt ) CG (C_prefix , F_prefix, J_prefix )},
  {S ("abs" , S_SFUNCT  , 1, 0.0 , v_abs    , d_abs  ) CG (C_abs    , F_prefix, J_prefix )},
  {S ("sgn" , S_SFUNCT  , 1, 0.0 , v_sgn    , d_const) CG (C_sgn    , F_sgn   , J_sgn    )},
  {S ("_END", S_NOSYMBOL, 0, 0.0 , v_dummy  , d_error) CG (C_error  , F_error , J_error )},
};
#undef S
#undef GC
#ifdef PARSER
static int s_uconst = 1;      
static int s_ufunct = 2;      
static int s_uident = 3;      
static int s_comma, s_minus, s_plus, s_mul, s_div, s_power;
static int _ros_start, _ros_end;    
static int _nas_start, _nas_end;    
static int _ans_start, _ans_end;    
static int _end;                    
#endif
