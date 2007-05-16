/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

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
  int    type;                   
  int    info;                   
  double val;                    
  double (*vcalc)(double l, double r);
  struct ftreenode *(*dcalc)(const struct ftreenode *node, int *error); 
  unsigned (*node2C)(struct unur_string *output, const struct ftreenode *node, const char *variable );
  unsigned (*node2F)(struct unur_string *output, const struct ftreenode *node, const char *variable );
  unsigned (*node2J)(struct unur_string *output, const struct ftreenode *node, const char *variable );
};
static struct symbols symbol[] = {   
  {""    , S_NOSYMBOL, 0, 0.0 , v_dummy  , d_error, C_error  , F_error , J_error  },
  {"UCONST",S_UCONST , 9, 0.0 , v_const  , d_const, C_const  , F_const , J_const  },
  {"UFUNCT",S_UFUNCT , 0, 0.0 , v_dummy  , d_error, C_error  , F_error , J_error  },
  {"VAR" , S_UIDENT  , 9, 0.0 , v_dummy  , d_var  , C_var    , F_var   , J_var    },
  {"_ROS", S_NOSYMBOL, 0, 0.0 , v_dummy  , d_error, C_error  , F_error , J_error  },
  {"<"   , S_REL_OP  , 1, 0.0 , v_less   , d_const, C_infix  , F_lt    , J_lt     },
  {"="   , S_REL_OP  , 1, 0.0 , v_equal  , d_const, C_equal  , F_eq    , J_eq     },
  {"=="  , S_REL_OP  , 1, 0.0 , v_equal  , d_const, C_equal  , F_eq    , J_eq     },
  {">"   , S_REL_OP  , 1, 0.0 , v_greater, d_const, C_infix  , F_gt    , J_gt     },
  {"<="  , S_REL_OP  , 1, 0.0 , v_less_or, d_const, C_infix  , F_le    , J_le     },
  {"<>"  , S_REL_OP  , 1, 0.0 , v_unequal, d_const, C_unequal, F_ne    , J_ne     },
  {"!="  , S_REL_OP  , 1, 0.0 , v_unequal, d_const, C_unequal, F_ne    , J_ne     },
  {">="  , S_REL_OP  , 1, 0.0 , v_grtr_or, d_const, C_infix  , F_ge    , J_ge     },
  {"_NAS", S_NOSYMBOL, 0, 0.0 , v_dummy  , d_error, C_error  , F_error , J_error  },
  {"("   , S_OTHERS  , 0, 0.0 , v_dummy  , d_error, C_error  , F_error , J_error  },
  {")"   , S_OTHERS  , 0, 0.0 , v_dummy  , d_error, C_error  , F_error , J_error  },
  {","   , S_OTHERS  , 0, 0.0 , v_dummy  , d_error, C_error  , F_error , J_error  },
  {"+"   , S_ADD_OP  , 2, 0.0 , v_plus   , d_add  , C_infix  , F_infix , J_infix  },
  {"-"   , S_ADD_OP  , 2, 0.0 , v_minus  , d_add  , C_minus  , F_minus , J_minus  },
  {"*"   , S_MUL_OP  , 4, 0.0 , v_mul    , d_mul  , C_infix  , F_infix , J_infix  },
  {"/"   , S_MUL_OP  , 4, 0.0 , v_div    , d_div  , C_infix  , F_infix , J_infix  },
  {"^"   , S_HPR_OP  , 5, 0.0 , v_power  , d_power, C_power  , F_power , J_power  },
  {"_ANS", S_NOSYMBOL, 0, 0.0 , v_dummy  , d_error, C_error  , F_error , J_error  },
  {"pi"  , S_SCONST  , 9, M_PI, v_const  , d_const, C_const  , F_const , J_const  },
  {"e"   , S_SCONST  , 9, M_E , v_const  , d_const, C_const  , F_const , J_const  },
  {"mod" , S_SFUNCT  , 2, 0.0 , v_mod    , d_const, C_mod    , F_prefix, J_mod    },
  {"exp" , S_SFUNCT  , 1, 0.0 , v_exp    , d_exp  , C_prefix , F_prefix, J_prefix },
  {"log" , S_SFUNCT  , 1, 0.0 , v_log    , d_log  , C_prefix , F_prefix, J_prefix },
  {"sin" , S_SFUNCT  , 1, 0.0 , v_sin    , d_sin  , C_prefix , F_prefix, J_prefix },
  {"cos" , S_SFUNCT  , 1, 0.0 , v_cos    , d_cos  , C_prefix , F_prefix, J_prefix },
  {"tan" , S_SFUNCT  , 1, 0.0 , v_tan    , d_tan  , C_prefix , F_prefix, J_prefix },
  {"sec" , S_SFUNCT  , 1, 0.0 , v_sec    , d_sec  , C_sec    , F_sec   , J_sec    },
  {"sqrt", S_SFUNCT  , 1, 0.0 , v_sqrt   , d_sqrt , C_prefix , F_prefix, J_prefix },
  {"abs" , S_SFUNCT  , 1, 0.0 , v_abs    , d_abs  , C_abs    , F_prefix, J_prefix },
  {"sgn" , S_SFUNCT  , 1, 0.0 , v_sgn    , d_const, C_sgn    , F_sgn   , J_sgn    },
  {"_END", S_NOSYMBOL, 0, 0.0 , v_dummy  , d_error, C_error  , F_error , J_error },
};
static int s_uconst = 1;      
static int s_ufunct = 2;      
static int s_uident = 3;      
static int s_comma, s_minus, s_plus, s_mul, s_div, s_power;
static int _ros_start, _ros_end;    
static int _nas_start, _nas_end;    
static int _ans_start, _ans_end;    
static int _end;                    
