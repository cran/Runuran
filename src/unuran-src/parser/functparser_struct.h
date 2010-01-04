/* Copyright (c) 2000-2010 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct ftreenode { 
  char            *symbol;      
  int             token;        
  int             type;         
  double          val;          
  struct ftreenode *left;       
  struct ftreenode *right;      
#ifdef UNUR_COOKIES
  unsigned cookie;              
#endif
}; 
