/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_slist {
  void **ptr;        
  int n_ptr;         
#ifdef UNUR_COOKIES
  unsigned cookie;   
#endif
};
