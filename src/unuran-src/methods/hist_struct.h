/* Copyright (c) 2000-2014 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_hist_par {
  int dummy;
};
struct unur_hist_gen {
  int     n_hist;       
  double *prob;         
  double *bins;         
  double  hmin, hmax;   
  double  hwidth;       
  double  sum;          
  double *cumpv;        
  int    *guide_table;  
};
