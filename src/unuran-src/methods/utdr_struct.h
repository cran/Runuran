/* Copyright (c) 2000-2024 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_utdr_par { 
  double  fm;               
  double  hm;               
  double  c_factor;         
  double  delta_factor;     
};
struct unur_utdr_gen { 
  double  il;               
  double  ir;               
  double  fm;               
  double  hm;               
  double  vollc,volcompl,voll,
    al,ar,col,cor,sal,sar,bl,br,ttlx,ttrx,
    brblvolc,drar,dlal,ooar2,ooal2;
  double  c_factor;         
  double  delta_factor;     
};
