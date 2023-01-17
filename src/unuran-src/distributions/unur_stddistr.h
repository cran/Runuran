/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef UNUR_STDDISTR_H_SEEN
#define UNUR_STDDISTR_H_SEEN
#define UNUR_DISTR_STD         (0x00000001u)  
enum {
  UNUR_DISTR_GENERIC          = 0x00000000u,
  UNUR_DISTR_CORDER           = 0x00000010u,  
  UNUR_DISTR_CXTRANS          = 0x00000020u,  
  UNUR_DISTR_CONDI            = 0x00000030u,  
  UNUR_DISTR_BETA             = 0x00000101u,  
  UNUR_DISTR_CAUCHY           = 0x00000201u,  
  UNUR_DISTR_CHI              = 0x00000301u,  
  UNUR_DISTR_CHISQUARE        = 0x00000401u,  
  UNUR_DISTR_EPANECHNIKOV     = 0x00000501u,  
  UNUR_DISTR_EXPONENTIAL      = 0x00000601u,  
  UNUR_DISTR_EXTREME_I        = 0x00000701u,  
  UNUR_DISTR_EXTREME_II       = 0x00000801u,  
  UNUR_DISTR_F                = 0x00000901u,  
  UNUR_DISTR_GAMMA            = 0x00000a01u,  
  UNUR_DISTR_GHYP             = 0x00002401u,  
  UNUR_DISTR_GIG              = 0x00000b01u,  
  UNUR_DISTR_GIG2             = 0x00002201u,  
  UNUR_DISTR_HYPERBOLIC       = 0x00002301u,  
  UNUR_DISTR_IG               = 0x00002101u,  
  UNUR_DISTR_LAPLACE          = 0x00000c01u,  
  UNUR_DISTR_LOGISTIC         = 0x00000d01u,  
  UNUR_DISTR_LOGNORMAL        = 0x00000e01u,  
  UNUR_DISTR_LOMAX            = 0x00000f01u,  
  UNUR_DISTR_MEIXNER          = 0x00002601u,  
  UNUR_DISTR_NORMAL           = 0x00001001u,  
   UNUR_DISTR_GAUSSIAN        = 0x00001001u,  
  UNUR_DISTR_PARETO           = 0x00001101u,  
  UNUR_DISTR_POWEREXPONENTIAL = 0x00001201u,  
  UNUR_DISTR_RAYLEIGH         = 0x00001301u,  
  UNUR_DISTR_SLASH            = 0x00001401u,  
  UNUR_DISTR_STUDENT          = 0x00001501u,  
  UNUR_DISTR_TRIANGULAR       = 0x00001601u,  
  UNUR_DISTR_UNIFORM          = 0x00002001u,  
   UNUR_DISTR_BOXCAR          = 0x00002001u,  
  UNUR_DISTR_VG               = 0x00002501u,  
  UNUR_DISTR_WEIBULL          = 0x00001801u,  
  UNUR_DISTR_BURR_I           = 0x0000b001u,  
  UNUR_DISTR_BURR_II          = 0x0000b101u,  
  UNUR_DISTR_BURR_III         = 0x0000b201u,  
  UNUR_DISTR_BURR_IV          = 0x0000b301u,  
  UNUR_DISTR_BURR_V           = 0x0000b401u,  
  UNUR_DISTR_BURR_VI          = 0x0000b501u,  
  UNUR_DISTR_BURR_VII         = 0x0000b601u,  
  UNUR_DISTR_BURR_VIII        = 0x0000b701u,  
  UNUR_DISTR_BURR_IX          = 0x0000b801u,  
  UNUR_DISTR_BURR_X           = 0x0000b901u,  
  UNUR_DISTR_BURR_XI          = 0x0000ba01u,  
  UNUR_DISTR_BURR_XII         = 0x0000bb01u,  
  UNUR_DISTR_BINOMIAL         = 0x00010001u,  
  UNUR_DISTR_GEOMETRIC        = 0x00020001u,  
  UNUR_DISTR_HYPERGEOMETRIC   = 0x00030001u,  
  UNUR_DISTR_LOGARITHMIC      = 0x00040001u,  
  UNUR_DISTR_NEGATIVEBINOMIAL = 0x00050001u,  
  UNUR_DISTR_POISSON          = 0x00060001u,  
  UNUR_DISTR_ZIPF             = 0x00070001u,  
  UNUR_DISTR_MCAUCHY          = 0x01000001u,  
  UNUR_DISTR_MNORMAL          = 0x02000001u,  
  UNUR_DISTR_MSTUDENT         = 0x03000001u,  
  UNUR_DISTR_MEXPONENTIAL     = 0x04000001u,  
  UNUR_DISTR_COPULA           = 0x05000001u,  
  UNUR_DISTR_MCORRELATION     = 0x10000001u   
};
#endif  
