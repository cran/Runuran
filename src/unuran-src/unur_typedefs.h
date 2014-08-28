/* Copyright (c) 2000-2014 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifndef UNUR_TYPEDEFS_H_SEEN
#define UNUR_TYPEDEFS_H_SEEN
struct unur_distr;                       
typedef struct unur_distr UNUR_DISTR;
struct unur_par;                         
typedef struct unur_par   UNUR_PAR;
struct unur_gen;                         
typedef struct unur_gen   UNUR_GEN;
struct unur_urng;                        
typedef struct unur_urng  UNUR_URNG;
#define UNUR_URNG_UNURAN 1
typedef double UNUR_FUNCT_CONT  (double x, const struct unur_distr *distr);
typedef double UNUR_FUNCT_DISCR (int x, const struct unur_distr *distr);
typedef int    UNUR_IFUNCT_DISCR(double x, const struct unur_distr *distr);
typedef double UNUR_FUNCT_CVEC (const double *x, struct unur_distr *distr);
typedef int    UNUR_VFUNCT_CVEC(double *result, const double *x, struct unur_distr *distr);
typedef double UNUR_FUNCTD_CVEC(const double *x, int coord, struct unur_distr *distr);
struct unur_slist;         
typedef void UNUR_ERROR_HANDLER( const char *objid, const char *file, int line, 
				 const char *errortype, int unur_errno, const char *reason );
#endif  
