/* Copyright (c) 2000-2022 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_distr_cont {
  UNUR_FUNCT_CONT *pdf;         
  UNUR_FUNCT_CONT *dpdf;        
  UNUR_FUNCT_CONT *cdf;         
  UNUR_FUNCT_CONT *invcdf;      
  UNUR_FUNCT_CONT *logpdf;      
  UNUR_FUNCT_CONT *dlogpdf;     
  UNUR_FUNCT_CONT *logcdf;      
  UNUR_FUNCT_CONT *hr;          
  double norm_constant;         
  double params[UNUR_DISTR_MAXPARAMS];  
  int    n_params;              
  double *param_vecs[UNUR_DISTR_MAXPARAMS];  
  int    n_param_vec[UNUR_DISTR_MAXPARAMS]; 
  double mode;                  
  double center;                
  double area;                  
  double domain[2];             
  double trunc[2];              
  struct ftreenode *pdftree;    
  struct ftreenode *dpdftree;   
  struct ftreenode *logpdftree; 
  struct ftreenode *dlogpdftree;
  struct ftreenode *cdftree;    
  struct ftreenode *logcdftree; 
  struct ftreenode *hrtree;     
  int (*set_params)(struct unur_distr *distr, const double *params, int n_params );
  int (*upd_mode)(struct unur_distr *distr);
  int (*upd_area)(struct unur_distr *distr);
  int  (*init)(struct unur_par *par,struct unur_gen *gen);
};
struct unur_distr_cvec {
  UNUR_FUNCT_CVEC *pdf;         
  UNUR_VFUNCT_CVEC *dpdf;       
  UNUR_FUNCTD_CVEC *pdpdf;      
  UNUR_FUNCT_CVEC *logpdf;      
  UNUR_VFUNCT_CVEC *dlogpdf;    
  UNUR_FUNCTD_CVEC *pdlogpdf;   
  double *mean;                 
  double *covar;                
  double *cholesky;             
  double *covar_inv;            
  double *rankcorr;             
  double *rk_cholesky;          
  struct unur_distr **marginals; 
  double params[UNUR_DISTR_MAXPARAMS];  
  int    n_params;              
  double *param_vecs[UNUR_DISTR_MAXPARAMS];  
  int    n_param_vec[UNUR_DISTR_MAXPARAMS]; 
  double norm_constant;         
  double *mode;                 
  double *center;               
  double volume;                
  double *domainrect;           
#ifdef USE_DEPRECATED_CODE
  struct unur_distr **stdmarginals; 
#endif
  int (*upd_mode)(struct unur_distr *distr);
  int (*upd_volume)(struct unur_distr *distr);
  int (*init)(struct unur_gen *gen);
};
struct unur_distr_discr {
  double *pv;                   
  int     n_pv;                 
  UNUR_FUNCT_DISCR  *pmf;       
  UNUR_FUNCT_DISCR  *cdf;       
  UNUR_IFUNCT_DISCR *invcdf;    
  double params[UNUR_DISTR_MAXPARAMS];  
  int    n_params;              
  double norm_constant;         
  int    mode;                  
  double sum;                   
  int (*set_params)(struct unur_distr *distr, const double *params, int n_params );
  int (*upd_mode)(struct unur_distr *distr);
  int (*upd_sum)(struct unur_distr *distr);
  int domain[2];                
  int trunc[2];                 
  struct ftreenode *pmftree;    
  struct ftreenode *cdftree;    
  int  (*init)(struct unur_par *par,struct unur_gen *gen);
};
struct unur_distr_cemp {
  int     n_sample;             
  double *sample;               
  int     n_hist;               
  double *hist_prob;            
  double  hmin, hmax;           
  double *hist_bins;            
};
struct unur_distr_cvemp {
  double *sample;              
  int    n_sample;             
};
struct unur_distr_matr {
  int n_rows;                   
  int n_cols;                   
  int  (*init)(struct unur_par *par,struct unur_gen *gen);
};
struct unur_distr {
  union {             
    struct unur_distr_cont  cont;   
    struct unur_distr_matr  matr;   
    struct unur_distr_cvec  cvec;   
    struct unur_distr_discr discr;  
    struct unur_distr_cemp  cemp;   
    struct unur_distr_cvemp cvemp;  
  } data;                           
  unsigned type;                    
  unsigned id;                      
  const char *name;                 
  char *name_str;                   
  int dim;                          
  unsigned set;                     
  const void *extobj;               
  struct unur_distr *base;          
  void (*destroy)(struct unur_distr *distr); 
  struct unur_distr* (*clone)(const struct unur_distr *distr ); 
#ifdef UNUR_COOKIES
  unsigned cookie;                  
#endif
};
