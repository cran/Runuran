/* Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_mvtdr_par { 
  int max_cones;                  
  int steps_min;                  
  double bound_splitting;         
};
typedef struct s_vertex           
{
  struct s_vertex *next;          
  int index;                      
  double *coord;                  
  double norm;                    
} VERTEX;
typedef struct s_cone             
{
  struct s_cone *next;            
  int level;                      
  VERTEX **v;                     
  double *center;                 
  double logdetf;                 
  double alpha;                   
  double beta;                    
  double *gv;                     
  double logai;                   
  double tp;                      
  double Hi;                      
  double Hsum;                    
  double Tfp;                     
  double height;                  
} CONE;
typedef struct s_edge_table       
{
  int  index[2];                  
  VERTEX *vertex;                 
  struct s_edge_table *next;      
} E_TABLE;
typedef struct s_tp_arg           
{
  double t;                       
  double logH;                    
  CONE *c;                        
  UNUR_GEN *gen;                  
  int status;                     
} TP_ARG;
enum {                            
  MVTDR_CONE_OK      = 0x000,     
  MVTDR_CONE_DOMAIN  = 0x001,     
  MVTDR_CONE_INVALID = 0x002      
};
struct unur_mvtdr_gen { 
  int  dim;                       
  int  has_domain;                
  double max_gamma;               
  const double *center;           
  CONE *cone;                     
  CONE *last_cone;                
  int n_cone;                     
  int max_cones;                  
  double bound_splitting;         
  VERTEX *vertex;                 
  VERTEX *last_vertex;            
  int n_vertex;                   
  E_TABLE **etable;               
  int etable_size;                
  CONE **guide;                   
  int guide_size;                 
  double *S;                      
  double *g;                      
  double *tp_coord;               
  double *tp_mcoord;              
  double *tp_Tgrad;               
  double Htot;                    
  int steps_min;                  
  int n_steps;                    
  double pdfcenter;               
};
