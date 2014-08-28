/* Copyright (c) 2000-2014 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#define _unur_init(par)               (par)->init(par)
#define _unur_sample_discr(gen)       (gen)->sample.discr(gen)
#define _unur_sample_cont(gen)        (gen)->sample.cont(gen)
#define _unur_sample_vec(gen,vector)  (gen)->sample.cvec(gen,vector)
#define _unur_free(gen)               do {if(gen) (gen)->destroy(gen);} while(0)
#define _unur_gen_is_discr(gen) ( ((gen)->distr->type == UNUR_DISTR_DISCR) ? 1 : 0 )
#define _unur_gen_is_cont(gen)  ( ((gen)->distr->type == UNUR_DISTR_CONT)  ? 1 : 0 )
#define _unur_gen_is_vec(gen)   ( ((gen)->distr->type == UNUR_DISTR_CVEC)  ? 1 : 0 )
int _unur_sample_discr_error( struct unur_gen *gen );
double _unur_sample_cont_error( struct unur_gen *gen );
int _unur_sample_cvec_error( struct unur_gen *gen, double *vec );
int _unur_sample_matr_error( struct unur_gen *gen, double *mat );
struct unur_par *_unur_par_new( size_t s );
struct unur_par *_unur_par_clone( const struct unur_par *par );
#define _unur_par_free(par)  do {free((par)->datap); free(par);} while(0)
struct unur_gen *_unur_generic_create( struct unur_par *par, size_t s );
struct unur_gen *_unur_generic_clone( const struct unur_gen *gen, const char *type );
#define _unur_gen_clone(gen)    ((gen)->clone(gen))
void _unur_generic_free( struct unur_gen *gen );
struct unur_gen **_unur_gen_list_set( struct unur_gen *gen, int n_gen_list );
struct unur_gen **_unur_gen_list_clone( struct unur_gen **gen_list, int n_gen_list );
void _unur_gen_list_free( struct unur_gen **gen_list, int n_gen_list );
