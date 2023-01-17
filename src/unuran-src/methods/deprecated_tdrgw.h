/* Copyright (c) 2000-2023 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#define unur_tdrgw_new( distr ) \
   (unur_ars_new( distr ))
#define unur_tdrgw_set_max_intervals( par,max_ivs ) \
   (unur_ars_set_max_intervals( (par),(max_ivs) ))
#define unur_tdrgw_set_cpoints( par,n_cpoints,cpoints ) \
   (unur_ars_set_cpoints( (par),(n_cpoints),(cpoints) ))
#define unur_tdrgw_set_reinit_percentiles( par,n_percentiles,percentiles ) \
   (unur_ars_set_reinit_percentiles( (par),(n_percentiles),(percentiles) ))
#define unur_tdrgw_chg_reinit_percentiles( gen,n_percentiles,percentiles ) \
   (unur_ars_chg_reinit_percentiles( (gen),(n_percentiles),(percentiles) ))
#define unur_tdrgw_set_reinit_ncpoints( par,ncpoints ) \
   (unur_ars_set_reinit_ncpoints( (par),(ncpoints) ))
#define unur_tdrgw_chg_reinit_ncpoints( gen,ncpoints ) \
   (unur_ars_chg_reinit_ncpoints( (gen),(ncpoints) ))
#define unur_tdrgw_set_verify( par,verify ) \
   (unur_ars_set_verify( (par),(verify) ))
#define unur_tdrgw_chg_verify( gen,verify ) \
   (unur_ars_chg_verify( (gen),(verify) ))
#define unur_tdrgw_set_pedantic( par,pedantic ) \
   (unur_ars_set_pedantic( (par),(pedantic) ))
#define unur_tdrgw_get_loghatarea( gen ) \
   (unur_ars_get_loghatarea( gen ))
#define unur_tdrgw_eval_invcdfhat( gen,u ) \
   (unur_ars_eval_invcdfhat( (gen),(u) ))
