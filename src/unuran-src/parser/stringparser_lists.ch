/* Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

struct unur_distr *
_unur_str_distr_new( char *distribution )
{
  struct unur_distr *distr = NULL;    
  char distr_unknown;
#ifdef UNUR_ENABLE_LOGGING
  char *name;                         
#endif
  char *params;                       
  double *darray = NULL;              
  int n_darray = 0;                   
#ifdef UNUR_ENABLE_LOGGING
  name = distribution;
#endif
  params = strchr(distribution,'(');
  if (params != NULL) {
    *params = '\0';                   
    ++params;                         
  }
  n_darray = _unur_parse_dlist(params, &darray );
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_distr(1,name,darray,n_darray);
#endif
  distr = (struct unur_distr *) &distr_unknown;
	 switch (*distribution) {
	 case 'b':
		 if ( !strcmp( distribution, "beta") ) {
			 distr = unur_distr_beta (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "binomial") ) {
			 distr = unur_distr_binomial (darray,n_darray);
			 break;
		 }
		 break;
	 case 'c':
		 if ( !strcmp( distribution, "cauchy") ) {
			 distr = unur_distr_cauchy (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "chi") ) {
			 distr = unur_distr_chi (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "chisquare") ) {
			 distr = unur_distr_chisquare (darray,n_darray);
			 break;
		 }
		 break;
	 case 'e':
		 if ( !strcmp( distribution, "exponential") ) {
			 distr = unur_distr_exponential (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "extremei") ) {
			 distr = unur_distr_extremeI (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "extremeii") ) {
			 distr = unur_distr_extremeII (darray,n_darray);
			 break;
		 }
		 break;
	 case 'f':
		 if ( !strcmp( distribution, "f") ) {
			 distr = unur_distr_F (darray,n_darray);
			 break;
		 }
		 break;
	 case 'g':
		 if ( !strcmp( distribution, "gamma") ) {
			 distr = unur_distr_gamma (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "geometric") ) {
			 distr = unur_distr_geometric (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "gig") ) {
			 distr = unur_distr_gig (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "gig2") ) {
			 distr = unur_distr_gig2 (darray,n_darray);
			 break;
		 }
		 break;
	 case 'h':
		 if ( !strcmp( distribution, "hyperbolic") ) {
			 distr = unur_distr_hyperbolic (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "hypergeometric") ) {
			 distr = unur_distr_hypergeometric (darray,n_darray);
			 break;
		 }
		 break;
	 case 'i':
		 if ( !strcmp( distribution, "ig") ) {
			 distr = unur_distr_ig (darray,n_darray);
			 break;
		 }
		 break;
	 case 'l':
		 if ( !strcmp( distribution, "laplace") ) {
			 distr = unur_distr_laplace (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "logarithmic") ) {
			 distr = unur_distr_logarithmic (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "logistic") ) {
			 distr = unur_distr_logistic (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "lognormal") ) {
			 distr = unur_distr_lognormal (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "lomax") ) {
			 distr = unur_distr_lomax (darray,n_darray);
			 break;
		 }
		 break;
	 case 'n':
		 if ( !strcmp( distribution, "negativebinomial") ) {
			 distr = unur_distr_negativebinomial (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "normal") ) {
			 distr = unur_distr_normal (darray,n_darray);
			 break;
		 }
		 break;
	 case 'p':
		 if ( !strcmp( distribution, "pareto") ) {
			 distr = unur_distr_pareto (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "poisson") ) {
			 distr = unur_distr_poisson (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "powerexponential") ) {
			 distr = unur_distr_powerexponential (darray,n_darray);
			 break;
		 }
		 break;
	 case 'r':
		 if ( !strcmp( distribution, "rayleigh") ) {
			 distr = unur_distr_rayleigh (darray,n_darray);
			 break;
		 }
		 break;
	 case 's':
		 if ( !strcmp( distribution, "slash") ) {
			 distr = unur_distr_slash (darray,n_darray);
			 break;
		 }
		 if ( !strcmp( distribution, "student") ) {
			 distr = unur_distr_student (darray,n_darray);
			 break;
		 }
		 break;
	 case 't':
		 if ( !strcmp( distribution, "triangular") ) {
			 distr = unur_distr_triangular (darray,n_darray);
			 break;
		 }
		 break;
	 case 'u':
		 if ( !strcmp( distribution, "uniform") ) {
			 distr = unur_distr_uniform (darray,n_darray);
			 break;
		 }
		 break;
	 case 'w':
		 if ( !strcmp( distribution, "weibull") ) {
			 distr = unur_distr_weibull (darray,n_darray);
			 break;
		 }
	 }
	 if (distr == (struct unur_distr *) &distr_unknown) { 
		 do {
			 if ( !strcmp( distribution, "cemp") ) {
				 distr = unur_distr_cemp_new();
				 break;
			 }
			 if ( !strcmp( distribution, "cont") ) {
				 distr = unur_distr_cont_new();
				 break;
			 }
			 if ( !strcmp( distribution, "discr") ) {
				 distr = unur_distr_discr_new();
				 break;
			 }
		 } while (0);
	 }
   if (distr == (struct unur_distr *) &distr_unknown) {
     _unur_error_unknown(distribution,"distribution");
     distr = NULL;
   }
   else if (distr == NULL) {
     _unur_error_invalid(distribution,"distribution");
   }
   if (darray) free(darray);
   return distr;
} 
int
_unur_str_distr_set( UNUR_DISTR **ptr_distr, const char *key, char *value )
{
  int result;      
  struct unur_distr *distr = *ptr_distr;
  char type_args[MAX_SET_ARGS+1];  
  char *args[MAX_SET_ARGS+1];      
  if (_unur_str_set_args( value, type_args, args, MAX_SET_ARGS ) < 0) {
    return UNUR_ERR_STR_SYNTAX;
  }
  result = UNUR_ERR_STR_UNKNOWN;
	 switch (distr->type) {
	 case UNUR_DISTR_CEMP:
		 switch (*key) {
		 case 'd':
			 if ( !strcmp(key, "data") ) {
				 result = _unur_str_distr_set_Di(distr,key,type_args,args,unur_distr_cemp_set_data);
				 break;
			 }
			 break;
		 case 'h':
			 if ( !strcmp(key, "hist_bins") ) {
				 result = _unur_str_distr_set_Di(distr,key,type_args,args,unur_distr_cemp_set_hist_bins);
				 break;
			 }
			 if ( !strcmp(key, "hist_domain") ) {
				 result = _unur_str_distr_set_dd(distr,key,type_args,args,unur_distr_cemp_set_hist_domain);
				 break;
			 }
			 if ( !strcmp(key, "hist_prob") ) {
				 result = _unur_str_distr_set_Di(distr,key,type_args,args,unur_distr_cemp_set_hist_prob);
				 break;
			 }
		 }
		 break;
	 case UNUR_DISTR_CONT:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "cdf") ) {
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_cont_set_cdfstr);
				 break;
			 }
			 if ( !strcmp(key, "cdfstr") ) {
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_cont_set_cdfstr);
				 break;
			 }
			 if ( !strcmp(key, "center") ) {
				 result = _unur_str_distr_set_d(distr,key,type_args,args,unur_distr_cont_set_center);
				 break;
			 }
			 break;
		 case 'd':
			 if ( !strcmp(key, "domain") ) {
				 result = _unur_str_distr_set_dd(distr,key,type_args,args,unur_distr_cont_set_domain);
				 break;
			 }
			 break;
		 case 'h':
			 if ( !strcmp(key, "hr") ) {
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_cont_set_hrstr);
				 break;
			 }
			 if ( !strcmp(key, "hrstr") ) {
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_cont_set_hrstr);
				 break;
			 }
			 break;
		 case 'l':
			 if ( !strcmp(key, "logcdf") ) {
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_cont_set_logcdfstr);
				 break;
			 }
			 if ( !strcmp(key, "logcdfstr") ) {
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_cont_set_logcdfstr);
				 break;
			 }
			 if ( !strcmp(key, "logpdf") ) {
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_cont_set_logpdfstr);
				 break;
			 }
			 if ( !strcmp(key, "logpdfstr") ) {
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_cont_set_logpdfstr);
				 break;
			 }
			 break;
		 case 'm':
			 if ( !strcmp(key, "mode") ) {
				 result = _unur_str_distr_set_d(distr,key,type_args,args,unur_distr_cont_set_mode);
				 break;
			 }
			 break;
		 case 'p':
			 if ( !strcmp(key, "pdf") ) {
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_cont_set_pdfstr);
				 break;
			 }
			 if ( !strcmp(key, "pdfarea") ) {
				 result = _unur_str_distr_set_d(distr,key,type_args,args,unur_distr_cont_set_pdfarea);
				 break;
			 }
			 if ( !strcmp(key, "pdfparams") ) {
				 result = _unur_str_distr_set_Di(distr,key,type_args,args,unur_distr_cont_set_pdfparams);
				 break;
			 }
			 if ( !strcmp(key, "pdfstr") ) {
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_cont_set_pdfstr);
				 break;
			 }
		 }
		 break;
	 case UNUR_DISTR_DISCR:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "cdf") ) {
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_discr_set_cdfstr);
				 break;
			 }
			 if ( !strcmp(key, "cdfstr") ) {
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_discr_set_cdfstr);
				 break;
			 }
			 break;
		 case 'd':
			 if ( !strcmp(key, "domain") ) {
				 result = _unur_str_distr_set_ii(distr,key,type_args,args,unur_distr_discr_set_domain);
				 break;
			 }
			 break;
		 case 'm':
			 if ( !strcmp(key, "mode") ) {
				 result = _unur_str_distr_set_i(distr,key,type_args,args,unur_distr_discr_set_mode);
				 break;
			 }
			 break;
		 case 'p':
			 if ( !strcmp(key, "pmf") ) {
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_discr_set_pmfstr);
				 break;
			 }
			 if ( !strcmp(key, "pmfparams") ) {
				 result = _unur_str_distr_set_Di(distr,key,type_args,args,unur_distr_discr_set_pmfparams);
				 break;
			 }
			 if ( !strcmp(key, "pmfstr") ) {
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_discr_set_pmfstr);
				 break;
			 }
			 if ( !strcmp(key, "pmfsum") ) {
				 result = _unur_str_distr_set_d(distr,key,type_args,args,unur_distr_discr_set_pmfsum);
				 break;
			 }
			 if ( !strcmp(key, "pv") ) {
				 result = _unur_str_distr_set_Di(distr,key,type_args,args,unur_distr_discr_set_pv);
				 break;
			 }
		 }
		 break;
	 }
	 if (result == UNUR_ERR_STR_UNKNOWN) {
		 switch (*key) {
		 case 'n':
			 if ( !strcmp(key, "name") ) {
				 result = _unur_str_distr_set_C(distr,key,type_args,args,unur_distr_set_name);
				 break;
			 }
		 }
	 }
   if (result == UNUR_ERR_STR_UNKNOWN)
     if (distr->type == UNUR_DISTR_CONT) {
       if ( !strcmp(key, "orderstatistics") ) {
	 *ptr_distr = _unur_str_distr_make_os (distr, key, type_args, args);
	 result = (*ptr_distr == NULL) ? UNUR_ERR_STR_SYNTAX : UNUR_SUCCESS;
       }
     }
  if (result == UNUR_ERR_STR_UNKNOWN) {
    _unur_error_unknown(key,"parameter for given distribution");
    return UNUR_ERR_STR_UNKNOWN;
  }
  else if (result != UNUR_SUCCESS) {
    _unur_error_invalid(key,"set call");
    return UNUR_ERR_STR_SYNTAX;
  }
  return UNUR_SUCCESS;
} 
struct unur_par *
_unur_str_par_new( const char *method, const UNUR_DISTR *distr )
{
  struct unur_par *par = NULL;
  char method_unknown;
#ifdef UNUR_ENABLE_LOGGING
  if (_unur_default_debugflag & UNUR_DEBUG_SETUP)
    _unur_str_debug_string(1,"method",method);
#endif
  par = (struct unur_par *) &method_unknown;
	 switch (*method) {
	 case 'a':
		 if ( !strcmp( method, "arou") ) {
			 par = unur_arou_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "ars") ) {
			 par = unur_ars_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "auto") ) {
			 par = unur_auto_new(distr);
			 break;
		 }
		 break;
	 case 'c':
		 if ( !strcmp( method, "cstd") ) {
			 par = unur_cstd_new(distr);
			 break;
		 }
		 break;
	 case 'd':
		 if ( !strcmp( method, "dari") ) {
			 par = unur_dari_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "dau") ) {
			 par = unur_dau_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "dgt") ) {
			 par = unur_dgt_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "dsrou") ) {
			 par = unur_dsrou_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "dss") ) {
			 par = unur_dss_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "dstd") ) {
			 par = unur_dstd_new(distr);
			 break;
		 }
		 break;
	 case 'e':
		 if ( !strcmp( method, "empk") ) {
			 par = unur_empk_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "empl") ) {
			 par = unur_empl_new(distr);
			 break;
		 }
		 break;
	 case 'g':
		 if ( !strcmp( method, "gibbs") ) {
			 par = unur_gibbs_new(distr);
			 break;
		 }
		 break;
	 case 'h':
		 if ( !strcmp( method, "hinv") ) {
			 par = unur_hinv_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "hist") ) {
			 par = unur_hist_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "hitro") ) {
			 par = unur_hitro_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "hrb") ) {
			 par = unur_hrb_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "hrd") ) {
			 par = unur_hrd_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "hri") ) {
			 par = unur_hri_new(distr);
			 break;
		 }
		 break;
	 case 'i':
		 if ( !strcmp( method, "itdr") ) {
			 par = unur_itdr_new(distr);
			 break;
		 }
		 break;
	 case 'm':
		 if ( !strcmp( method, "mcorr") ) {
			 par = unur_mcorr_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "mvstd") ) {
			 par = unur_mvstd_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "mvtdr") ) {
			 par = unur_mvtdr_new(distr);
			 break;
		 }
		 break;
	 case 'n':
		 if ( !strcmp( method, "ninv") ) {
			 par = unur_ninv_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "norta") ) {
			 par = unur_norta_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "nrou") ) {
			 par = unur_nrou_new(distr);
			 break;
		 }
		 break;
	 case 'p':
		 if ( !strcmp( method, "pinv") ) {
			 par = unur_pinv_new(distr);
			 break;
		 }
		 break;
	 case 's':
		 if ( !strcmp( method, "srou") ) {
			 par = unur_srou_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "ssr") ) {
			 par = unur_ssr_new(distr);
			 break;
		 }
		 break;
	 case 't':
		 if ( !strcmp( method, "tabl") ) {
			 par = unur_tabl_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "tdr") ) {
			 par = unur_tdr_new(distr);
			 break;
		 }
		 break;
	 case 'u':
		 if ( !strcmp( method, "unif") ) {
			 par = unur_unif_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "utdr") ) {
			 par = unur_utdr_new(distr);
			 break;
		 }
		 break;
	 case 'v':
		 if ( !strcmp( method, "vempk") ) {
			 par = unur_vempk_new(distr);
			 break;
		 }
		 if ( !strcmp( method, "vnrou") ) {
			 par = unur_vnrou_new(distr);
			 break;
		 }
	 }
   if (par == (struct unur_par *) &method_unknown) {
     _unur_error_unknown(method,"method");
     par = NULL;
   }
   else if (par == NULL) {
     _unur_error_invalid(method,"method");
   }
   return par;
} 
int
_unur_str_par_set( UNUR_PAR *par, const char *key, char *value, struct unur_slist *mlist )
{
  int result = 0;                  
  char type_args[MAX_SET_ARGS+1];  
  char *args[MAX_SET_ARGS+1];      
  if (_unur_str_set_args( value, type_args, args, MAX_SET_ARGS ) < 0) {
    return UNUR_ERR_STR_SYNTAX;
  }
  result = UNUR_ERR_STR_UNKNOWN;
	 switch (par->method) {
	 case UNUR_METH_AROU:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "cpoints") ) {
				 result = _unur_str_par_set_iD(par,key,type_args,args,unur_arou_set_cpoints,mlist);
				 break;
			 }
			 break;
		 case 'd':
			 if ( !strcmp(key, "darsfactor") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_arou_set_darsfactor);
				 break;
			 }
			 break;
		 case 'g':
			 if ( !strcmp(key, "guidefactor") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_arou_set_guidefactor);
				 break;
			 }
			 break;
		 case 'm':
			 if ( !strcmp(key, "max_segments") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_arou_set_max_segments);
				 break;
			 }
			 if ( !strcmp(key, "max_sqhratio") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_arou_set_max_sqhratio);
				 break;
			 }
			 break;
		 case 'p':
			 if ( !strcmp(key, "pedantic") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_arou_set_pedantic);
				 break;
			 }
			 break;
		 case 'u':
			 if ( !strcmp(key, "usecenter") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_arou_set_usecenter);
				 break;
			 }
			 if ( !strcmp(key, "usedars") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_arou_set_usedars);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_arou_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_ARS:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "cpoints") ) {
				 result = _unur_str_par_set_iD(par,key,type_args,args,unur_ars_set_cpoints,mlist);
				 break;
			 }
			 break;
		 case 'm':
			 if ( !strcmp(key, "max_intervals") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_ars_set_max_intervals);
				 break;
			 }
			 if ( !strcmp(key, "max_iter") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_ars_set_max_iter);
				 break;
			 }
			 break;
		 case 'p':
			 if ( !strcmp(key, "pedantic") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_ars_set_pedantic);
				 break;
			 }
			 break;
		 case 'r':
			 if ( !strcmp(key, "reinit_ncpoints") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_ars_set_reinit_ncpoints);
				 break;
			 }
			 if ( !strcmp(key, "reinit_percentiles") ) {
				 result = _unur_str_par_set_iD(par,key,type_args,args,unur_ars_set_reinit_percentiles,mlist);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_ars_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_AUTO:
		 switch (*key) {
		 case 'l':
			 if ( !strcmp(key, "logss") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_auto_set_logss);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_CSTD:
		 switch (*key) {
		 case 'v':
			 if ( !strcmp(key, "variant") ) {
				 result = _unur_str_par_set_u(par,key,type_args,args,unur_cstd_set_variant);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_DARI:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "cpfactor") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_dari_set_cpfactor);
				 break;
			 }
			 break;
		 case 's':
			 if ( !strcmp(key, "squeeze") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_dari_set_squeeze);
				 break;
			 }
			 break;
		 case 't':
			 if ( !strcmp(key, "tablesize") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_dari_set_tablesize);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_dari_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_DAU:
		 switch (*key) {
		 case 'u':
			 if ( !strcmp(key, "urnfactor") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_dau_set_urnfactor);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_DGT:
		 switch (*key) {
		 case 'g':
			 if ( !strcmp(key, "guidefactor") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_dgt_set_guidefactor);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "variant") ) {
				 result = _unur_str_par_set_u(par,key,type_args,args,unur_dgt_set_variant);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_DSROU:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "cdfatmode") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_dsrou_set_cdfatmode);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_dsrou_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_DSTD:
		 switch (*key) {
		 case 'v':
			 if ( !strcmp(key, "variant") ) {
				 result = _unur_str_par_set_u(par,key,type_args,args,unur_dstd_set_variant);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_EMPK:
		 switch (*key) {
		 case 'b':
			 if ( !strcmp(key, "beta") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_empk_set_beta);
				 break;
			 }
			 break;
		 case 'k':
			 if ( !strcmp(key, "kernel") ) {
				 result = _unur_str_par_set_u(par,key,type_args,args,unur_empk_set_kernel);
				 break;
			 }
			 break;
		 case 'p':
			 if ( !strcmp(key, "positive") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_empk_set_positive);
				 break;
			 }
			 break;
		 case 's':
			 if ( !strcmp(key, "smoothing") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_empk_set_smoothing);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "varcor") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_empk_set_varcor);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_GIBBS:
		 switch (*key) {
		 case 'b':
			 if ( !strcmp(key, "burnin") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_gibbs_set_burnin);
				 break;
			 }
			 break;
		 case 'c':
			 if ( !strcmp(key, "c") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_gibbs_set_c);
				 break;
			 }
			 break;
		 case 't':
			 if ( !strcmp(key, "thinning") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_gibbs_set_thinning);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "variant_coordinate") ) {
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_gibbs_set_variant_coordinate);
				 break;
			 }
			 if ( !strcmp(key, "variant_random_direction") ) {
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_gibbs_set_variant_random_direction);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_HINV:
		 switch (*key) {
		 case 'b':
			 if ( !strcmp(key, "boundary") ) {
				 result = _unur_str_par_set_dd(par,key,type_args,args,unur_hinv_set_boundary);
				 break;
			 }
			 break;
		 case 'c':
			 if ( !strcmp(key, "cpoints") ) {
				 result = _unur_str_par_set_Di(par,key,type_args,args,unur_hinv_set_cpoints,mlist);
				 break;
			 }
			 break;
		 case 'g':
			 if ( !strcmp(key, "guidefactor") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_hinv_set_guidefactor);
				 break;
			 }
			 break;
		 case 'm':
			 if ( !strcmp(key, "max_intervals") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_hinv_set_max_intervals);
				 break;
			 }
			 break;
		 case 'o':
			 if ( !strcmp(key, "order") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_hinv_set_order);
				 break;
			 }
			 break;
		 case 'u':
			 if ( !strcmp(key, "u_resolution") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_hinv_set_u_resolution);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_HITRO:
		 switch (*key) {
		 case 'a':
			 if ( !strcmp(key, "adaptive_multiplier") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_hitro_set_adaptive_multiplier);
				 break;
			 }
			 break;
		 case 'b':
			 if ( !strcmp(key, "burnin") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_hitro_set_burnin);
				 break;
			 }
			 break;
		 case 'r':
			 if ( !strcmp(key, "r") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_hitro_set_r);
				 break;
			 }
			 break;
		 case 't':
			 if ( !strcmp(key, "thinning") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_hitro_set_thinning);
				 break;
			 }
			 break;
		 case 'u':
			 if ( !strcmp(key, "use_adaptiveline") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_hitro_set_use_adaptiveline);
				 break;
			 }
			 if ( !strcmp(key, "use_adaptiverectangle") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_hitro_set_use_adaptiverectangle);
				 break;
			 }
			 if ( !strcmp(key, "use_boundingrectangle") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_hitro_set_use_boundingrectangle);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "v") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_hitro_set_v);
				 break;
			 }
			 if ( !strcmp(key, "variant_coordinate") ) {
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_hitro_set_variant_coordinate);
				 break;
			 }
			 if ( !strcmp(key, "variant_random_direction") ) {
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_hitro_set_variant_random_direction);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_HRB:
		 switch (*key) {
		 case 'u':
			 if ( !strcmp(key, "upperbound") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_hrb_set_upperbound);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_hrb_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_HRD:
		 switch (*key) {
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_hrd_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_HRI:
		 switch (*key) {
		 case 'p':
			 if ( !strcmp(key, "p0") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_hri_set_p0);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_hri_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_ITDR:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "cp") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_itdr_set_cp);
				 break;
			 }
			 if ( !strcmp(key, "ct") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_itdr_set_ct);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_itdr_set_verify);
				 break;
			 }
			 break;
		 case 'x':
			 if ( !strcmp(key, "xi") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_itdr_set_xi);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_MVTDR:
		 switch (*key) {
		 case 'b':
			 if ( !strcmp(key, "boundsplitting") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_mvtdr_set_boundsplitting);
				 break;
			 }
			 break;
		 case 'm':
			 if ( !strcmp(key, "maxcones") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_mvtdr_set_maxcones);
				 break;
			 }
			 break;
		 case 's':
			 if ( !strcmp(key, "stepsmin") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_mvtdr_set_stepsmin);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_mvtdr_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_NINV:
		 switch (*key) {
		 case 'm':
			 if ( !strcmp(key, "max_iter") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_ninv_set_max_iter);
				 break;
			 }
			 break;
		 case 's':
			 if ( !strcmp(key, "start") ) {
				 result = _unur_str_par_set_dd(par,key,type_args,args,unur_ninv_set_start);
				 break;
			 }
			 break;
		 case 't':
			 if ( !strcmp(key, "table") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_ninv_set_table);
				 break;
			 }
			 break;
		 case 'u':
			 if ( !strcmp(key, "u_resolution") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_ninv_set_u_resolution);
				 break;
			 }
			 if ( !strcmp(key, "usebisect") ) {
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_ninv_set_usebisect);
				 break;
			 }
			 if ( !strcmp(key, "usenewton") ) {
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_ninv_set_usenewton);
				 break;
			 }
			 if ( !strcmp(key, "useregula") ) {
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_ninv_set_useregula);
				 break;
			 }
			 break;
		 case 'x':
			 if ( !strcmp(key, "x_resolution") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_ninv_set_x_resolution);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_NROU:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "center") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_nrou_set_center);
				 break;
			 }
			 break;
		 case 'r':
			 if ( !strcmp(key, "r") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_nrou_set_r);
				 break;
			 }
			 break;
		 case 'u':
			 if ( !strcmp(key, "u") ) {
				 result = _unur_str_par_set_dd(par,key,type_args,args,unur_nrou_set_u);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "v") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_nrou_set_v);
				 break;
			 }
			 if ( !strcmp(key, "verify") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_nrou_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_PINV:
		 switch (*key) {
		 case 'b':
			 if ( !strcmp(key, "boundary") ) {
				 result = _unur_str_par_set_dd(par,key,type_args,args,unur_pinv_set_boundary);
				 break;
			 }
			 break;
		 case 'k':
			 if ( !strcmp(key, "keepcdf") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_pinv_set_keepcdf);
				 break;
			 }
			 break;
		 case 'm':
			 if ( !strcmp(key, "max_intervals") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_pinv_set_max_intervals);
				 break;
			 }
			 break;
		 case 'o':
			 if ( !strcmp(key, "order") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_pinv_set_order);
				 break;
			 }
			 break;
		 case 's':
			 if ( !strcmp(key, "searchboundary") ) {
				 result = _unur_str_par_set_ii(par,key,type_args,args,unur_pinv_set_searchboundary);
				 break;
			 }
			 if ( !strcmp(key, "smoothness") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_pinv_set_smoothness);
				 break;
			 }
			 break;
		 case 'u':
			 if ( !strcmp(key, "u_resolution") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_pinv_set_u_resolution);
				 break;
			 }
			 if ( !strcmp(key, "use_upoints") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_pinv_set_use_upoints);
				 break;
			 }
			 if ( !strcmp(key, "usecdf") ) {
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_pinv_set_usecdf);
				 break;
			 }
			 if ( !strcmp(key, "usepdf") ) {
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_pinv_set_usepdf);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_SROU:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "cdfatmode") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_srou_set_cdfatmode);
				 break;
			 }
			 break;
		 case 'p':
			 if ( !strcmp(key, "pdfatmode") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_srou_set_pdfatmode);
				 break;
			 }
			 break;
		 case 'r':
			 if ( !strcmp(key, "r") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_srou_set_r);
				 break;
			 }
			 break;
		 case 'u':
			 if ( !strcmp(key, "usemirror") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_srou_set_usemirror);
				 break;
			 }
			 if ( !strcmp(key, "usesqueeze") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_srou_set_usesqueeze);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_srou_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_SSR:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "cdfatmode") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_ssr_set_cdfatmode);
				 break;
			 }
			 break;
		 case 'p':
			 if ( !strcmp(key, "pdfatmode") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_ssr_set_pdfatmode);
				 break;
			 }
			 break;
		 case 'u':
			 if ( !strcmp(key, "usesqueeze") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_ssr_set_usesqueeze);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_ssr_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_TABL:
		 switch (*key) {
		 case 'a':
			 if ( !strcmp(key, "areafraction") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_tabl_set_areafraction);
				 break;
			 }
			 break;
		 case 'b':
			 if ( !strcmp(key, "boundary") ) {
				 result = _unur_str_par_set_dd(par,key,type_args,args,unur_tabl_set_boundary);
				 break;
			 }
			 break;
		 case 'c':
			 if ( !strcmp(key, "cpoints") ) {
				 result = _unur_str_par_set_iD(par,key,type_args,args,unur_tabl_set_cpoints,mlist);
				 break;
			 }
			 break;
		 case 'd':
			 if ( !strcmp(key, "darsfactor") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_tabl_set_darsfactor);
				 break;
			 }
			 break;
		 case 'g':
			 if ( !strcmp(key, "guidefactor") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_tabl_set_guidefactor);
				 break;
			 }
			 break;
		 case 'm':
			 if ( !strcmp(key, "max_intervals") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tabl_set_max_intervals);
				 break;
			 }
			 if ( !strcmp(key, "max_sqhratio") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_tabl_set_max_sqhratio);
				 break;
			 }
			 break;
		 case 'n':
			 if ( !strcmp(key, "nstp") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tabl_set_nstp);
				 break;
			 }
			 break;
		 case 'p':
			 if ( !strcmp(key, "pedantic") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tabl_set_pedantic);
				 break;
			 }
			 break;
		 case 's':
			 if ( !strcmp(key, "slopes") ) {
				 result = _unur_str_par_set_Di(par,key,type_args,args,unur_tabl_set_slopes,mlist);
				 break;
			 }
			 break;
		 case 'u':
			 if ( !strcmp(key, "usedars") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tabl_set_usedars);
				 break;
			 }
			 if ( !strcmp(key, "useear") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tabl_set_useear);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "variant_ia") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tabl_set_variant_ia);
				 break;
			 }
			 if ( !strcmp(key, "variant_splitmode") ) {
				 result = _unur_str_par_set_u(par,key,type_args,args,unur_tabl_set_variant_splitmode);
				 break;
			 }
			 if ( !strcmp(key, "verify") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tabl_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_TDR:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "c") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_tdr_set_c);
				 break;
			 }
			 if ( !strcmp(key, "cpoints") ) {
				 result = _unur_str_par_set_iD(par,key,type_args,args,unur_tdr_set_cpoints,mlist);
				 break;
			 }
			 break;
		 case 'd':
			 if ( !strcmp(key, "darsfactor") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_tdr_set_darsfactor);
				 break;
			 }
			 break;
		 case 'g':
			 if ( !strcmp(key, "guidefactor") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_tdr_set_guidefactor);
				 break;
			 }
			 break;
		 case 'm':
			 if ( !strcmp(key, "max_intervals") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tdr_set_max_intervals);
				 break;
			 }
			 if ( !strcmp(key, "max_sqhratio") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_tdr_set_max_sqhratio);
				 break;
			 }
			 break;
		 case 'p':
			 if ( !strcmp(key, "pedantic") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tdr_set_pedantic);
				 break;
			 }
			 break;
		 case 'r':
			 if ( !strcmp(key, "reinit_ncpoints") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tdr_set_reinit_ncpoints);
				 break;
			 }
			 if ( !strcmp(key, "reinit_percentiles") ) {
				 result = _unur_str_par_set_iD(par,key,type_args,args,unur_tdr_set_reinit_percentiles,mlist);
				 break;
			 }
			 break;
		 case 'u':
			 if ( !strcmp(key, "usecenter") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tdr_set_usecenter);
				 break;
			 }
			 if ( !strcmp(key, "usedars") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tdr_set_usedars);
				 break;
			 }
			 if ( !strcmp(key, "usemode") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tdr_set_usemode);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "variant_gw") ) {
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_tdr_set_variant_gw);
				 break;
			 }
			 if ( !strcmp(key, "variant_ia") ) {
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_tdr_set_variant_ia);
				 break;
			 }
			 if ( !strcmp(key, "variant_ps") ) {
				 result = _unur_str_par_set_void(par,key,type_args,args,unur_tdr_set_variant_ps);
				 break;
			 }
			 if ( !strcmp(key, "verify") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_tdr_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_UTDR:
		 switch (*key) {
		 case 'c':
			 if ( !strcmp(key, "cpfactor") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_utdr_set_cpfactor);
				 break;
			 }
			 break;
		 case 'd':
			 if ( !strcmp(key, "deltafactor") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_utdr_set_deltafactor);
				 break;
			 }
			 break;
		 case 'p':
			 if ( !strcmp(key, "pdfatmode") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_utdr_set_pdfatmode);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "verify") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_utdr_set_verify);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_VEMPK:
		 switch (*key) {
		 case 's':
			 if ( !strcmp(key, "smoothing") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_vempk_set_smoothing);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "varcor") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_vempk_set_varcor);
				 break;
			 }
		 }
		 break;
	 case UNUR_METH_VNROU:
		 switch (*key) {
		 case 'r':
			 if ( !strcmp(key, "r") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_vnrou_set_r);
				 break;
			 }
			 break;
		 case 'v':
			 if ( !strcmp(key, "v") ) {
				 result = _unur_str_par_set_d(par,key,type_args,args,unur_vnrou_set_v);
				 break;
			 }
			 if ( !strcmp(key, "verify") ) {
				 result = _unur_str_par_set_i(par,key,type_args,args,unur_vnrou_set_verify);
				 break;
			 }
		 }
		 break;
	 }
   if (result == UNUR_ERR_STR_UNKNOWN) {
     if ( !strcmp(key, "debug") ) {
       result = _unur_str_par_set_u(par,key,type_args,args,unur_set_debug);
     }
   }
  if (result == UNUR_ERR_STR_UNKNOWN) {
    _unur_error_unknown(key,"parameter for given method");
    return UNUR_ERR_STR_UNKNOWN;
  }
  else if (result != UNUR_SUCCESS) {
    _unur_error_invalid(key,"set call");
    return UNUR_ERR_STR_SYNTAX;
  }
  return UNUR_SUCCESS;
} 
