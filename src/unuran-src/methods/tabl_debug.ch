/* Copyright (c) 2000-2008 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifdef UNUR_ENABLE_LOGGING
#define empty_line() fprintf(log,"%s:\n",gen->genid);
void
_unur_tabl_debug_init_start( const struct unur_par *par, const struct unur_gen *gen )
{
  FILE *log;
  int i;
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_TABL_PAR,RETURN_VOID);
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);
  log = unur_get_stream();
  empty_line();
  fprintf(log,"%s: type    = continuous univariate random variates\n",gen->genid);
  fprintf(log,"%s: method  = rejection from piecewise constant hat\n",gen->genid);
  fprintf(log,"%s: variant = ",gen->genid);
  if (gen->variant & TABL_VARIANT_IA) 
    fprintf(log,"immediate acceptance ... IA\n");
  else
    fprintf(log,"acceptance/rejection ... RH\n");
  empty_line();
  _unur_distr_cont_debug( gen->distr, gen->genid );
  fprintf(log,"%s: sampling routine = _unur_tabl_",gen->genid);
  if (gen->variant & TABL_VARIANT_IA) 
    fprintf(log,"ia"); 
  else
    fprintf(log,"rh");
  if (gen->variant & TABL_VARFLAG_VERIFY)
    fprintf(log,"_sample_check()\n");
  else
    fprintf(log,"_sample()\n");
  fprintf(log,"%s: computation interval = (%g, %g)\n",gen->genid,GEN->bleft,GEN->bright);
  fprintf(log,"%s: maximum number of intervals = %d",gen->genid,GEN->max_ivs);
  _unur_print_if_default(gen,TABL_SET_MAX_IVS);
  fprintf(log,"\n");
  if (gen->variant & TABL_VARFLAG_USEEAR) {
    fprintf(log,"%s: use equal area rule\n",gen->genid);
    fprintf(log,"%s:\tarea fraction for equal area rule = %g ",gen->genid,PAR->area_fract);
    _unur_print_if_default(par,TABL_SET_AREAFRACTION);
    fprintf(log,"\n");
  }
  if (gen->variant & TABL_VARFLAG_USEDARS) {
    fprintf(log,"%s: Derandomized ARS enabled ",gen->genid);
    _unur_print_if_default(gen,TABL_SET_USE_DARS);
    fprintf(log,"\n%s:\tDARS factor = %g",gen->genid,GEN->darsfactor);
    _unur_print_if_default(gen,TABL_SET_DARS_FACTOR);
  }
  else {
    fprintf(log,"%s: Derandomized ARS disabled ",gen->genid);
    _unur_print_if_default(gen,TABL_SET_USE_DARS);
  }
  fprintf(log,"\n");
  fprintf(log,"%s: bound for ratio  Atotal / Asqueeze = %g%%",gen->genid,GEN->max_ratio*100.);
  _unur_print_if_default(gen,TABL_SET_MAX_SQHRATIO);
  fprintf(log,"\n");
  fprintf(log,"%s: split intervals at ",gen->genid);
  switch( gen->variant & TABL_VARMASK_SPLIT ) {
  case TABL_VARFLAG_SPLIT_MEAN:
    fprintf(log,"mean point");
    break;
  case TABL_VARFLAG_SPLIT_ARC:
    fprintf(log,"\"arcmean\" point");
    break;
  case TABL_VARFLAG_SPLIT_POINT:
  default:
    fprintf(log,"sample point");
    break;
  }
  fprintf(log," when using adaptive sampling.\n");
  empty_line();
  fprintf(log,"%s: sampling from list of intervals: indexed search (guide table method)\n",gen->genid);
  fprintf(log,"%s:    relative guide table size = %g%%",gen->genid,100.*GEN->guide_factor);
  _unur_print_if_default(gen,TABL_SET_GUIDEFACTOR);
  fprintf(log,"\n");
  empty_line();
  if (par->set & TABL_SET_SLOPES) {
    fprintf(log,"%s: slopes = %d\n",gen->genid,PAR->n_slopes);
    for (i=0; i<PAR->n_slopes; i++) {
      if ( PAR->slopes[2*i] > PAR->slopes[2*i+1] )
	fprintf(log,"%s:   (+)  ",gen->genid);
      else
	fprintf(log,"%s:   (-)  ",gen->genid);
      fprintf(log,"< %#g, %#g >\n", PAR->slopes[2*i], PAR->slopes[2*i+1] );
    }
  }
  else
    fprintf(log,"%s: no slopes given. compute from domain and mode.\n",gen->genid);
  if (PAR->cpoints != NULL) {
    fprintf(log,"%s: cpoints for slopes (%d):\n",gen->genid,PAR->n_cpoints);
    fprintf(log,"%s:\t %g",gen->genid,(PAR->cpoints)[0]);
    for (i=1; i<PAR->n_cpoints; i++)
      fprintf(log,", %g",(PAR->cpoints)[i]);
    fprintf(log,"\n");
  }
  fprintf(log,"%s: number of starting intervals (approx.) = %d",gen->genid,PAR->n_stp);
  _unur_print_if_default(par,TABL_SET_N_STP);
  fprintf(log,"\n");
  empty_line();
} 
void
_unur_tabl_debug_init_finished( const struct unur_gen *gen )
{
  FILE *log;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);
  log = unur_get_stream();
  _unur_tabl_debug_intervals(gen,"INIT completed",TRUE);
  fprintf(log,"%s: INIT completed **********************\n",gen->genid);
  empty_line();
} 
void 
_unur_tabl_debug_dars_start( const struct unur_par *par, const struct unur_gen *gen )
{
  FILE *log;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);
  CHECK_NULL(par,RETURN_VOID);  COOKIE_CHECK(par,CK_TABL_PAR,RETURN_VOID);
  log = unur_get_stream();
  fprintf(log,"%s: DARS started **********************\n",gen->genid);
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: DARS factor = %g",gen->genid,GEN->darsfactor);
  _unur_print_if_default(par,TABL_SET_DARS_FACTOR);
  fprintf(log,"\n%s:\n",gen->genid);
  fflush(log);
} 
void
_unur_tabl_debug_free( const struct unur_gen *gen )
{
  FILE *log;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);
  log = unur_get_stream();
  empty_line();
  if (gen->status == UNUR_SUCCESS) {
    fprintf(log,"%s: GENERATOR destroyed **********************\n",gen->genid);
    empty_line();
    _unur_tabl_debug_intervals(gen,NULL,TRUE);
  }
  else {
    fprintf(log,"%s: initialization of GENERATOR failed **********************\n",gen->genid);
  }
  empty_line();
  fflush(log);
} 
void
_unur_tabl_debug_intervals( const struct unur_gen *gen, const char *header, int print_areas )
{
  FILE *log;
  struct unur_tabl_interval *iv;
  double sAsqueeze, Atotal;
  int i;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TABL_GEN,RETURN_VOID);
  log = unur_get_stream();
  if (header)
    fprintf(log,"%s: %s\n",gen->genid,header);
  fprintf(log,"%s: intervals = %d\n",gen->genid,GEN->n_ivs);
  if (gen->debug & TABL_DEBUG_IV) {
    fprintf(log,"%s:             <   max       ,   min       >        f(max)          f(min) \n",gen->genid);
    fprintf(log,"%s:\n",gen->genid);
    for (iv = GEN->iv, i=0; iv!=NULL; iv=iv->next, i++) {
      COOKIE_CHECK(iv,CK_TABL_IV,RETURN_VOID);
      fprintf(log,"%s:[%3d]: (%s)   < %#-12.6g, %#-12.6g>   |  %#-12.6g    %#-12.6g  \n",
	      gen->genid, i, (iv->xmin <= iv->xmax)?"+":"-",
	      iv->xmax, iv->xmin, iv->fmax, iv->fmin);
    }
    empty_line();
  }
  if (!print_areas) return;
  if (GEN->Atotal <= 0.) {
    fprintf(log,"%s: Construction of hat function not successful\n",gen->genid);
    fprintf(log,"%s: Areas may be meaningless !!!!!!!!!!!!!!!!!!\n",gen->genid);
    empty_line();
    Atotal = -1.;   
  }
  else {
    Atotal = GEN->Atotal;
  }
  if (gen->debug & TABL_DEBUG_IV) {
    fprintf(log,"%s:Areas in intervals:\n",gen->genid);
    fprintf(log,"%s: Nr.\t below squeeze\t\t   below hat\t\t     cumulated\n",gen->genid);
    empty_line();
    sAsqueeze = 0.;
    for (iv = GEN->iv, i=0; iv!=NULL; iv=iv->next, i++) {
      COOKIE_CHECK(iv,CK_TABL_IV,RETURN_VOID); 
      sAsqueeze += iv->Asqueeze;
      fprintf(log,"%s:[%3d]: %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)\n",
	      gen->genid,i,
	      iv->Asqueeze, iv->Asqueeze * 100. / Atotal,
	      iv->Ahat, iv->Ahat * 100. / Atotal, 
	      iv->Acum, iv->Acum * 100. / Atotal);
    }
    fprintf(log,"%s:       ----------  ---------  +  ----------  ---------  +\n",gen->genid);
    fprintf(log,"%s: Sum : %-12.6g(%6.3f%%)     %-12.6g(100%%)\n",gen->genid,
	    sAsqueeze, sAsqueeze * 100. / Atotal, Atotal);
    empty_line();
  }
  fprintf(log,"%s: A(squeeze)     = %-12.6g  (%6.3f%%)\n",gen->genid,
	  GEN->Asqueeze, GEN->Asqueeze * 100./Atotal);
  fprintf(log,"%s: A(hat\\squeeze) = %-12.6g  (%6.3f%%)\n",gen->genid,
	  Atotal - GEN->Asqueeze, (Atotal - GEN->Asqueeze) * 100./Atotal);
  fprintf(log,"%s: A(total)       = %-12.6g\n",gen->genid, Atotal);
  empty_line();
} 
#endif   
