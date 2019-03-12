/* Copyright (c) 2000-2019 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifdef UNUR_ENABLE_INFO
void
_unur_mvtdr_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;
  int samplesize = 10000;
  double rc;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   dimension = %d\n",GEN->dim);
  _unur_string_append(info,"   functions = PDF dPDF\n");
  _unur_distr_cvec_info_domain(gen);
  if ( distr->set & UNUR_DISTR_SET_MODE ) {
    _unur_string_append(info,"   mode      = ");
    _unur_distr_info_vector( gen, DISTR.mode, GEN->dim);
  }
  _unur_string_append(info,"\n");
  _unur_string_append(info,"   center    = ");
  _unur_distr_info_vector( gen, GEN->center, GEN->dim);
  if ( !(distr->set & UNUR_DISTR_SET_CENTER) ) {
    if ( distr->set & UNUR_DISTR_SET_MODE )
      _unur_string_append(info,"  [= mode]");
    else
      _unur_string_append(info,"  [default]");
  }
  _unur_string_append(info,"\n\n");
  if (help) {
    if ( !(distr->set & UNUR_DISTR_SET_MODE) ) 
      _unur_string_append(info,"[ Hint: %s ]\n",
 			  "You can set the mode to improve the rejection constant.");
    _unur_string_append(info,"\n");
  } 
  _unur_string_append(info,"method: MVTDR (Multi-Variate Transformed Density Rejection)\n");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   volume(hat) = %g\n", GEN->Htot);
  _unur_string_append(info,"   rejection constant ");
  if (distr->set & UNUR_DISTR_SET_PDFVOLUME) {
    _unur_string_append(info,"= %g\n", GEN->Htot / DISTR.volume);
  }
  else {
    rc = unur_test_count_urn(gen,samplesize,0,NULL)/((1.+GEN->dim)*samplesize);
    _unur_string_append(info,"= %.2f  [approx.]\n", rc);
  }
  _unur_string_append(info,"   # cones = %d\n", GEN->n_cone);
  _unur_string_append(info,"   # vertices = %d\n", GEN->n_vertex);
  if (GEN->steps_min == GEN->n_steps)
    _unur_string_append(info,"   triangulation levels = %d\n", GEN->n_steps);
  else
    _unur_string_append(info,"   triangulation levels = %d-%d\n", GEN->steps_min, GEN->n_steps);
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters:\n");
    _unur_string_append(info,"   stepsmin = %d  %s\n", GEN->steps_min,
 			(gen->set & MVTDR_SET_STEPSMIN) ? "" : "[default]");
    _unur_string_append(info,"   maxcones = %d  %s\n", GEN->max_cones,
 			(gen->set & MVTDR_SET_MAXCONES) ? "" : "[default]");
    _unur_string_append(info,"   boundsplitting = %g  %s\n", GEN->bound_splitting,
 			(gen->set & MVTDR_SET_BOUNDSPLITTING) ? "" : "[default]");
    if (gen->variant & MVTDR_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");
    _unur_string_append(info,"\n");
  }
  if (help) {
    if ( !(gen->set & MVTDR_SET_STEPSMIN) )
      _unur_string_append(info,"[ Hint: %s ]\n",
 			  "You can increase \"stepsmin\" to improve the rejection constant." );
    if (GEN->max_cones <= GEN->n_cone)
      _unur_string_append(info,"[ Hint: %s ]\n",
 			  "You can increase \"maxcones\" to improve the rejection constant." );
    if ( !(gen->set & MVTDR_SET_BOUNDSPLITTING) )
      _unur_string_append(info,"[ Hint: %s ]\n",
 			  "You can change \"boundsplitting\" to change the creating of the hat function." );
    _unur_string_append(info,"\n");
  }
} 
#endif   
