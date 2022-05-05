/* Copyright (c) 2000-2022 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifdef UNUR_ENABLE_INFO
void
_unur_tdr_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = PDF dPDF\n");
  _unur_string_append(info,"   domain    = (%g, %g)", DISTR.trunc[0],DISTR.trunc[1]);
  if (gen->distr->set & UNUR_DISTR_SET_TRUNCATED) {
    _unur_string_append(info,"   [truncated from (%g, %g)]", DISTR.domain[0],DISTR.domain[1]);
  }
  _unur_string_append(info,"\n");
  _unur_string_append(info,"   center    = %g", unur_distr_cont_get_center(distr));
  if ( !(distr->set & UNUR_DISTR_SET_CENTER) ) {
    if ( distr->set & UNUR_DISTR_SET_MODE )
      _unur_string_append(info,"  [= mode]\n");
    else 
      _unur_string_append(info,"  [default]\n");
  }
  else {
    _unur_string_append(info,"\n");
  }
  if (help) {
    if ( !(distr->set & (UNUR_DISTR_SET_CENTER | UNUR_DISTR_SET_MODE )) ) 
      _unur_string_append(info,"\n[ Hint: %s ]\n",
			  "You may provide a point near the mode as \"center\"."); 
  }
  _unur_string_append(info,"\n");
  _unur_string_append(info,"method: TDR (Transformed Density Rejection)\n");
  _unur_string_append(info,"   variant   = ");
  switch (gen->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:
    _unur_string_append(info,"GW (original Gilks & Wild)\n"); break;
  case TDR_VARIANT_PS:
    _unur_string_append(info,"PS (proportional squeeze)\n"); break;
  case TDR_VARIANT_IA:
    _unur_string_append(info,"IA (immediate acceptance)\n"); break;
  }
  _unur_string_append(info,"   T_c(x)    = ");
  switch( gen->variant & TDR_VARMASK_T ) {
  case TDR_VAR_T_LOG:
    _unur_string_append(info,"log(x)  ... c = 0\n"); break;
  case TDR_VAR_T_SQRT:
    _unur_string_append(info,"-1/sqrt(x)  ... c = -1/2\n"); break;
  case TDR_VAR_T_POW:
    _unur_string_append(info,"-x^(%g)  ... c = %g\n",GEN->c_T,GEN->c_T); break;
  }
  _unur_string_append(info,"\n");
  _unur_string_append(info,"performance characteristics:\n");
  _unur_string_append(info,"   area(hat) = %g\n", GEN->Atotal);
  _unur_string_append(info,"   rejection constant ");
  if (distr->set & UNUR_DISTR_SET_PDFAREA)
    _unur_string_append(info,"= %g\n", GEN->Atotal/DISTR.area);
  else
    _unur_string_append(info,"<= %g\n", GEN->Atotal/GEN->Asqueeze);
  _unur_string_append(info,"   area ratio squeeze/hat = %g\n",
		      GEN->Asqueeze/GEN->Atotal);
  _unur_string_append(info,"   # intervals = %d\n", GEN->n_ivs);
  _unur_string_append(info,"\n");
  if (help) {
    _unur_string_append(info,"parameters:\n");
    switch (gen->variant & TDR_VARMASK_VARIANT) {
    case TDR_VARIANT_GW:
      _unur_string_append(info,"   variant_gw = on\n"); break;
    case TDR_VARIANT_PS:
      _unur_string_append(info,"   variant_ps = on  [default]\n"); break;
    case TDR_VARIANT_IA:
      _unur_string_append(info,"   variant_ia = on\n"); break;
    }
    _unur_string_append(info,"   c = %g  %s\n", GEN->c_T,
			(gen->set & TDR_SET_C) ? "" : "[default]");
    _unur_string_append(info,"   max_sqhratio = %g  %s\n", GEN->max_ratio,
			(gen->set & TDR_SET_MAX_SQHRATIO) ? "" : "[default]");
    _unur_string_append(info,"   max_intervals = %d  %s\n", GEN->max_ivs_info,
			(gen->set & TDR_SET_MAX_IVS) ? "" : "[default]");
    if (gen->variant & TDR_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");
    if (gen->variant & TDR_VARFLAG_PEDANTIC)
      _unur_string_append(info,"   pedantic = on\n");
    _unur_string_append(info,"\n");
  }
  if (help) {
    if ( (gen->variant & TDR_VARMASK_VARIANT) != TDR_VARIANT_IA) 
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You may use \"variant_ia\" for faster generation times."); 
    if ( !(gen->set & TDR_SET_MAX_SQHRATIO) )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can set \"max_sqhratio\" closer to 1 to decrease rejection constant." );
    if (GEN->Asqueeze/GEN->Atotal < GEN->max_ratio) 
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You should increase \"max_intervals\" to obtain the desired rejection constant." );
    _unur_string_append(info,"\n");
  }
} 
#endif
