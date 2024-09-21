/* Copyright (c) 2000-2024 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

#ifdef UNUR_ENABLE_INFO
void
_unur_tabl_info( struct unur_gen *gen, int help )
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;
  _unur_string_append(info,"generator ID: %s\n\n", gen->genid);
  _unur_string_append(info,"distribution:\n");
  _unur_distr_info_typename(gen);
  _unur_string_append(info,"   functions = PDF\n");
  _unur_string_append(info,"   domain    = (%g, %g)", DISTR.trunc[0],DISTR.trunc[1]);
  if (gen->distr->set & UNUR_DISTR_SET_TRUNCATED) {
    _unur_string_append(info,"   [truncated from (%g, %g)]", DISTR.domain[0],DISTR.domain[1]);
  }
  _unur_string_append(info,"\n");
  _unur_string_append(info,"   mode      = %g   %s\n", unur_distr_cont_get_mode(distr),
		      (distr->set & UNUR_DISTR_SET_MODE_APPROX) ? "[numeric.]" : "");
  _unur_string_append(info,"   area(PDF) = ");
  if (gen->distr->set & UNUR_DISTR_SET_PDFAREA)
    _unur_string_append(info,"%g\n", DISTR.area);
  else 
    _unur_string_append(info,"[not set: use 1.0]\n");
  _unur_string_append(info,"\n");
  _unur_string_append(info,"method: TABL (Ahrens' TABLe Method)\n");
  _unur_string_append(info,"   variant   = ");
  if (gen->variant & TABL_VARIANT_IA) 
    _unur_string_append(info,"immediate acceptance [ia = on]\n");
  else
    _unur_string_append(info,"acceptance/rejection [ia = off]\n");
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
    if (gen->variant & TABL_VARIANT_IA) 
      _unur_string_append(info,"   variant_ia = on  [default]\n");
    else
      _unur_string_append(info,"   variant_ia = off\n");
    _unur_string_append(info,"   max_sqhratio = %g  %s\n", GEN->max_ratio,
			(gen->set & TABL_SET_MAX_SQHRATIO) ? "" : "[default]");
    _unur_string_append(info,"   max_intervals = %d  %s\n", GEN->max_ivs_info,
			(gen->set & TABL_SET_MAX_IVS) ? "" : "[default]");
    if (gen->variant & TABL_VARFLAG_VERIFY)
      _unur_string_append(info,"   verify = on\n");
    if (gen->variant & TABL_VARFLAG_PEDANTIC)
      _unur_string_append(info,"   pedantic = on\n");
    _unur_string_append(info,"\n");
  }
  if (help) {
    if ( !(gen->set & TABL_SET_MAX_SQHRATIO) )
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You can set \"max_sqhratio\" closer to 1 to decrease rejection constant." );
    if (GEN->Asqueeze/GEN->Atotal < GEN->max_ratio)
      _unur_string_append(info,"[ Hint: %s ]\n",
			  "You should increase \"max_intervals\" to obtain the desired rejection constant." );
    _unur_string_append(info,"\n");
  }
} 
#endif
