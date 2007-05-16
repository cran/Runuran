/* Copyright (c) 2000-2007 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

void
_unur_tdr_ps_debug_intervals( const struct unur_gen *gen, int print_areas )
{
  FILE *log;
  struct unur_tdr_interval *iv;
  double sAsqueeze, sAhatl, sAhatr, Atotal;
  int i;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);
  log = unur_get_stream();
  fprintf(log,"%s:Intervals: %d\n",gen->genid,GEN->n_ivs);
  if (GEN->iv) {
    if (gen->debug & TDR_DEBUG_IV) {
      fprintf(log,"%s: Nr.       left ip           tp        f(tp)     T(f(tp))   d(T(f(tp)))       f(ip)   squ. ratio\n",gen->genid);
      for (iv=GEN->iv,i=0; iv->next; iv=iv->next, i++) {
	COOKIE_CHECK(iv,CK_TDR_IV,RETURN_VOID); 
	fprintf(log,"%s:[%3d]:%#12.6g %#12.6g %#12.6g %#12.6g %#12.6g %#12.6g %#12.6g\n", gen->genid, i,
		iv->ip, iv->x, iv->fx, iv->Tfx, iv->dTfx, iv->fip, iv->sq);
      }
      COOKIE_CHECK(iv,CK_TDR_IV,RETURN_VOID); 
      fprintf(log,"%s:[...]:%#12.6g\t\t\t\t\t\t       %#12.6g\n", gen->genid,
	      iv->ip, iv->fip);
    }
    fprintf(log,"%s:\n",gen->genid);
  }
  else
    fprintf(log,"%s: No intervals !\n",gen->genid);
  if (!print_areas || GEN->Atotal <= 0.) return;
  Atotal = GEN->Atotal;
  if (gen->debug & TDR_DEBUG_IV) {
    fprintf(log,"%s:Areas in intervals:\n",gen->genid);
    fprintf(log,"%s: Nr.\tbelow squeeze\t\t  below hat (left and right)\t\t  cumulated\n",gen->genid);
    sAsqueeze = sAhatl = sAhatr = 0.;
    if (GEN->iv) {
      for (iv=GEN->iv,i=0; iv->next; iv=iv->next, i++) {
	COOKIE_CHECK(iv,CK_TDR_IV,RETURN_VOID); 
	sAsqueeze += iv->Asqueeze;
	sAhatl += iv->Ahat - iv->Ahatr;
	sAhatr += iv->Ahatr;
	fprintf(log,"%s:[%3d]: %-12.6g(%6.3f%%)  |  %-12.6g+ %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)\n",
		gen->genid,i,
		iv->Asqueeze, iv->Asqueeze * 100. / Atotal,
		iv->Ahat-iv->Ahatr, iv->Ahatr, iv->Ahat * 100. / Atotal, 
		iv->Acum, iv->Acum * 100. / Atotal);
      }
      fprintf(log,"%s:       ----------  ---------  |  ------------------------  ---------  +\n",gen->genid);
      fprintf(log,"%s: Sum : %-12.6g(%6.3f%%)            %-12.6g      (%6.3f%%)\n",gen->genid,
	      sAsqueeze, sAsqueeze * 100. / Atotal,
	      sAhatl+sAhatr, (sAhatl+sAhatr) * 100. / Atotal);
      fprintf(log,"%s:\n",gen->genid);
    }
  }
  fprintf(log,"%s: A(squeeze)     = %-12.6g  (%6.3f%%)\n",gen->genid,
	  GEN->Asqueeze, GEN->Asqueeze * 100./Atotal);
  fprintf(log,"%s: A(hat\\squeeze) = %-12.6g  (%6.3f%%)\n",gen->genid,
	  Atotal - GEN->Asqueeze, (Atotal - GEN->Asqueeze) * 100./Atotal);
  fprintf(log,"%s: A(total)       = %-12.6g\n",gen->genid, Atotal);
  fprintf(log,"%s:\n",gen->genid);
} 
void
_unur_tdr_ps_debug_sample( const struct unur_gen *gen, 
			   const struct unur_tdr_interval *iv, 
			   double x, double fx, double hx, double sqx )
{
  FILE *log;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);
  CHECK_NULL(iv,RETURN_VOID);   COOKIE_CHECK(iv,CK_TDR_IV,RETURN_VOID);
  log = unur_get_stream();
  fprintf(log,"%s:\n",gen->genid);
  fprintf(log,"%s: construction point: x0 = %g\n",gen->genid,iv->x);
  fprintf(log,"%s: transformed hat Th(x) = %g + %g * (x - %g)\n",gen->genid,iv->Tfx,iv->dTfx,iv->x);
  fprintf(log,"%s: squeeze ratio = %g\n",gen->genid,iv->sq);
  fprintf(log,"%s: generated point: x = %g\n",gen->genid,x);
  fprintf(log,"%s:  h(x) = %.20g\n",gen->genid,hx);
  fprintf(log,"%s:  f(x) = %.20g\n",gen->genid,fx);
  fprintf(log,"%s:  s(x) = %.20g\n",gen->genid,sqx);
  fprintf(log,"%s:    h(x) - f(x) = %g",gen->genid,hx-fx);
  if (hx<fx) fprintf(log,"  <-- error\n");
  else       fprintf(log,"\n");
  fprintf(log,"%s:    f(x) - s(x) = %g",gen->genid,fx-sqx);
  if (fx<sqx) fprintf(log,"  <-- error\n");
  else       fprintf(log,"\n");
  fprintf(log,"%s:\n",gen->genid);
  fflush(log);
} 
void
_unur_tdr_ps_debug_split_start( const struct unur_gen *gen, 
				const struct unur_tdr_interval *iv_left, 
				const struct unur_tdr_interval *iv_right,
				double x, double fx )
{
  FILE *log;
  CHECK_NULL(gen,RETURN_VOID);      COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);
  log = unur_get_stream();
  fprintf(log,"%s: split interval at x = %g \t\tf(x) = %g\n",gen->genid,x,fx);
  fprintf(log,"%s: old intervals:\n",gen->genid);
  if (iv_left) {
    fprintf(log,"%s:   left boundary point      = %-12.6g\tf(x) = %-12.6g\n",gen->genid,iv_left->ip,iv_left->fip);
    fprintf(log,"%s:   left construction point  = %-12.6g\tf(x) = %-12.6g\n",gen->genid,iv_left->x,iv_left->fx);
  }
  fprintf(log,"%s:   middle boundary point    = %-12.6g\tf(x) = %-12.6g\n",gen->genid,iv_right->ip,iv_right->fip);
  if (iv_right->next) {
    fprintf(log,"%s:   right construction point = %-12.6g\tf(x) = %-12.6g\n",gen->genid,iv_right->x,iv_right->fx);
    fprintf(log,"%s:   right boundary point     = %-12.6g\tf(x) = %-12.6g\n",gen->genid,
	    iv_right->next->ip,iv_right->next->fip);
  }
  fprintf(log,"%s:   A(squeeze) =\n",gen->genid);
  if (iv_left)
    fprintf(log,"%s:\t%-12.6g\t(%6.3f%%)\n",gen->genid,
	    iv_left->Asqueeze,iv_left->Asqueeze*100./GEN->Atotal);
  if (iv_right->next)
    fprintf(log,"%s:\t%-12.6g\t(%6.3f%%)\n",gen->genid,
	    iv_right->Asqueeze,iv_right->Asqueeze*100./GEN->Atotal);
  fprintf(log,"%s:   A(hat\\squeeze) =\n",gen->genid);
  if (iv_left)
    fprintf(log,"%s:\t%-12.6g\t(%6.3f%%)\n",gen->genid,
	    (iv_left->Ahat - iv_left->Asqueeze),(iv_left->Ahat - iv_left->Asqueeze)*100./GEN->Atotal);
  if (iv_right->next)
    fprintf(log,"%s:\t%-12.6g\t(%6.3f%%)\n",gen->genid,
	  (iv_right->Ahat - iv_right->Asqueeze),(iv_right->Ahat - iv_right->Asqueeze)*100./GEN->Atotal);
  fprintf(log,"%s:   A(hat) =\n",gen->genid);
  if (iv_left)
    fprintf(log,"%s:\t%-12.6g\t(%6.3f%%)\n",gen->genid,
	    iv_left->Ahat, iv_left->Ahat*100./GEN->Atotal);
  if (iv_right->next)
    fprintf(log,"%s:\t%-12.6g\t(%6.3f%%)\n",gen->genid,
	    iv_right->Ahat, iv_right->Ahat*100./GEN->Atotal);
  fflush(log);
} 
void
_unur_tdr_ps_debug_split_stop( const struct unur_gen *gen, 
			       const struct unur_tdr_interval *iv_left, 
			       const struct unur_tdr_interval *iv_middle, 
			       const struct unur_tdr_interval *iv_right )
{
  FILE *log;
  CHECK_NULL(gen,RETURN_VOID);       COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);
  log = unur_get_stream();
  fprintf(log,"%s: new intervals:\n",gen->genid);
  if (iv_left) {
    fprintf(log,"%s:   left boundary point      = %-12.6g\tf(x) = %-12.6g\n",gen->genid,iv_left->ip,iv_left->fip);
    fprintf(log,"%s:   left construction point  = %-12.6g\tf(x) = %-12.6g\n",gen->genid,iv_left->x,iv_left->fx);
  }
  if (iv_middle) {
    fprintf(log,"%s:   middle boundary point    = %-12.6g\tf(x) = %-12.6g\n",gen->genid,iv_middle->ip,iv_middle->fip);
    fprintf(log,"%s:   middle construction point= %-12.6g\tf(x) = %-12.6g\n",gen->genid,iv_middle->x,iv_middle->fx);
  }
  fprintf(log,"%s:   middle boundary point    = %-12.6g\tf(x) = %-12.6g\n",gen->genid,iv_right->ip,iv_right->fip);
  if (iv_right->next) {
    fprintf(log,"%s:   right construction point = %-12.6g\tf(x) = %-12.6g\n",gen->genid,iv_right->x,iv_right->fx);
    fprintf(log,"%s:   right boundary point     = %-12.6g\tf(x) = %-12.6g\n",gen->genid,
	    iv_right->next->ip,iv_right->next->fip);
  }
  fprintf(log,"%s:   A(squeeze) =\n",gen->genid);
  if (iv_left)
    fprintf(log,"%s:\t%-12.6g\t(%6.3f%%)\n",gen->genid,
	    iv_left->Asqueeze,iv_left->Asqueeze*100./GEN->Atotal);
  if (iv_middle)
    fprintf(log,"%s:\t%-12.6g\t(%6.3f%%)\n",gen->genid,
	    iv_middle->Asqueeze,iv_middle->Asqueeze*100./GEN->Atotal);
  if (iv_right->next)
    fprintf(log,"%s:\t%-12.6g\t(%6.3f%%)\n",gen->genid,
	    iv_right->Asqueeze,iv_right->Asqueeze*100./GEN->Atotal);
  fprintf(log,"%s:   A(hat\\squeeze) =\n",gen->genid);
  if (iv_left)
    fprintf(log,"%s:\t%-12.6g\t(%6.3f%%)\n",gen->genid,
	    (iv_left->Ahat - iv_left->Asqueeze),(iv_left->Ahat - iv_left->Asqueeze)*100./GEN->Atotal);
  if (iv_middle)
    fprintf(log,"%s:\t%-12.6g\t(%6.3f%%)\n",gen->genid,
	    (iv_middle->Ahat - iv_middle->Asqueeze),(iv_middle->Ahat - iv_middle->Asqueeze)*100./GEN->Atotal);
  if (iv_right->next)
    fprintf(log,"%s:\t%-12.6g\t(%6.3f%%)\n",gen->genid,
	  (iv_right->Ahat - iv_right->Asqueeze),(iv_right->Ahat - iv_right->Asqueeze)*100./GEN->Atotal);
  fprintf(log,"%s:   A(hat) =\n",gen->genid);
  if (iv_left)
    fprintf(log,"%s:\t%-12.6g\t(%6.3f%%)\n",gen->genid,
	    iv_left->Ahat, iv_left->Ahat*100./GEN->Atotal);
  if (iv_middle)
    fprintf(log,"%s:\t%-12.6g\t(%6.3f%%)\n",gen->genid,
	    iv_middle->Ahat, iv_middle->Ahat*100./GEN->Atotal);
  if (iv_right->next)
    fprintf(log,"%s:\t%-12.6g\t(%6.3f%%)\n",gen->genid,
	    iv_right->Ahat, iv_right->Ahat*100./GEN->Atotal);
  fprintf(log,"%s: total areas:\n",gen->genid);
  fprintf(log,"%s:   A(squeeze)     = %-12.6g   (%6.3f%%)\n",gen->genid,
	  GEN->Asqueeze, GEN->Asqueeze * 100./GEN->Atotal);
  fprintf(log,"%s:   A(hat\\squeeze) = %-12.6g   (%6.3f%%)\n",gen->genid,
	  GEN->Atotal - GEN->Asqueeze, (GEN->Atotal - GEN->Asqueeze) * 100./GEN->Atotal);
  fprintf(log,"%s:   A(total)       = %-12.6g\n",gen->genid, GEN->Atotal);
  fprintf(log,"%s:\n",gen->genid);
  fflush(log);
} 
