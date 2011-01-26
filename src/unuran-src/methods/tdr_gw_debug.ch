/* Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold */
/* Department of Statistics and Mathematics, WU Wien, Austria  */

void
_unur_tdr_gw_debug_intervals( const struct unur_gen *gen, int print_areas )
{
  FILE *LOG;
  struct unur_tdr_interval *iv;
  double sAsqueeze, sAhatl, sAhatr, Atotal;
  int i;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:Intervals: %d\n",gen->genid,GEN->n_ivs);
  if (GEN->iv) {
    if (gen->debug & TDR_DEBUG_IV) {
      fprintf(LOG,"%s: Nr.            tp            ip          f(tp)      T(f(tp))    d(T(f(tp)))      squeeze\n",gen->genid);
      for (iv = GEN->iv, i=0; iv->next!=NULL; iv=iv->next, i++) {
	COOKIE_CHECK(iv,CK_TDR_IV,RETURN_VOID); 
	fprintf(LOG,"%s:[%3d]: %#12.6g  %#12.6g  %#12.6g  %#12.6g  %#12.6g  %#12.6g\n", gen->genid, i,
		iv->x, iv->ip, iv->fx, iv->Tfx, iv->dTfx, iv->sq);
      }
      COOKIE_CHECK(iv,CK_TDR_IV,RETURN_VOID); 
      fprintf(LOG,"%s:[...]: %#12.6g                %#12.6g  %#12.6g  %#12.6g\n", gen->genid,
	      iv->x, iv->fx, iv->Tfx, iv->dTfx);
    }
    fprintf(LOG,"%s:\n",gen->genid);
  }
  else
    fprintf(LOG,"%s: No intervals !\n",gen->genid);
  if (!print_areas || GEN->Atotal <= 0.) return;
  Atotal = GEN->Atotal;
  if (gen->debug & TDR_DEBUG_IV) {
    fprintf(LOG,"%s:Areas in intervals:\n",gen->genid);
    fprintf(LOG,"%s: Nr.\tbelow squeeze\t\t  below hat (left and right)\t\t  cumulated\n",gen->genid);
    sAsqueeze = sAhatl = sAhatr = 0.;
    if (GEN->iv) {
      for (iv = GEN->iv, i=0; iv->next!=NULL; iv=iv->next, i++) {
	COOKIE_CHECK(iv,CK_TDR_IV,RETURN_VOID); 
	sAsqueeze += iv->Asqueeze;
	sAhatl += iv->Ahat - iv->Ahatr;
	sAhatr += iv->Ahatr;
	fprintf(LOG,"%s:[%3d]: %-12.6g(%6.3f%%)  |  %-12.6g+ %-12.6g(%6.3f%%)  |  %-12.6g(%6.3f%%)\n",
		gen->genid,i,
		iv->Asqueeze, iv->Asqueeze * 100. / Atotal,
		iv->Ahat-iv->Ahatr, iv->Ahatr, iv->Ahat * 100. / Atotal, 
		iv->Acum, iv->Acum * 100. / Atotal);
      }
      fprintf(LOG,"%s:       ----------  ---------  |  ------------------------  ---------  +\n",gen->genid);
      fprintf(LOG,"%s: Sum : %-12.6g(%6.3f%%)            %-12.6g      (%6.3f%%)\n",gen->genid,
	      sAsqueeze, sAsqueeze * 100. / Atotal,
	      sAhatl+sAhatr, (sAhatl+sAhatr) * 100. / Atotal);
      fprintf(LOG,"%s:\n",gen->genid);
    }
  }
  fprintf(LOG,"%s: A(squeeze)     = %-12.6g  (%6.3f%%)\n",gen->genid,
	  GEN->Asqueeze, GEN->Asqueeze * 100./Atotal);
  fprintf(LOG,"%s: A(hat\\squeeze) = %-12.6g  (%6.3f%%)\n",gen->genid,
	  Atotal - GEN->Asqueeze, (Atotal - GEN->Asqueeze) * 100./Atotal);
  fprintf(LOG,"%s: A(total)       = %-12.6g\n",gen->genid, Atotal);
  fprintf(LOG,"%s:\n",gen->genid);
} 
void
_unur_tdr_gw_debug_sample( const struct unur_gen *gen, 
			   const struct unur_tdr_interval *iv, 
			   const struct unur_tdr_interval *pt, 
			   double x, double fx, double hx, double sqx )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);
  CHECK_NULL(iv,RETURN_VOID);   COOKIE_CHECK(iv,CK_TDR_IV,RETURN_VOID);
  CHECK_NULL(pt,RETURN_VOID);   COOKIE_CHECK(pt,CK_TDR_IV,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s:\n",gen->genid);
  if (iv == pt)
    fprintf(LOG,"%s: point generated in left part:\n",gen->genid);
  else
    fprintf(LOG,"%s: point generated in right part:\n",gen->genid);
  fprintf(LOG,"%s: construction point: x0 = %g\n",gen->genid,pt->x);
  fprintf(LOG,"%s: transformed hat     Th(x) = %g + %g * (x - %g)\n",gen->genid,pt->Tfx,pt->dTfx,pt->x);
  fprintf(LOG,"%s: transformed squeeze Ts(x) = %g + %g * (x - %g)\n",gen->genid,iv->Tfx,iv->sq,iv->x);
  fprintf(LOG,"%s: generated point: x = %g\n",gen->genid,x);
  fprintf(LOG,"%s:  h(x) = %.20g\n",gen->genid,hx);
  fprintf(LOG,"%s:  f(x) = %.20g\n",gen->genid,fx);
  fprintf(LOG,"%s:  s(x) = %.20g\n",gen->genid,sqx);
  fprintf(LOG,"%s:    hat: x - x0 = %g",gen->genid,x-pt->x);
  if (x < pt->x && iv == pt) fprintf(LOG,"  <-- error\n");
  else       fprintf(LOG,"\n");
  fprintf(LOG,"%s:    h(x) - f(x) = %g",gen->genid,hx-fx);
  if (hx<fx) fprintf(LOG,"  <-- error\n");
  else       fprintf(LOG,"\n");
  fprintf(LOG,"%s:    squeeze: x - x0 = %g",gen->genid,x-iv->x);
  if (x > pt->x && iv != pt) fprintf(LOG,"  <-- error\n");
  else       fprintf(LOG,"\n");
  fprintf(LOG,"%s:    f(x) - s(x) = %g",gen->genid,fx-sqx);
  if (fx<sqx) fprintf(LOG,"  <-- error\n");
  else       fprintf(LOG,"\n");
  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);
} 
void
_unur_tdr_gw_debug_split_start( const struct unur_gen *gen, 
				const struct unur_tdr_interval *iv,
				double x, double fx )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);  COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);
  CHECK_NULL(iv,RETURN_VOID);   COOKIE_CHECK(iv,CK_TDR_IV,RETURN_VOID);
  LOG = unur_get_stream();
  fprintf(LOG,"%s: split interval at x = %g \t\tf(x) = %g\n",gen->genid,x,fx);
  fprintf(LOG,"%s: old interval:\n",gen->genid);
  fprintf(LOG,"%s:   left  construction point = %-12.6g\tf(x) = %-12.6g\n",gen->genid,iv->x,iv->fx);
  fprintf(LOG,"%s:   right construction point = %-12.6g\tf(x) = %-12.6g\n",gen->genid,iv->next->x,iv->next->fx);
  fprintf(LOG,"%s:   A(squeeze)     = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  iv->Asqueeze,iv->Asqueeze*100./GEN->Atotal);
  fprintf(LOG,"%s:   A(hat\\squeeze) = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  (iv->Ahat - iv->Asqueeze),(iv->Ahat - iv->Asqueeze)*100./GEN->Atotal);
  fprintf(LOG,"%s:   A(hat)         = %-12.6g +  %-12.6g(%6.3f%%)\n",gen->genid,
	  iv->Ahat - iv->Ahatr, iv->Ahatr, iv->Ahat*100./GEN->Atotal);
  fflush(LOG);
} 
void
_unur_tdr_gw_debug_split_stop( const struct unur_gen *gen, 
			       const struct unur_tdr_interval *iv_left, 
			       const struct unur_tdr_interval *iv_right )
{
  FILE *LOG;
  CHECK_NULL(gen,RETURN_VOID);       COOKIE_CHECK(gen,CK_TDR_GEN,RETURN_VOID);
  CHECK_NULL(iv_left,RETURN_VOID);   COOKIE_CHECK(iv_left,CK_TDR_IV,RETURN_VOID);
  if (iv_right == NULL) iv_right = iv_left;
  LOG = unur_get_stream();
  fprintf(LOG,"%s: inserted point:\n",gen->genid);
  fprintf(LOG,"%s: x = %g, f(x) = %g, Tf(x)=%g, dTf(x) = %g, squeeze = %g:\n",
	  gen->genid, iv_right->x, iv_right->fx, iv_right->Tfx, iv_right->dTfx, iv_right->sq);
  fprintf(LOG,"%s: new intervals:\n",gen->genid);
  fprintf(LOG,"%s:   left   construction point = %g\n",gen->genid, iv_left->x);
  if (iv_left != iv_right)
    fprintf(LOG,"%s:   middle construction point = %g\n",gen->genid, iv_right->x);
  fprintf(LOG,"%s:   right  construction point = %g\n",gen->genid, iv_right->next->x);
  fprintf(LOG,"%s: left interval:\n",gen->genid);
  fprintf(LOG,"%s:   A(squeeze)     = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  iv_left->Asqueeze,
	  iv_left->Asqueeze*100./GEN->Atotal);
  fprintf(LOG,"%s:   A(hat\\squeeze) = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  (iv_left->Ahat - iv_left->Asqueeze),
	  (iv_left->Ahat - iv_left->Asqueeze) * 100./GEN->Atotal);
  fprintf(LOG,"%s:   A(hat)         = %-12.6g +  %-12.6g(%6.3f%%)\n",gen->genid,
	  iv_left->Ahat - iv_left->Ahatr,
	  iv_left->Ahatr,
	  iv_left->Ahat * 100./GEN->Atotal);
  if (iv_left == iv_right)
    fprintf(LOG,"%s: interval chopped.\n",gen->genid);
  else {
    fprintf(LOG,"%s: right interval:\n",gen->genid);
    fprintf(LOG,"%s:   A(squeeze)     = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	    iv_right->Asqueeze,
	    iv_right->Asqueeze*100./GEN->Atotal);
    fprintf(LOG,"%s:   A(hat\\squeeze) = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	    (iv_right->Ahat - iv_right->Asqueeze),
	    (iv_right->Ahat - iv_right->Asqueeze) * 100./GEN->Atotal);
    fprintf(LOG,"%s:   A(hat)         = %-12.6g +  %-12.6g(%6.3f%%)\n",gen->genid,
	    iv_right->Ahat - iv_right->Ahatr,
	    iv_right->Ahatr,
	    iv_right->Ahat * 100./GEN->Atotal);
  }
  fprintf(LOG,"%s: total areas:\n",gen->genid);
  fprintf(LOG,"%s:   A(squeeze)     = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  GEN->Asqueeze, GEN->Asqueeze * 100./GEN->Atotal);
  fprintf(LOG,"%s:   A(hat\\squeeze) = %-12.6g\t\t(%6.3f%%)\n",gen->genid,
	  GEN->Atotal - GEN->Asqueeze, (GEN->Atotal - GEN->Asqueeze) * 100./GEN->Atotal);
  fprintf(LOG,"%s:   A(total)       = %-12.6g\n",gen->genid, GEN->Atotal);
  fprintf(LOG,"%s:\n",gen->genid);
  fflush(LOG);
} 
