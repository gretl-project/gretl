/* print some window contents to a buffer suitable for
   text-to-speech
*/

#include "gretltypes.h"
#include "gretl_enums.h"

static void 
audioprint_model (MODEL *pmod, const DATAINFO *pdinfo, PRN *prn)
{
    if (pmod->ci != OLS) {
	pputs(prn, "Sorry, this model is not O.L.S.  I can't read it.\n");
    } else {
	pputs(prn, "O.L.S. model: not quite ready for speaking yet.\n");
    }
}

static void 
audioprint_summary (GRETLSUMMARY *summ, const DATAINFO *pdinfo,
		    PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN];
    int lo = summ->list[0], v, lv;

    ntodate(date1, pdinfo->t1, pdinfo);
    ntodate(date2, pdinfo->t2, pdinfo);

    if (lo == 1) {
	pprintf(prn, "Summary Statistics for the variable '%s' using the "
		"observations %s to %s.\n",
		pdinfo->varname[summ->list[1]], date1, date2);
    } else {
	pprintf(prn, "Summary Statistics, using the observations %s to %s.\n",
		date1, date2);
    }

    for (v=0; v<lo; v++) {
	lv = summ->list[v+1];
	if (lo > 1) {
	    pprintf(prn, "%s, ", pdinfo->varname[lv]);
	}
	pprintf(prn, "mean, %.4g,\n", summ->coeff[v]);
	pprintf(prn, "median, %.4g,\n", summ->xmedian[v]);
	pprintf(prn, "minimum, %.4g,\n", summ->xpx[v]);
	pprintf(prn, "maximum, %.4g,\n", summ->xpy[v]);
	pprintf(prn, "standard deviation, %.4g.\n", summ->sderr[v]);
    }
}

static void
audioprint_matrix (const double *vec, const int *list,
		   int t1, int t2, int n, int ci,
		   const DATAINFO *pdinfo, PRN *prn)
{
    int i, j, k;
    int lo = list[0];

    if (ci == CORR) {
	char date1[OBSLEN], date2[OBSLEN];

	ntodate(date1, t1, pdinfo);
	ntodate(date2, t2, pdinfo);

	pprintf(prn, "Correlation coefficients, using the observations "
		"%s to %s.\n", date1, date2);
	pprintf(prn, " The 5%% critical value (two-tailed) is %.3f.\n", 
		rhocrit95(n));
    } else {
	pputs(prn, "Coefficient covariance matrix.\n");
    }

    for (i=1; i<=lo; i++) {
	for (j=i; j<=lo; j++) {
	    k = ijton(i, j, lo);
	    if (i == j) {
		if (ci == CORR) continue;
		pprintf(prn, "%s, ", pdinfo->varname[list[i]]);
	    } else {
		pprintf(prn, "%s and %s, ", 
			pdinfo->varname[list[i]],
			pdinfo->varname[list[j]]);
	    }
	    if (ci == CORR) {
		pprintf(prn, "%.3f.\n", vec[k]);
	    } else {
		pprintf(prn, "%.4g.\n", vec[k]);
	    }
	}
    }
}

/* ........................................................... */

static void 
audioprint_corrmat (CORRMAT *corr,
		    const DATAINFO *pdinfo, 
		    PRN *prn)
{
    audioprint_matrix(corr->xpx, corr->list, corr->t1, corr->t2,
		      corr->n, CORR, pdinfo, prn);
}

/* .................................................................. */

static void 
audioprint_coeff_interval (const CONFINT *cf, const DATAINFO *pdinfo, 
			   int c, PRN *prn)
{
    pprintf(prn, "Variable %d, '%s', ", cf->list[c], pdinfo->varname[cf->list[c]]);
    pprintf(prn, "point estimate %.4g, 95%% confidence interval, ", cf->coeff[c-2]);

    if (isnan(cf->maxerr[c-2])) {
	pputs(prn, "undefined.\n");	
    } else {
	pprintf(prn, "%.4g to %.4g.\n", 
		cf->coeff[c-2] - cf->maxerr[c-2],
		cf->coeff[c-2] + cf->maxerr[c-2]);
    }
}

/* .................................................................. */

static void 
audioprint_confints (const CONFINT *cf, const DATAINFO *pdinfo, 
		     PRN *prn)
{
    int i, ncoeff = cf->list[0];

    for (i=2; i<=ncoeff; i++) {
	audioprint_coeff_interval(cf, pdinfo, i, prn);
    }
}

/* .................................................................. */

static void 
audioprint_vcv (const VCV *vcv, const DATAINFO *pdinfo, 
		PRN *prn)
{
    audioprint_matrix(vcv->vec, vcv->list, 0, 0, 0,
		      COVAR, pdinfo, prn);
}

/* .................................................................. */

#ifdef HAVE_FLITE

static int speak_buffer (const char *buf, int (*should_stop)())
{
    cst_voice *v;
    char line[128];

    flite_init();
    v = register_cmu_us_kal();

    bufgets(NULL, 0, buf);
    while (bufgets(line, 127, buf)) {
	if (should_stop()) {
	    flite_text_to_speech("OK, stopping", v, "play");
	    break;
	}
	flite_text_to_speech(line, v, "play");
    }
	
    return 0;
}

#else

static int speak_buffer (const char *buf, int (*should_stop)())
{
    ISpVoice *v = NULL;
    HRESULT hr;
    WCHAR *w;
    char line[128];

    hr = CoInitialize(NULL);
    if (!SUCCEEDED(hr)) return 1;
    hr = CoCreateInstance(&CLSID_SpVoice, 
                          NULL, 
                          CLSCTX_ALL, 
                          &IID_ISpVoice, 
                          (void **) &v);
    if (!SUCCEEDED(hr)) {
	CoUninitialize();
	return 1;
    }

    bufgets(NULL, 0, buf);
    while (bufgets(line, 127, buf)) {
	if (should_stop()) {
	    ISpVoice_Speak(v, L"OK, stopping", 0, NULL);
	    break;
	}
	w = wide_string(line);
	ISpVoice_Speak(v, w, 0, NULL);
	free(w);
    }

    ISpVoice_Release(v);
    CoUninitialize();

    return 0;
}

#endif

static int audio_print_special (int role, void *data, const DATAINFO *pdinfo,
				int (*should_stop)())
{
    PRN *prn;

    prn = gretl_print_new(GRETL_PRINT_BUFFER, NULL);
    if (prn == NULL) return 1;
    

    /* descriptive statistics */
    if (role == SUMMARY || role == VAR_SUMMARY) {
	GRETLSUMMARY *summ = (GRETLSUMMARY *) data;

	audioprint_summary(summ, pdinfo, prn);
    }
    else if (role == CORR) {
	CORRMAT *corr = (CORRMAT *) data;

	audioprint_corrmat(corr, pdinfo, prn);
    }
    else if (role == COVAR) {
	VCV *vcv = (VCV *) data;

	audioprint_vcv(vcv, pdinfo, prn);
    }
    else if (role == COEFFINT) {
	CONFINT *cf = (CONFINT *) data;

	audioprint_confints(cf, pdinfo, prn);
    }
    else if (role == VIEW_MODEL) {
	MODEL *pmod = (MODEL *) data;

	audioprint_model(pmod, pdinfo, prn);
    }

    speak_buffer(prn->buf, should_stop);
    gretl_print_destroy(prn);

    return 0;
}

#ifdef GLIB2

int read_window_text (windata_t *vwin, const DATAINFO *pdinfo,
		      int (*should_stop)())
{
    int err = 0;

    if (vwin->role == SUMMARY ||
	vwin->role == VAR_SUMMARY ||
	vwin->role == CORR ||
	vwin->role == COVAR ||
	vwin->role == COEFFINT ||
	vwin->role == VIEW_MODEL) {
	err = audio_print_special(vwin->role, vwin->data, pdinfo, should_stop);
    } else {
	GtkTextBuffer *tbuf;
	GtkTextIter start, end;
	gchar *window_text;

	tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));
	gtk_text_buffer_get_start_iter(tbuf, &start);
	gtk_text_buffer_get_end_iter(tbuf, &end);
	window_text = gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);

	err = speak_buffer(window_text, should_stop);
	g_free(window_text);
    }

    return err;
}

#else

int read_window_text (windata_t *vwin, const DATAINFO *pdinfo,
		      int (*should_stop)())
{
    int err = 0;

    if (vwin->role == SUMMARY ||
	vwin->role == VAR_SUMMARY ||
	vwin->role == CORR ||
	vwin->role == COVAR ||
	vwin->role == VIEW_MODEL) {
	err = audio_print_special(vwin->role, vwin->data, pdinfo, should_stop);
    } else {
	gchar *window_text;

	window_text = gtk_editable_get_chars(GTK_EDITABLE(vwin->w), 0, -1);
	err = speak_buffer(window_text, should_stop);
	g_free(window_text);
    }

    return err;
}

#endif /* GTK versions */

