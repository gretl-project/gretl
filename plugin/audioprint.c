/* print some window contents to a buffer suitable for
   text-to-speech
*/

enum {
    MATRIX_CORR,
    MATRIX_VCV
} matrix_types;

static void 
audioprint_summary (GRETLSUMMARY *summ, const DATAINFO *pdinfo,
		    PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN];
    int lo = summ->list[0], v, lv;

    ntodate(date1, pdinfo->t1, pdinfo);
    ntodate(date2, pdinfo->t2, pdinfo);

    pprintf(prn, "Summary Statistics, using the observations %s to %s\n",
	    date1, date2);

    for (v=0; v<lo; v++) {
	lv = summ->list[v+1];
	pprintf(prn, "%s: ", pdinfo->varname[lv]);
	pprintf(prn, "mean %.4g, ", summ->coeff[v]);
	pprintf(prn, "median %.4g, ", summ->xmedian[v]);
	pprintf(prn, "minimum %.4g, ", summ->xpx[v]);
	pprintf(prn, "maximum %.4g, ", summ->xpy[v]);
	pprintf(prn, "standard deviation %.4g.\n", summ->sderr[v]);
    }
}

static void
audioprint_matrix (const double *vec, const int *list,
		   int t1, int t2, int n, int code,
		   const DATAINFO *pdinfo, PRN *prn)
{
    int i, j, k;
    int lo = list[0];

    if (code == MATRIX_CORR) {
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
		if (code == MATRIX_CORR) continue;
		pprintf(prn, "%s: ", pdinfo->varname[list[i]]);
	    } else {
		pprintf(prn, "%s and %s: ", 
			pdinfo->varname[list[i]],
			pdinfo->varname[list[j]]);
	    }
	    if (code == MATRIX_CORR) {
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
		      corr->n, MATRIX_CORR, pdinfo, prn);
}

/* .................................................................. */

static void 
audioprint_coeff_interval (const CONFINT *cf, const DATAINFO *pdinfo, 
			   int c, PRN *prn)
{
    pprintf(prn, "%s (variable number %d): ", pdinfo->varname[cf->list[c]],
	    cf->list[c]);
    pprintf(prn, "point estimate %.4g, 95%% confidence interval ", cf->coeff[c-2]);

    if (isnan(cf->maxerr[c-2])) {
	pputs(prn, "undefined\n");	
    } else {
	pprintf(prn, "%.4g to %.4g\n", 
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
		      MATRIX_VCV, pdinfo, prn);
}

/* .................................................................. */

#ifdef HAVE_FLITE

static int speak_buffer (const char *buf)
{
    cst_voice *v;

    flite_init();
    v = register_cmu_us_kal();

    flite_text_to_speech(buf, v, "play");

    return 0;
}

#else

static int speak_buffer (const char *buf)
{
    ISpVoice *v = NULL;
    HRESULT hr;
    WCHAR *w;

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

    w = wide_string(dset->comments[i]);
    ISpVoice_Speak(v, w, 0, NULL);

    free(w);
    ISpVoice_Release(v);
    CoUninitialize();

    return 0;
}

#endif

#include "gretltypes.h"

#ifdef GLIB2

int read_window_text (windata_t *vwin, const DATAINFO *pdinfo)
{
    /* descriptive statistics */
    if (vwin->role == SUMMARY) {
	GRETLSUMMARY *summ = (GRETLSUMMARY *) vwin->data;
	PRN *prn;
	
	prn = gretl_print_new(GRETL_PRINT_BUFFER, NULL);
	if (prn == NULL) return 1;

	audioprint_summary(summ, pdinfo, prn);
	speak_buffer(prn->buf);
	gretl_print_destroy(prn);
    }

    /* correlation matrix */
    else if (vwin->role == CORR) {
	CORRMAT *corr = (CORRMAT *) vwin->data;
	PRN *prn;

	prn = gretl_print_new(GRETL_PRINT_BUFFER, NULL);
	if (prn == NULL) return 1;

	audioprint_corrmat(corr, pdinfo, prn);
	speak_buffer(prn->buf);
	gretl_print_destroy(prn);
    }

    else {
	GtkTextBuffer *tbuf;
	GtkTextIter start, end;
	gchar *window_text;

	tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(vwin->w));
	gtk_text_buffer_get_start_iter(tbuf, &start);
	gtk_text_buffer_get_end_iter(tbuf, &end);
	window_text = gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);

	speak_buffer(window_text);
	g_free(window_text);
    }

    return 0;
}

#else

int read_window_text (windata_t *vwin, const DATAINFO *pdinfo)
{
    gchar *window_text;

    window_text = gtk_editable_get_chars(GTK_EDITABLE(vwin->w), 0, -1);

    speak_buffer(window_text);

    g_free(window_text);

    return 0;
}

#endif /* GTK versions */

