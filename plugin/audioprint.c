/* print some window contents to a buffer suitable for
   text-to-speech
*/

#include "gretltypes.h"
#include "gretl_enums.h"

static int audioprint_coeff (const DATAINFO *pdinfo, const MODEL *pmod, 
			     int c, PRN *prn)
{
    pprintf(prn, "Variable %d, %s.\n", pmod->list[c], 
	    pdinfo->varname[pmod->list[c]]);

    /* print coeff value if well-defined */
    if (isnan(pmod->coeff[c-2]) || na(pmod->coeff[c-2])) {
	pputs(prn, "Coefficient is undefined.\n");
	return 1;
    } else {
	pprintf(prn, "Coefficient %g.\n", pmod->coeff[c-2]);
    }

    /* get out if std error is undefined */
    if (isnan(pmod->sderr[c-2]) || na(pmod->sderr[c-2])) {
	return 1;
    }

    /* std error is well-defined, but is it positive? */
    if (pmod->sderr[c-2] > 0.) {
	double t, pval;

	t = pmod->coeff[c-2] / pmod->sderr[c-2];
	pval = tprob(t, pmod->dfd);
	pprintf(prn, "P-value %.3g.\n", pval);
    } else { /* zero standard error */
	pputs(prn, "Standard error is zero.\n");
    }

    return 0;
}

static int audioprint_coefficients (const MODEL *pmod, const DATAINFO *pdinfo, 
				    PRN *prn)
{
    int i, err = 0, gotnan = 0;
    int n = pmod->ncoeff;

    for (i=0; i<n; i++) {
	err = audioprint_coeff(pdinfo, pmod, i + 2, prn);
	if (err) gotnan = 1;
    }

    return gotnan;
} 

static void audio_rsqline (const MODEL *pmod, PRN *prn)
{
    if (!na(pmod->rsq)) {
	pprintf(prn, "Unadjusted R-squared %.3f.\n", pmod->rsq);
    }
    if (!na(pmod->adjrsq)) {
	pprintf(prn, "Adjusted R-squared %.3f.\n", pmod->adjrsq);
    }
}

static void 
audioprint_model (MODEL *pmod, const DATAINFO *pdinfo, PRN *prn)
{
    char startdate[OBSLEN], enddate[OBSLEN];
    int t1 = pmod->t1, t2 = pmod->t2;

    if (pmod->ci != OLS) {
	pputs(prn, "Sorry, this model is not O.L.S.  I can't read it.\n");
	return;
    }

    if (pmod->data != NULL) {
	t2 += get_misscount(pmod);
    }

    ntodate(startdate, t1, pdinfo);
    ntodate(enddate, t2, pdinfo);

    pprintf(prn, "O.L.S estimates using the %d observations %s to %s.\n",
	    pmod->nobs, startdate, enddate);
    pprintf(prn, "Dependent variable %s.\n", pdinfo->varname[pmod->list[1]]);

    audioprint_coefficients(pmod, pdinfo, prn);  

    audio_rsqline(pmod, prn);
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

static int speak_line (const char *line)
{
    static cst_voice *v;

    if (v == NULL) {
	flite_init();
	v = register_cmu_us_kal();
    }

    if (v == NULL) return 1;

    flite_text_to_speech(line, v, "play");

    return 0;
}

#else

static ISpVoice *get_sapi_voice (void)
{
    ISpVoice *v;
    HRESULT hr;

    hr = CoInitialize(NULL);
    if (!SUCCEEDED(hr)) return NULL;

    hr = CoCreateInstance(&CLSID_SpVoice, 
                          NULL, 
                          CLSCTX_ALL, 
                          &IID_ISpVoice, 
                          (void **) &v);
    if (!SUCCEEDED(hr)) {
	CoUninitialize();
	return NULL;
    }

    return v;
}

static void release_sapi_voice (ISpVoice *v)
{
    ISpVoice_Release(v);
    CoUninitialize();
}

static int speak_buffer (const char *buf, int (*should_stop)())
{
    ISpVoice *v = NULL;
    WCHAR *w;
    char line[128];

    v = get_sapi_voice();
    if (v == NULL) return 1;

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

    release_sapi_voice(v);

    return 0;
}

static int speak_line (const char *line)
{
    static ISpVoice *v = NULL;

    if (line == NULL) {
	if (v != NULL) {
	    release_sapi_voice(v);
	}
	return 0;
    }

    if (v == NULL) {
	v = get_sapi_voice();
    }

    if (v != NULL) {
	WCHAR *w = wide_string(line);

	ISpVoice_Speak(v, w, 0, NULL);
	free(w);
    } else {
	return 1;
    }

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

static int read_listbox_content (windata_t *vwin, int (*should_stop)())
{
    GtkTreeModel *model;
    GtkTreeIter iter;
    gchar *line;
    gchar *tmpstr[3];
    int i, err = 0;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(vwin->listbox));
    gtk_tree_model_get_iter_first(model, &iter);

    err = speak_line("Contents of list box.\n");

    while (!err) {
	tmpstr[0] = tmpstr[1] = tmpstr[2] = NULL;

	gtk_tree_model_get(model, &iter, 
			   0, &tmpstr[0], 
			   1, &tmpstr[1],
			   2, &tmpstr[2],
			   -1);

	line = g_strdup_printf("%s. %s. %s.\n", tmpstr[0], 
			       ((tmpstr[1] != NULL && *tmpstr[1] != '\0')?
			       tmpstr[1] : "empty column"),
			       ((tmpstr[2] != NULL && *tmpstr[2] != '\0')?
			       tmpstr[2] : "empty column"));

	err = speak_line(line);

	for (i=0; i<3; i++) {
	    g_free(tmpstr[i]); 
	}
	g_free(line);

	if (!gtk_tree_model_iter_next(model, &iter) || should_stop()) {
	    break;
	}
    }

#ifdef G_OS_WIN32
    speak_line(NULL);
#endif

    return err;
}

int read_window_text (windata_t *vwin, const DATAINFO *pdinfo,
		      int (*should_stop)())
{
    int err = 0;

    if (pdinfo == NULL) {
	return read_listbox_content(vwin, should_stop);
    }

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

static int read_listbox_content (windata_t *vwin, int (*should_stop)())
{
    gchar *line;
    gchar *tmpstr[3];
    int i, err = 0;

    err = speak_line("Contents of list box.\n");

    for (i=0; !err; i++) {
	if (!gtk_clist_get_text(GTK_CLIST(vwin->listbox), i, 0, &tmpstr[0]) ||
	    !gtk_clist_get_text(GTK_CLIST(vwin->listbox), i, 1, &tmpstr[1]) ||
	    !gtk_clist_get_text(GTK_CLIST(vwin->listbox), i, 2, &tmpstr[2])) {
	    err = 1;
	}

	if (!err) {
	    line = g_strdup_printf("%s. %s. %s.\n", tmpstr[0], tmpstr[1], tmpstr[2]);
	    err = speak_line(line);
	    g_free(line);
	}

	if (!err && should_stop()) {
	    break;
	}
    }

    return err;
}

int read_window_text (windata_t *vwin, const DATAINFO *pdinfo,
		      int (*should_stop)())
{
    int err = 0;

    if (pdinfo == NULL) {
	return read_listbox_content(vwin, should_stop);
    }

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

