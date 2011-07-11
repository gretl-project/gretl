/* print some window contents to a buffer suitable for
   text-to-speech
*/

#include "gretl_enums.h"

static int audioprint_coeff (const DATASET *dset, const MODEL *pmod, 
			     int i, PRN *prn)
{
    char varname[24];

    gretl_model_get_param_name(pmod, dset, i, varname);

    pprintf(prn, "Variable %s.\n", varname);

    /* print coeff value if well-defined */
    if (isnan(pmod->coeff[i]) || na(pmod->coeff[i])) {
	pputs(prn, "Coefficient is undefined.\n");
	return 1;
    } else {
	pprintf(prn, "Coefficient %g.\n", pmod->coeff[i]);
    }

    /* get out if std error is undefined */
    if (isnan(pmod->sderr[i]) || na(pmod->sderr[i])) {
	return 1;
    }

    /* std error is well-defined, but is it positive? */
    if (pmod->sderr[i] > 0.) {
	double tval, pval;

	tval = pmod->coeff[i] / pmod->sderr[i];
	pval = coeff_pval(pmod->ci, tval, pmod->dfd);
	pprintf(prn, "P-value %.3g.\n", pval);
    } else { 
	/* zero standard error */
	pputs(prn, "Standard error is zero.\n");
    }

    return 0;
}

static int audioprint_coefficients (const MODEL *pmod, const DATASET *dset, 
				    PRN *prn)
{
    int i, err = 0, gotnan = 0;
    int n = pmod->ncoeff;

    for (i=0; i<n; i++) {
	err = audioprint_coeff(dset, pmod, i, prn);
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
audioprint_model (MODEL *pmod, const DATASET *dset, PRN *prn)
{
    char startdate[OBSLEN], enddate[OBSLEN];

    if (pmod->ci != OLS) {
	pputs(prn, "Sorry, this model is not O.L.S.  I can't read it.\n");
	return;
    }

    ntodate(startdate, pmod->t1, dset);
    ntodate(enddate, pmod->t2, dset);

    pprintf(prn, "O.L.S estimates using the %d observations %s to %s.\n",
	    pmod->nobs, startdate, enddate);
    pprintf(prn, "Dependent variable %s.\n", dset->varname[pmod->list[1]]);

    audioprint_coefficients(pmod, dset, prn);  

    audio_rsqline(pmod, prn);
}

static void 
audioprint_summary (Summary *summ, const DATASET *dset,
		    PRN *prn)
{
    char date1[OBSLEN], date2[OBSLEN];
    int lo = summ->list[0], i, vi;

    ntodate(date1, dset->t1, dset);
    ntodate(date2, dset->t2, dset);

    if (lo == 1) {
	pprintf(prn, "Summary Statistics for the variable '%s' using the "
		"observations %s to %s.\n",
		dset->varname[summ->list[1]], date1, date2);
    } else {
	pprintf(prn, "Summary Statistics, using the observations %s to %s.\n",
		date1, date2);
    }

    for (i=0; i<summ->list[0]; i++) {
	vi = summ->list[i + 1];
	if (lo > 1) {
	    pprintf(prn, "%s, ", dset->varname[vi]);
	}
	pprintf(prn, "mean, %.4g,\n", summ->mean[i]);
	pprintf(prn, "median, %.4g,\n", summ->median[i]);
	pprintf(prn, "minimum, %.4g,\n", summ->low[i]);
	pprintf(prn, "maximum, %.4g,\n", summ->high[i]);
	pprintf(prn, "standard deviation, %.4g.\n", summ->sd[i]);
    }
}

static void
audioprint_matrix (const VMatrix *vmat, const DATASET *dset,
		   PRN *prn)
{
    int i, j, k;
    int n = vmat->t2 - vmat->t1 + 1;
    int lo = vmat->dim;

    if (vmat->ci == CORR) {
	char date1[OBSLEN], date2[OBSLEN];

	ntodate(date1, vmat->t1, dset);
	ntodate(date2, vmat->t2, dset);

	pprintf(prn, "Correlation coefficients, using the observations "
		"%s to %s.\n", date1, date2);
	pprintf(prn, " The 5%% critical value (two-tailed) is %.3f.\n", 
		rhocrit95(n));
    } else {
	pputs(prn, "Coefficient covariance matrix.\n");
    }

    for (i=1; i<=lo; i++) {
	for (j=i; j<=lo; j++) {
	    k = ijton(i-1, j-1, lo);
	    if (i == j) {
		if (vmat->ci == CORR) {
		    continue;
		}
		pprintf(prn, "%s, ", vmat->names[i-1]);
	    } else {
		pprintf(prn, "%s and %s, ", vmat->names[i-1], 
			vmat->names[j-1]);
	    }
	    if (vmat->ci == CORR) {
		pprintf(prn, "%.3f.\n", vmat->vec[k]);
	    } else {
		pprintf(prn, "%.4g.\n", vmat->vec[k]);
	    }
	}
    }
}

static void 
audioprint_coeff_interval (const CoeffIntervals *cf, int i, PRN *prn)
{
    pprintf(prn, "Variable '%s', ", cf->names[i]);
    pprintf(prn, "point estimate %.4g, 95%% confidence interval, ", cf->coeff[i]);

    if (isnan(cf->maxerr[i])) {
	pputs(prn, "undefined.\n");	
    } else {
	pprintf(prn, "%.4g to %.4g.\n", 
		cf->coeff[i] - cf->maxerr[i],
		cf->coeff[i] + cf->maxerr[i]);
    }
}

static void 
audioprint_confints (const CoeffIntervals *cf, PRN *prn)
{
    int i;

    for (i=0; i<cf->ncoeff; i++) {
	audioprint_coeff_interval(cf, i, prn);
    }
}

#ifdef HAVE_FLITE

static int speak_buffer (const char *buf, int (*should_stop)())
{
    cst_voice *v;
    char line[2048];

    flite_init();
    v = register_cmu_us_kal();

    bufgets_init(buf);

    while (bufgets(line, sizeof line, buf)) {
	if (should_stop()) {
	    flite_text_to_speech("OK, stopping", v, "play");
	    break;
	}
	tailstrip(line);
	flite_text_to_speech(line, v, "play");
    }

    bufgets_finalize(buf);
	
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
                          (void *) &v);
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
    char line[2048];

    v = get_sapi_voice();
    if (v == NULL) return 1;

    bufgets_init(buf);

    while (bufgets(line, sizeof line, buf)) {
	if (should_stop()) {
	    ISpVoice_Speak(v, L"OK, stopping", 0, NULL);
	    break;
	}
	tailstrip(line);
	w = wide_string(line);
	ISpVoice_Speak(v, w, 0, NULL);
	free(w);
    }

    bufgets_finalize(buf);

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

static int audio_print_special (int role, void *data, const DATASET *dset,
				int (*should_stop)())
{
    PRN *prn;
    const char *buf;
    int err = 0;

    prn = gretl_print_new(GRETL_PRINT_BUFFER, &err);
    if (err) return err;
    
    /* descriptive statistics */
    if (role == SUMMARY || role == VAR_SUMMARY) {
	Summary *summ = (Summary *) data;

	audioprint_summary(summ, dset, prn);
    } else if (role == CORR || role == COVAR) {
	VMatrix *vmat = (VMatrix *) data;

	audioprint_matrix(vmat, dset, prn);
    } else if (role == COEFFINT) {
	CoeffIntervals *cf = (CoeffIntervals *) data;

	audioprint_confints(cf, prn);
    } else if (role == VIEW_MODEL) {
	MODEL *pmod = (MODEL *) data;

	audioprint_model(pmod, dset, prn);
    }

    buf = gretl_print_get_buffer(prn);

    speak_buffer(buf, should_stop);
    gretl_print_destroy(prn);

    return 0;
}

static int read_listbox_content (GtkWidget *listbox, int (*should_stop)())
{
    GtkTreeModel *model;
    GtkTreeIter iter;
    gchar *line;
    gchar *tmpstr[3];
    int i, err = 0;

    model = gtk_tree_view_get_model(GTK_TREE_VIEW(listbox));
    gtk_tree_model_get_iter_first(model, &iter);

    err = speak_line("Contents of list box.\n");

    while (!err) {
	tmpstr[0] = tmpstr[1] = tmpstr[2] = NULL;

	if (!GTK_IS_TREE_MODEL(model)) break;

	gtk_tree_model_get(model, &iter, 
			   0, &tmpstr[0], 
			   1, &tmpstr[1],
			   2, &tmpstr[2],
			   -1);

	line = g_strdup_printf("%s. %s. %s.\n", 
			       ((tmpstr[0] != NULL && *tmpstr[0] != '\0')?
				tmpstr[0] : "empty column"),
			       ((tmpstr[1] != NULL && *tmpstr[1] != '\0')?
				tmpstr[1] : "empty column"),
			       ((tmpstr[2] != NULL && *tmpstr[2] != '\0')?
				tmpstr[2] : "empty column"));

	err = speak_line(line);

	for (i=0; i<3; i++) {
	    g_free(tmpstr[i]); 
	}
	g_free(line);

	if (!GTK_IS_TREE_MODEL(model) ||
	    !gtk_tree_model_iter_next(model, &iter) || 
	    should_stop()) {
	    break;
	}
    }

#ifdef G_OS_WIN32
    speak_line(NULL);
#endif

    return err;
}

int read_window_text (GtkWidget *listbox, 
		      GtkWidget *text, 
		      int role,
		      gpointer data,
		      const DATASET *dset,
		      int (*should_stop)())
{
    int err = 0;

    if (dset == NULL) {
	return read_listbox_content(listbox, should_stop);
    }

    if (role == SUMMARY ||
	role == VAR_SUMMARY ||
	role == CORR ||
	role == COVAR ||
	role == COEFFINT ||
	role == VIEW_MODEL) {
	err = audio_print_special(role, data, dset, should_stop);
    } else {
	GtkTextBuffer *tbuf;
	GtkTextIter start, end;
	gchar *window_text;

	tbuf = gtk_text_view_get_buffer(GTK_TEXT_VIEW(text));
	gtk_text_buffer_get_start_iter(tbuf, &start);
	gtk_text_buffer_get_end_iter(tbuf, &end);
	window_text = gtk_text_buffer_get_text(tbuf, &start, &end, FALSE);

	err = speak_buffer(window_text, should_stop);
	g_free(window_text);
    }

    return err;
}

