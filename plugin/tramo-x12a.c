/* TRAMO/SEATS, X-12-ARIMA plugin for gretl */

#include "libgretl.h"
#include <gtk/gtk.h>

#ifdef OS_WIN32
# include <windows.h>
#endif

typedef enum {
    D11,     /* seasonally adjusted series */
    D12,     /* trend/cycle */
    D13,     /* irregular component */
    TRIGRAPH /* graph showing all of the above */
} x12a_objects;

const char *x12a_series_strings[] = {
    "d11", "d12", "d13"
};

const char *x12a_descrip_formats[] = {
    N_("seasonally adjusted %s"),
    N_("trend/cycle for %s"),
    N_("irregular component of %s")
};

const char *default_mdl = {
    "# ARIMA specifications that will be tried\n"
    "(0 1 1)(0 1 1) X\n"
    "(0 1 2)(0 1 1) X\n"
    "(2 1 0)(0 1 1) X\n"
    "(0 2 2)(0 1 1) X\n"
    "(2 1 2)(0 1 1)\n"
};  

typedef struct {
    GtkWidget *check;
    char save;
    unsigned short v;
} opt_info;

typedef struct {
    GtkWidget *dialog;
    opt_info opt[4];
    int savevars;
} x12a_request;

#if GTK_MAJOR_VERSION >= 2

static int x12a_dialog (x12a_request *request)
{
    GtkWidget *vbox, *tmp;
    gint ret;

    request->dialog = gtk_dialog_new_with_buttons ("X-12-ARIMA",
					  NULL,
					  GTK_DIALOG_MODAL | 
					  GTK_DIALOG_DESTROY_WITH_PARENT,
					  GTK_STOCK_OK,
					  GTK_RESPONSE_ACCEPT,
					  GTK_STOCK_CANCEL,
					  GTK_RESPONSE_REJECT,
					  NULL);

    tmp = gtk_label_new ("Save to data set");
    gtk_widget_show(tmp);
    vbox = gtk_vbox_new(FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 5);

    tmp = gtk_check_button_new_with_label(_("Seasonally adjusted series"));
    gtk_widget_show(tmp);
    gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 5);
    request->opt[D11].check = tmp;
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), TRUE);

    tmp = gtk_check_button_new_with_label(_("Trend/cycle"));
    gtk_widget_show(tmp);
    gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 5);
    request->opt[D12].check = tmp;
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), TRUE);

    tmp = gtk_check_button_new_with_label(_("Irregular"));
    gtk_widget_show(tmp);
    gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 5);
    request->opt[D13].check = tmp;
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), TRUE);

    tmp = gtk_hseparator_new();
    gtk_widget_show(tmp);
    gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 5);
    
    tmp = gtk_check_button_new_with_label(_("Generate graph"));
    gtk_widget_show(tmp);
    gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 5);
    request->opt[TRIGRAPH].check = tmp;
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), TRUE);

    gtk_widget_show(vbox);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(request->dialog)->vbox),
		       vbox, FALSE, FALSE, 5);

    ret = gtk_dialog_run (GTK_DIALOG(request->dialog));

    if (ret == GTK_RESPONSE_ACCEPT) return 1;

    return 0;
}

#endif

int write_tramo_data (char *fname, int varnum, 
		      double ***pZ, DATAINFO *pdinfo, 
		      PATHS *paths, int *graph,
		      const char *tramo,
		      const char *tramodir)
{
    int i, t;
    int tsamp = pdinfo->t2 - pdinfo->t1 + 1;
    double x;
    char tmp[9], varname[9], cmd[MAXLEN];
    int startyr, startper;
    FILE *fp = NULL;

    *gretl_errmsg = 0;

    if (!pdinfo->vector[varnum]) {
	sprintf(gretl_errmsg, "%s %s", pdinfo->varname[varnum], 
		_("is a scalar"));
	return 1;
    }

    sprintf(varname, pdinfo->varname[varnum]);
    lower(varname);
    sprintf(fname, "%s%c%s", tramodir, SLASH, varname);

    fp = fopen(fname, "w");
    if (fp == NULL) return 1;

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif   

    x = date(pdinfo->t1, pdinfo->pd, pdinfo->sd0);
    startyr = (int) x;
    sprintf(tmp, "%g", x);
    startper = atoi(strchr(tmp, '.') + 1);

    fprintf(fp, "%s\n", pdinfo->varname[varnum]);
    fprintf(fp, "%d %d %d %d\n", tsamp, startyr, startper, pdinfo->pd);

    for (t=pdinfo->t1; t<=pdinfo->t2; ) {
	for (i=0; i<pdinfo->pd; i++) {
	    if (t == pdinfo->t1) {
		i += startper - 1;
	    }
	    if (na((*pZ)[varnum][t])) {
		fprintf(fp, "-99999 ");
	    } else {
		fprintf(fp, "%g ", (*pZ)[varnum][t]);
	    }
	    t++;
	}
	fputs("\n", fp);
    }

    /* FIXME: make these values configurable */
    fprintf(fp, "$INPUT lam=-1,iatip=1,aio=2,va=3.3,noadmiss=1,seats=2,$\n");

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    if (fp != NULL) {
	fclose(fp);
    }

    /* testing */
#ifdef notdef
    sprintf(cmd, "cd %s && ./tramo -i %s >/dev/null && mv seats.itr serie && ./seats -OF %s", 
	    tramodir, varname, varname);
#endif
#ifdef OS_WIN32 /* FIXME */
    /* get into tramodir first */
    sprintf(cmd, "\"%s\" -i %s", tramo, varname);
    WinExec(cmd, SW_SHOWMINIMIZED);
#else
    sprintf(cmd, "cd \"%s\" && \"%s\" -i %s >/dev/null", tramodir, tramo, varname);
    system(cmd);
#endif

    strcpy(tmp, pdinfo->varname[varnum]);
    lower(tmp);
    sprintf(fname, "%s%coutput%c%s.out", tramodir, SLASH, SLASH, varname); /* .OUT for seats */

    return 0;
}

static void truncate (char *str, int n)
{
    int len = strlen(str);

    if (len > n) str[n] = 0;
}

static int graph_x12a_series (double **Z, DATAINFO *pdinfo, 
			      PATHS *paths)
{
    FILE *fp = NULL;
    int t;

    if (gnuplot_init(paths, &fp)) return E_FOPEN;

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    /* FIXME tics? */

    fputs("# X-12-ARIMA tri-graph (no auto-parse)\n", fp);

    fputs("set multiplot\nset size 1.0,0.32\n", fp);

    /* irregular component */
    fprintf(fp, "set bars 0\nset origin 0.0,0.0\n"
	    "plot '-' using :(1.0):($1-1.0) w yerrorbars title '%s', \\\n"
	    "1.0 notitle\n", I_("irregular"));
    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	fprintf(fp, "%g\n", Z[D13 + 1][t]);
    }
    fprintf(fp, "e\n");

    /* actual vs trend/cycle */
    fprintf(fp, "set origin 0.0,0.33\n"
	    "plot '-' using 1 w l title '%s', \\\n"
	    " '-' using 1 w l title '%s'\n",
	    pdinfo->varname[0], I_("trend/cycle"));
    for (t=pdinfo->t1; t<=pdinfo->t2; t++) 
	fprintf(fp, "%g\n", Z[0][t]);
    fprintf(fp, "e , \\\n");
    for (t=pdinfo->t1; t<=pdinfo->t2; t++) 
	fprintf(fp, "%g\n", Z[D12 + 1][t]);
    fprintf(fp, "e\n");

    /* actual vs seasonally adjusted */
    fprintf(fp, "set origin 0.0,0.66\n"
	    "plot '-' using 1 w l title '%s', \\\n"
	    " '-' using 1 w l title '%s'\n",
	    pdinfo->varname[0], I_("adjusted"));
    for (t=pdinfo->t1; t<=pdinfo->t2; t++) 
	fprintf(fp, "%g\n", Z[0][t]);
    fprintf(fp, "e\n");
    for (t=pdinfo->t1; t<=pdinfo->t2; t++) 
	fprintf(fp, "%g\n", Z[D11 + 1][t]);

    fprintf(fp, "unset multiplot\n");

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

#if defined(OS_WIN32) && !defined(GNUPLOT_PNG)
    fprintf(fp, "pause -1\n");
#endif

    fclose(fp);

    return 0;
}

static void copy_variable (double **targZ, DATAINFO *targinfo, int targv,
			   double **srcZ, DATAINFO *srcinfo, int srcv)
{
    int t;

    for (t=0; t<targinfo->n; t++) {
	targZ[targv][t] = srcZ[srcv][t];
    }

    strcpy(targinfo->varname[targv], srcinfo->varname[srcv]);
    strcpy(targinfo->label[targv], srcinfo->label[srcv]);
}

static int add_x12a_series (const char *fname, int code,
			    double **Z, DATAINFO *pdinfo,
			    int v)
{
    FILE *fp;
    char *p, line[128], varname[16], sfname[MAXLEN], date[8];
    double x;
    int d, yr, per, err = 0;
    int t;

    strcpy(sfname, fname);
    p = strrchr(sfname, '.');
    if (p != NULL) strcpy(p + 1, x12a_series_strings[code]);

    fp = fopen(sfname, "r");
    if (fp == NULL) {
	sprintf(gretl_errmsg, "%s %s", _("Couldn't open"), sfname);
	return 1;
    }

    /* formulate name of new variable to add */
    strcpy(varname, pdinfo->varname[0]);
    truncate(varname, 4);
    strcat(varname, "_");
    strcat(varname, x12a_series_strings[code]);

    /* copy varname and label into place */
    strcpy(pdinfo->varname[v], varname);
    sprintf(pdinfo->label[v], _(x12a_descrip_formats[code]), pdinfo->varname[0]);

    for (t=0; t<pdinfo->n; t++) Z[v][t] = NADBL;

    /* grab the data from the x12arima file */
    while (fgets(line, 127, fp)) {
	if (*line == 'd' || *line == '-') continue;
	if (sscanf(line, "%d %lf", &d, &x) != 2) {
	    err = 1; 
	    break;
	}
	yr = d / 100;
	per = d % 100;
	sprintf(date, "%d.%d", yr, per);
	t = dateton(date, pdinfo);
	/* fprintf(stderr, "date='%s', t=%d\n", date, t); */
	if (t < 0 || t >= pdinfo->n) {
	    err = 1;
	    break;
	}
	Z[v][t] = x;
    }

    fclose(fp);

    return err;
}

static void set_opts (x12a_request *request)
{
    int i;

    request->savevars = 0;

    for (i=0; i<4; i++) {
	if (GTK_TOGGLE_BUTTON(request->opt[i].check)->active) {
	    request->opt[i].save = 1;
	    if (i != TRIGRAPH) request->savevars++;
	} else {
	    request->opt[i].save = 0;
	} 
    }
}

static int write_spc_file (const char *fname, 
			   double **Z, DATAINFO *pdinfo, 
			   int varnum, int *varlist, 
			   x12a_request *request)
{
    double x;
    FILE *fp;
    int i, j, t;
    int startyr, startper;
    char tmp[8];

    fp = fopen(fname, "w");
    if (fp == NULL) return 1;    

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif 

    x = date(pdinfo->t1, pdinfo->pd, pdinfo->sd0);
    startyr = (int) x;
    sprintf(tmp, "%g", x);
    startper = atoi(strchr(tmp, '.') + 1);

    fprintf(fp, "series{\n period=%d\n title=\"%s\"\n", pdinfo->pd, 
	    pdinfo->varname[varnum]);
    fprintf(fp, " start=%d.%d\n", startyr, startper);
    fprintf(fp, " data=(\n");

    i = 0;
    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (na(Z[varnum][t])) {
	    fprintf(fp, "-99999 "); /* FIXME? */
	} else {
	    fprintf(fp, "%g ", Z[varnum][t]);
	}
	if ((i + 1) % 7 == 0) fputs("\n", fp);
	i++;
    }
    fputs(" )\n}\n", fp);

    /* FIXME: make these values configurable */
    fputs("automdl{}\nx11{", fp);

    /* construct list of vars to be grabbed */
    varlist[0] = 0;
    j = 1;
    for (i=0; i<TRIGRAPH; i++) {
	if (request->opt[TRIGRAPH].save || request->opt[i].save) {
	    varlist[0] += 1;
	    varlist[j++] = i;
	}
    }
    
    if (varlist[0] > 0) {
	if (varlist[0] == 1) {
	    fprintf(fp, " save=%s ", x12a_series_strings[varlist[1]]); 
	} else {
	    fputs(" save=( ", fp);
	    for (i=1; i<=varlist[0]; i++) {
		fprintf(fp, "%s ", x12a_series_strings[varlist[i]]);
	    }
	    fputs(") ", fp);
	}
    }

    fputs("}\n", fp);

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    fclose(fp);

    return 0;
}

static void copy_basic_data_info (DATAINFO *targ, DATAINFO *src)
{
    targ->sd0 = src->sd0;
    strcpy(targ->stobs, src->stobs); 
    targ->t1 = src->t1;
    targ->t2 = src->t2;
    targ->pd = src->pd;
    targ->time_series = src->time_series;
}

static int save_vars_to_dataset (double ***pZ, DATAINFO *pdinfo,
				 double **x12Z, DATAINFO *x12info,
				 int *varlist, x12a_request *request)
{
    int i, v, j, addvars = 0;

    /* how many vars are wanted, and new? */
    for (i=1; i<=varlist[0]; i++) {
	if (request->opt[varlist[i]].save && 
	    varindex(pdinfo, x12info->varname[i]) == pdinfo->v) {
	    addvars++;
	}
    }

    if (addvars > 0 && dataset_add_vars(addvars, pZ, pdinfo)) {
	strcpy(gretl_errmsg, _("Failed to allocate memory for new data"));
	return 1;
    }

    j = pdinfo->v - addvars;
    for (i=1; i<=varlist[0]; i++) {
	if (request->opt[varlist[i]].save) {
	    v = varindex(pdinfo, x12info->varname[i]);
	    if (v < pdinfo->v) {
		copy_variable(*pZ, pdinfo, v, x12Z, x12info, i);
	    } else {
		copy_variable(*pZ, pdinfo, j++, x12Z, x12info, i);
	    }
	}
    }

    return 0;
}

int write_x12a_data (char *fname, int varnum, 
		     double ***pZ, DATAINFO *pdinfo, 
		     PATHS *paths, int *graph,
		     const char *x12a, const char *x12adir)
{
    int i, doit, err = 0;
    char varname[9], cmd[MAXLEN];
    int varlist[4];
    FILE *fp = NULL;
    x12a_request request;
    double **x12Z;
    DATAINFO *x12info;

    /* sanity check */
    *gretl_errmsg = 0;
    if (!pdinfo->vector[varnum]) {
	sprintf(gretl_errmsg, "%s %s", pdinfo->varname[varnum], 
		_("is a scalar"));
	return 1;
    }

    /* show dialog and get option settings */
    doit = x12a_dialog(&request); 
    if (!doit) {
	gtk_widget_destroy(request.dialog);
	return 0;
    }
    set_opts(&request);
    gtk_widget_destroy(request.dialog);

    /* create little temporary dataset */
    x12info = create_new_dataset(&x12Z, 4, pdinfo->n, 0);
    if (x12info == NULL) return E_ALLOC;
    copy_basic_data_info(x12info, pdinfo);

    /* make a default x12a.mdl file if it doesn't already exist */
    sprintf(fname, "%s%cx12a.mdl", x12adir, SLASH);
    fp = fopen(fname, "r");
    if (fp == NULL) {
	fp = fopen(fname, "w");
	if (fp == NULL) return 1;
	fprintf(fp, "%s", default_mdl);
	fclose(fp);
    } else {
	fclose(fp);
    }

    sprintf(varname, pdinfo->varname[varnum]);

    /* write out the .spc file for x12a */
    sprintf(fname, "%s%c%s.spc", x12adir, SLASH, varname);
    write_spc_file(fname, *pZ, pdinfo, varnum, varlist, &request);

    /* run the x12a program -- FIXME win32 */
    sprintf(cmd, "cd %s && %s %s -r -p -q >/dev/null", x12adir, x12a, varname);
    err = system(cmd);
    
    if (!err) {
	sprintf(fname, "%s%c%s.out", x12adir, SLASH, varname); 

	/* save vars locally if needed; graph if wanted */
	if (varlist[0] > 0) {
	    copy_variable(x12Z, x12info, 0, *pZ, pdinfo, varnum);
	    for (i=1; i<=varlist[0]; i++) {
		err = add_x12a_series(fname, varlist[i], x12Z, x12info, i);
	    }
	    if (request.opt[TRIGRAPH].save) {
		err = graph_x12a_series(x12Z, x12info, paths);
		if (!err) *graph = 1;
	    }
	}

	/* now save the local vars to main dataset, if wanted */
	if (request.savevars > 0) {
	    err = save_vars_to_dataset(pZ, pdinfo, x12Z, x12info, varlist, &request);
	}
    }

    free_Z(x12Z, x12info);
    clear_datainfo(x12info, CLEAR_FULL);
    free(x12info);

    return err;
}

