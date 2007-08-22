/* 
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2001 Allin Cottrell and Riccardo "Jack" Lucchetti
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

/* TRAMO/SEATS, X-12-ARIMA plugin for gretl */

#include "libgretl.h"
#include <gtk/gtk.h>
#include "tramo_x12a.h"

#ifdef WIN32
# include <windows.h>
#else
# include <signal.h>
#endif

enum opt_codes {
    TRAMO_SEATS,
    TRAMO_ONLY,
    X12A
};

const char *x12a_series_strings[] = {
    "d11", "d12", "d13"
};

static int tramo_got_irfin;

const char *tramo_series_strings[] = {
    "safin.t", "trfin.t", "irfin.t", "irreg.t"
};

const char *tx_descrip_formats[] = {
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

#ifndef WIN32

#define SP_DEBUG 0

static int glib_spawn (const char *workdir, const char *fmt, ...)
{
    GError *error = NULL;
    gchar *sout = NULL;
    gchar *serr = NULL;
    gchar *argv[8];

    va_list ap;
    int i, ok, nargs;
    int status = 0, ret = 0;
    char *s;

    argv[0] = g_strdup(fmt);
    argv[1] = NULL;
    i = nargs = 1;

    va_start(ap, fmt);

    while ((s = va_arg(ap, char *))) {
	argv[i] = g_strdup(s);
	argv[++i] = NULL;
    }

    va_end(ap);

    nargs = i;

#if SP_DEBUG
    fputs("spawning the following:\n", stderr);
    for (i=0; i<nargs; i++) {
	fprintf(stderr, " argv[%d] = '%s'\n", i, argv[i]);
    }
#endif

    signal(SIGCHLD, SIG_DFL);

    ok = g_spawn_sync (workdir,
		       argv,
		       NULL,
		       G_SPAWN_SEARCH_PATH,
		       NULL,
		       NULL,
		       &sout,
		       &serr,
		       &status,
		       &error);

    if (!ok) {
	fprintf(stderr, "spawn: '%s'\n", error->message);
	g_error_free(error);
	ret = 1;
    } else if (serr && *serr) {
	fprintf(stderr, "stderr: '%s'\n", serr);
	ret = 1;
    } else if (status != 0) {
	fprintf(stderr, "status=%d: stdout: '%s'\n", status, sout);
	ret = 1;
    }

    if (serr != NULL) g_free(serr);
    if (sout != NULL) g_free(sout);

    for (i=0; i<nargs; i++) {
	if (ret != 0) {
	    if (i == 0) {
		fputc(' ', stderr);
	    }
	    fprintf(stderr, "%s ", argv[i]);
	    if (i == nargs - 1) {
		fputc('\n', stderr);
	    }
	}
	free(argv[i]);
    }

    return ret;
}

#endif /* !WIN32 */

static int tx_dialog (tx_request *request)
{
    GtkWidget *hbox, *vbox, *tmp;
    gint i, ret = 0;

    for (i=0; i<N_COMMON_OPTS; i++) {
	request->opt[i].check = NULL;
    }

    request->dialog = 
	gtk_dialog_new_with_buttons (
				     (request->code == TRAMO_SEATS)?
				     "TRAMO/SEATS" : "X-12-ARIMA",
				     NULL,
				     GTK_DIALOG_MODAL | 
				     GTK_DIALOG_DESTROY_WITH_PARENT,
				     GTK_STOCK_CANCEL,
				     GTK_RESPONSE_REJECT,
				     GTK_STOCK_OK,
				     GTK_RESPONSE_ACCEPT,
				     NULL);

    vbox = gtk_vbox_new(FALSE, 0);    

    if (request->code == TRAMO_SEATS) {
	gtk_dialog_set_has_separator(GTK_DIALOG(request->dialog), FALSE);
	show_tramo_options(request, vbox);
    } else {
	/* X-12-ARIMA */
	tmp = gtk_label_new (_("Save data"));
	gtk_widget_show(tmp);
	gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 5);

	tmp = gtk_check_button_new_with_label(_("Seasonally adjusted series"));
	gtk_widget_show(tmp);
	gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 5);
	request->opt[D11].check = tmp;
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), FALSE);

	tmp = gtk_check_button_new_with_label(_("Trend/cycle"));
	gtk_widget_show(tmp);
	gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 5);
	request->opt[D12].check = tmp;
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), FALSE);

	tmp = gtk_check_button_new_with_label(_("Irregular"));
	gtk_widget_show(tmp);
	gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 5);
	request->opt[D13].check = tmp;
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), FALSE);

	tmp = gtk_hseparator_new();
	gtk_widget_show(tmp);
	gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 5);
    
	tmp = gtk_check_button_new_with_label(_("Generate graph"));
	gtk_widget_show(tmp);
	gtk_box_pack_start(GTK_BOX(vbox), tmp, FALSE, FALSE, 5);
	request->opt[TRIGRAPH].check = tmp;
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), TRUE);
    }

    gtk_widget_show(vbox);
    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 5);
    gtk_widget_show(hbox);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(request->dialog)->vbox),
		       hbox, FALSE, FALSE, 5);

    ret = gtk_dialog_run(GTK_DIALOG(request->dialog));

    if (ret == GTK_RESPONSE_ACCEPT) ret = 1;
    else ret = 0;

    return ret;
}

static void get_seats_command (char *seats, const char *tramo)
{
    char *p;

    strcpy(seats, tramo);
    p = strrchr(seats, SLASH);
    if (p != NULL) {
	strcpy(p + 1, "seats");
    } else {
	strcpy(seats, "seats");
    }
}

static int graph_series (double **Z, DATAINFO *pdinfo, int opt)
{
    FILE *fp = NULL;
    const double *obs;
    char title[32];
    int t;

    obs = gretl_plotx(pdinfo);
    if (obs == NULL) {
	return E_ALLOC;
    }

    if (gnuplot_init(PLOT_TRI_GRAPH, &fp)) {
	return E_FOPEN;
    }

    gretl_push_c_numeric_locale();

    if (pdinfo->pd == 4) {
	if ((pdinfo->t2 - pdinfo->t1) / 4 < 8) {
	    fputs("set xtics nomirror 0,1\n", fp); 
	    fputs("set mxtics 4\n", fp);
	}
    } else if (pdinfo->pd == 12) {
	if ((pdinfo->t2 - pdinfo->t1) / 12 < 8) {
	    fputs("set xtics nomirror 0,1\n", fp); 
	    fputs("set mxtics 12\n", fp);
	}
    }

    fputs("set size 1.0,1.0\nset multiplot\nset size 1.0,0.32\n", fp);

    /* irregular component */
    if (opt == TRAMO_SEATS && !tramo_got_irfin) {
	sprintf(title, "%s", G_("irregular"));
    } else {
	sprintf(title, "%s - 1", G_("irregular"));
    }

    fprintf(fp, "set bars 0\n"
	    "set origin 0.0,0.0\n"
	    "plot '-' using 1:($2-1.0) title '%s' w impulses\n",
	    title);

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	double y = Z[D13 + 1][t];

	fprintf(fp, "%g %g\n", obs[t], 
		(opt == TRAMO_SEATS && tramo_got_irfin)? y / 100.0 : y);
    }
    fputs("e\n", fp);

    /* actual vs trend/cycle */

    fprintf(fp, "set origin 0.0,0.33\n"
	    "plot '-' using 1:2 title '%s' w l, \\\n"
	    " '-' using 1:2 title '%s' w l\n",
	    pdinfo->varname[0], G_("trend/cycle"));

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) { 
	fprintf(fp, "%g %g\n", obs[t], Z[0][t]);
    }
    fputs("e , \\\n", fp);

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) { 
	fprintf(fp, "%g %g\n", obs[t], Z[D12 + 1][t]);
    }
    fputs("e\n", fp);

    /* actual vs seasonally adjusted */

    fprintf(fp, "set origin 0.0,0.66\n"
	    "plot '-' using 1:2 title '%s' w l, \\\n"
	    " '-' using 1:2 title '%s' w l\n",
	    pdinfo->varname[0], G_("adjusted"));

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	fprintf(fp, "%g %g\n", obs[t], Z[0][t]);
    }
    fputs("e\n", fp);

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	fprintf(fp, "%g %g\n", obs[t], Z[D11 + 1][t]);
    }
    fputs("e\n", fp);

    fputs("set nomultiplot\n", fp);

    gretl_pop_c_numeric_locale();

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
    strcpy(VARLABEL(targinfo, targv), VARLABEL(srcinfo, srcv));
}

static void clear_tramo_files (const char *tpath, const char *varname)
{
    char tfname[MAXLEN];
    int i;

    for (i=D11; i<=TRIGRAPH; i++) {
	sprintf(tfname, "%s%cgraph%cseries%c%s", tpath, SLASH, SLASH, SLASH,
		tramo_series_strings[i]);
	remove(tfname);
    }

    sprintf(tfname, "%s%coutput%c%s.out", tpath, SLASH, SLASH, varname);
    remove(tfname);
}

static int add_series_from_file (const char *fname, int code,
				 double **Z, DATAINFO *pdinfo,
				 int v, int opt, char *errmsg)
{
    FILE *fp;
    char *p, line[128], sfname[MAXLEN], date[8];
    char varname[VNAMELEN];
    double x;
    int d, yr, per, err = 0;
    int t;

    if (opt == TRAMO_SEATS) {
	tramo_got_irfin = 1;
	sprintf(sfname, "%s%cgraph%cseries%c%s", fname, SLASH, SLASH, SLASH,
		tramo_series_strings[code]);
    } else {
	strcpy(sfname, fname);
	p = strrchr(sfname, '.');
	if (p != NULL) strcpy(p + 1, x12a_series_strings[code]);
    }

    fp = gretl_fopen(sfname, "r");

    if (fp == NULL) {
	int gotit = 0;

	/* this is a bit of a pest: under some configurations, tramo/seats
	   outputs a series "irfin"; sometimes that is not created, but
	   we do get an "irreg".  So if we can't find the one, try looking
	   for the other.
	*/
	if (opt == TRAMO_SEATS && code == D13) { 
	    sprintf(sfname, "%s%cgraph%cseries%c%s", fname, SLASH, SLASH, SLASH,
		    tramo_series_strings[code + 1]);
	    fp = gretl_fopen(sfname, "r");
	    if (fp != NULL) {
		gotit = 1;
	    }
	    tramo_got_irfin = 0;
	}

	if (!gotit) {
	    fprintf(stderr, "Couldn't open %s\n", sfname);
	    sprintf(errmsg, _("Couldn't open %s"), sfname);
	    return 1;
	}
    }

    /* formulate name of new variable to add */
    if (opt == TRAMO_SEATS) {
	sprintf(varname, "%.5s_%.2s", pdinfo->varname[0], 
		tramo_series_strings[code]);
    } else {
	sprintf(varname, "%.4s_%s", pdinfo->varname[0], 
		x12a_series_strings[code]);
    }

    /* copy varname and label into place */
    strcpy(pdinfo->varname[v], varname);
    sprintf(VARLABEL(pdinfo, v), _(tx_descrip_formats[code]), pdinfo->varname[0]);

    if (opt == TRAMO_SEATS) {
	strcat(VARLABEL(pdinfo, v), " (TRAMO/SEATS)");
    } else {
	strcat(VARLABEL(pdinfo, v), " (X-12-ARIMA)");
    }	

    for (t=0; t<pdinfo->n; t++) {
	Z[v][t] = NADBL;
    }

    gretl_push_c_numeric_locale();

    if (opt == TRAMO_SEATS) {
	int i = 0;

	t = pdinfo->t1;
	while (fgets(line, 127, fp)) {
	    i++;
	    if (i >= 7 && sscanf(line, " %lf", &x) == 1) {
		if (t >= pdinfo->n) {
		    fprintf(stderr, "t = %d >= pdinfo->n = %d\n", t, pdinfo->n);
		    err = 1;
		    break;
		}		
		Z[v][t++] = x;
	    }
	}
    } else {
	/* grab the data from the x12arima file */
	while (fgets(line, 127, fp)) {
	    if (*line == 'd' || *line == '-') {
		continue;
	    }
	    if (sscanf(line, "%d %lf", &d, &x) != 2) {
		err = 1; 
		break;
	    }
	    yr = d / 100;
	    per = d % 100;
	    sprintf(date, "%d.%d", yr, per);
	    t = dateton(date, pdinfo);
	    if (t < 0 || t >= pdinfo->n) {
		err = 1;
		break;
	    }
	    Z[v][t] = x;
	}
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return err;
}

static void set_opts (tx_request *request)
{
    int i;

    request->savevars = 0;

    for (i=0; i<N_COMMON_OPTS; i++) {
	if (request->opt[i].check != NULL && 
	    GTK_TOGGLE_BUTTON(request->opt[i].check)->active) {
	    request->opt[i].save = 1;
	    if (i != TRIGRAPH) {
		request->savevars++;
	    }
	} else {
	    request->opt[i].save = 0;
	} 
    }
}

static void cancel_savevars (tx_request *request)
{
    int i;

    request->savevars = 0;

    for (i=0; i<N_COMMON_OPTS; i++) {
	request->opt[i].save = 0;
    } 
}

static int write_tramo_file (const char *fname, 
			     double **Z, DATAINFO *pdinfo,
			     int v, tx_request *request) 
{
    int startyr, startper, tsamp = pdinfo->t2 - pdinfo->t1 + 1; 
    char *p, tmp[8];
    double x;
    FILE *fp;
    int t;

    fp = gretl_fopen(fname, "w");
    if (fp == NULL) {
	return 1;
    }

    gretl_push_c_numeric_locale();

    x = date(pdinfo->t1, pdinfo->pd, pdinfo->sd0);
    startyr = (int) x;
    sprintf(tmp, "%g", x);
    p = strchr(tmp, '.');
    if (p != NULL) {
	startper = atoi(p + 1);
    } else {
	startper = 1;
    }

    fprintf(fp, "%s\n", pdinfo->varname[v]);
    fprintf(fp, "%d %d %d %d\n", tsamp, startyr, startper, pdinfo->pd);

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (t && t % pdinfo->pd == 0) fputc('\n', fp);
	if (na(Z[v][t])) {
	    fputs("-99999 ", fp);
	} else {
	    fprintf(fp, "%g ", Z[v][t]);
	}
    }
    fputc('\n', fp);

    if (print_tramo_options(request, fp) == 0) {
	request->code = TRAMO_ONLY; /* not running SEATS */
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return 0;
}

static int write_spc_file (const char *fname, 
			   double **Z, DATAINFO *pdinfo, 
			   int v, int *varlist) 
{
    int startyr, startper;
    char *p, tmp[8];
    double x;
    FILE *fp;
    int i, t;

    fp = gretl_fopen(fname, "w");
    if (fp == NULL) {
	return 1;
    }

    gretl_push_c_numeric_locale();

    x = date(pdinfo->t1, pdinfo->pd, pdinfo->sd0);
    startyr = (int) x;
    sprintf(tmp, "%g", x);
    p = strchr(tmp, '.');
    if (p != NULL) {
	startper = atoi(p + 1);
    } else {
	startper = 1;
    }

    fprintf(fp, "series{\n period=%d\n title=\"%s\"\n", pdinfo->pd, 
	    pdinfo->varname[v]);
    fprintf(fp, " start=%d.%d\n", startyr, startper);
    fputs(" data=(\n", fp);

    i = 0;
    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (na(Z[v][t])) {
	    fputs("-99999 ", fp); /* FIXME? */
	} else {
	    fprintf(fp, "%g ", Z[v][t]);
	}
	if ((i + 1) % 7 == 0) {
	    fputc('\n', fp);
	}
	i++;
    }
    fputs(" )\n}\n", fp);

    /* FIXME: make these values configurable */
    fputs("automdl{}\nx11{", fp);

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

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return 0;
}

static void form_varlist (int *varlist, tx_request *request)
{
    int i, j = 1;

    varlist[0] = 0;

    for (i=0; i<TRIGRAPH; i++) {
	if (request->opt[TRIGRAPH].save || request->opt[i].save) {
	    varlist[0] += 1;
	    varlist[j++] = i;
	}
    }
}

static void copy_basic_data_info (DATAINFO *targ, DATAINFO *src)
{
    targ->sd0 = src->sd0;
    strcpy(targ->stobs, src->stobs); 
    targ->t1 = src->t1;
    targ->t2 = src->t2;
    targ->pd = src->pd;
    targ->structure = src->structure;
}

static int save_vars_to_dataset (double ***pZ, DATAINFO *pdinfo,
				 double **tmpZ, DATAINFO *tmpinfo,
				 int *varlist, tx_request *request,
				 char *errmsg)
{
    int i, v, j, addvars = 0;

    /* how many vars are wanted, and new? */
    for (i=1; i<=varlist[0]; i++) {
	if (request->opt[varlist[i]].save && 
	    varindex(pdinfo, tmpinfo->varname[i]) == pdinfo->v) {
	    addvars++;
	}
    }

    if (addvars > 0 && dataset_add_series(addvars, pZ, pdinfo)) {
	strcpy(errmsg, _("Failed to allocate memory for new data"));
	return 1;
    }

    j = pdinfo->v - addvars;

    for (i=1; i<=varlist[0]; i++) {
	if (request->opt[varlist[i]].save) {
	    v = varindex(pdinfo, tmpinfo->varname[i]);
	    if (v < pdinfo->v) {
		copy_variable(*pZ, pdinfo, v, tmpZ, tmpinfo, i);
	    } else {
		copy_variable(*pZ, pdinfo, j++, tmpZ, tmpinfo, i);
	    }
	}
    }

    return 0;
}

#ifdef WIN32

static int helper_spawn (const char *prog, const char *vname,
			 const char *workdir, int code)
{
    char *cmd = NULL;
    int ret;

    if (code == TRAMO_ONLY) {
	cmd = g_strdup_printf("\"%s\" -i %s -k serie", prog, vname);
    } else if (code == TRAMO_SEATS) {
	cmd = g_strdup_printf("\"%s\" -OF %s", prog, vname);
    } else if (code == X12A) {
	cmd = g_strdup_printf("\"%s\" %s -r -p -q", prog, vname);
    } else {
	return 1;
    }

    if (cmd == NULL) {
	ret = E_ALLOC;
    } else {
	ret = winfork(cmd, workdir, SW_SHOWMINIMIZED, 
		      CREATE_NEW_CONSOLE | HIGH_PRIORITY_CLASS);
	g_free(cmd);
    }

    return ret;
}

#else

static int helper_spawn (const char *prog, const char *vname,
			 const char *workdir, int code)
{
    int err;

    if (code == TRAMO_ONLY) {
	err = glib_spawn(workdir, prog, "-i", vname, "-k", "serie", NULL);
    } else if (code == TRAMO_SEATS) {
	err = glib_spawn(workdir, prog, "-OF", vname, NULL);
    } else if (code == X12A) {
	err = glib_spawn(workdir, prog, vname, "-r", "-p", "-q", NULL);
    } else {
	err = 1;
    }

    return err;
}

#endif

int write_tx_data (char *fname, int varnum, 
		   double ***pZ, DATAINFO *pdinfo, int *graph, 
		   const char *prog, const char *workdir,
		   char *errmsg)
{
    int i, doit, err = 0;
    char varname[VNAMELEN];
    int varlist[4];
    FILE *fp = NULL;
    tx_request request;
    double **tmpZ;
    DATAINFO *tmpinfo;

    *errmsg = 0;

    /* sanity check */
    if (var_is_scalar(pdinfo, varnum)) {
	sprintf(errmsg, "%s %s", pdinfo->varname[varnum], 
		_("is a scalar"));
	return 1;
    }

    /* figure out which program we're using */
    if (strstr(prog, "tramo") != NULL) {
	request.code = TRAMO_SEATS;
    } else {
	request.code = X12A;
    }

    if (request.code == TRAMO_SEATS && (pdinfo->t2 - pdinfo->t1) > 599) {
	strcpy(errmsg, _("TRAMO can't handle more than 600 observations.\n"
			 "Please select a smaller sample."));
	return 1;
    } else if (request.code == X12A && (pdinfo->t2 - pdinfo->t1) > 719) {
	strcpy(errmsg, _("X-12-ARIMA can't handle more than 720 observations.\n"
			 "Please select a smaller sample."));
	return 1;
    }	

    request.pd = pdinfo->pd;

    /* show dialog and get option settings */
    doit = tx_dialog(&request); 
    if (!doit) {
	gtk_widget_destroy(request.dialog);
	return 0;
    }
    set_opts(&request);
    gtk_widget_destroy(request.dialog);

#if 0
    if (request.code == TRAMO_SEATS) {
	print_tramo_options(&request, stderr);
	return 1;
    }
#endif

    /* create little temporary dataset */
    tmpinfo = create_new_dataset(&tmpZ, 4, pdinfo->n, 0);
    if (tmpinfo == NULL) {
	return E_ALLOC;
    }

    copy_basic_data_info(tmpinfo, pdinfo);

    if (request.code == X12A) { 
	/* make a default x12a.mdl file if it doesn't already exist */
	sprintf(fname, "%s%cx12a.mdl", workdir, SLASH);
	fp = gretl_fopen(fname, "r");
	if (fp == NULL) {
	    fp = gretl_fopen(fname, "w");
	    if (fp == NULL) return 1;
	    fprintf(fp, "%s", default_mdl);
	    fclose(fp);
	} else {
	    fclose(fp);
	}
    } 

    sprintf(varname, pdinfo->varname[varnum]);
    form_varlist(varlist, &request);

    if (request.code == X12A) { 
	/* write out the .spc file for x12a */
	sprintf(fname, "%s%c%s.spc", workdir, SLASH, varname);
	write_spc_file(fname, *pZ, pdinfo, varnum, varlist);
    } else { 
	/* TRAMO, possibly plus SEATS */
	lower(varname);
	sprintf(fname, "%s%c%s", workdir, SLASH, varname);
	/* next line: this also sets request->code = TRAMO_ONLY if
	   seats is not to be run */
	write_tramo_file(fname, *pZ, pdinfo, varnum, &request);
	if (request.code == TRAMO_ONLY) {
	    cancel_savevars(&request); /* FIXME later */
	    varlist[0] = 0;
	}
    }

    /* now run the program(s) */

    if (request.code == X12A) {
	err = helper_spawn(prog, varname, workdir, X12A);
    } else { 
	char seats[MAXLEN];

	/* ensure any stale files get deleted first, just in case */
	clear_tramo_files(workdir, varname);

	err = helper_spawn(prog, varname, workdir, TRAMO_ONLY);

	if (!err && request.code == TRAMO_SEATS) {
	    get_seats_command(seats, prog);
	    err = helper_spawn(seats, varname, workdir, TRAMO_SEATS);
	}
    }

    if (!err) {
	if (request.code == X12A) {
	    sprintf(fname, "%s%c%s.out", workdir, SLASH, varname); 
	} else {
	    sprintf(fname, "%s%coutput%c%s.out", workdir, SLASH, SLASH, varname);
	} 

	/* save vars locally if needed; graph if wanted */
	if (varlist[0] > 0) {
	    copy_variable(tmpZ, tmpinfo, 0, *pZ, pdinfo, varnum);

	    for (i=1; i<=varlist[0]; i++) {
		err = add_series_from_file((request.code == X12A)? fname : workdir, 
					   varlist[i], tmpZ, tmpinfo, i, 
					   request.code, errmsg);
		if (err) {
		    fprintf(stderr, "add_series_from_file() failed\n");
		}
	    }

	    if (request.opt[TRIGRAPH].save) {
		err = graph_series(tmpZ, tmpinfo, request.code);
		if (err) {
		    fprintf(stderr, "graph_series() failed\n");
		} else {
		    *graph = 1;
		}
	    }
	}

	/* now save the local vars to main dataset, if wanted */
	if (request.savevars > 0) {
	    err = save_vars_to_dataset(pZ, pdinfo, tmpZ, tmpinfo, varlist, 
				       &request, errmsg);
	}
    }

    destroy_dataset(tmpZ, tmpinfo);

    return err;
}

