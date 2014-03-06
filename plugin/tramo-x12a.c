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
#include "version.h"
#include "estim_private.h"

#include <gtk/gtk.h>
#include "tramo_x12a.h"

#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 18)
# include "gtk_compat.h"
#endif

#define button_is_active(b) (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(b)))

#ifdef WIN32
# include <windows.h>
#endif

enum prog_codes {
    TRAMO_SEATS,
    TRAMO_ONLY,
    X12A
};

const char *x12a_save_strings[] = {
    "d11", /* seasonally adjusted */
    "d12", /* trend/cycle */
    "d13", /* irregular */
    NULL
};

static int tramo_got_irfin;

const char *tramo_save_strings[] = {
    "safin.t", /* final seasonally adjusted series */
    "trfin.t", /* final trend */
    "irfin.t", /* final irregular factor (component) */
    "irreg.t", /* irregular component (logs) */
    NULL
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

const char *get_tramo_save_string (int i)
{
    return tramo_save_strings[i];
}

#ifndef WIN32

#define SP_DEBUG 0

static int glib_spawn (const char *workdir, const char *fmt, ...)
{
    GError *gerr = NULL;
    gchar *sout = NULL;
    gchar *serr = NULL;
    gchar *argv[8];
    char *s;
    va_list ap;
    int i, ok, nargs;
    int status = 0;
    int err = 0;

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

    gretl_error_clear();

    ok = g_spawn_sync(workdir,
		      argv,
		      NULL,
		      G_SPAWN_SEARCH_PATH,
		      NULL,
		      NULL,
		      &sout,
		      &serr,
		      &status,
		      &gerr);

    if (!ok) {
	gretl_errmsg_set(gerr->message);
	fprintf(stderr, "spawn failed: '%s'\n", gerr->message);
	g_error_free(gerr);
	err = E_EXTERNAL;
    } else if (status != 0) {
	if (sout && *sout) {
	    gretl_errmsg_set(sout);
	    fprintf(stderr, "spawn: status = %d: '%s'\n", status, sout);
	} else {
	    gretl_errmsg_set(_("Command failed"));
	    fprintf(stderr, "spawn: status = %d\n", status);
	}
	err = E_DATA;
    } else if (serr && *serr) {
	fprintf(stderr, "stderr: '%s'\n", serr);
    }

    if (serr != NULL) g_free(serr);
    if (sout != NULL) g_free(sout);

    for (i=0; i<nargs; i++) {
	if (err) {
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

    return err;
}

#endif /* !WIN32 */

static void toggle_outliers (GtkToggleButton *b, tx_request *request)
{
    request->xopt.outliers = gtk_toggle_button_get_active(b);
}

static void toggle_trading_days (GtkToggleButton *b, tx_request *request)
{
    request->xopt.trdays = gtk_toggle_button_get_active(b);
}

static void set_logtrans (GtkButton *b, tx_request *request)
{
    gpointer p = g_object_get_data(G_OBJECT(b), "transval");

    request->xopt.logtrans = GPOINTER_TO_INT(p);
}

static void toggle_edit_script (GtkToggleButton *b, tx_request *request)
{
    GtkWidget **chk = g_object_get_data(G_OBJECT(b), "checks");
    gboolean s = gtk_toggle_button_get_active(b);
    int i;

    if (s) {
	*request->popt |= OPT_S;
    } else {
	*request->popt &= ~OPT_S;
    }

    for (i=0; i<4; i++) {
	gtk_widget_set_sensitive(chk[i], !s);
    }
}

void sensitize_tx_entry (GtkToggleButton *b, GtkWidget *w)
{
    gtk_widget_set_sensitive(w, button_is_active(b));
}

void update_tx_savename (GtkEntry *entry, char *name)
{
    strcpy(name, gtk_entry_get_text(entry));
}

static void add_x12a_options (tx_request *request, GtkBox *vbox)
{
    const gchar *save_strs[] = {
	N_("Seasonally adjusted series"),
	N_("Trend/cycle"),
	N_("Irregular")
    };
    gretlopt save_opts[] = {
	OPT_A, OPT_B, OPT_C
    };
    int save_codes[] = {
	TX_SA, TX_TR, TX_IR
    };
    GtkWidget *tmp, *b[3], *chk[4];
    GtkWidget *tbl;
    GSList *group;
    int i;

    tmp = gtk_label_new(_("X-12-ARIMA options"));
    gtk_box_pack_start(vbox, tmp, TRUE, TRUE, 5);

    tmp = gtk_check_button_new_with_label(_("Detect and correct for outliers"));
    gtk_box_pack_start(vbox, tmp, FALSE, FALSE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), request->xopt.outliers);
    g_signal_connect(GTK_TOGGLE_BUTTON(tmp), "toggled",
		     G_CALLBACK(toggle_outliers), request);

    tmp = gtk_check_button_new_with_label(_("Correct for trading days effect"));
    gtk_box_pack_start(vbox, tmp, FALSE, FALSE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), request->xopt.trdays);
    g_signal_connect(GTK_TOGGLE_BUTTON(tmp), "toggled",
		     G_CALLBACK(toggle_trading_days), request);

    tmp = gtk_hseparator_new();
    gtk_box_pack_start(vbox, tmp, FALSE, FALSE, 5);

    b[0] = gtk_radio_button_new_with_label(NULL, _("Log transformation"));
    gtk_box_pack_start(vbox, b[0], FALSE, FALSE, 0);
    g_signal_connect(GTK_TOGGLE_BUTTON(b[0]), "toggled",
		     G_CALLBACK(set_logtrans), request);
    g_object_set_data(G_OBJECT(b[0]), "transval", GINT_TO_POINTER(1));

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b[0]));
    b[1] = gtk_radio_button_new_with_label(group, _("No log transformation"));
    gtk_box_pack_start(vbox, b[1], FALSE, FALSE, 0);
    g_signal_connect(GTK_TOGGLE_BUTTON(b[1]), "toggled",
		     G_CALLBACK(set_logtrans), request);
    g_object_set_data(G_OBJECT(b[1]), "transval", GINT_TO_POINTER(2));

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b[1]));
    b[2] = gtk_radio_button_new_with_label(group, _("Automatic"));
    gtk_box_pack_start(vbox, b[2], FALSE, FALSE, 0);
    g_signal_connect(GTK_TOGGLE_BUTTON(b[2]), "toggled",
		     G_CALLBACK(set_logtrans), request);
    g_object_set_data(G_OBJECT(b[2]), "transval", GINT_TO_POINTER(3));

    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b[request->xopt.logtrans - 1]), 
				 TRUE); 

    tmp = gtk_hseparator_new();
    gtk_box_pack_start(vbox, tmp, FALSE, FALSE, 5);

    tmp = gtk_label_new(_("Save data"));
    gtk_box_pack_start(vbox, tmp, FALSE, FALSE, 0);

    tbl = gtk_table_new(3, 2, FALSE);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);

    for (i=0; i<3; i++) {
	/* buttons plus entries for saving series */
	gboolean active = (*request->popt & save_opts[i])? TRUE : FALSE;
	GtkWidget *entry;
	int idx = save_codes[i];

	chk[i] = gtk_check_button_new_with_label(_(save_strs[i]));
	request->opts[idx].check = chk[i];
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(chk[i]), active);
	gtk_table_attach_defaults(GTK_TABLE(tbl), chk[i], 0, 1, i, i+1);
	entry = gtk_entry_new();
	gtk_widget_set_sensitive(entry, active);
	gtk_entry_set_max_length(GTK_ENTRY(entry), VNAMELEN-1);
	gtk_entry_set_width_chars(GTK_ENTRY(entry), VNAMELEN-1);
	sprintf(request->opts[idx].savename, "%.8s_%s", request->yname, 
		x12a_save_strings[i]);
	gtk_entry_set_text(GTK_ENTRY(entry), request->opts[idx].savename);
	gtk_table_attach_defaults(GTK_TABLE(tbl), entry, 1, 2, i, i+1);
	g_object_set_data(G_OBJECT(chk[i]), "entry", entry);
	g_signal_connect(G_OBJECT(chk[i]), "toggled",
			 G_CALLBACK(sensitize_tx_entry), entry);
	g_signal_connect(G_OBJECT(GTK_EDITABLE(entry)), "changed",
			 G_CALLBACK(update_tx_savename), 
			 request->opts[idx].savename);
    }

    gtk_box_pack_start(vbox, tbl, FALSE, FALSE, 5);

    tmp = gtk_hseparator_new();
    gtk_box_pack_start(vbox, tmp, FALSE, FALSE, 5);
    
    chk[3] = gtk_check_button_new_with_label(_("Generate graph"));
    gtk_box_pack_start(vbox, chk[3], FALSE, FALSE, 0);
    request->opts[TRIGRAPH].check = chk[3];
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(chk[3]), 
				 (*request->popt & OPT_G)? TRUE : FALSE);

    tmp = gtk_check_button_new_with_label(_("Show full output"));
    gtk_box_pack_start(vbox, tmp, FALSE, FALSE, 0);
    request->opts[TEXTOUT].check = tmp;
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp),
				 (*request->popt & OPT_Q)? FALSE : TRUE);

    tmp = gtk_hseparator_new();
    gtk_box_pack_start(vbox, tmp, FALSE, FALSE, 5);

    b[0] = gtk_radio_button_new_with_label(NULL, _("Execute X-12-ARIMA directly"));
    gtk_box_pack_start(vbox, b[0], FALSE, FALSE, 0);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b[0]));
    b[1] = gtk_radio_button_new_with_label(group, _("Make X-12-ARIMA command file"));
    gtk_box_pack_start(vbox, b[1], FALSE, FALSE, 0);
    g_object_set_data(G_OBJECT(b[1]), "checks", chk);
    g_signal_connect(GTK_TOGGLE_BUTTON(b[1]), "toggled",
		     G_CALLBACK(toggle_edit_script), request);

    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b[0]), TRUE); 
}

static GtkWidget *x12a_help_button (GtkWidget *hbox, tx_request *request)
{
    GtkWidget *button;

    button = gtk_button_new_from_stock(GTK_STOCK_HELP);
    gtk_container_add(GTK_CONTAINER(hbox), button);
    gtk_button_box_set_child_secondary(GTK_BUTTON_BOX(hbox),
				       button, TRUE);
#if GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION >= 2
    gtk_button_box_set_child_non_homogeneous(GTK_BUTTON_BOX(hbox),
					     button, TRUE);
#endif
    g_signal_connect(G_OBJECT(button), "clicked", 
		     G_CALLBACK(request->helpfunc), 
		     NULL);
    
    return button;
}

static void tx_errbox (tx_request *request)
{
    GtkWidget *w = gtk_message_dialog_new(GTK_WINDOW(request->dialog),
					  GTK_DIALOG_DESTROY_WITH_PARENT,
					  GTK_MESSAGE_ERROR,
					  GTK_BUTTONS_CLOSE,
					  _("Expected a valid variable name"));

    gtk_dialog_run(GTK_DIALOG(w));
    gtk_widget_destroy(w);
}

static int check_savevars (tx_request *request)
{
    int i, err = 0;

    for (i=0; i<=TX_IR && !err; i++) {
	GtkWidget *w = request->opts[i].check;

	if (w != NULL && gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	    const char *vname = request->opts[i].savename;

	    if (*vname == '\0') {
		err = 1;
	    } else {
		err = check_varname(vname);
	    }
	    if (err) {
		GtkWidget *entry, *book;

		entry = g_object_get_data(G_OBJECT(w), "entry");
		book = g_object_get_data(G_OBJECT(entry), "book");
		if (book != NULL) {
		    gint pg = GPOINTER_TO_INT
			(g_object_get_data(G_OBJECT(entry), "output-page"));

		    gtk_notebook_set_current_page(GTK_NOTEBOOK(book), pg);
		}
		tx_errbox(request);
		gtk_widget_grab_focus(entry);
	    }
	}
    }

    return err;
}

static void tx_dialog_callback (GtkDialog *dlg, gint id, int *ret)
{
    int err = 0;

    if (id == GTK_RESPONSE_REJECT || id == GTK_RESPONSE_ACCEPT) {
	*ret = id;
    } else if (id == GTK_RESPONSE_DELETE_EVENT) {
	*ret = GTK_RESPONSE_REJECT;
    }

    if (*ret == GTK_RESPONSE_ACCEPT) {
	tx_request *request = g_object_get_data(G_OBJECT(dlg), "request");
	
	err = check_savevars(request);
    }

    if (!err) {
	gtk_main_quit();
    }
} 

static void nullify_request_dialog (GtkWidget *dlg, 
				    tx_request *request)
{
    request->dialog = NULL;
}

static int tx_dialog (tx_request *request, GtkWindow *parent)
{
    GtkDialogFlags dflags = GTK_DIALOG_DESTROY_WITH_PARENT;
    GtkWidget *hbox, *vbox;
    gint i, ret = 0;

    for (i=0; i<TX_MAXOPT; i++) {
	request->opts[i].check = NULL;
    }

    if (request->prog != X12A) {
	dflags |= GTK_DIALOG_MODAL;
    }

    request->dialog = 
	gtk_dialog_new_with_buttons((request->prog == TRAMO_SEATS)?
				    "TRAMO/SEATS" : "X-12-ARIMA",
				    parent,
				    dflags,
				    GTK_STOCK_CANCEL,
				    GTK_RESPONSE_REJECT,
				    GTK_STOCK_OK,
				    GTK_RESPONSE_ACCEPT,
				    NULL);

    g_signal_connect(G_OBJECT(request->dialog), "destroy",
		     G_CALLBACK(nullify_request_dialog),
		     request);
    g_object_set_data(G_OBJECT(request->dialog), "request",
		      request);

    vbox = gtk_vbox_new(FALSE, 0);    

    if (request->prog == TRAMO_SEATS) {
#if GTK_MAJOR_VERSION < 3    
	gtk_dialog_set_has_separator(GTK_DIALOG(request->dialog), FALSE);
#endif	
	add_tramo_options(request, vbox);
    } else {
	add_x12a_options(request, GTK_BOX(vbox));
    }

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 5);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(request->dialog));
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    if (request->prog == X12A) {
	hbox = gtk_dialog_get_action_area(GTK_DIALOG(request->dialog));
	x12a_help_button(hbox, request);
    }

    g_signal_connect(G_OBJECT(request->dialog), "response", 
		     G_CALLBACK(tx_dialog_callback), &ret);
    gtk_widget_show_all(request->dialog);
    gtk_main(); /* note: block */

    return (ret == GTK_RESPONSE_ACCEPT)? 1 : 0;
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

/* try to avoid collision of date and graph key, which
   defaults to right top
*/

static void set_keypos (const double *x, int t1, int t2,
			FILE *fp)
{
    int T = t2 - t1 + 1;

    if (T <= 12) {
	if (x[t2] > x[t1]) {
	    fputs("set key left top\n", fp);
	}
    } else {
	double m1, m2;
	int r = T / 6;

	m1 = gretl_mean(t1, t1 + r, x);
	m2 = gretl_mean(t2 - r, t2, x);

	if (m2 > m1) {
	    fputs("set key left top\n", fp);
	}
    }
}

static int graph_series (const DATASET *dset, tx_request *req)
{
    FILE *fp = NULL;
    const double *obs;
    int v_sa = TX_SA + 1;
    int v_tr = TX_TR + 1;
    int v_ir = TX_IR + 1;
    double irbar, irmax;
    int sub1 = 0;
    double f1;
    char title[32];
    int t, err = 0;

    obs = gretl_plotx(dset, OPT_NONE);
    if (obs == NULL) {
	return E_ALLOC;
    }

    fp = open_plot_input_file(PLOT_TRI_GRAPH, &err);
    if (err) {
	return err;
    }

    gretl_push_c_numeric_locale();

    if (dset->pd == 4) {
	if ((dset->t2 - dset->t1) / 4 < 8) {
	    fputs("set xtics nomirror 0,1\n", fp); 
	    fputs("set mxtics 4\n", fp);
	}
    } else if (dset->pd == 12) {
	if ((dset->t2 - dset->t1) / 12 < 8) {
	    fputs("set xtics nomirror 0,1\n", fp); 
	    fputs("set mxtics 12\n", fp);
	}
    }

    if (req->seasonal_ok) {
	f1 = 0.33;
	fputs("set size 1.0,1.0\nset multiplot\nset size 1.0,0.32\n", fp);
    } else {
	f1 = 0.5;
	fputs("set size 1.0,1.0\nset multiplot\nset size 1.0,0.48\n", fp);
	tramo_got_irfin = 0; /* I _think_ this may be right */
    }

    if (req->prog == TRAMO_SEATS && tramo_got_irfin) {
	/* need to divide by 100? */
	irmax = 10.0;
    } else {
	irmax = 0.5;
    }

    irbar = gretl_mean(dset->t1, dset->t2, dset->Z[v_ir]);
    if (irbar > irmax) {
	sub1 = 1;
    }

    /* irregular component */
    if (sub1) {
	sprintf(title, "%s - 1", _("irregular"));
    } else {
	sprintf(title, "%s", _("irregular"));
    }

    fprintf(fp, "set bars 0\n"
	    "set origin 0.0,0.0\n"
	    "set xzeroaxis\n"
	    "plot '-' using 1:%s title '%s' w impulses\n",
	    (sub1)? "($2-1.0)" : "2", title);

    for (t=dset->t1; t<=dset->t2; t++) {
	double yt = dset->Z[v_ir][t];

	if (req->prog == TRAMO_SEATS && tramo_got_irfin) {
	    yt /= 100.0;
	}

	fprintf(fp, "%.10g %.10g\n", obs[t], yt);
    }
    fputs("e\n", fp);

    set_keypos(dset->Z[0], dset->t1, dset->t2, fp);

    /* actual (in var 0) vs trend/cycle */

    fprintf(fp, "set origin 0.0,%.2f\n"
	    "plot '-' using 1:2 title '%s' w l, \\\n"
	    " '-' using 1:2 title '%s' w l\n",
	    f1, dset->varname[0], _("trend/cycle"));

    for (t=dset->t1; t<=dset->t2; t++) { 
	fprintf(fp, "%.10g %.10g\n", obs[t], dset->Z[0][t]);
    }
    fputs("e , \\\n", fp);

    for (t=dset->t1; t<=dset->t2; t++) { 
	fprintf(fp, "%.10g %.10g\n", obs[t], dset->Z[v_tr][t]);
    }
    fputs("e\n", fp);

    if (req->seasonal_ok) {
	/* actual vs seasonally adjusted */
	fprintf(fp, "set origin 0.0,0.66\n"
		"plot '-' using 1:2 title '%s' w l, \\\n"
		" '-' using 1:2 title '%s' w l\n",
		dset->varname[0], _("adjusted"));

	for (t=dset->t1; t<=dset->t2; t++) {
	    fprintf(fp, "%.10g %.10g\n", obs[t], dset->Z[0][t]);
	}
	fputs("e\n", fp);

	for (t=dset->t1; t<=dset->t2; t++) {
	    fprintf(fp, "%.10g %.10g\n", obs[t], dset->Z[v_sa][t]);
	}
	fputs("e\n", fp);
    }

    fputs("unset multiplot\n", fp);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

static void copy_variable (DATASET *targ, int targv,
			   DATASET *src, int srcv)
{
    int t;

    for (t=0; t<targ->n; t++) {
	targ->Z[targv][t] = src->Z[srcv][t];
    }

    strcpy(targ->varname[targv], src->varname[srcv]);
    series_set_label(targ, targv, series_get_label(src, srcv));
}

static void clear_tramo_files (const char *path, const char *vname)
{
    char fname[MAXLEN];
    int i;

    for (i=0; tramo_save_strings[i] != NULL; i++) {
	sprintf(fname, "%s%cgraph%cseries%c%s", path, SLASH, SLASH, SLASH,
		tramo_save_strings[i]);
	gretl_remove(fname);
    }

    sprintf(fname, "%s%coutput%c%s.out", path, SLASH, SLASH, vname);
    gretl_remove(fname);
}

static void clear_x12a_files (const char *path, const char *vname)
{
    char fname[MAXLEN];
    int i;

    for (i=0; x12a_save_strings[i] != NULL; i++) {
	sprintf(fname, "%s%c%s.%s", path, SLASH, vname,
		x12a_save_strings[i]);
	gretl_remove(fname);
    }    

    sprintf(fname, "%s%c%s.out", path, SLASH, vname);
    gretl_remove(fname);

    sprintf(fname, "%s%c%s.err", path, SLASH, vname);
    gretl_remove(fname);
}

static int add_series_from_file (const char *path, int src,
				 DATASET *dset, int targv, 
				 tx_request *request)
{
    FILE *fp;
    char line[128], sfname[MAXLEN];
    char varname[VNAMELEN], date[8];
    char label[MAXLABEL];
    double x;
    int d, yr, per, err = 0;
    int t;

    if (request->prog == TRAMO_SEATS) {
	tramo_got_irfin = 1;
	sprintf(sfname, "%s%cgraph%cseries%c%s", path, SLASH, SLASH, SLASH,
		tramo_save_strings[src]);
    } else {
	char *p;

	strcpy(sfname, path);
	p = strrchr(sfname, '.');
	if (p != NULL) {
	    strcpy(p + 1, x12a_save_strings[src]);
	}
    }

    fp = gretl_fopen(sfname, "r");

    if (fp == NULL) {
	/* couldn't open the file we wanted */
	int gotit = 0;

	/* This is a bit of a pest: under some configurations, tramo/seats
	   outputs a series "irfin"; sometimes that is not created, but
	   we do get an "irreg".  So if we can't find the one, try looking
	   for the other.  Also, the seasonally adjusted series "safin"
	   is not always available.
	*/
	if (request->prog == TRAMO_SEATS) {
	    if (src == TX_IR) { 
		/* try "irreg" */
		sprintf(sfname, "%s%cgraph%cseries%c%s", path, SLASH, SLASH, SLASH,
			tramo_save_strings[src + 1]);
		fp = gretl_fopen(sfname, "r");
		if (fp != NULL) {
		    gotit = 1;
		}
		tramo_got_irfin = 0;
	    } else if (src == TX_SA) {
		/* scrub all use of seasonal series */
		request->seasonal_ok = 0;
		if (request->opts[src].save) {
		    request->opts[src].save = 0;
		    request->savevars -= 1;
		}
		return 0;
	    }
	}

	if (!gotit) {
	    gretl_errmsg_sprintf(_("Couldn't open %s"), sfname);
	    return 1;
	}
    }

    /* formulate name of new variable to add */
    strcpy(varname, request->opts[src].savename);
    if (*varname == '\0') {
	if (request->prog == TRAMO_SEATS) {
	    sprintf(varname, "%.8s_%.2s", dset->varname[0], 
		    tramo_save_strings[src]);
	} else {
	    sprintf(varname, "%.8s_%s", dset->varname[0], 
		    x12a_save_strings[src]);
	}
    }

    /* copy varname and label into place */
    strcpy(dset->varname[targv], varname);
    sprintf(label, _(tx_descrip_formats[src]), dset->varname[0]);
    if (request->prog == TRAMO_SEATS) {
	strcat(label, " (TRAMO/SEATS)");
    } else {
	strcat(label, " (X-12-ARIMA)");
    }
    series_set_label(dset, targv, label);

    for (t=0; t<dset->n; t++) {
	dset->Z[targv][t] = NADBL;
    }

    gretl_push_c_numeric_locale();

    if (request->prog == TRAMO_SEATS) {
	int i = 0;

	t = dset->t1;
	while (fgets(line, 127, fp)) {
	    i++;
	    if (i >= 7 && sscanf(line, " %lf", &x) == 1) {
		if (t >= dset->n) {
		    fprintf(stderr, "t = %d >= dset->n = %d\n", t, dset->n);
		    err = 1;
		    break;
		}		
		dset->Z[targv][t++] = x;
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
	    t = dateton(date, dset);
	    if (t < 0 || t >= dset->n) {
		err = 1;
		break;
	    }
	    dset->Z[targv][t] = x;
	}
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return err;
}

static int grab_deseasonal_series (double *y, const DATASET *dset,
				   int prog, const char *path)
{
    FILE *fp;
    char line[128], sfname[MAXLEN], date[8];
    double yt;
    int d, yr, per, err = 0;
    int t;

    if (prog == TRAMO_SEATS) {
	sprintf(sfname, "%s%cgraph%cseries%c%s", path, SLASH, SLASH, SLASH,
		tramo_save_strings[TX_SA]);
    } else {
	char *p;

	strcpy(sfname, path);
	p = strrchr(sfname, '.');
	if (p != NULL) {
	    strcpy(p + 1, x12a_save_strings[TX_SA]);
	}
    }

    fp = gretl_fopen(sfname, "r");
    if (fp == NULL) {
	return E_FOPEN;
    }

    gretl_push_c_numeric_locale();

    if (prog == TRAMO_SEATS) {
	int i = 0;

	t = dset->t1;
	while (fgets(line, 127, fp)) {
	    i++;
	    if (i >= 7 && sscanf(line, " %lf", &yt) == 1) {
		if (t >= dset->n) {
		    fprintf(stderr, "t = %d >= dset->n = %d\n", t, dset->n);
		    err = E_DATA;
		    break;
		}		
		y[t++] = yt;
	    }
	}
    } else {
	/* grab the data from the x12arima file */
	while (fgets(line, 127, fp)) {
	    if (*line == 'd' || *line == '-') {
		continue;
	    }
	    if (sscanf(line, "%d %lf", &d, &yt) != 2) {
		err = 1; 
		break;
	    }
	    yr = d / 100;
	    per = d % 100;
	    sprintf(date, "%d.%d", yr, per);
	    t = dateton(date, dset);
	    if (t < 0 || t >= dset->n) {
		err = E_DATA;
		break;
	    }
	    y[t] = yt;
	}
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return err;
}

static void request_opts_init (tx_request *request, const DATASET *dset,
			       int varnum, void (*helpfunc))
{
    int i;

    request->savevars = 0;
    request->helpfunc = helpfunc;

    strcpy(request->yname, dset->varname[varnum]);

    request->xopt.logtrans = 3; /* x12a: automatic logs or not */
    request->xopt.outliers = 1; /* x12a: detect outliers */
    request->xopt.trdays = 0;   /* x12a: trading days correction */

    for (i=0; i<TX_MAXOPT; i++) {
	request->opts[i].save = 0;
	request->opts[i].savename[0] = '\0';
    }

    request->seasonal_ok = 1;
}

static void set_opts (tx_request *request)
{
    GtkWidget *w;
    int i;

    request->savevars = 0;

    *request->popt &= ~(OPT_A | OPT_B | OPT_C | OPT_G);

    for (i=0; i<TX_MAXOPT; i++) {
	w = request->opts[i].check;
	if (w != NULL && gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	    request->opts[i].save = 1;
	    if (i < TRIGRAPH) {
		request->savevars++;
		if (i == 0) {
		    *request->popt |= OPT_A;
		} else if (i == 1) {
		    *request->popt |= OPT_B;
		} else if (i == 2) {
		    *request->popt |= OPT_C;
		}
	    } else if (i == TRIGRAPH) {
		*request->popt |= OPT_G;
	    }
	} else {
	    request->opts[i].save = 0;
	} 
    }
}

static void cancel_savevars (tx_request *request)
{
    int i;

    request->savevars = 0;

    for (i=0; i<TX_MAXOPT; i++) {
	request->opts[i].save = 0;
    } 
}

static int write_tramo_file (const char *fname, 
			     const double *y, 
			     const char *vname, 
			     const DATASET *dset,
			     tx_request *request) 
{
    int startyr, startper;
    int T = dset->t2 - dset->t1 + 1; 
    char *p, tmp[8];
    double x;
    FILE *fp;
    int t;

    fp = gretl_fopen(fname, "w");
    if (fp == NULL) {
	return 1;
    }

    gretl_push_c_numeric_locale();

    x = date_as_double(dset->t1, dset->pd, dset->sd0);
    startyr = (int) x;
    sprintf(tmp, "%g", x);
    p = strchr(tmp, '.');
    if (p != NULL) {
	startper = atoi(p + 1);
    } else {
	startper = 1;
    }

    fprintf(fp, "%s\n", vname);
    fprintf(fp, "%d %d %d %d\n", T, startyr, startper, dset->pd);

    for (t=dset->t1; t<=dset->t2; t++) {
	if (t && t % dset->pd == 0) fputc('\n', fp);
	if (na(y[t])) {
	    fputs("-99999 ", fp);
	} else {
	    fprintf(fp, "%g ", y[t]);
	}
    }
    fputc('\n', fp);

    if (request == NULL) {
	fputs("$INPUT rsa=3,out=2,$END\n", fp);
    } else if (print_tramo_options(request, fp) == 0) {
	/* not running SEATS */
	request->prog = TRAMO_ONLY; 
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return 0;
}

static int x12_get_subperiod (double x, const DATASET *dset)
{
    int i, d = ceil(log10(dset->pd));
    int ret;

    x -= floor(x);
    for (i=0; i<d; i++) {
	x *= 10;
    }

    ret = (x-floor(x)) >.5 ? ceil(x) : floor(x);

    return ret;
}

static int write_spc_file (const char *fname, const double *y,
			   const char *vname,
			   const DATASET *dset, 
			   const int *savelist,
			   x12a_opts *xopt)
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

    x = date_as_double(dset->t1, dset->pd, dset->sd0);
    startyr = (int) x;
    sprintf(tmp, "%g", x);
    p = strchr(tmp, '.');
    if (p != NULL) {
	startper = x12_get_subperiod(x, dset);
    } else {
	startper = 1;
    }

    fprintf(fp, "series{\n period=%d\n title=\"%s\"\n", dset->pd, vname);
    fprintf(fp, " start=%d.%d\n", startyr, startper);

    for (t=dset->t1; t<=dset->t2; t++) {
	if (na(y[t])) {
	    fputs(" missingcode=-99999\n", fp);
	    break;
	}
    }

    fputs(" data=(\n", fp);

    i = 0;
    for (t=dset->t1; t<=dset->t2; t++) {
	if (na(y[t])) {
	    fputs("-99999 ", fp);
	} else {
	    fprintf(fp, "%g ", y[t]);
	}
	if ((i + 1) % 7 == 0) {
	    fputc('\n', fp);
	}
	i++;
    }
    fputs(" )\n}\n", fp);

    if (xopt->logtrans == 1) {
	fputs("transform{function=log}\n", fp);
    } else if (xopt->logtrans == 2) {
	fputs("transform{function=none}\n", fp);
    } else {
	fputs("transform{function=auto}\n", fp);
    }

    if (xopt->trdays) {
	fputs("regression{variables = td}\n", fp);
    }

    if (xopt->outliers) {
	fputs("outlier{}\n", fp);
    }

    fputs("automdl{}\n", fp); 

    fputs("x11{", fp);

    if (savelist[0] > 0) {
	if (savelist[0] == 1) {
	    fprintf(fp, " save=%s ", x12a_save_strings[savelist[1]]); 
	} else {
	    fputs(" save=( ", fp);
	    for (i=1; i<=savelist[0]; i++) {
		fprintf(fp, "%s ", x12a_save_strings[savelist[i]]);
	    }
	    fputs(") ", fp);
	}
    }

    fputs("}\n", fp);

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return 0;
}

static void form_savelist (int *list, tx_request *request)
{
    int i, j = 1;

    list[0] = 0;

    for (i=0; i<TRIGRAPH; i++) {
	if (request->opts[TRIGRAPH].save || request->opts[i].save) {
	    list[0] += 1;
	    list[j++] = i;
	}
    }
}

static void copy_basic_data_info (DATASET *targ, DATASET *src)
{
    targ->sd0 = src->sd0;
    strcpy(targ->stobs, src->stobs); 
    targ->t1 = src->t1;
    targ->t2 = src->t2;
    targ->pd = src->pd;
    targ->structure = src->structure;
}

static int save_vars_to_dataset (DATASET *dset,
				 DATASET *tmpset,
				 int *varlist, 
				 tx_request *request)
{
    int i, v, j, addvars = 0;

    /* how many vars are wanted, and how many are new? */
    for (i=1; i<=varlist[0]; i++) {
	if (request->opts[varlist[i]].save && 
	    series_index(dset, tmpset->varname[i]) == dset->v) {
	    addvars++;
	}
    }

    if (addvars > 0 && dataset_add_series(dset, addvars)) {
	return E_ALLOC;
    }

    j = dset->v - addvars;

    for (i=1; i<=varlist[0]; i++) {
	if (request->opts[varlist[i]].save) {
	    v = series_index(dset, tmpset->varname[i]);
	    if (v < dset->v) {
		copy_variable(dset, v, tmpset, i);
	    } else {
		copy_variable(dset, j++, tmpset, i);
	    }
	}
    }

    return 0;
}

#ifdef WIN32

static int helper_spawn (const char *path, const char *vname,
			 const char *workdir, int prog)
{
    char *cmd = NULL;
    int err = 0;

    if (prog == TRAMO_ONLY) {
	cmd = g_strdup_printf("\"%s\" -i %s -k serie", path, vname);
    } else if (prog == TRAMO_SEATS) {
	cmd = g_strdup_printf("\"%s\" -OF %s", path, vname);
    } else if (prog == X12A) {
	cmd = g_strdup_printf("\"%s\" %s -r -p -q", path, vname);
    } else {
	return E_EXTERNAL;
    }

    if (cmd == NULL) {
	err = E_ALLOC;
    } else {
	err = win_run_sync(cmd, workdir);
	g_free(cmd);
    }

    return err;
}

#else

static int helper_spawn (const char *path, const char *vname,
			 const char *workdir, int prog)
{
    int err;

    if (prog == TRAMO_ONLY) {
	err = glib_spawn(workdir, path, "-i", vname, "-k", "serie", NULL);
    } else if (prog == TRAMO_SEATS) {
	err = glib_spawn(workdir, path, "-OF", vname, NULL);
    } else if (prog == X12A) {
	err = glib_spawn(workdir, path, vname, "-r", "-p", "-q", NULL);
    } else {
	err = E_EXTERNAL;
    }

    return err;
}

#endif

/* The x12a .err file is always produced, but often contains no
   warning or error messages: in that case it contains 2 lines
   of boiler-plate text followed by two blank lines. Here we
   check if we got more than 4 lines of output.
*/

static int got_x12a_warning (const char *fname)
{
    FILE *fp = gretl_fopen(fname, "r");
    int ret = 0;

    if (fp != NULL) {
	char line[128];
	int n = 0;

	while (fgets(line, sizeof line, fp)) {
	    if (++n > 4 && !string_is_blank(line)) {
		ret = 1;
		break;
	    }
	}

	fclose(fp);
    }

    return ret;
}

/* make a default x12a.mdl file if it doesn't already exist */

static int check_x12a_model_file (const char *workdir)
{
    gchar *fname;
    FILE *fp;
    int err = 0;

    fname = g_strdup_printf("%s%cx12a.mdl", workdir, SLASH);
    fp = gretl_fopen(fname, "r");

    if (fp != NULL) {
	fclose(fp); /* assume we're OK */
    } else {
	fp = gretl_fopen(fname, "w");
	if (fp == NULL) {
	    err = E_FOPEN;
	} else {
	    fprintf(fp, "%s", default_mdl);
	    fclose(fp);
	}
    } 

    g_free(fname);

    return err;
}

static int check_sample_bound (int prog, const DATASET *dset)
{
    int T = dset->t2 - dset->t1 + 1;

    if (prog == TRAMO_SEATS && T > 600) {
	gretl_errmsg_set(_("TRAMO can't handle more than 600 observations.\n"
			 "Please select a smaller sample."));
	return E_EXTERNAL;
    } else if (prog == X12A) {
	int pdmax = get_x12a_maxpd();

	if (T > 50 * pdmax) {
	    gretl_errmsg_sprintf(_("X-12-ARIMA can't handle more than %d observations.\n"
				   "Please select a smaller sample."), 50 * pdmax);
	    return E_EXTERNAL;
	}
    }

    return 0;
}

/* Callback for the "Run" button in a gretl editor window displaying
   an x12a command file: run the commands (which are provided in @buf),
   and fill @outname with the name of the x12a output file.
*/

int exec_tx_script (char *outname, const gchar *buf)
{
    const char *exepath;
    const char *workdir = NULL;
    const char *tmpname = "x12atmp";
    gchar *tmppath;
    FILE *fp;
    int err = 0;

    *outname = '\0';

    exepath = gretl_x12_arima();
    workdir = gretl_x12_arima_dir();

    tmppath = g_strdup_printf("%s%c%s.spc", workdir, SLASH, tmpname);
    fp = gretl_fopen(tmppath, "w");
    if (fp == NULL) {
	g_free(tmppath);
	return E_FOPEN;
    }

    fputs(buf, fp);
    fclose(fp);
    g_free(tmppath);
	
    clear_x12a_files(workdir, tmpname);
    err = helper_spawn(exepath, tmpname, workdir, X12A);

    if (err == E_EXTERNAL) {
	; /* fatal: couldn't run program */
    } else if (err) {
	sprintf(outname, "%s%c%s.err", workdir, SLASH, tmpname);
    } else {
	/* set the output filename */
	sprintf(outname, "%s%c%s.out", workdir, SLASH, tmpname); 
    }

    return err;
}

/* Driver for access to TRAMO/SEATS (if @tramo is non-zero) or
   X-12-ARIMA, for analysis of the series with ID number @varnum. The
   @opt pointer is filled out based on selections from a dialog box.

   The @fname argument is filled out to tell the caller the name of
   the output file from the third-party program, if available; it will
   be an empty string if the user cancels, or if we fail to execute
   the program in question.

   At present the @help_func argument is used only for X-12-ARIMA,
   for which we offer a special option, namely writing a spec file
   for display in an editor window rather than directly executing
   x12a commands. If that option is selected we signal the caller
   by putting OPT_S into @opt, and we return the name of the gretl-
   generated command file in @fname.
*/

int write_tx_data (char *fname, int varnum, 
		   DATASET *dset, 
		   gretlopt *opt, int tramo,
		   GtkWindow *mainwin,
		   void (*help_func))
{
    const char *exepath;
    const char *workdir;
    char vname[VNAMELEN];
    int savelist[4];
    tx_request request;
    DATASET *tmpset = NULL;
    int savescript = 0;
    int i, doit;
    int err = 0;

    if (tramo) {
	request.prog = TRAMO_SEATS;
	exepath = gretl_tramo();
	workdir = gretl_tramo_dir();
    } else {
	request.prog = X12A;
	exepath = gretl_x12_arima();
	workdir = gretl_x12_arima_dir();
    }	
	
    request_opts_init(&request, dset, varnum, help_func);

    err = check_sample_bound(request.prog, dset);
    if (err) {
	return err;
    }

    request.pd = dset->pd;
    request.popt = opt;

    /* show dialog and get option settings */
    doit = tx_dialog(&request, mainwin); 
    if (doit) {
	set_opts(&request);
    } 
    if (request.dialog != NULL) {
	gtk_widget_destroy(request.dialog);
    }
    if (!doit) {
	*fname = '\0';
	return 0;
    }

#if 0
    if (request.prog == TRAMO_SEATS) {
	print_tramo_options(&request, stderr);
	return 1;
    }
#endif

    if (*opt & OPT_S) {
	savescript = 1;
    } else {
	/* create little temporary dataset */
	tmpset = create_auxiliary_dataset(4, dset->n, 0);
	if (tmpset == NULL) {
	    return E_ALLOC;
	}

	copy_basic_data_info(tmpset, dset);

	if (request.prog == X12A) { 
	    err = check_x12a_model_file(workdir);
	    if (err) {
		goto bailout;
	    }
	}
    } 

    strcpy(vname, dset->varname[varnum]);
    form_savelist(savelist, &request);

    if (request.prog == X12A) { 
	/* write out the .spc file for x12a */
	sprintf(fname, "%s%c%s.spc", workdir, SLASH, vname);
	write_spc_file(fname, dset->Z[varnum], vname, dset, savelist,
		       &request.xopt);
    } else { 
	/* TRAMO, possibly plus SEATS */
	gretl_lower(vname);
	gretl_trunc(vname, 8);
	sprintf(fname, "%s%c%s", workdir, SLASH, vname);
	/* next line: this also sets request->prog = TRAMO_ONLY if
	   SEATS is not to be run */
	write_tramo_file(fname, dset->Z[varnum], vname, dset, &request);
	if (request.prog == TRAMO_ONLY) {
	    cancel_savevars(&request); /* FIXME later */
	    savelist[0] = 0;
	}
    }

    if (savescript) {
	/* x12a file written; we're done */
	return 0;
    }

    /* now run the program(s): we try to ensure that any
       old output files get deleted first 
    */

    if (request.prog == X12A) {
	clear_x12a_files(workdir, vname);
	err = helper_spawn(exepath, vname, workdir, X12A);
    } else { 
	char seats[MAXLEN];

	clear_tramo_files(workdir, vname);
	err = helper_spawn(exepath, vname, workdir, TRAMO_ONLY);

	if (!err && request.prog == TRAMO_SEATS) {
	    get_seats_command(seats, exepath);
	    err = helper_spawn(seats, vname, workdir, TRAMO_SEATS);
	}
    }

    if (err == E_EXTERNAL) {
	/* fatal: couldn't run program */
	*fname = '\0';
    } else if (err) {
	if (request.prog == X12A) {
	    sprintf(fname, "%s%c%s.err", workdir, SLASH, vname);
	} else {
	    sprintf(fname, "%s%coutput%c%s.out", workdir, SLASH, SLASH, vname);
	}
    } else {
	if (request.prog == X12A) {
	    /* first check the .err file for warnings */
	    sprintf(fname, "%s%c%s.err", workdir, SLASH, vname);
	    if (got_x12a_warning(fname)) {
		*opt |= OPT_W;
	    }
	    /* then set the output filename */
	    sprintf(fname, "%s%c%s.out", workdir, SLASH, vname); 
	} else {
	    sprintf(fname, "%s%coutput%c%s.out", workdir, SLASH, SLASH, vname);
	    if (request.prog == TRAMO_ONLY) {
		/* no graph offered */
		request.opts[TRIGRAPH].save = 0;
		*opt |= OPT_T;
	    }
	} 

	/* save vars locally if needed; graph if wanted */
	if (savelist[0] > 0) {
	    const char *path = (request.prog == X12A)? fname : workdir;

	    copy_variable(tmpset, 0, dset, varnum);

	    for (i=1; i<=savelist[0]; i++) {
		err = add_series_from_file(path, savelist[i], tmpset,
					   i, &request);
		if (err) {
		    fprintf(stderr, "i = %d: add_series_from_file() failed\n", i);
		    if (request.prog == X12A) {
			/* switch to X12A error file */
			sprintf(fname, "%s%c%s.err", workdir, SLASH, vname);
		    }
		    break;
		} 
	    }

	    if (!err) {
		if (request.opts[TRIGRAPH].save) {
		    err = graph_series(tmpset, &request);
		    if (err) {
			fprintf(stderr, "graph_series() failed\n");
		    } else {
			*opt |= OPT_G;
		    }
		} else {
		    *opt &= ~OPT_G;
		}
	    }

	    if (request.prog == X12A) {
		if (request.opts[TEXTOUT].save) {
		    *opt &= ~OPT_Q;
		} else {
		    *opt |= OPT_Q;
		}
	    }
	}

	/* now save the local vars to main dataset, if wanted */
	if (!err && request.savevars > 0) {
	    err = save_vars_to_dataset(dset, tmpset, savelist, 
				       &request);
	}
    }

 bailout:

    destroy_dataset(tmpset);

    return err;
}

int adjust_series (const double *x, double *y, const DATASET *dset, 
		   int tramo, int use_log)
{
    int prog = (tramo)? TRAMO_SEATS : X12A;
    int savelist[2] = {1, TX_SA};
    const char *vname = "x";
    const char *exepath;
    const char *workdir;
    char fname[MAXLEN];
    int err = 0;

    if (prog == X12A) { 
	exepath = gretl_x12_arima();
	workdir = gretl_x12_arima_dir();
    } else { 
	exepath = gretl_tramo();
	workdir = gretl_tramo_dir();
    }    

    if (prog == X12A) { 
	err = check_x12a_model_file(workdir);
	if (err) {
	    return err;
	}
    } 

    if (prog == X12A) { 
	x12a_opts xopt = { 2, /* log transformation flag (2 == no) */
			   0, /* don't correct for outliers */
			   0  /* trading days correction */
	};

	if (use_log) {
	    xopt.logtrans = 1;
	}
	if (dset->pd == 12) {
	    xopt.trdays = 1;
	}
	sprintf(fname, "%s%c%s.spc", workdir, SLASH, vname);
	write_spc_file(fname, x, vname, dset, savelist, &xopt);
    } else { 
	sprintf(fname, "%s%c%s", workdir, SLASH, vname);
	write_tramo_file(fname, x, vname, dset, NULL); 
    }

    if (prog == X12A) {
	clear_x12a_files(workdir, vname);
	err = helper_spawn(exepath, vname, workdir, X12A);
    } else { 
	char seats[MAXLEN];

	clear_tramo_files(workdir, vname);
	err = helper_spawn(exepath, vname, workdir, TRAMO_ONLY);
	if (!err) {
	    get_seats_command(seats, exepath);
	    err = helper_spawn(seats, vname, workdir, TRAMO_SEATS);
	}
    }

    if (!err) {
	const char *path = (prog == X12A)? fname : workdir;

	err = grab_deseasonal_series(y, dset, prog, path);
    }

    return err;
}
