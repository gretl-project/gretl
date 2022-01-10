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

/* TRAMO/SEATS, X-13ARIMA plugin for gretl */

#include "libgretl.h"
#include "version.h"
#include "estim_private.h"

#ifdef WIN32
# include "gretl_win32.h"
#endif

#include <gtk/gtk.h>
#include "tramo_x12a.h"

#define DSDEBUG 0

#define button_is_active(b) (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(b)))

#ifdef WIN32
# include <windows.h>
#endif

static void display_x13a_output (char *fname, int err, PRN *prn);

enum prog_codes {
    TRAMO_SEATS,
    TRAMO_ONLY,
    X13A
};

const char *x11_save_strings[] = {
    "d11", /* seasonally adjusted */
    "d12", /* trend/cycle */
    "d13", /* irregular */
    NULL
};

const char *x13_seats_save_strings[] = {
    "s11", /* seasonally adjusted */
    "s12", /* trend/cycle */
    "s13", /* irregular */
    NULL
};

static int tramo_got_irfin;

const char *tramo_save_strings[] = {
    "safin.t", /* final seasonally adjusted series */
    "trfin.t", /* final trend */
    "irfin.t", /* final irregular factor (component) */
    "xlin.t",  /* linearized series */
    NULL
};

const char *tx_descrip_formats[] = {
    N_("seasonally adjusted %s"),
    N_("trend/cycle for %s"),
    N_("irregular component of %s"),
    N_("linearized %s")
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

static void toggle_seats (GtkToggleButton *b, tx_request *request)
{
    request->xopt.seats = gtk_toggle_button_get_active(b);
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

static void add_x13a_options (tx_request *request, GtkBox *vbox)
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
    GtkWidget *tmp, *hb, *b[3], *chk[4];
    GtkWidget *rx, *rs;
    GtkWidget *tbl;
    GSList *group;
    int i;

    tmp = gtk_label_new(_("X-13ARIMA options"));
    gtk_box_pack_start(vbox, tmp, TRUE, TRUE, 5);

    hb = gtk_hbox_new(FALSE, 0);
    tmp = gtk_label_new(_("Seasonal adjustment algrithm:"));
    gtk_box_pack_start(GTK_BOX(hb), tmp, 0, 0, 0);
    gtk_box_pack_start(vbox, hb, FALSE, FALSE, 0);

    /* radio buttons for X11 vs SEATS */
    rx = gtk_radio_button_new_with_label(NULL, "X11");
    gtk_box_pack_start(vbox, rx, FALSE, FALSE, 0);
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(rx));
    rs = gtk_radio_button_new_with_label(group, "SEATS");
    gtk_box_pack_start(vbox, rs, FALSE, FALSE, 0);
    g_signal_connect(GTK_TOGGLE_BUTTON(rs), "toggled",
		     G_CALLBACK(toggle_seats), request);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(rx), TRUE);

    tmp = gtk_hseparator_new();
    gtk_box_pack_start(vbox, tmp, FALSE, FALSE, 5);

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

    hb = gtk_hbox_new(FALSE, 0);
    tmp = gtk_label_new(_("Save data"));
    gtk_box_pack_start(GTK_BOX(hb), tmp, 0, 0, 0);
    gtk_box_pack_start(vbox, hb, FALSE, FALSE, 0);

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
		x11_save_strings[i]);
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

    b[0] = gtk_radio_button_new_with_label(NULL, _("Execute X-13ARIMA directly"));
    gtk_box_pack_start(vbox, b[0], FALSE, FALSE, 0);

    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b[0]));
    b[1] = gtk_radio_button_new_with_label(group, _("Make X-13ARIMA command file"));
    gtk_box_pack_start(vbox, b[1], FALSE, FALSE, 0);
    g_object_set_data(G_OBJECT(b[1]), "checks", chk);
    g_signal_connect(GTK_TOGGLE_BUTTON(b[1]), "toggled",
		     G_CALLBACK(toggle_edit_script), request);

    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b[0]), TRUE);
}

static GtkWidget *x13a_help_button (GtkWidget *hbox, tx_request *request)
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
    int imax = request->prog == X13A ? TX_IR : TX_LN;
    int i, err = 0;

    for (i=0; i<=imax && !err; i++) {
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

    if (request->prog != X13A) {
	dflags |= GTK_DIALOG_MODAL;
    }

    request->dialog =
	gtk_dialog_new_with_buttons((request->prog == TRAMO_SEATS)?
				    "TRAMO/SEATS" : "X-13ARIMA",
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
	add_x13a_options(request, GTK_BOX(vbox));
    }

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), vbox, FALSE, FALSE, 5);

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(request->dialog));
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    if (request->prog == X13A) {
	hbox = gtk_dialog_get_action_area(GTK_DIALOG(request->dialog));
	x13a_help_button(hbox, request);
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
    p = strrslash(seats);
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
    gchar *title;
    int t, err = 0;

    obs = gretl_plotx(dset, OPT_NONE);
    if (obs == NULL) {
	return E_ALLOC;
    }

    fp = open_plot_input_file(PLOT_TRI_GRAPH, 0, &err);
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
	title = g_strdup_printf("%s - 1", _("irregular"));
    } else {
	title = g_strdup(_("irregular"));
    }

    fprintf(fp, "set bars 0\n"
	    "set origin 0.0,0.0\n"
	    "set xzeroaxis\n"
	    "plot '-' using 1:%s title '%s' w impulses\n",
	    (sub1)? "($2-1.0)" : "2", title);
    g_free(title);

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
    const char *vlabel;
    int t;

    for (t=0; t<targ->n; t++) {
	targ->Z[targv][t] = src->Z[srcv][t];
    }

    strcpy(targ->varname[targv], src->varname[srcv]);
    vlabel = series_get_label(src, srcv);
    if (vlabel != NULL && *vlabel != '\0') {
	series_set_label(targ, targv, vlabel);
    }
}

static void clear_tramo_files (const char *path, const char *vname)
{
    char fname[MAXLEN];
    int i;

    for (i=0; tramo_save_strings[i] != NULL; i++) {
	gretl_build_path(fname, path, "graph", "series",
			     tramo_save_strings[i], NULL);
	gretl_remove(fname);
    }

    /* just in case, clear "irreg" too */
    gretl_build_path(fname, path, "graph", "series", "irreg.t", NULL);
    gretl_remove(fname);

    /* clear @vname.out */
    gretl_build_path(fname, path, "output", vname, NULL);
    strcat(fname, ".out");
    gretl_remove(fname);

    /* clear the generic "summary.txt" */
    gretl_build_path(fname, path, "output", "summary.txt", NULL);
    gretl_remove(fname);
}

static void clear_x13a_files (const char *path, const char *vname)
{
    char fname[MAXLEN];
    int i;

    gretl_build_path(fname, path, vname, NULL);

    for (i=0; x11_save_strings[i] != NULL; i++) {
	switch_ext_in_place(fname, x11_save_strings[i]);
	gretl_remove(fname);
    }
    for (i=0; x13_seats_save_strings[i] != NULL; i++) {
	switch_ext_in_place(fname, x13_seats_save_strings[i]);
	gretl_remove(fname);
    }

    switch_ext_in_place(fname, "out");
    gretl_remove(fname);

    switch_ext_in_place(fname, "err");
    gretl_remove(fname);
}

/* Peek into output/x.out and see if it contains evidence
   that SEATS worked OK but did not find any evidence of
   seasonality, in which case return 1. If we can't open
   x.out, or can't find such evidence, return 0.
*/

static int seats_no_seasonal (const char *path)
{
    char *p, line[MAXLEN];
    char outname[MAXLEN];
    FILE *fp;
    int ret = 0;

    gretl_build_path(outname, path, "output", "x.out", NULL);

    fp = gretl_fopen(outname, "r");
    if (fp == NULL) {
	return 0;
    }

    while (fgets(line, sizeof line, fp)) {
	if (strstr(line, "NO DECOMPOSITION IS PERFORMED")) {
	    ret = 1;
	    break;
	} else if ((p = strstr(line, "SEASONAL")) != NULL) {
	    ret = string_is_blank(p + 8);
	    break;
	}
    }

    fclose(fp);

    return ret;
}

static const char *addstr (tx_request *request)
{
    if (request->prog == X13A) {
	return "(X-13ARIMA)";
    } else if (request->prog == TRAMO_SEATS) {
	return "(TRAMO/SEATS)";
    } else {
	return "(TRAMO)";
    }
}

static int add_series_from_file (const char *path, int src,
				 DATASET *dset, int targv,
				 tx_request *request)
{
    FILE *fp;
    char line[128], sfname[MAXLEN];
    char varname[VNAMELEN], date[16];
    gchar *label, *tmp;
    double x;
    int d, yr, per, err = 0;
    int t;

    if (request->prog == X13A) {
	const char **save_strings;
	int seats = request->xopt.seats;
	char *p;

	save_strings = seats ? x13_seats_save_strings : x11_save_strings;
	strcpy(sfname, path);
	p = strrchr(sfname, '.');
	if (p != NULL) {
	    strcpy(p + 1, save_strings[src]);
	}
    } else {
	tramo_got_irfin = 1;
	gretl_build_path(sfname, path, "graph", "series",
			     tramo_save_strings[src], NULL);
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
	if (request->prog == TRAMO_SEATS || request->prog == TRAMO_ONLY) {
	    if (src == TX_IR) {
		/* try "irreg" */
		gretl_build_path(sfname, path, "graph", "series",
				     "irreg.t", NULL);
		fp = gretl_fopen(sfname, "r");
		if (fp != NULL) {
		    gotit = 1;
		}
		tramo_got_irfin = 0;
	    } else if (src == TX_LN) {
		/* maybe no linearization was required? */
		gretl_build_path(sfname, path, "graph", "series",
				 "xorigt.t", NULL);
		fp = gretl_fopen(sfname, "r");
		if (fp != NULL) {
		    gotit = 1;
		}
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
	if (request->prog == X13A) {
	    sprintf(varname, "%.8s_%s", dset->varname[0],
		    x11_save_strings[src]);
	} else {
	    sprintf(varname, "%.8s_%.2s", dset->varname[0],
		    tramo_save_strings[src]);
	}
    }

    /* copy varname and label into place */
    strcpy(dset->varname[targv], varname);
    tmp = g_strdup_printf(_(tx_descrip_formats[src]), dset->varname[0]);
    label = g_strdup_printf("%s %s", tmp, addstr(request));
    series_set_label(dset, targv, label);
    g_free(label);
    g_free(tmp);

    for (t=0; t<dset->n; t++) {
	dset->Z[targv][t] = NADBL;
    }

    gretl_push_c_numeric_locale();

    if (request->prog == TRAMO_SEATS || request->prog == TRAMO_ONLY) {
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
	/* grab the data from the x13as file */
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

static int grab_adjusted_series (double *y, const double *x,
				 const DATASET *dset,
				 int prog, x13a_opts *xopt,
				 const char *path, PRN *prn)
{
    FILE *fp;
    char line[128], sfname[MAXLEN], date[16];
    double yt;
    int no_seas = 0;
    int t, d, yr, per;
    int err = 0;

    if (prog == TRAMO_SEATS) {
	gretl_build_path(sfname, path, "graph", "series",
			 tramo_save_strings[TX_SA], NULL);
    } else {
	/* x13as: @path should end with ".spc" */
	const char **save_strings;
	char *p;

	save_strings = xopt->seats ? x13_seats_save_strings : x11_save_strings;
	strcpy(sfname, path);
	p = strrchr(sfname, '.');
	if (p != NULL) {
	    strcpy(p + 1, save_strings[xopt->output]);
	}
    }

    fp = gretl_fopen(sfname, "r");

    if (fp == NULL) {
	err = E_FOPEN;
	if (prog == TRAMO_SEATS) {
	    if (seats_no_seasonal(path)) {
		gretl_warnmsg_set(_("no seasonality was detected"));
		no_seas = 1;
		err = 0;
	    }
	} else if (prn != NULL) {
	    display_x13a_output(sfname, 1, prn);
	}
	if (err) {
	    return err;
	}
    }

    gretl_push_c_numeric_locale();

    if (no_seas) {
	/* give back the original series? */
	for (t=dset->t1; t<=dset->t2; t++) {
	    y[t] = x[t];
	}
    } else if (prog == TRAMO_SEATS) {
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
	/* grab the data from the x13as file */
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

    if (fp != NULL) {
	fclose(fp);
    }

    return err;
}

static int grab_linearized_series (double *y, const double *x,
				   const DATASET *dset,
				   const char *path,
				   const char *vname)
{
    FILE *fp;
    char line[128], sfname[MAXLEN];
    double yt;
    int err = 0;
    int i, t;

    gretl_build_path(sfname, path, "graph", "series",
			 tramo_save_strings[TX_LN], NULL);
    fp = gretl_fopen(sfname, "r");

    /* The linearized series may not have been produced: is this
       because an error has occurred, or simply because TRAMO
       judges that the series doesn't need linearizing?
    */

    if (fp == NULL) {
	gretl_build_path(sfname, path, "output", vname, NULL);
	strcat(sfname, ".out");
	fp = gretl_fopen(sfname, "r");
	if (fp != NULL) {
	    fclose(fp); /* OK ? */
	    gretl_build_path(sfname, path, "output", "summary.txt", NULL);
	    fp = gretl_fopen(sfname, "r");
	    if (fp != NULL) {
		fclose(fp); /* again, OK ? */
		/* is use of xorigt.t correct? */
		gretl_build_path(sfname, path, "graph", "series",
				     "xorigt.t", NULL);
		fp = gretl_fopen(sfname, "r");
	    }
	}
    }

    if (fp == NULL) {
	return E_FOPEN;
    }

    gretl_push_c_numeric_locale();

    i = 0;
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

    request->xopt.logtrans = 3; /* x13a: automatic logs or not */
    request->xopt.outliers = 1; /* x13a: detect outliers */
    request->xopt.trdays = 0;   /* x13a: trading days correction */
    request->xopt.wdays = 0;    /* x13a: working days correction */
    request->xopt.easter = 0;   /* x13a: Easter effect */
    request->xopt.seats = 0;    /* x13a: use SEATS rather than X11 */
    request->xopt.airline = 0;  /* x13a: force "airline" ARIMA spec */
    request->xopt.critical = NADBL; /* for use with outliers */

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

    *request->popt &= ~(OPT_A | OPT_B | OPT_C | OPT_D | OPT_G);

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
		} else if (i == 3) {
		    *request->popt |= OPT_D;
		}
	    } else if (i == TRIGRAPH) {
		*request->popt |= OPT_G;
	    }
	} else {
	    request->opts[i].save = 0;
	}
    }
}

static void adjust_savevars (tx_request *request,
			     int *savelist)
{
    int i;

    request->savevars = 0;

    for (i=0; i<TX_MAXOPT; i++) {
	if (i == TX_LN && request->opts[i].save) {
	    request->savevars = 1;
	} else {
	    request->opts[i].save = 0;
	}
    }

    if (request->savevars == 1) {
	savelist[0] = 1;
	savelist[1] = TX_LN;
    } else {
	savelist[0] = 0;
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
    char *p, tmp[16];
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
	if (na(y[t])) {
	    fputs("-99999\n", fp);
	} else {
	    fprintf(fp, "%.12g\n", y[t]);
	}
    }

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

static int x13_get_subperiod (double x, const DATASET *dset)
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

static int write_spc_file (const char *fname,
			   const double *y,
			   const char *vname,
			   const DATASET *dset,
			   const int *savelist,
			   x13a_opts *xopt)
{
    const char **save_strings;
    int startyr, startper;
    char *p, tmp[16];
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
	startper = x13_get_subperiod(x, dset);
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
	    fprintf(fp, "%.12g ", y[t]);
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
	if (xopt->easter) {
	    if (xopt->trdays == 2) {
		fprintf(fp, "regression{aictest = (td easter)}\n");
	    } else if (xopt->trdays) {
		fputs("regression{variables = (td easter[8])}\n", fp);
	    }
	} else {
	    if (xopt->trdays == 2) {
		fprintf(fp, "regression{aictest = (%s)}\n", "td");
	    } else if (xopt->trdays) {
		fputs("regression{variables = td}\n", fp);
	    }
	}
    } else if (xopt->wdays) {
    	if (xopt->easter) {
	    if (xopt->wdays == 2) {
		fprintf(fp, "regression{aictest = (td1coef easter)}\n");
	    } else if (xopt->wdays) {
		fputs("regression{variables = (td1coef easter[8])}\n", fp);
	    }
	} else {
	    if (xopt->wdays == 2) {
		fprintf(fp, "regression{aictest = (%s)}\n", "td1coef");
	    } else if (xopt->wdays) {
		fputs("regression{variables = td1coef}\n", fp);
	    }
	}
    }

    if (xopt->outliers) {
	if (!na(xopt->critical)) {
	    fprintf(fp, "outlier{critical = %g}\n", xopt->critical);
	} else {
	    fputs("outlier{}\n", fp);
	}
    }
    if (xopt->airline) {
	fputs("arima {model=(0,1,1)(0,1,1)}\n", fp);
    } else {
	fputs("automdl{}\n", fp);
    }

#if DSDEBUG
    fprintf(stderr, "x13as: using %s\n", xopt->seats ? "seats" : "x11");
#endif

    if (xopt->seats) {
	save_strings = x13_seats_save_strings;
	fputs("seats{", fp);
    } else {
	save_strings = x11_save_strings;
	fputs("x11{", fp);
    }

    if (savelist[0] > 0) {
	if (savelist[0] == 1) {
	    fprintf(fp, " save=%s ", save_strings[savelist[1]]);
	} else {
	    fputs(" save=( ", fp);
	    for (i=1; i<=savelist[0]; i++) {
		fprintf(fp, "%s ", save_strings[savelist[i]]);
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
    int imax = request->prog == X13A ? TX_IR : TX_LN;
    int i, j = 1;

    list[0] = 0;

    for (i=0; i<=imax; i++) {
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

static int helper_spawn (const char *path,
			 const char *vname,
			 const char *workdir,
			 int prog)
{
    char *cmd = NULL;
    int err = 0;

    if (prog == TRAMO_ONLY) {
	cmd = g_strdup_printf("\"%s\" -i %s -k serie", path, vname);
    } else if (prog == TRAMO_SEATS) {
	cmd = g_strdup_printf("\"%s\" -OF %s", path, vname);
    } else if (prog == X13A) {
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

static int helper_spawn (const char *path,
			 const char *vname,
			 const char *workdir,
			 int prog)
{
    int err;

    if (prog == TRAMO_ONLY) {
	err = glib_spawn(workdir, path, "-i", vname, "-k", "serie", NULL);
    } else if (prog == TRAMO_SEATS) {
	err = glib_spawn(workdir, path, "-OF", vname, NULL);
    } else if (prog == X13A) {
	err = glib_spawn(workdir, path, vname, "-r", "-p", "-q", NULL);
    } else {
	err = E_EXTERNAL;
    }

    return err;
}

#endif

/* The x13a .err file is always produced, but often contains no
   warning or error messages: in that case it contains 2 lines
   of boiler-plate text followed by two blank lines. Here we
   check if we got more than 4 lines of output.
*/

static int got_x13a_warning (const char *fname, int *err)
{
    FILE *fp = gretl_fopen(fname, "r");
    int ret = 0;

    if (fp != NULL) {
	char line[128];
	int n = 0;

	while (fgets(line, sizeof line, fp)) {
	    if (strstr(line, "ERROR:") != NULL) {
		*err = ret = 1;
	    }
	    if (++n > 4 && !string_is_blank(line)) {
		ret = 1;
	    }
	}
	fclose(fp);
    }

    return ret;
}

/* make a default x13a.mdl file if it doesn't already exist */

static int check_x13a_model_file (const char *workdir)
{
    gchar *fname;
    FILE *fp;
    int err = 0;

    fname = g_strdup_printf("%s%cx13a.mdl", workdir, SLASH);
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
    } else if (prog == X13A) {
	int pdmax = get_x13as_maxpd();

	if (T > 50 * pdmax) {
	    gretl_errmsg_sprintf(_("X-13ARIMA can't handle more than %d observations.\n"
				   "Please select a smaller sample."), 50 * pdmax);
	    return E_EXTERNAL;
	}
    }

    return 0;
}

/* Callback for the "Run" button in a gretl editor window displaying
   an x13a command file: run the commands (which are provided in @buf),
   and fill @outname with the name of the x13a output file.
*/

int exec_tx_script (char *outname, const gchar *buf)
{
    const char *exepath;
    const char *workdir = NULL;
    const char *tmpname = "x13atmp";
    FILE *fp;
    int err = 0;

    *outname = '\0';
    exepath = gretl_x12_arima();
    workdir = gretl_x12_arima_dir();

    gretl_build_path(outname, workdir, tmpname, NULL);
    strcat(outname, ".spc");
    fp = gretl_fopen(outname, "w");
    *outname = '\0';

    if (fp == NULL) {
	return E_FOPEN;
    }

    fputs(buf, fp);
    fclose(fp);

    clear_x13a_files(workdir, tmpname);
    err = helper_spawn(exepath, tmpname, workdir, X13A);

    if (err == E_EXTERNAL) {
	; /* fatal: couldn't run program */
    } else if (err) {
	gretl_build_path(outname, workdir, tmpname, NULL);
	strcat(outname, ".err");
    } else {
	/* set the output filename */
	gretl_build_path(outname, workdir, tmpname, NULL);
	strcat(outname, ".out");
    }

    return err;
}

/* Driver for access to TRAMO/SEATS (if @tramo is non-zero) or
   X-13ARIMA, for analysis of the series with ID number @varnum. The
   @opt pointer is filled out based on selections from a dialog box.

   The @fname argument is filled out to tell the caller the name of
   the output file from the third-party program, if available; it will
   be an empty string if the user cancels, or if we fail to execute
   the program in question.

   At present the @help_func argument is used only for X-13ARIMA,
   for which we offer a special option, namely writing a spec file
   for display in an editor window rather than directly executing
   x13a commands. If that option is selected we signal the caller
   by putting OPT_S into @opt, and we return the name of the gretl-
   generated command file in @fname.
*/

int write_tx_data (char *fname,
		   int varnum,
		   DATASET *dset,
		   gretlopt *opt, int tramo,
		   int *warning,
		   GtkWindow *mainwin,
		   void (*help_func))
{
    const char *exepath;
    const char *workdir;
    char vname[VNAMELEN];
    int savelist[5];
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
	request.prog = X13A;
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
	tmpset = create_auxiliary_dataset(5, dset->n, 0);
	if (tmpset == NULL) {
	    return E_ALLOC;
	}

	copy_basic_data_info(tmpset, dset);

	if (request.prog == X13A) {
	    err = check_x13a_model_file(workdir);
	    if (err) {
		goto bailout;
	    }
	}
    }

    strcpy(vname, dset->varname[varnum]);
    form_savelist(savelist, &request);

    if (request.prog == X13A) {
	/* write out the .spc file for x13a */
	gretl_build_path(fname, workdir, vname, NULL);
	strcat(fname, ".spc");
	write_spc_file(fname, dset->Z[varnum], vname, dset, savelist,
		       &request.xopt);
    } else {
	/* TRAMO, possibly plus SEATS */
	gretl_lower(vname);
	gretl_trunc(vname, 8);
	gretl_build_path(fname, workdir, vname, NULL);
	/* next line: this also sets request->prog = TRAMO_ONLY if
	   SEATS is not to be run */
	write_tramo_file(fname, dset->Z[varnum], vname, dset, &request);
	if (request.prog == TRAMO_ONLY) {
	    adjust_savevars(&request, savelist);
	}
    }

    if (savescript) {
	/* x13a file written; we're done */
	return 0;
    }

    /* now run the program(s): we try to ensure that any
       old output files get deleted first
    */

    if (request.prog == X13A) {
	clear_x13a_files(workdir, vname);
	err = helper_spawn(exepath, vname, workdir, X13A);
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
	goto bailout;
    }

    if (request.prog == X13A) {
	/* see if we got a warning -- and if so, whether it
	   should count as an error */
	gretl_build_path(fname, workdir, vname, NULL);
	strcat(fname, ".err");
	*warning = got_x13a_warning(fname, &err);
	if (!err) {
	    /* switch @fname to the .out file */
	    switch_ext_in_place(fname, "out");
	}
    } else {
	gretl_build_path(fname, workdir, "output", vname, NULL);
	strcat(fname, ".out");
	if (request.prog == TRAMO_ONLY) {
	    /* no graph offered */
	    request.opts[TRIGRAPH].save = 0;
	    *opt |= OPT_T;
	}
    }

    if (err) {
	goto bailout;
    }

    /* save vars locally if needed; graph if wanted */
    if (savelist[0] > 0) {
	const char *path = request.prog == X13A ? fname : workdir;

	copy_variable(tmpset, 0, dset, varnum);

	for (i=1; i<=savelist[0]; i++) {
	    err = add_series_from_file(path, savelist[i], tmpset,
				       i, &request);
	    if (err) {
		fprintf(stderr, "i = %d: add_series_from_file() failed\n", i);
		if (request.prog == X13A) {
		    /* switch @fname to point to X13A error file */
		    switch_ext_in_place(fname, "err");
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

	if (request.prog == X13A) {
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

 bailout:

    destroy_dataset(tmpset);

    return err;
}

static int parse_deseas_bundle (x13a_opts *xopt, gretl_bundle *b,
				PRN *prn)
{
    const char *trival_strs[] = {
	"no", "yes", "auto"
    };
    int lt = 2; /* log transformation */
    int td = 2; /* trading days */
    int wd = 0; /* working days */
    int got_td_spec = 0;
    int err = 0;

    xopt->outliers = gretl_bundle_get_bool(b, "outlier_correction", 0);
    xopt->seats    = gretl_bundle_get_bool(b, "seats", 0);
    xopt->airline  = gretl_bundle_get_bool(b, "airline", 0);
    xopt->easter   = gretl_bundle_get_bool(b, "easter", 0);

    if (gretl_bundle_has_key(b, "logtrans")) {
	lt = gretl_bundle_get_int(b, "logtrans", &err);
	if (!err) {
	    if (lt == 0) {
		xopt->logtrans = 2; /* no log */
	    } else if (lt == 1) {
		xopt->logtrans = 1; /* force log */
	    } else if (lt == 2) {
		xopt->logtrans = 3; /* automatic */
	    } else {
		err = E_INVARG;
	    }
	}
    }

    if (gretl_bundle_has_key(b, "trading_days")) {
	got_td_spec = 1;
	td = gretl_bundle_get_int(b, "trading_days", &err);
	if (!err && (td < 0 || td > 2)) {
	    err = E_INVARG;
	}
	if (!err) {
	    xopt->trdays = td;
	}
    }

    if (gretl_bundle_has_key(b, "working_days") && (!td || !got_td_spec)) {
	/* cannot kick in if trading days explicitly set to non-zero */
    	wd = gretl_bundle_get_int(b, "working_days", &err);
    	if (!err && (wd < 0 || wd > 2)) {
	    err = E_INVARG;
	}
	if (!err) {
	    /* working days and trading days should not coexist */
	    xopt->wdays = wd;
	    xopt->trdays = 0;
	}
    }

    if (gretl_bundle_has_key(b, "critical")) {
	double crit = gretl_bundle_get_scalar(b, "critical", &err);

	if (!err && (crit < 2 || crit > 10)) {
	    err = E_INVARG;
	}
	if (!err) {
	    xopt->critical = crit;
	}
    }

    if (gretl_bundle_has_key(b, "output")) {
	int output = gretl_bundle_get_int(b, "output", &err);

	if (!err && (output < 1 || output > 3)) {
	    err = E_INVARG;
	}
	if (!err) {
	    /* zero-based */
	    xopt->output = output - 1;
	}
    }

    if (!err) {
	xopt->verbose = gretl_bundle_get_int(b, "verbose", NULL);
    }

    if (xopt->verbose > 0) {
	const char **ostrs = xopt->seats ? x13_seats_save_strings :
	    x11_save_strings;

	/* FIXME translations */

	pprintf(prn, "x13as options:\n");
	pprintf(prn, "  adjustment algorithm:    %s\n", xopt->seats ? "SEATS" : "X11");
	if (xopt->outliers) {
	    if (na(xopt->critical)) {
		pprintf(prn, "  outlier correction:      %s\n", "yes");
	    } else {
		pprintf(prn, "  outlier correction:      %s (critical value %g)\n",
			"yes", xopt->critical);
	    }
	} else {
	    pprintf(prn, "  outlier correction:      %s\n", "no");
	}
	pprintf(prn, "  trading days correction: %s\n", trival_strs[td]);
	pprintf(prn, "  working days correction: %s\n", trival_strs[wd]);
	pprintf(prn, "  easter effect:           %s\n", xopt->easter ? "yes" : "no");
	pprintf(prn, "  log transformation:      %s\n", trival_strs[lt]);
	pprintf(prn, "  force 'airline' model:   %s\n", xopt->airline ? "yes" : "no");
	pprintf(prn, "  output series:           %s\n", ostrs[xopt->output]);
	pputc(prn, '\n');
    }

    return err;
}

static void display_x13a_output (char *fname, int err, PRN *prn)
{
    FILE *fp;

    if (err) {
	switch_ext_in_place(fname, "err");
    } else {
	switch_ext_in_place(fname, "out");
    }

    fp = fopen(fname, "r");

    if (fp != NULL) {
	char line[1024];
	gchar *tmp;

	while (fgets(line, sizeof line, fp)) {
	    if (g_utf8_validate(line, -1, NULL)) {
		pputs(prn, line);
	    } else {
		tmp = g_convert(line, -1, "UTF-8", "ISO-8859-1",
				NULL, NULL, NULL);
		if (tmp != NULL) {
		    pputs(prn, tmp);
		    g_free(tmp);
		}
	    }
	}
	fclose(fp);
    }
}

/* implements the deseas() function */

int adjust_series (const double *x, double *y,
		   const char *vname, const DATASET *dset,
		   int tramo, gretl_bundle *opts,
		   PRN *prn)
{
    int prog = (tramo)? TRAMO_SEATS : X13A;
    int savelist[2] = {1, TX_SA};
    x13a_opts xopt = {3, 0, 0, 0, 0, 0, 0, 0, 0, NADBL};
    const char *exepath;
    const char *workdir;
    char fname[MAXLEN];
    int err = 0;

#if DSDEBUG
    fprintf(stderr, "deseas: using %s\n", prog == X13A ? "x13as" : "tramo");
#endif

    if (vname == NULL) {
	vname = "x";
    }

    if (prog == X13A) {
	exepath = gretl_x12_arima();
	workdir = gretl_x12_arima_dir();
    } else {
	exepath = gretl_tramo();
	workdir = gretl_tramo_dir();
    }

    if (prog == X13A) {
	err = check_x13a_model_file(workdir);
	if (err) {
	    return err;
	}
    }

    if (prog == X13A) {
	if (opts != NULL) {
	    err = parse_deseas_bundle(&xopt, opts, prn);
	    savelist[1] = xopt.output;
	}
	if (!err) {
	    gretl_build_path(fname, workdir, vname, NULL);
	    strcat(fname, ".spc");
	    write_spc_file(fname, x, vname, dset, savelist, &xopt);
	}
    } else {
	gretl_build_path(fname, workdir, vname, NULL);
	write_tramo_file(fname, x, vname, dset, NULL);
    }

    if (!err && prog == X13A) {
	clear_x13a_files(workdir, vname);
	err = helper_spawn(exepath, vname, workdir, X13A);
    } else if (!err) {
	char seats[MAXLEN];

	clear_tramo_files(workdir, vname);
	err = helper_spawn(exepath, vname, workdir, TRAMO_ONLY);
	if (!err) {
	    get_seats_command(seats, exepath);
	    err = helper_spawn(seats, vname, workdir, TRAMO_SEATS);
	}
    }

    if (!err) {
	const char *path = (prog == X13A)? fname : workdir;

	err = grab_adjusted_series(y, x, dset, prog, &xopt, path, prn);
    }

    if (!err && xopt.verbose > 1) {
	display_x13a_output(fname, 0, prn);
    }

    return err;
}

/* implements the linearize() function */

int linearize_series (const double *x, double *y, const DATASET *dset)
{
    const char *vname = "x";
    const char *exepath;
    const char *workdir;
    char fname[MAXLEN];
    int err = 0;

    exepath = gretl_tramo();
    workdir = gretl_tramo_dir();

    gretl_build_path(fname, workdir, vname, NULL);
    write_tramo_file(fname, x, vname, dset, NULL);
    clear_tramo_files(workdir, vname);
    err = helper_spawn(exepath, vname, workdir, TRAMO_ONLY);

    if (!err) {
	err = grab_linearized_series(y, x, dset, workdir, vname);
    }

    return err;
}
