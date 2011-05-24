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

#include "libgretl.h"
#include <gtk/gtk.h>
#include "tramo_x12a.h"

#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 18)
# include "gtk_compat.h"
#endif

#define button_is_active(b) (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(b)))

typedef struct _tramo_options tramo_options;

struct _tramo_options {
    int rsa;      /* for standard auto setting, rsa = 3 (guide.pdf) */

    int iatip;    /* correction for outliers? 0=no, 1=auto */
    int aio;      /* sorts of outliers recognized (codes 1, 2, or 3;
		     0 is also acceptable in case only tramo is run) */
    float va;     /* critical value for outliers (or 0.0 for auto) */

    /* widgets for dealing with outliers */
    GtkWidget *iatip_button;
    GtkWidget *aio_transitory_button;
    GtkWidget *aio_shift_button;
    GtkWidget *aio_innov_button;
    GtkWidget *va_button;
    GtkWidget *va_spinner;
    GtkWidget *aio_label;
    GtkWidget *va_label;

    int lam;      /* log transformation? (0=logs, 1=levels, -1=auto) */
    int imean;    /* mean correction? (0=no, 1=yes) */

    int inic;     /* controls automatic model identification */
    int idif;     /* ditto */

    int auto_arima;  /* 1 = ARIMA spec will be left up to TRAMO; 
			0 for manual specification */
    int d;        /* number of non-seasonal differences */
    int bd;       /* number of seasonal differences */
    int p;        /* number of non-seasonal AR terms */
    int bp;       /* number of seasonal AR terms */
    int q;        /* number of non-seasonal MA terms */
    int bq;       /* number of seasonal MA terms */

    /* widgets for dealing with ARIMA terms */
    GtkWidget *d_spin;
    GtkWidget *bd_spin;
    GtkWidget *p_spin;
    GtkWidget *bp_spin;
    GtkWidget *q_spin;
    GtkWidget *bq_spin;

    int mq;       /* periodicity (12 = monthly; acceptable values are
		     1, 2, 3, 4, 5, 6 and 12 */
    int noadmiss; /* when model "does not accept an admissible decomposition",
		     use an approximation (1) or not (0) */
    int seats;    /* 0: no SEATS input file created 
		     1: create SEATS file; estimation re-done in SEATS
		     2: create SEATS file; estimation not re-done in SEATS
		  */
    int out;      /* verbosity level: 
		     0 = detailed output file per series
		     1 = reduced output file per series
		     2 = very brief summary
		     3 = no output file
		  */

    tx_request *request; /* generic tramo/x-12-arima request struct */
};

#define option_widgets_shown(p)  (p->va_spinner != NULL)

static void tramo_options_set_defaults (tramo_options *opts, int pd)
{
    opts->rsa = 3;           /* standard auto analysis */
    opts->iatip = 1;         /* detect outliers */
    opts->aio = 2;           /* both transitory changes and level shifts */
    opts->va = 0.0;          /* let critical value be decided by tramo */
    opts->lam = -1;          /* leave log/level decision to tramo */
    opts->imean = 1;         /* mean correction */
    opts->auto_arima = 1;    /* leave ARIMA spec to tramo */
    opts->inic = 3;
    opts->idif = 3;
    opts->d = opts->bd = 1;  
    opts->p = opts->bp = 0;
    opts->q = opts->bq = 1;
    opts->mq = pd;
    opts->noadmiss = 1;      /* use approximation if needed */
    opts->seats = (pd > 1);  /* make SEATS file and re-do estimation? */
    opts->out = 0;           /* verbose */
}

static void 
tramo_custom_tabs_set_sensitive (GtkWidget *notebook, gboolean s)
{
    gint i;
    GtkWidget *page;

    for (i=2; i<5; i++) {
	page = gtk_notebook_get_nth_page(GTK_NOTEBOOK(notebook), i);
	gtk_widget_set_sensitive(page, s);
    }
}

static void main_auto_callback (GtkWidget *w, GtkWidget *notebook)
{
    tramo_options *opts = g_object_get_data(G_OBJECT(notebook), "opts");

    if (w == NULL || button_is_active(w)) {
	tramo_custom_tabs_set_sensitive(notebook, FALSE);
	opts->rsa = 3;
    } else {
	tramo_custom_tabs_set_sensitive(notebook, TRUE);
	opts->rsa = 0;
    }
}

static void va_spinner_set_state (tramo_options *opts)
{
    if (!option_widgets_shown(opts)) return;
    
    gtk_widget_set_sensitive(opts->va_spinner, 
			     gtk_widget_is_sensitive(opts->va_label) &&
			     !button_is_active(opts->va_button));
}

static void outlier_options_set_sensitive (tramo_options *opts, gboolean s)
{
    if (!option_widgets_shown(opts)) return;

    gtk_widget_set_sensitive(opts->aio_label, s);
    gtk_widget_set_sensitive(opts->aio_transitory_button, s);
    gtk_widget_set_sensitive(opts->aio_shift_button, s);
    gtk_widget_set_sensitive(opts->aio_innov_button, s && opts->seats == 0);
    gtk_widget_set_sensitive(opts->va_label, s);
    gtk_widget_set_sensitive(opts->va_button, s);
    va_spinner_set_state(opts);
}

static void flip_iatip (GtkWidget *w, tramo_options *opts)
{
    if (!option_widgets_shown(opts)) return;
    
    if (button_is_active(w)) {
	outlier_options_set_sensitive(opts, TRUE);
	opts->iatip = 1;
    } else {
	outlier_options_set_sensitive(opts, FALSE);
	opts->iatip = 0;
    }	
}

static void flip_auto_va (GtkWidget *w, tramo_options *opts)
{
    if (!option_widgets_shown(opts)) return;  
  
    if (button_is_active(w)) {
	gtk_widget_set_sensitive(opts->va_spinner, FALSE);
	opts->va = 0.0;
    } else {
	gtk_widget_set_sensitive(opts->va_spinner, TRUE);
    }
}

static void arima_options_set_sensitive (tramo_options *opts, gboolean s)
{
    gtk_widget_set_sensitive(opts->d_spin, s);
    gtk_widget_set_sensitive(opts->p_spin, s);
    gtk_widget_set_sensitive(opts->q_spin, s);

    if (opts->request->pd > 1) {
	gtk_widget_set_sensitive(opts->bd_spin, s);
	gtk_widget_set_sensitive(opts->bp_spin, s);
	gtk_widget_set_sensitive(opts->bq_spin, s);
    }
}

static void flip_auto_arima (GtkWidget *w, tramo_options *opts)
{
    if (!option_widgets_shown(opts)) return;  
  
    if (button_is_active(w)) {
	arima_options_set_sensitive(opts, FALSE);
	opts->auto_arima = 1;
    } else {
	arima_options_set_sensitive(opts, TRUE);
	opts->auto_arima = 0;
    }
}

static void tramo_innov_callback (GtkWidget *w, tramo_options *opts)
{
    if (button_is_active(w)) {
	gtk_toggle_button_set_active
	    (GTK_TOGGLE_BUTTON(opts->aio_transitory_button), TRUE);
	gtk_toggle_button_set_active
	    (GTK_TOGGLE_BUTTON(opts->aio_shift_button), TRUE);
	opts->aio = 0;
	opts->seats = 0;
    }
}

static void tramo_aio_callback (GtkWidget *w, tramo_options *opts)
{
    GtkWidget *other_button, *this_button = w;

    if (!option_widgets_shown(opts)) return;  

    if (this_button == opts->aio_transitory_button) {
	other_button = opts->aio_shift_button;
    } else {
	other_button = opts->aio_transitory_button;
    }

    /* must have at least one box checked */
    if (!button_is_active(this_button) && !button_is_active(other_button)) {
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(other_button), TRUE);
    }

    if (button_is_active(opts->aio_transitory_button)) {
	if (button_is_active(opts->aio_shift_button)) {
	    /* both transitory and level-shifts */
	    opts->aio = 2;
	} else {
	    /* only transitory */
	    opts->aio = 1;
	}
    } else {
	/* only level-shifts */
	opts->aio = 3;
    }
}

static void 
seats_specific_widgets_set_sensitive (tramo_options *opts,
				      gboolean s)
{
    tx_request *request = opts->request;
    int i;

    if (opts->aio_innov_button != NULL) {
	gtk_widget_set_sensitive(opts->aio_innov_button, !s); 
    }
    
    for (i=0; i<TX_MAXOPT; i++) {
	if (request->opts[i].check != NULL) {
	    gtk_widget_set_sensitive(request->opts[i].check, s);
	}
    }
}				      

static void real_set_seats (tramo_options *opts, gint val)
{
    if (val > 0) {
	seats_specific_widgets_set_sensitive(opts, TRUE);
    } else {
	seats_specific_widgets_set_sensitive(opts, FALSE);
    }

    opts->seats = val;
}

static void set_seats (GtkWidget *w, tramo_options *opts)
{
    real_set_seats(opts, 1);
}

static void set_no_seats (GtkWidget *w, tramo_options *opts)
{
    real_set_seats(opts, 0);
}


static void set_out (GtkWidget *w, tramo_options *opts)
{
    opts->out = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), 
						  "out_value"));
}

static void set_lam (GtkWidget *w, tramo_options *opts)
{
    opts->lam = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), 
						  "lam_value"));
}

static void set_imean (GtkWidget *w, tramo_options *opts)
{
    opts->imean = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), 
						    "imean_value"));
}

static GtkWidget *make_notebook_page_table (GtkWidget *notebook, 
					    const gchar *tab_title,
					    gint rows, gint cols,
					    gint *page)
{
    GtkWidget *box, *tmp, *tbl;
    gint pg;
 
    box = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(box), 10);
    gtk_widget_show(box);

    tmp = gtk_label_new(tab_title);
    gtk_widget_show(tmp);

    pg = gtk_notebook_append_page(GTK_NOTEBOOK(notebook), box, tmp);
    if (page != NULL) {
	*page = pg;
    }

    tbl = gtk_table_new(rows, cols, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(box), tbl, FALSE, FALSE, 0);
    gtk_widget_show(tbl); 

    return tbl;
}

static void tramo_tab_general (GtkWidget *notebook, tx_request *request)
{
    GtkWidget *tbl, *tmp;
    int tbl_len = 4, row = 0;
    GSList *group = NULL;
    const char *radio_labels[] = {
	N_("Time-series model plus seasonal adjustment"),
	N_("Time-series model only")
    };
   
    tbl = make_notebook_page_table(notebook, _("General"), tbl_len, 2, NULL);

    /* checkbox for standard default run */
    tmp = gtk_check_button_new_with_label(_("Standard automatic analysis"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 1, row, row + 1);
    row++;
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), TRUE);
    g_object_set_data(G_OBJECT(notebook), "opts", request->gui);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(main_auto_callback), 
		     notebook);

    /* horizontal separator */
    tmp = gtk_hseparator_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;

    /* TRAMO + SEATS option */
    tmp = gtk_radio_button_new_with_label(NULL, _(radio_labels[0]));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), (request->pd > 1));
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(tmp));
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;

    if (request->pd == 1) {
	gtk_widget_set_sensitive(tmp, FALSE);
    } else {    
	g_signal_connect(G_OBJECT(tmp), "clicked",
			 G_CALLBACK(set_seats), 
			 request->gui);
    }

    /* TRAMO-only option */
    tmp = gtk_radio_button_new_with_label(group, _(radio_labels[1]));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), (request->pd == 1));
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(tmp));
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;

    if (request->pd == 1) {
	gtk_widget_set_sensitive(tmp, FALSE);
    } else {   
	g_signal_connect(G_OBJECT(tmp), "clicked",
			 G_CALLBACK(set_no_seats), 
			 request->gui);
    }
}

static void tramo_tab_output (GtkWidget *notebook, tx_request *request)
{
    const gchar *save_strs[] = {
	N_("Seasonally adjusted series"),
	N_("Trend/cycle"),
	N_("Irregular")
    };
    int save_codes[] = {
	TX_SA, TX_TR, TX_IR
    };
    GtkWidget *tbl, *tmp;
    int tbl_len = 10, row = 0;
    GSList *group = NULL;
    int i, pg;

    if (request->pd == 1) {
	tbl_len -= 2;
    }

    tbl = make_notebook_page_table(notebook, _("Output"), tbl_len, 2, &pg);

    /* label for output window detail */
    tmp = gtk_label_new(_("Output window:"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;

    /* full detail option */
    tmp = gtk_radio_button_new_with_label(NULL, _("Full details"));
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), TRUE);    
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(tmp));
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    g_signal_connect(G_OBJECT(tmp), "clicked", G_CALLBACK(set_out), 
		     request->gui);
    g_object_set_data(G_OBJECT(tmp), "out_value", GINT_TO_POINTER(0));
    row++;

    /* reduced output option */
    tmp = gtk_radio_button_new_with_label(group, _("Reduced output"));
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(tmp));
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    g_signal_connect(G_OBJECT(tmp), "clicked", G_CALLBACK(set_out), 
		     request->gui);
    g_object_set_data(G_OBJECT(tmp), "out_value", GINT_TO_POINTER(1));
    row++;

    tmp = gtk_hseparator_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;

    /* label pertaining to saving series */
    tmp = gtk_label_new(_("Save to data set:"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;

    for (i=0; i<3; i++) {
	/* buttons plus entries for saving series */
	gboolean active = FALSE;
	GtkWidget *entry;
	int idx = save_codes[i];

	if (request->pd == 1 && i == 0) {
	    request->opts[idx].check = NULL;
	    continue;
	}

	tmp = gtk_check_button_new_with_label(_(save_strs[i]));
	request->opts[idx].check = tmp;
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), active);
	gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 1, row, row+1);
	entry = gtk_entry_new();
	gtk_widget_set_sensitive(entry, active);
	gtk_entry_set_max_length(GTK_ENTRY(entry), VNAMELEN-1);
	gtk_entry_set_width_chars(GTK_ENTRY(entry), VNAMELEN-1);
	g_object_set_data(G_OBJECT(tmp), "entry", entry);
	g_object_set_data(G_OBJECT(entry), "book", notebook);
	g_object_set_data(G_OBJECT(entry), "output-page", 
			  GINT_TO_POINTER(pg));
	sprintf(request->opts[idx].savename, "%.8s_%.2s", request->yname,
		get_tramo_save_string(i));
	gtk_entry_set_text(GTK_ENTRY(entry), request->opts[idx].savename);
	gtk_table_attach_defaults(GTK_TABLE(tbl), entry, 1, 2, row, row+1);
	g_signal_connect(G_OBJECT(tmp), "toggled",
			 G_CALLBACK(sensitize_tx_entry), entry);
	g_signal_connect(G_OBJECT(GTK_EDITABLE(entry)), "changed",
			 G_CALLBACK(update_tx_savename), 
			 request->opts[idx].savename);
	row++;
    }

    tmp = gtk_hseparator_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 1, row, row + 1);
    row++;

    if (request->pd > 1) {
	tmp = gtk_check_button_new_with_label(_("Generate graph"));
	gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 1, row, row + 1);
	request->opts[TRIGRAPH].check = tmp;
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), TRUE);
    } else {
	request->opts[TRIGRAPH].check = NULL;
    }
}

static void tramo_tab_transform (GtkWidget *notebook, tramo_options *opts)
{
    GtkWidget *b1, *b2, *b3, *tbl, *tmp;
    GSList *log_group = NULL, *mean_group = NULL;
    int tbl_len = 6, row = 0;

    tbl = make_notebook_page_table(notebook, _("Transformations"), 
				   tbl_len, 2, NULL);

    /* logs versus levels: radio group */

    /* logs option */
    b1 = gtk_radio_button_new_with_label(NULL, _("Log transformation"));
    log_group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    gtk_table_attach_defaults(GTK_TABLE(tbl), b1, 0, 2, row, row + 1);
    g_signal_connect(G_OBJECT(b1), "clicked", G_CALLBACK(set_lam), opts);
    g_object_set_data(G_OBJECT(b1), "lam_value", GINT_TO_POINTER(0));
    row++;

    /* no logs option */
    b2 = gtk_radio_button_new_with_label(log_group, _("No log transformation"));
    log_group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b2));
    gtk_table_attach_defaults(GTK_TABLE(tbl), b2, 0, 2, row, row + 1);
    g_signal_connect(G_OBJECT(b2), "clicked", G_CALLBACK(set_lam), opts);
    g_object_set_data(G_OBJECT(b2), "lam_value", GINT_TO_POINTER(1));
    row++;

    /* automatic log/level option */
    b3 = gtk_radio_button_new_with_label(log_group, _("Automatic"));
    log_group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b3));
    gtk_table_attach_defaults(GTK_TABLE(tbl), b3, 0, 2, row, row + 1);
    g_signal_connect(G_OBJECT(b3), "clicked", G_CALLBACK(set_lam), opts);
    g_object_set_data(G_OBJECT(b3), "lam_value", GINT_TO_POINTER(-1));
    row++;

    if (opts->lam == 0) {
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b1), TRUE);
    } else if (opts->lam == 1) {
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b2), TRUE);
    } else if (opts->lam == -1) {
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b3), TRUE); 
    }

    /* horizontal separator */
    tmp = gtk_hseparator_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;

    /* mean correction: radio group */
    b1 = gtk_radio_button_new_with_label(NULL, _("Mean correction"));
    mean_group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    gtk_table_attach_defaults(GTK_TABLE(tbl), b1, 0, 2, row, row + 1);
    g_signal_connect(G_OBJECT(b1), "clicked", G_CALLBACK(set_imean), 
		     opts);
    g_object_set_data(G_OBJECT(b1), "imean_value", GINT_TO_POINTER(1));
    row++;

    b2 = gtk_radio_button_new_with_label(mean_group, _("No mean correction"));
    mean_group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b2));
    gtk_table_attach_defaults(GTK_TABLE(tbl), b2, 0, 2, row, row + 1);
    g_signal_connect(G_OBJECT(b2), "clicked", G_CALLBACK(set_imean), 
		     opts);
    g_object_set_data(G_OBJECT(b2), "imean_value", GINT_TO_POINTER(0));

    if (opts->imean == 0) {
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b2), TRUE);
    } else if (opts->imean == 1) {
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b1), TRUE);
    }
}

static void get_va_value (GtkSpinButton *sb, tramo_options *opts)
{
    opts->va = (float) gtk_spin_button_get_value(sb);
}

static void tramo_tab_outliers (GtkWidget *notebook, tramo_options *opts)
{
    GtkWidget *tbl, *tmp;
    GtkAdjustment *adj;
    int tbl_len = 9, row = 0;

    tbl = make_notebook_page_table(notebook, _("Outliers"), tbl_len, 2, NULL);

    /* Overall choice: outlier correction or no? */
    tmp = gtk_check_button_new_with_label(_("Detect and correct for outliers"));
    opts->iatip_button = tmp;
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), (opts->iatip != 0));
    g_signal_connect(G_OBJECT(tmp), "clicked", G_CALLBACK(flip_iatip), 
		     opts);
    row++;

    /* horizontal separator */
    tmp = gtk_hseparator_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;
    
    /* label pertaining to transitory and level-shift buttons */
    tmp = gtk_label_new(_("Besides additive outliers, allow for:"));
    opts->aio_label = tmp;
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;

    /* "transitory" button */
    tmp = gtk_check_button_new_with_label(_("transitory changes"));
    opts->aio_transitory_button = tmp;
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), (opts->aio < 3));
    g_signal_connect(G_OBJECT(tmp), "clicked", G_CALLBACK(tramo_aio_callback), 
		     opts);
    row++;

    /* level-shift button */
    tmp = gtk_check_button_new_with_label(_("shifts of level"));
    opts->aio_shift_button = tmp;
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), (opts->aio > 1));
    g_signal_connect(G_OBJECT(tmp), "clicked", G_CALLBACK(tramo_aio_callback), 
		     opts);
    row++;

    /* innovationals button */
    tmp = gtk_check_button_new_with_label(_("innovational outliers"));
    opts->aio_innov_button = tmp;
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), (opts->aio == 0));
    gtk_widget_set_sensitive(tmp, opts->seats == 0);
    g_signal_connect(G_OBJECT(tmp), "clicked", G_CALLBACK(tramo_innov_callback), 
		     opts);
    row++;
    
    /* horizontal separator */
    tmp = gtk_hseparator_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;

    /* label pertaining to critical value */
    tmp = gtk_label_new(_("Critical value for outliers:"));
    opts->va_label = tmp;
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;

    /* check button for auto critical value */
    tmp = gtk_check_button_new_with_label(_("Automatic"));
    opts->va_button = tmp;
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), (opts->va == 0.0));
    g_signal_connect(G_OBJECT(tmp), "clicked", G_CALLBACK(flip_auto_va), 
		     opts);
    row++;

    /* spinner for manual critical value */
    adj = (GtkAdjustment *) gtk_adjustment_new((opts->va == 0.0)? 3.3 : opts->va, 
					       2.1, 6.0, 0.1, 1.0, 0);
    tmp = gtk_spin_button_new(adj, 0.1, 1);
    opts->va_spinner = tmp;
    gtk_table_attach(GTK_TABLE(tbl), tmp, 0, 1, row, row + 1,
		     0, 0, 0, 0);
    gtk_widget_set_sensitive(tmp, (opts->va != 0.0));
    g_signal_connect(G_OBJECT(tmp), "value-changed",
		     G_CALLBACK(get_va_value), opts);
}

static void arima_spin_callback (GtkSpinButton *sb, gint *var)
{
    *var = (int) gtk_spin_button_get_value(sb);
}

static GtkWidget *arima_spinner (const gchar *label, 
				 GtkWidget *tbl, gint row,
				 gint imax, gint *var)
{
    GtkWidget *w;
    GtkAdjustment *adj;

    w = gtk_label_new(label);
    gtk_label_set_justify(GTK_LABEL(w), GTK_JUSTIFY_RIGHT);
    gtk_table_attach(GTK_TABLE(tbl), w, 0, 1, row, row + 1,
		     0, 0, 0, 0);
    gtk_widget_show(w);

    adj = (GtkAdjustment *) gtk_adjustment_new(*var, 0, imax, 1, 1, 0);
    w = gtk_spin_button_new(adj, 1, 0);
    gtk_table_attach(GTK_TABLE(tbl), w, 1, 2, row, row + 1,
		     0, 0, 0, 0);
    g_signal_connect(G_OBJECT(w), "value-changed",
		     G_CALLBACK(arima_spin_callback), 
		     var);
    gtk_widget_show(w);

    return w;
}

static void tramo_tab_arima (GtkWidget *notebook, tramo_options *opts, int pd)
{
    GtkWidget *tbl, *tmp;
    int tbl_len = (pd == 1)? 7 : 10;
    int row = 0;

    tbl = make_notebook_page_table(notebook, _("ARIMA"), tbl_len, 2, NULL);
    gtk_table_set_homogeneous(GTK_TABLE(tbl), FALSE);

    /* auto versus manual button */
    tmp = gtk_check_button_new_with_label(_("Automatic"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;
    gtk_widget_show(tmp);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), opts->auto_arima);

    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(flip_auto_arima), 
		     opts);

    /* difference terms */
    opts->d_spin = arima_spinner(_("Non-seasonal differences:"), tbl, row++,
				 3, &opts->d);

    if (pd > 1) {
	opts->bd_spin = arima_spinner(_("Seasonal differences:"), tbl, row++,
				      2, &opts->bd);
    } else {
	opts->bd_spin = NULL;
    }
	    
    /* horizontal separator */
    tmp = gtk_hseparator_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;
    gtk_widget_show(tmp);

    /* AR terms */
    opts->p_spin = arima_spinner(_("Non-seasonal AR terms:"), tbl, row++,
				 10, &opts->p);

    if (pd > 1) {
	opts->bp_spin = arima_spinner(_("Seasonal AR terms:"), tbl, row++,
				      2, &opts->bp);
    } else {
	opts->bp_spin = NULL;
    }

    /* horizontal separator */
    tmp = gtk_hseparator_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;
    gtk_widget_show(tmp);

    /* MA terms */
    opts->q_spin = arima_spinner(_("Non-seasonal MA terms:"), tbl, row++,
				 10, &opts->q);
    
    if (pd > 1) {
	opts->bq_spin = arima_spinner(_("Seasonal MA terms:"), tbl, row,
				      2, &opts->bq);
    } else {
	opts->bq_spin = NULL;
    }

    arima_options_set_sensitive(opts, (opts->auto_arima == 0));
}

static void real_show_tramo_options (tx_request *request, GtkWidget *vbox) 
{
    tramo_options *opts = (tramo_options *) request->gui;
    GtkWidget *book;

    book = gtk_notebook_new();
    gtk_box_pack_start(GTK_BOX(vbox), book, TRUE, TRUE, 0);
    gtk_widget_show(book);

    tramo_tab_general(book, request);    
    tramo_tab_output(book, request);
    tramo_tab_outliers(book, opts);
    tramo_tab_transform(book, opts);
    tramo_tab_arima(book, opts, request->pd);

    if (opts->rsa == 3) {
	main_auto_callback(NULL, book);
    }
}

static tramo_options *tramo_options_new (int pd)
{
    tramo_options *opts;

    opts = malloc(sizeof *opts);
    if (opts == NULL) return NULL;

    if (pd == 4 || pd == 12) {
	tramo_options_set_defaults(opts, pd);
    } else {
	tramo_options_set_defaults(opts, 0);
    }

    opts->iatip_button = NULL;
    opts->aio_transitory_button = NULL;
    opts->aio_shift_button = NULL;
    opts->va_button = NULL;
    opts->va_spinner = NULL;
    opts->aio_label = NULL;
    opts->va_label = NULL;    

    return opts;
}

int add_tramo_options (tx_request *request, GtkWidget *vbox)
{
    tramo_options *opts;

    opts = tramo_options_new(request->pd);
    if (opts == NULL) {
	return 1;
    }

    /* mutual pointer hook-up */
    request->gui = opts;
    opts->request = request;

    real_show_tramo_options(request, vbox);

    return 0;
}

/* Below: print then free the tramo options structure.
   Return an indication of whether SEATS will be run (1) or not (0)
*/

int print_tramo_options (tx_request *request, FILE *fp)
{
    tramo_options *opts;

    if (request->gui == NULL) {
	return 0;
    }

    opts = request->gui;

    fputs("$INPUT ", fp);

    if (opts->rsa == 3) {
	fputs("rsa=3,", fp);
	goto set_out;
    }

    /* note: if values are at their TRAMO defaults, don't bother
       printing them */

    if (opts->lam != -1) {
	fprintf(fp, "lam=%d,", opts->lam);
    }

    if (opts->imean != 1) {
	fprintf(fp, "imean=%d,", opts->imean);
    }

    fprintf(fp, "iatip=%d,", opts->iatip);

    if (opts->iatip == 1) {
	if (opts->aio != 2) {
	    fprintf(fp, "aio=%d,", opts->aio);
	}
	if (opts->va != 0.0) {
	    fprintf(fp, "va=%.1f,", opts->va);
	}
    }

    if (!opts->auto_arima) {
	fprintf(fp, "D=%d,BD=%d,", opts->d, opts->bd);
	fprintf(fp, "P=%d,BP=%d,", opts->p, opts->bp);
	fprintf(fp, "Q=%d,BQ=%d,", opts->q, opts->bq);
    } else {
	fprintf(fp, "inic=%d,", opts->inic);
	fprintf(fp, "idif=%d,", opts->idif);
    }

    if (opts->mq > 0) {
	fprintf(fp, "mq=%d,", opts->mq);
    }

    if (opts->noadmiss != 1) {
	fprintf(fp, "noadmiss=%d,", opts->noadmiss);
    }

    fprintf(fp, "seats=%d,", opts->seats);

 set_out:
    if (opts->out != 0) {
	fprintf(fp, "out=%d,", opts->out);
    }

    fputs("$END\n", fp);

    free(opts);
    request->gui = NULL;

    return (opts->seats > 0);
}

