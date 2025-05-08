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

/* filters.c : for gretl */

#include "gretl.h"
#include "gui_utils.h"
#include "dlgutils.h"
#include "filters.h"
#include "libset.h"
#include "plotspec.h"
#include "uservar.h"
#include "usermat.h"
#include "cmdstack.h"
#include "ssheet.h"

enum {
    FILTER_SMA = 1,
    FILTER_EMA,
    FILTER_HP,
    FILTER_BK,
    FILTER_BW,
    FILTER_POLY,
    FILTER_FD
};

#define FILTER_SAVE_NONE 0
#define FILTER_SAVE_TREND OPT_A
#define FILTER_SAVE_CYCLE OPT_B

#define FILTER_GRAPH_NONE 0
#define FILTER_GRAPH_TREND OPT_A
#define FILTER_GRAPH_CYCLE OPT_B
#define FILTER_GRAPH_RESP OPT_C

typedef struct filter_info_ filter_info;

struct filter_info_ {
    int ftype;            /* type of filter */
    int vnum;             /* ID number of series to filter */
    const char *vname;    /* name of series */
    int t1;               /* starting observation */
    int t2;               /* ending observation */
    int k;                /* no. of terms in MA, or no. of obs for mean */
    int center;           /* for simple MA: center or not */
    int oneside;          /* for HP filter: one-sided variant? */
    double lambda;        /* general purpose parameter */
    double fx0;           /* pre-sample EMA value */
    int bkk;              /* Baxter-King parameter */
    int bkl;              /* Baxter-King lower value */
    int bku;              /* Baxter-King upper value */
    int order;            /* Butterworth or polynomial order */
    int cutoff;           /* Butterworth cut-off */
    double midfrac;       /* for weighted polynomial fit */
    gretlopt wopt;        /* weighting-type option */
    gretlopt graph_opt;
    gretlopt save_opt;
    char save_t[VNAMELEN];
    char save_c[VNAMELEN];
    char label_t[MAXLABEL];
    char label_c[MAXLABEL];
    GtkWidget *dlg;
    GtkWidget *entry1;
    GtkWidget *entry2;
    GtkWidget *kspin;
    GtkWidget *fx0spin;
    GtkWidget *spin1;
    GtkWidget *spin2;
    GtkWidget *button;
    GtkWidget *label;
};

static int calculate_filter (filter_info *finfo);
static void weights_shape_graph (GtkWidget *button, filter_info *finfo);
static void butterworth_poles_graph (GtkWidget *button, filter_info *finfo);

static const char *filter_get_title (int ftype)
{
    if (ftype == FILTER_SMA) {
	return N_("Simple moving average");
    } else if (ftype == FILTER_EMA) {
	return N_("Exponential moving average");
    } else if (ftype == FILTER_HP) {
	return N_("Hodrick-Prescott filter");
    } else if (ftype == FILTER_BK) {
	return N_("Baxter-King Band-pass filter");
    } else if (ftype == FILTER_BW) {
	return N_("Butterworth filter");
    } else if (ftype == FILTER_POLY) {
	return N_("Polynomial trend");
    } else if (ftype == FILTER_FD) {
	return N_("Fractional difference");
    } else {
	return "";
    }
}

static void filter_info_init (filter_info *finfo, int ftype, int v,
			      int t1, int t2)
{
    int pd = dataset->pd;

    finfo->ftype = ftype;
    finfo->vnum = v;
    finfo->vname = dataset->varname[v];
    finfo->t1 = t1;
    finfo->t2 = t2;
    finfo->k = 0;
    finfo->center = 0;
    finfo->oneside = 0;
    finfo->lambda = 0.0;
    finfo->fx0 = NADBL;
    finfo->bkk = 0;
    finfo->bkl = 0;
    finfo->bku = 0;

    if (ftype == FILTER_BK) {
	finfo->graph_opt = FILTER_GRAPH_CYCLE;
    } else {
	finfo->graph_opt = FILTER_GRAPH_TREND;
    }

    finfo->save_opt = FILTER_SAVE_NONE;

    *finfo->save_t = 0;
    *finfo->save_c = 0;
    *finfo->label_t = 0;
    *finfo->label_c = 0;

    finfo->dlg = NULL;
    finfo->entry1 = NULL;
    finfo->entry2 = NULL;
    finfo->kspin = NULL;
    finfo->fx0spin = NULL;
    finfo->spin1 = NULL;
    finfo->spin2 = NULL;
    finfo->button = NULL;
    finfo->label = NULL;

    if (ftype == FILTER_SMA) {
	finfo->center = 0;
	finfo->k = (pd == 1 || pd == 10)? 3 : pd;
    } else if (ftype == FILTER_EMA) {
	finfo->lambda = 0.2;
    } else if (ftype == FILTER_HP) {
	finfo->lambda = 100 * pd * pd;
    } else if (ftype == FILTER_BK) {
	finfo->bkk = get_bkbp_k(dataset);
	get_bkbp_periods(dataset, &finfo->bkl, &finfo->bku);
    } else if (ftype == FILTER_BW) {
	finfo->order = 8;
	finfo->cutoff = 67;
    } else if (ftype == FILTER_POLY) {
	finfo->order = 3;
	finfo->wopt = OPT_NONE;
	finfo->lambda = 1.0;
	finfo->midfrac = 0.0;
    } else if (ftype == FILTER_FD) {
	finfo->lambda = 0.5;
    }
}

static void filter_make_savename (filter_info *finfo, int i)
{
    char *targ = (i == 0)? finfo->save_t : finfo->save_c;
    int len;

    if (finfo->ftype == FILTER_SMA) {
	strcpy(targ, (i == 0)? "ma_" : "mc_");
    } else if (finfo->ftype == FILTER_EMA) {
	strcpy(targ, (i == 0)? "ema_" : "emc_");
    } else if (finfo->ftype == FILTER_HP) {
	strcpy(targ, (i == 0)? "hpt_" : "hp_");
    } else if (finfo->ftype == FILTER_BK) {
	strcpy(targ, "bk_");
    } else if (finfo->ftype == FILTER_BW) {
	strcpy(targ, (i == 0)? "bt_" : "bc_");
    } else if (finfo->ftype == FILTER_POLY) {
	strcpy(targ, (i == 0)? "pt_" : "pc_");
    } else if (finfo->ftype == FILTER_FD) {
	strcpy(targ, "fd_");
    }

    len = strlen(targ);
    strncat(targ, finfo->vname, VNAMELEN - len - 1);
}

static void filter_make_varlabel (filter_info *finfo, int v, int i)
{
    gchar *s = NULL;

    if (finfo->ftype == FILTER_SMA) {
	if (i == FILTER_SAVE_TREND) {
	    if (finfo->center) {
		s = g_strdup_printf(_("Centered %d-period moving average of %s"),
				    finfo->k, finfo->vname);
	    } else {
		s = g_strdup_printf(_("%d-period moving average of %s"),
				    finfo->k, finfo->vname);
	    }
	} else {
	    s = g_strdup_printf(_("Residual from %d-period MA of %s"),
				finfo->k, finfo->vname);
	}
    } else if (finfo->ftype == FILTER_EMA) {
	if (i == FILTER_SAVE_TREND) {
	    s = g_strdup_printf(_("Exponential moving average of %s (current weight %g)"),
				finfo->vname, finfo->lambda);
	} else {
	    s = g_strdup_printf(_("Residual from EMA of %s (current weight %g)"),
				finfo->vname, finfo->lambda);
	}
    } else if (finfo->ftype == FILTER_HP) {
	if (finfo->oneside) {
	    if (i == FILTER_SAVE_TREND) {
		s = g_strdup_printf(_("Filtered %s: one-sided Hodrick-Prescott trend (lambda = %g)"),
				    finfo->vname, finfo->lambda);
	    } else {
		s = g_strdup_printf(_("Filtered %s: one-sided Hodrick-Prescott cycle (lambda = %g)"),
				    finfo->vname, finfo->lambda);
	    }
	} else {
	    if (i == FILTER_SAVE_TREND) {
		s = g_strdup_printf(_("Filtered %s: Hodrick-Prescott trend (lambda = %g)"),
				    finfo->vname, finfo->lambda);
	    } else {
		s = g_strdup_printf(_("Filtered %s: Hodrick-Prescott cycle (lambda = %g)"),
				    finfo->vname, finfo->lambda);
	    }
	}
    } else if (finfo->ftype == FILTER_BK) {
	s = g_strdup_printf(_("Filtered %s: Baxter-King, frequency %d to %d"),
			    finfo->vname, finfo->bkl, finfo->bku);
    } else if (finfo->ftype == FILTER_BW) {
	if (i == FILTER_SAVE_TREND) {
	    s = g_strdup_printf(_("Filtered %s: Butterworth low-pass (n=%d, cutoff=%d)"),
				finfo->vname, finfo->order, finfo->cutoff);
	} else {
	    s = g_strdup_printf(_("Filtered %s: Butterworth high-pass (n=%d, cutoff=%d)"),
				finfo->vname, finfo->order, finfo->cutoff);
	}
    } else if (finfo->ftype == FILTER_POLY) {
	if (i == FILTER_SAVE_TREND) {
	    s = g_strdup_printf(_("Filtered %s: polynomial of order %d"),
				finfo->vname, finfo->order);
	} else {
	    s = g_strdup_printf(_("Filtered %s: residual from polynomial of order %d"),
				finfo->vname, finfo->order);
	}
    } else if (finfo->ftype == FILTER_FD) {
	s = g_strdup_printf("fracdiff(%s, %g)", finfo->vname, finfo->lambda);
    }

    if (s != NULL) {
	gretl_utf8_truncate(s, MAXLABEL-1);
    } else {
	s = g_strdup("");
    }

    if (i == FILTER_SAVE_TREND) {
	strcpy(finfo->label_t, s);
    } else {
	strcpy(finfo->label_c, s);
    }

    series_set_label(dataset, v, s);
    g_free(s);
}

static void check_bk_limits1 (GtkWidget *s1, GtkWidget *s2)
{
    int n1, n2;

    n1 = (int) gtk_spin_button_get_value(GTK_SPIN_BUTTON(s1));
    n2 = (int) gtk_spin_button_get_value(GTK_SPIN_BUTTON(s2));

    if (n1 >= n2) {
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(s2), n1 + 1);
    }
}

static void check_bk_limits2 (GtkWidget *s2, GtkWidget *s1)
{
    int n1, n2;

    n1 = (int) gtk_spin_button_get_value(GTK_SPIN_BUTTON(s1));
    n2 = (int) gtk_spin_button_get_value(GTK_SPIN_BUTTON(s2));

    if (n1 >= n2) {
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(s1), n2 - 1);
    }
}

static void binary_option_callback (GtkWidget *w, int *o)
{
    *o = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));
}

static int varname_error (filter_info *finfo, int i)
{
    const char *s = (i == 1)? finfo->save_t : finfo->save_c;

    if (*s == 0) {
	errbox(_("Variable name is missing"));
	return 1;
    } else if (gui_validate_varname(s, GRETL_TYPE_SERIES, finfo->dlg)) {
	return 1;
    } else if (i == 2 && (finfo->save_opt & FILTER_SAVE_TREND)) {
	if (!strcmp(s, finfo->save_t)) {
	    errbox(_("Conflicting variable names"));
	    return 1;
	}
    }

    return 0;
}

static void filter_dialog_hsep (GtkWidget *dlg)
{
    GtkWidget *vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));
    GtkWidget *hs = gtk_hseparator_new();

    gtk_box_pack_start(GTK_BOX(vbox), hs, TRUE, TRUE, 0);
}

static void filter_graph_check_button (GtkWidget *dlg, filter_info *finfo,
				       const char *txt, gretlopt g)
{
    GtkWidget *vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));
    GtkWidget *hbox, *w;

    hbox = gtk_hbox_new(FALSE, 5);
    w = gretl_option_check_button(txt, &finfo->graph_opt, g);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), finfo->graph_opt & g);
    gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
}

static void filter_save_check_buttons (GtkWidget *dlg, filter_info *finfo)
{
    GtkWidget *vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));
    GtkWidget *hbox, *tab, *w;
    int i, imax = 2;

    if (finfo->ftype == FILTER_FD) {
	imax = 1;
    }

    hbox = gtk_hbox_new(FALSE, 5);
    tab = gtk_table_new(2, 2, FALSE);
    gtk_table_set_col_spacing(GTK_TABLE(tab), 0, 5);

    for (i=0; i<imax; i++) {
	if (finfo->ftype == FILTER_FD) {
	    w = gretl_option_check_button(_("Save output as"),
					  &finfo->save_opt,
					  FILTER_SAVE_TREND);
	} else if (i == 0) {
	    w = gretl_option_check_button(_("Save smoothed series as"),
					  &finfo->save_opt,
					  FILTER_SAVE_TREND);
	} else {
	    w = gretl_option_check_button(_("Save cyclical component as"),
					  &finfo->save_opt,
					  FILTER_SAVE_CYCLE);
	}
	gtk_table_attach_defaults(GTK_TABLE(tab), w, 0, 1, i, i+1);

	w = gtk_entry_new();
	gtk_entry_set_max_length(GTK_ENTRY(w), VNAMELEN - 1);
	gtk_entry_set_width_chars(GTK_ENTRY(w), VNAMELEN + 3);
	filter_make_savename(finfo, i);
	gtk_entry_set_text(GTK_ENTRY(w), (i == 0)? finfo->save_t : finfo->save_c);
	gtk_entry_set_activates_default(GTK_ENTRY(w), TRUE);
	gtk_table_attach_defaults(GTK_TABLE(tab), w, 1, 2, i, i+1);
	if (i == 0) {
	    finfo->entry1 = w;
	} else {
	    finfo->entry2 = w;
	}
    }

    gtk_box_pack_start(GTK_BOX(hbox), tab, TRUE, TRUE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
}

static void bkbp_frequencies_table (GtkWidget *dlg, filter_info *finfo)
{
    GtkWidget *vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));
    GtkWidget *tab, *hbox, *w;
    GtkWidget *s1, *s2;

    hbox = gtk_hbox_new(FALSE, 5);
    tab = gtk_table_new(2, 2, FALSE);
    gtk_table_set_col_spacing(GTK_TABLE(tab), 0, 5);

    int T = finfo->t2 - finfo->t1 + 1;

    /* lower limit */
    w = gtk_label_new(_("Lower frequency bound:"));
    gtk_table_attach_defaults(GTK_TABLE(tab), w, 0, 1, 0, 1);
    s1 = gtk_spin_button_new_with_range(1, 64, 1);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(s1), finfo->bkl);
    g_signal_connect(G_OBJECT(s1), "value-changed",
		     G_CALLBACK(set_int_from_spin), &finfo->bkl);
    gtk_table_attach_defaults(GTK_TABLE(tab), s1, 1, 2, 0, 1);

    /* upper limit */
    w = gtk_label_new(_("Upper frequency bound:"));
    gtk_table_attach_defaults(GTK_TABLE(tab), w, 0, 1, 1, 2);
    s2 = gtk_spin_button_new_with_range(2, (double) T, 1);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(s2), finfo->bku);
    g_signal_connect(G_OBJECT(s2), "value-changed",
		     G_CALLBACK(set_int_from_spin), &finfo->bku);
    gtk_table_attach_defaults(GTK_TABLE(tab), s2, 1, 2, 1, 2);

    /* inter-connect the lower and upper spin buttons */
    g_signal_connect(G_OBJECT(s1), "value-changed",
		     G_CALLBACK(check_bk_limits1), s2);
    g_signal_connect(G_OBJECT(s2), "value-changed",
		     G_CALLBACK(check_bk_limits2), s1);

    /* pack everything */
    gtk_box_pack_start(GTK_BOX(hbox), tab, TRUE, TRUE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
}

static void set_ema_init (GtkWidget *b, filter_info *finfo)
{
    gboolean t = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(b));
    int i = widget_get_int(b, "idx");

    if (i == 1) {
	gtk_widget_set_sensitive(finfo->kspin, t);
	gtk_widget_set_sensitive(finfo->fx0spin, !t);
    } else if (i == 2) {
	gtk_widget_set_sensitive(finfo->kspin, !t);
	gtk_widget_set_sensitive(finfo->fx0spin, !t);
    } else if (i == 3) {
	gtk_widget_set_sensitive(finfo->kspin, !t);
	gtk_widget_set_sensitive(finfo->fx0spin, t);
    }
}

static void ema_obs_radios (GtkWidget *dlg, filter_info *finfo)
{
    GtkWidget *vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));
    GtkWidget *tab, *hbox, *w;
    GtkWidget *r1, *r2, *r3;
    GSList *group;

    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_label_new(_("The initial EMA value is"));
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    hbox = gtk_hbox_new(FALSE, 5);
    tab = gtk_table_new(3, 2, FALSE);
    gtk_table_set_col_spacing(GTK_TABLE(tab), 0, 5);

    /* default */
    r1 = gtk_radio_button_new_with_label(NULL, _("the mean of the first n observations"));
    gtk_table_attach_defaults(GTK_TABLE(tab), r1, 0, 1, 0, 1);
    w = gtk_spin_button_new_with_range(1, (finfo->t2 - finfo->t1 + 1) / 2, 1);
    gtk_table_attach_defaults(GTK_TABLE(tab), w, 1, 2, 0, 1);
    finfo->kspin = w;
    widget_set_int(r1, "idx", 1);

    /* alternate */
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(r1));
    r2 = gtk_radio_button_new_with_label(group, _("the mean of the whole series"));
    gtk_table_attach_defaults(GTK_TABLE(tab), r2, 0, 1, 1, 2);
    widget_set_int(r2, "idx", 2);

    /* further alternate */
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(r2));
    r3 = gtk_radio_button_new_with_label(group, _("a user-specified value"));
    gtk_table_attach_defaults(GTK_TABLE(tab), r3, 0, 1, 2, 3);
    w = gtk_spin_button_new_with_range(-1.0e6, 1.0e6, 0.001);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(w), 0.0);
    gtk_widget_set_sensitive(w, FALSE);
    gtk_table_attach_defaults(GTK_TABLE(tab), w, 1, 2, 2, 3);
    finfo->fx0spin = w;
    widget_set_int(r3, "idx", 3);

    /* connect switching signals */
    g_signal_connect(G_OBJECT(r1), "toggled", G_CALLBACK(set_ema_init), finfo);
    g_signal_connect(G_OBJECT(r2), "toggled", G_CALLBACK(set_ema_init), finfo);
    g_signal_connect(G_OBJECT(r3), "toggled", G_CALLBACK(set_ema_init), finfo);

    /* pack everything */
    gtk_box_pack_start(GTK_BOX(hbox), tab, TRUE, TRUE, 10);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
}

static void weight_scheme_changed (GtkComboBox *combo, filter_info *finfo)
{
    int i = gtk_combo_box_get_active(combo);

    gtk_widget_set_sensitive(finfo->spin1, (i > 0));
    gtk_widget_set_sensitive(finfo->spin2, (i > 0));
    gtk_widget_set_sensitive(finfo->button, (i > 0));

    if (i == 0) {
	finfo->wopt = OPT_NONE;
    } else if (i == 1) {
	finfo->wopt = OPT_Q;
    } else if (i == 2) {
	finfo->wopt = OPT_B;
    } else {
	finfo->wopt = OPT_C;
    }
}

static void add_poly_weights_options (GtkWidget *vbox, filter_info *finfo)
{
    GtkWidget *hbox, *label, *combo;
    const char *wstrs[] = {
	N_("uniform"),
	N_("quadratic"),
	N_("cosine-bell"),
	N_("steps")
    };
    int i;

    /* first row: weighting scheme */
    hbox = gtk_hbox_new(FALSE, 5);
    label = gtk_label_new(_("Weights:"));
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    combo = gtk_combo_box_text_new();
    for (i=0; i<4; i++) {
	combo_box_append_text(combo, _(wstrs[i]));
    }
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo), 0);
    g_signal_connect(GTK_COMBO_BOX(combo), "changed",
		     G_CALLBACK(weight_scheme_changed), finfo);
    gtk_box_pack_start(GTK_BOX(hbox), combo, FALSE, FALSE, 0);
    label = gtk_label_new(_("with maximum:"));
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    finfo->spin1 = gtk_spin_button_new_with_range(1.0, 10.0, 0.1);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(finfo->spin1), 2.0);
    gtk_widget_set_sensitive(finfo->spin1, FALSE);
    gtk_box_pack_start(GTK_BOX(hbox), finfo->spin1, FALSE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

    /* second row: size of central section plus show button */
    hbox = gtk_hbox_new(FALSE, 5);
    label = gtk_label_new(_("Central fraction:"));
    gtk_box_pack_start(GTK_BOX(hbox), label, FALSE, FALSE, 5);
    finfo->spin2 = gtk_spin_button_new_with_range(0.0, 1.0, 0.1);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(finfo->spin2), 0.4);
    gtk_widget_set_sensitive(finfo->spin2, FALSE);
    gtk_box_pack_start(GTK_BOX(hbox), finfo->spin2, FALSE, FALSE, 0);
    finfo->button = gtk_button_new_with_label(_("Show weights"));
    gtk_widget_set_sensitive(finfo->button, FALSE);
    g_signal_connect(GTK_BUTTON(finfo->button), "clicked",
		     G_CALLBACK(weights_shape_graph), finfo);
    gtk_box_pack_end(GTK_BOX(hbox), finfo->button, FALSE, FALSE, 15);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
}

static void filter_dialog_ok (GtkWidget *w, filter_info *finfo)
{
    if (finfo->save_opt & FILTER_SAVE_TREND) {
	strcpy(finfo->save_t, gtk_entry_get_text(GTK_ENTRY(finfo->entry1)));
	if (varname_error(finfo, 1)) {
	    return;
	}
    }

    if (finfo->save_opt & FILTER_SAVE_CYCLE) {
	strcpy(finfo->save_c, gtk_entry_get_text(GTK_ENTRY(finfo->entry2)));
	if (varname_error(finfo, 2)) {
	    return;
	}
    }

    if (finfo->kspin != NULL) {
	if (gtk_widget_is_sensitive(finfo->kspin)) {
	    finfo->k = (int)
		gtk_spin_button_get_value(GTK_SPIN_BUTTON(finfo->kspin));
	} else {
	    finfo->k = 0;
	}
    }

    if (finfo->fx0spin != NULL && gtk_widget_is_sensitive(finfo->fx0spin)) {
	finfo->fx0 = gtk_spin_button_get_value(GTK_SPIN_BUTTON(finfo->fx0spin));
    }

    if (finfo->spin1 != NULL && finfo->spin2 != NULL) {
	if (gtk_widget_is_sensitive(finfo->spin1)) {
	    finfo->lambda = gtk_spin_button_get_value(GTK_SPIN_BUTTON(finfo->spin1));
	    finfo->midfrac = gtk_spin_button_get_value(GTK_SPIN_BUTTON(finfo->spin2));
	}
    }

    if (finfo->graph_opt != FILTER_GRAPH_NONE ||
	finfo->save_opt != FILTER_SAVE_NONE) {
	int err = calculate_filter(finfo);

	if (err) {
	    gui_errmsg(err);
	}
    }
}

#if 0
static gboolean check_bw_feasibility (GtkSpinButton *b,
				      filter_info *finfo)
{
    if (finfo->spin1 != NULL && finfo->spin2 != NULL) {
	const double b0 = 0.131981;
	const double b1 = 0.0126241;
	const double b2 = 0.0658865;
	const double b3 = 0.883492;
	double s, lrc;
	int c, n;

	n = gtk_spin_button_get_value(GTK_SPIN_BUTTON(finfo->spin1));
	c = gtk_spin_button_get_value(GTK_SPIN_BUTTON(finfo->spin2));

	s = log(sin(c*M_PI/360));
	lrc = b0 + b1*n + b2*s + b3*n*s;

	fprintf(stderr, "lrc = %g\n", lrc);

	gtk_widget_set_sensitive(finfo->label, lrc < -12);
    }

    return FALSE;
}
#else
static gboolean check_bw_feasibility (GtkSpinButton *b,
				      filter_info *finfo)
{
    if (finfo->spin1 != NULL && finfo->spin2 != NULL) {
	double b0 = -8.56644923648263;
	double b1 = -0.66992007228155;
	double b2 =  2.82851981714539;
	double x;
	int c, n;

	n = gtk_spin_button_get_value(GTK_SPIN_BUTTON(finfo->spin1));
	c = gtk_spin_button_get_value(GTK_SPIN_BUTTON(finfo->spin2));

	/* let's be conservative here */
	c--;
	n++;

	x = b0 + b1*c + b2*n;
	x = 1/(1+exp(-x));

	gtk_widget_set_sensitive(finfo->label, x > 0.47971759156793);
    }

    return FALSE;
}
#endif

static void filter_dialog_quit (GtkWidget *w, gpointer data)
{
    set_dataset_locked(FALSE);
}

static void filter_dialog (filter_info *finfo)
{
    GtkWidget *dlg;
    GtkWidget *vbox;
    GtkWidget *hbox;
    GtkWidget *w;

    dlg = gretl_dialog_new(_("gretl: time-series filter"), mdata->main,
			   GRETL_DLG_BLOCK);
    finfo->dlg = dlg;

    g_signal_connect(G_OBJECT(dlg), "destroy",
		     G_CALLBACK(filter_dialog_quit), NULL);

    /* box title */
    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_label_new(_(filter_get_title(finfo->ftype)));
    gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 5);
    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dlg));
    gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 5);

    if (finfo->ftype == FILTER_SMA) {
	GtkWidget *kspin;
	int T = finfo->t2 - finfo->t1 + 1;
	int nmax = T / 2;

	/* select number of observations */
	hbox = gtk_hbox_new(FALSE, 5);
	w = gtk_label_new(_("Number of observations in average:"));
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	kspin = gtk_spin_button_new_with_range(2, (gdouble) nmax, 1);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(kspin), (gdouble) finfo->k);
	g_signal_connect(G_OBJECT(kspin), "value-changed",
			 G_CALLBACK(set_int_from_spin), &finfo->k);
	gtk_box_pack_start(GTK_BOX(hbox), kspin, TRUE, TRUE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

	/* select centered or not */
	hbox = gtk_hbox_new(FALSE, 5);
	w = gtk_check_button_new_with_label(_("Centered"));
	if (finfo->center) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), TRUE);
	}
	g_signal_connect(G_OBJECT(w), "toggled",
			 G_CALLBACK(binary_option_callback), &finfo->center);
	gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    } else if (finfo->ftype == FILTER_EMA) {
	/* set weight on current observation */
	hbox = gtk_hbox_new(FALSE, 5);
	w = gtk_label_new(_("Weight on current observation:"));
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	w = gtk_spin_button_new_with_range(0.001, 0.999, 0.001);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(w), finfo->lambda);
	g_signal_connect(G_OBJECT(w), "value-changed",
			 G_CALLBACK(set_double_from_spin), &finfo->lambda);
	gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
	/* set calculation of initial EMA value */
	ema_obs_radios(dlg, finfo);
    } else if (finfo->ftype == FILTER_HP) {
	/* set H-P "lambda" parameter */
	hbox = gtk_hbox_new(FALSE, 5);
	w = gtk_label_new(_("lambda (higher values -> smoother trend):"));
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	w = gtk_spin_button_new_with_range(1.0, 999999, 1);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(w), finfo->lambda);
	g_signal_connect(G_OBJECT(w), "value-changed",
			 G_CALLBACK(set_double_from_spin), &finfo->lambda);
	gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);

	/* select one-sided or not */
	hbox = gtk_hbox_new(FALSE, 5);
	w = gtk_check_button_new_with_label(_("One-sided"));
	if (finfo->oneside) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), TRUE);
	}
	g_signal_connect(G_OBJECT(w), "toggled",
			 G_CALLBACK(binary_option_callback), &finfo->oneside);
	gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    } else if (finfo->ftype == FILTER_BK) {
	/* set Baxter-King "k" */
	hbox = gtk_hbox_new(FALSE, 5);
	w = gtk_label_new(_("k (higher values -> better approximation):"));
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	w = gtk_spin_button_new_with_range(1.0, 4 * finfo->bkk, 1);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(w), finfo->bkk);
	g_signal_connect(G_OBJECT(w), "value-changed",
			 G_CALLBACK(set_int_from_spin), &finfo->bkk);
	gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
	/* set periods */
	bkbp_frequencies_table(dlg, finfo);
    } else if (finfo->ftype == FILTER_BW) {
	/* set Butterworth order */
	hbox = gtk_hbox_new(FALSE, 5);
	w = gtk_label_new(_("n (higher values -> better approximation):"));
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	finfo->spin1 = w = gtk_spin_button_new_with_range(1.0, 16, 1);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(w), finfo->order);
	g_signal_connect(G_OBJECT(w), "value-changed",
			 G_CALLBACK(set_int_from_spin), &finfo->order);
	g_signal_connect(G_OBJECT(w), "value-changed",
			 G_CALLBACK(check_bw_feasibility), finfo);
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
	/* set cutoff */
	hbox = gtk_hbox_new(FALSE, 5);
	w = gtk_label_new(_("frequency cutoff (degrees):"));
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
        finfo->spin2 = w = gtk_spin_button_new_with_range(1, 179, 1);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(w), finfo->cutoff);
	g_signal_connect(G_OBJECT(w), "value-changed",
			 G_CALLBACK(set_int_from_spin), &finfo->cutoff);
	g_signal_connect(G_OBJECT(w), "value-changed",
			 G_CALLBACK(check_bw_feasibility), finfo);
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
	hbox = gtk_hbox_new(FALSE, 5);
	w = gtk_button_new_with_label(_("show poles"));
	g_signal_connect(G_OBJECT(w), "clicked",
			 G_CALLBACK(butterworth_poles_graph), finfo);
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	finfo->label = w = gtk_label_new(_("filter may be numerically unstable"));
	gtk_widget_set_sensitive(w, FALSE);
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    } else if (finfo->ftype == FILTER_POLY) {
	/* set polynomial order */
	hbox = gtk_hbox_new(FALSE, 5);
	w = gtk_label_new(_("Order:"));
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	w = gtk_spin_button_new_with_range(1.0, 15, 1);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(w), finfo->order);
	g_signal_connect(G_OBJECT(w), "value-changed",
			 G_CALLBACK(set_int_from_spin), &finfo->order);
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
	/* weights settings (optional) */
	filter_dialog_hsep(dlg);
	add_poly_weights_options(vbox, finfo);
    } else if (finfo->ftype == FILTER_FD) {
	/* set "d" */
	hbox = gtk_hbox_new(FALSE, 5);
	w = gtk_label_new(_("differencing parameter:"));
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	w = gtk_spin_button_new_with_range(-10, 10, 0.001);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(w), finfo->lambda);
	g_signal_connect(G_OBJECT(w), "value-changed",
			 G_CALLBACK(set_double_from_spin), &finfo->lambda);
	gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 5);
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    }

    filter_dialog_hsep(dlg);

    if (finfo->ftype == FILTER_BK) {
	/* graphical output? */
	filter_graph_check_button(dlg, finfo, _("Graph filtered series"),
				  FILTER_GRAPH_CYCLE);
	/* save to dataset? */
	hbox = gtk_hbox_new(FALSE, 5);
	w = gretl_option_check_button(_("Save filtered series as"),
				      &finfo->save_opt, FILTER_SAVE_CYCLE);
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	w = gtk_entry_new();
	gtk_entry_set_max_length(GTK_ENTRY(w), VNAMELEN - 1);
	gtk_entry_set_width_chars(GTK_ENTRY(w), VNAMELEN + 3);
	filter_make_savename(finfo, 1);
	gtk_entry_set_text(GTK_ENTRY(w), finfo->save_c);
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	finfo->entry2 = w;
	gtk_box_pack_start(GTK_BOX(vbox), hbox, FALSE, FALSE, 0);
    } else if (finfo->ftype == FILTER_FD) {
	/* graphical output? */
	filter_graph_check_button(dlg, finfo, _("Graph differenced series"),
				  FILTER_GRAPH_TREND);
	/* save to dataset? */
	filter_save_check_buttons(dlg, finfo);
    } else {
	/* graphical output? */
	filter_graph_check_button(dlg, finfo, _("Graph original and smoothed series"),
				  FILTER_GRAPH_TREND);
	filter_graph_check_button(dlg, finfo, _("Graph residual or cycle series"),
				  FILTER_GRAPH_CYCLE);
	if (finfo->ftype == FILTER_HP || finfo->ftype == FILTER_BW) {
	    filter_graph_check_button(dlg, finfo, _("Graph frequency response"),
				      FILTER_GRAPH_RESP);
	}
	/* save to dataset? */
	filter_dialog_hsep(dlg);
	filter_save_check_buttons(dlg, finfo);
    }

    /* convenience pointer for button box */
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dlg));

    /* Close button */
    w = close_button(hbox);
    g_signal_connect_swapped(G_OBJECT(w), "clicked",
			     G_CALLBACK(gtk_widget_destroy), dlg);

    /* "OK" button */
    w = ok_button(hbox);
    g_signal_connect(G_OBJECT(w), "clicked",
		     G_CALLBACK(filter_dialog_ok), finfo);
    gtk_widget_grab_default(w);

    if (finfo->ftype == FILTER_BW) {
	context_help_button(hbox, BWFILTER);
    } else if (finfo->ftype == FILTER_POLY) {
	context_help_button(hbox, POLYWEIGHTS);
    } else if (finfo->ftype == FILTER_EMA) {
	context_help_button(hbox, EMAFILTER);
    }

    set_dataset_locked(TRUE);

    gtk_widget_show_all(dlg);
}

static void print_gp_filter_data (filter_info *finfo,
				  const double *obs,
				  const double *x,
				  FILE *fp)
{
    int t;

    for (t=finfo->t1; t<=finfo->t2; t++) {
	/* same print format as in function printvars() in
	   lib/src/graphing.c
	*/
	if (na(x[t])) {
	    fprintf(fp, "%.10g ?\n", obs[t]);
	} else {
	    fprintf(fp, "%.10g %.10g\n", obs[t], x[t]);
	}
    }
}

static int
do_filter_graph (filter_info *finfo, const double *fx, const double *u)
{
    gchar *xtitle = NULL;
    gchar *ztitle = NULL;
    gchar *title = NULL;
    int twoplot = 0;
    int zkeypos = 'R';
    FILE *fp = NULL;
    const double *obs;
    int v, err = 0;

    obs = gretl_plotx(dataset, OPT_NONE);
    if (obs == NULL) {
	return E_ALLOC;
    }

    if ((finfo->graph_opt & FILTER_GRAPH_TREND) &&
	(finfo->graph_opt & FILTER_GRAPH_CYCLE)) {
	twoplot = 1;
    }

    if (twoplot) {
	fp = open_plot_input_file(PLOT_TRI_GRAPH, 0, &err);
    } else {
	fp = open_plot_input_file(PLOT_REGULAR, GPT_LETTERBOX, &err);
    }
    if (err) {
	return err;
    }

    if (!twoplot) {
	fprintf(fp, "# timeseries %d (letterbox)\n", dataset->pd);
    }

    if (dataset->pd == 4) {
	if ((finfo->t2 - finfo->t1) / 4 < 8) {
	    fputs("set xtics nomirror 0,1\n", fp);
	    fputs("set mxtics 4\n", fp);
	}
    } else if (dataset->pd == 12) {
	if ((finfo->t2 - finfo->t1) / 12 < 8) {
	    fputs("set xtics nomirror 0,1\n", fp);
	    fputs("set mxtics 12\n", fp);
	}
    }

    v = finfo->vnum;

    if (finfo->graph_opt & FILTER_GRAPH_TREND) {
	if (dataset->Z[v][finfo->t2] > dataset->Z[v][finfo->t1]) {
	    zkeypos = 'L';
	}
    }

    gretl_push_c_numeric_locale();

    fprintf(fp, "set xrange [%g:%g]\n", floor(obs[dataset->t1]),
	    ceil(obs[dataset->t2]));

    if (twoplot) {
	fputs("set size 1.0,1.0\nset multiplot\nset size 1.0,0.60\n", fp);

	fputs("set origin 0.0,0.4\n", fp);
	if (zkeypos == 'L') {
	    print_keypos_string(GP_KEY_LEFT | GP_KEY_TOP, fp);
	}
	xtitle = g_strdup_printf(_("%s (original data)"), finfo->vname);
	ztitle = g_strdup_printf(_("%s (smoothed)"), finfo->vname);
	fprintf(fp, "plot '-' using 1:2 title '%s' w lines, \\\n"
		" '-' using 1:2 title '%s' w lines\n", xtitle, ztitle);
	print_gp_filter_data(finfo, obs, dataset->Z[v], fp);
	fputs("e , \\\n", fp);
	print_gp_filter_data(finfo, obs, fx, fp);
	fputs("e\n", fp);

	fputs("set size 1.0,0.38\n", fp);
	fputs("set origin 0.0,0.0\n", fp);
	fputs("set xzeroaxis\n", fp);
	title = g_strdup_printf(_("Cyclical component of %s"), finfo->vname);
	fprintf(fp, "plot '-' using 1:2 title '%s' w lines\n", title);
	print_gp_filter_data(finfo, obs, u, fp);
	fputs("e\n", fp);
	fputs("unset multiplot\n", fp);
    } else if (finfo->ftype == FILTER_FD) {
	ztitle = g_strdup_printf("fracdiff(%s, %g)", finfo->vname, finfo->lambda);
	fprintf(fp, "plot '-' using 1:2 title '%s' w lines\n", ztitle);
	print_gp_filter_data(finfo, obs, fx, fp);
	fputs("e\n", fp);
     } else if (finfo->graph_opt & FILTER_GRAPH_TREND) {
	if (zkeypos == 'L') {
	    print_keypos_string(GP_KEY_LEFT | GP_KEY_TOP, fp);
	}
	xtitle = g_strdup_printf(_("%s (original data)"), finfo->vname);
	ztitle = g_strdup_printf(_("%s (smoothed)"), finfo->vname);
	fprintf(fp, "plot '-' using 1:2 title '%s' w lines, \\\n"
		" '-' using 1:2 title '%s' w lines\n", xtitle, ztitle);
	print_gp_filter_data(finfo, obs, dataset->Z[v], fp);
	fputs("e , \\\n", fp);
	print_gp_filter_data(finfo, obs, fx, fp);
	fputs("e\n", fp);
    } else if (finfo->graph_opt & FILTER_GRAPH_CYCLE) {
	if (finfo->ftype == FILTER_BK) {
	    title =
		g_strdup_printf(_("Baxter-King component of %s at frequency %d to %d"),
				finfo->vname, finfo->bkl, finfo->bku);
	} else {
	    title = g_strdup_printf(_("Cyclical component of %s"), finfo->vname);
	}
	fprintf(fp, "set title '%s'\n", title);
	fputs("set xzeroaxis\n", fp);
	fprintf(fp, "plot '-' using 1:2 title '' w lines\n");
	print_gp_filter_data(finfo, obs, u, fp);
	fputs("e\n", fp);
    }

    gretl_pop_c_numeric_locale();

    g_free(xtitle);
    g_free(ztitle);
    g_free(title);

    err = finalize_plot_input_file(fp);
    gui_graph_handler(err);

    return err;
}

static void butterworth_poles_graph (GtkWidget *button, filter_info *finfo)
{
    const char *fmt = N_("Butterworth poles (n = %d, cutoff = %d degrees)");
    gchar *title;
    FILE *fp;
    double cut, theta, delta;
    double x, y, px, py;
    int i, n = finfo->order;
    int err = 0;

    fp = open_plot_input_file(PLOT_ROOTS, 0, &err);
    if (err) {
	gui_errmsg(err);
	return;
    }

    cut = finfo->cutoff * M_PI / 180;
    cut = tan(cut / 2);

    title = g_strdup_printf(_(fmt), finfo->order, finfo->cutoff);
    fprintf(fp, "set title \"%s\"\n", title);
    g_free(title);

    fputs("unset border\n", fp);
    fputs("unset key\n", fp);
    fputs("set xzeroaxis\n", fp);
    fputs("set yzeroaxis\n", fp);
    fputs("unset xtics\n", fp);
    fputs("unset ytics\n", fp);
    fputs("set size square\n", fp);
    fputs("set polar\n", fp);
    fputs("plot 1.0 w lines, \\\n"
	  "'-' w points pt 7\n", fp);

    gretl_push_c_numeric_locale();

    for (i=1; i<=n; i++) {
	theta = M_PI * (2 * i - 1) / (2 * n);
	delta = 1 + cut * (2 * sin(theta) + cut);
	x = (1 - cut * cut) / delta;
	y = 2 * cut * cos(theta) / delta;
	px = atan2(y, x);
	py = sqrt(x * x + y * y);
	fprintf(fp, "%g %g # %.4f,%.4f\n", px, py, x, y);
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    err = finalize_plot_input_file(fp);
    gui_graph_handler(err);
}

static void weights_shape_graph (GtkWidget *button, filter_info *finfo)
{
    gretl_matrix *W;
    int T = sample_size(dataset);
    int err = 0;

    W = gretl_matrix_alloc(T, 2);

    if (W == NULL) {
	err = E_ALLOC;
    } else {
	int t, list[3] = {2, 1, 2};

	err = private_matrix_add(W, "$Weights");

	if (!err) {
	    err = umatrix_set_names_from_string(W, "weight t", 0);
	}

	if (!err) {
	    char rstr[64];

	    for (t=0; t<T; t++) {
		W->val[T+t] = t;
	    }

	    finfo->lambda = gtk_spin_button_get_value(GTK_SPIN_BUTTON(finfo->spin1));
	    finfo->midfrac = gtk_spin_button_get_value(GTK_SPIN_BUTTON(finfo->spin2));
	    poly_weights(W->val, T, finfo->lambda, finfo->midfrac, finfo->wopt);
	    gretl_push_c_numeric_locale();
	    sprintf(rstr, "set xrange [0:%d] ; set yrange [0.8:%g] ;", T - 1,
		    1.1 * finfo->lambda);
	    gretl_pop_c_numeric_locale();
	    err = matrix_plot(W, list, rstr, OPT_O);
	}

	destroy_private_matrices();
    }

    gui_graph_handler(err);
}

static int do_filter_response_graph (filter_info *finfo)
{
    gretl_matrix *G;
    FILE *fp = NULL;
    gchar *title = NULL;
    double omega_star = 0.0;
    int i, nlit, err = 0;

    if (finfo->ftype == FILTER_BW) {
	omega_star = finfo->cutoff * M_PI / 180;
	G = butterworth_gain(finfo->order, omega_star, 0);
	nlit = 2;
    } else {
	G = hp_gain(finfo->lambda, 0);
	nlit = 1;
    }

    if (G == NULL) {
	gui_errmsg(E_ALLOC);
	return E_ALLOC;
    }

    fp = open_plot_input_file(PLOT_REGULAR, 0, &err);
    if (err) {
	gui_errmsg(err);
	gretl_matrix_free(G);
	return err;
    }

    fputs("set xrange [0:3.1416]\n", fp);
    fputs("set yrange [0:1.1]\n", fp);

    if (finfo->ftype == FILTER_BW) {
	title = g_strdup_printf("%s (n = %d, %s %.2fπ)", _("Gain for Butterworth filter"),
				finfo->order, _("nominal cutoff"),
				(double) finfo->cutoff / 180);
    } else {
	title = g_strdup_printf(_("Gain for H-P filter (lambda = %g)"), finfo->lambda);
    }
    fprintf(fp, "set title \"%s\"\n", title);
    g_free(title);

    fprintf(fp, "# literal lines = %d\n", nlit);
    fputs("set xtics (\"0\" 0, \"π/4\" pi/4, \"π/2\" pi/2, \"3π/4\" 3*pi/4, "
	  "\"π\" pi)\n", fp);

    gretl_push_c_numeric_locale();

    if (finfo->ftype == FILTER_BW) {
	fprintf(fp, "set arrow from %g,0 to %g,1.1 nohead\n",
		omega_star, omega_star);
    }

    fputs("plot \\\n", fp);
    fputs("'-' using 1:2 notitle w lines\n", fp);
    for (i=0; i<G->rows; i++) {
	fprintf(fp, "%.5f %.5f\n", gretl_matrix_get(G, i, 0),
		gretl_matrix_get(G, i, 1));
    }
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    gretl_matrix_free(G);

    err = finalize_plot_input_file(fp);
    gui_graph_handler(err);

    return err;
}

static int save_filtered_var (filter_info *finfo, double *x, int i,
			      int *saved)
{
    const char *vname = (i == FILTER_SAVE_TREND)? finfo->save_t :
	finfo->save_c;
    int v = series_index(dataset, vname);
    int err = 0;

    if (v == dataset->v) {
	err = dataset_add_allocated_series(dataset, x);
    } else {
	free(dataset->Z[v]);
	dataset->Z[v] = x;
	series_set_discrete(dataset, v, 0);
    }

    if (!err) {
	strcpy(dataset->varname[v], vname);
	filter_make_varlabel(finfo, v, i);
	*saved = 1;
    }

    return err;
}

static void record_filter_command (filter_info *finfo)
{
    int trend = finfo->save_opt & FILTER_SAVE_TREND;
    int cycle = finfo->save_opt & FILTER_SAVE_CYCLE;
    GString *s = g_string_new(NULL);

    if (finfo->ftype == FILTER_BK) {
	g_string_append_printf(s, "series %s = ", finfo->save_c);
    } else if (finfo->ftype == FILTER_HP) {
	if (trend) {
	    g_string_append_printf(s, "series %s = %s - ", finfo->save_t, finfo->vname);
	} else if (cycle) {
	    g_string_append_printf(s, "series %s = ", finfo->save_c);
	}
    } else if (trend) {
	g_string_append_printf(s, "series %s = ", finfo->save_t);
    } else if (cycle) {
	g_string_append_printf(s, "series %s = %s - ", finfo->save_c, finfo->vname);
    }

    gretl_push_c_numeric_locale();

    if (finfo->ftype == FILTER_SMA) {
	g_string_append_printf(s, "movavg(%s, %d, %d)\n", finfo->vname,
                               finfo->k, finfo->center);
    } else if (finfo->ftype == FILTER_EMA) {
	g_string_append_printf(s, "movavg(%s, %g, %d)\n", finfo->vname,
                               finfo->lambda, finfo->k);
    } else if (finfo->ftype == FILTER_HP) {
	g_string_append_printf(s, "hpfilt(%s, %g, %d)\n", finfo->vname,
                               finfo->lambda, finfo->oneside);
    } else if (finfo->ftype == FILTER_BK) {
	g_string_append_printf(s, "bkfilt(%s, %d, %d, %d)\n", finfo->vname,
                               finfo->bkl, finfo->bku, finfo->bkk);
    } else if (finfo->ftype == FILTER_BW) {
	g_string_append_printf(s, "bwfilt(%s, %d, %d)\n", finfo->vname,
                               finfo->order, finfo->cutoff);
    } else if (finfo->ftype == FILTER_POLY) {
	g_string_append_printf(s, "polyfit(%s, %d)\n", finfo->vname, finfo->order);
    } else if (finfo->ftype == FILTER_FD) {
	g_string_append_printf(s, "fracdiff(%s, %g)\n", finfo->vname, finfo->lambda);
    }

    gretl_pop_c_numeric_locale();

    if (trend) {
	g_string_append_printf(s, "setinfo %s --description=\"%s\"\n",
                               finfo->save_t, finfo->label_t);
    } else if (cycle) {
	g_string_append_printf(s, "setinfo %s --description=\"%s\"\n",
                               finfo->save_c, finfo->label_c);
    }

    if (trend && cycle) {
	g_string_append_printf(s, "series %s = %s - %s\n",
                               finfo->save_c, finfo->vname, finfo->save_t);
	g_string_append_printf(s, "setinfo %s --description=\"%s\"\n",
                               finfo->save_c, finfo->label_c);
    }

    add_command_to_stack(s->str, 0);
    g_string_free(s, TRUE);
}

static int calculate_filter (filter_info *finfo)
{
    int save_t1 = dataset->t1;
    int save_t2 = dataset->t2;
    const double *x = dataset->Z[finfo->vnum];
    double *fx = NULL;
    double *u = NULL;
    int saved = 0;
    int t, err = 0;

    if (finfo->ftype != FILTER_BK) {
	fx = malloc(dataset->n * sizeof *fx);
	if (fx == NULL) {
	    return E_ALLOC;
	}
	for (t=0; t<dataset->n; t++) {
	    fx[t] = NADBL;
	}
    }

    if ((finfo->graph_opt & FILTER_GRAPH_CYCLE) ||
	finfo->save_opt & FILTER_SAVE_CYCLE) {
	u = malloc(dataset->n * sizeof *u);
	if (u == NULL) {
	    free(fx);
	    return E_ALLOC;
	}
    }

    dataset->t1 = finfo->t1;
    dataset->t2 = finfo->t2;

    if (finfo->ftype == FILTER_SMA) {
	/* simple moving average */
	movavg_series(x, fx, dataset, finfo->k, finfo->center);
    } else if (finfo->ftype == FILTER_EMA) {
	/* exponential moving average */
	exponential_movavg_series(x, fx, dataset, finfo->lambda,
				  finfo->k, finfo->fx0);
    } else if (finfo->ftype == FILTER_HP) {
	/* Hodrick-Prescott */
	if (finfo->oneside) {
	    err = oshp_filter(x, fx, dataset, finfo->lambda, OPT_T);
	} else {
	    err = hp_filter(x, fx, dataset, finfo->lambda, OPT_T);
	}
    } else if (finfo->ftype == FILTER_BK) {
	/* Baxter and King bandpass */
	err = bkbp_filter(x, u, dataset, finfo->bkl, finfo->bku, finfo->bkk);
    } else if (finfo->ftype == FILTER_BW) {
	/* Butterworth */
	err = butterworth_filter(x, fx, dataset, finfo->order, finfo->cutoff);
    } else if (finfo->ftype == FILTER_POLY) {
	if (finfo->wopt != OPT_NONE) {
	    err = weighted_poly_trend(x, fx, dataset, finfo->order,
				      finfo->wopt, finfo->lambda,
				      finfo->midfrac);
	} else {
	    err = poly_trend(x, fx, dataset, finfo->order);
	}
    } else if (finfo->ftype == FILTER_FD) {
	/* fractional differencing */
	err = fracdiff_series(x, fx, finfo->lambda, 1, -1, dataset);
    }

    dataset->t1 = save_t1;
    dataset->t2 = save_t2;

    if (!err && fx != NULL && u != NULL) {
	for (t=0; t<dataset->n; t++) {
	    if (na(x[t]) || na(fx[t])) {
		u[t] = NADBL;
	    } else {
		u[t] = x[t] - fx[t];
	    }
	}
    }

    if ((finfo->graph_opt & FILTER_GRAPH_TREND) ||
	(finfo->graph_opt & FILTER_GRAPH_CYCLE)) {
	do_filter_graph(finfo, fx, u);
    }

    if (finfo->graph_opt & FILTER_GRAPH_RESP) {
	do_filter_response_graph(finfo);
    }

    if (finfo->save_opt & FILTER_SAVE_TREND) {
	err = save_filtered_var(finfo, fx, FILTER_SAVE_TREND, &saved);
    } else {
	free(fx);
    }

    if (!err && finfo->save_opt & FILTER_SAVE_CYCLE) {
	err = save_filtered_var(finfo, u, FILTER_SAVE_CYCLE, &saved);
    } else {
	free(u);
    }

    if (saved) {
	populate_varlist();
	mark_dataset_as_modified();
    }

    if (!err && finfo->save_opt) {
	record_filter_command(finfo);
    }

    return err;
}

static int filter_code (GtkAction *action)
{
    const gchar *s = gtk_action_get_name(action);

    if (!strcmp(s, "FilterSMA"))
	return FILTER_SMA;
    else if (!strcmp(s, "FilterEMA"))
	return FILTER_EMA;
    else if (!strcmp(s, "FilterHP"))
	return FILTER_HP;
    else if (!strcmp(s, "FilterBK"))
	return FILTER_BK;
    else if (!strcmp(s, "FilterBW"))
	return FILTER_BW;
    else if (!strcmp(s, "FilterPoly"))
	return FILTER_POLY;
    else if (!strcmp(s, "FilterFD"))
	return FILTER_FD;
    else
	return FILTER_SMA;
}

void filter_callback (GtkAction *action)
{
    filter_info finfo;
    int code = filter_code(action);
    int v = mdata_active_var();
    int t1 = dataset->t1;
    int t2 = dataset->t2;
    int err = 0;

    err = series_adjust_sample(dataset->Z[v], &t1, &t2);

    if (!err && t2 - t1 + 1 < 4) {
	err = E_TOOFEW;
    }

    if (err) {
	gui_errmsg(err);
    } else {
	filter_info_init(&finfo, code, v, t1, t2);
	filter_dialog(&finfo);
    }
}
