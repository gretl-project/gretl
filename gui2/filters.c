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
#include "dlgutils.h"
#include "filters.h"
#include "libset.h"

enum {
    FILTER_SAVE_NONE = 0,
    FILTER_SAVE_TREND = 1 << 0,
    FILTER_SAVE_CYCLE = 1 << 1
};

enum {
    FILTER_GRAPH_NONE = 0,
    FILTER_GRAPH_TREND = 1 << 0,
    FILTER_GRAPH_CYCLE = 1 << 1
};

typedef struct filter_info_ filter_info;

struct filter_info_ {
    int ftype;
    int vnum;
    int t1;
    int t2;
    int nobs;
    int center;
    double lambda;
    int k;
    int bkl;
    int bku;
    int graph_opt;
    int save_opt;
    char save1[VNAMELEN];
    char save2[VNAMELEN];
    GtkWidget *dlg;
    GtkWidget *entry1;
    GtkWidget *entry2;
    GtkWidget *nspin;
};

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
    } else {
	return "";
    }
}

static void filter_info_init (filter_info *finfo, int ftype, int v,
			      int t1, int t2)
{
    finfo->ftype = ftype;
    finfo->vnum = v;
    finfo->t1 = t1;
    finfo->t2 = t2;
    finfo->nobs = 0;
    finfo->center = 0;
    finfo->lambda = 0.0;
    finfo->k = 0;
    finfo->bkl = 0;
    finfo->bku = 0;

    if (ftype == FILTER_BK) {
	finfo->graph_opt = FILTER_GRAPH_CYCLE;
    } else {
	finfo->graph_opt = FILTER_GRAPH_TREND;
    }

    finfo->save_opt = FILTER_SAVE_NONE;

    *finfo->save1 = 0;
    *finfo->save2 = 0;

    finfo->dlg = NULL;
    finfo->entry1 = NULL;
    finfo->entry2 = NULL;
    finfo->nspin = NULL;

    if (ftype == FILTER_SMA) {
	finfo->center = 1;
	finfo->nobs = (datainfo->pd == 1 || datainfo->pd == 10)? 3 :
	    datainfo->pd;
    } else if (ftype == FILTER_EMA) {
	finfo->lambda = 0.1;
    } else if (ftype == FILTER_HP) {
	finfo->lambda = libset_get_double("hp_lambda");
	if (na(finfo->lambda)) {
	    finfo->lambda = 100 * datainfo->pd * datainfo->pd;
	}
    } else if (ftype == FILTER_BK) {
	finfo->k = get_bkbp_k(datainfo);
	get_bkbp_periods(datainfo, &finfo->bkl, &finfo->bku);
	if (finfo->bku <= finfo->bkl) {
	    finfo->bku = finfo->bkl + 1;
	}
    }
}

static void filter_make_savename (filter_info *finfo, int i)
{
    char *targ = (i == 0)? finfo->save1 : finfo->save2;
    int len;

    if (finfo->ftype == FILTER_SMA) {
	strcpy(targ, (i == 0)? "ma_" : "mc_");
    } else if (finfo->ftype == FILTER_EMA) {
	strcpy(targ, (i == 0)? "ema_" : "emc_");
    } else if (finfo->ftype == FILTER_HP) {
	strcpy(targ, (i == 0)? "hpt_" : "hp_");
    } else if (finfo->ftype == FILTER_BK) {
	strcpy(targ, "bk_");
    }

    len = strlen(targ);
    strncat(targ, datainfo->varname[finfo->vnum], VNAMELEN - len - 1);
}

static void filter_make_varlabel (filter_info *finfo, int v, int i)
{
    char *targ = VARLABEL(datainfo, v);

    if (finfo->ftype == FILTER_SMA) {
	if (i == FILTER_SAVE_TREND) {
	    sprintf(targ, "%s%d-period moving average of %s", 
		    finfo->center? "Centered " : "",
		    finfo->nobs, datainfo->varname[finfo->vnum]);
	} else {
	    sprintf(targ, "Residual from %s%d-period MA of %s", 
		    finfo->center? "Centered " : "",
		    finfo->nobs, datainfo->varname[finfo->vnum]);
	}
    } else if (finfo->ftype == FILTER_EMA) {
	if (i == FILTER_SAVE_TREND) {
	    sprintf(targ, "Exponential moving average of %s (current weight %g)",
		    datainfo->varname[finfo->vnum], 1.0 - finfo->lambda);
	} else {
	    sprintf(targ, "Residual from EMA of %s (current weight %g)",
		    datainfo->varname[finfo->vnum], 1.0 - finfo->lambda);
	}	    
    } else if (finfo->ftype == FILTER_HP) {
	if (i == FILTER_SAVE_TREND) {
	    sprintf(targ, "Filtered %s: Hodrick-Prescott trend (lambda = %g)", 
		    datainfo->varname[finfo->vnum], finfo->lambda);
	} else {
	    sprintf(targ, "Filtered %s: Hodrick-Prescott cycle (lambda = %g)", 
		    datainfo->varname[finfo->vnum], finfo->lambda);
	}	    
    } else if (finfo->ftype == FILTER_BK) {
	sprintf(targ, "Filtered %s: Baxter-King cycle", 
		datainfo->varname[finfo->vnum]);
    }
}

static void spinner_set_int (GtkWidget *s, int *n)
{
    *n = (int) gtk_spin_button_get_value(GTK_SPIN_BUTTON(s));
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

static void spinner_set_double (GtkWidget *b, double *a)
{
    *a = gtk_spin_button_get_value(GTK_SPIN_BUTTON(b));
}

static void sma_center_callback (GtkWidget *w, int *c)
{
    *c = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));
}

static void set_trend_graph (GtkWidget *w, int *opt)
{
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	*opt |= FILTER_GRAPH_TREND;
    } else {
	*opt &= ~FILTER_GRAPH_TREND;
    }
}

static void set_cycle_graph (GtkWidget *w, int *opt)
{
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	*opt |= FILTER_GRAPH_CYCLE;
    } else {
	*opt &= ~FILTER_GRAPH_CYCLE;
    }
}

static void set_trend_save (GtkWidget *w, int *opt)
{
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	*opt |= FILTER_SAVE_TREND;
    } else {
	*opt &= ~FILTER_SAVE_TREND;
    }
}

static void set_cycle_save (GtkWidget *w, int *opt)
{
    if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w))) {
	*opt |= FILTER_SAVE_CYCLE;
    } else {
	*opt &= ~FILTER_SAVE_CYCLE;
    }
}

static int varname_error (filter_info *finfo, int i) 
{
    const char *s = (i == 1)? finfo->save1 : finfo->save2;

    if (*s == 0) {
	errbox(_("Variable name is missing"));
	return 1;
    } else if (validate_varname(s)) {
	return 1;
    } else if (i == 2 && (finfo->save_opt & FILTER_SAVE_TREND)) {
	if (!strcmp(s, finfo->save1)) {
	    errbox(_("Conflicting variable names"));
	    return 1;
	}
    }

    return 0;
}

static void filter_dialog_ok (GtkWidget *w, filter_info *finfo)
{
    if (finfo->save_opt & FILTER_SAVE_TREND) {
	strcpy(finfo->save1, gtk_entry_get_text(GTK_ENTRY(finfo->entry1)));
	if (varname_error(finfo, 1)) {
	    return;
	}
    }

    if (finfo->save_opt & FILTER_SAVE_CYCLE) {
	strcpy(finfo->save2, gtk_entry_get_text(GTK_ENTRY(finfo->entry2)));
	if (varname_error(finfo, 2)) {
	    return;
	} 
    } 

    if (finfo->nspin != NULL) {
	if (GTK_WIDGET_SENSITIVE(finfo->nspin)) {
	    finfo->nobs = (int) 
		gtk_spin_button_get_value(GTK_SPIN_BUTTON(finfo->nspin));
	} else {
	    finfo->nobs = 0;
	}
    }
    

    gtk_widget_destroy(finfo->dlg);
}

static void filter_dialog_hsep (GtkWidget *dlg)
{
    GtkWidget *hs = gtk_hseparator_new();
    
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), hs, TRUE, TRUE, 0);
}

static void filter_graph_check_button (GtkWidget *dlg, filter_info *finfo,
				       const char *txt, int g)
{
    GtkWidget *hbox, *w;

    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_check_button_new_with_label(txt);
    g_signal_connect(G_OBJECT(w), "toggled", 
		     (g == FILTER_GRAPH_TREND)? G_CALLBACK(set_trend_graph) :
		     G_CALLBACK(set_cycle_graph), &finfo->graph_opt);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), finfo->graph_opt & g);
    gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 5);    
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), hbox, FALSE, FALSE, 0);
}

static void filter_save_check_buttons (GtkWidget *dlg, filter_info *finfo)
{
    GtkWidget *hbox, *tab, *w;
    int i;

    hbox = gtk_hbox_new(FALSE, 5);
    tab = gtk_table_new(2, 2, FALSE);
    gtk_table_set_col_spacing(GTK_TABLE(tab), 0, 5);

    for (i=0; i<2; i++) {
	w = gtk_check_button_new_with_label((i == 0)? _("Save smoothed series as") : 
					    _("Save cyclical component as"));
	g_signal_connect(G_OBJECT(w), "toggled", 
			 (i == 0)? G_CALLBACK(set_trend_save) :
			 G_CALLBACK(set_cycle_save), &finfo->save_opt);
	gtk_table_attach_defaults(GTK_TABLE(tab), w, 0, 1, i, i+1);

	w = gtk_entry_new();
	gtk_entry_set_max_length(GTK_ENTRY(w), VNAMELEN - 1);
	gtk_entry_set_width_chars(GTK_ENTRY(w), VNAMELEN + 3);
	filter_make_savename(finfo, i);
	gtk_entry_set_text(GTK_ENTRY(w), (i == 0)? finfo->save1 : finfo->save2);
	gtk_table_attach_defaults(GTK_TABLE(tab), w, 1, 2, i, i+1);
	if (i == 0) {
	    finfo->entry1 = w;
	} else {
	    finfo->entry2 = w;
	}
    }

    gtk_box_pack_start(GTK_BOX(hbox), tab, TRUE, TRUE, 5);  
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), hbox, FALSE, FALSE, 0);
}

static void bkbp_frequencies_table (GtkWidget *dlg, filter_info *finfo)
{
    GtkWidget *tab, *hbox, *w;
    GtkWidget *s1, *s2;

    hbox = gtk_hbox_new(FALSE, 5);
    tab = gtk_table_new(2, 2, FALSE);
    gtk_table_set_col_spacing(GTK_TABLE(tab), 0, 5);

    /* lower limit */
    w = gtk_label_new(_("Lower frequency bound:"));
    gtk_table_attach_defaults(GTK_TABLE(tab), w, 0, 1, 0, 1);
    s1 = gtk_spin_button_new_with_range(1, 32, 1);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(s1), finfo->bkl);
    g_signal_connect(G_OBJECT(s1), "value-changed",
		     G_CALLBACK(spinner_set_int), &finfo->bkl);
    gtk_table_attach_defaults(GTK_TABLE(tab), s1, 1, 2, 0, 1);
	
    /* upper limit */
    w = gtk_label_new(_("Upper frequency bound:"));
    gtk_table_attach_defaults(GTK_TABLE(tab), w, 0, 1, 1, 2);
    s2 = gtk_spin_button_new_with_range(4, 64, 1);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(s2), finfo->bku);
    g_signal_connect(G_OBJECT(s2), "value-changed",
		     G_CALLBACK(spinner_set_int), &finfo->bku);
    gtk_table_attach_defaults(GTK_TABLE(tab), s2, 1, 2, 1, 2);

    /* inter-connect the lower and upper spinners */
    g_signal_connect(G_OBJECT(s1), "value-changed",
		     G_CALLBACK(check_bk_limits1), s2);
    g_signal_connect(G_OBJECT(s2), "value-changed",
		     G_CALLBACK(check_bk_limits2), s1);

    /* pack everything */
    gtk_box_pack_start(GTK_BOX(hbox), tab, TRUE, TRUE, 5);  
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), hbox, FALSE, FALSE, 0);
}

static void use_submean (GtkWidget *w, GtkWidget *s)
{
    gtk_widget_set_sensitive(s, GTK_TOGGLE_BUTTON(w)->active);
}

static void ema_obs_radios (GtkWidget *dlg, filter_info *finfo)
{
    GtkWidget *tab, *hbox, *w;
    GtkWidget *r1, *r2;
    GSList *group;

    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_label_new(_("The first EMA value is"));
    gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), hbox, FALSE, FALSE, 0);

    hbox = gtk_hbox_new(FALSE, 5);
    tab = gtk_table_new(2, 2, FALSE);
    gtk_table_set_col_spacing(GTK_TABLE(tab), 0, 5);

    /* default */
    r1 = gtk_radio_button_new_with_label(NULL, _("the mean of the whole series"));
    gtk_table_attach_defaults(GTK_TABLE(tab), r1, 0, 1, 0, 1);

    /* alternate */
    group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(r1));
    r2 = gtk_radio_button_new_with_label(group, _("the mean of the first n observations"));
    gtk_table_attach_defaults(GTK_TABLE(tab), r2, 0, 1, 1, 2);
    w = gtk_spin_button_new_with_range(1, (finfo->t2 - finfo->t1 + 1) / 4, 1);
    gtk_table_attach_defaults(GTK_TABLE(tab), w, 1, 2, 1, 2);
    gtk_widget_set_sensitive(w, FALSE);
    finfo->nspin = w;

    /* connect */
    g_signal_connect(G_OBJECT(r2), "clicked", G_CALLBACK(use_submean), w);

    /* pack everything */
    gtk_box_pack_start(GTK_BOX(hbox), tab, TRUE, TRUE, 10);  
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), hbox, FALSE, FALSE, 0);
}

static int filter_dialog (filter_info *finfo)
{
    GtkWidget *dlg;
    GtkWidget *hbox;
    GtkWidget *w;
    int ret = 0;

    dlg = gretl_dialog_new(_("gretl: time-series filter"), mdata->w,
			   GRETL_DLG_BLOCK);
    finfo->dlg = dlg;

    /* box title */
    hbox = gtk_hbox_new(FALSE, 5);
    w = gtk_label_new(_(filter_get_title(finfo->ftype)));
    gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 5);    
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), hbox, FALSE, FALSE, 5);

    if (finfo->ftype == FILTER_SMA) {
	GtkWidget *nspin;
	int T = finfo->t2 - finfo->t1 + 1;
	int nmax = T / 2;

	/* select number of observations */
	hbox = gtk_hbox_new(FALSE, 5);
	w = gtk_label_new(_("Number of observations in average:"));
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	nspin = gtk_spin_button_new_with_range(2, (gdouble) nmax, 1);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(nspin), (gdouble) finfo->nobs);
	g_signal_connect(G_OBJECT(nspin), "value-changed",
			 G_CALLBACK(spinner_set_int), &finfo->nobs);
	gtk_box_pack_start(GTK_BOX(hbox), nspin, TRUE, TRUE, 5);    
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), hbox, FALSE, FALSE, 0);

	/* select centered or not */
	hbox = gtk_hbox_new(FALSE, 5);
	w = gtk_check_button_new_with_label(_("Centered"));
	if (finfo->center) {
	    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(w), TRUE);
	}
	g_signal_connect(G_OBJECT(w), "toggled", 
			 G_CALLBACK(sma_center_callback), &finfo->center);
	gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 5);    
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), hbox, FALSE, FALSE, 0);
    } else if (finfo->ftype == FILTER_EMA) {
	/* set weight on current observation */
	hbox = gtk_hbox_new(FALSE, 5);
	w = gtk_label_new(_("Weight on current observation:"));
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	w = gtk_spin_button_new_with_range(0.001, 0.999, 0.001);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(w), finfo->lambda);
	g_signal_connect(G_OBJECT(w), "value-changed",
			 G_CALLBACK(spinner_set_double), &finfo->lambda);
	gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 5);    
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), hbox, FALSE, FALSE, 0);

	/* set calculation of initial EMA value */
	ema_obs_radios(dlg, finfo);
    } else if (finfo->ftype == FILTER_HP) {
	/* set H-P "lambda" parameter */
	hbox = gtk_hbox_new(FALSE, 5);
	w = gtk_label_new(_("lambda (higher values -> smoother trend):"));
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	w = gtk_spin_button_new_with_range(50.0, 2.0 * finfo->lambda, 1);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(w), finfo->lambda);
	g_signal_connect(G_OBJECT(w), "value-changed",
			 G_CALLBACK(spinner_set_double), &finfo->lambda);
	gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 5);    
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), hbox, FALSE, FALSE, 0);
    } else if (finfo->ftype == FILTER_BK) {
	/* set "k" */
	hbox = gtk_hbox_new(FALSE, 5);
	w = gtk_label_new(_("k (higher values -> better approximation):"));
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	w = gtk_spin_button_new_with_range(1.0, 4 * finfo->k, 1);
	gtk_spin_button_set_value(GTK_SPIN_BUTTON(w), finfo->k);
	g_signal_connect(G_OBJECT(w), "value-changed",
			 G_CALLBACK(spinner_set_int), &finfo->k);
	gtk_box_pack_start(GTK_BOX(hbox), w, TRUE, TRUE, 5);    
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), hbox, FALSE, FALSE, 0);

	/* set periods */
	bkbp_frequencies_table(dlg, finfo);
    }

    filter_dialog_hsep(dlg);

    if (finfo->ftype == FILTER_BK) {
	/* graphical output? */
	filter_graph_check_button(dlg, finfo, _("Graph filtered series"),
				  FILTER_GRAPH_CYCLE);
	/* save to dataset? */
	hbox = gtk_hbox_new(FALSE, 5);
	w = gtk_check_button_new_with_label(_("Save cyclical series as"));
	g_signal_connect(G_OBJECT(w), "toggled", 
			 G_CALLBACK(set_cycle_save), &finfo->save_opt);
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	w = gtk_entry_new();
	gtk_entry_set_max_length(GTK_ENTRY(w), VNAMELEN - 1);
	gtk_entry_set_width_chars(GTK_ENTRY(w), VNAMELEN + 3);
	filter_make_savename(finfo, 1);
	gtk_entry_set_text(GTK_ENTRY(w), finfo->save2);
	gtk_box_pack_start(GTK_BOX(hbox), w, FALSE, FALSE, 5);
	finfo->entry2 = w;
	gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dlg)->vbox), hbox, FALSE, FALSE, 0);
    } else {
	/* graphical output? */
	filter_graph_check_button(dlg, finfo, _("Graph original and smoothed series"),
				  FILTER_GRAPH_TREND);
	filter_graph_check_button(dlg, finfo, _("Graph residual or cycle series"),
				  FILTER_GRAPH_CYCLE);
	/* save to dataset? */
	filter_dialog_hsep(dlg);
	filter_save_check_buttons(dlg, finfo);
    }

    /* convenience pointer for button box */
    hbox = GTK_DIALOG(dlg)->action_area;

    /* Cancel button */
    w = cancel_delete_button(hbox, dlg, &ret);

    /* "OK" button */
    w = ok_button(hbox);
    g_signal_connect(G_OBJECT(w), "clicked", 
		     G_CALLBACK(filter_dialog_ok), finfo);
    gtk_widget_grab_default(w);

    gtk_widget_show_all(dlg);

    return ret;
}

static void print_gp_data (filter_info *finfo, const double *obs,
			   const double *x, FILE *fp)
{
    int t;

    for (t=finfo->t1; t<=finfo->t2; t++) {
	if (na(x[t])) {
	    fprintf(fp, "%g ?\n", obs[t]);
	} else {
	    fprintf(fp, "%g %g\n", obs[t], x[t]);
	}
    }
}

static int 
do_filter_graph (filter_info *finfo, const double *fx, const double *u)
{
    int twoplot = 0;
    int zkeypos = 'R';
    FILE *fp = NULL;
    const double *obs;
    char xtitle[48];
    char ztitle[48];
    char title[128];
    int v;

    obs = gretl_plotx(datainfo);
    if (obs == NULL) {
	return E_ALLOC;
    }

    if ((finfo->graph_opt & FILTER_GRAPH_TREND) &&
	(finfo->graph_opt & FILTER_GRAPH_CYCLE)) {
	twoplot = 1;
    }

    if (gnuplot_init((twoplot)? PLOT_TRI_GRAPH : PLOT_REGULAR, &fp)) { 
	return E_FOPEN;
    }

    if (datainfo->pd == 4) {
	if ((finfo->t2 - finfo->t1) / 4 < 8) {
	    fputs("set xtics nomirror 0,1\n", fp); 
	    fputs("set mxtics 4\n", fp);
	}
    } else if (datainfo->pd == 12) {
	if ((finfo->t2 - finfo->t1) / 12 < 8) {
	    fputs("set xtics nomirror 0,1\n", fp); 
	    fputs("set mxtics 12\n", fp);
	}
    }

    v = finfo->vnum;

    if (finfo->graph_opt & FILTER_GRAPH_TREND) {
	if (Z[v][finfo->t2] > Z[v][finfo->t1]) {
	    zkeypos = 'L';
	}
    }

    gretl_push_c_numeric_locale();

    if (twoplot) {
	fputs("set size 1.0,1.0\nset multiplot\nset size 1.0,0.60\n", fp);

	fputs("set origin 0.0,0.4\n", fp);
	if (zkeypos == 'L') {
	    fputs("set key top left\n", fp);
	}
	sprintf(xtitle, I_("%s (original data)"), datainfo->varname[v]);
	sprintf(ztitle, I_("%s (smoothed)"), datainfo->varname[v]);
	fprintf(fp, "plot '-' using 1:2 title '%s' w lines, \\\n"
		" '-' using 1:2 title '%s' w lines\n", xtitle, ztitle);
	print_gp_data(finfo, obs, Z[v], fp);
	fputs("e , \\\n", fp);
	print_gp_data(finfo, obs, fx, fp);
	fputs("e\n", fp);

	fputs("set size 1.0,0.38\n", fp);
	fputs("set origin 0.0,0.0\n", fp);
	fputs("set xzeroaxis\n", fp);
	sprintf(title, I_("Cyclical component of %s"), datainfo->varname[v]);
	fprintf(fp, "plot '-' using 1:2 title '%s' w lines\n", title);
	print_gp_data(finfo, obs, u, fp);
	fputs("e\n", fp);
	fputs("set nomultiplot\n", fp);
    } else if (finfo->graph_opt & FILTER_GRAPH_TREND) {
	if (zkeypos == 'L') {
	    fputs("set key top left\n", fp);
	}
	sprintf(xtitle, I_("%s (original data)"), datainfo->varname[v]);
	sprintf(ztitle, I_("%s (smoothed)"), datainfo->varname[v]);
	fprintf(fp, "plot '-' using 1:2 title '%s' w lines, \\\n"
		" '-' using 1:2 title '%s' w lines\n", xtitle, ztitle);
	print_gp_data(finfo, obs, Z[v], fp);
	fputs("e , \\\n", fp);
	print_gp_data(finfo, obs, fx, fp);
	fputs("e\n", fp);
    } else if (finfo->graph_opt & FILTER_GRAPH_CYCLE) {
	if (finfo->ftype == FILTER_BK) {
	    sprintf(title, I_("Baxter-King component of %s at frequency %d to %d"), 
		    datainfo->varname[v], finfo->bkl, finfo->bku);
	} else {
	    sprintf(title, I_("Cyclical component of %s"), datainfo->varname[v]);
	}
	fprintf(fp, "set title '%s'\n", title); 
	fputs("set xzeroaxis\n", fp);
	fprintf(fp, "plot '-' using 1:2 title '' w lines\n");
	print_gp_data(finfo, obs, u, fp);
	fputs("e\n", fp);
    }	

    gretl_pop_c_numeric_locale();

    fclose(fp);

    if (gnuplot_make_graph()) {
	errbox(_("gnuplot command failed"));
	return 1;
    } 

    register_graph();

    return 0;
}

static int save_filtered_var (filter_info *finfo, double *x, int i,
			      int *saved)
{
    const char *vname = (i == FILTER_SAVE_TREND)? finfo->save1 :
	finfo->save2;
    int v = varindex(datainfo, vname);
    int err = 0;

    if (v == datainfo->v) {
	err = dataset_add_allocated_series(x, &Z, datainfo);
    } else {
	free(Z[v]);
	Z[v] = x;
	set_var_scalar(datainfo, v, 0);
	set_var_discrete(datainfo, v, 0);
    }

    if (!err) {
	strcpy(datainfo->varname[v], vname);
	filter_make_varlabel(finfo, v, i);
	*saved = 1;
    }

    return err;
}

/* centered MA for even number of terms in average */

static int sma_special (const filter_info *finfo, double *fx,
			const double *x)
{
    int offset = finfo->nobs / 2;
    double *tmp;
    int t, t1, t2;
    int i, k, n;

    t1 = finfo->t1 + offset;
    t2 = finfo->t2 - offset + 1;

    n = t2 - t1 + 1;

    tmp = malloc(n * sizeof *tmp);
    if (tmp == NULL) {
	return E_ALLOC;
    }

    for (t=t1; t<=t2; t++) {
	fx[t] = 0.0;
	for (i=0; i<finfo->nobs; i++) {
	    k = i - offset + 1;
	    fx[t] += x[t-k];
	}
	fx[t] /= finfo->nobs;
    }

    for (t=0; t<n; t++) {
	tmp[t] = fx[t+t1];
    }

    for (t=t1; t<t2; t++) {
	k = t - t1;
	fx[t] = (tmp[k] + tmp[k+1]) / 2.0;
    }
    fx[t2] = NADBL;

    free(tmp);

    return 0;
}

static int calculate_filter (filter_info *finfo)
{
    const double *x = Z[finfo->vnum];
    double *fx = NULL;
    double *u = NULL;
    int saved = 0;
    int t, err = 0;

    if (finfo->ftype != FILTER_BK) {
	fx = malloc(datainfo->n * sizeof *fx);
	if (fx == NULL) {
	    return E_ALLOC;
	}
	for (t=0; t<datainfo->n; t++) {
	    fx[t] = NADBL;
	}
    }

    if ((finfo->graph_opt & FILTER_GRAPH_CYCLE) ||
	finfo->save_opt & FILTER_SAVE_CYCLE) {
	u = malloc(datainfo->n * sizeof *u);
	if (u == NULL) {
	    free(fx);
	    return E_ALLOC;
	}
    }

    if (finfo->ftype == FILTER_SMA) {
	if (finfo->center && finfo->nobs % 2 == 0) {
	    err = sma_special(finfo, fx, x);
	} else {
	    int offset = finfo->center ? finfo->nobs / 2 : 0;
	    int i, k, t1, t2;

	    t1 = finfo->center ? finfo->t1 + offset : finfo->t1 + finfo->nobs - 1;
	    t2 = finfo->center ? finfo->t2 - offset : finfo->t2;

	    for (t=t1; t<=t2; t++) {
		fx[t] = 0.0;
		for (i=0; i<finfo->nobs; i++) {
		    k = i - offset;
		    fx[t] += x[t-k];
		}
		fx[t] /= finfo->nobs;
	    }
	}
    } else if (finfo->ftype == FILTER_EMA) {
	int t1 = finfo->t1 + 1;

	if (finfo->nobs == 0) {
	    fx[finfo->t1] = gretl_mean(finfo->t1, finfo->t2, x);
	} else {
	    t1 = finfo->t1 + finfo->nobs - 1;
	    fx[t1] = gretl_mean(finfo->t1, t1, x);
	    t1++;
	}
	for (t=t1; t<=finfo->t2; t++) {
	    fx[t] = finfo->lambda * x[t] + (1.0 - finfo->lambda) * fx[t-1];
	}
    } else if (finfo->ftype == FILTER_HP) {
	double l = libset_get_double("hp_lambda");

	set_hp_lambda(finfo->lambda);
	err = hp_filter(x, fx, datainfo, OPT_T);
	set_hp_lambda(l);
    } else if (finfo->ftype == FILTER_BK) {
	set_bkbp_k(finfo->k);
	set_bkbp_periods(finfo->bkl, finfo->bku);
	err = bkbp_filter(x, u, datainfo);
	unset_bkbp_k();
	unset_bkbp_periods();
    }
	
    if (!err && fx != NULL && u != NULL) {
	for (t=0; t<datainfo->n; t++) {
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
    }

    return err;
}

void filter_callback (gpointer p, guint code, GtkWidget *w)
{
    filter_info finfo;
    int v = mdata_active_var();
    int t1 = datainfo->t1;
    int t2 = datainfo->t2;
    int cancel = 0;
    int err = 0;

    if (var_is_scalar(datainfo, v)) {
	errbox(_("This variable is a scalar"));
	return;
    }

    err = array_adjust_t1t2(Z[v], &t1, &t2);
    if (err) {
	gui_errmsg(E_MISSDATA);
	return;
    }

    if (t2 - t1 + 1 < 4) {
	gui_errmsg(E_MISSDATA);
	return;
    }	

    filter_info_init(&finfo, code, v, t1, t2);

    cancel = filter_dialog(&finfo);
    if (cancel) {
	return;
    }

    if (finfo.graph_opt == FILTER_GRAPH_NONE &&
	finfo.save_opt == FILTER_SAVE_NONE) {
	return;
    }

    err = calculate_filter(&finfo);
    if (err) {
	gui_errmsg(err);
    }    
}
