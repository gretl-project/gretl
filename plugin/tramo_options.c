/*
 *  Copyright (c) by Allin Cottrell
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */

#include "libgretl.h"
#include <gtk/gtk.h>
#include "tramo_x12a.h"

typedef struct _tramo_options tramo_options;

struct _tramo_options {
    int iatip;    /* correction for outliers? 0=no, 1=auto */
    int aio;      /* sorts of outliers recognized (codes 1, 2, or 3) */
    float va;     /* critical value for outliers (or 0.0 for auto) */
    GtkWidget *iatip_button;
    GtkWidget *aio_transitory_button;
    GtkWidget *aio_shift_button;
    GtkWidget *va_button;
    GtkWidget *va_spinner;
    GtkWidget *aio_label;
    GtkWidget *va_label;

    int lam;      /* log transformation? (0=yes, 1=no, -1=auto) FIXME */
    int imean;    /* mean correction? (0=yes, 1=no) FIXME */

    int d;        /* no. of non-seasonal differences */
    int bd;       /* no. of seasonal differences */
    int p;        /* no of non-seasonal AR terms */
    int bp;       /* no of seasonal AR terms */
    int q;        /* no of non-seasonal MA terms */
    int bq;       /* no of seasonal MA terms */
};

#define option_widgets_shown(p)  (p->va_spinner != NULL)

static void tramo_options_set_defaults (tramo_options *opts)
{
    opts->iatip = 1;         /* detect outliers */
    opts->aio = 2;           /* both transitory changes and level shifts */
    opts->va = 0.0;          /* let critical value be decided by tramo */
    opts->lam = -1;          /* leave log/level decision to tramo */
    opts->imean = 1;         /* no mean correction */
    opts->d = opts->bd = 1;  
    opts->p = opts->bp = 0;
    opts->q = opts->bq = 1;
}

static void va_spinner_set_state (tramo_options *opts)
{
    if (!option_widgets_shown(opts)) return;
    
    gtk_widget_set_sensitive(opts->va_spinner, 
			     GTK_WIDGET_IS_SENSITIVE(opts->va_label) &&
			     !GTK_TOGGLE_BUTTON(opts->va_button)->active);
}

static void outlier_options_set_sensitive (tramo_options *opts, gboolean s)
{
    if (!option_widgets_shown(opts)) return;

    gtk_widget_set_sensitive(opts->aio_label, s);
    gtk_widget_set_sensitive(opts->aio_transitory_button, s);
    gtk_widget_set_sensitive(opts->aio_shift_button, s);
    gtk_widget_set_sensitive(opts->va_label, s);
    gtk_widget_set_sensitive(opts->va_button, s);
    va_spinner_set_state(opts);
}

static void flip_iatip (GtkWidget *w, tramo_options *opts)
{
    if (!option_widgets_shown(opts)) return;
    
    if (GTK_TOGGLE_BUTTON(w)->active) {
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
  
    if (GTK_TOGGLE_BUTTON(w)->active) {
	gtk_widget_set_sensitive(opts->va_spinner, FALSE);
	opts->va = 0.0;
    } else {
	gtk_widget_set_sensitive(opts->va_spinner, TRUE);
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
    if (!GTK_TOGGLE_BUTTON(this_button)->active && 
	!GTK_TOGGLE_BUTTON(other_button)->active) {
	gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(other_button), TRUE);
    }

    if (GTK_TOGGLE_BUTTON(opts->aio_transitory_button)->active) {
	if (GTK_TOGGLE_BUTTON(opts->aio_shift_button)->active) {
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

static void set_lam (GtkWidget *w, tramo_options *opts)
{
    opts->lam = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "lam_value"));
}

static void set_imean (GtkWidget *w, tramo_options *opts)
{
    opts->imean = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "imean_value"));
}

static GtkWidget *make_notebook_page_table (GtkWidget *notebook, 
					    const gchar *tab_title,
					    gint rows, gint cols)
{
    GtkWidget *box, *tmp, *tbl;

    box = gtk_vbox_new(FALSE, 0);
    gtk_container_set_border_width(GTK_CONTAINER(box), 10);
    gtk_widget_show(box);

    tmp = gtk_label_new(tab_title);
    gtk_widget_show(tmp);

    gtk_notebook_append_page(GTK_NOTEBOOK(notebook), box, tmp);  

    tbl = gtk_table_new(rows, cols, FALSE);
    gtk_table_set_row_spacings(GTK_TABLE(tbl), 5);
    gtk_table_set_col_spacings(GTK_TABLE(tbl), 5);
    gtk_box_pack_start(GTK_BOX(box), tbl, FALSE, FALSE, 0);
    gtk_widget_show(tbl); 

    return tbl;
}

static void tramo_tab_output (GtkWidget *notebook, tx_request *request)
{
    GtkWidget *tbl, *tmp;
    int tbl_len = 6, row = 0;

    tbl = make_notebook_page_table(notebook, _("Output"), tbl_len, 2);

    /* label pertaining saving series */
    tmp = gtk_label_new(_("Save to data set:"));
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 1, row, row + 1);
    row++;
    gtk_widget_show(tmp);

    tmp = gtk_check_button_new_with_label(_("Seasonally adjusted series"));
    gtk_widget_show(tmp);
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 1, row, row + 1);
    row++;
    request->opt[D11].check = tmp;
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), TRUE);

    tmp = gtk_check_button_new_with_label(_("Trend/cycle"));
    gtk_widget_show(tmp);
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 1, row, row + 1);
    row++;
    request->opt[D12].check = tmp;
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), FALSE);

    tmp = gtk_check_button_new_with_label(_("Irregular"));
    gtk_widget_show(tmp);
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 1, row, row + 1);
    row++;
    request->opt[D13].check = tmp;
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), FALSE);

    tmp = gtk_hseparator_new();
    gtk_widget_show(tmp);
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 1, row, row + 1);
    row++;
    
    tmp = gtk_check_button_new_with_label(_("Generate graph"));
    gtk_widget_show(tmp);
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 1, row, row + 1);
    request->opt[TRIGRAPH].check = tmp;
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), TRUE);
}

static void tramo_tab_transform (GtkWidget *notebook, tramo_options *opts)
{
    GtkWidget *b1, *b2, *b3, *tbl, *tmp;
    GSList *log_group = NULL, *mean_group = NULL;
    int tbl_len = 6, row = 0;

    tbl = make_notebook_page_table(notebook, _("Transformations"), tbl_len, 2);

    /* logs versus levels: radio group */

    /* logs option */
    b1 = gtk_radio_button_new_with_label(NULL, _("Log transformation"));
    log_group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    gtk_widget_show(b1);
    gtk_table_attach_defaults(GTK_TABLE(tbl), b1, 0, 2, row, row + 1);
    row++;
    g_signal_connect(G_OBJECT(b1), "clicked",
		     G_CALLBACK(set_lam), 
		     opts);
    g_object_set_data(G_OBJECT(b1), "lam_value", 
                      GINT_TO_POINTER(0));

    /* no logs option */
    b2 = gtk_radio_button_new_with_label(log_group, _("No log transformation"));
    log_group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b2));
    gtk_widget_show(b2);
    gtk_table_attach_defaults(GTK_TABLE(tbl), b2, 0, 2, row, row + 1);
    row++;
    g_signal_connect(G_OBJECT(b2), "clicked",
		     G_CALLBACK(set_lam), 
		     opts);
    g_object_set_data(G_OBJECT(b2), "lam_value", 
                      GINT_TO_POINTER(1));

    /* automatic log/level option */
    b3 = gtk_radio_button_new_with_label(log_group, _("Automatic"));
    log_group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b3));
    gtk_widget_show(b3);
    gtk_table_attach_defaults(GTK_TABLE(tbl), b3, 0, 2, row, row + 1);
    row++;
    g_signal_connect(G_OBJECT(b3), "clicked",
		     G_CALLBACK(set_lam), 
		     opts);
    g_object_set_data(G_OBJECT(b3), "lam_value", 
                      GINT_TO_POINTER(-1));

    switch (opts->lam) {
    case  0: gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b1), TRUE); break;
    case  1: gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b2), TRUE); break;
    case -1: gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b3), TRUE); break;
    default: break;
    }

    /* horizontal separator */
    tmp = gtk_hseparator_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;
    gtk_widget_show(tmp);

    /* mean correction: radio group */
    b1 = gtk_radio_button_new_with_label(NULL, _("Mean correction"));
    mean_group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b1));
    gtk_widget_show(b1);
    gtk_table_attach_defaults(GTK_TABLE(tbl), b1, 0, 2, row, row + 1);
    row++;    
    g_signal_connect(G_OBJECT(b1), "clicked",
		     G_CALLBACK(set_imean), 
		     opts);
    g_object_set_data(G_OBJECT(b1), "imean_value", 
                      GINT_TO_POINTER(0));

    b2 = gtk_radio_button_new_with_label(mean_group, _("No mean correction"));
    mean_group = gtk_radio_button_get_group(GTK_RADIO_BUTTON(b2));
    gtk_widget_show(b2);
    gtk_table_attach_defaults(GTK_TABLE(tbl), b2, 0, 2, row, row + 1);
    g_signal_connect(G_OBJECT(b2), "clicked",
		     G_CALLBACK(set_imean), 
		     opts);
    g_object_set_data(G_OBJECT(b2), "imean_value", 
                      GINT_TO_POINTER(1));

    switch (opts->imean) {
    case  0: gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b1), TRUE); break;
    case  1: gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(b2), TRUE); break;
    default: break;
    }
}

static void get_va_value (GtkSpinButton *sb, tramo_options *opts)
{
    opts->va = (float) gtk_spin_button_get_value(sb);
}

static void tramo_tab_outliers (GtkWidget *notebook, tramo_options *opts)
{
    GtkWidget *tbl, *tmp;
    GtkObject *adj;
    int tbl_len = 9, row = 0;

    tbl = make_notebook_page_table(notebook, _("Outliers"), tbl_len, 2);

    /* Overall choice: outlier correction or no? */
    tmp = gtk_check_button_new_with_label(_("Detect and correct for outliers"));
    opts->iatip_button = tmp;
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;
    gtk_widget_show(tmp);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), (opts->iatip != 0));
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(flip_iatip), 
		     opts);

    /* horizontal separator */
    tmp = gtk_hseparator_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;
    gtk_widget_show(tmp);
    
    /* label pertaining to transitory and level-shift buttons */
    tmp = gtk_label_new(_("Besides additive outliers, allow for:"));
    opts->aio_label = tmp;
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;
    gtk_widget_show(tmp);

    /* "transitory" button */
    tmp = gtk_check_button_new_with_label(_("transitory changes"));
    opts->aio_transitory_button = tmp;
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;
    gtk_widget_show(tmp);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), (opts->aio < 3));
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(tramo_aio_callback), 
		     opts);

    /* level-shift button */
    tmp = gtk_check_button_new_with_label(_("shifts of level"));
    opts->aio_shift_button = tmp;
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;
    gtk_widget_show(tmp);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), (opts->aio > 1));
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(tramo_aio_callback), 
		     opts);
    
    /* horizontal separator */
    tmp = gtk_hseparator_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;
    gtk_widget_show(tmp);

    /* label pertaining to critical value */
    tmp = gtk_label_new(_("Critical value for outliers:"));
    opts->va_label = tmp;
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;
    gtk_widget_show(tmp);

    /* check button for auto critical value */
    tmp = gtk_check_button_new_with_label(_("Automatic"));
    opts->va_button = tmp;
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;
    gtk_widget_show(tmp);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(tmp), (opts->va == 0.0));
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(flip_auto_va), 
		     opts);

    /* spinner for manual critical value */
    adj = gtk_adjustment_new((opts->va == 0.0)? 3.3 : opts->va, 
			     2.1, 6.0, 0.1, 1.0, 1.0);
    tmp = gtk_spin_button_new(GTK_ADJUSTMENT(adj), 0.1, 1);
    opts->va_spinner = tmp;
    gtk_table_attach(GTK_TABLE(tbl), tmp, 0, 1, row, row + 1,
		     0, 0, 0, 0);
    gtk_widget_show(tmp);
    gtk_widget_set_sensitive(tmp, (opts->va != 0.0));
    g_signal_connect(G_OBJECT(tmp), "value-changed",
		     G_CALLBACK(get_va_value), 
		     opts);
}

static void tramo_arima_callback (GtkWidget *w, gint *var)
{
    GtkWidget *entry = g_object_get_data(G_OBJECT(w), "entry");

    *var = atoi(gtk_entry_get_text(GTK_ENTRY(entry)));
}

static GtkWidget *make_labeled_combo (const gchar *label, 
				      GtkWidget *tbl, gint row,
				      GList *list, gint *var)
{
    GtkWidget *w;
    char numstr[2];

    w = gtk_label_new(label);
    gtk_label_set_justify(GTK_LABEL(w), GTK_JUSTIFY_RIGHT);
    gtk_table_attach(GTK_TABLE(tbl), w, 0, 1, row, row + 1,
		     0, 0, 0, 0);
    gtk_widget_show(w);
    w = gtk_combo_new();
    gtk_combo_set_popdown_strings(GTK_COMBO(w), list); 
    sprintf(numstr, "%d", *var);
    gtk_entry_set_text(GTK_ENTRY(GTK_COMBO(w)->entry), numstr);
    gtk_widget_set_size_request(w, 48, -1);
    gtk_table_attach(GTK_TABLE(tbl), w, 1, 2, row, row + 1,
		     0, 0, 0, 0);
    g_object_set_data(G_OBJECT(GTK_COMBO(w)->list), "entry", 
		      GTK_COMBO(w)->entry);
    g_signal_connect(G_OBJECT(GTK_COMBO(w)->list), "selection-changed",
		     G_CALLBACK(tramo_arima_callback), 
		     var);

    return w;
}

static void tramo_tab_arima (GtkWidget *notebook, tramo_options *opts)
{
    GtkWidget *tbl, *tmp;
    int i, tbl_len = 9, row = 0;
    GList *onelist = NULL, *twolist = NULL, *threelist = NULL;
    gchar *intvals[] = {
	"0", "1", "2", "3"
    };

    for (i=0; i<2; i++) {
	onelist = g_list_append(onelist, intvals[i]);
    }
    for (i=0; i<3; i++) {
	twolist = g_list_append(twolist, intvals[i]);
    }
    for (i=0; i<4; i++) {
	threelist = g_list_append(threelist, intvals[i]);
    }

    tbl = make_notebook_page_table(notebook, _("ARIMA"), tbl_len, 2);
    gtk_table_set_homogeneous(GTK_TABLE(tbl), FALSE);

    /* difference terms */
    tmp = make_labeled_combo(_("Non-seasonal differences:"), tbl, row,
			     twolist, &opts->d);
    gtk_widget_show(tmp);
    row++;
    tmp = make_labeled_combo(_("Seasonal differences:"), tbl, row,
			     onelist, &opts->bd);
    gtk_widget_show(tmp);
    row++;
	    
    /* horizontal separator */
    tmp = gtk_hseparator_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;
    gtk_widget_show(tmp);

    /* AR terms */
    tmp = make_labeled_combo(_("Non-seasonal AR terms:"), tbl, row,
			     threelist, &opts->p);
    gtk_widget_show(tmp);
    row++;
    tmp = make_labeled_combo(_("Seasonal AR terms:"), tbl, row,
			     onelist, &opts->bp);
    gtk_widget_show(tmp);
    row++;

    /* horizontal separator */
    tmp = gtk_hseparator_new();
    gtk_table_attach_defaults(GTK_TABLE(tbl), tmp, 0, 2, row, row + 1);
    row++;
    gtk_widget_show(tmp);

    /* MA terms */
    tmp = make_labeled_combo(_("Non-seasonal MA terms:"), tbl, row,
			     threelist, &opts->q);
    gtk_widget_show(tmp);
    row++;
    tmp = make_labeled_combo(_("Seasonal MA terms:"), tbl, row,
			     onelist, &opts->bq);
    gtk_widget_show(tmp);
}

static void real_show_tramo_options (tx_request *request, GtkWidget *vbox) 
{
    GtkWidget *notebook;

    notebook = gtk_notebook_new();
    gtk_box_pack_start(GTK_BOX(vbox), notebook, TRUE, TRUE, 0);
    gtk_widget_show(notebook);

    tramo_tab_output(notebook, request);
    tramo_tab_outliers(notebook, request->opts);
    tramo_tab_transform(notebook, request->opts);
    tramo_tab_arima(notebook, request->opts);
}

static tramo_options *tramo_options_new (void)
{
    tramo_options *opts;

    opts = malloc(sizeof *opts);
    if (opts == NULL) return NULL;

    tramo_options_set_defaults(opts);

    opts->iatip_button = NULL;
    opts->aio_transitory_button = NULL;
    opts->aio_shift_button = NULL;
    opts->va_button = NULL;
    opts->va_spinner = NULL;
    opts->aio_label = NULL;
    opts->va_label = NULL;    

    return opts;
}

int show_tramo_options (tx_request *request, GtkWidget *vbox)
{
    tramo_options *opts;

    opts = tramo_options_new();
    if (opts == NULL) return 1;

    request->opts = opts;

    real_show_tramo_options(request, vbox);

    return 0;
}

void print_tramo_options (tx_request *request, FILE *fp)
{
    tramo_options *opts;

    if (request->opts == NULL) return;

    opts = request->opts;

    fputs("$INPUT ", fp);
    fprintf(fp, "lam=%d,", opts->lam);
    fprintf(fp, "imean=%d,", opts->imean);
    fprintf(fp, "iatip=%d,", opts->iatip);
    if (opts->iatip == 1) {
	fprintf(fp, "aio=%d,", opts->aio);
	if (opts->va != 0.0) {
	    fprintf(fp, "va=%.1f,", opts->va);
	}
    }

#if 1
    fprintf(fp, "D=%d,BD=%d,", opts->d, opts->bd);
    fprintf(fp, "P=%d,BP=%d,", opts->p, opts->bp);
    fprintf(fp, "Q=%d,BQ=%d,", opts->q, opts->bq);
#endif

    fputs("noadmiss=1,seats=2,$\n", fp);

    free(opts);
    request->opts = NULL;
}


