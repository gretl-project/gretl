/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "libgretl.h"
#include "internal.h"

#include <gtk/gtk.h>

#if GTK_MAJOR_VERSION < 2
#define OLD_GTK
#endif

struct lmax_opt {
    GtkWidget *dlg;
    GtkWidget *entry;
    double *lmax;
};

static void lmax_opt_free (GtkWidget *w, struct lmax_opt *opt)
{
    free(opt);
    gtk_main_quit();
}

static void lmax_opt_finalize (GtkWidget *w, struct lmax_opt *opt)
{
    const gchar *numstr;
    char *test;
    double x;

    numstr = gtk_entry_get_text(GTK_ENTRY(opt->entry));
    x = strtod(numstr, &test);
    if (*test != 0 || x < 0.0) {
	gretl_errmsg_set(_("Invalid value for the maximum of the "
			   "dependent variable"));
	*opt->lmax = NADBL;
    } else {
	*opt->lmax = x;
    }

    gtk_widget_destroy(opt->dlg);
}

static void lmax_opt_cancel (GtkWidget *w, struct lmax_opt *opt)
{
    *opt->lmax = 0.0;
    gtk_widget_destroy(opt->dlg);
}

static void lmax_dialog (double *lmax)
{
    GtkWidget *tmp, *hbox;
    gchar *numstr;
    struct lmax_opt *opt;

    opt = malloc(sizeof *opt);
    if (opt == NULL) return;

    opt->dlg = gtk_dialog_new();
    opt->lmax = lmax;

    gtk_window_set_title(GTK_WINDOW(opt->dlg), _("Logistic model"));
    gtk_container_set_border_width (GTK_CONTAINER 
				    (GTK_DIALOG (opt->dlg)->vbox), 10);
    gtk_container_set_border_width (GTK_CONTAINER 
				    (GTK_DIALOG (opt->dlg)->action_area), 5); 
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (opt->dlg)->vbox), 5);
    gtk_window_set_position (GTK_WINDOW (opt->dlg), GTK_WIN_POS_MOUSE);

#ifndef OLD_GTK
    g_signal_connect (G_OBJECT(opt->dlg), "destroy", 
		      G_CALLBACK(lmax_opt_free), opt);
#else
    gtk_signal_connect (GTK_OBJECT(opt->dlg), "destroy", 
			GTK_SIGNAL_FUNC(lmax_opt_free), opt);
#endif

    /* label */
    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new (_("Maximum (asymptote) for the\n"
			   "dependent variable"));
    gtk_label_set_justify(GTK_LABEL(tmp), GTK_JUSTIFY_CENTER);
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(opt->dlg)->vbox), 
		       hbox, FALSE, FALSE, 5);

    /* lmax entry */
    hbox = gtk_hbox_new(FALSE, 5);
    opt->entry = gtk_entry_new();
    gtk_entry_set_max_length(GTK_ENTRY(opt->entry), 12);
    numstr = g_strdup_printf("%g", *lmax);
    gtk_entry_set_text(GTK_ENTRY(opt->entry), numstr);
    g_free(numstr);
#ifndef OLD_GTK
    gtk_entry_set_width_chars(GTK_ENTRY(opt->entry), 6);
    gtk_entry_set_activates_default(GTK_ENTRY(opt->entry), TRUE);
#endif
    gtk_editable_select_region(GTK_EDITABLE(opt->entry), 0, -1);
    gtk_box_pack_start(GTK_BOX(hbox), opt->entry, TRUE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(opt->dlg)->vbox), 
		       hbox, FALSE, FALSE, 5);

    /* Create the "OK" button */
#ifndef OLD_GTK
    tmp = gtk_button_new_from_stock(GTK_STOCK_OK);
#else
    tmp = gtk_button_new_with_label(_("OK"));
#endif
    GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX(GTK_DIALOG (opt->dlg)->action_area), 
			tmp, TRUE, TRUE, 0);
    gtk_widget_grab_default (tmp);
#ifndef OLD_GTK
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(lmax_opt_finalize), opt);
#else
    gtk_signal_connect(GTK_OBJECT(tmp), "clicked",
		       GTK_SIGNAL_FUNC(lmax_opt_finalize), opt);
#endif

    /* And a Cancel button */
#ifndef OLD_GTK
    tmp = gtk_button_new_from_stock(GTK_STOCK_CANCEL);
#else
    tmp = gtk_button_new_with_label(_("Cancel"));
#endif
    gtk_box_pack_start (GTK_BOX(GTK_DIALOG(opt->dlg)->action_area), 
			tmp, TRUE, TRUE, 0);
#ifndef OLD_GTK
    g_signal_connect (G_OBJECT (tmp), "clicked", 
		      G_CALLBACK(lmax_opt_cancel), opt);
#else
    gtk_signal_connect(GTK_OBJECT(tmp), "clicked",
		       GTK_SIGNAL_FUNC(lmax_opt_cancel), opt);
#endif

    gtk_widget_show_all(opt->dlg);

    gtk_main();
}

static double get_lmax (const double *y, const DATAINFO *pdinfo,
			const char *lmstr)
{
    double lmax, ymax = 0.0;
    int t, lmax_auto = 1;

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (y[t] <= 0.0) {
	    gretl_errmsg_set(_("Illegal non-positive value of the "
			       "dependent variable"));
	    return NADBL;
	}
	if (y[t] > ymax) {
	    ymax = y[t];
	}
    }

    if (lmstr != NULL && *lmstr != '\0') {
	lmax = atof(lmstr + 5);
	lmax_auto = 0;
    }	

    if (lmax_auto) {
	if (ymax < 1.0) lmax = 1.0;
	else if (ymax < 100.0) lmax = 100.0;
	else lmax = 1.1 * ymax;
    } else if (lmax <= ymax) {
	gretl_errmsg_set(_("Invalid value for the maximum of the "
			   "dependent variable"));
	lmax = NADBL;
    }

    if (lmstr == NULL) {
	lmax_dialog(&lmax);
    }
	    
    return lmax;
}

static int make_logistic_depvar (double ***pZ, DATAINFO *pdinfo, 
				 int dv, double lmax)
{
    int t, v = pdinfo->v;
    int err;

    err = dataset_add_vars(1, pZ, pdinfo);
    if (err) return 1;

    for (t=0; t<pdinfo->n; t++) {
	double p = (*pZ)[dv][t];

	if (na(p)) continue;
	(*pZ)[v][t] = log(p / (lmax - p));
    }

    return 0;
}

static int rewrite_logistic_stats (const double **Z, const DATAINFO *pdinfo,
				   MODEL *pmod, int dv, double lmax)
{
    int t;
    double x;

    pmod->ybar = _esl_mean(pmod->t1, pmod->t2, Z[dv]);
    pmod->sdy = _esl_stddev(pmod->t1, pmod->t2, Z[dv]);

    /* make the VCV matrix before messing with the model stats */
    makevcv(pmod);

    pmod->ess = 0.0;
    for (t=0; t<pdinfo->n; t++) {
	x = pmod->yhat[t];
	if (na(x)) continue;
	pmod->yhat[t] = lmax / (1.0 + exp(-x));
	pmod->uhat[t] = Z[dv][t] - pmod->yhat[t];
	pmod->ess += pmod->uhat[t] * pmod->uhat[t];
    }

    pmod->sigma = sqrt(pmod->ess / pmod->dfd);

    pmod->tss = 0.0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	pmod->tss += (Z[dv][t] - pmod->ybar) * (Z[dv][t] - pmod->ybar);
    }

    pmod->fstt = pmod->dfd * (pmod->tss - pmod->ess) / (pmod->dfn * pmod->ess);

    pmod->rsq = pmod->adjrsq = NADBL;

    if (pmod->tss > 0) {
	pmod->rsq = 1.0 - (pmod->ess / pmod->tss);
	if (pmod->dfd > 0) {
	    double den = pmod->tss * pmod->dfd;

	    pmod->adjrsq = 1.0 - (pmod->ess * (pmod->nobs - 1) / den);
	}
    }

    pmod->list[1] = dv;

    gretl_model_set_double(pmod, "lmax", lmax);

    pmod->ci = LOGISTIC;

    gretl_aic_etc(pmod);

    return 0;
}

MODEL logistic_estimate (int *list, double ***pZ, DATAINFO *pdinfo,
			 const char *param) 
{
    double lmax;
    int dv = list[1];
    MODEL lmod;

    gretl_model_init(&lmod, pdinfo); 

    lmax = get_lmax((*pZ)[dv], pdinfo, param);

    if (na(lmax)) {
	lmod.errcode = E_DATA;
	return lmod;
    } else if (lmax == 0.0) {
	lmod.errcode = E_CANCEL;
	return lmod;
    }

    if (make_logistic_depvar(pZ, pdinfo, dv, lmax)) {
	lmod.errcode = E_ALLOC;	
	return lmod;
    }

    list[1] = pdinfo->v - 1;

    lmod = lsq(list, pZ, pdinfo, OLS, OPT_D, 0.0);
    if (!lmod.errcode) {
	rewrite_logistic_stats((const double **) *pZ, pdinfo, &lmod,
			       dv, lmax);
    }

    dataset_drop_vars(1, pZ, pdinfo);
    list[1] = dv;
    
    return lmod;
}



    

    

    
    
