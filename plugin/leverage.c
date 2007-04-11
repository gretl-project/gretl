/*
 *  Copyright (c) 2003 by Allin Cottrell
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
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 */

#include "libgretl.h"
#include "gretl_matrix.h"

#include "f2c.h"
#include "clapack_double.h"

#include <gtk/gtk.h>

struct flag_info {
    GtkWidget *dialog;
    GtkWidget *levcheck;
    GtkWidget *infcheck;
    GtkWidget *dffcheck;
    unsigned char *flag;
};

static gboolean destroy_save_dialog (GtkWidget *w, struct flag_info *finfo)
{
    free(finfo);
    gtk_main_quit();
    return FALSE;
}

static gboolean update_save_flag (GtkWidget *w, struct flag_info *finfo)
{
    int checked = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(w));

    if (w == finfo->levcheck) {
	if (checked) {
	    *finfo->flag |= SAVE_LEVERAGE;
	} else {
	    *finfo->flag &= ~SAVE_LEVERAGE;
	}
    } else if (w == finfo->infcheck) {
	if (checked) {
	    *finfo->flag |= SAVE_INFLUENCE;
	} else {
	    *finfo->flag &= ~SAVE_INFLUENCE;
	}
    } else if (w == finfo->dffcheck) {
	if (checked) {
	    *finfo->flag |= SAVE_DFFITS;
	} else {
	    *finfo->flag &= ~SAVE_DFFITS;
	}
    }    

    return FALSE;
}

static gboolean cancel_set_flag (GtkWidget *w, struct flag_info *finfo)
{
    *(finfo->flag) = 0;
    gtk_widget_destroy(finfo->dialog);
    return FALSE;
}

static gboolean save_dialog_finalize (GtkWidget *w, struct flag_info *finfo)
{
    gtk_widget_destroy(finfo->dialog);
    return FALSE;
}

unsigned char leverage_data_dialog (void)
{
    struct flag_info *finfo;
    GtkWidget *dialog, *tmp, *button, *hbox;
    GtkWidget *internal_vbox;
    unsigned char flag = SAVE_LEVERAGE | SAVE_INFLUENCE | SAVE_DFFITS;

    finfo = malloc(sizeof *finfo);
    if (finfo == NULL) return 0;

    dialog = gtk_dialog_new();

    finfo->dialog = dialog;
    finfo->flag = &flag;
    
    gtk_window_set_title(GTK_WINDOW (dialog), _("gretl: save data")); 
    gtk_window_set_resizable(GTK_WINDOW(dialog), FALSE);
    gtk_container_set_border_width(GTK_CONTAINER 
				   (GTK_DIALOG(dialog)->vbox), 10);
    gtk_container_set_border_width(GTK_CONTAINER 
				   (GTK_DIALOG(dialog)->action_area), 5);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG (dialog)->vbox), 5);

    gtk_window_set_position(GTK_WINDOW(dialog), GTK_WIN_POS_MOUSE);

    g_signal_connect(G_OBJECT(dialog), "destroy", 
		     G_CALLBACK(destroy_save_dialog), finfo);

    internal_vbox = gtk_vbox_new(FALSE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Variables to save:"));
    gtk_box_pack_start (GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
    gtk_widget_show(tmp);
    gtk_box_pack_start(GTK_BOX(internal_vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show(hbox); 

    /* Leverage */
    button = gtk_check_button_new_with_label(_("leverage"));
    gtk_box_pack_start(GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(update_save_flag), finfo);
    gtk_widget_show (button);
    finfo->levcheck = button;

    /* Influence */
    button = gtk_check_button_new_with_label(_("influence"));
    gtk_box_pack_start (GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(update_save_flag), finfo);
    gtk_widget_show (button);
    finfo->infcheck = button;

    /* DFFITS */
    button = gtk_check_button_new_with_label(_("DFFITS"));
    gtk_box_pack_start (GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(update_save_flag), finfo);
    gtk_widget_show (button);
    finfo->dffcheck = button;

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), internal_vbox, TRUE, TRUE, 5);
    gtk_widget_show(hbox);

    gtk_widget_show(internal_vbox);
    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show(hbox);

    hbox = GTK_DIALOG(dialog)->action_area;
    gtk_button_box_set_layout(GTK_BUTTON_BOX(hbox), GTK_BUTTONBOX_END);
    gtk_button_box_set_spacing(GTK_BUTTON_BOX(hbox), 10);

    /* Cancel button */
    tmp = gtk_button_new_from_stock(GTK_STOCK_CANCEL);
    gtk_container_add(GTK_CONTAINER(hbox), tmp);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(cancel_set_flag), finfo);
    gtk_widget_show(tmp);    

    /* "OK" button */
    tmp = gtk_button_new_from_stock(GTK_STOCK_OK);
    gtk_container_add(GTK_CONTAINER(hbox), tmp);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(save_dialog_finalize), finfo);
    GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
    gtk_widget_grab_default(tmp);
    gtk_widget_show(tmp);

    gtk_widget_show(dialog);

    gtk_main();

    return flag;
}

static void 
leverage_x_range (int t1, int t2, const double *x, FILE *fp)
{
    double xrange, xmin;
    double xmin0 = x[t1];
    double xmax = x[t2];

    xrange = xmax - xmin0;
    xmin = xmin0 - xrange * .025;
    xmax += xrange * .025;

    if (xmin0 >= 0.0 && xmin < 0.0) {
	xmin = 0.0;
    }

    fprintf(fp, "set xrange [%.7g:%.7g]\n", xmin, xmax);
}

static int leverage_plot (const MODEL *pmod, gretl_matrix *S,
			  double ***pZ, DATAINFO *pdinfo)
{
    FILE *fp = NULL;
    const double *obs = NULL;
    int t;

    if (gnuplot_init(PLOT_LEVERAGE, &fp)) {
	return E_FOPEN;
    }

    if (dataset_is_time_series(pdinfo)) { 
	obs = gretl_plotx(pdinfo);
	if (obs == NULL) {
	    if (fp != NULL) {
		fclose(fp);
	    }
	    return 1;
	}
    }

    gretl_push_c_numeric_locale();

    fputs("set size 1.0,1.0\nset multiplot\nset size 1.0,0.48\n", fp);
    fputs("set xzeroaxis\n", fp);
    fputs("set nokey\n", fp); 

    if (obs != NULL) {
	leverage_x_range(pmod->t1, pmod->t2, obs, fp);
    } else {
	fprintf(fp, "set xrange [%g:%g]\n", 
		pmod->t1 + 0.5, pmod->t2 + 1.5);
    } 

    /* upper plot: leverage factor */
    fputs("set origin 0.0,0.50\n", fp);
    fputs("set missing '?'\n", fp);
    fputs("set yrange [0:1]\n", fp);
    fprintf(fp, "set title '%s'\n", I_("leverage"));
    fputs("plot \\\n'-' using 1:2 w impulses\n", fp);

    for (t=pmod->t1; t<=pmod->t2; t++) {
	double h = gretl_matrix_get(S, t - pmod->t1, 0);

	if (na(h)) {
	    if (obs != NULL) {
		fprintf(fp, "%g ?\n", obs[t]);
	    } else { 
		fprintf(fp, "%d ?\n", t + 1);
	    }
	} else {
	    if (obs != NULL) {
		fprintf(fp, "%g %g\n", obs[t], h);
	    } else { 
		fprintf(fp, "%d %g\n", t + 1, h);
	    }
	} 
    }
    fputs("e\n", fp);

    /* lower plot: influence factor */
    fputs("set origin 0.0,0.0\n", fp);
    fputs("set missing '?'\n", fp);
    fputs("set yrange [*:*]\n", fp);
    fprintf(fp, "set title '%s'\n", I_("influence")); 
    fputs("plot \\\n'-' using 1:2 w impulses\n", fp);

    for (t=pmod->t1; t<=pmod->t2; t++) {
	double f = gretl_matrix_get(S, t - pmod->t1, 1);

	if (na(f)) {
	    if (obs != NULL) {
		fprintf(fp, "%g ?\n", obs[t]);
	    } else {
		fprintf(fp, "%d ?\n", t + 1);
	    }
	} else {
	    if (obs != NULL) {
		fprintf(fp, "%g %g\n", obs[t], f);
	    } else {
		fprintf(fp, "%d %g\n", t + 1, f);
	    }
	}
    }
    fputs("e\n", fp);

    fputs("set nomultiplot\n", fp);

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return 0;
}

static int studentized_residuals (const MODEL *pmod, double ***pZ, 
				  DATAINFO *pdinfo, gretl_matrix *S)
{
    double *dum = NULL;
    int *slist = NULL;
    MODEL smod;  
    int orig_v = pdinfo->v;
    int err = 0;
    int i, t, ts, k;

    /* create a full-length dummy variable */
    dum = malloc(pdinfo->n * sizeof *dum);
    if (dum == NULL) {
	return E_ALLOC;
    }

    /* allocate regression list */
    slist = malloc((pmod->list[0] + 2) * sizeof *slist);
    if (slist == NULL) {
	free(dum);
	return E_ALLOC;
    }

    if (dataset_add_allocated_series(dum, pZ, pdinfo)) {
	free(dum);
	free(slist);
	return E_ALLOC;	
    }

    /* zero out the dummy */
    for (t=0; t<pdinfo->n; t++) {
	dum[t] = 0.0;
    }

    slist[0] = pmod->list[0] + 1;
    for (i=1; i<=pmod->list[0]; i++) {
	slist[i] = pmod->list[i];
    }
    slist[slist[0]] = pdinfo->v - 1; /* last var added */  
    k = slist[0] - 2;

    for (t=pmod->t1; t<=pmod->t2 && !err; t++) {
	ts = t - pmod->t1;

	if (model_missing(pmod, t)) {
	    gretl_matrix_set(S, ts, 2, NADBL);
	    dum[t-1] = 0.0;
	    continue;
	}	
	dum[t] = 1.0;
	if (t > pmod->t1) {
	    dum[t-1] = 0.0;
	}
	smod = lsq(slist, pZ, pdinfo, OLS, OPT_A);
	if (smod.errcode) {
	    err = smod.errcode;
	} else {
	    gretl_matrix_set(S, ts, 2, smod.coeff[k] / smod.sderr[k]);
	}
	clear_model(&smod);
    }

    if (err) {
	int modn = pmod->t2 - pmod->t1 + 1;

	for (t=0; t<modn; t++) {
	    gretl_matrix_set(S, t, 2, NADBL);
	}
    }

    free(slist);

    dataset_drop_last_variables(pdinfo->v - orig_v, pZ, pdinfo);

    return err;
}

/* In fortran arrays, column entries are contiguous.
   Columns of data matrix X hold variables, rows hold observations.
   So in a fortran array, entries for a given variable are
   contiguous.
*/

gretl_matrix *model_leverage (const MODEL *pmod, double ***pZ, 
			      DATAINFO *pdinfo, gretlopt opt,
			      PRN *prn, int *err)
{
    integer info, lwork;
    integer m, n, lda;
    gretl_matrix *Q, *S = NULL;
    doublereal *tau, *work;
    double lp;
    int i, j, t, tq, tmod;
    int serr = 0, gotlp = 0;
    /* allow for missing obs in model range */
    int modn = pmod->t2 - pmod->t1 + 1;

    m = pmod->nobs;              /* # of rows = # of observations */
    lda = m;                     /* leading dimension of Q */
    n = pmod->list[0] - 1;       /* # of cols = # of variables */

    Q = gretl_matrix_alloc(m, n);

    /* dim of tau is min (m, n) */
    tau = malloc(n * sizeof *tau);
    work = malloc(sizeof *work);

    if (Q == NULL || tau == NULL || work == NULL) {
	*err = E_ALLOC;
	goto qr_cleanup;
    }

    /* copy independent var values into Q */
    j = 0;
    for (i=2; i<=pmod->list[0]; i++) {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (!model_missing(pmod, t)) {
		Q->val[j++] = (*pZ)[pmod->list[i]][t];
	    }
	}
    }

    /* do a workspace size query */
    lwork = -1;
    info = 0;
    dgeqrf_(&m, &n, Q->val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
	fprintf(stderr, "dgeqrf: info = %d\n", (int) info);
	*err = 1;
	goto qr_cleanup;
    }

    /* set up optimally sized work array */
    lwork = (integer) work[0];
    work = realloc(work, (size_t) lwork * sizeof *work);
    if (work == NULL) {
	*err = E_ALLOC;
	goto qr_cleanup;
    }

    /* run actual QR factorization */
    dgeqrf_(&m, &n, Q->val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
	fprintf(stderr, "dgeqrf: info = %d\n", (int) info);
	*err = 1;
	goto qr_cleanup;
    }

    /* obtain the real "Q" matrix */
    dorgqr_(&m, &n, &n, Q->val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
	*err = 1;
	goto qr_cleanup;
    }

    free(tau);
    tau = NULL;
    free(work);
    work = NULL;

    S = gretl_matrix_alloc(modn, 3);
    if (S == NULL) {
	*err = E_ALLOC;
	goto qr_cleanup;
    }	

    pputs(prn, "        ");
    pprintf(prn, "%*s", UTF_WIDTH(_("residual"), 16), _("residual"));
    pprintf(prn, "%*s", UTF_WIDTH(_("leverage"), 16), _("leverage"));
    pprintf(prn, "%*s", UTF_WIDTH(_("influence"), 16), _("influence"));
    pprintf(prn, "%*s", UTF_WIDTH(_("DFFITS"), 14), _("DFFITS"));
    pputs(prn, "\n        ");
    pputs(prn, "            u          0<=h<=1         u*h/(1-h)\n\n");

    /* do the "h" calculations */
    tq = 0;
    for (t=0; t<modn; t++) {
	double q, h = 0.0;

	if (model_missing(pmod, t + pmod->t1)) {
	    gretl_matrix_set(S, t, 0, NADBL);
	} else {
	    for (i=0; i<n; i++) {
		q = gretl_matrix_get(Q, tq, i);
		h += q * q;
	    }
	    tq++;
	    gretl_matrix_set(S, t, 0, h);
	}
    }

    /* put studentized resids into S[2] */
    serr = studentized_residuals(pmod, pZ, pdinfo, S);

    lp = 2.0 * n / m;

    obs_marker_init(pdinfo);

    /* print the results */
    for (t=0; t<modn; t++) {
	double h, s, d, f = NADBL;
	char fstr[32];

	tmod = t + pmod->t1;

	if (model_missing(pmod, tmod)) {
	    print_obs_marker(tmod, pdinfo, prn);
	    gretl_matrix_set(S, t, 1, f);
	    pputc(prn, '\n');
	    continue;
	}

	h = gretl_matrix_get(S, t, 0);
	if (h > lp) {
	    gotlp = 1;
	}

	if (h < 1.0) {
	    f = pmod->uhat[tmod] * h / (1.0 - h);
	    sprintf(fstr, "%15.5g", f);
	} else {
	    f = NADBL;
	    sprintf(fstr, "%15s", _("undefined"));
	}

	print_obs_marker(tmod, pdinfo, prn);

	if (!serr) {
	    s = gretl_matrix_get(S, t, 2);
	    d = s * sqrt(h / (1.0 - h));
	    pprintf(prn, "%14.5g %14.3f%s %s %14.3f\n", pmod->uhat[tmod], h, 
		    (h > lp)? "*" : " ", fstr, d);

	} else {
	    pprintf(prn, "%14.5g %14.3f%s %s\n", pmod->uhat[tmod], h, 
		    (h > lp)? "*" : " ", fstr);
	}
	gretl_matrix_set(S, t, 1, f);
    }

    if (gotlp) {
	pprintf(prn, "\n%s\n\n", _("('*' indicates a leverage point)"));
    } else {
	pprintf(prn, "\n%s\n\n", _("No leverage points were found"));
    }

    if (opt & OPT_P) {
	leverage_plot(pmod, S, pZ, pdinfo);
    }

 qr_cleanup:

    if (Q != NULL) {
	gretl_matrix_free(Q);
    }
    if (tau != NULL) {
	free(tau); 
    }
    if (work != NULL) {
	free(work);
    }

    if (S != NULL) {
	gretl_matrix_set_t1(S, pmod->t1);
    }

    return S;    
}

