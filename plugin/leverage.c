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
#include "gretl_matrix_private.h"

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
	if (checked) *finfo->flag |= SAVE_LEVERAGE;
	else *finfo->flag &= ~SAVE_LEVERAGE;
    }
    else if (w == finfo->infcheck) {
	if (checked) *finfo->flag |= SAVE_INFLUENCE;
	else *finfo->flag &= ~SAVE_INFLUENCE;
    }
    else if (w == finfo->dffcheck) {
	if (checked) *finfo->flag |= SAVE_DFFITS;
	else *finfo->flag &= ~SAVE_DFFITS;
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
    GtkWidget *dialog, *tempwid, *button, *hbox;
    GtkWidget *internal_vbox;
    unsigned char flag = SAVE_LEVERAGE | SAVE_INFLUENCE | SAVE_DFFITS;

    finfo = malloc(sizeof *finfo);
    if (finfo == NULL) return 0;

    dialog = gtk_dialog_new();

    finfo->dialog = dialog;
    finfo->flag = &flag;
    
    gtk_window_set_title (GTK_WINDOW (dialog), _("gretl: save data")); 
#if GTK_MAJOR_VERSION >= 2
    gtk_window_set_resizable (GTK_WINDOW (dialog), FALSE);
#endif
    gtk_container_set_border_width (GTK_CONTAINER 
				    (GTK_DIALOG (dialog)->vbox), 10);
    gtk_container_set_border_width (GTK_CONTAINER 
				    (GTK_DIALOG (dialog)->action_area), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (dialog)->vbox), 5);

    gtk_window_set_position (GTK_WINDOW (dialog), GTK_WIN_POS_MOUSE);

#if GTK_MAJOR_VERSION >= 2
    g_signal_connect (G_OBJECT(dialog), "destroy", 
		      G_CALLBACK(destroy_save_dialog), finfo);
#else
    gtk_signal_connect (GTK_OBJECT(dialog), "destroy", 
			GTK_SIGNAL_FUNC(destroy_save_dialog), finfo);
#endif

    internal_vbox = gtk_vbox_new (FALSE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    tempwid = gtk_label_new (_("Variables to save:"));
    gtk_box_pack_start (GTK_BOX(hbox), tempwid, TRUE, TRUE, 5);
    gtk_widget_show(tempwid);
    gtk_box_pack_start (GTK_BOX(internal_vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show(hbox); 

    /* Leverage */
    button = gtk_check_button_new_with_label(_("leverage"));
    gtk_box_pack_start (GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
#if GTK_MAJOR_VERSION >= 2
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(update_save_flag), finfo);
#else
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
		       GTK_SIGNAL_FUNC(update_save_flag), finfo);
#endif   
    gtk_widget_show (button);
    finfo->levcheck = button;

    /* Influence */
    button = gtk_check_button_new_with_label(_("influence"));
    gtk_box_pack_start (GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
#if GTK_MAJOR_VERSION >= 2
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(update_save_flag), finfo);
#else
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
		       GTK_SIGNAL_FUNC(update_save_flag), finfo);
#endif
    gtk_widget_show (button);
    finfo->infcheck = button;

    /* DFFITS */
    button = gtk_check_button_new_with_label(_("DFFITS"));
    gtk_box_pack_start (GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
#if GTK_MAJOR_VERSION >= 2
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(update_save_flag), finfo);
#else
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
		       GTK_SIGNAL_FUNC(update_save_flag), finfo);
#endif
    gtk_widget_show (button);
    finfo->dffcheck = button;

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), internal_vbox, TRUE, TRUE, 5);
    gtk_widget_show (hbox);

    gtk_widget_show (internal_vbox);

    gtk_box_pack_start(GTK_BOX(GTK_DIALOG(dialog)->vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show (hbox);

    /* Create the "OK" button */
#if GTK_MAJOR_VERSION >= 2
    tempwid = gtk_button_new_from_stock (GTK_STOCK_OK);
    g_signal_connect(G_OBJECT(tempwid), "clicked",
		     G_CALLBACK(save_dialog_finalize), finfo);
#else
    tempwid = gtk_button_new_with_label(_("OK"));
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
		       GTK_SIGNAL_FUNC(save_dialog_finalize), finfo);
#endif
    gtk_box_pack_start (GTK_BOX(GTK_DIALOG (dialog)->action_area), 
			tempwid, TRUE, TRUE, 0);
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    /* "Cancel" button */
#if GTK_MAJOR_VERSION >= 2
    tempwid = gtk_button_new_from_stock (GTK_STOCK_CANCEL);
    g_signal_connect(G_OBJECT(tempwid), "clicked",
		     G_CALLBACK(cancel_set_flag), finfo);
#else
    tempwid = gtk_button_new_with_label(_("Cancel"));
    gtk_signal_connect(GTK_OBJECT(tempwid), "clicked",
		       GTK_SIGNAL_FUNC(cancel_set_flag), finfo);
#endif    
    gtk_box_pack_start (GTK_BOX(GTK_DIALOG (dialog)->action_area), 
			tempwid, TRUE, TRUE, 0);
    gtk_widget_show (tempwid);

    gtk_widget_show(dialog);

    gtk_main();

    return flag;
}

static int leverage_plot (int n, int tstart, gretl_matrix *S,
			  double ***pZ, DATAINFO *pdinfo, 
			  PATHS *ppaths)
{
    FILE *fp = NULL;
    int t, tmod;
    int timeplot = 0;

    if (gnuplot_init(ppaths, PLOT_LEVERAGE, &fp)) return 1;

    if (dataset_is_time_series(pdinfo) && 
	(pdinfo->pd == 1 || pdinfo->pd == 4 || pdinfo->pd == 12)) {
	char per[8];

	if (pdinfo->pd == 1) strcpy(per, "annual");
	else if (pdinfo->pd == 4) strcpy(per, "qtrs");
	else if (pdinfo->pd == 12) strcpy(per, "months");
	timeplot = plotvar(pZ, pdinfo, per);
	if (timeplot < 0) {
	    if (fp != NULL) fclose(fp);
	    return 1;
	}
    }

    fputs("# leverage/influence plot\n", fp);
    fputs("set size 1.0,1.0\nset multiplot\nset size 1.0,0.48\n", fp);
    fputs("set xzeroaxis\n", fp);
    if (!timeplot) { 
	fprintf(fp, "set xrange [%g:%g]\n", 
		tstart + 0.5, tstart + n + 0.5);
    }
    fputs("set nokey\n", fp); 

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "C");
#endif

    /* upper plot: leverage factor */
    fputs("set origin 0.0,0.50\n", fp);
    fputs("set yrange [0:1]\n", fp);
    fprintf(fp, "set title '%s'\n", I_("leverage"));
    fputs("plot \\\n'-' using 1:2 w impulses\n", fp);
    for (t=0; t<n; t++) {
	double h = gretl_matrix_get(S, t, 0);

	if (timeplot) {
	    tmod = t + tstart;
	    fprintf(fp, "%g %g\n", (*pZ)[timeplot][tmod], h);
	} else { 
	    fprintf(fp, "%d %g\n", t+tstart+1, h);
	}
    }
    fputs("e\n", fp);

    /* lower plot: influence factor */
    fputs("set origin 0.0,0.0\n", fp);
    fputs("set missing '?'\n", fp);
    fputs("set yrange [*:*]\n", fp);
    fprintf(fp, "set title '%s'\n", I_("influence")); 
    fputs("plot \\\n'-' using 1:2 w impulses\n", fp);
    for (t=0; t<n; t++) {
	double f = gretl_matrix_get(S, t, 1);

	tmod = t + tstart;
	if (!na(f)) {
	    if (timeplot) {
		fprintf(fp, "%g %g\n", (*pZ)[timeplot][tmod], f);
	    } else {
		fprintf(fp, "%d %g\n", tmod + 1, f);
	    }
	} else {
	    if (timeplot) {
		fprintf(fp, "%g ?\n", (*pZ)[timeplot][tmod]);
	    } else {
		fprintf(fp, "%d ?\n", tmod + 1);
	    }
	}
    }
    fputs("e\n", fp);
    fputs("set nomultiplot\n", fp);

#ifdef ENABLE_NLS
    setlocale(LC_NUMERIC, "");
#endif

    fclose(fp);

    return 0;
}

static int studentized_residuals (const MODEL *pmod, double ***pZ, 
				  DATAINFO *pdinfo, gretl_matrix *S)
{
    double *dum = NULL, *suhat = NULL;
    int *slist = NULL;
    MODEL smod;  
    int orig_v = pdinfo->v;
    int err = 0;
    int i, t, k;

    dum = malloc(pdinfo->n * sizeof *dum);
    if (dum == NULL) {
	return E_ALLOC;
    }

    suhat = malloc(pdinfo->n * sizeof *suhat);
    if (suhat == NULL) {
	free(dum);
	return E_ALLOC;
    }

    slist = malloc((pmod->list[0] + 2) * sizeof *slist);
    if (slist == NULL) {
	free(dum);
	free(suhat);
	return E_ALLOC;
    }

    if (dataset_add_allocated_var(dum, pZ, pdinfo)) {
	free(dum);
	free(suhat);
	free(slist);
	return E_ALLOC;	
    }

    for (t=0; t<pdinfo->n; t++) {
	dum[t] = 0.0;
    }

    slist[0] = pmod->list[0] + 1;
    for (i=1; i<=pmod->list[0]; i++) {
	slist[i] = pmod->list[i];
    }
    slist[slist[0]] = pdinfo->v - 1; /* last var added */  
    k = slist[0] - 2;

    gretl_model_init(&smod, NULL);

    for (t=pmod->t1; t<=pmod->t2 && !err; t++) {
	dum[t] = 1.0;
	if (t > pmod->t1) dum[t-1] = 0.0;
	smod = lsq(slist, pZ, pdinfo, OLS, OPT_A, 0.0);
	if (smod.errcode) {
	    err = smod.errcode;
	} else {
	    suhat[t] = smod.coeff[k] / smod.sderr[k];
	}
	clear_model(&smod, NULL);
    }

    if (!err) {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    gretl_matrix_set(S, t - pmod->t1, 2, suhat[t]);
	}
    } else {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    gretl_matrix_set(S, t - pmod->t1, 2, NADBL);
	}
    }	

    free(suhat);
    free(slist);
    dataset_drop_vars(pdinfo->v - orig_v, pZ, pdinfo);

    return err;
}

/* In fortran arrays, column entries are contiguous.
   Columns of data matrix X hold variables, rows hold observations.
   So in a fortran array, entries for a given variable are
   contiguous.
*/

gretl_matrix *model_leverage (const MODEL *pmod, double ***pZ, 
			      DATAINFO *pdinfo, PRN *prn,
			      PATHS *ppaths)
{
    integer info, lwork;
    integer m, n, lda;
    gretl_matrix *Q, *S = NULL;
    doublereal *tau, *work;
    double lp;
    int i, j, t;
    int err = 0, serr = 0, gotlp = 0;

    m = pmod->t2 - pmod->t1 + 1; /* # of rows = # of observations */
    lda = m;                     /* leading dimension of Q */
    n = pmod->list[0] - 1;       /* # of cols = # of variables */

    Q = gretl_matrix_alloc(m, n);
    /* dim of tau is min (m, n) */
    tau = malloc(n * sizeof *tau);
    work = malloc(sizeof *work);

    if (Q == NULL || tau == NULL || work == NULL) {
	err = E_ALLOC;
	goto qr_cleanup;
    }

    /* copy independent var values into Q */
    j = 0;
    for (i=2; i<=pmod->list[0]; i++) {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    Q->val[j++] = (*pZ)[pmod->list[i]][t];
	}
    }

    /* do a workspace size query */
    lwork = -1;
    info = 0;
    dgeqrf_(&m, &n, Q->val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
	fprintf(stderr, "dgeqrf: info = %d\n", (int) info);
	err = 1;
	goto qr_cleanup;
    }

    /* set up optimally sized work array */
    lwork = (integer) work[0];
    work = realloc(work, (size_t) lwork * sizeof *work);
    if (work == NULL) {
	err = E_ALLOC;
	goto qr_cleanup;
    }

    /* run actual QR factorization */
    dgeqrf_(&m, &n, Q->val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
	fprintf(stderr, "dgeqrf: info = %d\n", (int) info);
	err = 1;
	goto qr_cleanup;
    }

    /* obtain the real "Q" matrix */
    dorgqr_(&m, &n, &n, Q->val, &lda, tau, work, &lwork, &info);
    if (info != 0) {
	err = 1;
	goto qr_cleanup;
    }

    free(tau);
    tau = NULL;
    free(work);
    work = NULL;

    S = gretl_matrix_alloc(m, 3);
    if (S == NULL) {
	err = 1;
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
    for (t=0; t<m; t++) {
	double q, h = 0.0;

	for (i=0; i<n; i++) {
	    q = gretl_matrix_get(Q, t, i);
	    h += q * q;
	}
	gretl_matrix_set(S, t, 0, h);
    }

    /* put studentized resids into S[2] */
    serr = studentized_residuals(pmod, pZ, pdinfo, S);

    lp = 2.0 * n / m;

    /* print the results */
    for (t=0; t<m; t++) {
	double f, h, s, d;
	int tmod = t + pmod->t1;
	char fstr[24];

	h = gretl_matrix_get(S, t, 0);
	if (h > lp) gotlp = 1;
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

    if (ppaths != NULL) {
	leverage_plot(m, pmod->t1, S, pZ, pdinfo, ppaths);
    }

 qr_cleanup:

    if (Q != NULL) gretl_matrix_free(Q);
    if (tau != NULL) free(tau); 
    if (work != NULL) free(work);

    if (S != NULL) gretl_matrix_set_int(S, pmod->t1);

    return S;    
}

