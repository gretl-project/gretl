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
#include "version.h"
#include "gretl_matrix.h"

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
    *finfo->flag = 0;
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
    GtkWidget *dialog, *tmp, *button, *vbox, *hbox;
    GtkWidget *internal_vbox;
    unsigned char flag = SAVE_LEVERAGE | SAVE_INFLUENCE | SAVE_DFFITS;

    finfo = malloc(sizeof *finfo);
    if (finfo == NULL) return 0;

    dialog = gtk_dialog_new();

    finfo->dialog = dialog;
    finfo->flag = &flag;

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));

    gtk_window_set_title(GTK_WINDOW(dialog), _("gretl: save data"));
    gtk_window_set_resizable(GTK_WINDOW(dialog), FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 10);
    gtk_container_set_border_width(GTK_CONTAINER(hbox), 5);
    gtk_box_set_spacing(GTK_BOX(vbox), 5);

    gtk_window_set_position(GTK_WINDOW(dialog), GTK_WIN_POS_MOUSE);

    g_signal_connect(G_OBJECT(dialog), "destroy",
		     G_CALLBACK(destroy_save_dialog), finfo);

    internal_vbox = gtk_vbox_new(FALSE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new(_("Variables to save:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
    gtk_box_pack_start(GTK_BOX(internal_vbox), hbox, TRUE, TRUE, 5);

    /* Leverage */
    button = gtk_check_button_new_with_label(_("leverage"));
    gtk_box_pack_start(GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(update_save_flag), finfo);
    finfo->levcheck = button;

    /* Influence */
    button = gtk_check_button_new_with_label(_("influence"));
    gtk_box_pack_start(GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(update_save_flag), finfo);
    finfo->infcheck = button;

    /* DFFITS */
    button = gtk_check_button_new_with_label(_("DFFITS"));
    gtk_box_pack_start(GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON (button), TRUE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(update_save_flag), finfo);
    finfo->dffcheck = button;

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), internal_vbox, TRUE, TRUE, 5);

    gtk_widget_show(internal_vbox);
    gtk_box_pack_start(GTK_BOX(vbox), hbox, TRUE, TRUE, 5);

    hbox = gtk_dialog_get_action_area(GTK_DIALOG(dialog));
    gtk_button_box_set_layout(GTK_BUTTON_BOX(hbox), GTK_BUTTONBOX_END);
    gtk_box_set_spacing(GTK_BOX(hbox), 10);

    /* Cancel button */
    tmp = gtk_button_new_from_stock(GTK_STOCK_CANCEL);
    gtk_container_add(GTK_CONTAINER(hbox), tmp);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(cancel_set_flag), finfo);

    /* "OK" button */
    tmp = gtk_button_new_from_stock(GTK_STOCK_OK);
    gtk_container_add(GTK_CONTAINER(hbox), tmp);
    g_signal_connect(G_OBJECT(tmp), "clicked",
		     G_CALLBACK(save_dialog_finalize), finfo);
    gtk_widget_set_can_default(tmp, TRUE);
    gtk_widget_grab_default(tmp);

    gtk_widget_show_all(dialog);

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
			  DATASET *dset)
{
    FILE *fp;
    const double *obs = NULL;
    int t, err = 0;

    fp = open_plot_input_file(PLOT_LEVERAGE, 0, &err);
    if (err) {
	return err;
    }

    if (dataset_is_time_series(dset)) {
	obs = gretl_plotx(dset, OPT_NONE);
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
    gnuplot_missval_string(fp);
    fputs("set yrange [0:1]\n", fp);
    fprintf(fp, "set title '%s'\n", _("leverage"));
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
    gnuplot_missval_string(fp);
    fputs("set yrange [*:*]\n", fp);
    fprintf(fp, "set title '%s'\n", _("influence"));
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

    fputs("unset multiplot\n", fp);

    gretl_pop_c_numeric_locale();

    return finalize_plot_input_file(fp);
}

static void leverage_print (const MODEL *pmod,
			    gretl_matrix *S,
			    double Xvalcrit,
			    DATASET *dset,
			    PRN *prn)
{
    double lp = 2.0 * pmod->ncoeff / pmod->nobs;
    int obslen = max_obs_marker_length(dset);
    int t, j, gotlp = 0;

    if (obslen < 8) {
	obslen = 8;
    }

    bufspace(obslen, prn);
    pprintf(prn, "%*s", UTF_WIDTH(_("residual"), 16), _("residual"));
    pprintf(prn, "%*s", UTF_WIDTH(_("leverage"), 16), _("leverage"));
    pprintf(prn, "%*s", UTF_WIDTH(_("influence"), 16), _("influence"));
    pprintf(prn, "%*s", UTF_WIDTH(_("DFFITS"), 14), _("DFFITS"));
    pputc(prn, '\n');
    bufspace(obslen, prn);
    pputs(prn, "            u          0<=h<=1         u*h/(1-h)\n\n");

    for (t=pmod->t1, j=0; t<=pmod->t2; t++, j++) {
	double h, st, d, f;
	gchar *fstr = NULL;

	if (na(pmod->uhat[t])) {
	    print_obs_marker(t, dset, obslen, prn);
	    pputc(prn, '\n');
	    continue;
	}

	h = gretl_matrix_get(S, j, 0);
	if (h > lp) {
	    gotlp = 1;
	}

	f = gretl_matrix_get(S, j, 1);
	if (na(f)) {
	    fstr = g_strdup(_("undefined"));
	    gretl_utf8_truncate(fstr, 15);
	} else {
	    fstr = g_strdup_printf("%15.5g", f);
	}

	print_obs_marker(t, dset, obslen, prn);

	st = gretl_matrix_get(S, j, 2);
	d = st * sqrt(h / (1.0 - h));
	pprintf(prn, "%14.5g %14.3f%s %s %14.3f\n", pmod->uhat[t], h,
		(h > lp)? "*" : " ", fstr, d);
	g_free(fstr);
    }

    if (gotlp) {
	pprintf(prn, "\n%s\n", _("('*' indicates a leverage point)"));
    } else {
	pprintf(prn, "\n%s\n", _("No leverage points were found"));
    }

    pprintf(prn, "\n%s = %g\n\n", _("Cross-validation criterion"), Xvalcrit);
}

gretl_matrix *model_leverage (const MODEL *pmod, DATASET *dset,
			      gretlopt opt, PRN *prn, int *err)
{
    gretl_matrix *Q, *S = NULL;
    char **cnames = NULL;
    double Xvalcrit, s21;
    int i, j, s, t, vi, df;
    /* allow for missing obs in model range */
    int modn = pmod->t2 - pmod->t1 + 1;
    int k = pmod->list[0] - 1;

    Q = gretl_matrix_alloc(pmod->nobs, k);
    if (Q == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* copy regressors into Q, skipping missing obs */
    j = 0;
    for (i=2; i<=pmod->list[0]; i++) {
	vi = pmod->list[i];
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (!na(pmod->uhat[t])) {
		Q->val[j++] = dset->Z[vi][t];
	    }
	}
    }

    *err = gretl_matrix_QR_decomp(Q, NULL);
    if (*err) {
	goto qr_cleanup;
    }

    /* allocate the results matrix */
    S = gretl_matrix_alloc(modn, 3);
    if (S == NULL) {
	*err = E_ALLOC;
	goto qr_cleanup;
    }

    cnames = malloc(3 * sizeof *cnames);
    /* maybe mark cnames as translatable? */
    cnames[0] = gretl_strdup("leverage");
    cnames[1] = gretl_strdup("influence");
    cnames[2] = gretl_strdup("studres");
    gretl_matrix_set_colnames(S, cnames);
    gretl_matrix_set_t1(S, pmod->t1);

    /* revised df and first component of variance for the
       studentized residuals
    */
    df = pmod->dfd - 1;
    s21 = pmod->ess / df;

    /* initialize cross-validation criterion */
    Xvalcrit = 0.0;

    /* do the "h" calculations, etc. */
    s = 0;
    for (t=pmod->t1, i=0; t<=pmod->t2; t++, i++) {
	double s22, f = NADBL, d = NADBL;
	double q, h, et = pmod->uhat[t];

	if (na(et)) {
	    h = NADBL;
	} else {
	    h = 0.0;
	    for (j=0; j<k; j++) {
		q = gretl_matrix_get(Q, s, j);
		h += q * q;
	    }
	    if (h < 1.0) {
		f = et / (1 - h);
		Xvalcrit += f * f;
		f -= et;
		/* studentized residual (note: agrees with R, Stata) */
		s22 = et * et / (df * (1 - h));
		d = et / sqrt((s21 - s22) * (1 - h));
	    }
	    s++;
	}
	gretl_matrix_set(S, i, 0, h);
	gretl_matrix_set(S, i, 1, f);
	gretl_matrix_set(S, i, 2, d);
    }

    record_test_result(Xvalcrit, NADBL);

    /* print the results, unless in quiet mode */
    if (!(opt & OPT_Q)) {
	leverage_print(pmod, S, Xvalcrit, dset, prn);
	if (gnuplot_graph_wanted(PLOT_LEVERAGE, opt, err)) {
	    *err = leverage_plot(pmod, S, dset);
	}
    }

 qr_cleanup:

    gretl_matrix_free(Q);

    return S;
}
