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

#include <gtk/gtk.h>

#undef PCA_DEBUG

struct flag_info {
    GtkWidget *dialog;
    gint *flag;
};

enum pca_flags {
    PCA_SAVE_NONE,
    PCA_SAVE_MAIN,
    PCA_SAVE_ALL
};

static gboolean destroy_pca_dialog (GtkWidget *w, struct flag_info *finfo)
{
    free(finfo);
    gtk_main_quit();
    return FALSE;
}

static gboolean set_pca_flag (GtkWidget *w, struct flag_info *finfo)
{
    gint opt = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(w), "opt"));

    *(finfo->flag) = opt;
    return FALSE;
}

static gboolean cancel_set_flag (GtkWidget *w, struct flag_info *finfo)
{
    *(finfo->flag) = PCA_SAVE_NONE;
    gtk_widget_destroy(finfo->dialog);
    return FALSE;
}

static gboolean pca_dialog_finalize (GtkWidget *w, struct flag_info *finfo)
{
    gtk_widget_destroy(finfo->dialog);
    return FALSE;
}

static gretlopt pca_flag_dialog (void)
{
    struct flag_info *finfo;
    GtkWidget *dialog, *tmp, *button, *hbox;
    GtkWidget *internal_vbox;
    GSList *group;
    gint flag = PCA_SAVE_MAIN;

    finfo = malloc(sizeof *finfo);
    if (finfo == NULL) return 0;

    dialog = gtk_dialog_new();

    finfo->dialog = dialog;
    finfo->flag = &flag;
    
    gtk_window_set_title(GTK_WINDOW(dialog), _("gretl: save data")); 
    gtk_window_set_resizable(GTK_WINDOW(dialog), FALSE);
    gtk_container_set_border_width(GTK_CONTAINER 
				   (GTK_DIALOG(dialog)->vbox), 10);
    gtk_box_set_spacing(GTK_BOX(GTK_DIALOG(dialog)->vbox), 5);
    gtk_window_set_position(GTK_WINDOW(dialog), GTK_WIN_POS_MOUSE);

    g_signal_connect(G_OBJECT(dialog), "destroy", 
		     G_CALLBACK(destroy_pca_dialog), finfo);

    internal_vbox = gtk_vbox_new(FALSE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new (_("Variables to save:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
    gtk_widget_show(tmp);
    gtk_box_pack_start(GTK_BOX(internal_vbox), hbox, TRUE, TRUE, 5);
    gtk_widget_show(hbox); 

    /* Only those with eigenvalues > 1.0 */
    button = gtk_radio_button_new_with_label(NULL, 
					     _("Components with eigenvalues > 1.0"));
    gtk_box_pack_start(GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_pca_flag), finfo);
    g_object_set_data(G_OBJECT(button), "opt", GINT_TO_POINTER(PCA_SAVE_MAIN)); 
    gtk_widget_show (button);   

    /* All components */
    group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("All components"));
    gtk_box_pack_start (GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), FALSE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_pca_flag), finfo);
    g_object_set_data(G_OBJECT(button), "opt", GINT_TO_POINTER(PCA_SAVE_ALL)); 
    gtk_widget_show (button);

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
		     G_CALLBACK(pca_dialog_finalize), finfo);
    GTK_WIDGET_SET_FLAGS(tmp, GTK_CAN_DEFAULT);
    gtk_widget_grab_default(tmp);
    gtk_widget_show(tmp);

    gtk_widget_show(dialog);

    gtk_main();

    if (flag == PCA_SAVE_MAIN) return OPT_O;
    if (flag == PCA_SAVE_ALL) return OPT_A;

    return 0L;
}

static void pca_print (VMatrix *vmat, gretl_matrix *m,
		       gretl_matrix *evals, PRN *prn)
{
    double x, y;
    int n = vmat->dim;
    int i, j, cols;

    pprintf(prn, "%s\n\n", _("Principal Components Analysis"));
    pprintf(prn, "%s\n\n", _("Eigenanalysis of the Correlation Matrix"));

    pputs(prn, _("Component  Eigenvalue  Proportion   Cumulative\n"));

    x = 0.0;
    y = 0.0;

    for (i=n-1; i>=0; i--) {
	y += evals->val[i] / n;
	pprintf(prn, "%5d%13.4f%13.4f%13.4f\n", n - i,
		evals->val[i], evals->val[i] / n, y);
	x += evals->val[i];
    }
    pputc(prn, '\n');

#ifdef PCA_DEBUG
    fprintf(stderr, "check: sum of evals = %g\n", x);
#endif

    pprintf(prn, "%s\n\n", _("Eigenvectors (component loadings)"));

    cols = n;
    while (cols > 0) {
	int colsdone = 0;

	pprintf(prn, "%-16s", _("Variable"));
	for (i=n-cols; i<n-cols+7 && i<n; i++) {
	    char pcname[8];

	    sprintf(pcname, "PC%d", i + 1);
	    pprintf(prn, "%9s", pcname);
	    colsdone++;
	}
	pputc(prn, '\n');
	for (i=0; i<n; i++) {
	    pprintf(prn, "%-16s", vmat->names[i]);
	    for (j=cols-1; j>cols-8 && j>=0; j--) {
		pprintf(prn, "%9.3f", gretl_matrix_get(m, i, j));
	    }
	    pputc(prn, '\n');
	}
	cols -= colsdone;
	pputc(prn, '\n');
    }
}

static int standardize (double *y, const double *x, int n)
{
    double xbar, sd;
    int i, err;

    err = gretl_moments(0, n-1, x, &xbar, &sd, NULL, NULL, 1);
    if (err) {
	return err;
    }

    for (i=0; i<n; i++) {
	if (na(x[i])) {
	    y[i] = NADBL;
	} else {
	    y[i] = (x[i] - xbar) / sd;
	}
    }

    return 0;
}

int pca_from_corrmat (VMatrix *corrmat, double ***pZ,
		      DATAINFO *pdinfo, gretlopt *popt,
		      PRN *prn)
{
    gretl_matrix *C;
    gretl_matrix *evals = NULL;
    int k = corrmat->dim;
    int i, j, t, vi, idx;
    double x;
    gretlopt opt = OPT_NONE;
    int err = 0;

    if (popt != NULL) {
	opt = *popt;
    }

    if (opt & OPT_D) { 
	opt = pca_flag_dialog();
	if (!opt) {
	    /* canceled */
	    *popt = OPT_NONE;
	    return 0; 
	}
    }    

    C = gretl_matrix_alloc(k, k);
    if (C == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<k; i++) {
	for (j=0; j<k; j++) {
	    idx = ijton(i, j, k);
	    x = corrmat->vec[idx];
	    gretl_matrix_set(C, i, j, x);
	}
    }

    evals = gretl_symmetric_matrix_eigenvals(C, 1, &err);
    if (err) {
	gretl_matrix_free(C);
	return err;
    }

    if (prn != NULL) {
	pca_print(corrmat, C, evals, prn);
    }

    if (opt) {
	/* add PCs to the dataset */
	double **sZ = NULL;
	int *plist = NULL;
	int m, v = pdinfo->v;

	if (opt & OPT_A) {
	    m = k;
	} else {
	    m = 0;
	    for (i=0; i<k; i++) {
		if (evals->val[i] > 1.0) {
		    m++;
		}
	    }
	}

	plist = gretl_list_new(m);
	if (plist == NULL) {
	    err = E_ALLOC;
	}

	if (!err) {
	    /* build list of PCs */
	    j = 1;
	    for (i=k-1; i>=0; i--) {
		if ((opt & OPT_A) || evals->val[i] > 1.0) {
		    plist[j++] = i;
		}
	    }
	    err = dataset_add_series(m, pZ, pdinfo);
	}

	if (!err) {
	    /* construct standardized versions of variables */
	    sZ = doubles_array_new(k, pdinfo->n); 
	    if (sZ == NULL) {
		err = E_ALLOC;
	    } else {
		for (i=0; i<k && !err; i++) {
		    vi = corrmat->list[i+1];
		    err = standardize(sZ[i], (const double *) (*pZ)[vi], 
				      pdinfo->n);
		}
	    }
	}

	if (!err) {
	    for (i=1; i<=plist[0]; i++) {
		int pi = plist[i];
		double load;

		vi = v + i - 1;
		sprintf(pdinfo->varname[vi], "PC%d", i);
		make_varname_unique(pdinfo->varname[vi], vi, pdinfo);
		sprintf(VARLABEL(pdinfo, vi), "Component with "
			"eigenvalue = %.4f", evals->val[pi]);

		for (t=0; t<pdinfo->n; t++) {
		    (*pZ)[vi][t] = 0.0;
		    for (j=0; j<k; j++) {
			x = sZ[j][t];
			if (na(x)) {
			    (*pZ)[vi][t] = NADBL;
			    break;
			} else {
			    load = gretl_matrix_get(C, j, pi);
			    (*pZ)[vi][t] += load * x;
			}
		    }
		}
	    }
	}

	free(plist);
	doubles_array_free(sZ, k);

    } /* end opt (save PCs) conditional */

    gretl_matrix_free(evals);
    gretl_matrix_free(C);

    if (popt != NULL) {
	*popt = opt;
    }

    return 0;
}
