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

#if (GTK_MAJOR_VERSION == 2 && GTK_MINOR_VERSION < 18)
# include "gtk_compat.h"
#endif

#define PCA_DEBUG 0

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
    GtkWidget *vbox, *internal_vbox;
    GSList *group;
    gint flag = PCA_SAVE_MAIN;

    finfo = malloc(sizeof *finfo);
    if (finfo == NULL) return 0;

    dialog = gtk_dialog_new();

    finfo->dialog = dialog;
    finfo->flag = &flag;

    vbox = gtk_dialog_get_content_area(GTK_DIALOG(dialog));
    
    gtk_window_set_title(GTK_WINDOW(dialog), _("gretl: save data")); 
    gtk_window_set_resizable(GTK_WINDOW(dialog), FALSE);
    gtk_container_set_border_width(GTK_CONTAINER(vbox), 10);
    gtk_box_set_spacing(GTK_BOX(vbox), 5);
    gtk_window_set_position(GTK_WINDOW(dialog), GTK_WIN_POS_MOUSE);

    g_signal_connect(G_OBJECT(dialog), "destroy", 
		     G_CALLBACK(destroy_pca_dialog), finfo);

    internal_vbox = gtk_vbox_new(FALSE, 5);

    hbox = gtk_hbox_new(FALSE, 5);
    tmp = gtk_label_new (_("Variables to save:"));
    gtk_box_pack_start(GTK_BOX(hbox), tmp, TRUE, TRUE, 5);
    gtk_box_pack_start(GTK_BOX(internal_vbox), hbox, TRUE, TRUE, 5);

    /* Only those with eigenvalues > 1.0 */
    button = gtk_radio_button_new_with_label(NULL, 
					     _("Components with eigenvalues > mean"));
    gtk_box_pack_start(GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(button), TRUE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_pca_flag), finfo);
    g_object_set_data(G_OBJECT(button), "opt", GINT_TO_POINTER(PCA_SAVE_MAIN)); 

    /* All components */
    group = gtk_radio_button_get_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, _("All components"));
    gtk_box_pack_start (GTK_BOX(internal_vbox), button, TRUE, TRUE, 0);
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), FALSE);
    g_signal_connect(G_OBJECT(button), "clicked",
		     G_CALLBACK(set_pca_flag), finfo);
    g_object_set_data(G_OBJECT(button), "opt", GINT_TO_POINTER(PCA_SAVE_ALL)); 

    hbox = gtk_hbox_new(FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox), internal_vbox, TRUE, TRUE, 5);

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
		     G_CALLBACK(pca_dialog_finalize), finfo);
    gtk_widget_set_can_default(tmp, TRUE);
    gtk_widget_grab_default(tmp);

    gtk_widget_show_all(dialog);

    gtk_main();

    if (flag == PCA_SAVE_MAIN) {
	return OPT_O;
    } else if (flag == PCA_SAVE_ALL) {
	return OPT_A;
    } else {
	return OPT_NONE;
    }
}

#define PCA_COLS 7

static void pca_print (VMatrix *cmat, gretl_matrix *E,
		       gretl_matrix *C, PRN *prn)
{
    double cum, esum;
    char pcname[8];
    int nl, namelen = 8;
    int n = cmat->dim;
    int done, todo;
    int i, j;

    pprintf(prn, "%s\n\n", _("Principal Components Analysis"));

    if (cmat->ci == CORR) {
	pprintf(prn, "%s\n\n", _("Eigenanalysis of the Correlation Matrix"));
    } else {
	pprintf(prn, "%s\n\n", _("Eigenanalysis of the Covariance Matrix"));
    }

    pputs(prn, _("Component  Eigenvalue  Proportion   Cumulative\n"));

    if (cmat->ci == CORR) {
	esum = n;
    } else {
	esum = 0.0;
	for (i=0; i<n; i++) {
	    esum += E->val[i];
	}
    }

    cum = 0.0;
    for (i=0; i<n; i++) {
	cum += E->val[i] / esum;
	pprintf(prn, "%5d%13.4f%13.4f%13.4f\n", i + 1,
		E->val[i], E->val[i] / esum, cum);
	nl = strlen(cmat->names[i]);
	if (nl > namelen) {
	    namelen = nl;
	}
    }
    pputc(prn, '\n');

    pprintf(prn, "%s\n\n", _("Eigenvectors (component loadings)"));

    done = 0;
    todo = n;

    while (todo > 0) {
	int ncols = todo > PCA_COLS ? PCA_COLS : todo; 

	pprintf(prn, "%-*s", namelen + 1, " ");

	for (j=0; j<ncols; j++) {
	    sprintf(pcname, "PC%d", done + j + 1);
	    pprintf(prn, "%9s", pcname);
	}
	pputc(prn, '\n');

	for (i=0; i<n; i++) {
	    pprintf(prn, "%-*s", namelen + 1, cmat->names[i]);
	    for (j=0; j<ncols; j++) {
		pprintf(prn, "%9.3f", gretl_matrix_get(C, i, done + j));
	    }
	    pputc(prn, '\n');
	}
	pputc(prn, '\n');

	todo -= ncols;
	done += ncols;
    }
}

static int standardize (double *y, const double *x, int n)
{
    double xbar, sd;
    int i, err;

    err = gretl_moments(0, n-1, x, NULL, &xbar, &sd, NULL, NULL, 1);

    if (!err) {
	for (i=0; i<n; i++) {
	    y[i] = na(x[i]) ? NADBL : (x[i] - xbar) / sd;
	}
    }

    return err;
}

/* Add components to the dataset, either "major" ones (eigenvalues
   greater than 1.0), or a specified number (if @nsave > 0), or 
   all of them.
*/

static int pca_save_components (VMatrix *cmat, 
				gretl_matrix *E, 
				gretl_matrix *C, 
				DATASET *dset,
				int nsave,
				gretlopt opt)
{
    int save_all = (opt & OPT_A);
    double x, **sZ = NULL;
    int m = 0, v = dset->v;
    int k = cmat->dim;
    int i, j, t, vi;
    int err = 0;

    if (save_all) {
	m = k;
    } else if (nsave > 0) {
	m = nsave > k ? k : nsave;
    } else {
	for (i=0; E->val[i] > 1.0; i++) {
	    m++;
	} 
    }

    err = dataset_add_series(dset, m);

    if (!err) {
	/* construct standardized versions of all variables */
	sZ = doubles_array_new(k, dset->n); 
	if (sZ == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0; i<k && !err; i++) {
		vi = cmat->list[i+1];
		err = standardize(sZ[i], dset->Z[vi], dset->n);
	    }
	}
    }

    if (!err) {
	gchar *label;
	double load;

	for (i=0; i<m; i++) {
	    vi = v + i;
	    sprintf(dset->varname[vi], "PC%d", i+1);
	    make_varname_unique(dset->varname[vi], vi, dset);
	    label = g_strdup_printf(_("Component with eigenvalue = %.4f"),
				    E->val[i]);
	    series_set_label(dset, vi, label);
	    g_free(label);
	    for (t=0; t<dset->n; t++) {
		if (t < dset->t1 || t > dset->t2) {
		    dset->Z[vi][t] = NADBL;
		    continue;
		}
		dset->Z[vi][t] = 0.0;
		for (j=0; j<k; j++) {
		    x = sZ[j][t];
		    if (na(x)) {
			dset->Z[vi][t] = NADBL;
			break;
		    } else {
			load = gretl_matrix_get(C, j, i);
			dset->Z[vi][t] += load * x;
		    }
		}
	    }
	}
    }

    doubles_array_free(sZ, k);

    return err;
}

/* The incoming options here: 

   CLI: the option may be OPT_O (save the first XX components) or
   OPT_A (save all the components). OPT_Q suppresses printing of
   the results.

   GUI: no option (simply display the results) or OPT_D. The latter
   means that we should not display the results, but should put up a 
   dialog box allowing the user to decide what to save. 

   Note that depending on the original option supplied to the pca
   command, the incoming matrix may be either a correlation matrix
   (the default) or a covariance matrix.  This prior option is not
   included in the @opt passed here; rather it's encoded in the ci
   member of the VMatrix struct.
*/

int pca_from_cmatrix (VMatrix *cmat, DATASET *dset, 
		      gretlopt opt, PRN *prn)
{
    gretl_matrix *C;
    gretl_matrix *evals = NULL;
    gretlopt saveopt = opt;
    int k = cmat->dim;
    int nsave = 0;
    int i, j, idx;
    double x;
    int err = 0;

    if (opt & OPT_D) { 
	saveopt = pca_flag_dialog();
	if (saveopt == OPT_NONE) {
	    /* canceled */
	    return 0; 
	}
    } else if (opt & OPT_O) {
	nsave = get_optval_int(PCA, OPT_O, &err);
	if (err) {
	    return err;
	}
    }

    C = gretl_matrix_alloc(k, k);

    if (C == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<k; i++) {
	for (j=0; j<k; j++) {
	    idx = ijton(i, j, k);
	    x = cmat->vec[idx];
	    gretl_matrix_set(C, i, j, x);
	}
    }

#if PCA_DEBUG
    gretl_matrix_print(C, "original C, in pca");
#endif

    evals = gretl_symm_matrix_eigenvals_descending(C, 1, &err);

#if PCA_DEBUG
    gretl_matrix_print(C, "revised C (eigenvecs)");
#endif

    if (!err && !(opt & OPT_Q) && prn != NULL) {
	pca_print(cmat, evals, C, prn);
    }

    if (!err && saveopt) {
	err = pca_save_components(cmat, evals, C, dset, nsave,
				  saveopt);
    }

    gretl_matrix_free(evals);
    gretl_matrix_free(C);

    return err;
}
