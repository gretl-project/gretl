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
#include "gretl_matrix.h"

#include <gtk/gtk.h>

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
					     _("Components with eigenvalues > mean"));
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

    return OPT_NONE;
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

    nl = g_utf8_strlen(_("Variable"), -1);
    if (nl > namelen) {
	namelen = nl;
    }

    done = 0;
    todo = n;

    while (todo > 0) {
	int ncols = (todo > PCA_COLS)? PCA_COLS : todo; 

	pprintf(prn, "%-*s", namelen + 1, _("Variable"));

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

/* add PCs to the dataset, either "major" ones or all */

static int pca_save_components (VMatrix *cmat, 
				gretl_matrix *E, gretl_matrix *C, 
				double ***pZ, DATAINFO *pdinfo,
				gretlopt opt)
{
    int save_all = (opt & OPT_A);
    double **sZ = NULL;
    double x;
    int *plist = NULL;
    int m, v = pdinfo->v;
    int k = cmat->dim;
    int i, j, t, vi;
    int err = 0;

    if (save_all) {
	m = k;
    } else {
	m = 0;
	for (i=0; i<k; i++) {
	    if (E->val[i] > 1.0) {
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
	for (i=0; i<k; i++) {
	    if (save_all || E->val[i] > 1.0) {
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
		vi = cmat->list[i+1];
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
		    "eigenvalue = %.4f", E->val[pi]);

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

    return err;
}

/* The incoming option here: When this function is called from the
   CLI, the option may be OPT_O (save the first component), or OPT_A
   (save all the components), or none.  The results are printed in all
   cases.  As for the GUI, we either get no option (simply display the
   results) or OPT_D.  The latter means that we should not display the
   results, but should put up a dialog box allowing the user to decide
   what to save.

   Note that depending on the original option supplied to the "pca"
   command, the incoming matrix may be either a correlation matrix
   (the default) or a covariance matrix.  This prior option is not
   included in the gretlopt passed here; it's encoded in the "ci"
   member of the VMatrix struct.
 */

int pca_from_cmatrix (VMatrix *cmat, double ***pZ,
		      DATAINFO *pdinfo, gretlopt opt,
		      PRN *prn)
{
    gretl_matrix *C;
    gretl_matrix *evals = NULL;
    gretlopt saveopt = opt;
    int k = cmat->dim;
    int i, j, idx;
    double x;
    int err = 0;

    if (opt & OPT_D) { 
	saveopt = pca_flag_dialog();
	if (saveopt == OPT_NONE) {
	    /* canceled */
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
	    x = cmat->vec[idx];
	    gretl_matrix_set(C, i, j, x);
	}
    }

#if PCA_DEBUG
    gretl_matrix_print(C, "original C, in pca");
#endif

    evals = gretl_symmetric_matrix_eigenvals(C, 1, &err);

    if (!err) {
	err = gretl_symmetric_eigen_sort(evals, C, 0);
    }

    if (!err && prn != NULL) {
	pca_print(cmat, evals, C, prn);
    }

    if (!err && saveopt) {
	err = pca_save_components(cmat, evals, C, pZ, pdinfo, saveopt);
    }

    gretl_matrix_free(evals);
    gretl_matrix_free(C);

    return err;
}
