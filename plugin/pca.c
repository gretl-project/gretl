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
		       double *evals, PRN *prn)
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
	y += evals[i] / n;
	pprintf(prn, "%5d%13.4f%13.4f%13.4f\n", n - i,
		evals[i], evals[i] / n, y);
	x += evals[i];
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

static double *standardize (const double *x, int n)
{
    double *sx;
    double xbar, sd;
    int i, err;

    err = gretl_moments(0, n-1, x, &xbar, &sd, NULL, NULL, 1);
    if (err) {
	return NULL;
    }

    sx = malloc(n * sizeof *sx);
    if (sx == NULL) {
	return NULL;
    }

    for (i=0; i<n; i++) {
	if (na(x[i])) {
	    sx[i] = NADBL;
	} else {
	    sx[i] = (x[i] - xbar) / sd;
	}
    }

    return sx;
}

int pca_from_corrmat (VMatrix *corrmat, double ***pZ,
		      DATAINFO *pdinfo, gretlopt *pflag,
		      PRN *prn)
{
    gretl_matrix *m;
    double x;
    int i, j, idx, n = corrmat->dim;
    double *evals;
    gretlopt oflag = 0L;
    int err = 0;

    if (pflag != NULL) oflag = *pflag;

    if (oflag & OPT_D) { 
	oflag = pca_flag_dialog();
	if (!oflag) {
	    /* canceled */
	    *pflag = 0L;
	    return 0; 
	}
    }    

    m = gretl_matrix_alloc(n, n);
    if (m == NULL) return E_ALLOC;

    for (i=0; i<n; i++) {
	for (j=0; j<n; j++) {
	    idx = ijton(i, j, n);
	    x = corrmat->vec[idx];
	    gretl_matrix_set(m, i, j, x);
	}
    }

    evals = gretl_symmetric_matrix_eigenvals(m, 1, &err);
    if (err) {
	gretl_matrix_free(m);
	return err;
    }

    if (prn != NULL) {
	pca_print(corrmat, m, evals, prn);
    }

    if (oflag) {
	/* add components with eigenvalues > 1 to the dataset */
	int v = pdinfo->v;
	int nc = 0, err = 0;
	double **sZ = NULL;
	int add_all = (oflag == OPT_A);
	int *plist;

	if (add_all) {
	    nc = n;
	} else {
	    for (i=0; i<n; i++) {
		if (evals[i] > 1.0) nc++;
	    }
	}

	plist = malloc((nc + 1) * sizeof *plist);
	if (plist == NULL) err = E_ALLOC;

	if (!err) {
	    /* build list of PCs (with eigenvals > 1?) */
	    plist[0] = nc;
	    j = 1;
	    for (i=n-1; i>=0; i--) {
		if (add_all || evals[i] > 1.0) {
		    plist[j++] = i;
		}
	    }
#ifdef PCA_DEBUG
	    printlist(plist, "pclist");
#endif
	    err = dataset_add_series(nc, pZ, pdinfo);
	}

	if (!err) {
	    /* construct standardized versions of variables */
	    sZ = malloc(n * sizeof *sZ);
	    if (sZ == NULL) {
		err = E_ALLOC;
	    } else {
		for (i=0; i<n; i++) {
		    sZ[i] = NULL;
		}
		for (i=0; i<n; i++) {
		    int oldv = corrmat->list[i+1];

#ifdef PCA_DEBUG
		    fprintf(stderr, "Getting standardized version of "
			    "var %d\n", oldv);
#endif
		    sZ[i] = standardize((const double *) (*pZ)[oldv], 
					pdinfo->n);
		    if (sZ[i] == NULL) {
			err = E_ALLOC;
			break;
		    }
		}
		if (err) {
		    for (i=0; i<n; i++) {
			free(sZ[i]);
		    }
		    free(sZ);
		    sZ = NULL;
		}
	    }
	}

	if (!err) {
	    for (i=1; i<=plist[0]; i++) {
		int newv = v + i - 1;
		int pcnum = plist[i];
		int t;

		sprintf(pdinfo->varname[newv], "PC%d", i);
		make_varname_unique(pdinfo->varname[newv], newv, pdinfo);
		sprintf(VARLABEL(pdinfo, newv), "Component with "
			"eigenvalue = %.4f", evals[pcnum]);
		for (t=0; t<pdinfo->n; t++) {
#ifdef PCA_DEBUG
		    fprintf(stderr, "Obs %d\n", t);
#endif
		    (*pZ)[newv][t] = 0.0;
		    for (j=0; j<n; j++) {
			double load = gretl_matrix_get(m, j, pcnum);
			double val = sZ[j][t];

#ifdef PCA_DEBUG
			fprintf(stderr, "j=%d,pcnum=%d,load=%g,val=%g\n",
				j,pcnum,load,val);
#endif

			if (na(val)) {
			    (*pZ)[newv][t] = NADBL;
			    break;
			} else {
			    (*pZ)[newv][t] += load * val;
			}
		    }
		} /* end loop over observations */
	    } /* end loop over components */
	} /* end !err conditional */

	free(plist);
	if (sZ != NULL) {
	    for (i=0; i<n; i++) {
		free(sZ[i]);
	    }
	    free(sZ);
	}

    } /* end oflag conditional */

    free(evals);
    gretl_matrix_free(m);

    if (pflag != NULL) *pflag = oflag;

    return 0;
}
