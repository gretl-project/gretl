/*
 *  Copyright (c) by Ramu Ramanathan and Allin Cottrell
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

/* panel data plugin for gretl */

#include "libgretl.h"
#include <gtk/gtk.h>
#include "../gui/gretltypes.h"

typedef struct {
    int ns;
    double sigma_e;
    double H;
    double *bdiff;
    double *sigma;
} hausman_t;

/* .................................................................. */

static int get_panel_structure (DATAINFO *pdinfo, int *nunits, int *T)
{
    int err = 0;

    if (pdinfo->time_series == 2) {
	*nunits = pdinfo->n / pdinfo->pd;
	*T = pdinfo->pd;
    } 
    else if (pdinfo->time_series == 3) {
	char Tstr[8];

	if (sscanf(pdinfo->endobs, "%[^.].%d", Tstr, nunits) != 2)
	    err = 1;
	else 
	    *T = atoi(Tstr);
    } else err = 1;

    return err;
}

/* .................................................................. */

static void print_panel_const (MODEL *panelmod, print_t *prn)
{
    char numstr[18];
    int i = panelmod->list[0] - 1;

    sprintf(numstr, "(%.5g)", panelmod->sderr[i]);
    pprintf(prn, " constant: %14.5g %15s\n", 
	    panelmod->coeff[i], numstr);
}

/* .................................................................. */

static void print_panel_coeff (MODEL *pmod, MODEL *panelmod,
			       DATAINFO *pdinfo, int i, 
			       print_t *prn)
{
    char numstr[18];

    sprintf(numstr, "(%.5g)", panelmod->sderr[i]);
    pprintf(prn, "%9s: %14.5g %15s\n", 
	    pdinfo->varname[pmod->list[i+1]],
	    panelmod->coeff[i], numstr);
}

/* .................................................................. */

static double group_means_variance (MODEL *pmod, 
				    double *Z, DATAINFO *pdinfo,
				    double **groupZ, int nunits, int T) 
{
    int i, j, t, start, *list;
    double xx;
    MODEL meanmod;
    DATAINFO *ginfo;

    ginfo = create_new_dataset(groupZ, pmod->list[0], nunits, 0);
    if (ginfo == NULL) return NADBL;

    list = malloc((pmod->list[0] + 1) * sizeof *list);
    if (list == NULL) {
	clear_datainfo(ginfo, 1);
	free(ginfo);
	return NADBL;
    }

    list[0] = pmod->list[0];
    list[list[0]] = 0;
    for (j=1; j<list[0]; j++) { /* the variables */
	list[j] = j;
	start = 0;
	for (i=0; i<nunits; i++) { /* the observations */
	    xx = 0.0;
	    if (pdinfo->time_series == 2) {
		for (t=start; t<start+T; t++) 
		    xx += Z[pdinfo->n * pmod->list[j] + t];
		start += T;
	    } else {
		for (t=start; t<pdinfo->n; t += nunits) 
		    xx += Z[pdinfo->n * pmod->list[j] + t];
		start++;
	    }
	    xx /= (double) T;
	    (*groupZ)[nunits*j + i] = xx;
	}
    }

    meanmod = lsq(list, *groupZ, ginfo, OLS, 0, 0.0);
    if (meanmod.errcode) xx = NADBL;
    else xx = meanmod.sigma * meanmod.sigma;
    clear_model(&meanmod, NULL, NULL);
    clear_datainfo(ginfo, 1);
    free(ginfo);
    free(list);
    return xx;
}

/* .................................................................. */

static void vcv_slopes (hausman_t *haus, MODEL *pmod, int nunits, 
			int subt)
{
    int i, j, k = 0, min = 0, max;

    for (j=0; j<haus->ns; j++) {
	max = min + haus->ns - j;
	for (i=min; i<max; i++) {
	    /*  printf("vcv[%d] = %g\n", k, pmod->vcv[i]); */
	    if (subt) haus->sigma[k++] -= pmod->vcv[i];
	    else haus->sigma[k++] = pmod->vcv[i];
	}
	if (subt) min = max + 1;
	else min = max + nunits;
    }
}

/* .......................................................... */

#define TINY 1.0e-20

static int lu_decomp (double **a, int n, int *idx)
{
    int i, j, k, imax;
    double big, dum, sum, tmp;
    double *xx;

    xx = malloc((n + 1) * sizeof *xx);
    if (xx == NULL) return 1;
    for (i=0; i<=n; i++) xx[i] = 1.0;

    for (i=1; i<=n; i++) {
	big = 0.0;
	for (j=1; j<=n; j++) 
	    if ((tmp = fabs(a[i][j])) > big) big = tmp;
	if (floateq(big, 0.0)) {
	    free(xx);
	    return 1;
	}
	xx[i] = 1.0/big;
    }
    for (j=1; j<=n; j++) {
	for (i=1; i<j; i++) {
	    sum = a[i][j];
	    for (k=1; k<i; k++) sum -= a[i][k] * a[k][j];
	    a[i][j] = sum;
	}
	big = 0.0;
	for (i=j; i<=n; i++) {
	    sum = a[i][j];
	    for (k=1; k<j; k++) sum -= a[i][k] * a[k][j];
	    a[i][j] = sum;
	    if ((dum = xx[i] * fabs(sum)) >= big) {
		big = dum;
		imax = i;
	    }
	}
	if (j != imax) {
	    for (k=1; k<=n; k++) {
		dum = a[imax][k];
		a[imax][k] = a[j][k];
		a[j][k] = dum;
	    }
	    xx[imax] = xx[j];
	}
	idx[j] = imax;
	if (floateq(a[j][j], 0.0)) a[j][j] = TINY;
	if (j != n) {
	    dum = 1.0/a[j][j];
	    for (i=j+1; i<=n; i++) a[i][j] *= dum;
	}
    }
    free(xx);
    return 0;
}

/* .................................................................. */

static void lu_backsub (double **a, int n, int *idx, double *b)
{
    int i, k = 0, ip, j;
    double sum;

    for (i=1; i<=n; i++) {
	ip = idx[i];
	sum = b[ip];
	b[ip] = b[i];
	if (k) 
	    for (j=k; j<=i-1; j++) sum -= a[i][j] * b[j];
	else if (floatneq(sum, 0.0)) k = i;
	b[i] = sum;
    }
    for (i=n; i>=1; i--) {
	sum = b[i];
	for (j=i+1; j<=n; j++) sum -= a[i][j] * b[j];
	b[i] = sum / a[i][i];
    }
}

/* .................................................................. */

static double bXb (double *b, double **X, int n)
{
    int i, j;
    double row, xx = 0.0;

    for (i=1; i<=n; i++) {
	row = 0.0;
	for (j=1; j<=n; j++) row += b[j] * X[j][i];
	xx += b[i] * row;
    }
    return xx;
}

/* .................................................................. */

static int haus_invert (hausman_t *haus)
{
    double **a, **y, *col;
    int i, j, k, *idx;
    int err = 0, n = haus->ns;

    a = malloc((n + 1) * sizeof *a);
    if (a == NULL) return 1;
    for (i=1; i<=n; i++) {
	a[i] = malloc((n + 1) * sizeof **a);
	if (a[i] == NULL) return 1;
    }
    y = malloc((n + 1) * sizeof *y);
    if (y == NULL) return 1;
    for (i=1; i<=n; i++) {
	y[i] = malloc((n + 1) * sizeof **y);
	if (y[i] == NULL) return 1;
    }
    col = malloc((n + 1) * sizeof *col);
    if (col == NULL) return 1;
    idx = malloc((n + 1) * sizeof *idx);
    if (idx == NULL) return 1;

    k = 0;
    for (i=1; i<=n; i++) {
	for (j=i; j<=n; j++) {
	    a[i][j] = haus->sigma[k];
            if (i != j) 
	       a[j][i] = haus->sigma[k];
            k++;
	}
    }

    err = lu_decomp(a, n, idx);
    if (!err) {
	for (j=1; j<=n; j++) {
	    for (i=1; i<=n; i++) col[i] = 0.0;
	    col[j] = 1.0;
	    lu_backsub(a, n, idx, col);
	    for (i=1; i<=n; i++) y[i][j] = col[i];
	}
	haus->H = bXb(haus->bdiff, y, n);
    }

    for (i=1; i<=n; i++) {
	free(a[i]);
	free(y[i]);
    }
    free(a);
    free(y);
    free(col);
    free(idx);
    return err;
}

/* .................................................................. */

static double LSDV (MODEL *pmod, double **pZ, DATAINFO *pdinfo,
		    int nunits, int T, hausman_t *haus, print_t *prn) 
{
    int i, t, oldv = pdinfo->v, start;
    int *dvlist;
    double var, F;
    MODEL lsdv;

    dvlist = malloc((pmod->list[0] + nunits) * sizeof *dvlist);
    if (dvlist == NULL) return NADBL;
    if (dataset_add_vars(nunits - 1, pZ, pdinfo)) {
	free(dvlist);
	return NADBL;
    }

    start = 0;
    for (i=0; i<nunits-1; i++) {
	for (t=0; t<pdinfo->n; t++) 
	    (*pZ)[(oldv+i)*pdinfo->n + t] = 0.0;
	if (pdinfo->time_series == 2) {
	    for (t=start; t<start+T; t++) 
		(*pZ)[(oldv+i)*pdinfo->n + t] = 1.0;
	    start += T;
	} else {
	    for (t=start; t<pdinfo->n; t += nunits) 
		(*pZ)[(oldv+i)*pdinfo->n + t] = 1.0;
	    start++;
	}
    }

    dvlist[0] = pmod->list[0] + nunits - 1;
    for (i=1; i<=pmod->list[0]; i++) 
	dvlist[i] = pmod->list[i];
    for (i=1; i<nunits; i++) 
	dvlist[pmod->list[0] + i] = oldv + i - 1;

    lsdv = lsq(dvlist, *pZ, pdinfo, OLS, 0, 0.0);
    if (lsdv.errcode) {
	var = NADBL;
	pprintf(prn, "Error estimating fixed effects model\n");
	errmsg(lsdv.errcode, lsdv.errmsg, prn);
    } else {
	haus->sigma_e = lsdv.sigma;
	var = lsdv.sigma * lsdv.sigma;
	pprintf(prn, 
		"                          Fixed effects estimator\n"
		"          allows for differing intercepts by cross-sectional "
		"unit\n"
		"         (slope standard errors in parentheses, a_i = "
		"intercepts)\n\n");
	for (i=1; i<pmod->list[0] - 1; i++) {
	    print_panel_coeff(&lsdv, &lsdv, pdinfo, i, prn);
	    haus->bdiff[i] = lsdv.coeff[i];
	} for (i=pmod->list[0]; i<=dvlist[0]; i++) {
	    if (i < dvlist[0]) 
		lsdv.coeff[i-1] += lsdv.coeff[dvlist[0] - 1];
	    pprintf(prn, "      a_%d: %14.4g\n", 
		    i - pmod->list[0] + 1, lsdv.coeff[i-1]);
	}
	pprintf(prn, "\nResidual variance: %g/(%d - %d) = %g\n", 
		lsdv.ess, pdinfo->n, lsdv.ncoeff, var);
	F = (pmod->ess - lsdv.ess) * lsdv.dfd /
	    (lsdv.ess * (nunits - 1.0));
	pprintf(prn, "Joint significance of unit dummy variables:\n"
		" F(%d, %d) = %g with p-value %g\n", nunits - 1,
		lsdv.dfd, F, fdist(F, nunits - 1, lsdv.dfd));
	pprintf(prn, "(A low p-value counts against the null hypothesis that "
		"the pooled OLS model\nis adequate, in favor of the fixed "
		"effects alternative.)\n\n");
	makevcv(&lsdv);
	vcv_slopes(haus, &lsdv, nunits, 0);
    }
    clear_model(&lsdv, NULL, NULL);
    dataset_drop_vars(nunits - 1, pZ, pdinfo);
    free(dvlist);
    return var;
}

/* .................................................................. */

static int random_effects (MODEL *pmod, double *Z, DATAINFO *pdinfo, 
			   double *groupZ, double theta, int nunits, int T, 
			   hausman_t *haus, print_t *prn)
{
    double *reZ;
    DATAINFO *reinfo;
    MODEL remod;
    int *relist;
    int i, j, t, err = 0;

    reinfo = create_new_dataset(&reZ, pmod->list[0], pdinfo->n, 0);
    if (reinfo == NULL) return E_ALLOC;

    relist = malloc((pmod->list[0] + 1) * sizeof *relist);
    if (relist == NULL) {
	clear_datainfo(reinfo, 1);
	free(reinfo);
	free(reZ);
	return E_ALLOC;
    }

    relist[0] = pmod->list[0];
    relist[relist[0]] = 0;
    /* create transformed variables */
    for (i=1; i<relist[0]; i++) {
	relist[i] = i;
	j = 0;
	if (pdinfo->time_series == 2) { /* stacked time series */
	    for (t=0; t<pdinfo->n; t++) {
		if (t && (t % T == 0)) j++; 
		reZ[i*reinfo->n + t] = Z[pmod->list[i]*pdinfo->n + t] 
		    - theta * groupZ[i*nunits + j];
	    }
	} else { /* stacked cross sections */
	    for (t=0; t<pdinfo->n; t++) {
		if (t && t % nunits == 0) j = 0; /* FIXME ?? */
		reZ[i*reinfo->n + t] = Z[pmod->list[i]*pdinfo->n + t] 
		    - theta * groupZ[i*nunits + j];
		j++;
	    }
	}
    }
    for (t=0; t<pdinfo->n; t++) reZ[t] = 1.0 - theta;

    remod = lsq(relist, reZ, reinfo, OLS, 0, 0.0);
    if ((err = remod.errcode)) {
	pprintf(prn, "Error estimating random effects model\n");
	errmsg(err, remod.errmsg, prn);
    } else {
	pprintf(prn,
		"                         Random effects estimator\n"
		"           allows for a unit-specific component to the "
		"error term\n"
		"                     (standard errors in parentheses)\n\n");
	print_panel_const(&remod, prn);
	for (i=1; i<relist[0] - 1; i++) {
	    print_panel_coeff(pmod, &remod, pdinfo, i, prn);
	    haus->bdiff[i] -= remod.coeff[i];
	}
	makevcv(&remod);
	vcv_slopes(haus, &remod, nunits, 1);
    }
    clear_model(&remod, NULL, NULL);
    free(reZ);
    clear_datainfo(reinfo, 1);
    free(reinfo);
    free(relist);    

    return err;
}

/* .................................................................. */

int breusch_pagan_LM (MODEL *pmod, DATAINFO *pdinfo, 
		      int nunits, int T, print_t *prn)
{
    double *ubar, LM, eprime = 0.0;
    int i, t, start = 0;

    ubar = malloc(nunits * sizeof *ubar);
    if (ubar == NULL) return E_ALLOC;

    for (i=0; i<nunits; i++) {
	ubar[i] = 0.0;
	if (pdinfo->time_series == 2) {
	    for (t=start; t<start+T; t++) 
		ubar[i] += pmod->uhat[t];
	    start += T;
	} else {
	    for (t=start; t<pdinfo->n; t += nunits) 
		ubar[i] += pmod->uhat[t];
	    start++;
	}
	ubar[i] /= (double) T;
	eprime += ubar[i] * ubar[i];
    }

    pprintf(prn, "\nMeans of pooled OLS residuals for cross-sectional "
	    "units:\n\n");
    for (i=0; i<nunits; i++) {
	pprintf(prn, " unit %2d: %13.5g\n", 
		i + 1, ubar[i]);
    }
    free(ubar);

    LM = (double) pdinfo->n/(2.0*(T - 1.0)) * 
	pow((T * T * eprime/pmod->ess) - 1.0, 2);
    pprintf(prn, "\nBreusch-Pagan test statistic:\n"
	    " LM = %g with p-value = prob(chi-square(1) > %g) = %g\n", 
	    LM, LM, chisq(LM, 1));
    pprintf(prn, "(A low p-value counts against the null hypothesis that "
	    "the pooled OLS model\nis adequate, in favor of the random "
	    "effects alternative.)\n\n");
    return 0;
}

/* .................................................................. */

static int hausman_test (hausman_t *haus, print_t *prn)
{
/*      int i, ns = haus->ns; */
/*      int nterms = (ns * ns + ns) / 2; */

/*      for (i=1; i<=ns; i++)  */
/*  	pprintf(prn, "b%d_FE - beta%d_RE = %g\n", i, i, haus->bdiff[i]); */
/*      pprintf(prn, "\n"); */

/*      for (i=0; i<nterms; i++)  */
/*  	pprintf(prn, "vcv_diff[%d] = %g\n", i, haus->sigma[i]); */

    if (haus_invert(haus)) { 
	pprintf(prn, "Error attempting to invert vcv difference matrix\n");
	return 1;
    }
    if (haus->H < 0) 
	pprintf(prn, "\nHausman test matrix is not positive definite (this "
		"result may be treated as\n\"fail to reject\" the random effects "
		"specification).\n");
    else {
	pprintf(prn, "\nHausman test statistic:\n"
		" H = %g with p-value = prob(chi-square(%d) > %g) = %g\n",
		haus->H, haus->ns, haus->H, chisq(haus->H, haus->ns));
	pprintf(prn, "(A low p-value counts against the null hypothesis that "
		"the random effects\nmodel is consistent, in favor of the fixed "
		"effects model.)\n");
    }

    return 0;
}

/* .................................................................. */

int panel_diagnostics (MODEL *pmod, double **pZ, DATAINFO *pdinfo, 
		       print_t *prn)
{
    int nunits, ns, T;
    double var1, var2, theta;
    double *groupZ;
    hausman_t haus;

    if (get_panel_structure(pdinfo, &nunits, &T))
	return 1;

    if (nunits > pmod->ncoeff) {
	ns = haus.ns = pmod->ncoeff - 1;
	haus.bdiff = malloc(pmod->ncoeff * sizeof(double));
	if (haus.bdiff == NULL) return E_ALLOC;
	haus.sigma = malloc(((ns * ns + ns) / 2) * sizeof(double));
	if (haus.sigma == NULL) return E_ALLOC; 
    }   
    
/*      fprintf(stderr, "diagnostics: printing to buffer\n"); */
    pprintf(prn, "      Diagnostics: assuming a balanced panel with %d "
	    "cross-sectional units\n "
	    "                        observed over %d periods\n\n", 
	    nunits, T);

    var2 = LSDV(pmod, pZ, pdinfo, nunits, T, &haus, prn);

    breusch_pagan_LM(pmod, pdinfo, nunits, T, prn);
    
    if (nunits > pmod->ncoeff && var2 > 0) {
	var1 = group_means_variance(pmod, *pZ, pdinfo, &groupZ, nunits, T);
	if (var1 < 0) 
	    pprintf(prn, "Couldn't estimate group means regression\n");
	else {
	    pprintf(prn, "Residual variance for group means "
		    "regression: %g\n\n", var1);    
	    theta = 1.0 - sqrt(var2 / (T * var1));
	    random_effects(pmod, *pZ, pdinfo, groupZ, theta, nunits, T, 
			   &haus, prn);
	    hausman_test(&haus, prn);
	}
	free(groupZ);
	free(haus.bdiff);
	free(haus.sigma);
    }

    return 0;
}

/* .................................................................. */

static void set_panel_structure (GtkWidget *w, gpointer data)
{
    gint i;

    if (GTK_TOGGLE_BUTTON (w)->active) {
	DATAINFO *pdinfo = (DATAINFO *) data;
	i = GPOINTER_TO_INT(gtk_object_get_data(GTK_OBJECT(w), "action"));
	if (dataset_is_panel(pdinfo))
	    pdinfo->time_series = (short) i;
    }
}

/* .................................................................. */

void panel_structure_dialog (DATAINFO *pdinfo, GtkWidget *w,
			     void (*cleanfun)(), void (*helpfun)())
{
    dialog_t *d, *cancel_d;
    GtkWidget *button;
    GtkWidget *tempwid;
    GSList *group;

    d = malloc(sizeof *d);
    if (d == NULL) return;
    cancel_d = malloc(sizeof *cancel_d);
    if (cancel_d == NULL) {
	free(d);
	return;
    }
    
    d->data = cancel_d->data = NULL;
    cancel_d->all_buttons = d->all_buttons = NULL;

    d->dialog = gtk_dialog_new();
    w = d->dialog;

    gtk_window_set_title (GTK_WINDOW (d->dialog), "gretl: panel data structure");
    gtk_window_set_policy (GTK_WINDOW (d->dialog), FALSE, FALSE, FALSE);
    gtk_container_border_width (GTK_CONTAINER 
				(GTK_DIALOG (d->dialog)->vbox), 10);
    gtk_container_border_width (GTK_CONTAINER 
				(GTK_DIALOG (d->dialog)->action_area), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 5);
    gtk_box_set_spacing (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 15);
    gtk_box_set_homogeneous (GTK_BOX 
			     (GTK_DIALOG (d->dialog)->action_area), TRUE);
    gtk_window_set_position (GTK_WINDOW (d->dialog), GTK_WIN_POS_MOUSE);

    gtk_signal_connect (GTK_OBJECT (d->dialog), "destroy", 
			GTK_SIGNAL_FUNC (cleanfun), 
			cancel_d);

    button = gtk_radio_button_new_with_label (NULL, "Stacked time series");
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			button, TRUE, TRUE, FALSE);
    if (pdinfo->time_series == 2)
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
                       GTK_SIGNAL_FUNC(set_panel_structure), pdinfo);
    gtk_object_set_data(GTK_OBJECT(button), "action", GINT_TO_POINTER(2));
    gtk_widget_show (button);

    group = gtk_radio_button_group (GTK_RADIO_BUTTON (button));
    button = gtk_radio_button_new_with_label(group, "Stacked cross sections");
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->vbox), 
			button, TRUE, TRUE, FALSE);
    if (pdinfo->time_series == 3)
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (button), TRUE);
    gtk_signal_connect(GTK_OBJECT(button), "clicked",
                       GTK_SIGNAL_FUNC(set_panel_structure), pdinfo);
    gtk_object_set_data(GTK_OBJECT(button), "action", GINT_TO_POINTER(3));
    gtk_widget_show (button);

    /* Create the "OK" button */
    tempwid = gtk_button_new_with_label ("OK");
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 
			tempwid, TRUE, TRUE, FALSE);
    gtk_signal_connect_object (GTK_OBJECT (tempwid), "clicked", 
			       GTK_SIGNAL_FUNC (gtk_widget_destroy), 
			       GTK_OBJECT (d->dialog));
    gtk_widget_grab_default (tempwid);
    gtk_widget_show (tempwid);

    /* Create the "Cancel" button */
    tempwid = gtk_button_new_with_label ("Cancel");
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 
			tempwid, TRUE, TRUE, FALSE);
    gtk_signal_connect_object (GTK_OBJECT (tempwid), "clicked", 
			       GTK_SIGNAL_FUNC (gtk_widget_destroy), 
			       GTK_OBJECT (d->dialog));
    gtk_widget_show (tempwid);

    /* Create a "Help" button */
    tempwid = gtk_button_new_with_label ("Help");
    GTK_WIDGET_SET_FLAGS (tempwid, GTK_CAN_DEFAULT);
    gtk_box_pack_start (GTK_BOX (GTK_DIALOG (d->dialog)->action_area), 
			tempwid, TRUE, TRUE, FALSE);
    gtk_signal_connect (GTK_OBJECT (tempwid), "clicked", 
			GTK_SIGNAL_FUNC (helpfun), 
			GINT_TO_POINTER (PANEL));
    gtk_widget_show (tempwid);

    gtk_widget_show (d->dialog);
    gtk_main();
}




