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

#include "gretl_f2c.h"
#include "clapack_double.h"

static double packed_matrix_norm (const double *x, int k)
{
    double csum, cmax = 0.0;
    int i, j;

    for (j=0; j<k; j++) {
	csum = 0.0;
	for (i=0; i<k; i++) {
	    csum += fabs(x[ijton(i, j, k)]);
	}
	if (csum > cmax) {
	    cmax = csum;
	}
    }

    return cmax;
}

/* Get 1-norm, determinant and reciprocal condition number using
   Cholesky */

static int 
decomp_etc (double *xpx, int k, double *xnorm, double *det, double *rcond)
{
    char uplo = 'L';
    integer n = k;
    integer info = 0;
    integer *iwork = NULL;
    double *work = NULL;
    int i, err = 0;

    work = malloc((3 * n) * sizeof *work);
    iwork = malloc(n * sizeof *iwork);

    if (work == NULL || iwork == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    *xnorm = packed_matrix_norm(xpx, k);

    dpptrf_(&uplo, &n, xpx, &info);

    if (info != 0) {
	err = 1;
    } else {
	double d = 1.0;

	for (i=0; i<k; i++) {
	    d *= xpx[ijton(i,i,k)];
	}

	*det = d * d;
	dppcon_(&uplo, &n, xpx, xnorm, rcond, work, iwork, &info);
	if (info != 0) {
	    err = 1;
	}
    } 

 bailout:

    free(work);
    free(iwork);

    return err;
}

static int XTX_properties (const MODEL *pmod, DATASET *dset,
			   PRN *prn)
{
    double *xpx = NULL;
    int k = pmod->ncoeff;
    double xnorm, rcond, det = 1;
    int err = 0;

    xpx = gretl_XTX(pmod, dset, &err);

    if (!err) {
	err = decomp_etc(xpx, k, &xnorm, &det, &rcond);
    }

    if (!err) {
	pprintf(prn, "\n%s:\n\n", _("Properties of matrix X'X"));
	pprintf(prn, " %s = %.8g\n", _("1-norm"), xnorm);
	pprintf(prn, " %s = %.8g\n", _("Determinant"), det);
	pprintf(prn, " %s = %.8g\n", _("Reciprocal condition number"), rcond);
	pputc(prn, '\n');
    }

    free(xpx);

    return err;
}

/* run the vif regression for regressor k */

static double get_vif (MODEL *mod, const int *xlist, 
		       int *vlist, int k,
		       DATASET *dset,
		       int *err)
{
    double vk = NADBL;
    int i, j;

    vlist[1] = xlist[k]; /* dep. var. is regressor k */
    /* position 2 in vlist holds 0 = const */
    j = 3;
    for (i=1; i<=xlist[0]; i++) {
	if (i != k) {
	    vlist[j++] = xlist[i];
	}
    }

    *mod = lsq(vlist, dset, OLS, OPT_A); 
    *err = mod->errcode;

    if (!*err && !xna(mod->rsq) && mod->rsq != 1.0) {
	vk = 1.0 / (1.0 - mod->rsq);
    }

    clear_model(mod);

    return vk;
}

/* run regressions of each x_i on the other x_j's */

static double *model_vif_vector (MODEL *pmod, const int *xlist,
				 DATASET *dset, int *err)
{
    MODEL tmpmod;
    double *vif = NULL;
    int *vlist = NULL;
    int nvif = xlist[0];
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int i;

    if (nvif <= 1) {
	gretl_errmsg_set(_("The statistic you requested is not meaningful "
			   "for this model"));
	return NULL;
    }

    vif = malloc(nvif * sizeof *vif);
    if (vif == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* vlist is the list for the vif regressions:
       allow space for the constant */
    vlist = gretl_list_new(nvif + 1);
    if (vlist == NULL) {
	*err = E_ALLOC;
	free(vif);
	return NULL;
    }

    /* impose original model sample */
    dset->t1 = pmod->t1;
    dset->t2 = pmod->t2;

    for (i=1; i<=xlist[0] && !*err; i++) {
	vif[i-1] = get_vif(&tmpmod, xlist, vlist, i, dset, err);
    }

    /* reinstate sample */
    dset->t1 = save_t1;
    dset->t2 = save_t2;

    free(vlist);

    if (*err) {
	free(vif);
	vif = NULL;
    }

    return vif;
}

#define xtx_ok(c) (c == OLS || c == AR1 || c == WLS) 

int print_vifs (MODEL *pmod, DATASET *dset, PRN *prn)
{
    double *vif = NULL;
    int *xlist;
    double vj;
    int vi, i, n;
    int maxlen = 0;
    int err = 0;

    /* fetch list of regressors */
    xlist = gretl_model_get_x_list(pmod);
    if (xlist == NULL) {
	return E_DATA;
    }

    /* drop the constant if present in xlist */
    for (i=1; i<=xlist[0]; i++) {
	if (xlist[i] == 0) {
	    gretl_list_delete_at_pos(xlist, i);
	    break;
	}
    }

    vif = model_vif_vector(pmod, xlist, dset, &err);
    if (err) {
	return err;
    }

    pprintf(prn, "%s\n\n", _("Variance Inflation Factors"));

    pprintf(prn, "%s\n", _("Minimum possible value = 1.0"));
    pprintf(prn, "%s\n", _("Values > 10.0 may indicate a collinearity problem"));
    pputc(prn, '\n');

    for (i=1; i<=xlist[0]; i++) {
	vi = xlist[i];
	vj = vif[i-1];
	if (!na(vj)) {
	    n = strlen(dset->varname[vi]);
	    if (n > maxlen) {
		maxlen = n;
	    }
	}
    }

    maxlen = maxlen < 12 ? 12 : maxlen;

    for (i=1; i<=xlist[0]; i++) {
	vi = xlist[i];
	vj = vif[i-1];
	if (!na(vj)) {
	    pprintf(prn, "%*s %8.3f\n", maxlen, dset->varname[vi], vj);
	}
    }
    pputc(prn, '\n');

    pputs(prn, _("VIF(j) = 1/(1 - R(j)^2), where R(j) is the "
		 "multiple correlation coefficient\nbetween "
		 "variable j and the other independent variables"));
    pputc(prn, '\n');

    if (xtx_ok(pmod->ci)) {
	XTX_properties(pmod, dset, prn);
    }

    free(vif);
    free(xlist);

    return 0;
}
