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

static int XTX_properties (const MODEL *pmod, const double **Z,
			   PRN *prn)
{
    double *xpx = NULL;
    int k = pmod->ncoeff;
    double xnorm, rcond, det = 1;
    int err = 0;

    xpx = gretl_XTX(pmod, Z, &err);

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

static double get_vif (const MODEL *pmod, double ***pZ, 
		       DATAINFO *pdinfo, int *vlist, int k,
		       int *err)
{
    MODEL tmpmod;
    double x = NADBL;
    int i, j;

    vlist[1] = pmod->list[k];
    j = 2;
    for (i=2; i<=pmod->list[0]; i++) {
	if (i != k) {
	    vlist[j++] = pmod->list[i];
	}
    }

    tmpmod = lsq(vlist, pZ, pdinfo, OLS, OPT_A); 
    *err = tmpmod.errcode;

    if (!*err && !na(tmpmod.rsq) && tmpmod.rsq != 1.0) {
	x = 1.0 / (1.0 - tmpmod.rsq);
    }

    clear_model(&tmpmod);

    return x;
}

static int testlist (const int *list)
{
    int i;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    return 1;
	}
    }

    return 0;
}

static double *
model_vif_vector (MODEL *pmod, double ***pZ, DATAINFO *pdinfo,
		  int *err)
{
    double *vif = NULL;
    int *vlist = NULL;
    int t1 = pdinfo->t1, t2 = pdinfo->t2;
    int nvif = pmod->ncoeff - pmod->ifc;
    int m = pmod->list[0] - 1;
    int i, j;

    if (nvif <= 1) {
	gretl_errmsg_set(_("The statistic you requested is not meaningful "
			   "for this model"));
	return NULL;
    }

    if (testlist(pmod->list)) {
	*err = E_DATA;
	return NULL;
    }

    vif = malloc(nvif * sizeof *vif);
    if (vif == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    vlist = gretl_list_new(m);
    if (vlist == NULL) {
	*err = E_ALLOC;
	free(vif);
	return NULL;
    }

    /* impose original model sample */
    pdinfo->t1 = pmod->t1;
    pdinfo->t2 = pmod->t2;

    j = 0;
    for (i=2; i<=pmod->list[0] && !*err; i++) {
	if (pmod->list[i] != 0) {
	    vif[j++] = get_vif(pmod, pZ, pdinfo, vlist, i, err);
	}
    }

    /* reinstate sample */
    pdinfo->t1 = t1;
    pdinfo->t2 = t2;

    free(vlist);

    if (*err) {
	free(vif);
	vif = NULL;
    }

    return vif;
}

#define xtx_ok(c) (c == OLS || c == AR1 || \
		   c == WLS || c == HSK) 

int print_vifs (MODEL *pmod, double ***pZ, DATAINFO *pdinfo, 
		PRN *prn)
{
    double *vif = NULL;
    int v, i, j;
    int err = 0;

    vif = model_vif_vector(pmod, pZ, pdinfo, &err);
    if (err) {
	return err;
    }

    pprintf(prn, "%s\n\n", _("Variance Inflation Factors"));

    pprintf(prn, " %s\n", _("Minimum possible value = 1.0"));
    pprintf(prn, " %s\n", _("Values > 10.0 may indicate a collinearity problem"));
    pputc(prn, '\n');

    j = 0;
    for (i=2; i<=pmod->list[0]; i++) {
	v = pmod->list[i];
	if (v != 0) {
	    pprintf(prn, " %3d) %15s %8.3f\n", v, pdinfo->varname[v], vif[j++]);
	}
    }
    pputc(prn, '\n');

    pputs(prn, _("VIF(j) = 1/(1 - R(j)^2), where R(j) is the "
		 "multiple correlation coefficient\nbetween "
		 "variable j and the other independent variables"));
    pputc(prn, '\n');

    if (xtx_ok(pmod->ci)) {
	XTX_properties(pmod, (const double **) *pZ, prn);
    }

    free(vif);

    return 0;
}
