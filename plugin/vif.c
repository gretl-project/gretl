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
#include "matrix_extra.h"
#include "version.h"

#define NAMEWIDTH 9

static void BKW_print (gretl_matrix *B, int namelen, PRN *prn);

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

    if (!*err && !na(mod->rsq) && mod->rsq != 1.0) {
	vk = 1.0 / (1.0 - mod->rsq);
    }

    clear_model(mod);

    return vk;
}

/* run regressions of each x_i on the other x_j's */

static gretl_vector *model_vif_vector (MODEL *pmod, const int *xlist,
				       DATASET *dset, int *err)
{
    MODEL tmpmod;
    gretl_vector *vif = NULL;
    int *vlist = NULL;
    int nvif = xlist[0];
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int i;

    vif = gretl_column_vector_alloc(nvif);
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
	vif->val[i-1] = get_vif(&tmpmod, xlist, vlist, i, dset, err);
    }

    /* reinstate sample */
    dset->t1 = save_t1;
    dset->t2 = save_t2;

    free(vlist);

    if (*err) {
	gretl_matrix_free(vif);
	vif = NULL;
    }

    return vif;
}

static int bkw_add_colnames (gretl_matrix *BKW,
			     const gretl_matrix *VCV,
			     gretl_array *pnames)
{
    char **S = strings_array_new(BKW->cols);
    const char **S0 = NULL;
    int maxlen = 0;

    if (pnames == NULL) {
	S0 = gretl_matrix_get_colnames(VCV);
    }

    if (S != NULL) {
	int i, len, k = BKW->cols - 2;
	const char *si;
	char tmp[16];

	S[0] = gretl_strdup("lambda");
	S[1] = gretl_strdup("cond");

	for (i=0; i<k; i++) {
	    if (pnames != NULL || S0 != NULL) {
		si = pnames != NULL ? gretl_array_get_data(pnames, i) : S0[i];
		if (strlen(si) > NAMEWIDTH) {
		    *tmp = '\0';
		    strncat(tmp, si, NAMEWIDTH - 1);
		    strcat(tmp, "~");
		    S[i+2] = gretl_strdup(tmp);
		} else if (pnames != NULL) {
		    S[i+2] = gretl_array_get_data(pnames, i);
		    gretl_array_set_data(pnames, i, NULL);
		} else {
		    S[i+2] = gretl_strdup(si);
		}
	    } else {
		sprintf(tmp, "x%d", i+1);
		S[i+2] = gretl_strdup(tmp);
	    }
	    len = strlen(S[i+2]);
	    if (len > maxlen) {
		maxlen = len;
	    }
	}
	gretl_matrix_set_colnames(BKW, S);
    }

    return maxlen;
}

/* note: we're assuming in bkw_matrix() that the array argument
   @pnames is disposable: we pull out its entries and set them
   to NULL, but it's still up to the caller to destroy the
   array itself.
*/

gretl_matrix *bkw_matrix (const gretl_matrix *VCV,
			  gretl_array *pnames,
			  PRN *prn, int *err)
{
    gretl_matrix *Vi = NULL;
    gretl_matrix *S = NULL;
    gretl_matrix *Q = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *lambda = NULL;
    gretl_matrix *BKW = NULL;
    double x, y;
    int k = VCV->rows;
    int namelen = 0;
    int i, j;

    if (pnames != NULL) {
	if (gretl_array_get_length(pnames) != k) {
	    fprintf(stderr, "bkw_matrix: expected %d names but got %d\n",
		    k, gretl_array_get_length(pnames));
	    *err = E_INVARG;
	    return NULL;
	}
    }

    /* copy the covariance matrix */
    Vi = gretl_matrix_copy(VCV);
    if (Vi == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* and invert it */
    *err = gretl_invert_symmetric_matrix(Vi);
    if (*err) {
	goto bailout;
    }

    /* allocate workspace */
    S = gretl_identity_matrix_new(k);
    Q = gretl_matrix_alloc(k, k);
    BKW = gretl_matrix_alloc(k, k+2);

    if (S == NULL || Q == NULL || BKW == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }

    for (i=0; i<k; i++) {
	x = gretl_matrix_get(Vi, i, i);
	gretl_matrix_set(S, i, i, 1/sqrt(x));
    }

    *err = gretl_matrix_qform(S, GRETL_MOD_TRANSPOSE,
			      Vi, Q, GRETL_MOD_NONE);

    if (!*err) {
	*err = gretl_matrix_SVD(Q, NULL, &lambda, &V, 0);
    }

    if (*err) {
	goto bailout;
    }

    /* S = (1/lambda) ** ones(k, 1) */
    for (j=0; j<k; j++) {
	x = lambda->val[j];
	for (i=0; i<k; i++) {
	    gretl_matrix_set(S, i, j, 1/x);
	}
    }

    for (i=0; i<k; i++) {
	for (j=0; j<k; j++) {
	    x = gretl_matrix_get(V, j, i);
	    y = gretl_matrix_get(S, i, j);
	    gretl_matrix_set(Q, i, j, x * x * y);
	}
    }

    for (i=0; i<k; i++) {
	/* compute row sums */
	y = 0.0;
	for (j=0; j<k; j++) {
	    y += gretl_matrix_get(Q, i, j);
	}
	for (j=0; j<k; j++) {
	    x = gretl_matrix_get(Q, i, j);
	    gretl_matrix_set(V, j, i, x/y);
	}
    }

    y = lambda->val[0];

    /* assemble the matrix to return */
    for (i=0; i<k; i++) {
	x = lambda->val[i];
	gretl_matrix_set(BKW, i, 0, x);
	gretl_matrix_set(BKW, i, 1, sqrt(y / x));
	for (j=0; j<k; j++) {
	    x = gretl_matrix_get(V, i, j);
	    gretl_matrix_set(BKW, i, j+2, x);
	}
    }

    namelen = bkw_add_colnames(BKW, VCV, pnames);

 bailout:

    gretl_matrix_free(Vi);
    gretl_matrix_free(S);
    gretl_matrix_free(Q);
    gretl_matrix_free(V);
    gretl_matrix_free(lambda);

    if (*err) {
	gretl_matrix_free(BKW);
	BKW = NULL;
    } else if (prn != NULL) {
	BKW_print(BKW, namelen, prn);
    }

    return BKW;
}

static int do_proportion_sums (const gretl_matrix *B,
			       const char **bnames,
			       const char *label,
			       double cval, PRN *prn)
{
    gretl_matrix *P;
    char **cnames;
    double x;
    int np = B->cols - 2;
    int len, namelen = 0;
    int ngp5 = 0;
    int i, j, k;

    cnames = strings_array_new(np);
    if (cnames == NULL) {
	return E_ALLOC;
    }

    P = gretl_zero_matrix_new(1, np);
    if (P == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<B->rows; i++) {
	if (gretl_matrix_get(B, i, 1) >= cval) {
	    for (j=2; j<B->cols; j++) {
		x = 0;
		for (k=i; k<B->rows; k++) {
		    x += gretl_matrix_get(B, k, j);
		}
		if (x >= 0.5) {
		    P->val[ngp5] = x;
		    cnames[ngp5] = gretl_strdup(bnames[j]);
		    len = strlen(cnames[ngp5]);
		    if (len > namelen) {
			namelen = len;
		    }
		    ngp5++;
		}
	    }
	    break;
	}
    }

    if (ngp5 > 0) {
	char fmt[16];
	int len = 8;

	if (len < namelen + 1) {
	    len = namelen + 1;
	}
	sprintf(fmt, "%%%d.3f", len);
	P->cols = ngp5;
	gretl_matrix_set_colnames(P, cnames);
	pprintf(prn, "%s:\n\n", _(label));
	gretl_matrix_print_with_format(P, fmt, 0, 0, prn);
	pputc(prn, '\n');
    } else {
	pprintf(prn, "%s: 0\n\n", _(label));
	strings_array_free(cnames, np);
    }

    gretl_matrix_free(P);

    return 0;
}

static int BKW_analyse (gretl_matrix *B, double maxcond,
			const char *fmt, PRN *prn)
{
    const char *labels[] = {
	N_("Count of condition indices >= 30"),
	N_("Variance proportions >= 0.5 associated with cond >= 30"),
	N_("Count of condition indices >= 10"),
        N_("Variance proportions >= 0.5 associated with cond >= 10"),
	N_("No evidence of excessive collinearity")
    };
    const char **bnames = NULL;
    int rows = B->rows;
    int ngc10 = 0, ngc30 = 0;
    int i, err = 0;

    pputs(prn, _("According to BKW, cond >= 30 indicates \"strong\" near linear dependence,\n"
		 "and cond between 10 and 30 \"moderately strong\".  Parameter estimates whose\n"
		 "variance is mostly associated with problematic cond values may themselves\n"
		 "be considered problematic."));
    pputs(prn, "\n\n");

    /* count rows with condition index >= thresholds */
    for (i=rows-1; i>=0; i--) {
	if (gretl_matrix_get(B, i, 1) >= 30) {
	    ngc30++;
	}
	if (gretl_matrix_get(B, i, 1) >= 10) {
	    ngc10++;
	} else {
	    break;
	}
    }

    if (ngc10 > 0) {
	bnames = gretl_matrix_get_colnames(B);
    }

    pprintf(prn, "%s: %d\n", _(labels[0]), ngc30);
    if (ngc30 == 0 && ngc10 > 0) {
	pputc(prn, '\n');
    }

    if (ngc30 > 0) {
	/* variance proportion sums, cond >= 30 */
	do_proportion_sums(B, bnames, _(labels[1]), 30, prn);
    }

    pprintf(prn, "%s: %d\n", _(labels[2]), ngc10);

    if (ngc10 > ngc30) {
	/* variance proportion sums, cond >= 10 */
	do_proportion_sums(B, bnames, _(labels[3]), 10, prn);
    } else if (ngc10 == 0) {
	pprintf(prn, "\n%s\n", _(labels[4]));
    }

    if (!err) {
	pputc(prn, '\n');
    }

    return err;
}

static void BKW_print (gretl_matrix *B, int namelen, PRN *prn)
{
    const char *strs[] = {
	N_("Belsley-Kuh-Welsch collinearity diagnostics"),
	N_("variance proportions"),
	N_("eigenvalues of inverse covariance matrix"),
	N_("condition index"),
	N_("note: variance proportions columns sum to 1.0")
    };
    double maxcond = gretl_matrix_get(B, B->rows - 1, 1);
    int sp = maxcond >= 10000 ? 10 : maxcond >= 1000 ? 9 : 8;
    int maxcols = 9;
    char fmt[16];

    if (sp < namelen + 1) {
	sp = namelen + 1;
    }

    sprintf(fmt, "%%%d.3f", sp);
    pprintf(prn, "\n%s:\n\n", _(strs[0]));
    bufspace(2, prn);
    pprintf(prn, "%s\n\n", _(strs[1]));

    if (B->cols > maxcols) {
	int err = 0;
	gretl_array *a = gretl_matrix_col_split(B, 2, maxcols, &err);

	if (a != NULL) {
	    int i, n = gretl_array_get_length(a);
	    gretl_matrix *ai;

	    for (i=0; i<n; i++) {
		ai = gretl_array_get_data(a, i);
		gretl_matrix_print_with_format(ai, fmt, 0, 0, prn);
		if (i < n-1) {
		    pputc(prn, '\n');
		}
	    }
	    gretl_array_destroy(a);
	}
    } else {
	gretl_matrix_print_with_format(B, fmt, 0, 0, prn);
    }

    pprintf(prn, "\n  lambda = %s (smallest is %g)\n", _(strs[2]),
	    gretl_matrix_get(B, B->rows - 1, 0));
    pprintf(prn, "  cond   = %s\n", _(strs[3]));
    pprintf(prn, "  %s\n\n", _(strs[4]));

    BKW_analyse(B, maxcond, fmt, prn);
}

static gretl_array *BKW_pnames (MODEL *pmod, DATASET *dset)
{
    gretl_array *pnames;
    char pname[VNAMELEN];
    int i, k = pmod->ncoeff;
    int err = 0;

    pnames = gretl_array_new(GRETL_TYPE_STRINGS, k, &err);

    if (pnames != NULL) {
	for (i=0; i<pmod->ncoeff; i++) {
	    gretl_model_get_param_name(pmod, dset, i, pname);
	    gretl_array_set_string(pnames, i, pname, 1);
	}
    }

    return pnames;
}

int compute_vifs (MODEL *pmod, DATASET *dset,
		  gretlopt opt, PRN *prn)
{
    gretl_vector *vif = NULL;
    int *xlist;
    int quiet = (opt & OPT_Q);
    int i, err = 0;

    /* fetch list of regressors */
    xlist = gretl_model_get_x_list(pmod);
    if (xlist == NULL) {
	return E_DATA;
    }

    /* drop the constant if present in xlist */
    for (i=xlist[0]; i>0; i--) {
	if (xlist[i] == 0) {
	    gretl_list_delete_at_pos(xlist, i);
	    break;
	}
    }

    if (xlist[0] > 1) {
	vif = model_vif_vector(pmod, xlist, dset, &err);
    }

    if (vif != NULL && !quiet) {
	int vlen = gretl_vector_get_length(vif);
	int vi, n, maxlen = 0;
	double vj;

	pprintf(prn, "\n%s\n", _("Variance Inflation Factors"));
	pprintf(prn, "%s\n", _("Minimum possible value = 1.0"));
	pprintf(prn, "%s\n", _("Values > 10.0 may indicate a collinearity problem"));
	pputc(prn, '\n');

	for (i=0; i<vlen; i++) {
	    vi = xlist[i+1];
	    vj = vif->val[i];
	    if (!na(vj)) {
		n = strlen(dset->varname[vi]);
		if (n > maxlen) {
		    maxlen = n;
		}
	    }
	}

	maxlen = maxlen < 12 ? 12 : maxlen;

	for (i=0; i<vlen; i++) {
	    vi = xlist[i+1];
	    vj = vif->val[i];
	    if (!quiet && !na(vj)) {
		pprintf(prn, "%*s %8.3f\n", maxlen, dset->varname[vi], vj);
	    }
	}

	pputc(prn, '\n');
	pputs(prn, _("VIF(j) = 1/(1 - R(j)^2), where R(j) is the "
		     "multiple correlation coefficient\nbetween "
		     "variable j and the other independent variables"));
	pputc(prn, '\n');
    }

    if (!err && !(opt & OPT_G)) {
	set_last_result_data(vif, GRETL_TYPE_MATRIX);
    } else {
	gretl_matrix_free(vif);
    }

    free(xlist);

    return err;
}

int compute_bkw (MODEL *pmod, DATASET *dset,
		 gretlopt opt, PRN *prn)
{
    gretl_matrix *V, *BKW = NULL;
    int quiet = (opt & OPT_Q);
    int err = 0;

    V = gretl_vcv_matrix_from_model(pmod, NULL, &err);

    if (!err) {
	gretl_array *pnames = BKW_pnames(pmod, dset);
	PRN *vprn = quiet ? NULL : prn;

	BKW = bkw_matrix(V, pnames, vprn, &err);
	gretl_array_destroy(pnames);
	gretl_matrix_free(V);
    }

    if (!err && !(opt & OPT_G)) {
	set_last_result_data(BKW, GRETL_TYPE_MATRIX);
    } else {
	gretl_matrix_free(BKW);
    }

    return err;
}
