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

static void BKW_print (gretl_matrix *B, PRN *prn);

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

static void bkw_add_colnames (gretl_matrix *BKW,
			      gretl_array *pnames)
{
    char **S = strings_array_new(BKW->cols);

    if (S != NULL) {
	int i, k = BKW->cols - 2;
	char tmp[16];

	S[0] = gretl_strdup("lambda");
	S[1] = gretl_strdup("cond");
	for (i=0; i<k; i++) {
	    if (pnames != NULL) {
		S[i+2] = gretl_array_get_data(pnames, i);
		gretl_array_set_data(pnames, i, NULL);
	    } else {
		sprintf(tmp, "x%d", i+1);
		S[i+2] = gretl_strdup(tmp);
	    }
	}
	gretl_matrix_set_colnames(BKW, S);
    }
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
    int i, j;

    if (pnames != NULL && gretl_array_get_length(pnames) != k) {
	fprintf(stderr, "bkw_matrix: expected %d names but got %d\n",
		k, gretl_array_get_length(pnames));
	*err = E_INVARG;
	return NULL;
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
	*err = gretl_matrix_SVD(Q, NULL, &lambda, &V);
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

    bkw_add_colnames(BKW, pnames);

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
	BKW_print(BKW, prn);
    }

    return BKW;
}

static int BKW_analyse (gretl_matrix *B, double maxcond,
			const char *fmt, PRN *prn)
{
    const char *labels[] = {
        N_("Potentially problematic rows, with condition index >= 10"),
	N_("Sums of variance proportions attributable to the affected rows"),
	N_("Problematic parameter estimates (proportion >= 0.5)"),
	N_("No evidence of excessive collinearity")
    };
    const char **bnames = NULL;
    char **cnames = NULL;
    gretl_matrix *P = NULL;
    gretl_matrix *S = NULL;
    double x;
    int rows = B->rows;
    int np = B->cols - 2;
    int nbigc = 0, ngp5 = 0;
    int bigc = 10;
    int i, j, k;
    int err = 0;

    /* count rows with condition index >= bigc */
    for (i=rows-1; i>=0; i--) {
	if (gretl_matrix_get(B, i, 1) >= bigc) {
	    nbigc++;
	} else {
	    break;
	}
    }

    if (nbigc > 0) {
	/* construct matrix showing problematic rows */
	P = gretl_matrix_alloc(nbigc, B->cols - 1);
	if (P == NULL) {
	    err = E_ALLOC;
	} else {
	    for (j=1; j<B->cols; j++) {
		k = 0;
		for (i=rows-nbigc; i<rows; i++) {
		    x = gretl_matrix_get(B, i, j);
		    gretl_matrix_set(P, k++, j-1, x);
		}
	    }
	    bnames = gretl_matrix_get_colnames(B);
	    cnames = strings_array_dup((char **) bnames + 1, B->cols - 1);
	    gretl_matrix_set_colnames(P, cnames);
	    pprintf(prn, "%s:\n\n", _(labels[0]));
	    gretl_matrix_print_with_format(P, fmt, 0, 0, prn);
	}
    }

    if (nbigc > 0 && !err) {
	/* row vector showing sums of variance proportions */
	S = gretl_matrix_alloc(1, np);
	if (S == NULL) {
	    err = E_ALLOC;
	} else {
	    for (j=1; j<P->cols; j++) {
		x = 0;
		for (i=0; i<P->rows; i++) {
		    x += gretl_matrix_get(P, i, j);
		}
		S->val[j-1] = x;
		if (x >= 0.5) {
		    ngp5++;
		}
	    }
	    cnames = strings_array_dup((char **) bnames + 2, B->cols - 2);
	    gretl_matrix_set_colnames(S, cnames);
	    pprintf(prn, "\n%s:\n\n", _(labels[1]));
	    gretl_matrix_print_with_format(S, "%8.3f", 0, 0, prn);
	}
    }

    if (ngp5 > 1 && !err) {
	/* row vector showing problematic estimates */
	gretl_matrix_free(P);
	P = gretl_matrix_alloc(1, ngp5);
	if (P == NULL) {
	    err = E_ALLOC;
	} else {
	    cnames = strings_array_new(ngp5);
	    k = 0;
	    for (j=0; j<S->cols; j++) {
		if (S->val[j] >= 0.5) {
		    P->val[k] = S->val[j];
		    cnames[k] = gretl_strdup(bnames[j+2]);
		    k++;
		}
	    }
	    gretl_matrix_set_colnames(P, cnames);
	    pprintf(prn, "\n%s:\n\n", _(labels[2]));
	    gretl_matrix_print_with_format(P, "%8.3f", 0, 0, prn);
	}
    }

    if (!err) {
	if (ngp5 < 1) {
	    pprintf(prn, "%s\n", _(labels[3]));
	} else {
	    pputc(prn, '\n');
	}
    }

    gretl_matrix_free(P);
    gretl_matrix_free(S);

    return err;
}

static void BKW_print (gretl_matrix *B, PRN *prn)
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
    char fmt[8];

    sprintf(fmt, "%%%d.3f", sp);
    pprintf(prn, "\n%s:\n\n", _(strs[0]));
    bufspace(25, prn);
    pprintf(prn, "--- %s ---\n", _(strs[1]));
    gretl_matrix_print_with_format(B, fmt, 0, 0, prn);
    pprintf(prn, "\n  lambda = %s\n", _(strs[2]));
    pprintf(prn, "  cond   = %s\n", _(strs[3]));
    pprintf(prn, "  %s\n\n", _(strs[4]));

    BKW_analyse(B, maxcond, fmt, prn);
}

static void maybe_truncate_param_name (char *s)
{
    int n = strlen(s);

    if (n > 9) {
	char tmp[VNAMELEN];

	tmp[0] = '\0';
	strncat(tmp, s, 8);
	strcat(tmp, "~");
	strcpy(s, tmp);
    }
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
	    maybe_truncate_param_name(pname);
	    gretl_array_set_string(pnames, i, pname, 1);
	}
    }

    return pnames;
}

int compute_vifs (MODEL *pmod, DATASET *dset,
		  gretlopt opt, PRN *prn)
{
    gretl_matrix *BKW = NULL;
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
	if (err) {
	    return err;
	}
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

    if (1) {
	/* this should get hived off to a separate function,
	   but for now...
	*/
	gretl_matrix *V;

	V = gretl_vcv_matrix_from_model(pmod, NULL, &err);

	if (!err) {
	    gretl_array *pnames = BKW_pnames(pmod, dset);
	    PRN *vprn = quiet ? NULL : prn;

	    BKW = bkw_matrix(V, pnames, vprn, &err);
	    gretl_array_destroy(pnames);
	}
	gretl_matrix_free(V);
    }

    if (!err) {
	gretl_bundle *b = gretl_bundle_new();

	if (b != NULL && (vif != NULL || BKW != NULL)) {
	    if (vif != NULL) {
		gretl_bundle_donate_data(b, "vif", vif, GRETL_TYPE_MATRIX, 0);
		vif = NULL;
	    }
	    if (BKW != NULL) {
		gretl_bundle_donate_data(b, "BKW", BKW, GRETL_TYPE_MATRIX, 0);
		BKW = NULL;
	    }
	    set_last_result_data(b, GRETL_TYPE_BUNDLE);
	}
    }

    gretl_matrix_free(vif);
    gretl_matrix_free(BKW);
    free(xlist);

    return 0;
}

int gui_bkw (MODEL *pmod, const DATASET *dset, PRN *prn)
{
    gretl_matrix *V = NULL, *B = NULL;
    gretl_array *pnames = NULL;
    int err = 0;

    V = gretl_vcv_matrix_from_model(pmod, NULL, &err);

    if (!err) {
	pnames = gretl_model_get_param_names(pmod, dset, &err);
    }

    if (!err) {
	B = bkw_matrix(V, pnames, prn, &err);
    }

    gretl_matrix_free(V);
    gretl_matrix_free(B);
    gretl_array_destroy(pnames);

    return err;
}
