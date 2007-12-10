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

/* tsls.c - Two-Stage Least Squares for gretl */

#include "libgretl.h"
#include "qr_estimate.h"
#include "libset.h"
#include "estim_private.h"
#include "matrix_extra.h"

#define TDEBUG 0

static void tsls_omitzero (int *list, const double **Z, int t1, int t2)
{
    int i, v;

    for (i=2; i<=list[0]; i++) {
        v = list[i];
        if (gretl_iszero(t1, t2, Z[v])) {
	    gretl_list_delete_at_pos(list, i--);
	}
    }
}

void tsls_free_data (const MODEL *pmod)
{
    const char *endog = gretl_model_get_data(pmod, "endog");
    double **X = gretl_model_get_data(pmod, "tslsX");
    int i, m = 0;

    if (endog != NULL && X != NULL) {
	for (i=0; i<pmod->ncoeff; i++) {
	    if (endog[i]) m++;
	}
	for (i=0; i<m; i++) {
	    free(X[i]);
	}
    }
}

static int 
tsls_save_data (MODEL *pmod, const int *hatlist, const int *exolist,
		double **Z, DATAINFO *pdinfo)
{
    double **X = NULL;
    char *endog = NULL;
    int i, j, v, err = 0;
    size_t Xsize, esize;
    size_t xs_old, es_old;
    int recycle_X = 0;
    int recycle_e = 0;

    esize = pmod->ncoeff * sizeof *endog;
    Xsize = hatlist[0] * sizeof *X;

    /* re-use old pointers if applicable */

    X = gretl_model_get_data_and_size(pmod, "tslsX", &xs_old);
    if (X != NULL) {
	if (Xsize == xs_old) {
	    recycle_X = 1;
	} else {
	    tsls_free_data(pmod);
	    gretl_model_detach_data_item(pmod, "tslsX");
	    free(X);
	    X = NULL;
	} 
    }

    if (!recycle_X && Xsize > 0) {
	X = malloc(Xsize);
	if (X == NULL) {
	    return E_ALLOC;
	}
    }

    endog = gretl_model_get_data_and_size(pmod, "endog", &es_old);
    if (endog != NULL) {
	if (esize == es_old) {
	    recycle_e = 1;
	} else {
	    gretl_model_detach_data_item(pmod, "endog");
	    free(endog);
	    endog = NULL;
	}
    }

    if (!recycle_e && esize > 0) {
	endog = malloc(esize);
	if (endog == NULL) {
	    free(X);
	    return E_ALLOC;
	}
    }

    /* steal the appropriate columns from Z */
    for (i=1; i<=hatlist[0]; i++) {
	v = hatlist[i];
	X[i-1] = Z[v];
	Z[v] = NULL;
    }

    for (i=0; i<pmod->ncoeff; i++) {
	v = pmod->list[i+2];
	endog[i] = 1;
	for (j=1; j<=exolist[0]; j++) {
	    if (exolist[j] == v) {
		endog[i] = 0;
		break;
	    }
	}
    }

    /* now attach X and endog to the model */
    if (!recycle_X) {
	gretl_model_set_data(pmod, "tslsX", X, GRETL_TYPE_DOUBLE_ARRAY,
			     Xsize);
    }
    if (!recycle_e) {
	gretl_model_set_data(pmod, "endog", endog, GRETL_TYPE_CHAR_ARRAY,
			     esize);
    }

    return err;
}

double *tsls_get_Xi (const MODEL *pmod, double **Z, int i)
{
    const char *endog = gretl_model_get_data(pmod, "endog");
    double **X;
    double *ret;

    if (endog == NULL) {
	return NULL;
    }

    X = gretl_model_get_data(pmod, "tslsX");

    if (!endog[i]) {
	ret = Z[pmod->list[i+2]];
    } else if (X == NULL) {
	ret = NULL;
    } else {
	int j, k = 0;

	for (j=0; j<i; j++) {
	    if (endog[j]) k++;
	}
	ret = X[k];
    }

    return ret;
}

static void add_tsls_var (int *list, int v, gretlopt opt)
{
    int pos = gretl_list_separator_position(list);
    int i;

    if (opt & OPT_T) {
	/* add as instrument */
	list[0] += 1;
	list[list[0]] = v;
    } else if (opt & OPT_B) {
	/* add as exogenous regressor */
	list[0] += 2;
	for (i=list[0]-1; i>pos; i--) {
	    list[i] = list[i-1];
	}
	list[pos] = v;
	list[list[0]] = v;
    } else {
	/* add as endogenous regressor */
	list[0] += 1;
	for (i=list[0]; i>pos; i--) {
	    list[i] = list[i-1];
	}
	list[pos] = v;
    }
}

static void delete_tsls_var (int *list, int v, gretlopt opt)
{
    int pos = gretl_list_separator_position(list);
    int i, j;

    if (opt & OPT_T) {
	/* delete as instrument */
	for (i=pos+1; i<=list[0]; i++) {
	    if (list[i] == v) {
		for (j=i; j<list[0]; j++) {
		    list[j] = list[j + 1]; 
		}
		list[0] -= 1;
		break;
	    }
	}
    } else if (opt & OPT_B) {	
	/* delete from both sub-lists */
	for (i=2; i<=list[0]; i++) {
	    if (list[i] == v) {
		for (j=i; j<list[0]; j++) {
		    list[j] = list[j + 1]; 
		}
		list[0] -= 1;
	    }
	}
    } else {
	/* delete from regressor list */
	for (i=2; i<pos; i++) {
	    if (list[i] == v) {
		for (j=i; j<list[0]; j++) {
		    list[j] = list[j + 1]; 
		}
		list[0] -= 1;
		break;
	    }
	}
    }
}

static int in_tsls_list (const int *list, int v, gretlopt opt)
{
    int pos = gretl_list_separator_position(list);
    int i, imin, imax;

    if (pos == 0) {
	return 0;
    }

    if (opt & OPT_T) {
	imin = pos + 1;
	imax = list[0];
    } else if (opt & OPT_B) {
	imin = 2;
	imax = list[0];
    } else {
	imin = 2;
	imax = pos - 1;
    }

    for (i=imin; i<=imax; i++) {
	if (list[i] == v) {
	    return i;
	}
    }

    return 0;
}

/**
 * tsls_list_omit:
 * @orig: list specifying original TSLS model.
 * @drop: list of variables to be omitted.
 * @opt: may include %OPT_T (omit variable only as instrument)
 * or %OPT_B (omit both as regressor and instrument).  The default
 * is to omit (only) as regressor.
 * @err: pointer to receive error code.
 *
 * Creates a new TSLS specification list, by first copying @orig
 * then deleting from the copy the variables found in @drop.
 *
 * Returns: the new, reduced list or %NULL on error.
 */

int *
tsls_list_omit (const int *orig, const int *drop, gretlopt opt, int *err)
{
    int *newlist;
    int i;

    if ((opt & OPT_T) && (opt & OPT_B)) {
	/* incoherent options */
	*err = E_BADOPT;
	return NULL;
    }

    newlist = gretl_list_copy(orig);

    for (i=1; i<=drop[0]; i++) {
	if (in_tsls_list(orig, drop[i], opt)) {
	    delete_tsls_var(newlist, drop[i], opt);
	} else {
	    *err = E_UNSPEC;
	}
    }

    if (*err) {
	free(newlist);
	newlist = NULL;
    }

    return newlist;
}

/**
 * tsls_list_add:
 * @orig: list specifying original TSLS model.
 * @add: list of variables to be added.
 * @opt: may include %OPT_T (add variable only as instrument)
 * or %OPT_B (add both as regressor and instrument).  The default
 * is to add as an endogenous regressor.
 * @err: pointer to receive error code.
 *
 * Creates a new TSLS specification list, by copying @orig
 * and adding the variables found in @add.
 *
 * Returns: the new, augmented list or %NULL on error.
 */

int *
tsls_list_add (const int *orig, const int *add, gretlopt opt, int *err)
{
    int norig = orig[0];
    int nadd = add[0];
    int *newlist;
    int i;

    if ((opt & OPT_T) && (opt & OPT_B)) {
	/* incoherent options */
	*err = E_BADOPT;
	return NULL;
    }

    if (opt & OPT_B) {
	nadd *= 2;
    }

    newlist = gretl_list_new(norig + nadd);

    for (i=0; i<=norig; i++) {
	newlist[i] = orig[i];
    }

    for (i=1; i<=add[0]; i++) {
	if (in_tsls_list(orig, add[i], opt)) {
	    *err = E_ADDDUP;
	} else {
	    add_tsls_var(newlist, add[i], opt);
	} 
    }

    if (*err) {
	free(newlist);
	newlist = NULL;
    }

    return newlist;
}

/*
  tsls_make_hatlist: determines which variables in reglist, when
  checked against the predetermined and exogenous vars in instlist,
  need to be instrumented, and populates hatlist accordingly
*/

static int 
tsls_make_hatlist (const int *reglist, int *instlist, int *hatlist)
{
    int i, j, k = 0;
    int endog, addconst = 0;

    for (i=2; i<=reglist[0]; i++) {  
	endog = 1;
	for (j=1; j<=instlist[0]; j++) {
	    if (instlist[j] == reglist[i]) {
		endog = 0;
		break;
	    }
	}
	if (endog) {
	    if (reglist[i] == 0) {
		/* found const in reglist but not instlist: needs fixing */
		addconst = 1;
	    } else {
		hatlist[++k] = reglist[i];
	    }
	} 
    }

    hatlist[0] = k;

    if (addconst) {
	/* add constant to list of instruments */
	instlist[0] += 1;
	for (i=instlist[0]; i>=2; i--) {
	    instlist[i] = instlist[i - 1];
	}
	instlist[1] = 0;
    }

    return addconst;
}

/* perform the Sargan overidentification test for a model
   estimated via TSLS */

static int 
tsls_sargan_test (MODEL *tsls_model, int Orank, int *instlist, 
		  double ***pZ, DATAINFO *pdinfo)
{
    int t1 = tsls_model->t1;
    int t2 = tsls_model->t2;
    int ninst = instlist[0];
    int i, t, err = 0;

    MODEL smod;
    int *OT_list = NULL;
    int nv = pdinfo->v;

    if (Orank == 0) {
	return 0;
    }

    if (tsls_model->list[0] == 2 && tsls_model->list[2] == 0) {
	/* degenerate model with const only */
	return 0;
    }

    err = dataset_add_series(1, pZ, pdinfo);
    if (err) {
	return err;
    }

    for (t=t1; t<=t2; t++) {
	(*pZ)[nv][t] = tsls_model->uhat[t];
    }

    OT_list = gretl_list_new(ninst + 1);
    if (OT_list == NULL) {
	dataset_drop_last_variables(1, pZ, pdinfo);
	return E_ALLOC;
    }

    OT_list[1] = nv;
    for (i=2; i<=OT_list[0]; i++) {
	OT_list[i] = instlist[i-1];
    }

    smod = lsq(OT_list, pZ, pdinfo, OLS, OPT_A);
    if (smod.errcode) {
	fprintf(stderr, "tsls_sargan_test: smod.errcode = %d\n", smod.errcode);
	err = smod.errcode;
    } else {
	ModelTest *test = model_test_new(GRETL_TEST_SARGAN);
	double OTest = smod.rsq * smod.nobs;

	if (test != NULL) {
	    model_test_set_teststat(test, GRETL_STAT_TR2);
	    model_test_set_dfn(test, Orank);
	    model_test_set_value(test, OTest);
	    model_test_set_pvalue(test, chisq_cdf_comp(OTest, Orank));
	    maybe_add_test_to_model(tsls_model, test);
	}	
    }

    clear_model(&smod);
    dataset_drop_last_variables(1, pZ, pdinfo);
    free(OT_list);

    return err;
}

static int 
tsls_hausman_test (MODEL *tmod, int *reglist, int *hatlist,
		   double ***pZ, DATAINFO *pdinfo)
{
    double RRSS;
    int df, err = 0;

    MODEL hmod;
    int *HT_list = NULL;

    hmod = lsq(reglist, pZ, pdinfo, OLS, OPT_A);
    if (hmod.errcode) {
	fprintf(stderr, "tsls_hausman_test: hmod.errcode (R) = %d\n", hmod.errcode);
	err = hmod.errcode;
	goto bailout;
    } 

    RRSS = hmod.ess;
    clear_model(&hmod);

    HT_list = gretl_list_add(reglist, hatlist, &err);
    if (err) {
	fprintf(stderr, "tsls_hausman_test: gretl_list_add: err = %d\n", err);
	goto bailout;
    } 

    hmod = lsq(HT_list, pZ, pdinfo, OLS, OPT_A);

    if (hmod.errcode) {
	fprintf(stderr, "tsls_hausman_test: hmod.errcode (U) = %d\n", hmod.errcode);
	err = hmod.errcode;
	goto bailout;
    } else if ((df = hmod.list[0] - reglist[0]) > 0) {
	/* Degrees of freedom for the Hausman test need not be equal to the
	   number of added regressors, since some of them may not really be 
	   endogenous.
	*/
	double URSS = hmod.ess;
	double HTest = (RRSS/URSS - 1) * hmod.nobs;
	ModelTest *test = model_test_new(GRETL_TEST_TSLS_HAUSMAN);

	if (test != NULL) {
	    model_test_set_teststat(test, GRETL_STAT_WALD_CHISQ);
	    model_test_set_dfn(test, df);
	    model_test_set_value(test, HTest);
	    model_test_set_pvalue(test, chisq_cdf_comp(HTest, df));
	    maybe_add_test_to_model(tmod, test);
	}
    }

 bailout:

    clear_model(&hmod);
    free(HT_list);

    return err;
}

/* form matrix of instruments and perform QR decomposition
   of this matrix; return Q */

static gretl_matrix *tsls_Q (int *instlist, int *reglist, int **pdlist,
			     const double **Z, int t1, int t2, char *mask,
			     int *err)
{
    gretl_matrix *Q = NULL;
    gretl_matrix *R = NULL;
    int rank, ndrop = 0;
    int *droplist = NULL;
    double test;
    int i, j, k;

    Q = gretl_matrix_data_subset(instlist, Z, t1, t2, mask, err);
    if (*err) {
	return NULL;
    }

    k = gretl_matrix_cols(Q);

    R = gretl_matrix_alloc(k, k);
    if (R == NULL) {
	*err = E_ALLOC;
	goto bailout;
    }    

    *err = gretl_matrix_QR_decomp(Q, R);
    if (*err) {
	goto bailout;
    }

    rank = gretl_check_QR_rank(R, err);

#if TDEBUG
    fprintf(stderr, "tsls_Q: k = Q->cols = %d, rank = %d\n", k, rank);
    if (rank == 0) {
	gretl_matrix_print(R, "R");
    }
#endif

    if (*err) {
	goto bailout;
    } else if (rank < k) {
	fprintf(stderr, "tsls_Q: k = %d, rank = %d\n", k, rank);
	ndrop = k - rank;
    }

    if (ndrop > 0) {
	droplist = gretl_list_new(ndrop);
	if (droplist != NULL) {
	    droplist[0] = 0;
	}

	j = 1;
	for (i=0; i<k; i++) {
	    test = gretl_matrix_get(R, i, i);
	    if (fabs(test) < R_DIAG_MIN) {
		if (droplist != NULL) {
		    droplist[0] += 1;
		    droplist[droplist[0]] = instlist[j];
		}	    
		fprintf(stderr, "tsls_Q: dropping redundant instrument %d\n", 
			instlist[j]);
		gretl_list_delete_at_pos(instlist, j--);
	    }
	    j++;
	}

	k = instlist[0];
	gretl_matrix_free(Q);
	Q = gretl_matrix_data_subset(instlist, Z, t1, t2, mask, err);
	if (!*err) {
	    R = gretl_matrix_reuse(R, k, k);
	    *err = gretl_matrix_QR_decomp(Q, R);
	}
    }

 bailout:

    gretl_matrix_free(R);

    if (*err) {
	free(droplist);
	gretl_matrix_free(Q);
	Q = NULL;
    } else {
	*pdlist = droplist;
    }

    return Q;
}

static int tsls_form_xhat (gretl_matrix *Q, double *x, double *xhat,
			   DATAINFO *pdinfo, const char *mask)
{
    double *r;
    int k = gretl_matrix_cols(Q);
    int i, t, s = 0;

    r = malloc(k * sizeof *r);
    if (r == NULL) {
	return E_ALLOC;
    }

#if TDEBUG
    fprintf(stderr, "tsls_form_xhat: t1=%d, t2=%d, mask=%p\n",
	    pdinfo->t1, pdinfo->t2, (void *) mask);
#endif

    /* form r = Q'y */
    for (i=0; i<k; i++) {
	r[i] = 0.0;
	s = 0;
	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    if (mask != NULL && mask[t - pdinfo->t1]) {
		continue;
	    }
	    r[i] += gretl_matrix_get(Q, s++, i) * x[t];
	}
    }

    /* form Qr = QQ'y */
    s = 0;
    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (mask != NULL && mask[t - pdinfo->t1]) {
	    xhat[t] = NADBL;
	    continue;
	}
	xhat[t] = 0.0;
	for (i=0; i<k; i++) {
	    xhat[t] += gretl_matrix_get(Q, s, i) * r[i];
	}
	s++;
    }

    free(r);

    return 0;
}

static void tsls_residuals (MODEL *pmod, const int *reglist,
			    const double **Z, gretlopt opt)
{
    int den = (opt & OPT_N)? pmod->nobs : pmod->dfd;
    int yno = pmod->list[1];
    double sigma0 = pmod->sigma;
    double yh;
    int i, t;

    pmod->ess = 0.0;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (model_missing(pmod, t)) {
	    continue;
	}
	yh = 0.0;
	for (i=0; i<pmod->ncoeff; i++) {
	    yh += pmod->coeff[i] * Z[reglist[i+2]][t];
	}
	pmod->yhat[t] = yh; 
	pmod->uhat[t] = Z[yno][t] - yh;
	pmod->ess += pmod->uhat[t] * pmod->uhat[t];
    }

    if (pmod->ess <= 0.0) {
	pmod->sigma = 0.0;
    } else {
	pmod->sigma = sqrt(pmod->ess / den);
    }

    if (sigma0 > 0.0) {
	double corr = pmod->sigma / sigma0;

	for (i=0; i<pmod->ncoeff; i++) {
	    pmod->sderr[i] *= corr;
	}
    }
}

static void tsls_extra_stats (MODEL *pmod, const double **Z)
{
    int yno = pmod->list[1];
    double r;

    pmod->rsq = gretl_corr_rsq(pmod->t1, pmod->t2, Z[yno], pmod->yhat);
    r = 1.0 - pmod->rsq;

    pmod->adjrsq = 1.0 - (r * (pmod->nobs - 1.0) / pmod->dfd);

    pmod->fstt = pmod->rsq * pmod->dfd / (pmod->dfn * r);

    ls_criteria(pmod);

    if (pmod->missmask == NULL) {
	/* no missing obs within sample range */
	pmod->rho = rhohat(1, pmod->t1, pmod->t2, pmod->uhat);
	pmod->dw = dwstat(1, pmod, Z);
    } else {
	pmod->rho = pmod->dw = NADBL;
    }
}

static void tsls_recreate_full_list (MODEL *pmod, const int *reglist,
				     const int *instlist)
{
    int len = reglist[0] + instlist[0] + 1;
    int *full_list;

    full_list = gretl_list_new(len);

    if (full_list == NULL) {
	pmod->errcode = E_ALLOC;
    } else {
	int i, j = 1;

	for (i=1; i<=reglist[0]; i++) {
	    full_list[j++] = reglist[i];
	}

	full_list[j++] = LISTSEP;

	for (i=1; i<=instlist[0]; i++) {
	    full_list[j++] = instlist[i];
	}

	free(pmod->list);
	pmod->list = full_list;
    }
}

static void replace_list_element (int *list, int targ, int repl)
{
    int i;

    for (i=2; i<=list[0]; i++) {
	if (list[i] == targ) {
	    list[i] = repl;
	    break;
	}
    }
}

static int 
reglist_remove_redundant_vars (const MODEL *tmod, int *s2list, int *reglist)
{
    int *dlist = gretl_model_get_data(tmod, "droplist");
    int i, pos;

    if (dlist == NULL) {
	return 1;
    }

    for (i=1; i<=dlist[0]; i++) {
	pos = in_gretl_list(s2list, dlist[i]);
	if (pos > 1) {
	    gretl_list_delete_at_pos(reglist, pos);
	}
    }

    return 0;
}

static int 
compute_first_stage_F (MODEL *pmod, int v, const int *reglist, 
		       const int *instlist, double ***pZ, 
		       DATAINFO *pdinfo, gretlopt opt)
{
    MODEL mod1;
    gretlopt myopt;
    int *list1 = NULL;
    int *flist = NULL;
    double F;
    int i, vi;
    int err = 0;

    gretl_model_init(&mod1);

    list1 = gretl_list_append_term(&list1, v);
    if (list1 == NULL) {
	return E_ALLOC;
    }

    for (i=1; i<=instlist[0]; i++) {
	vi = instlist[i];
	list1 = gretl_list_append_term(&list1, vi);
	if (list1 == NULL) {
	    err = E_ALLOC;
	    break;
	}
	if (!in_gretl_list(reglist, vi)) {
	    flist = gretl_list_append_term(&flist, vi);
	    if (flist == NULL) {
		err = E_ALLOC;
		break;
	    }
	}
    }

    if (!err) {
	myopt = (opt & OPT_R)? (OPT_A | OPT_R) : OPT_A;
	mod1 = lsq(list1, pZ, pdinfo, OLS, myopt);
	err = mod1.errcode;
    }

    if (!err) {
	F = robust_omit_F(flist, &mod1);
	if (na(F)) {
	    err = 1;
	}
    }

    if (!err) {
	gretl_model_set_double(pmod, "stage1-F", F);
	gretl_model_set_int(pmod, "stage1-dfn", flist[0]);
	gretl_model_set_int(pmod, "stage1-dfd", mod1.dfd);
    }
    
    clear_model(&mod1);
    free(list1);
    free(flist);

    return err;
}

static int 
tsls_adjust_sample (const int *list, int *t1, int *t2, 
		    const double **Z, char **pmask)
{
    int i, t, t1min = *t1, t2max = *t2;
    char *mask = NULL;
    int T, vi, missobs;
    int err = 0;

    /* advance start of sample range to skip missing obs? */
    for (t=t1min; t<t2max; t++) {
	missobs = 0;
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    if (vi == 0 || vi == LISTSEP) {
		continue;
	    }
	    if (na(Z[vi][t])) {
		missobs = 1;
		break;
	    }
	}
	if (missobs) {
	    t1min++;
	} else {
	    break;
	}
    }

    /* retard end of sample range to skip missing obs? */
    for (t=t2max; t>t1min; t--) {
	missobs = 0;
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    if (vi == 0 || vi == LISTSEP) {
		continue;
	    }
	    if (na(Z[vi][t])) {
		missobs = 1;
		break;
	    }
	}
	if (missobs) {
	    t2max--;
	} else {
	    break;
	}	
    }

    /* count missing values in mid-range of data */
    missobs = 0;
    for (t=t1min; t<=t2max; t++) {
	for (i=1; i<=list[0]; i++) {
	    vi = list[i];
	    if (vi == 0 || vi == LISTSEP) {
		continue;
	    }
	    if (na(Z[vi][t])) {
		missobs++;
		break;
	    }
	}
    }
    
    T = t2max - t1min + 1;

    if (missobs == T) {
	err = E_MISSDATA;
    } else if (missobs > 0) {
	mask = calloc(T, 1); /* all NUL bytes */
	if (mask == NULL) {
	    err = E_ALLOC;
	} else {
	    for (t=t1min; t<=t2max; t++) {
		for (i=1; i<=list[0]; i++) {
		    vi = list[i];
		    if (vi == 0 || vi == LISTSEP) {
			continue;
		    }
		    if (na(Z[vi][t])) {
			mask[t - t1min] = 1;
			break;
		    }
		}
	    }
	}
    }

#if TDEBUG
    fprintf(stderr, "tsls_adjust_sample: t1=%d, t2=%d, missobs=%d, ok obs=%d\n", 
	    t1min, t2max, missobs, t2max - t1min + 1 - missobs);
#endif

    *t1 = t1min; 
    *t2 = t2max;
    *pmask = mask;

    return err;
}

/**
 * tsls_func:
 * @list: dependent variable plus list of regressors.
 * @ci: %TSLS for regular two-stage least squares, or
 * %SYSTEM if the 2SLS estimates are wanted in the context of
 * estimation of a system of equations.
 * @pZ: pointer to data matrix.
 * @pdinfo: information on the data set.
 * @opt: may contain %OPT_O for printing covariance matrix, 
 * %OPT_R for robust VCV, %OPT_N for no df correction.
 *
 * Estimate the model given in @list by means of Two-Stage Least
 * Squares.  If @ci is %SYSTEM, fitted values from the first-stage
 * regressions are saved with the model: see tsls_get_Xi().
 * 
 * Returns: a #MODEL struct, containing the estimates.
 */

MODEL tsls_func (const int *list, int ci, double ***pZ, DATAINFO *pdinfo,
		 gretlopt opt)
{
    gretl_matrix *Q = NULL;
    char *missmask = NULL;
    int orig_t1 = pdinfo->t1, orig_t2 = pdinfo->t2;
    int *reglist = NULL, *instlist = NULL;
    int *hatlist = NULL, *s2list = NULL;
    int *exolist = NULL;
    int *droplist = NULL;
    int addconst = 0;
    int pos, orig_nvar = pdinfo->v;
    int rlen, ninst, ev = 0;
    int OverIdRank = 0;
    int i;

    MODEL tsls;

    /* initialize model in case we bail out early on error */
    gretl_model_init(&tsls);

    gretl_error_clear();

    pos = gretl_list_separator_position(list);
    if (pos == 0) {
	tsls.errcode = E_PARSE;
	return tsls;
    }

    /* number of elements in regression list (vars preceding the
       separator at "pos") */
    rlen = pos - 1;

    /* number of instruments (vars following the separator at
       "pos") */
    ninst = list[0] - pos;

    if (rlen < 2 || ninst <= 0) {
	tsls.errcode = E_DATA;
	return tsls;
    }

    /* adjust sample range for missing observations */
    tsls.errcode = tsls_adjust_sample(list, &pdinfo->t1, &pdinfo->t2, 
				      (const double **) *pZ, &missmask); 
    if (tsls.errcode) {
	return tsls;
    }

    /* 
       reglist: dependent var plus list of regressors
       instlist: list of instruments
       s2list: regression list for second-stage regression
       hatlist: list of vars to be replaced by fitted values
                for the second-stage regression
    */
    reglist = gretl_list_new(rlen);
    instlist = gretl_list_new(ninst + 1); /* allow for adding const */
    s2list = gretl_list_new(rlen);
    hatlist = gretl_list_new(rlen);

    if (reglist == NULL || instlist == NULL || 
	s2list == NULL || hatlist == NULL) {
	tsls.errcode = E_ALLOC;
	goto bailout;
    }

    /* dependent variable plus regressors: first portion of input
       list */
    for (i=1; i<=rlen; i++) {
	reglist[i] = list[i];
    }    

    /* list of instruments: second portion of input list */
    instlist[0] = ninst;
    for (i=1; i<=instlist[0]; i++) {
	instlist[i] = list[i + pos];
	if (instlist[i] == list[1]) {
	    strcpy(gretl_errmsg, "You can't use the dependent variable as an instrument");
	    tsls.errcode = E_DATA;
	    goto bailout;
	}
    }	

    /* in case we drop any instruments as redundant */
    if (ci == SYSTEM) {
	exolist = gretl_list_copy(instlist);
    }

    /* drop any vars that are all zero, and reshuffle the constant
       into first position among the independent vars 
    */
    tsls_omitzero(reglist, (const double **) *pZ, pdinfo->t1, pdinfo->t2);
    reglist_check_for_const(reglist, (const double **) *pZ, pdinfo);
    tsls_omitzero(instlist, (const double **) *pZ, pdinfo->t1, pdinfo->t2);

    /* determine the list of variables (hatlist) for which we need to
       obtain fitted values in the first stage 
    */
    addconst = tsls_make_hatlist(reglist, instlist, hatlist);

    Q = tsls_Q(instlist, reglist, &droplist,
	       (const double **) *pZ, pdinfo->t1, pdinfo->t2,
	       missmask, &tsls.errcode);
    if (tsls.errcode) {
	goto bailout;
    } 

    /* check for order condition */
    OverIdRank = instlist[0] - reglist[0] + 1;
    if (OverIdRank < 0) {
        sprintf(gretl_errmsg, 
		_("The order condition for identification is not satisfied.\n"
		  "At least %d more instruments are needed."), -OverIdRank);
	tsls.errcode = E_UNSPEC; 
	goto bailout;
    }

    /* hatlist[0] holds the number of fitted vars to create */
    if (dataset_add_series(hatlist[0], pZ, pdinfo)) {
	tsls.errcode = E_ALLOC;
	goto bailout;
    } 

    /* initial composition of second-stage regression list (will be
       modified in the loop below) 
    */
    for (i=0; i<=reglist[0]; i++) {
	s2list[i] = reglist[i];
    }

    /* 
       Deal with the variables for which instruments are needed: loop
       across the list of variables to be instrumented (hatlist),
       form the fitted values as QQ'x_i, and add these fitted values
       into the data matrix Z.
    */

    /* prepare for first-stage F-test, if only one endogenous regressor */
    if (hatlist[0] == 1) {
	ev = hatlist[1];
    }    

    for (i=1; i<=hatlist[0]; i++) {
	int orig = hatlist[i];
	int newv = orig_nvar + i - 1;

	/* form xhat_i = QQ'x_i */
	tsls.errcode = tsls_form_xhat(Q, (*pZ)[orig], (*pZ)[newv], pdinfo, 
				      missmask);
	if (tsls.errcode) {
	    goto bailout;
	}

	/* give the fitted series the same name as the original */
	strcpy(pdinfo->varname[newv], pdinfo->varname[orig]);

	/* substitute newv into the right place in the second-stage
	   regression list */
	replace_list_element(s2list, orig, newv);

	/* update hatlist */
	hatlist[i] = newv;
    }

    /* second-stage regression */
    tsls = lsq(s2list, pZ, pdinfo, OLS, (ci == TSLS)? OPT_NONE : OPT_Z);
    if (tsls.errcode) {
	fprintf(stderr, "tsls_func, stage 2: tsls.errcode = %d\n", tsls.errcode);
	goto bailout;
    }

    if (ci == TSLS && hatlist[0] == 1) {
	compute_first_stage_F(&tsls, ev, reglist, instlist, pZ, pdinfo, opt);
    } 

    if (tsls.list[0] < s2list[0]) {
	/* collinear regressors dropped? If so, adjustments are needed */
	OverIdRank += s2list[0] - tsls.list[0];
	reglist_remove_redundant_vars(&tsls, s2list, reglist);
    }

    /* special: we need to use the original RHS vars to compute
       residuals and associated statistics */
    tsls_residuals(&tsls, reglist, (const double **) *pZ, opt);

    /* compute standard errors */
    if ((opt & OPT_R) || libset_get_bool(USE_QR)) {
	/* QR decomp in force, or robust standard errors called for */
	qr_tsls_vcv(&tsls, (const double **) *pZ, pdinfo, opt);
	if (tsls.errcode) {
	    goto bailout;
	}	
    } 

    /* compute additional statistics (R^2, F, etc.) */
    tsls_extra_stats(&tsls, (const double **) *pZ);

    if (ci == TSLS) {
	tsls_hausman_test(&tsls, reglist, hatlist, pZ, pdinfo);
	if (OverIdRank > 0) {
	    tsls_sargan_test(&tsls, OverIdRank, instlist, pZ, pdinfo);
	}
    }

    /* set command code on the model */
    tsls.ci = TSLS;

    /* write the full 2SLS list (dep. var. and regressors, followed by
       instruments, possibly purged of redundant elements) into the model
       for future reference
    */
    tsls_recreate_full_list(&tsls, reglist, instlist);

    if (droplist != NULL) {
	gretl_model_set_list_as_data(&tsls, "inst_droplist", droplist); 
    }

#if 0
    if (addconst) {
	gretl_model_set_int(&tsls, "addconst", addconst);
    }
#endif

 bailout:

    gretl_matrix_free(Q);
    free(missmask);

    /* save first-stage fitted values, if wanted */
    if (ci == SYSTEM && tsls.errcode == 0) {
	tsls_save_data(&tsls, hatlist, exolist, *pZ, pdinfo);
    }

    /* delete first-stage fitted values from dataset */
    dataset_drop_last_variables(pdinfo->v - orig_nvar, pZ, pdinfo);

    free(hatlist);
    free(reglist); 
    free(instlist);
    free(exolist);
    free(s2list);

    if (tsls.errcode) {
	model_count_minus();
    }

    /* restore original sample range */
    pdinfo->t1 = orig_t1;
    pdinfo->t2 = orig_t2;

    return tsls;
}


