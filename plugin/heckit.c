/* 
 * Copyright (C) Riccardo "Jack" Lucchetti and Allin Cottrell
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "libgretl.h"
#include "matrix_extra.h"
#include "libset.h"
#include "missing_private.h"

#define HDEBUG 0

static int transcribe_heckit_params (MODEL *hm, MODEL *pm, DATAINFO *pdinfo)
{
    double *fullcoeff;
    int ko = hm->ncoeff;
    int kp = pm->ncoeff;
    int k = ko + kp;
    int i, err = 0;

    fullcoeff = malloc(k * sizeof *fullcoeff);
    if (fullcoeff == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	err = gretl_model_allocate_params(hm, k);
    }

    if (!err) {
	for (i=0; i<ko; i++) {
	    strcpy(hm->params[i], pdinfo->varname[hm->list[i+2]]);
	    fullcoeff[i] = hm->coeff[i];
	}
	for (i=0; i<kp; i++) {
	    strcpy(hm->params[i+ko], pdinfo->varname[pm->list[i+2]]);
	    fullcoeff[i+ko] = pm->coeff[i];
	}
    }

    if (!err) {
	free(hm->coeff);
	hm->coeff = fullcoeff;
	hm->ncoeff = k;
	gretl_model_set_coeff_separator(hm, N_("Selection equation"), ko);
	hm->list[hm->list[0]] = 0;
    }
    
    return err;
}

static int transcribe_heckit_vcv (MODEL *pmod, const gretl_matrix *S, 
				  const gretl_matrix *Vp)
{
    int m = S->rows;
    int n = Vp->rows;
    int nvc, tot = m + n;
    int i, j, k = 0;
    double vij;

    if (pmod->vcv != NULL) {
	free(pmod->vcv);
    }

    if (pmod->sderr != NULL) {
	free(pmod->sderr);
    }
    
    nvc = (tot * tot + tot) / 2;

    pmod->vcv = malloc(nvc * sizeof *pmod->vcv);
    pmod->sderr = malloc(tot * sizeof *pmod->sderr);

    if (pmod->vcv == NULL || pmod->sderr == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<tot; i++) {
	for (j=i; j<tot; j++) {
	    if (i < m) {
		vij = (j < m)? gretl_matrix_get(S, i, j) : NADBL;
	    } else {
		vij = gretl_matrix_get(Vp, i-m, j-m);
	    }
	    pmod->vcv[k++] = vij;
	    if (i == j) {
		pmod->sderr[i] = sqrt(vij);
	    }
	}
    }

    return 0;
}

static int make_heckit_vcv (const int *Xl, const int *Zl, int vdelta, 
			    const double **Z, DATAINFO *pdinfo, 
			    MODEL *olsmod, gretl_matrix *Vp)
{
    gretl_matrix *X = NULL;
    gretl_matrix *w = NULL;
    gretl_matrix *Xw = NULL;
    gretl_matrix *W = NULL;

    gretl_matrix *XX = NULL;
    gretl_matrix *XXi = NULL;
    gretl_matrix *XXw = NULL;
    gretl_matrix *XwZ = NULL;
    gretl_matrix *S = NULL;

    int *Xlist = NULL;
    int *Zlist = NULL;
    char *mmask = NULL;

    int i, v, n;
    int nX, nZ;
    int dlist[2] = { 1 , vdelta };
    int err = 0;

    double delta, s2, sigma, bmills, rho, mdelta = 0;

    n = 0;
    for (i=pdinfo->t1; i<=pdinfo->t2; i++) {
	delta = Z[vdelta][i]; 
	if (!na(delta)) {
	    mdelta += delta;
	    n++;
	}
    }
    mdelta /= n;

    bmills = olsmod->coeff[olsmod->ncoeff-1];
    s2 = olsmod->ess / olsmod->nobs + mdelta * bmills * bmills;
    sigma = sqrt(s2);
    rho = bmills / sigma;

    mmask = malloc(pdinfo->t2 - pdinfo->t1 + 1);
    if (mmask == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    v = Xl[0];
    for (i=pdinfo->t1; i<=pdinfo->t2; i++) {
	mmask[i] = na(Z[Xl[v]][i]);
    }

    nX = Xl[0] - 1;
    nZ = Zl[0] - 1;

    Xlist = gretl_list_new(nX);
    Zlist = gretl_list_new(nZ);
    if (Xlist == NULL || Zlist == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    for (i=1; i<=Xlist[0]; i++) {
	Xlist[i] = Xl[i+1];
    }

    for (i=1; i<=Zlist[0]; i++) {
	Zlist[i] = Zl[i+1];
    }

    X = gretl_matrix_data_subset(Xlist, Z, pdinfo->t1, pdinfo->t2, mmask);
    w = gretl_matrix_data_subset(dlist, Z, pdinfo->t1, pdinfo->t2, mmask);
    W = gretl_matrix_data_subset(Zlist, Z, pdinfo->t1, pdinfo->t2, mmask);
    Xw = gretl_matrix_dot_op(w, X, '*', &err);

    XX = gretl_matrix_alloc(nX, nX);
    XXi = gretl_matrix_alloc(nX, nX);
    XXw = gretl_matrix_alloc(nX, nX);
    XwZ = gretl_matrix_alloc(nX, nZ);
    S = gretl_matrix_alloc(nX, nX);

    if (X == NULL || w == NULL || W == NULL || Xw == NULL ||
	XX == NULL || XXi == NULL || XXw == NULL ||
	XwZ == NULL || S == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE, 
			      X, GRETL_MOD_NONE, 
			      XX, GRETL_MOD_NONE);
    gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE, 
			      Xw, GRETL_MOD_NONE, 
			      XXw, GRETL_MOD_NONE);
    gretl_matrix_multiply_mod(Xw, GRETL_MOD_TRANSPOSE, 
			      W, GRETL_MOD_NONE, 
			      XwZ, GRETL_MOD_NONE);

    gretl_matrix_copy_values(XXi, XX);
    err = gretl_invert_symmetric_matrix(XXi); 
    if (err) {
	fprintf(stderr, "make_heckit_vcv: error inverting X'X\n");
	goto bailout;
    }

    gretl_matrix_qform(XwZ, GRETL_MOD_NONE,
		       Vp, S, GRETL_MOD_NONE);
    gretl_matrix_subtract_from(XXw, S);
    gretl_matrix_multiply_by_scalar(XXw, rho * rho);
    gretl_matrix_subtract_from(XX, XXw);
    gretl_matrix_qform(XXi, GRETL_MOD_NONE,
		       XX, S, GRETL_MOD_NONE);
    gretl_matrix_multiply_by_scalar(S, s2);

    olsmod->sigma = sigma;
    olsmod->rho = rho;

    err = transcribe_heckit_vcv(olsmod, S, Vp);

 bailout:

    free(Xlist);
    free(Zlist);
    free(mmask);
    gretl_matrix_free(X);
    gretl_matrix_free(w);
    gretl_matrix_free(W);
    gretl_matrix_free(Xw);

    gretl_matrix_free(XX);
    gretl_matrix_free(XXi);
    gretl_matrix_free(XXw);
    gretl_matrix_free(XwZ);
    gretl_matrix_free(S);

    return err;
}

static char *heckit_suppress_mask (int sel, const int *list, const double **Z,
				   const DATAINFO *pdinfo, int *err)
{
    char *mask = NULL;
    int T = pdinfo->t2 - pdinfo->t1 + 1;
    int i, t, s = 0;

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	if (Z[sel][t]) {
	    for (i=1; i<=list[0]; i++) {
		if (na(Z[list[i]][t])) {
		    s++;
		    break;
		}
	    }
	}
    }

    if (s == T) {
	*err = E_DATA;
	return NULL;
    } else if (s > 0) {
	mask = malloc(pdinfo->n + 1);
	if (mask == NULL) {
	    *err = E_ALLOC;
	    return NULL;
	}

	memset(mask, '0', pdinfo->n);
	mask[pdinfo->n] = 0;

	for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	    if (Z[sel][t]) {
		for (i=1; i<=list[0]; i++) {
		    if (na(Z[list[i]][t])) {
			mask[t] = '1';
			break;
		    }
		}
	    }
	}
    }

    return mask;
}

/* the driver function for the plugin */

MODEL heckit_2step (int *list, double ***pZ, DATAINFO *pdinfo,
		      PRN *prn) 
{
    MODEL hm;
    MODEL probmod;
    int *Xlist = NULL;
    int *Zlist = NULL;
    gretl_matrix *Vp = NULL;
    char *mask = NULL;
    int v = pdinfo->v;
    int N, T, t, sel, nsel;
    double u, xb, delta;
    int err = 0;

    gretl_model_init(&hm);
    gretl_model_init(&probmod);

    err = gretl_list_split_on_separator(list, &Xlist, &Zlist);
    if (err) {
	hm.errcode = err;
	goto bailout;
    } 

    mask = heckit_suppress_mask(Zlist[1], Xlist, (const double **) *pZ, 
				pdinfo, &err);
    if (err) {
	hm.errcode = err;
	goto bailout;
    }  

    if (mask != NULL) {
	copy_to_reference_missmask(mask);
    }

    /* run initial probit */
    probmod = logit_probit(Zlist, pZ, pdinfo, PROBIT, OPT_NONE, prn);
    if (probmod.errcode) {
	free(Xlist);
	free(Zlist);
	goto bailout;
    }

    T = probmod.nobs;
    if (prn != NULL) {
	printmodel(&probmod, pdinfo, OPT_NONE, prn);
    }

    Vp = gretl_vcv_matrix_from_model(&probmod, NULL);

    /* add the auxiliary series to the dataset */

    err = dataset_add_series(2, pZ, pdinfo);
    if (err) {
	hm.errcode = err;
	clear_model(&probmod);
	goto bailout;
    } 

    for (t=pdinfo->t1; t<=pdinfo->t2; t++) {
	sel = ((*pZ)[Zlist[1]][t] == 1.0);
	if (sel && (mask == NULL || mask[t] == '0')) {
	    nsel++;
	    u = probmod.uhat[t];
	    xb = normal_cdf_inverse(probmod.yhat[t]);
	    delta = u * (u + xb);
	    (*pZ)[v][t] = u;
	    (*pZ)[v+1][t] = delta;
	} else {
	    (*pZ)[v][t] = (*pZ)[v+1][t] = NADBL;
	}
    }

    strcpy(pdinfo->varname[v], "lambda");
    strcpy(pdinfo->varname[v+1], "delta");
    gretl_list_append_term(&Xlist, v);

    /* run OLS */
    hm = lsq(Xlist, pZ, pdinfo, OLS, OPT_A);
    hm.ci = HECKIT;
    N = hm.nobs;
    gretl_model_set_int(&hm, "totobs", T);

    /* compute appropriate correction to covariances */
    err = make_heckit_vcv(Xlist, Zlist, v+1, (const double **) *pZ, 
			  pdinfo, &hm, Vp);

    if (err) {
	hm.errcode = err;
    } else {
	err = transcribe_heckit_params(&hm, &probmod, pdinfo);
    }

    if (err) {
	hm.errcode = err;
    }
	
    clear_model(&probmod);
    dataset_drop_last_variables(2, pZ, pdinfo);

 bailout:

    free(Xlist);
    free(Zlist);
    free(mask);
    gretl_matrix_free(Vp);

    return hm;
}
