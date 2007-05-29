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
#define VCV 1

#if VCV
static int make_heckit_vcv (int *Xl, int *Zl, int vdelta, 
			    const double **Z, DATAINFO *pdinfo, 
			    MODEL *olsmod, gretl_matrix *Vp, gretl_matrix *V)
{
    int i, j, k, v, n;
    int t1 = pdinfo->t1;
    int t2 = pdinfo->t2;
    int *Xlist = NULL;
    int *Zlist = NULL;
    int nX, nZ;
    int dlist[2] = { 1 , 1 };
    dlist[1] = vdelta;
    char *mmask;
    int err = 0;

    double delta, s2, sigma, bmills, rho, rho2, mdelta = 0;

    n = 0;
    for (i=pdinfo->t1; i<=pdinfo->t2; i++) {
	delta = Z[vdelta][i]; 
	if (!na(delta)) {
	    mdelta += delta;
	    n++;
	}
    }
    mdelta = mdelta / n;

    bmills = olsmod->coeff[olsmod->ncoeff-1];
    s2 = olsmod->ess/olsmod->nobs + mdelta*(bmills*bmills);
    sigma = sqrt(s2);
    rho = bmills/sigma;
    rho2 = rho*rho;

    v = Xl[0];
    mmask = malloc((t2-t1+1)* sizeof *mmask);
    for (i=t1; i<=t2; i++) {
	mmask[i] = na(Z[Xl[v]][i]);
    }

    nX = Xl[0] - 1;
    nZ = Zl[0] - 1;

    Xlist = gretl_list_new(nX);
    Zlist = gretl_list_new(nZ);

    for (i=1; i<=Xlist[0]; i++) {
	Xlist[i] = Xl[i+1];
    }

    for (i=1; i<=Zlist[0]; i++) {
	Zlist[i] = Zl[i+1];
    }

    gretl_matrix *X = NULL;
    gretl_matrix *w = NULL;
    gretl_matrix *Xw = NULL;
    gretl_matrix *W = NULL;

    X = gretl_matrix_data_subset(Xlist, Z, t1, t2, mmask);
    w = gretl_matrix_data_subset(dlist, Z, t1, t2, mmask);
    W = gretl_matrix_data_subset(Zlist, Z, t1, t2, mmask);
    Xw = gretl_matrix_dot_op (w, X, '*', &err);

    gretl_matrix *XX = gretl_matrix_alloc(nX, nX);
    gretl_matrix *XXi = gretl_matrix_alloc(nX, nX);
    gretl_matrix *XXw = gretl_matrix_alloc(nX, nX);
    gretl_matrix *XwZ = gretl_matrix_alloc(nX, nZ);
    gretl_matrix *tmp = gretl_matrix_alloc(nX, nX);
    gretl_matrix *S = gretl_matrix_alloc(nX, nX);

    err = gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE, 
				    X, GRETL_MOD_NONE, 
				    XX, GRETL_MOD_NONE);
    err = gretl_matrix_multiply_mod(X, GRETL_MOD_TRANSPOSE, 
				    Xw, GRETL_MOD_NONE, 
				    XXw, GRETL_MOD_NONE);
    err = gretl_matrix_multiply_mod(Xw, GRETL_MOD_TRANSPOSE, 
				    W, GRETL_MOD_NONE, 
				    XwZ, GRETL_MOD_NONE);

    gretl_matrix_free(X);
    gretl_matrix_free(w);
    gretl_matrix_free(W);
    gretl_matrix_free(Xw);

    XXi = gretl_matrix_copy(XX);
    err = gretl_invert_symmetric_matrix(XXi); 
    err = gretl_matrix_qform (XwZ, GRETL_MOD_NONE,
			      Vp, tmp, GRETL_MOD_NONE);
    err = gretl_matrix_subtract_from(XXw, tmp);
    gretl_matrix_multiply_by_scalar (XXw, rho2);
    err = gretl_matrix_subtract_from(XX, XXw);
    err = gretl_matrix_qform (XXi, GRETL_MOD_NONE,
			      XX, S, GRETL_MOD_NONE);
    gretl_matrix_multiply_by_scalar (S, s2);

    olsmod->sigma = sigma;
    olsmod->rho = rho;

    k = 0;
    makevcv(olsmod, sigma);
    for (i=0; i<nX; i++) {
	olsmod->sderr[i] = sqrt(gretl_matrix_get(S,i,i));
	for (j=i; j<nX; j++) {
	    olsmod->vcv[k++] = gretl_matrix_get(S,i,j);
	}
    }

    gretl_matrix_free(XX);
    gretl_matrix_free(XXi);
    gretl_matrix_free(XXw);
    gretl_matrix_free(XwZ);
    gretl_matrix_free(tmp);
    gretl_matrix_free(S);

    return err;
}
#endif

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

MODEL heckit_estimate (const int *list, double ***pZ, DATAINFO *pdinfo,
		      PRN *prn) 
{
    MODEL hm;
    MODEL probmod;
    int *Xlist = NULL;
    int *Zlist = NULL;
    gretl_matrix *V = NULL;
    gretl_matrix *Vp = NULL;

    char *mask = NULL;
    int v = pdinfo->v;
    int i, N, T, t, k, sel, nsel;
    double u, xb, delta;
    int err = 0;

    gretl_model_init(&hm);

    /* adjust sample --- TODO */

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
    fprintf(stderr, "Probit obs = %d\n", T);
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

    clear_model(&probmod);

    strcpy(pdinfo->varname[v], "Mills factor");
    strcpy(pdinfo->varname[v+1], "delta");
    gretl_list_append_term(&Xlist, v);

    /* run OLS */
    hm = lsq(Xlist, pZ, pdinfo, OLS, OPT_A);
    hm.ci = HECKIT;
    k = hm.ncoeff;
    N = hm.nobs;
    gretl_model_set_int(&hm, "totobs", T);

    /* compute  appropriate correction to covariances */
    make_heckit_vcv(Xlist, Zlist, v+1, (const double **) *pZ, 
		    pdinfo, &hm, Vp, V);

    gretl_model_allocate_params(&hm, k);

    for (i=0; i<k-1; i++) {
	strcpy(hm.params[i], pdinfo->varname[hm.list[i+2]]);
    }
    strcpy(hm.params[k-1], "lambda");
	
    dataset_drop_last_variables(2, pZ, pdinfo);

 bailout:

    free(Xlist);
    free(Zlist);
    free(mask);

    return hm;
}
