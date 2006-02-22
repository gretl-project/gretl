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

#include "libgretl.h"
#include "objstack.h"
#include "modelspec.h"

#include "glib.h"

#define MDEBUG 0

struct model_data_item_ {
    char *key;
    void *ptr;
    size_t size;
    void (*destructor) (void *);
};

struct ModelTest_ {
    int type;
    int order;
    char *param;
    unsigned char teststat;
    int dfn, dfd;
    double value;
    double pvalue;
};

static void free_model_data_item (model_data_item *item)
{
    if (item->destructor != NULL) {
	(*item->destructor)(item->ptr);
    } else {
	free(item->ptr);
    }
    free(item->key);
    free(item);
}

static model_data_item *create_data_item (const char *key, void *ptr, size_t size,
					  void (*destructor) (void *))
{
    model_data_item *item;

    item = malloc(sizeof *item);
    if (item != NULL) {
	item->key = gretl_strdup(key);
	if (item->key == NULL) {
	    free(item);
	    item = NULL;
	} else {
	    item->ptr = ptr;
	    item->size = size;
	    item->destructor = destructor;
	}
    }

    return item;
}

/**
 * gretl_model_set_data_with_destructor:
 * @pmod: pointer to #MODEL.
 * @key: key string for data, used in retrieval.
 * @ptr: data-pointer to be attached to model.
 * @size: size of data in bytes.
 * @destructor: pointer to function that should be used to free
 * the data-pointer in question.
 *
 * Attaches data to @pmod: the data can be retrieved later using
 * gretl_model_get_data().  Note that the data are not "physically"
 * copied to the model; simply, @ptr is recorded on the model.
 * This means that the data referenced by the pointer now in 
 * effect belong to @pmod.  When @pmod is cleared with clear_model(),
 * @destructor will be invoked with @ptr as its single argument.
 * If a simple "free" is OK for freeing the data, you can use
 * gretl_model_set_data() instead.
 *
 * The @size is needed in case the model is copied with
 * copy_model(), in which case the target of the copying
 * operation receives a newly allocated copy of the data in
 * question.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_data_with_destructor (MODEL *pmod, const char *key, void *ptr, 
					  size_t size, void (*destructor) (void *))
{
    model_data_item **items;
    model_data_item *item;
    int n_items = pmod->n_data_items + 1;

    items = realloc(pmod->data_items, n_items * sizeof *items);
    if (items == NULL) return 1;

    pmod->data_items = items;

    item = create_data_item(key, ptr, size, destructor);
    if (item == NULL) return 1;

    pmod->data_items[n_items - 1] = item;
    pmod->n_data_items += 1;

    return 0;
}

/**
 * gretl_model_set_data:
 * @pmod: pointer to #MODEL.
 * @key: key string for data, used in retrieval.
 * @ptr: data-pointer to be attached to model.
 * @size: size of data in bytes.
 *
 * Attaches data to @pmod: the data can be retrieved later using
 * gretl_model_get_data().  Note that the data are not "physically"
 * copied to the model; simply, @ptr is recorded on the model.
 * This means that the data referenced by the pointer now in 
 * effect belong to @pmod.  The data pointer will be freed when 
 * @pmod is cleared with clear_model().  If the data has deep
 * structure that requires special treatment on freeing, use
 * gretl_model_set_data_with_destructor() instead.
 *
 * The @size is needed in case the model is copied with
 * copy_model(), in which case the target of the copying
 * operation receives a newly allocated copy of the data in
 * question.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_data (MODEL *pmod, const char *key, void *ptr, size_t size)
{
    return gretl_model_set_data_with_destructor(pmod, key, ptr, size, NULL);
}

/**
 * gretl_model_set_list_as_data:
 * @pmod: pointer to #MODEL.
 * @key: key string, used in retrieval.
 * @list: list to attach.
 *
 * Attaches @list to @pmod as data, recoverable via the key @key 
 * using gretl_model_get_data().
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_list_as_data (MODEL *pmod, const char *key, int *list)
{
    size_t size = (list[0] + 1) * sizeof *list;

    return gretl_model_set_data_with_destructor(pmod, key, (void *) list, 
						size, NULL);
}


/**
 * gretl_model_set_int:
 * @pmod: pointer to #MODEL.
 * @key: key string, used in retrieval.
 * @val: integer value to set.
 *
 * Records an integer value on a model: the value can be retrieved 
 * later using gretl_model_get_int(), using the appropriate @key.  
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_int (MODEL *pmod, const char *key, int val)
{
    int *valp;
    int err;

    /* if value is already set, reset it */
    valp = gretl_model_get_data(pmod, key);
    if (valp != NULL) {
	*valp = val;
	return 0;
    }

    valp = malloc(sizeof *valp);
    if (valp == NULL) return 1;

    *valp = val;

    err = gretl_model_set_data(pmod, key, valp, sizeof(int));
    if (err) free(valp);

    return err;
}

/**
 * gretl_model_set_double:
 * @pmod: pointer to model.
 * @key: key string, used in retrieval.
 * @val: double-precision value to set.
 *
 * Records a floating-point value on @pmod: the value can be 
 * retrieved later using gretl_model_get_double() with the
 * appropriate @key. 
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_double (MODEL *pmod, const char *key, double val)
{
    double *valp;
    int err;

    /* if value is already set, reset it */
    valp = gretl_model_get_data(pmod, key);
    if (valp != NULL) {
	*valp = val;
	return 0;
    }

    valp = malloc(sizeof *valp);
    if (valp == NULL) return 1;

    *valp = val;

    err = gretl_model_set_data(pmod, key, valp, sizeof(double));
    if (err) free(valp);

    return err;
}

/**
 * gretl_model_get_data_and_size:
 * @pmod: pointer to model.
 * @key: key string.
 * @sz: location to receive the size of the data.
 *
 * Returns: the data pointer identified by @key, or %NULL on failure.
 */

void *gretl_model_get_data_and_size (const MODEL *pmod, const char *key,
				     size_t *sz)
{
    int i;

    for (i=0; i<pmod->n_data_items; i++) {
	if (!strcmp(key, pmod->data_items[i]->key)) {
	    if (sz != NULL) {
		*sz = pmod->data_items[i]->size;
	    }
	    return pmod->data_items[i]->ptr;
	}
    }

    return NULL;
}

/**
 * gretl_model_get_data:
 * @pmod: pointer to model.
 * @key: key string.
 *
 * Returns: the data pointer identified by @key, or %NULL on failure.
 */

void *gretl_model_get_data (const MODEL *pmod, const char *key)
{
    return gretl_model_get_data_and_size(pmod, key, NULL);
}

/**
 * gretl_model_get_int:
 * @pmod: pointer to model.
 * @key: key string.
 *
 * Returns: the integer value identified by @key, or 0 on failure.
 */

int gretl_model_get_int (const MODEL *pmod, const char *key)
{
    int *valp = NULL;
    int i;

    for (i=0; i<pmod->n_data_items; i++) {
	if (!strcmp(key, (pmod->data_items[i])->key)) {
	    valp = (int *) (pmod->data_items[i])->ptr;
	    return *valp;
	}
    }

    return 0;
}

/**
 * gretl_model_get_double:
 * @pmod: pointer to model.
 * @key: key string.
 *
 * Returns: the double-precision value identified by @key, or 
 * #NADBL on failure.
 */

double gretl_model_get_double (const MODEL *pmod, const char *key)
{
    double *valp = NULL;
    int i;

    for (i=0; i<pmod->n_data_items; i++) {
	if (!strcmp(key, (pmod->data_items[i])->key)) {
	    valp = (double *) (pmod->data_items[i])->ptr;
	    return *valp;
	}
    }

    return NADBL;
}

static void adjust_vecm_name (const char *orig, char *cname)
{
    int cnum;
    char cc;

    if (sscanf(orig, "EC%d%c", &cnum, &cc) == 2) {
	sprintf(cname, "EC%d", cnum);
    } else {
	strcpy(cname, orig);
    }
}

static void make_cname (const char *orig, char *cname)
{
    char *p;

    if (orig == NULL || *orig == 0) {
	return;
    }

    p = strrchr(orig, '_');

    if (p == NULL) {
	strcpy(cname, orig);
    } else {
	unsigned char c = (unsigned char) *(p + 1);

	if (isdigit(c)) {
	    int lag = atoi(++p);

	    sprintf(cname, "ut^2(-%d)", lag);
	}
    }
}

/**
 * gretl_model_get_param_name:
 * @pmod: pointer to model.
 * @pdinfo: dataset information.
 * @i: index number for parameter, zero-based, corresponding
 * to position in the %coeff array in @pmod.
 * @targ: string into which to write param name.
 *
 * Writes the appropriate parameter name into @targ, which
 * should be at least #VNAMELEN bytes long.  Usually this is
 * the name of a variable in the dataset, but sometimes it is
 * a special string (e.g. for nonlinear models).
 *
 * Returns: @targ.
 */

char *gretl_model_get_param_name (const MODEL *pmod, const DATAINFO *pdinfo,
				  int i, char *targ)
{
    *targ = '\0';

    if (pmod != NULL) {
	/* special treatment for ARCH, ARMA, GARCH, NLS, MLE */
	if (pmod->aux == AUX_ARCH) {
	    make_cname(pdinfo->varname[pmod->list[i + 2]], targ);
	} else if (pmod->ci == NLS || pmod->ci == MLE ||
		   pmod->ci == ARMA || pmod->ci == GARCH) {
	    strcpy(targ, pmod->params[i + 1]);
	} else if (pmod->aux == AUX_VECM) {
	    adjust_vecm_name(pdinfo->varname[pmod->list[i + 2]], targ);
	} else {
	    strcpy(targ, pdinfo->varname[pmod->list[i + 2]]);
	}
    }

    return targ;
}

void free_coeff_intervals (CoeffIntervals *cf)
{
    int i;

    free(cf->coeff);
    free(cf->maxerr);

    if (cf->names != NULL) {
	for (i=0; i<cf->ncoeff; i++) {
	    free(cf->names[i]);
	}
	free(cf->names);
    }

    free(cf);
}

/**
 * gretl_model_get_coeff_intervals:
 * @pmod: pointer to gretl model.
 * @pdinfo: dataset information.
 *
 * Save the 95 percent confidence intervals for the parameter
 * estimates in @pmod.
 * 
 * Returns: pointer to #CONFINT struct containing the results.
 */

CoeffIntervals *
gretl_model_get_coeff_intervals (const MODEL *pmod, 
				 const DATAINFO *pdinfo)
{
    CoeffIntervals *cf;
    double t = tcrit95(pmod->dfd);
    char pname[24];
    int i, err = 0;

    cf = malloc(sizeof *cf);
    if (cf == NULL) {
	return NULL;
    }

    cf->ncoeff = pmod->ncoeff;
    cf->df = pmod->dfd;
    cf->ifc = pmod->ifc;

    cf->coeff = NULL;
    cf->maxerr = NULL;
    cf->names = NULL;

    cf->coeff = malloc(cf->ncoeff * sizeof *cf->coeff);
    if (cf->coeff == NULL) {
	err = 1;
	goto bailout;
    }

    cf->maxerr = malloc(cf->ncoeff * sizeof *cf->maxerr);
    if (cf->maxerr == NULL) {
	err = 1;
	goto bailout;
    }

    cf->names = malloc(cf->ncoeff * sizeof *cf->names);
    if (cf->names == NULL) {
	err = 1;
	goto bailout;
    }    

    for (i=0; i<cf->ncoeff && !err; i++) { 
	cf->coeff[i] = pmod->coeff[i];
	cf->maxerr[i] = (pmod->sderr[i] > 0)? t * pmod->sderr[i] : 0.0;
	gretl_model_get_param_name(pmod, pdinfo, i, pname);
	cf->names[i] = gretl_strdup(pname);
	if (cf->names[i] == NULL) {
	    int j;
	    
	    for (j=0; j<i; j++) {
		free(cf->names[i]);
	    }
	    free(cf->names);
	    cf->names = NULL;
	    err = 1;
	}
    }

 bailout:

    if (err) {
	free_coeff_intervals(cf);
	cf = NULL;
    }

    return cf;
}

static int gretl_is_arima_model (const MODEL *pmod)
{
    int d = (int) gretl_model_get_data(pmod, "arima_d");
    int D = (int) gretl_model_get_data(pmod, "arima_D");

    return (d > 0 || D > 0);
}

/**
 * gretl_arma_model_nonseasonal_AR_order:
 * @pmod: pointer to gretl model.
 *
 * Returns: the non-seasonal autoregressive order of @pmod, or 0 if
 * @pmod is not an ARMA model.
 */

int gretl_arma_model_nonseasonal_AR_order (const MODEL *pmod)
{
    int p = 0;

    if (pmod->ci == ARMA) {
	p = pmod->list[1];
    }

    return p;
}

/**
 * gretl_arma_model_nonseasonal_MA_order:
 * @pmod: pointer to gretl model.
 *
 * Returns: the non-seasonal moving-average order of @pmod, or 0 if
 * @pmod is not an ARMA model.
 */

int gretl_arma_model_nonseasonal_MA_order (const MODEL *pmod)
{
    int q = 0;

    if (pmod->ci == ARMA) {
	if (gretl_is_arima_model(pmod)) {
	    q = pmod->list[3];
	} else {
	    q = pmod->list[2];
	}
    }

    return q;
}

/**
 * gretl_arma_model_max_AR_lag:
 * @pmod: pointer to gretl model.
 *
 * Returns: the maximum autoregressive lag in @pmod, or 0 if
 * @pmod is not an ARMA model.
 */

int gretl_arma_model_max_AR_lag (const MODEL *pmod)
{
    int pmax = 0;

    if (pmod->ci == ARMA) {
	int p, P, pd;

	p = gretl_arma_model_nonseasonal_AR_order(pmod);
	P = gretl_model_get_int(pmod, "arma_P");

	if (P == 0) {
	    pmax = p;
	} else {
	    pd = gretl_model_get_int(pmod, "arma_pd");
	    pmax = P * pd + p;
	}
    }

    return pmax;
}

/**
 * gretl_arma_model_max_MA_lag:
 * @pmod: pointer to gretl model.
 *
 * Returns: the maximum moving-average lag in @pmod, or 0 if
 * @pmod is not an ARMA model.
 */

int gretl_arma_model_max_MA_lag (const MODEL *pmod)
{
    int qmax = 0;

    if (pmod->ci == ARMA) {
	int q, Q, pd;

	q = gretl_arma_model_nonseasonal_MA_order(pmod);
	Q = gretl_model_get_int(pmod, "arma_Q");
	
	if (Q == 0) {
	    qmax = q;
	} else {
	    pd = gretl_model_get_int(pmod, "arma_pd");
	    qmax = Q * pd + q;
	}
    }

    return qmax;
}

/**
 * gretl_arma_model_get_AR_MA_coeffs:
 * @pmod: pointer to gretl model.
 * @arvec: pointer to receive AR coeff vector.
 * @mavec: pointer to receive MA coeff vector.
 *
 * Creates allocated copies of the AR and MA coefficient vectors from
 * @pmod.  If @pmod includes seasonal ARMA terms, the coefficient vectors
 * are suitably expanded, and include the interactions, if any, between 
 * seasonal and non-seasonal terms.  The length of these vectors can
 * be determined using gretl_arma_model_get_max_AR_lag() and
 * gretl_arma_model_get_max_MA_lag() respectively.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_arma_model_get_AR_MA_coeffs (const MODEL *pmod,
				       double **arvec,
				       double **mavec)
{
    double *ac = NULL;
    double *mc = NULL;
    int err = 0;

    if (pmod->ci != ARMA) {
	err = 1;
    } else {
	const double *phi = NULL, *Phi = NULL;
	const double *theta = NULL, *Theta = NULL;

	int p = gretl_arma_model_nonseasonal_AR_order(pmod);
	int q = gretl_arma_model_nonseasonal_MA_order(pmod);
	int P = gretl_model_get_int(pmod, "arma_P");
	int Q = gretl_model_get_int(pmod, "arma_Q");
	int pd = gretl_model_get_int(pmod, "arma_pd");
	int pmax, qmax;
	int i, j, k;

	pmax = (P > 0)? P * pd + p : p;
	qmax = (Q > 0)? Q * pd + q : q;
	
	if (pmax > 0) {
	    ac = malloc(pmax * sizeof *ac);
	    if (ac == NULL) {
		err = E_ALLOC;
	    }
	}

	if (!err && qmax > 0) {
	    mc = malloc(qmax * sizeof *ac);
	    if (mc == NULL) {
		free(ac);
		ac = NULL;
		err = E_ALLOC;
	    }
	}

	if (!err) {
	    phi = pmod->coeff + pmod->ifc; /* non-seasonal AR coeffs */
	    Phi = phi + p;                 /* seasonal AR coeffs */
	    theta = Phi + P;               /* non-seasonal MA coeffs */
	    Theta = theta + q;             /* seasonal MA coeffs */
	}

	if (ac != NULL) {
	    for (i=0; i<p; i++) {
		ac[i] = phi[i];
	    }	    
	    if (P > 0) {
		for (i=p; i<pmax; i++) {
		    ac[i] = 0.0;
		}
		for (i=0; i<P; i++) {
		    k = pd * (i+1) - 1;
		    ac[k] += Phi[i];
		    for (j=0; j<p; j++) {
			ac[k+j+1] += Phi[i] * phi[j];
		    }
		}
	    }		
	}

	if (mc != NULL) {
	    for (i=0; i<q; i++) {
		mc[i] = theta[i];
	    }	    
	    if (Q > 0) {
		for (i=q; i<qmax; i++) {
		    mc[i] = 0.0;
		}
		for (i=0; i<Q; i++) {
		    k = pd * (i+1) - 1;
		    mc[k] += Theta[i];
		    for (j=0; j<q; j++) {
			mc[k+j+1] += Theta[i] * theta[j];
		    }
		}
	    }		
	}
    }

    if (!err) {
	*arvec = ac;
	*mavec = mc;
    }

    return err;
}

/**
 * gretl_arma_model_get_x_coeffs:
 * @pmod: pointer to gretl model.
 *
 * Returns: pointer to the array of coefficients on the exogenous
 * regressors in @pmod, or %NULL if the model is not ARMA or if there
 * are no such regressors.
 */

const double *gretl_arma_model_get_x_coeffs (const MODEL *pmod)
{
    const double *xc = NULL;

    if (pmod->ci == ARMA && gretl_model_get_int(pmod, "armax")) {
	xc = pmod->coeff;
	xc += pmod->ifc;
	xc += gretl_arma_model_nonseasonal_AR_order(pmod);
	xc += gretl_arma_model_nonseasonal_MA_order(pmod);
	xc += gretl_model_get_int(pmod, "arma_P");
	xc += gretl_model_get_int(pmod, "arma_Q");
    }

    return xc;
}

static int arma_depvar_pos (const MODEL *pmod)
{
    int seasonal = 0;
    int arima = 0;
    int dvpos;

    if (gretl_model_get_int(pmod, "arma_P") ||
	gretl_model_get_int(pmod, "arima_D") ||
	gretl_model_get_int(pmod, "arma_Q")) {
	seasonal = 1;
    }

    if (gretl_model_get_int(pmod, "arima_d") ||
	gretl_model_get_int(pmod, "arima_D")) {
	arima = 1;
    }

    if (arima) {
	dvpos = (seasonal)? 9 : 5;
    } else {
	dvpos = (seasonal)? 7 : 4;
    }

    return dvpos;
}
    
/**
 * gretl_model_get_depvar:
 * @pmod: pointer to gretl model.
 *
 * Returns: the ID number of the dependent variable in @pmod.
 */

int gretl_model_get_depvar (const MODEL *pmod)
{
    int dv = 0;

    if (pmod != NULL && pmod->list != NULL) {
	if (pmod->ci == GARCH) {
	    dv = pmod->list[4];
	} else if (pmod->ci == ARMA) {
	    dv = pmod->list[arma_depvar_pos(pmod)];
	} else {
	    dv = pmod->list[1];
	}
    }

    return dv;
}

/**
 * gretl_model_get_x_list:
 * @pmod: pointer to gretl model.
 *
 * Returns: an allocated copy of the list of independent
 * variables included in @pmod, or %NULL on failure. 
 */

int *gretl_model_get_x_list (const MODEL *pmod)
{
    int *list = NULL;
    int i, nx;

    if (pmod->ci == ARMA) {
	int start = arma_depvar_pos(pmod);

	nx = pmod->list[0] - start + pmod->ifc;
	if (nx > 0) {
	    list = gretl_list_new(nx);
	    if (list != NULL) {
		if (pmod->ifc) {
		    list[1] = 0;
		    for (i=2; i<=list[0]; i++) {
			list[i] = pmod->list[i + start - 1];
		    }
		} else {
		    for (i=1; i<=list[0]; i++) {
			list[i] = pmod->list[i + start];
		    }
		}		    
	    }
	}
    } else if (pmod->ci == GARCH) {
	nx = pmod->list[0] - 4;
	if (nx > 0) {
	    list = gretl_list_new(nx);
	    if (list != NULL) {
		for (i=1; i<=list[0]; i++) {
		    list[i] = pmod->list[i + 4];
		}
	    }
	}
    } else if (pmod->ci != NLS && pmod->ci != MLE) {
	nx = pmod->ncoeff;
	list = gretl_list_new(nx);
	if (list != NULL) {
	    for (i=1; i<=list[0]; i++) {
		list[i] = pmod->list[i + 1];
	    }
	}
    }

    return list;
}

/**
 * gretl_model_new_vcv:
 * @pmod: pointer to model.
 * @nelem: pointer to receive number of elements in
 * the packed array, or %NULL;
 * 
 * Allocates space for a packed coefficient covariance matrix
 * in @pmod (if such space is not already allocated).  Sets
 * all entries in the array to zero.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int gretl_model_new_vcv (MODEL *pmod, int *nelem)
{
    int nv = pmod->ncoeff;
    int nxpx = (nv * nv + nv) / 2; 
    int i, err = 0;

    if (pmod->vcv == NULL) {
	/* not already allocated */
	pmod->vcv = malloc(nxpx * sizeof *pmod->vcv);
	if (pmod->vcv == NULL) {
	    err = E_ALLOC;
	} 
    }

    if (pmod->vcv != NULL) {
	for (i=0; i<nxpx; i++) {
	    pmod->vcv[i] = 0.0;
	}
	if (nelem != NULL) {
	    *nelem = nxpx;
	}
    }	

    return err;
}

/**
 * gretl_model_get_vcv:
 * @pmod: pointer to model.
 * @pdinfo: dataset information.
 * 
 * Supplies the caller with a copy of the variance-covariance 
 * matrix for the parameter estimates in @pmod.  See also
 * free_vcv().  To get the covariance matrix in gretl_matrix
 * format, see gretl_vcv_matrix_from_model().
 *
 * Returns: #VMatrix struct or %NULL on error.
 */

VMatrix *gretl_model_get_vcv (MODEL *pmod, const DATAINFO *pdinfo)
{
    char varname[VNAMELEN];
    int i, nt, nc = pmod->ncoeff;
    VMatrix *vcv;

    vcv = vmatrix_new();

    if (vcv == NULL) {
	return NULL;
    }

    vcv->names = create_strings_array(nc);
    if (vcv->names == NULL) {
	free(vcv);
	return NULL;
    }

    for (i=0; i<nc; i++) {
	gretl_model_get_param_name(pmod, pdinfo, i, varname);
	vcv->names[i] = gretl_strdup(varname);
	if (vcv->names[i] == NULL) {
	    free_vmatrix(vcv);
	    return NULL;
	}
    }

    if (pmod->vcv == NULL && makevcv(pmod)) {
	free_vmatrix(vcv);
	return NULL;
    }

    /* calculate number of elements in vcv */
    nt = (nc * nc + nc) / 2;

    /* copy vcv */
    vcv->vec = copyvec(pmod->vcv, nt);
    if (vcv->vec == NULL) {
	free_vmatrix(vcv);
	return NULL;
    }

    vcv->ci = pmod->ci;
    vcv->dim = nc;
    vcv->t1 = pmod->t1;
    vcv->t2 = pmod->t2;
    
    return vcv;
}

/**
 * impose_model_smpl:
 * @pmod: pointer to model.
 * @pdinfo: dataset information.
 *
 * Sets on @pdinfo the sample range (starting and ending
 * observations) that was in effect when @pmod was estimated.
 * This is not always the same as the data range over which
 * @pmod was actually estimated (e.g. in case of 
 * autoregressive models, where observations are dropped
 * to allow for lags).
 */

void impose_model_smpl (const MODEL *pmod, DATAINFO *pdinfo)
{
    pdinfo->t1 = pmod->smpl.t1;
    pdinfo->t2 = pmod->smpl.t2;
}

#ifdef HAVE_X12A

static void maybe_delete_x12_file (const MODEL *pmod)
{
    char *fname = NULL;

    fname = gretl_model_get_data(pmod, "x12a_output");
    if (fname != NULL) remove(fname);
}

#endif

static void destroy_all_data_items (MODEL *pmod)
{
    int i;

    if (pmod->n_data_items == 0) {
	return;
    }

#ifdef HAVE_X12A
    maybe_delete_x12_file(pmod);
#endif

    for (i=0; i<pmod->n_data_items; i++) {
	free_model_data_item(pmod->data_items[i]);
    }

    free(pmod->data_items);
    pmod->data_items = NULL;
}

static int discard_model_data_item (MODEL *pmod, const char *key,
				    int free_data)
{
    model_data_item *junk = NULL;
    int i, targ = 0;
    int err = 0;

    for (i=0; i<pmod->n_data_items; i++) {
	if (!strcmp(key, pmod->data_items[i]->key)) {
	    junk = pmod->data_items[i];
	    targ = i;
	    break;
	}
    }

    if (junk == NULL) {
	err = 1;
    } else {
	model_data_item **items;
	int n_items = pmod->n_data_items - 1;

	for (i=targ; i<n_items; i++) {
	    pmod->data_items[i] = pmod->data_items[i+1];
	}

	items = realloc(pmod->data_items, n_items * sizeof *items);
	if (items != NULL) {
	    pmod->data_items = items;
	}

	pmod->n_data_items -= 1;

	if (free_data) {
	    /* deep free the data item */
	    free_model_data_item(junk);
	} else {
	    /* just free the item, not the actual data */
	    free(junk->key);
	    free(junk);
	}
    }

    return err;
}

/**
 * gretl_model_destroy_data_item:
 * @pmod: pointer to model.
 * @key: key string.
 *
 * Looks up the data pointer, attached to @pmod, that is
 * identified by @key, and if a pointer is found, frees
 * it (or applies the destructor function that was set for
 * the item, if any) and removes it from the model's list of 
 * data items.  If you want to remove the item from the
 * model's list without freeing the underlying data pointer,
 * use gretl_model_detach_data_item().
 *
 * Returns: 0 on success, 1 on failure (pointer not found).
 */

int gretl_model_destroy_data_item (MODEL *pmod, const char *key)
{
    return discard_model_data_item(pmod, key, 1);
}

/**
 * gretl_model_detach_data_item:
 * @pmod: pointer to model.
 * @key: key string.
 *
 * Looks up the data item, attached to @pmod, that is
 * identified by @key, and if an item is found, removes
 * it from the model's list of such items.  The data
 * pointer associated with @key is not touched.  If you
 * want the underlying resources associated with @key to be 
 * freed, use gretl_model_destroy_data_item().
 *
 * Returns: 0 on success, 1 on failure (key not found).
 */

int gretl_model_detach_data_item (MODEL *pmod, const char *key)
{
    return discard_model_data_item(pmod, key, 0);
}

/**
 * gretl_model_set_auxiliary:
 * @pmod: pointer to model.
 * @aux: code indicating a model's function in an
 * auxiliary role (typically, in relation to a hypothesis
 * test on another model).
 *
 * Sets an auxiliary code on @pmod, which may be relevant
 * for how the model is printed.
 */

void gretl_model_set_auxiliary (MODEL *pmod, ModelAuxCode aux)
{
    pmod->aux = aux;
}

static void gretl_model_init_pointers (MODEL *pmod)
{
    pmod->list = NULL;
    pmod->submask = NULL;
    pmod->missmask = NULL;
    pmod->coeff = NULL;
    pmod->sderr = NULL;
    pmod->yhat = NULL;
    pmod->uhat = NULL;
    pmod->xpx = NULL;
    pmod->vcv = NULL;
    pmod->arinfo = NULL;
    pmod->name = NULL;
    pmod->params = NULL;
    pmod->tests = NULL;
    pmod->data = NULL;
    pmod->dataset = NULL;
    pmod->data_items = NULL;
}

/**
 * gretl_model_init:
 * @pmod: pointer to model.
 *
 * Initializes a gretl #MODEL, including setting its pointer members
 * to %NULL. This initialization should be done if the caller has
 * declared a #MODEL struct directly, rather than obtaining a pointer to
 * #MODEL using gretl_model_new() (in which case the initialization is
 * done automatically).
 */

void gretl_model_init (MODEL *pmod)
{
    int i;

    if (pmod == NULL) return;

#if MDEBUG
    fprintf(stderr, "gretl_model_init: pmod at %p\n", (void *) pmod);
#endif

    pmod->ID = 0;
    pmod->refcount = 1;
    pmod->full_n = 0;

    pmod->ntests = 0;
    pmod->nparams = 0;
    pmod->errcode = 0;
    pmod->ci = 0;
    pmod->ifc = 0;
    pmod->aux = AUX_NONE;

    for (i=0; i<C_MAX; i++) {
	pmod->criterion[i] = NADBL;
    }

    gretl_model_init_pointers(pmod);
    pmod->n_data_items = 0;

    *gretl_msg = '\0';
}

/**
 * gretl_model_smpl_init:
 * @pmod: pointer to model.
 * @pdinfo: dataset information.
 *
 * Records the start and end of the current sample range in
 * the model @pmod, which may be necessary for future reference
 * if a hypothesis test is to be performed.  Note that this
 * sample range may not be the same as the data range over
 * which the model is actually estimated (for example, in the
 * case of autoregressive models where observations have to
 * be dropped to allow for lags).
 */

void gretl_model_smpl_init (MODEL *pmod, const DATAINFO *pdinfo)
{
    pmod->smpl.t1 = pdinfo->t1;
    pmod->smpl.t2 = pdinfo->t2;
}

static MODEL *real_gretl_model_new (int protected)
{
    MODEL *pmod = malloc(sizeof *pmod);

#if MDEBUG
    fprintf(stderr, "real_gretl_model_new: model at %p\n", (void *) pmod);
#endif

    if (pmod != NULL) {
	gretl_model_init(pmod);
    }

    if (protected) {
#if MDEBUG
	fprintf(stderr, " protecting this model\n");
#endif
	gretl_model_protect(pmod);
    }

    return pmod;
}

/**
 * gretl_model_new:
 * 
 * Allocates memory for a gretl #MODEL struct and initializes the struct,
 * using gretl_model_init().
 *
 * Returns: pointer to model (or %NULL if allocation fails).
 */

MODEL *gretl_model_new (void)
{
    return real_gretl_model_new(0);
}

/**
 * gretl_model_new_protected:
 * 
 * Allocates memory for a gretl #MODEL struct and initializes the struct,
 * using gretl_model_init().  Unlike plain gretl_model_new(), this function
 * creates a model which will not be freed by the mechanisms in objstack.c.
 *
 * Returns: pointer to model (or %NULL if allocation fails).
 */

MODEL *gretl_model_new_protected (void)
{
    return real_gretl_model_new(1);
}

/**
 * gretl_model_array_new:
 * @n: number of models in array.
 * 
 * Allocates memory for an array of @n gretl #MODEL structs and 
 * initializes each model using gretl_model_init().
 *
 * Returns: pointer to models array (or %NULL if allocation fails).
 */

MODEL **gretl_model_array_new (int n)
{
    MODEL **models;
    int i, j;

    models = malloc(n * sizeof *models);

    if (models != NULL) {
	for (i=0; i<n; i++) {
	    models[i] = gretl_model_new();
	    if (models[i] == NULL) {
		for (j=0; j<i; j++) {
		    free(models[i]);
		}
		free(models);
		models = NULL;
		break;
	    }
	}
    } 

    return models;
}

/**
 * gretl_model_array_destroy:
 * @models: array of gretl models.
 * @n: number of models in array.
 * 
 * Frees all resources associated with an array of models, which
 * should have been obtained via gretl_model_array_new().
 */

void gretl_model_array_destroy (MODEL **models, int n)
{
    int i;

    if (models != NULL) {
	for (i=0; i<n; i++) {
	    clear_model(models[i]);
	    free(models[i]);
	}
	free(models);
    }
}

static void clear_ar_info (MODEL *pmod)
{
    if (pmod->arinfo->arlist) {
	free(pmod->arinfo->arlist);
    }
    if (pmod->arinfo->rho) {
	free(pmod->arinfo->rho);
    }
    if (pmod->arinfo->sderr) {
	free(pmod->arinfo->sderr);
    }

    free(pmod->arinfo);
    pmod->arinfo = NULL;
}

#if MDEBUG > 1

static void 
debug_print_model_info (const MODEL *pmod, const char *msg)
{
    fprintf(stderr, "%s:\n"
	    " pmod = %p\n"
	    " pmod->list = %p\n"
	    " pmod->submask = %p\n"
	    " pmod->missmask = %p\n"
	    " pmod->coeff = %p\n"
	    " pmod->sderr = %p\n"
	    " pmod->yhat = %p\n"
	    " pmod->uhat = %p\n"
	    " pmod->xpx = %p\n"
	    " pmod->vcv = %p\n"
	    " pmod->name = %p\n"
	    " pmod->params = %p\n"
	    " pmod->arinfo = %p\n"
	    " pmod->tests = %p\n"
	    " pmod->data = %p\n", msg,
	    (void *) pmod, (void *) pmod->list, 
	    (void *) pmod->submask, (void *) pmod->missmask, 
	    (void *) pmod->coeff, (void *) pmod->sderr, 
	    (void *) pmod->yhat, (void *) pmod->uhat, 
	    (void *) pmod->xpx, (void *) pmod->vcv, 
	    (void *) pmod->name, (void *) pmod->params, 
	    (void *) pmod->arinfo, (void *) pmod->tests, 
	    (void *) pmod->data);
}

#endif /* MDEBUG */

/**
 * gretl_model_destroy_tests:
 * @pmod: pointer to model.
 *
 * Clears any hypothesis test structs that have been attached
 * to @pmod.
 */

void gretl_model_destroy_tests (MODEL *pmod)
{
    if (pmod != NULL && pmod->ntests > 0) {
	int i;

	for (i=0; i<pmod->ntests; i++) {
	    if (pmod->tests[i].param != NULL) {
		free(pmod->tests[i].param);
	    }
	}
	free(pmod->tests);
	pmod->tests = NULL;
	pmod->ntests = 0;
    }
}

/**
 * clear_model:
 * @pmod: pointer to model.
 *
 * Clears a gretl #MODEL, freeing all allocated storage and setting
 * pointer members to %NULL.  Also frees any data pointers attached
 * via gretl_model_set_data().  The model pointer itself is not
 * freed, so this function may be called on a #MODEL which has been 
 * declared directly by the caller (by passing the address of the
 * #MODEL).  
 */

void clear_model (MODEL *pmod)
{
    if (pmod != NULL) {
	int i;

#if MDEBUG
	fprintf(stderr, "clear model: model at %p\n", (void *) pmod);
#endif

#if MDEBUG > 1
	debug_print_model_info(pmod, "Doing clear_model");
#endif
	if (pmod->list) free(pmod->list);
	if (pmod->submask) free(pmod->submask);
	if (pmod->missmask) free(pmod->missmask);
	if (pmod->coeff) free(pmod->coeff);
	if (pmod->sderr) free(pmod->sderr);
	if (pmod->yhat) free(pmod->yhat);
	if (pmod->uhat) free(pmod->uhat);
	if (pmod->xpx) free(pmod->xpx);
	if (pmod->vcv) free(pmod->vcv);
	if (pmod->name) free(pmod->name);
	if (pmod->arinfo) {
	    clear_ar_info(pmod);
	}
	if (pmod->params) {
	    for (i=0; i<pmod->nparams; i++) {
		free(pmod->params[i]);
	    }
	    free(pmod->params);
	}
	if (pmod->dataset) {
	    free_model_dataset(pmod);
	}
	gretl_model_destroy_tests(pmod);
	destroy_all_data_items(pmod);
    }

    /* this may be redundant */
#if MDEBUG
    fprintf(stderr, "clear_model, calling gretl_model_init\n");
#endif
#if 1
    gretl_model_init(pmod);
#endif
}

/**
 * gretl_model_free:
 * @pmod: pointer to #MODEL.
 *
 * Free allocated content of @pmod then the pointer itself,
 * if the model's reference count has reached zero.
 */

void gretl_model_free (MODEL *pmod)
{
    if (pmod != NULL) {
#if MDEBUG
	fprintf(stderr, "gretl_model_free: pmod at %p, incoming refcount = %d\n",
		(void *) pmod, pmod->refcount);
#endif
	pmod->refcount -= 1;
	if (pmod->refcount <= 0) {
	    clear_model(pmod);
	    free(pmod);
	}
    }
}

/**
 * gretl_model_free_on_exit:
 * @pmod: pointer to #MODEL.
 *
 * Free allocated content of @pmod then the pointer itself,
 * without regard to the model's reference count.
 */

void gretl_model_free_on_exit (MODEL *pmod)
{
    if (pmod != NULL) {
	remove_model_from_stack(pmod);
	clear_model(pmod);
	free(pmod);
    }
}

static void copy_test (ModelTest *targ, const ModelTest *src)
{
    targ->type = src->type;
    
    if (src->param != NULL) {
	targ->param = gretl_strdup(src->param);
    } else {
	targ->param = NULL;
    }

    targ->teststat = src->teststat;
    targ->dfn = src->dfn;
    targ->dfd = src->dfd;
    targ->order = src->order;
    targ->value = src->value;
    targ->pvalue = src->pvalue;
}

static int copy_model_tests (MODEL *targ, const MODEL *src)
{
    int i, n = src->ntests;

    if (n <= 0 || src->tests == NULL) {
	return 0;
    }

    targ->tests = malloc(n * sizeof *targ->tests);
    if (targ->tests == NULL) return 1;

    for (i=0; i<n; i++) {
	copy_test(&targ->tests[i], &src->tests[i]);
    }

    return 0;
}

static void gretl_test_init (ModelTest *test, ModelTestType ttype)
{
    test->type = ttype;
    test->order = 0;
    test->param = NULL;
    test->teststat = GRETL_STAT_NONE;
    test->dfn = test->dfd = 0;
    test->value = test->pvalue = NADBL;
}

static int model_tests_differ (ModelTest *mt1, ModelTest *mt2)
{
    int ret = 0;

    if (mt1->type != mt2->type) {
	ret = 1;
    } else if (mt1->order != mt2->order) {
	ret = 1;
    } else if (mt1->param != NULL && mt2->param != NULL &&
	       strcmp(mt1->param, mt2->param)) {
	ret = 1;
    } else if (mt1->teststat != mt2->teststat) {
	ret = 1;
    } else if (mt1->value != mt2->value) {
	ret = 1;
    }

    return ret;
}

/**
 * model_test_new:
 * @ttype: type of test to add.
 *
 * Returns: new #ModelTest pointer, or %NULL on failure.
 */

ModelTest *model_test_new (ModelTestType ttype)
{
    ModelTest *test = malloc(sizeof *test);

    if (test != NULL) {
	gretl_test_init(test, ttype);
    }

    return test;
}

/**
 * maybe_add_test_to_model:
 * @pmod: pointer to model.
 * @test: model test to be added.
 *
 * Adds a #ModelTest to @pmod, if the test in question has
 * not already been performed and recorded.
 *
 * Returns: 1 if the test was added, otherwise 0.
 */

int maybe_add_test_to_model (MODEL *pmod, ModelTest *test)
{
    ModelTest *tests = NULL;
    int i, nt = pmod->ntests;
    int done = 0, add = 0;

    if (test == NULL || test->teststat == GRETL_STAT_NONE) {
	return 0;
    }

    for (i=0; i<nt; i++) {
	if (!model_tests_differ(test, &pmod->tests[i])) {
	    done = 1;
	}
    } 

    if (!done) {
	tests = realloc(pmod->tests, (nt + 1) * sizeof *tests);
    }

    if (tests != NULL) {
	pmod->tests = tests;
	pmod->ntests += 1;
	copy_test(&pmod->tests[nt], test);
	add = 1;
    }

    free(test->param);
    free(test);

    return add;
}

void model_test_set_teststat (ModelTest *test, unsigned char ts)
{
    test->teststat = ts;
}

void model_test_set_order (ModelTest *test, int order)
{
    test->order = order;
}

void model_test_set_dfn (ModelTest *test, int df)
{
    test->dfn = df;
}

void model_test_set_dfd (ModelTest *test, int df)
{
    test->dfd = df;
}

void model_test_set_value (ModelTest *test, double val)
{
    test->value = val;
}

void model_test_set_pvalue (ModelTest *test, double pval)
{
    test->pvalue = pval;
}

void model_test_set_param (ModelTest *test, const char *s)
{
    test->param = gretl_strdup(s);
}

void model_test_set_allocated_param (ModelTest *test, char *s)
{
    test->param = s;
}

static int gretl_test_print_string (const ModelTest *test, PRN *prn)
{
    const char *test_strs[] = {
	N_("Test for addition of variables"),
	N_("Test for ARCH of order %s"),
	N_("LM test for autocorrelation up to order %s"),
	N_("Chow test for structural break at observation %s"),
	N_("CUSUM test for parameter stability"),
	N_("Likelihood ratio test for groupwise heteroskedasticity"),
	N_("Non-linearity test (logs)"),
	N_("Test for normality of residual"),
	N_("Test for omission of variables"),
	N_("RESET test for specification"),
	N_("Non-linearity test (squares)"),
	N_("White's test for heteroskedasticity"),
	N_("Sargan over-identification test"),
	N_("Hausman test")
    };
    char ordstr[16];
    char *param = NULL;

    if (test->type >= GRETL_TEST_MAX) {
	return 1;
    }

    if (test->order > 0) {
	sprintf(ordstr, "%d", test->order);
	param = ordstr;
    } else if (test->type == GRETL_TEST_CHOW) {
	param = test->param;
    }

    if (param != NULL) {
	if (plain_format(prn)) {
	    pprintf(prn, _(test_strs[test->type]), param);
	} else {
	    pprintf(prn, I_(test_strs[test->type]), param);
	}
    } else {
	if (plain_format(prn)) {
	    pputs(prn, _(test_strs[test->type]));
	} else {
	    pputs(prn, I_(test_strs[test->type]));
	}
    }

    return 0;
}

static int print_add_omit_varnames (const char *s, PRN *prn)
{
    const char *endings[] = {
	"\n    ",
	"\\\\\n \\qquad ",
	"\\par\n    "
    };
    const char *ending = NULL;

    if (s == NULL || *s == '\0') {
	return 1;
    }

    if (plain_format(prn)) {
	ending = endings[0];
    } else if (tex_format(prn)) {
	ending = endings[1];
    } else if (rtf_format(prn)) {
	ending = endings[2];
    } else {
	return 1;
    }

    pputs(prn, ending);

    while (*s) {
	if (*s == ' ') {
	    pputs(prn, ending);
	} else {
	    pputc(prn, *s);
	}
	s++;
    }

    return 0;
}

static int gretl_test_print_h_0 (const ModelTest *test, PRN *prn)
{
    const char *h_0_strs[] = {
	N_("parameters are zero for the variables"),
	N_("no ARCH effect is present"),
	N_("no autocorrelation"),
	N_("no structural break"),
	N_("no change in parameters"),
	N_("the units have a common error variance"),
	N_("relationship is linear"),
	N_("error is normally distributed"),
	N_("parameters are zero for the variables"),
	N_("specification is adequate"),
	N_("relationship is linear"),
	N_("heteroskedasticity not present"),
	N_("all instruments are valid"),
	N_("OLS estimates are consistent")
    };

    if (test->type >= GRETL_TEST_MAX) {
	return 1;
    }  

    if (plain_format(prn)) {
	pputs(prn, _(h_0_strs[test->type]));
    } else {
	pputs(prn, I_(h_0_strs[test->type]));
    }  

    if (test->type == GRETL_TEST_ADD || test->type == GRETL_TEST_OMIT) {
	print_add_omit_varnames(test->param, prn);
    } 

    return 0;
}

static void 
get_test_stat_string (const ModelTest *test, char *str, PRN *prn)
{
    int tex = tex_format(prn);

    switch (test->teststat) {
    case GRETL_STAT_TR2:
	if (tex) {
	    sprintf(str, "$TR^2$ = %g", test->value);
	} else if (rtf_format(prn)) {
	    sprintf(str, "TR{\\super 2} = %g", test->value);
	} else {
	    sprintf(str, "TR^2 = %g", test->value);
	}
	break;
    case GRETL_STAT_F:
    case GRETL_STAT_RESET:
	if (tex) {
	    sprintf(str, "$F(%d, %d)$ = %g", test->dfn, test->dfd, test->value);
	} else {
	    sprintf(str, "F(%d, %d) = %g", test->dfn, test->dfd, test->value);
	}
	break;
    case GRETL_STAT_LMF:
	sprintf(str, "LMF = %g", test->value);
	break;
    case GRETL_STAT_HARVEY_COLLIER:
	if (tex) {
	    sprintf(str, "Harvey--Collier $t(%d)$ = %g", test->dfn, test->value);
	} else {
	    sprintf(str, "Harvey-Collier t(%d) = %g", test->dfn, test->value);
	}
	break;
    case GRETL_STAT_NORMAL_CHISQ:
	if (tex) {
	    sprintf(str, "$\\chi^2_2$ = %g", test->value); 
	} else {
	    sprintf(str, "%s(2) = %g", _("Chi-square"), test->value);
	}
	break;
    case GRETL_STAT_LR:
    case GRETL_STAT_WALD_CHISQ:
	if (tex) {
	    sprintf(str, "$\\chi^2_%d$ = %g", test->dfn, test->value); 
	} else {
	    sprintf(str, "%s(%d) = %g", _("Chi-square"), test->dfn, test->value);
	}
	break;
    default:
	*str = 0;
    }
}

static void 
get_test_pval_string (const ModelTest *test, char *str, PRN *prn)
{
    int tex = tex_format(prn);

    switch (test->teststat) {
    case GRETL_STAT_TR2:
	if (tex) sprintf(str, "$P$($\\chi^2_{%d} >$ %g) = %g", 
			 test->dfn, test->value, test->pvalue);
	else sprintf(str, "P(Chi-Square(%d) > %g) = %g", 
		     test->dfn, test->value, test->pvalue);
	break;
    case GRETL_STAT_F:
    case GRETL_STAT_RESET:
	if (tex) {
	    sprintf(str, "$P$($F(%d, %d) >$ %g) = %g", 
		    test->dfn, test->dfd, test->value, test->pvalue);
	} else {
	    sprintf(str, "P(F(%d, %d) > %g) = %g", 
		    test->dfn, test->dfd, test->value, test->pvalue);
	}
	break;
    case GRETL_STAT_LMF:
	if (tex) {
	    sprintf(str, "$P$($F(%d, %d) >$ %g) = %g", 
		    test->dfn, test->dfd, test->value, test->pvalue);
	} else {
	    sprintf(str, "P(F(%d,%d) > %g) = %g", 
		    test->dfn, test->dfd, test->value, test->pvalue);
	}
	break;
    case GRETL_STAT_HARVEY_COLLIER:
	if (tex) {
	    sprintf(str, "$P$($t_{%d} >$ %g)  = %g", 
		    test->dfn, test->value, test->pvalue);
	} else {
	    sprintf(str, "P(t(%d) > %g) = %g", 
		    test->dfn, test->value, test->pvalue);
	}
	break;
    case GRETL_STAT_NORMAL_CHISQ:
    case GRETL_STAT_LR:
    case GRETL_STAT_WALD_CHISQ:
	sprintf(str, "%g", test->pvalue);
	break;
    default:
	*str = 0;
    }
}

void gretl_model_test_print (const MODEL *pmod, int i, PRN *prn)
{
    const ModelTest *test;
    char test_str[128], pval_str[128];
    const char *tstat;

    if (i >= pmod->ntests) {
	return;
    }

    test = &pmod->tests[i];

    if (test->teststat == GRETL_STAT_WALD_CHISQ) {
	tstat = N_("Asymptotic test statistic");
    } else {
	tstat = N_("Test statistic");
    }

    get_test_stat_string(test, test_str, prn);
    get_test_pval_string(test, pval_str, prn);

    if (plain_format(prn)) {
	gretl_test_print_string(test, prn);
	pprintf(prn, " -\n  %s: ", _("Null hypothesis"));
	gretl_test_print_h_0(test, prn);
	pprintf(prn, "\n  %s: %s\n"
		"  %s = %s\n\n",
		_(tstat), test_str, 
		_("with p-value"), pval_str);
    } else if (tex_format(prn)) {
	gretl_test_print_string(test, prn);
	pprintf(prn, " --\\\\\n\\quad %s: ", I_("Null hypothesis"));
	gretl_test_print_h_0(test, prn);
	pprintf(prn, "\\\\\n\\quad %s: %s\\\\\n"
		"\\quad %s = %s\\\\\n",
		I_(tstat), test_str, 
		I_("with p-value"), pval_str);
    } else if (rtf_format(prn)) {
	pputs(prn, "\\par \\ql ");
	gretl_test_print_string(test, prn);
	pprintf(prn, " -\\par\n %s: ", I_("Null hypothesis"));
	gretl_test_print_h_0(test, prn);
	pprintf(prn, "\\par\n %s: %s\\par\n"
		" %s = %s\\par\n\n",
		I_(tstat), test_str, 
		I_("with p-value"), pval_str);
    }
}

void gretl_model_print_last_test (const MODEL *pmod, PRN *prn)
{
    gretl_model_test_print(pmod, pmod->ntests - 1, prn);
}

static ARINFO *copy_ar_info (const ARINFO *src)
{
    ARINFO *targ;

    targ = malloc(sizeof *targ);
    if (targ == NULL) return NULL;

    if (src->arlist != NULL) {
	int i, m = src->arlist[0];

      	if ((targ->rho = copyvec(src->rho, m)) == NULL) {
	    free(targ);
	    targ = NULL;
	    return NULL; 
	}

      	if ((targ->sderr = copyvec(src->sderr, m)) == NULL) { 
	    free(targ->rho);
	    free(targ);
	    targ = NULL;
	    return NULL; 
	}

	targ->arlist = malloc((m + 1) * sizeof(int));
	if (targ->arlist == NULL) {
	    free(targ->rho);
	    free(targ->sderr);
	    free(targ);
	    targ = NULL;
	    return NULL; 
	}
	for (i=0; i<=m; i++) {
	    targ->arlist[i] = src->arlist[i];
	}
    }    

    return targ;
}

static int copy_model_params (MODEL *targ, const MODEL *src)
{
    int i, j, n = src->nparams;

    targ->params = malloc(n * sizeof *targ->params);
    if (targ->params == NULL) return 1;

    targ->nparams = n;

    for (i=0; i<n; i++) {
	targ->params[i] = gretl_strdup(src->params[i]);
	if (targ->params[i] == NULL) {
	    for (j=0; j<i; j++) {
		free(targ->params[j]);
	    }
	    free(targ->params);
	    targ->params = NULL;
	    targ->nparams = 0;
	    return 1;
	}
    }

    return 0;
}

static int copy_model_data_items (MODEL *targ, const MODEL *src)
{
    int i, n = src->n_data_items;
    int err = 0;

    targ->data_items = malloc(n * sizeof *targ->data_items);
    if (targ->data_items == NULL) return 1;

    for (i=0; i<src->n_data_items; i++) {
	targ->data_items[i] = NULL;
    }

    for (i=0; i<src->n_data_items; i++) {
	model_data_item *targitem, *srcitem;

	targitem = malloc(sizeof *targitem);
	if (targitem == NULL) {
	    err = 1;
	    break;
	}

	srcitem = src->data_items[i];
	targ->data_items[i] = targitem;
	
	targitem->key = gretl_strdup(srcitem->key);
	if (targitem->key == NULL) {
	    err = 1;
	    break;
	}

	targitem->ptr = malloc(srcitem->size);
	if (targitem->ptr == NULL) {
	    free(targitem->key);
	    err = 1;
	    break;
	}

	memcpy(targitem->ptr, srcitem->ptr, srcitem->size);
	targitem->size = srcitem->size;
	targitem->destructor = srcitem->destructor;
    }

    if (err) {
	for (i=0; i<targ->n_data_items; i++) {
	    free(targ->data_items[i]);
	}
	free(targ->data_items);
	targ->data_items = NULL;
	targ->n_data_items = 0;
    }

    targ->n_data_items = src->n_data_items;

    return err;
}

static char *copy_missmask (const MODEL *pmod)
{
    char *mask;
    int n = pmod->t2 - pmod->t1 + 1;

    mask = malloc(n);
    if (mask == NULL) {
	return NULL;
    }

    memcpy(mask, pmod->missmask, n);

    return mask;
}

static int copy_model (MODEL *targ, const MODEL *src)
{
    int k = src->ncoeff;
    int m = k * (k + 1) / 2;

    /* monolithic copy of structure */
    *targ = *src;

    /* temporarily zero various count members just in case
       copying of allocated info fails */
    targ->ntests = 0;
    targ->nparams = 0;
    targ->n_data_items = 0;

    /* now work on pointer members */
    gretl_model_init_pointers(targ);

    if (src->coeff != NULL &&
	(targ->coeff = copyvec(src->coeff, src->ncoeff)) == NULL) {
	return 1;
    }

    if (src->sderr != NULL &&
	(targ->sderr = copyvec(src->sderr, src->ncoeff)) == NULL) {  
	return 1;
    }

    if (src->uhat != NULL && 
	(targ->uhat = copyvec(src->uhat, src->full_n)) == NULL) {
	return 1;
    }

    if (src->yhat != NULL && 
	(targ->yhat = copyvec(src->yhat, src->full_n)) == NULL) {
	return 1;
    }

    if (src->submask != NULL && 
	(targ->submask = copy_subsample_mask(src->submask)) == NULL) { 
	return 1;
    }

    if (src->missmask != NULL && 
	(targ->missmask = copy_missmask(src)) == NULL) { 
	return 1;
    }

    if (src->xpx != NULL && (targ->xpx = copyvec(src->xpx, m)) == NULL) {
	return 1;
    }

    if (src->vcv != NULL && (targ->vcv = copyvec(src->vcv, m)) == NULL) {
	return 1;
    }

    if (src->arinfo != NULL && 
	(targ->arinfo = copy_ar_info(src->arinfo)) == NULL) {
	return 1; 
    }

    if (src->ntests > 0 && src->tests != NULL) {
	copy_model_tests(targ, src);
	if (targ->tests == NULL) {
	    return 1;
	}
	targ->ntests = src->ntests;
    }

    if (src->nparams > 0 && src->params != NULL) {
	copy_model_params(targ, src);
	if (targ->params == NULL) {
	    return 1;
	}
	targ->nparams = src->nparams;
    }    

    if (src->n_data_items > 0) {
	copy_model_data_items(targ, src);
	if (targ->data_items == NULL) {
	    return 1;
	}
	targ->n_data_items = src->n_data_items;
    }

    if (src->list != NULL && 
	(targ->list = gretl_list_copy(src->list)) == NULL) {
	return 1;
    } 

    return 0;
}

/**
 * gretl_model_copy:
 * @pmod: pointer to #MODEL to copy.
 *
 * Does a deep copy of @pmod: allocates a new #MODEL pointer
 * which has its own allocated copies of all the pointer
 * members of @pmod.  The only feature of @pmod that is
 * not duplicated is the reference count, which is set
 * to zero in the copy. 
 *
 * Returns: the copied model, or %NULL on failure.
 */

MODEL *gretl_model_copy (const MODEL *pmod)
{
    MODEL *new = malloc(sizeof *new);

    if (new != NULL) {
	int err = copy_model(new, pmod);

	if (err) {
	    clear_model(new);
	    free(new);
	    new = NULL;
	} else {
	    /* nota bene */
	    new->refcount = 0;
	}
    }

    return new;
}

/**
 * swap_models:
 * @targ: pointer to pointer to #MODEL.
 * @src: pointer to pointer to #MODEL.
 *
 * Swaps the model pointers.
 *
 * Returns: 0 on success.
 */

int swap_models (MODEL **targ, MODEL **src)
{
    MODEL *tmp = *targ;

#if MDEBUG
    fprintf(stderr, "swap_models: %p <-> %p\n", *targ, *src);
#endif

    *targ = *src;
    *src = tmp;

    /* handle the case where we have a model that is currently
       defined as the "last model" (objstack.c) */
    maybe_swap_into_last_model(*targ, tmp);

    return 0;
}

int is_model_cmd (const char *s)
{
    int ret = 0;

    /* FIXME mle? */

    if (s != NULL && *s != '\0') {
	if (!strncmp(s, "ols", 3)  ||
	    !strncmp(s, "corc", 4) ||
	    !strncmp(s, "hilu", 4) ||
	    !strncmp(s, "wls", 3)  ||
	    !strncmp(s, "pwe", 3)  ||
	    !strncmp(s, "pooled", 6)  ||
	    !strncmp(s, "hccm", 4) ||
	    !strncmp(s, "hsk", 3)  ||
	    !strncmp(s, "add ", 4)  ||
	    !strncmp(s, "lad", 3)  ||
	    !strncmp(s, "omit", 4) ||
	    !strncmp(s, "tsls", 4) ||
	    !strncmp(s, "logit", 5)  ||
	    !strncmp(s, "probit", 6) ||
	    !strncmp(s, "tobit", 5) ||
	    !strncmp(s, "poisson", 7) ||
	    !strncmp(s, "garch", 5) ||
	    !strncmp(s, "logistic", 8) ||
	    !strncmp(s, "end nls", 7) ||
	    !strncmp(s, "arma", 4) ||
	    !strncmp(s, "arima", 5) ||
	    !strncmp(s, "ar ", 3) ||
	    !strcmp(s, "ar")) {
	    ret = 1;
	}
    }

    return ret;
}

int is_quiet_model_test (int ci, gretlopt opt)
{
    int ret = 0;

    if ((opt & OPT_Q) && (ci == OMIT || ci == ADD ||
			  ci == OMITFROM || ci == ADDTO)) {
	ret = 1;
    }

    return ret;
}

/**
 * command_ok_for_model:
 * @test_ci:  index of command to be tested.
 * @model_ci: command index of a gretl model (for example,
 * %OLS, %HCCM or %CORC).
 *
 * Returns: 1 if the model-related command in question is
 * meaningful and acceptable in the context of the specific
 * sort of model indentified by @model_ci, otherwise 0.
 */

int command_ok_for_model (int test_ci, int model_ci)
{
    int ok = 1;

    if (model_ci == MLE) {
	return 0;
    }

    switch (test_ci) {
    case ADD:
    case ADDTO:
    case OMIT:
    case OMITFROM:
	if (model_ci == NLS || model_ci == ARMA || 
	    model_ci == GARCH) {
	    ok = 0;
	}
	break;

    case COEFFSUM:
    case VIF:
	if (model_ci == NLS || model_ci == TSLS ||
	    model_ci == ARMA || model_ci == GARCH) {
	    ok = 0;
	}
	break;

    case EQNPRINT:
	if (model_ci == ARMA || model_ci == NLS) {
	    ok = 0; 
	}
	break;

    case LMTEST:
	if (model_ci != OLS && model_ci != POOLED) ok = 0;
	break;

    case ARCH:
    case CHOW:
    case CUSUM:
    case LEVERAGE:
    case RESET:
	if (model_ci != OLS) ok = 0;
	break;

    case HAUSMAN:
	if (model_ci != POOLED) ok = 0;
	break;

    case RESTRICT:
	if (model_ci == LAD || model_ci == NLS) ok = 0;
	break;

    case TESTUHAT:
	/* do we really need to exclude garch? */
	if (model_ci == TOBIT || model_ci == GARCH) ok = 0;
	break;

    default:
	break;
    }

    return ok;
}

static int gretl_model_count;

int get_model_count (void)
{
    return gretl_model_count;
}

void reset_model_count (void)
{
    gretl_model_count = 0;
}

int model_count_plus (void)
{
    return ++gretl_model_count;
}

void model_count_minus (void)
{
    --gretl_model_count;
}

void set_model_id (MODEL *pmod)
{
    if (pmod->errcode == 0) {
	pmod->ID = ++gretl_model_count;
    }
}

void model_list_to_string (int *list, char *buf)
{
    char numstr[8];
    int i;

    for (i=1; i<=list[0]; i++) {
	if (list[i] == LISTSEP) {
	    strcat(buf, "; ");
	} else {
	    sprintf(numstr, "%d ", list[i]);
	    strcat(buf, numstr);
	}
    }
}

int highest_numbered_var_in_model (const MODEL *pmod, 
				   const DATAINFO *pdinfo)
{
    int i, v, vmax = 0;
    int gotsep = 0;

    if (pmod->ci == MLE) {
	return 0;
    }

    for (i=1; i<=pmod->list[0]; i++) {
	v = pmod->list[i];
	if (v == LISTSEP) {
	    gotsep = 1;
	    continue;
	}
	if (v >= pdinfo->v) {
	    /* temporary variables, already gone? */
	    continue;
	}
	if ((pmod->ci == ARMA || pmod->ci == GARCH) && !gotsep) {
	    /* real vars start after LISTSEP */
	    continue;
	}
#if 0
	fprintf(stderr, "highest numbered... checking var %d\n", v);
#endif
	if (v > vmax) {
	    vmax = v;
	}
	if (pmod->ci == NLS) {
	    /* only the dependent var can be tested */
	    break;
	}
    }

    return vmax;
}

int mle_criteria (MODEL *pmod, int addk)
{
    int err = 0;

    if (na(pmod->lnL)) {
	pmod->criterion[C_AIC] = NADBL;
	pmod->criterion[C_BIC] = NADBL;
	pmod->criterion[C_HQC] = NADBL;
	err = 1;
    } else {
	int k = pmod->ncoeff + addk;

	pmod->criterion[C_AIC] = -2.0 * pmod->lnL + 2.0 * k;
	pmod->criterion[C_BIC] = -2.0 * pmod->lnL + k * log(pmod->nobs);
	pmod->criterion[C_HQC] = -2.0 * pmod->lnL + 2 * k * log(log(pmod->nobs));
    }

    return err;
}

/* As of 2005-02-19, not ready for this yet: need to adjust
   the text in the model printout.
*/

double coeff_pval (const MODEL *pmod, double x, int df)
{
    if (pmod->ci == MLE) {
	return normal_pvalue_2(x);
    } else if (0 && ML_ESTIMATOR(pmod->ci)) {
	/* not ready */
	return normal_pvalue_2(x);
    } else {
	return t_pvalue_2(x, df);
    }
}

/**
 * gretl_model_add_arma_varnames:
 * @pmod: pointer to target model.
 * @pdinfo: dataset information.
 * @yno: ID number of dependent variable.
 * @p: non-seasonal AR order.
 * @q: non-seasonal MA order.
 * @P: seasonal AR order.
 * @Q: seasonal MA order.
 * @r: number of exogenous regressors (excluding the constant).
 * 
 * Composes a set of names to be given to the regressors in an
 * ARMA model.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_model_add_arma_varnames (MODEL *pmod, const DATAINFO *pdinfo,
				   int yno, int p, int q, int P, int Q, 
				   int r)
{
    int np = p + P + q + Q + r + pmod->ifc + 1;
    int xstart;
    int i, j;

    pmod->params = malloc(np * sizeof pmod->params);
    if (pmod->params == NULL) {
	pmod->errcode = E_ALLOC;
	return 1;
    }

    pmod->nparams = np;

    for (i=0; i<np; i++) {
	pmod->params[i] = malloc(VNAMELEN);
	if (pmod->params[i] == NULL) {
	    for (j=0; j<i; j++) {
		free(pmod->params[j]);
	    }
	    free(pmod->params);
	    pmod->params = NULL;
	    pmod->nparams = 0;
	    pmod->errcode = E_ALLOC;
	    return 1;
	}
    }

    strcpy(pmod->params[0], pdinfo->varname[yno]);

    if (pmod->ifc) {
	strcpy(pmod->params[1], pdinfo->varname[0]);
	j = 2;
    } else {
	j = 1;
    }

    for (i=0; i<p; i++) {
	sprintf(pmod->params[j++], "phi_%d", i + 1);
    }

    for (i=0; i<P; i++) {
	sprintf(pmod->params[j++], "Phi_%d", i + 1);
    }   

    for (i=0; i<q; i++) {
	sprintf(pmod->params[j++], "theta_%d", i + 1);
    }

    for (i=0; i<Q; i++) {
	sprintf(pmod->params[j++], "Theta_%d", i + 1);
    }       

    xstart = (P || Q)? 8 : 5;

    for (i=0; i<r; i++) {
	strcpy(pmod->params[j++], pdinfo->varname[pmod->list[xstart+i]]); 
    }   

    return 0;
}

/* try to tell if an OLS model with two independent variables is
   actually a quadratic model (x_2 = x_1^2). */

static int model_is_quadratic (const MODEL *pmod, const double **Z,
			       const DATAINFO *pdinfo)
{
    const double *x1, *x2;
    int t, ret = 1;

    x1 = Z[pmod->list[3]];
    x2 = Z[pmod->list[4]];

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(x1[t]) && x2[t] != x1[t] * x1[t]) {
	    ret = 0;
	    break;
	}
    }

    return ret;
}

/**
 * gretl_model_get_fitted_formula:
 * @pmod: pointer to target model.
 * @xvar: ID number of variable that _may_ be "x" in the model.
 * @Z: data array.
 * @pdinfo: dataset information.
 * 
 * If @pmod is a simple linear, quadratic or logistic model, 
 * and if @xvar is in fact the "x" variable from the model, 
 * returns a string representing the formula for generating the 
 * fitted values as a function of x.  This formula may be used 
 * in the context of a fitted versus actual plot.
 *
 * Returns: formula for fitted values, or %NULL if this is
 * not available.
 */

char *gretl_model_get_fitted_formula (const MODEL *pmod, int xvar, 
				      const double **Z,
				      const DATAINFO *pdinfo)
{
    const double **mZ;
    const DATAINFO *mdinfo;
    char *ret = NULL;

    /* only OLS and logistic are handled */
    if (xvar == 0 || (pmod->ci != OLS && pmod->ci != LOGISTIC)) {
	return NULL;
    }

    if (pmod->dataset != NULL) {
	mZ = (const double **) pmod->dataset->Z;
	mdinfo = pmod->dataset->dinfo;
    } else {
	mZ = Z;
	mdinfo = pdinfo;
    }

    gretl_push_c_numeric_locale();

    if (pmod->ci == LOGISTIC) {
	if (pmod->ifc && pmod->ncoeff == 2 && xvar == pmod->list[3]) {
	    double lmax = gretl_model_get_double(pmod, "lmax");

	    if (!na(lmax)) {
		ret = g_strdup_printf("yformula: %g/(1.0+exp(-(%g%s%g*x)))",
				      lmax, pmod->coeff[0], 
				      (pmod->coeff[1] >= 0)? "+" : "",
				      pmod->coeff[1]);
	    }
	}
    } else if (!pmod->ifc && pmod->ncoeff == 1 && xvar == pmod->list[2]) {
	ret = g_strdup_printf("yformula: %g*x", pmod->coeff[0]);
    } else if (pmod->ifc && pmod->ncoeff == 2 && xvar == pmod->list[3]) {
	ret = g_strdup_printf("yformula: %g%s%g*x", pmod->coeff[0], 
			      (pmod->coeff[1] >= 0)? "+" : "",
			      pmod->coeff[1]);
    } else if (pmod->ifc && pmod->ncoeff == 3 && xvar == pmod->list[3]) {
	if (model_is_quadratic(pmod, mZ, mdinfo)) {
	    ret = g_strdup_printf("yformula: %g%s%g*x%s%g*x**2", pmod->coeff[0], 
				  (pmod->coeff[1] >= 0)? "+" : "",
				  pmod->coeff[1], 
				  (pmod->coeff[2] >= 0)? "+" : "",
				  pmod->coeff[2]);
	}
    }

    gretl_pop_c_numeric_locale();
	
    return ret;
}

void gretl_model_set_name (MODEL *pmod, const char *name)
{
    if (pmod->name != NULL) {
	free(pmod->name);
    }

    pmod->name = gretl_strdup(name);
}

const char *gretl_model_get_name (const MODEL *pmod)
{
    if (pmod != NULL) {
	return pmod->name;
    }

    return NULL;
}

int gretl_model_data_index (const char *s)
{
    char test[8] = {0};
    int ssel = 0;
    int msel = 0;

    strncat(test, s, 7);
    lower(test);

    /* scalar values */
    if (!strcmp(test, "$ess"))  
	return M_ESS;
    if (!strcmp(test, "$t")) 
	return M_T;
    if (!strcmp(test, "$rsq"))  
	return M_RSQ;
    if (!strcmp(test, "$sigma"))  
	return M_SIGMA;
    if (!strcmp(test, "$df"))   
	return M_DF;
    if (!strcmp(test, "$ncoeff"))   
	return M_NCOEFF;
    if (!strcmp(test, "$lnl"))   
	return M_LNL;
    if (!strcmp(test, "$aic"))   
	return M_AIC;
    if (!strcmp(test, "$bic"))   
	return M_BIC;
    if (!strcmp(test, "$hqc"))   
	return M_HQC;
    if (!strcmp(test, "$nrsq") || 
	!strcmp(test, "$trsq")) 
	return M_TRSQ;

    sscanf(s, "%7[^[( ]", test);
    lower(test);

    if (strchr(s, '[')) {
	/* selecting submatrix? */
	msel = 1;
    }

    /* series or matrices */
    if (!strcmp(test, "$uhat")) 
	return M_UHAT;
    if (!strcmp(test, "$yhat"))
	return M_YHAT;
    if (!strcmp(test, "$h"))
	return M_H;

    if (!msel && strchr(s, '(')) {
	/* selecting scalar element from array? */
	ssel = 1;
    } 

    /* matrices, not series */
    if (!strcmp(test, "$coeff"))  
	return (ssel)? M_COEFF_S : M_COEFF;
    if (!strcmp(test, "$stderr"))  
	return (ssel)? M_SE_S : M_SE;
    if (!strcmp(test, "$vcv"))   
	return (ssel)? M_VCV_S : M_VCV;
    if (!strcmp(test, "$rho"))   
	return (ssel)? M_RHO_S : M_RHO;

    return 0;
}

double gretl_model_get_scalar (const MODEL *pmod, int idx, int *err)
{
    double x = NADBL;

    if (pmod == NULL) {
	*err = E_BADSTAT;
	return x;
    }

    switch (idx) {  
    case M_ESS:
	x = pmod->ess;
	break;
    case M_RSQ:
	x = pmod->rsq;
	break;
    case M_LNL:
	x = pmod->lnL;
	break;
    case M_AIC:
	x = pmod->criterion[C_AIC];
	break;
    case M_BIC:
	x = pmod->criterion[C_BIC];
	break;
    case M_HQC:
	x = pmod->criterion[C_HQC];
	break;
    case M_SIGMA:
	if (pmod->nwt) x = pmod->sigma_wt;
	else x = pmod->sigma;
	break;
    case M_TRSQ:
	if (!na(pmod->rsq)) {
	    x = pmod->nobs * pmod->rsq;
	}
	break;
    case M_DF:
	x = (double) pmod->dfd;
	break;
    case M_NCOEFF:
	x = (double) pmod->ncoeff; /* is this always available? */
	break;
    case M_T:
	x = (double) pmod->nobs;
	break;	
    default:
	break;
    }

    if (na(x)) {
	*err = E_BADSTAT;
    }

    return x;
}

double *
gretl_model_get_series (const MODEL *pmod, const DATAINFO *pdinfo, 
			int idx, int *err)
{
    double *x = NULL;
    double *garch_h = NULL;
    int t;

    if (pmod->t2 - pmod->t1 + 1 > pdinfo->n || 
	model_sample_issue(pmod, NULL, 0, pdinfo)) {
	strcpy(gretl_errmsg, 
	       (idx == M_UHAT)? 
	       _("Can't retrieve uhat: data set has changed") :
	       (idx == M_YHAT)?
	       _("Can't retrieve yhat: data set has changed") :
	       _("Can't retrieve ht: data set has changed"));
	*err = 1;
	return NULL;
    }   

    if ((idx == M_UHAT && pmod->uhat == NULL) ||
	(idx == M_YHAT && pmod->yhat == NULL)) {
	*err = 1;
	return NULL;
    }

    if (idx == M_H) {
	garch_h = gretl_model_get_data(pmod, "garch_h");
	if (garch_h == NULL) {
	    strcpy(gretl_errmsg, _("Can't retrieve error variance"));
	    *err = 1;
	    return NULL;
	}
    }

    x = malloc(pdinfo->n * sizeof *x);
    if (x == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (t=0; t<pdinfo->n; t++) {
	if (t < pmod->t1 || t > pmod->t2) {
	    x[t] = NADBL;
	} else {
	    if (idx == M_UHAT) {
		x[t] = pmod->uhat[t];
	    } else if (idx == M_YHAT) {
		x[t] = pmod->yhat[t];
	    } else if (idx == M_H) {
		x[t] = garch_h[t];
	    }
	}
    }
	    
    return x;
}

static gretl_matrix *
model_get_estvec (const MODEL *pmod, int idx, int *err)
{
    gretl_matrix *v = NULL;
    double x;
    int i;

    v = gretl_column_vector_alloc(pmod->ncoeff);
    if (v == NULL) {
	*err = E_ALLOC;
    } else {
	for (i=0; i<pmod->ncoeff; i++) {
	    x = (idx == M_COEFF)? pmod->coeff[i] : pmod->sderr[i];
	    gretl_vector_set(v, i, x);
	}
    }

    return v;
}

static gretl_matrix *
model_get_hatvec (const MODEL *pmod, int idx, int *err)
{
    gretl_matrix *v = NULL;
    double x;
    int t;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(pmod->uhat[t])) {
	    *err = E_MISSDATA;
	    break;
	}
    }

    if (!*err) {
	v = gretl_column_vector_alloc(pmod->t2 - pmod->t1 + 1);
	if (v == NULL) {
	    *err = E_ALLOC;
	} else {
	    for (t=pmod->t1; t<=pmod->t2; t++) {
		x = (idx == M_UHAT)? pmod->uhat[t] : pmod->yhat[t];
		gretl_vector_set(v, t - pmod->t1, x);
	    }
	}
    }

    return v;
}

static gretl_matrix *model_get_hvec (const MODEL *pmod, int *err)
{
    double *garch_h = gretl_model_get_data(pmod, "garch_h");
    gretl_matrix *v = NULL;
    int t;

    if (garch_h == NULL) {
	*err = E_BADSTAT;
    }

    v = gretl_column_vector_alloc(pmod->t2 - pmod->t1 + 1);
    if (v == NULL) {
	*err = E_ALLOC;
    } else {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    /* FIXME: is indexation right? */
	    gretl_vector_set(v, t - pmod->t1, garch_h[t]);
	}
    }

    return v;
}

static gretl_matrix *model_get_rhovec (const MODEL *pmod, int *err)
{
    gretl_matrix *r = NULL;

    if (AR1_MODEL(pmod->ci)) {
	double x = gretl_model_get_double(pmod, "rho_in");

	r = gretl_matrix_from_scalar(x);
    } else if (pmod->ci != AR) {
	r = gretl_matrix_from_scalar(pmod->rho);
    } else if (pmod->arinfo == NULL || 
	       pmod->arinfo->arlist == NULL || 
	       pmod->arinfo->rho == NULL) {
	*err = E_INVARG;
    } else {
	int i, l0 = pmod->arinfo->arlist[0];
	int lmax = pmod->arinfo->arlist[l0];

	r = gretl_vector_alloc(lmax);
	if (r != NULL) {
	    gretl_matrix_zero(r);
	    for (i=1; i<=l0; i++) {
		gretl_vector_set(r, pmod->arinfo->arlist[i] - 1,
				 pmod->arinfo->rho[i-1]);
	    }
	}
    }

    return r;
}

gretl_matrix *gretl_model_get_matrix (MODEL *pmod, int idx, int *err)
{
    gretl_matrix *M = NULL;

    if (pmod == NULL) {
	*err = E_BADSTAT;
	return M;
    }

    if (*err) return M;

    switch (idx) {  
    case M_UHAT:
    case M_YHAT:
	M = model_get_hatvec(pmod, idx, err);
	break;
    case M_COEFF:
    case M_SE:
	M = model_get_estvec(pmod, idx, err);
	break;
    case M_VCV:
	M = gretl_vcv_matrix_from_model(pmod, NULL);
	break;
    case M_H:
	if (pmod->ci != GARCH) {
	    *err = E_BADSTAT;
	} else {
	    M = model_get_hvec(pmod, err);
	}
	break;
    case M_RHO:
	M = model_get_rhovec(pmod, err);
	break;
    default:
	*err = E_BADSTAT;
	break;
    }

    if (M == NULL && !*err) {
	*err = E_ALLOC;
    }

    return M;
}
