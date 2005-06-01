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

#undef MODEL_DEBUG

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
	/* special treatment for ARCH, ARMA, GARCH, NLS */
	if (pmod->aux == AUX_ARCH) {
	    make_cname(pdinfo->varname[pmod->list[i + 2]], targ);
	} else if (pmod->ci == NLS || pmod->ci == ARMA || pmod->ci == GARCH) {
	    strcpy(targ, pmod->params[i + 1]);
	} else {
	    strcpy(targ, pdinfo->varname[pmod->list[i + 2]]);
	}
    }

    return targ;
}

/**
 * free_vcv:
 * @vcv: pointer to covariance matrix struct.
 * 
 * Frees the resources associated with @vcv, then frees the
 * pointer itself.
 */

void free_vcv (VCV *vcv)
{
    free(vcv->vec);
    free(vcv->list);
    free(vcv);
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
 * 
 * Supplies the caller with a copy of the variance-covariance 
 * matrix for the parameter estimates in @pmod.  See also
 * free_vcv().  To get the covariance matrix in gretl_matrix
 * format, see gretl_vcv_matrix_from_model().
 *
 * Returns: #VCV struct or %NULL on error.
 */

VCV *gretl_model_get_vcv (MODEL *pmod)
{
    int i, nv = pmod->ncoeff;
    VCV *vcv;

    vcv = malloc(sizeof *vcv);
    if (vcv == NULL) return NULL;

    vcv->list = malloc((nv + 1) * sizeof *vcv->list);
    if (vcv->list == NULL) {
	free(vcv);
	return NULL;
    }

    vcv->list[0] = nv;
    for (i=1; i<=nv; i++) {
	vcv->list[i] = pmod->list[i+1];
    }

    if (pmod->vcv == NULL && makevcv(pmod)) {
	free(vcv->list);
	free(vcv);
	return NULL;
    }

    /* calculate number of elements in vcv */
    nv = (nv * nv + nv) / 2;

    /* copy vcv */
    vcv->vec = copyvec(pmod->vcv, nv + 1);
    if (vcv->vec == NULL) {
	free(vcv->list);
	free(vcv);
	return NULL;
    }

    vcv->ci = pmod->ci;
    
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

    pmod->ID = 0;

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
    MODEL *pmod = malloc(sizeof *pmod);

    gretl_model_init(pmod);
    return pmod;
}

/* .......................................................... */

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

#ifdef MODEL_DEBUG

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

#endif /* MODEL_DEBUG */

static void gretl_test_free (ModelTest *test)
{
    if (test->param != NULL) {
	free(test->param);
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

#ifdef MODEL_DEBUG
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
	if (pmod->ntests) {
	    for (i=0; i<pmod->ntests; i++) {
		gretl_test_free(&pmod->tests[i]);
	    }
	    free(pmod->tests);
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
	destroy_all_data_items(pmod);
    }

    /* this may be redundant */
    gretl_model_init(pmod);
}

/* ........................................................... */

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
    test->teststat = 0;
    test->dfn = test->dfd = 0;
    test->value = test->pvalue = NADBL;
}

/**
 * new_test_on_model:
 * @pmod: pointer to model.
 * @ttype: type of test to add.
 *
 * Adds a #ModelTest to @pmod, if a test of the given type has
 * not already been performed and recorded.
 *
 * Returns: model test pointer, or %NULL if the test is
 * already present or on failure.
 */

ModelTest *new_test_on_model (MODEL *pmod, ModelTestType ttype)
{
    ModelTest *tests = NULL;
    ModelTest *ret = NULL;
    int i, nt = pmod->ntests;
    int done = 0;

    for (i=0; i<nt; i++) {
	if (ttype == pmod->tests[i].type) {
	    /* already done */
	    done = 1;
	}
    }

    if (!done) {
	tests = realloc(pmod->tests, (nt + 1) * sizeof *tests);
    }

    if (tests != NULL) {
	pmod->tests = tests;
	pmod->ntests += 1;
	ret = &pmod->tests[nt];
	gretl_test_init(ret, ttype);
    }

    return ret;
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
	N_("White's test for heteroskedasticity")
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
	N_("heteroskedasticity not present")
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

/**
 * copy_model:
 * @targ: pointer to #MODEL to copy to.
 * @src: pointer to #MODEL to copy from.
 * @pdinfo: pointer to dataset information.
 *
 * Does a deep copy of @src to @targ.  That is, @targ ends up with
 * its own allocated copies of all the pointer members of @src.
 *
 * Returns: 0 on success, 1 on failure.
 */

int copy_model (MODEL *targ, const MODEL *src, const DATAINFO *pdinfo)
{
    int i = src->list[0] - 1;
    int m = i * (i + 1) / 2;

    /* monolithic copy of structure */
    *targ = *src;

    /* now work on pointer members */
    gretl_model_init_pointers(targ);
    if ((targ->coeff = copyvec(src->coeff, src->ncoeff)) == NULL) 
	return 1;
    if ((targ->sderr = copyvec(src->sderr, src->ncoeff))  == NULL)   
	return 1;
    if ((targ->uhat = copyvec(src->uhat, pdinfo->n)) == NULL) 
	return 1;
    if ((targ->yhat = copyvec(src->yhat, pdinfo->n)) == NULL) 
	return 1;
    if (src->submask != NULL && 
	(targ->submask = copy_subsample_mask(src->submask, pdinfo->n)) == NULL) 
	return 1;
    if (src->missmask != NULL && 
	(targ->missmask = copy_missmask(src)) == NULL) 
	return 1;

    if (src->xpx != NULL &&
	(targ->xpx = copyvec(src->xpx, m)) == NULL) return 1;
    if (src->vcv != NULL && 
	(targ->vcv = copyvec(src->vcv, m)) == NULL) return 1;

    if (src->arinfo != NULL) {
	targ->arinfo = copy_ar_info(src->arinfo);
	if (targ->arinfo == NULL) 
	    return 1; 
    }

    if (src->ntests > 0 && src->tests != NULL) {
	copy_model_tests(targ, src);
	if (targ->tests == NULL) {
	    return 1;
	}
    }

    if (src->nparams > 0 && src->params != NULL) {
	copy_model_params(targ, src);
	if (targ->params == NULL) {
	    return 1;
	}
    }    

    if (src->n_data_items > 0) {
	copy_model_data_items(targ, src);
	if (targ->data_items == NULL) {
	    return 1;
	}
    }

    m = src->list[0];
    targ->list = malloc((m + 1) * sizeof *targ->list);
    if (targ->list == NULL) return 1;
    for (i=0; i<=m; i++) targ->list[i] = src->list[i];    

    return 0;
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

    *targ = *src;
    *src = tmp;
    return 0;
}

int is_model_cmd (const char *s)
{
    int ret = 0;

    if (s != NULL && *s != '\0') {
	if (!strncmp(s, "ols", 3)  ||
	    !strncmp(s, "corc", 4) ||
	    !strncmp(s, "hilu", 4) ||
	    !strncmp(s, "wls", 3)  ||
	    !strncmp(s, "pwe", 3)  ||
	    !strncmp(s, "pooled", 6)  ||
	    !strncmp(s, "hccm", 4) ||
	    !strncmp(s, "hsk", 3)  ||
	    !strncmp(s, "add", 3)  ||
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

    switch (test_ci) {
    case ADD:
    case ADDTO:
    case COEFFSUM:
    case VIF:
	if (model_ci == NLS || model_ci == TSLS ||
	    model_ci == ARMA || model_ci == GARCH) ok = 0;
	break;

    case OMIT:
    case OMITFROM:
	if (model_ci == NLS || model_ci == ARMA || 
	    model_ci == GARCH) ok = 0;
	break;

    case EQNPRINT:
	if (model_ci != OLS) ok = 0; /* FIXME: unduly restrictive? */
	break;

    case FCAST:
    case FIT:
	break;

    case FCASTERR:
	if (AR_MODEL(model_ci)) ok = 0;
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
    int i;
    char numstr[5];

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

int mle_aic_bic (MODEL *pmod, int addk)
{
    int err = 0;

    if (na(pmod->lnL)) {
	pmod->criterion[C_AIC] = NADBL;
	pmod->criterion[C_BIC] = NADBL;
	err = 1;
    } else {
	int k = pmod->ncoeff + addk;

	pmod->criterion[C_AIC] = -2.0 * pmod->lnL + 2.0 * k;
	pmod->criterion[C_BIC] = -2.0 * pmod->lnL + k * log(pmod->nobs);
    }

    return err;
}

/* As of 2005-02-19, not ready for this yet: need to adjust
   the text in the model printout.
*/

double coeff_pval (const MODEL *pmod, double x, int df)
{
    if (0 && ML_ESTIMATOR(pmod->ci)) {
	return normal_pvalue_2(x);
    } else {
	return t_pvalue_2(x, df);
    }
}



