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
#include "gretl_private.h"

#undef MODEL_DEBUG

struct model_data_item_ {
    char *key;
    void *ptr;
    size_t size;
};

static void free_item_data (const char *key, void *ptr)
{
    /* FIXME? (case of structs) */
    free(ptr);
}

static model_data_item *create_data_item (const char *key, void *ptr, size_t size)
{
    model_data_item *item;

    item = malloc(sizeof *item);
    if (item != NULL) {
	item->key = malloc(strlen(key) + 1);
	if (item->key == NULL) {
	    free(item);
	    item = NULL;
	} else {
	    strcpy(item->key, key);
	    item->ptr = ptr;
	    item->size = size;
	}
    }

    return item;
}

/**
 * gretl_model_set_data:
 * @pmod: pointer to #MODEL.
 * @key: key string for data, used in retrieval.
 * @ptr: data-pointer to be attached to model.
 * @size: size of data in bytes.
 *
 * Attaches data to a model: the data can be retrieved later using
 * gretl_model_get_data().  Note that the data are not "physically"
 * copied to the model; simply, the pointer is recorded.  The 
 * size is needed in case the model is copied.  The data pointer
 * will be freed when the model is cleared, with clear_model().
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_data (MODEL *pmod, const char *key, void *ptr, size_t size)
{
    model_data_item **items;
    model_data_item *item;
    int n_items = pmod->n_data_items + 1;

    items = realloc(pmod->data_items, n_items * sizeof *items);
    if (items == NULL) return 1;

    pmod->data_items = items;

    item = create_data_item(key, ptr, size);
    if (item == NULL) return 1;

    pmod->data_items[n_items - 1] = item;
    pmod->n_data_items += 1;

    return 0;
}

/**
 * gretl_model_set_int:
 * @pmod: pointer to #MODEL.
 * @key: key string, used in retrieval.
 * @val: integer value to set.
 *
 * Records an integer value on a model: the value can be retrieved 
 * later using gretl_model_get_int().  
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
 * @pmod: pointer to #MODEL.
 * @key: key string, used in retrieval.
 * @val: double-precision value to set.
 *
 * Records a floating-point value on a model: the value can be 
 * retrieved later using gretl_model_get_double().  
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
 * @pmod: pointer to #MODEL.
 * @key: key string.
 * @sz: pointer to receive the size of the data
 *
 * Returns the data pointer identified by @key, or %NULL on failure.
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
 * @pmod: pointer to #MODEL.
 * @key: key string.
 *
 * Returns the data pointer identified by @key, or %NULL on failure.
 */

void *gretl_model_get_data (const MODEL *pmod, const char *key)
{
    return gretl_model_get_data_and_size(pmod, key, NULL);
}

/**
 * gretl_model_get_int:
 * @pmod: pointer to #MODEL.
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
 * @pmod: pointer to #MODEL.
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
    model_data_item *item;
    int i;

    if (pmod->n_data_items == 0) return;

#ifdef HAVE_X12A
    maybe_delete_x12_file(pmod);
#endif

    for (i=0; i<pmod->n_data_items; i++) {
	item = pmod->data_items[i];
	free_item_data(item->key, item->ptr);
	free(item->key);
	free(item);
    }

    free(pmod->data_items);
    pmod->data_items = NULL;
}

int gretl_model_destroy_data_item (MODEL *pmod, const char *key)
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

	free(junk->key);
	free(junk);
    }

    return err;
}

/* .......................................................... */

void gretl_model_set_auxiliary (MODEL *pmod, int aux)
{
    pmod->aux = aux;
}

static void gretl_model_init_pointers (MODEL *pmod)
{
    pmod->list = NULL;
    pmod->subdum = NULL;
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
 * @pmod: pointer to #MODEL.
 *
 * Initializes a gretl #MODEL, including setting its pointer members
 * to %NULL.
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
 * @pmod: pointer to #MODEL.
 * @pdinfo: dataset information.
 *
 * Sets the start and end of the model's sample to the current dataset
 * values.
 */

void gretl_model_smpl_init (MODEL *pmod, const DATAINFO *pdinfo)
{
    pmod->smpl.t1 = pdinfo->t1;
    pmod->smpl.t2 = pdinfo->t2;
}

/**
 * gretl_model_new:
 * 
 * Allocates memory for a gretl MODEL struct and initializes the struct,
 * using gretl_model_init().
 *
 * Returns: pointer to #MODEL (or %NULL if allocation fails).
 */

MODEL *gretl_model_new (void)
{
    MODEL *pmod = malloc(sizeof *pmod);

    gretl_model_init(pmod);
    return pmod;
}

/**
 * exchange_smpl:
 * @pmod: pointer to #MODEL.
 * @pdinfo: pointer to data information struct.
 * 
 * Swaps the starting and ending values for the data sample
 * between @pmod and @pdinfo.
 */

void exchange_smpl (MODEL *pmod, DATAINFO *pdinfo)
{
    int t1 = pdinfo->t1, t2 = pdinfo->t2;

    pdinfo->t1 = pmod->smpl.t1;
    pdinfo->t2 = pmod->smpl.t2;

    pmod->smpl.t1 = t1;
    pmod->smpl.t2 = t2;
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
	    " pmod->subdum = %p\n"
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
	    (void *) pmod->subdum, (void *) pmod->missmask, 
	    (void *) pmod->coeff, (void *) pmod->sderr, 
	    (void *) pmod->yhat, (void *) pmod->uhat, 
	    (void *) pmod->xpx, (void *) pmod->vcv, 
	    (void *) pmod->name, (void *) pmod->params, 
	    (void *) pmod->arinfo, (void *) pmod->tests, 
	    (void *) pmod->data);
}

#endif /* MODEL_DEBUG */

/**
 * clear_model:
 * @pmod: pointer to #MODEL.
 *
 * Clears a gretl #MODEL, freeing all allocated storage and setting
 * pointer members to %NULL.  Also frees any data pointers attached
 * via gretl_model_set_data().
 */

void clear_model (MODEL *pmod)
{
    if (pmod != NULL) {
#ifdef MODEL_DEBUG
	debug_print_model_info(pmod, "Doing clear_model");
#endif
	if (pmod->list) free(pmod->list);
	if (pmod->subdum) free(pmod->subdum);
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
	if (pmod->ntests) free(pmod->tests);
	if (pmod->params) {
	    int i;

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

static void copy_test (GRETLTEST *targ, const GRETLTEST *src)
{
    strcpy(targ->type, src->type);
    strcpy(targ->h_0, src->h_0);
    strcpy(targ->param, src->param);
    targ->teststat = src->teststat;
    targ->dfn = src->dfn;
    targ->dfd = src->dfd;
    targ->value = src->value;
    targ->pvalue = src->pvalue;
}

static int copy_model_tests (MODEL *targ, const MODEL *src)
{
    int i, n = src->ntests;

    if (n <= 0 || src->tests == NULL) return 0;

    targ->tests = malloc(n * sizeof *targ->tests);
    if (targ->tests == NULL) return 1;

    for (i=0; i<n; i++) {
	copy_test(&targ->tests[i], &src->tests[i]);
    }

    return 0;
}

void gretl_test_init (GRETLTEST *test)
{
    test->type[0] = 0;
    test->h_0[0] = 0;
    test->param[0] = 0;
    test->teststat = 0;
    test->dfn = test->dfd = 0;
    test->value = test->pvalue = NADBL;
}

int add_test_to_model (MODEL *pmod, const GRETLTEST *test)
{
    GRETLTEST *tests;
    int i, nt = pmod->ntests;

    for (i=0; i<nt; i++) {
	if (!strcmp(test->type, pmod->tests[i].type)) {
	    /* already done */
	    return -1;
	}
    }

    tests = realloc(pmod->tests, (nt + 1) * sizeof *tests);
    if (tests == NULL) {
	return 1;
    }

    pmod->tests = tests;

    strcpy(pmod->tests[nt].type, test->type);
    strcpy(pmod->tests[nt].h_0, test->h_0);
    strcpy(pmod->tests[nt].param, test->param);

    pmod->tests[nt].teststat = test->teststat;
    pmod->tests[nt].value = test->value;
    pmod->tests[nt].dfn = test->dfn;
    pmod->tests[nt].dfd = test->dfd;
    pmod->tests[nt].pvalue = test->pvalue;

    pmod->ntests += 1;

    return 0;
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
    if (src->subdum != NULL && 
	(targ->subdum = copy_subdum(src->subdum, pdinfo->n)) == NULL) 
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

/**
 * command_ok_for_model:
 * @test_ci:  index of command to be tested.
 * @model_ci: command index of a gretl #MODEL (for example,
 * OLS, HCCM or CORC).
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
    case OMIT:
    case OMITFROM:
    case COEFFSUM:
    case VIF:
	if (model_ci == TSLS || model_ci == NLS || 
	    model_ci == ARMA || model_ci == GARCH) ok = 0;
	break;

    case EQNPRINT:
	if (model_ci != OLS) ok = 0; /* unduly restrictive? */
	break;

    case FCAST:
    case FIT:
	break;

    case FCASTERR:
	if (model_ci != OLS) ok = 0;
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
	/* need to exclude garch? */
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



