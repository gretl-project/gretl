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

/* .......................................................... */

struct _model_data_item {
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

void *gretl_model_get_data (const MODEL *pmod, const char *key)
{
    int i;

    for (i=0; i<pmod->n_data_items; i++) {
	if (!strcmp(key, (pmod->data_items[i])->key)) {
	    return (pmod->data_items[i])->ptr;
	}
    }

    return NULL;
}

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

/* .......................................................... */

void gretl_model_set_auxiliary (MODEL *pmod, int aux)
{
    pmod->aux = aux;
}

static void gretl_model_init_pointers (MODEL *pmod)
{
    pmod->list = NULL;
    pmod->subdum = NULL;
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

void gretl_model_init (MODEL *pmod, const DATAINFO *pdinfo)
{
    int i;

    if (pmod == NULL) return;

    pmod->ID = 0;
    if (pdinfo != NULL) {
	pmod->smpl.t1 = pdinfo->t1;
	pmod->smpl.t2 = pdinfo->t2;
    }

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
 * gretl_model_new:
 * @pdinfo: pointer to data information struct.
 * 
 * Allocates memory for a gretl MODEL struct and initializes the struct.
 *
 * Returns: pointer to #MODEL (or NULL if allocation fails).
 */

MODEL *gretl_model_new (const DATAINFO *pdinfo)
{
    MODEL *pmod = malloc(sizeof *pmod);

    gretl_model_init(pmod, pdinfo);
    return pmod;
}

/* ........................................................... */

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
    if (pmod->arinfo->arlist) free(pmod->arinfo->arlist);
    if (pmod->arinfo->rho) free(pmod->arinfo->rho);
    if (pmod->arinfo->sderr) free(pmod->arinfo->sderr);

    free(pmod->arinfo);
    pmod->arinfo = NULL;
}

#if 0

static void 
debug_print_model_info (const MODEL *pmod, const char *msg)
{
    fprintf(stderr, "%s:\n"
	    " pmod = %p\n"
	    " pmod->list = %p\n"
	    " pmod->subdum = %p\n"
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
	    (void *) pmod->subdum, (void *) pmod->coeff,
	    (void *) pmod->sderr, (void *) pmod->yhat, 
	    (void *) pmod->uhat, (void *) pmod->xpx,
	    (void *) pmod->vcv, (void *) pmod->name, 
	    (void *) pmod->params, (void *) pmod->arinfo, 
	    (void *) pmod->tests, (void *) pmod->data);
}

#endif

/* .......................................................... */

void clear_model (MODEL *pmod, const DATAINFO *pdinfo)
{
    if (pmod != NULL) {
#if 0
	debug_print_model_info(pmod, "Doing clear_model");
#endif
	if (pmod->list) free(pmod->list);
	if (pmod->subdum) free(pmod->subdum);
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
	if (pmod->data) {
	    MISSOBS *mobs = (MISSOBS *) pmod->data;
	    free(mobs->missvec);
	    free(pmod->data);
	}
	if (pmod->dataset) {
	    free_model_dataset(pmod);
	}
	destroy_all_data_items(pmod);
    }

    gretl_model_init(pmod, pdinfo);
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

static char *my_strdup (const char *src)
{
    char *targ = malloc(strlen(src) + 1);

    if (src != NULL) {
	strcpy(targ, src);
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
	targ->params[i] = my_strdup(src->params[i]);
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
	
	targitem->key = my_strdup(srcitem->key);
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

/* ........................................................... */

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

/* ........................................................... */

int swap_models (MODEL **targ, MODEL **src)
{
    MODEL *tmp = *targ;

    *targ = *src;
    *src = tmp;
    return 0;
}

/* ........................................................... */

int command_ok_for_model (int test_ci, int model_ci)
{
    int ok = 1;

    switch (test_ci) {
    case ADD:
    case ADDTO:
    case OMIT:
    case OMITFROM:
    case COEFFSUM:	
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

    case ARCH:
    case CHOW:
    case CUSUM:
    case LMTEST:
    case LEVERAGE:
    case RESET:
	if (model_ci != OLS) ok = 0;
	break;

    case HAUSMAN:
	if (model_ci != POOLED) ok = 0;
	break;
    case RESTRICT:
	if (model_ci == LAD) ok = 0;
	break;
    case TESTUHAT:
	if (model_ci == TOBIT || model_ci == GARCH) ok = 0; /* garch? */
	break;

    default:
	break;
    }

    return ok;
}

int model_ci_from_modelspec (MODELSPEC *spec)
{
    char mword[9];
    int ci;

    if (!sscanf(spec->cmd, "%8s", mword)) {
	ci = -1;
    } else {
	ci = gretl_command_number(mword);
    }

    return ci;
}
