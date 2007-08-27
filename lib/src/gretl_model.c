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
#include "modelspec.h"
#include "gretl_xml.h"
#include "matrix_extra.h"

#define MDEBUG 0

struct model_data_item_ {
    char *key;
    void *ptr;
    int type;
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
    double crit;
    double alpha;
};

typedef struct CoeffSep_ CoeffSep;

#define CSLEN 64

struct CoeffSep_ {
    char str[CSLEN];
    int pos;
};

#if 0 /* not ready yet: this apparatus may be used for printing
	 "matrix models" at some future point */

typedef struct gretl_matrix_model_ gretl_matrix_model;

struct gretl_matrix_model_ {
    gretl_matrix *b;
    gretl_matrix *se;
    char **bnames;
    int nobs;
    double R2;
};

gretl_matrix_model *gretl_matrix_model_new (void)
{
    gretl_matrix_model *m = malloc(sizeof *m);
    
    if (m != NULL) {
	m->b = NULL;
	m->se = NULL;
	m->bnames = NULL;
	m->nobs = 0;
	m->R2 = NADBL;
    }

    return m;
}

void gretl_matrix_model_destroy (gretl_matrix_model *m)
{
    if (m != NULL) {
	if (m->b != NULL) {
	    int k = gretl_vector_get_length(m->b);

	    free_strings_array(m->bnames, k);
	    gretl_matrix_free(m->b);
	}
	gretl_matrix_free(m->se);
	free(m);
    }
}

void gretl_matrix_model_set_data (gretl_matrix_model *m, gretl_matrix *b,
				  gretl_matrix *se, char **bnames,
				  int nobs, double R2)
{
    m->b = b;
    m->se = se;
    m->bnames = bnames;
    m->nobs = nobs;
    m->R2 = R2;
}

void gretl_matrix_model_get_data (gretl_matrix_model *m, gretl_matrix **b,
				  gretl_matrix **se, char ***bnames,
				  int *nobs, double *R2)
{
    *b = m->b;
    *se = m->se;
    *bnames = m->bnames;
    *nobs = m->nobs;
    *R2 = m->R2;
}

#endif

static void gretl_test_init (ModelTest *test, ModelTestType ttype)
{
    test->type = ttype;
    test->order = 0;
    test->param = NULL;
    test->teststat = GRETL_STAT_NONE;
    test->dfn = test->dfd = 0;
    test->value = test->pvalue = NADBL;
    test->crit = test->alpha = NADBL;
}

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

static model_data_item *create_data_item (const char *key, void *ptr, 
					  ModelDataType type, size_t size,
					  void (*destructor) (void *))
{
    model_data_item *item = malloc(sizeof *item);

    if (item != NULL) {
	item->key = gretl_strdup(key);
	if (item->key == NULL) {
	    free(item);
	    item = NULL;
	} else {
	    item->ptr = ptr;
	    item->type = type;
	    item->size = size;
	    item->destructor = destructor;
	}
    }

    return item;
}

static model_data_item *
replicate_data_item (const model_data_item *orig)
{
    model_data_item *item = malloc(sizeof *item);

    if (item != NULL) {
	item->key = gretl_strdup(orig->key);
	if (item->key == NULL) {
	    free(item);
	    item = NULL;
	}
    }

    if (item != NULL) {
	item->ptr = malloc(orig->size);
	if (item->ptr == NULL) {
	    free(item->key);
	    free(item);
	    item = NULL;
	}
    }
    
    if (item != NULL) {
	memcpy(item->ptr, orig->ptr, orig->size);
	item->type = orig->type;
	item->size = orig->size;
	item->destructor = orig->destructor;
    }

    return item;
}

/**
 * gretl_model_set_data_with_destructor:
 * @pmod: pointer to #MODEL.
 * @key: key string for data, used in retrieval.
 * @ptr: data-pointer to be attached to model.
 * @type: type of data to set.
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
					  ModelDataType type, size_t size, 
					  void (*destructor) (void *))
{
    model_data_item **items;
    model_data_item *item;
    int i, n;

    for (i=0; i<pmod->n_data_items; i++) {
	item = pmod->data_items[i];
	if (!strcmp(key, item->key)) {
	    /* there's a pre-existing item of this name */
	    if (item->destructor != NULL) {
		(*item->destructor)(item->ptr);
	    } else {
		free(item->ptr);
	    }
	    item->type = type;
	    item->ptr = ptr;
	    item->size = size;
	    item->destructor = destructor;
	    /* handled */
	    return 0;
	}
    } 

    n = pmod->n_data_items + 1;

    items = realloc(pmod->data_items, n * sizeof *items);
    if (items == NULL) {
	return 1;
    }

    pmod->data_items = items;

    item = create_data_item(key, ptr, type, size, destructor);
    if (item == NULL) {
	return 1;
    }

    pmod->data_items[n - 1] = item;
    pmod->n_data_items += 1;

    return 0;
}

/**
 * gretl_model_set_data:
 * @pmod: pointer to #MODEL.
 * @key: key string for data, used in retrieval.
 * @ptr: data-pointer to be attached to model.
 * @type: type of the data to set.
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

int gretl_model_set_data (MODEL *pmod, const char *key, void *ptr, 
			  ModelDataType type, size_t size)
{
    return gretl_model_set_data_with_destructor(pmod, key, ptr, type, 
						size, NULL);
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
						MODEL_DATA_LIST, size, 
						NULL);
}

/**
 * gretl_model_set_string_as_data:
 * @pmod: pointer to #MODEL.
 * @key: key string, used in retrieval.
 * @str: string to attach.
 *
 * Attaches @str to @pmod as data, recoverable via the key @key 
 * using gretl_model_get_data().
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_string_as_data (MODEL *pmod, const char *key, char *str)
{
    size_t size = strlen(str) + 1;

    return gretl_model_set_data_with_destructor(pmod, key, (void *) str, 
						MODEL_DATA_STRING, size, 
						NULL);
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

    err = gretl_model_set_data(pmod, key, valp, MODEL_DATA_INT, 
			       sizeof(int));
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

    err = gretl_model_set_data(pmod, key, valp, MODEL_DATA_DOUBLE,
			       sizeof(double));
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
	if (pmod->data_items[i]->type != MODEL_DATA_INT) {
	    continue;
	}
	if (!strcmp(key, pmod->data_items[i]->key)) {
	    valp = (int *) pmod->data_items[i]->ptr;
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
	if (pmod->data_items[i]->type != MODEL_DATA_DOUBLE) {
	    continue;
	}
	if (!strcmp(key, pmod->data_items[i]->key)) {
	    valp = (double *) pmod->data_items[i]->ptr;
	    return *valp;
	}
    }

    return NADBL;
}

/**
 * gretl_model_get_list:
 * @pmod: pointer to model.
 * @key: key string.
 *
 * Returns: the list of integers identified by @key, or 
 * %NULL on failure.
 */

int *gretl_model_get_list (const MODEL *pmod, const char *key)
{
    int *list = NULL;
    int i;

    for (i=0; i<pmod->n_data_items; i++) {
	if (pmod->data_items[i]->type != MODEL_DATA_LIST) {
	    continue;
	}
	if (!strcmp(key, pmod->data_items[i]->key)) {
	    list = (int *) pmod->data_items[i]->ptr;
	    break;
	}
    }

    return list;
}

/**
 * gretl_model_set_coeff_separator:
 * @pmod: pointer to model.
 * @s: informative string (or %NULL).
 * @pos: position in the array of coefficients.
 *
 * Arranges for the insertion of the given string (or a blank line if
 * @s is %NULL) at the given position in the array of coefficients,
 * when the model is printed.  The extra line is printed before 
 * coefficient @pos, where @pos is interpreted as a zero-based
 * index.
 *
 * Returns: 0 on success, %E_ALLOC on failure.
 */

int gretl_model_set_coeff_separator (MODEL *pmod, const char *s, int pos)
{
    CoeffSep *cs = malloc(sizeof *cs);
    int err = 0;

    if (cs == NULL) {
	return E_ALLOC;
    }

    cs->str[0] = '\0';
    if (s != NULL) {
	strncat(cs->str, s, CSLEN - 1);
    }
    cs->pos = pos;

    err = gretl_model_set_data(pmod, "coeffsep", cs, MODEL_DATA_STRUCT, 
			       sizeof *cs);
    if (err) {
	free(cs);
    }

    return err;
}

/**
 * gretl_model_get_coeff_separator:
 * @pmod: pointer to model.
 * @ps: location to receive string, if any.
 * @ppos: location to receive position, if any.
 *
 * Retrieves information that has been set on @pmod regarding the
 * insertion of an extra line when printing the coefficients, if any.
 *
 * Returns: 1 if such information is present, 0 otherwise.
 */

int gretl_model_get_coeff_separator (const MODEL *pmod, const char **ps, int *ppos)
{
    CoeffSep *cs = gretl_model_get_data(pmod, "coeffsep");

    if (cs == NULL) {
	return 0;
    }

    *ps = cs->str;
    *ppos = cs->pos;

    return 1;
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
	int j = i + 2;
	int k = -1;

	if (pmod->aux == AUX_ARCH) {
	    make_cname(pdinfo->varname[pmod->list[j]], targ);
	} else if (pmod->ci == NLS || pmod->ci == MLE || pmod->ci == GMM ||
		   pmod->ci == ARMA || pmod->ci == PANEL ||
		   pmod->ci == ARBOND || pmod->ci == GARCH) {
	    k = i;
	} else if (pmod->ci == MPOLS && pmod->params != NULL) {
	    k = i;
	} else if ((pmod->ci == PROBIT || pmod->ci == LOGIT ||
		    pmod->ci == HECKIT) && pmod->params != NULL) {
	    k = i;
	} else if (pmod->list == NULL || j > pmod->list[0]) {
	    k = i;
	} else {
	    strcpy(targ, pdinfo->varname[pmod->list[j]]);
	}

	if (k >= 0) {
	    if (pmod->params != NULL) {
		strcpy(targ, pmod->params[k]);
	    } else {
		strcpy(targ, "unknown");
	    }
	}

    }

    return targ;
}

/**
 * gretl_model_get_param_number:
 * @pmod: pointer to model.
 * @pdinfo: dataset information.
 * @s: name of model parameter.
 *
 * Returns the zero-based index of the coefficient in @pmod
 * corresponding to @s, or -1 if @s is not the name
 * of a parameter.
 */

int 
gretl_model_get_param_number (const MODEL *pmod, const DATAINFO *pdinfo,
			      const char *s)
{
    char pname[16];
    int idx = -1;

    if (pmod == NULL) {
	return -1;
    }

    if (!strcmp(s, "0")) {
	strcpy(pname, "const");
    } else {
	strcpy(pname, s);
    }

    if (pmod->params != NULL) {
	int i;

	for (i=0; i<pmod->nparams; i++) {
	    if (!strcmp(pname, pmod->params[i])) {
		idx = i;
		break;
	    }
	}
    } else if (pmod->list != NULL) {
	int v = varindex(pdinfo, pname);

	if (v < pdinfo->v) {
	    idx = gretl_list_position(v, pmod->list);
	    if (idx > 1) {
		idx -= 2;
	    } else {
		idx = -1;
	    }
	}
    }

    return idx;
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

    if (ASYMPTOTIC_MODEL(pmod->ci)) {
	cf->asy = 1;
	cf->t = normal_cdf_inverse(0.975);
    } else {
	cf->asy = 0;
	cf->t = tcrit95(pmod->dfd);
    }

    for (i=0; i<cf->ncoeff && !err; i++) { 
	cf->coeff[i] = pmod->coeff[i];
	cf->maxerr[i] = (pmod->sderr[i] > 0)? cf->t * pmod->sderr[i] : 0.0;
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

int gretl_is_arima_model (const MODEL *pmod)
{
    int d = gretl_model_get_int(pmod, "arima_d");
    int D = gretl_model_get_int(pmod, "arima_D");

    return (d > 0 || D > 0);
}

/**
 * arma_model_nonseasonal_AR_order:
 * @pmod: pointer to gretl model.
 *
 * Returns: the non-seasonal autoregressive order of @pmod, or 0 if
 * @pmod is not an ARMA model.
 */

int arma_model_nonseasonal_AR_order (const MODEL *pmod)
{
    int p = 0;

    if (pmod->ci == ARMA) {
	p = pmod->list[1];
    }

    return p;
}

/**
 * arma_model_nonseasonal_MA_order:
 * @pmod: pointer to gretl model.
 *
 * Returns: the non-seasonal moving-average order of @pmod, or 0 if
 * @pmod is not an ARMA model.
 */

int arma_model_nonseasonal_MA_order (const MODEL *pmod)
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
 * arma_model_max_AR_lag:
 * @pmod: pointer to gretl model.
 *
 * Returns: the maximum autoregressive lag in @pmod, or 0 if
 * @pmod is not an ARMA model. The maximum AR lag takes into
 * account any differencing (seasonal and/or non-seasonal) in
 * an ARIMA model.
 */

int arma_model_max_AR_lag (const MODEL *pmod)
{
    int pmax = 0;

    if (pmod->ci == ARMA) {
	int p = arma_model_nonseasonal_AR_order(pmod);
	int P = gretl_model_get_int(pmod, "arma_P");
	int s = gretl_model_get_int(pmod, "arma_pd");
	int d = gretl_model_get_int(pmod, "arima_d");
	int D = gretl_model_get_int(pmod, "arima_D");

	pmax = p + s * P;
	pmax += d + s * D;
    }

    return pmax;
}

/**
 * arma_model_max_MA_lag:
 * @pmod: pointer to gretl model.
 *
 * Returns: the maximum moving-average lag in @pmod, or 0 if
 * @pmod is not an ARMA model.
 */

int arma_model_max_MA_lag (const MODEL *pmod)
{
    int qmax = 0;

    if (pmod->ci == ARMA) {
	int q, Q, pd;

	q = arma_model_nonseasonal_MA_order(pmod);
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

/* From Box and Jenkins, 1976, pp 506-7, "Program 4": a clever
   algorithm for "unscrambling the coefficients", or in effect
   producing reduced-form AR coefficients that take into account any
   differencing, seasonal and/or non-seasonal.
*/

static int ar_coeff_integrate (double *c0, int d, int D, int s, int pmax)
{
    int pstar = pmax + d + s * D;
    double *c1;
    int i, j, pp;

    c1 = malloc((pstar + 1) * sizeof *c1);
    if (c1 == NULL) {
	return E_ALLOC;
    }

    for (j=0; j<=pstar; j++) {
	c1[j] = 0.0;
    }

    pp = pmax;

    for (i=0; D>0 && i<D; i++) {
	for (j=0; j<=pstar; j++) {
	    if (j < s) {
		c1[j] = c0[j];
	    } else if (j <= pp) {
		c1[j] = c0[j] - c0[j-s];
	    } else if (j <= pp + s) {
		c1[j] = -c0[j-s];
	    }
	}
	pp += s;
	for (j=0; j<=pstar; j++) {
	    c0[j] = c1[j];
	} 
    }

    for (i=0; d>0 && i<d; i++) {
	for (j=0; j<=pstar; j++) {
	    if (j < 1) {
		c1[j] = c0[j];
	    } else if (j <= pp) {
		c1[j] = c0[j] - c0[j-1];
	    } else if (j <= pp + 1) {
		c1[j] = -c0[j-1];
	    }		
	}
	pp += 1;
	for (j=0; j<=pstar; j++) {
	    c0[j] = c1[j];
	}
    }

    free(c1);

    return 0;
}

/**
 * arma_model_integrated_AR_MA_coeffs:
 * @pmod: pointer to gretl model.
 * @phi_star: pointer to receive AR coeff vector.
 * @theta_star: pointer to receive MA coeff vector.
 *
 * Creates consolidated versions of the AR and MA coefficient vectors
 * from @pmod.  If @pmod includes seasonal ARMA terms, the vectors are
 * suitably expanded to include the interactions between seasonal and
 * non-seasonal terms.  If the dependent variable has been
 * differenced, the AR coefficients are integrated to account for the
 * differencing.  These are the \Phi^* and \Theta^* as used by Box and
 * Jenkins for forecasting.
 *
 * The length of these vectors can be determined using 
 * gretl_arma_model_get_max_AR_lag() and 
 * gretl_arma_model_get_max_MA_lag() respectively.
 *
 * Returns: 0 on success, non-zero on error.
 */

int arma_model_integrated_AR_MA_coeffs (const MODEL *pmod,
					double **phi_star,
					double **theta_star)
{
    double *ac = NULL;
    double *mc = NULL;
    int err = 0;

    if (pmod->ci != ARMA) {
	err = 1;
    } else {
	const double *phi = NULL, *Phi = NULL;
	const double *theta = NULL, *Theta = NULL;
	double x, y;

	int p = arma_model_nonseasonal_AR_order(pmod);
	int q = arma_model_nonseasonal_MA_order(pmod);
	int P = gretl_model_get_int(pmod, "arma_P");
	int Q = gretl_model_get_int(pmod, "arma_Q");
	int d = gretl_model_get_int(pmod, "arima_d");
	int D = gretl_model_get_int(pmod, "arima_D");
	int s = gretl_model_get_int(pmod, "arma_pd");
	int pmax, pstar, qmax;
	int i, j, k;

	pmax = p + s * P;
	pstar = pmax + d + s * D;
	qmax = q + s * Q;
	
	if (pstar > 0) {
	    ac = malloc((pstar + 1) * sizeof *ac);
	    if (ac == NULL) {
		err = E_ALLOC;
	    }
	}

	if (!err && qmax > 0) {
	    mc = malloc((qmax + 1) * sizeof *ac);
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
	    for (i=0; i<=pstar; i++) {
		ac[i] = 0.0;
	    }
	    for (i=0; i<=P; i++) {
		x = (i == 0)? -1 : Phi[i-1];
		for (j=0; j<=p; j++) {
		    y = (j == 0)? -1 : phi[j-1];
		    k = j + s * i;
		    ac[k] -= x * y;
		}
	    }
	    if (D > 0 || d > 0) {
		ar_coeff_integrate(ac, d, D, s, pmax);
	    }
	}	

	if (mc != NULL) {
	    for (i=0; i<=qmax; i++) {
		mc[i] = 0.0;
	    }
	    for (i=0; i<=Q; i++) {
		x = (i == 0)? -1 : Theta[i-1];
		for (j=0; j<=q; j++) {
		    y = (j == 0)? -1 : theta[j-1];
		    k = j + s * i;
		    mc[k] -= x * y;
		}
	    }
	}
    }

    if (!err) {
	*phi_star = ac;
	*theta_star = mc;
    }

    return err;
}

/**
 * regarma_model_AR_coeffs:
 * @pmod: pointer to gretl model.
 * @phi0: pointer to receive AR coeff vector.
 * @pp: pointer to receive length of @phi0.
 *
 * Creates a consolidated version of the AR coefficients from @pmod.
 * If @pmod includes seasonal AR terms the vector is suitably expanded
 * to include the interactions between seasonal and non-seasonal
 * terms, but it is not integrated with respect to any differencing of
 * the dependent variable.
 *
 * Returns: 0 on success, non-zero on error.
 */

int regarma_model_AR_coeffs (const MODEL *pmod,
			     double **phi0,
			     int *pp)
{
    int p = arma_model_nonseasonal_AR_order(pmod);
    int P = gretl_model_get_int(pmod, "arma_P");
    int s = gretl_model_get_int(pmod, "arma_pd");
    const double *phi = NULL, *Phi = NULL;
    double *ac = NULL;
    double x, y;
    int i, j, k, pmax;
    int err = 0;

    pmax = p + s * P;

    if (pmax == 0) {
	*pp = 0;
	return 0;
    }
	
    ac = malloc((pmax + 1) * sizeof *ac);
    if (ac == NULL) {
	return E_ALLOC;
    }

    phi = pmod->coeff + pmod->ifc;
    Phi = phi + p;

    for (i=0; i<=pmax; i++) {
	ac[i] = 0.0;
    }

    for (i=0; i<=P; i++) {
	x = (i == 0)? -1 : Phi[i-1];
	for (j=0; j<=p; j++) {
	    y = (j == 0)? -1 : phi[j-1];
	    k = j + s * i;
	    ac[k] -= x * y;
	}
    }

    *phi0 = ac;
    *pp = pmax;

    return err;
}

/**
 * arma_model_get_x_coeffs:
 * @pmod: pointer to gretl model.
 *
 * Returns: pointer to the array of coefficients on the exogenous
 * regressors in @pmod, or %NULL if the model is not ARMA or if there
 * are no such regressors.
 */

const double *arma_model_get_x_coeffs (const MODEL *pmod)
{
    const double *xc = NULL;

    if (pmod->ci == ARMA && gretl_model_get_int(pmod, "armax")) {
	xc = pmod->coeff;
	xc += pmod->ifc;
	xc += arma_model_nonseasonal_AR_order(pmod);
	xc += arma_model_nonseasonal_MA_order(pmod);
	xc += gretl_model_get_int(pmod, "arma_P");
	xc += gretl_model_get_int(pmod, "arma_Q");
    }

    return xc;
}

static int arbond_get_depvar (const MODEL *pmod)
{
    int i;

    for (i=1; i<=pmod->list[0]; i++) {
	if (pmod->list[i] == LISTSEP && i < pmod->list[0]) {
	    return pmod->list[i+1];
	}
    }

    return 0;
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

    /* safety: should be impossible */
    if (dvpos == LISTSEP) {
	fprintf(stderr, "internal error in arma_depvar_pos\n");
	dvpos = 0;
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
    int dv = gretl_model_get_int(pmod, "yno");

    if (dv > 0) {
	return dv;
    }

    if (pmod != NULL && pmod->list != NULL) {
	if (pmod->ci == GARCH) {
	    dv = pmod->list[4];
	} else if (pmod->ci == ARMA) {
	    dv = pmod->list[arma_depvar_pos(pmod)];
	} else if (pmod->ci == ARBOND) {
	    dv = arbond_get_depvar(pmod);
	} else {
	    dv = pmod->list[1];
	}
    }

    return dv;
}

/**
 * gretl_model_get_depvar_name:
 * @pmod: pointer to gretl model.
 * @pdinfo: dataset information.
 *
 * Returns: the name of the dependent variable in @pmod.
 */

const char *gretl_model_get_depvar_name (const MODEL *pmod,
					 const DATAINFO *pdinfo)
{
    int dv;

    if (pmod->depvar != NULL) {
	return pmod->depvar;
    }

    dv = gretl_model_get_int(pmod, "yno");

    if (dv == 0) {
	if (pmod != NULL && pmod->list != NULL) {
	    if (pmod->ci == GARCH) {
		dv = pmod->list[4];
	    } else if (pmod->ci == ARMA) {
		dv = pmod->list[arma_depvar_pos(pmod)];
	    } else if (pmod->ci == ARBOND) {
		dv = arbond_get_depvar(pmod);
	    } else {
		dv = pmod->list[1];
	    }
	}
    }

    return pdinfo->varname[dv];
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
    } else if (pmod->ci == PANEL) {
	nx = pmod->list[0] - 1;
	list = gretl_list_new(nx);
	if (list != NULL) {
	    for (i=1; i<=list[0]; i++) {
		list[i] = pmod->list[i + 1];
	    }
	}
    } else if (pmod->ci == ARBOND) {
	int sep = 0;

	nx = 0;
	for (i=2; i<=pmod->list[0]; i++) {
	    if (pmod->list[i] == LISTSEP) {
		sep++;
		if (sep == 1) {
		    i += 2;
		} else if (sep == 2) {
		    break;
		}
	    }
	    if (sep == 1 && i <= pmod->list[0]) {
		list = gretl_list_append_term(&list, pmod->list[i]);
		if (list == NULL) {
		    return NULL;
		}
	    }
	}
    } else if (pmod->ci == HECKIT) {
	nx = gretl_model_get_int(pmod, "base-coeffs");
	if (nx > 0) {
	    list = gretl_list_new(nx);
	    if (list != NULL) {
		for (i=1; i<=list[0]; i++) {
		    list[i] = pmod->list[i + 1];
		}
	    }
	}	    
    } else if (pmod->ci != NLS && pmod->ci != MLE && pmod->ci != GMM) {
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
 * gretl_model_allocate_storage:
 * @pmod: pointer to model.
 * 
 * Allocates space for coefficients and standard errors,
 * residuals and fitted values in @pmod. The sizes of
 * the arrays are based on the %ncoeff and %full_n
 * members of @pmod, which must be set first. The
 * residuals and fitted values are initialized to
 * gretl's missing value.
 *
 * Returns: 0 on success, %E_ALLOC on error.
 */

int gretl_model_allocate_storage (MODEL *pmod)
{
    int k = pmod->ncoeff;
    int T = pmod->full_n;

    if (k > 0) {
	pmod->coeff = malloc(k * sizeof *pmod->coeff);
	if (pmod->coeff == NULL) {
	    return E_ALLOC;
	}
	pmod->sderr = malloc(k * sizeof *pmod->sderr);
	if (pmod->sderr == NULL) {
	    return E_ALLOC;
	}
    }

    if (T > 0) {
	int t;

	pmod->uhat = malloc(T * sizeof *pmod->uhat);
	if (pmod->uhat == NULL) {
	    return E_ALLOC;
	}    
	pmod->yhat = malloc(T * sizeof *pmod->yhat);
	if (pmod->yhat == NULL) {
	    return E_ALLOC;
	}
	for (t=0; t<T; t++) {
	    pmod->uhat[t] = pmod->yhat[t] = NADBL;
	}
    }

    return 0;
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

static double *copy_vcv_subset (const MODEL *pmod)
{
    double *V;
    int nc = pmod->ncoeff;
    int k = pmod->list[0] - 1;
    int n = k * (k + 1) / 2;
    int i, j;

    V = malloc(n * sizeof *V);
    if (V == NULL) {
	return NULL;
    }

    for (i=0; i<k; i++) {
	for (j=0; j<=i; j++) {
	    V[ijton(i, j, k)] = pmod->vcv[ijton(i, j, nc)];
	}
    }

    return V;
}

/**
 * gretl_model_get_vcv:
 * @pmod: pointer to model.
 * @pdinfo: dataset information.
 * 
 * Supplies the caller with a copy of the variance-covariance 
 * matrix for the parameter estimates in @pmod, in a format
 * suitable for printing.  See also free_vcv().  To get the 
 * covariance matrix as a gretl_matrix, see 
 * gretl_vcv_matrix_from_model().
 *
 * Returns: #VMatrix struct or %NULL on error.
 */

VMatrix *gretl_model_get_vcv (MODEL *pmod, const DATAINFO *pdinfo)
{
    char varname[VNAMELEN];
    int i, k = pmod->ncoeff;
    int special = 0;
    VMatrix *vcv;

    vcv = vmatrix_new();

    if (vcv == NULL) {
	return NULL;
    }

    /* special for fixed effects panel model: strip out the
       vcv elements for per-unit dummies, if present
    */
    if (pmod->ci == PANEL) {
	int k2 = pmod->list[0] - 1;

	if (k > k2) {
	    k = k2;
	    special = 1;
	}
    } 

    vcv->names = strings_array_new(k);
    if (vcv->names == NULL) {
	free(vcv);
	return NULL;
    }

    for (i=0; i<k; i++) {
	gretl_model_get_param_name(pmod, pdinfo, i, varname);
	vcv->names[i] = gretl_strdup(varname);
	if (vcv->names[i] == NULL) {
	    free_vmatrix(vcv);
	    return NULL;
	}
    }

    if (pmod->vcv == NULL && makevcv(pmod, pmod->sigma)) {
	free_vmatrix(vcv);
	return NULL;
    }

    if (special) {
	/* copy subset of vcv */
	vcv->vec = copy_vcv_subset(pmod);
    } else {
	/* copy full vcv */
	vcv->vec = copyvec(pmod->vcv, k * (k + 1)  / 2);
    }

    if (vcv->vec == NULL) {
	free_vmatrix(vcv);
	return NULL;
    }

    vcv->ci = pmod->ci;
    vcv->dim = k;
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
#if MDEBUG
    fprintf(stderr, "impose_model_smpl: set t1=%d, t2=%d\n", 
	    pdinfo->t1, pdinfo->t2);
#endif
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

/**
 * gretl_model_add_arinfo:
 * @pmod: pointer to model.
 * @nterms: number of autoregressive coefficients.
 *
 * Performs initial setup for structure to hold info 
 * on autoregressive coefficients.
 * 
 * Returns: 0 on success, 1 on error.
 */

int gretl_model_add_arinfo (MODEL *pmod, int nterms)
{
    int i;

    pmod->arinfo = malloc(sizeof *pmod->arinfo);
    if (pmod->arinfo == NULL) {
	return 1;
    }

    if (nterms == 0) {
	pmod->arinfo->arlist = NULL;
	pmod->arinfo->rho = NULL;
	pmod->arinfo->sderr = NULL;
	return 0;
    }

    pmod->arinfo->arlist = gretl_list_new(nterms);
    if (pmod->arinfo->arlist == NULL) {
	free(pmod->arinfo);
	pmod->arinfo = NULL;
	return 1; 
    }

    pmod->arinfo->rho = malloc(nterms * sizeof *pmod->arinfo->rho);
    if (pmod->arinfo->rho == NULL) {
	free(pmod->arinfo->arlist);
	free(pmod->arinfo);
	pmod->arinfo = NULL;
	return 1; 
    }

    pmod->arinfo->sderr = malloc(nterms * sizeof *pmod->arinfo->sderr);
    if (pmod->arinfo->sderr == NULL) {
	free(pmod->arinfo->arlist);
	free(pmod->arinfo->rho);
	free(pmod->arinfo);
	pmod->arinfo = NULL;
	return 1; 
    }

    for (i=0; i<nterms; i++) {
	pmod->arinfo->sderr[i] = pmod->arinfo->rho[i] = NADBL;
    }

    return 0;
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
    pmod->depvar = NULL;
    pmod->params = NULL;
    pmod->tests = NULL;
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
    pmod->refcount = 0;
    pmod->full_n = 0;
    pmod->t1 = 0;
    pmod->t2 = 0;
    pmod->nobs = 0;

    pmod->smpl.t1 = 0;
    pmod->smpl.t2 = 0;

    pmod->ncoeff = 0;
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
#if MDEBUG
    fprintf(stderr, "gretl_model_smpl_init: set t1=%d, t2=%d\n",
	    pdinfo->t1, pdinfo->t2);
#endif
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

#if MDEBUG
    fprintf(stderr, "gretl_model_new: model at %p\n", (void *) pmod);
    if (pmod != NULL) fprintf(stderr, " calling gretl_model_init\n");
#endif

    if (pmod != NULL) {
	gretl_model_init(pmod);
    }

    return pmod;
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
 * allocate_working_models:
 * @n: number of models in array.
 * 
 * Allocates memory for an array of @n gretl #MODEL structs and 
 * initializes each model using gretl_model_init().  The models
 * are "protected" against deletion.
 *
 * Returns: pointer to models array (or %NULL if allocation fails).
 */

MODEL **allocate_working_models (int n)
{
    MODEL **models;
    int i, err = 0;

    models = gretl_model_array_new(n);
    if (models == NULL) {
	return NULL;
    }

    for (i=0; i<n && !err; i++) {
	err = gretl_model_protect(models[i]);
    }

    if (err) {
	gretl_model_array_destroy(models, n);
	models = NULL;
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

void destroy_working_models (MODEL **models, int n)
{
    int i;

    if (models == NULL) {
	return;
    }

    for (i=0; i<n; i++) {
	gretl_model_free_on_exit(models[i]);
    }

    free(models);
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
	    " pmod->depvar = %p\n"
	    " pmod->params = %p\n"
	    " pmod->arinfo = %p\n"
	    " pmod->tests = %p\n" 
	    " pmod->data_items = %p\n", msg,
	    (void *) pmod, (void *) pmod->list, 
	    (void *) pmod->submask, (void *) pmod->missmask, 
	    (void *) pmod->coeff, (void *) pmod->sderr, 
	    (void *) pmod->yhat, (void *) pmod->uhat, 
	    (void *) pmod->xpx, (void *) pmod->vcv, 
	    (void *) pmod->name, (void *) pmod->depvar,
	    (void *) pmod->params, (void *) pmod->arinfo, 
	    (void *) pmod->tests, (void *) pmod->data_items);
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
	if (pmod->list != NULL) free(pmod->list);
	if (pmod->submask != NULL) free(pmod->submask);
	if (pmod->missmask != NULL) free(pmod->missmask);
	if (pmod->coeff != NULL) free(pmod->coeff);
	if (pmod->sderr != NULL) free(pmod->sderr);
	if (pmod->yhat != NULL) free(pmod->yhat);
	if (pmod->uhat != NULL) free(pmod->uhat);
	if (pmod->xpx != NULL) free(pmod->xpx);
	if (pmod->vcv != NULL) free(pmod->vcv);
	if (pmod->name != NULL) free(pmod->name);
	if (pmod->depvar != NULL) free(pmod->depvar);

	if (pmod->arinfo != NULL) {
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
 * Free allocated content of @pmod then the pointer itself.
 */

void gretl_model_free (MODEL *pmod)
{
    if (pmod != NULL) {
	clear_model(pmod);
#if MDEBUG
	fprintf(stderr, "gretl_model_free: freeing at %p\n", (void *) pmod);
#endif
	free(pmod);
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
#if MDEBUG
    fprintf(stderr, "gretl_model_free_on_exit: pmod at %p\n", (void *) pmod);
#endif
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
    targ->crit = src->crit;
    targ->alpha = src->alpha;
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

static int real_add_test_to_model (MODEL *pmod, ModelTest *test)
{
    ModelTest *tests = NULL;
    int nt = pmod->ntests;
    int err = 0;

    tests = realloc(pmod->tests, (nt + 1) * sizeof *tests);

    if (tests == NULL) {
	err = E_ALLOC;
    } else {
	pmod->tests = tests;
	pmod->ntests += 1;
	copy_test(&pmod->tests[nt], test);
    }

    return err;
}

static int attach_model_params_from_xml (xmlNodePtr node, xmlDocPtr doc,
					 MODEL *pmod)
{
    char **S;
    int np = 0;
    int err = 0;

    S = gretl_xml_get_strings_array(node, doc, &np, &err);
    if (!err) {
	pmod->params = S;
	pmod->nparams = np;
    }

    return err;
}

int attach_model_tests_from_xml (MODEL *pmod, xmlNodePtr node)
{
    ModelTest test;
    xmlNodePtr cur = node->xmlChildrenNode;
    int got, err = 0;

    gretl_test_init(&test, 0);

    while (cur != NULL && !err) {
	got = 0;
	got += gretl_xml_get_prop_as_int(cur, "type", &test.type);
	got += gretl_xml_get_prop_as_uchar(cur, "teststat", &test.teststat);
	got += gretl_xml_get_prop_as_int(cur, "dfn", &test.dfn);
	got += gretl_xml_get_prop_as_int(cur, "dfd", &test.dfd);
	got += gretl_xml_get_prop_as_int(cur, "order", &test.order);
	got += gretl_xml_get_prop_as_double(cur, "value", &test.value);
	got += gretl_xml_get_prop_as_double(cur, "pvalue", &test.pvalue);
	got += gretl_xml_get_prop_as_string(cur, "param", &test.param);
	got += gretl_xml_get_prop_as_double(cur, "crit", &test.crit);
	got += gretl_xml_get_prop_as_double(cur, "alpha", &test.alpha);
	if (got < 7) {
	    err = E_DATA;
	} else {
	    err = real_add_test_to_model(pmod, &test);
	}
	free(test.param);
	cur = cur->next;
    }

    return err;
}

static void serialize_test (const ModelTest *src, FILE *fp)
{
    fprintf(fp, "<test type=\"%d\" ", src->type);
    
    if (src->param != NULL) {
	fprintf(fp, "param=\"%s\" ", src->param);
    }

    fprintf(fp, "teststat=\"%d\" ", (int) src->teststat);
    fprintf(fp, "dfn=\"%d\" ", src->dfn);
    fprintf(fp, "dfd=\"%d\" ", src->dfd);
    fprintf(fp, "order=\"%d\" ", src->order);
    fprintf(fp, "value=\"%.15g\" ", src->value);
    fprintf(fp, "pvalue=\"%.15g\" ", src->pvalue);
    
    if (!na(src->crit)) {
	fprintf(fp, "crit=\"%g\" ", src->crit);
	fprintf(fp, "alpha=\"%g\" ", src->alpha);
    }

    fputs("/>\n", fp);
}

static int serialize_model_tests (const MODEL *pmod, FILE *fp)
{
    int i, n = pmod->ntests;

    if (n <= 0 || pmod->tests == NULL) {
	return 0;
    }

    fprintf(fp, "<tests count=\"%d\">\n", pmod->ntests);

    for (i=0; i<n; i++) {
	serialize_test(&pmod->tests[i], fp);
    }

    fputs("</tests>\n", fp);

    return 0;
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
 * model_test_free:
 * @test: object to free.
 */

void model_test_free (ModelTest *test)
{
    free(test->param);
    free(test);
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
 * not already been performed and recorded.  Note that this
 * function takes care of freeing @test.
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

void model_test_set_crit_and_alpha (ModelTest *test, 
				    double crit,
				    double alpha)
{
    test->crit = crit;
    test->alpha = alpha;
}

void model_test_set_param (ModelTest *test, const char *s)
{
    test->param = gretl_strdup(s);
}

void model_test_set_allocated_param (ModelTest *test, char *s)
{
    test->param = s;
}

struct test_strings {
    ModelTestType ID;
    const char *descrip;
    const char *H0;
};

static struct test_strings tstrings[] = {
    { GRETL_TEST_ADD,
      N_("Test for addition of variables"),
      N_("parameters are zero for the variables") },
    { GRETL_TEST_ARCH,
      N_("Test for ARCH of order %s"),
      N_("no ARCH effect is present") },
    { GRETL_TEST_AUTOCORR,
      N_("LM test for autocorrelation up to order %s"),
      N_("no autocorrelation") },
    { GRETL_TEST_CHOW,
      N_("Chow test for structural break at observation %s"),
      N_("no structural break") },
    { GRETL_TEST_CUSUM,
      N_("CUSUM test for parameter stability"),
      N_("no change in parameters") },
    { GRETL_TEST_QLR,
      N_("QLR test for structural break"),
      N_("no structural break") },
    { GRETL_TEST_GROUPWISE,
      N_("Likelihood ratio test for groupwise heteroskedasticity"),
      N_("the units have a common error variance") },
    { GRETL_TEST_LOGS,
      N_("Non-linearity test (logs)"),
      N_("relationship is linear") },
    { GRETL_TEST_NORMAL,
      N_("Test for normality of residual"),
      N_("error is normally distributed") },
    { GRETL_TEST_OMIT,
      N_("Test for omission of variables"),
      N_("parameters are zero for the variables") },
    { GRETL_TEST_RESET,
      N_("RESET test for specification"),
      N_("specification is adequate") },
    { GRETL_TEST_SQUARES,
      N_("Non-linearity test (squares)"),
      N_("relationship is linear") },
    { GRETL_TEST_WHITES,
      N_("White's test for heteroskedasticity"),
      N_("heteroskedasticity not present") },
    { GRETL_TEST_SARGAN,
      N_("Sargan over-identification test"),
      N_("all instruments are valid") },
    { GRETL_TEST_TSLS_HAUSMAN,
      N_("Hausman test"),
      N_("OLS estimates are consistent") },
    { GRETL_TEST_PANEL_HAUSMAN,
      N_("Hausman test"),
      N_("GLS estimates are consistent") },
    { GRETL_TEST_PANEL_F,
      N_("Test for differing group intercepts"),
      N_("The groups have a common intercept") },
    { GRETL_TEST_PANEL_BP,
      N_("Breusch-Pagan test"),
      N_("Variance of the unit-specific error = 0") },
    { GRETL_TEST_PANEL_TIMEDUM,
      N_("Wald test for joint significance of time dummies"),
      NULL },
    { GRETL_TEST_MAX, NULL, NULL }
};   

static int gretl_test_print_heading (const ModelTest *test, PRN *prn)
{
    const char *descrip = NULL;
    char *param = NULL;
    char ordstr[16];
    int i;

    for (i=0; tstrings[i].ID < GRETL_TEST_MAX; i++) {
	if (test->type == tstrings[i].ID) {
	    descrip = tstrings[i].descrip;
	}
    }

    if (descrip == NULL) {
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
	    pprintf(prn, _(descrip), param);
	} else {
	    pprintf(prn, I_(descrip), param);
	}
    } else {
	if (plain_format(prn)) {
	    pputs(prn, _(descrip));
	} else {
	    pputs(prn, I_(descrip));
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

static void gretl_test_print_h_0 (const ModelTest *test, int heading,
				  PRN *prn)
{
    const char *H0 = NULL;
    int i;

    for (i=0; tstrings[i].ID < GRETL_TEST_MAX; i++) {
	if (test->type == tstrings[i].ID) {
	    H0 = tstrings[i].H0;
	}
    }

    if (H0 == NULL) {
	return;
    }    

    if (heading) {
	pputs(prn, " -");
	if (tex_format(prn)) {
	    pputc(prn, '-');
	}
    }

    if (plain_format(prn)) {
	pprintf(prn, "\n  %s: ", _("Null hypothesis"));
	pputs(prn, _(H0));
    } else if (tex_format(prn)) {
	pprintf(prn, "\\\\\n\\quad %s: ", I_("Null hypothesis"));
	pputs(prn, I_(H0));
    } else if (rtf_format(prn)) {
	pprintf(prn, "\\par\n %s: ", I_("Null hypothesis"));
	pputs(prn, I_(H0));
    }

    if (test->type == GRETL_TEST_ADD || test->type == GRETL_TEST_OMIT) {
	print_add_omit_varnames(test->param, prn);
    } 
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
    case GRETL_STAT_SUP_WALD:	
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

void gretl_model_test_print_direct (const ModelTest *test, int heading, PRN *prn)
{
    char buf[128];
    const char *tstat;

    if (test->teststat == GRETL_STAT_WALD_CHISQ) {
	tstat = N_("Asymptotic test statistic");
    } else {
	tstat = N_("Test statistic");
    }

    if (rtf_format(prn)) {
	pputs(prn, "\\par \\ql ");
    }

    get_test_stat_string(test, buf, prn);

    if (heading) {
	gretl_test_print_heading(test, prn);
    }

    gretl_test_print_h_0(test, heading, prn);

    if (plain_format(prn)) {
	pprintf(prn, "\n  %s: %s\n", _(tstat), buf);
    } else if (tex_format(prn)) {
	pprintf(prn, "\\\\\n\\quad %s: %s\\\\\n", I_(tstat), buf);
    } else if (rtf_format(prn)) {
	pprintf(prn, "\\par\n %s: %s\\par\n", I_(tstat), buf);
    }

    get_test_pval_string(test, buf, prn);

    if (*buf) {
	if (plain_format(prn)) {
	    pprintf(prn, "  %s = %s\n\n", _("with p-value"), buf);
	} else if (tex_format(prn)) {
	    pprintf(prn, "\\quad %s = %s\\\\\n", I_("with p-value"), buf);
	} else if (rtf_format(prn)) {
	    pprintf(prn, " %s = %s\\par\n\n", I_("with p-value"), buf);
	}
    } else if (!na(test->crit) && !na(test->alpha)) {
	double a = test->alpha * 100.0;

	if (plain_format(prn)) {
	    sprintf(buf, _("%g percent critical value"), a);
	    pprintf(prn, "  (%s = %.2f)\n\n", buf, test->crit);
	} else if (tex_format(prn)) {
	    sprintf(buf, I_("%g percent critical value"), a);
	    pprintf(prn, "\\quad (%s = %.2f)\\\\\n", buf, test->crit);
	} else if (rtf_format(prn)) {
	    sprintf(buf, I_("%g percent critical value"), a);
	    pprintf(prn, " (%s = %.2f)\\par\n\n", buf, test->crit);
	}
    }
}

void gretl_model_test_print (const MODEL *pmod, int i, PRN *prn)
{
    const ModelTest *test;

    if (i < pmod->ntests) {
	test = &pmod->tests[i];
	gretl_model_test_print_direct(test, 1, prn);
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
	int m = src->arlist[0];

      	if ((targ->rho = copyvec(src->rho, m)) == NULL) {
	    free(targ);
	    return NULL; 
	}

      	if ((targ->sderr = copyvec(src->sderr, m)) == NULL) { 
	    free(targ->rho);
	    free(targ);
	    return NULL; 
	}

	targ->arlist = gretl_list_copy(src->arlist);
	if (targ->arlist == NULL) {
	    free(targ->rho);
	    free(targ->sderr);
	    free(targ);
	    return NULL; 
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

static void print_model_coeff_sep (model_data_item *item, FILE *fp)
{
    CoeffSep *cs = (CoeffSep *) item->ptr;

    fprintf(fp, " pos=\"%d\"", cs->pos);
    if (cs->str != NULL) {
	fputs(" string=\"", fp);
	gretl_xml_put_raw_string(cs->str, fp);
	fputc('"', fp);
    }
    fputs("/>\n", fp);
}

/* FIXME updating and placement of these function */

struct type_mapper {
    int type;
    const char *name;
    const char *compat;
};

static struct type_mapper mapper[] = {
    { MODEL_DATA_INT,          "int",         "1" },
    { MODEL_DATA_LIST,         "list",        "2" },
    { MODEL_DATA_DOUBLE,       "double",      "3" },
    { MODEL_DATA_INT_ARRAY,    "intarray",    "4" },
    { MODEL_DATA_DOUBLE_ARRAY, "doublearray", "5" },
    { MODEL_DATA_STRING,       "string",      "6" },
    { MODEL_DATA_CHAR_ARRAY,   "chararray",   "7" },
    { MODEL_DATA_CMPLX_ARRAY,  "cmplxarray",  "8" },
    { MODEL_DATA_STRUCT,       "struct" ,     "9" },
    { MODEL_DATA_NONE,         "none",        "0" },
};

static int type_from_type_string (const char *s)
{
    int i;

    if (isdigit(*s)) {
	/* backward compatibility */
	for (i=0; mapper[i].type != MODEL_DATA_NONE; i++) {
	    if (*s == mapper[i].compat[0]) {
		return mapper[i].type;
	    }
	}
    } else {
	for (i=0; mapper[i].type != MODEL_DATA_NONE; i++) {
	    if (!strcmp(s, mapper[i].name)) {
		return mapper[i].type;
	    }
	}
    }

    return MODEL_DATA_NONE;
}

static const char *gretl_type_name (int t)
{
    int i;

    for (i=0; mapper[i].type != MODEL_DATA_NONE; i++) {
	if (t == mapper[i].type) {
	    return mapper[i].name;
	}
    }

    return "unknown";
}
	
static void serialize_model_data_items (const MODEL *pmod, FILE *fp)
{
    model_data_item *item;
    int i, j, nelem;

    fprintf(fp, "<data-items count=\"%d\">\n", pmod->n_data_items);

    for (i=0; i<pmod->n_data_items; i++) {
	item = pmod->data_items[i];
	nelem = 0;

	fprintf(fp, "<data-item key=\"%s\" type=\"%s\"", 
		item->key, gretl_type_name(item->type));

	if (!strcmp(item->key, "coeffsep")) {
	    print_model_coeff_sep(item, fp);
	    continue;
	}

	if (item->type == MODEL_DATA_INT_ARRAY) {
	    nelem = item->size / sizeof(int);
	} else if (item->type == MODEL_DATA_DOUBLE_ARRAY) {
	    nelem = item->size / sizeof(double);
	} else if (item->type == MODEL_DATA_CMPLX_ARRAY) {
	    nelem = item->size / sizeof(cmplx);
	}

	if (nelem > 0) {
	    fprintf(fp, " count=\"%d\">\n", nelem);
	} else {
	    fputs(">\n", fp);
	}

	if (item->type == MODEL_DATA_INT) {
	    fprintf(fp, "%d", *(int *) item->ptr);
	} else if (item->type == MODEL_DATA_DOUBLE) {
	    fprintf(fp, "%.15g", *(double *) item->ptr);
	} else if (item->type == MODEL_DATA_INT_ARRAY) {
	    int *vals = (int *) item->ptr;

	    for (j=0; j<nelem; j++) {
		fprintf(fp, "%d ", vals[i]);
	    }
	} else if (item->type == MODEL_DATA_DOUBLE_ARRAY) {
	    double *vals = (double *) item->ptr;

	    for (j=0; j<nelem; j++) {
		fprintf(fp, "%.15g ", vals[i]);
	    }	    
	} else if (item->type == MODEL_DATA_CMPLX_ARRAY) {
	    cmplx *vals = (cmplx *) item->ptr;
	    
	    for (j=0; j<nelem; j++) {
		fprintf(fp, "%.15g %.15g ", vals[j].r, vals[j].i);
	    }	    
	} else if (item->type == MODEL_DATA_LIST) {
	    int *list = (int *) item->ptr;

	    for (j=0; j<=list[0]; j++) {
		fprintf(fp, "%d ", list[j]);
	    }
	} else if (item->type == MODEL_DATA_STRING) {
	    fprintf(fp, "%s", (char *) item->ptr);
	} else {
	    ; /* no-op: not handled */
	}

	fputs("</data-item>\n", fp);
    }    

    fputs("</data-items>\n", fp);
}

static int copy_model_data_items (MODEL *targ, const MODEL *src)
{
    int i, n = src->n_data_items;
    int err = 0;

    targ->data_items = malloc(n * sizeof *targ->data_items);
    if (targ->data_items == NULL) {
	return 1;
    }

    for (i=0; i<n; i++) {
	targ->data_items[i] = NULL;
    }

    for (i=0; i<n; i++) {
	targ->data_items[i] = replicate_data_item(src->data_items[i]);
	if (targ->data_items[i] == NULL) {
	    err = 1;
	    break;
	}
    }

    if (err) {
	for (i=0; i<n; i++) {
	    free(targ->data_items[i]);
	}
	free(targ->data_items);
	targ->data_items = NULL;
	targ->n_data_items = 0;
    } else {
	targ->n_data_items = n;
    }

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

    if (src->name != NULL) {
	targ->name = gretl_strdup(src->name);
	if (targ->name == NULL) {
	    return 1;
	}
    }    

    if (src->depvar != NULL) {
	targ->depvar = gretl_strdup(src->depvar);
	if (targ->depvar == NULL) {
	    return 1;
	}
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

    /* src->dataset? */

    return 0;
}

int gretl_model_serialize (const MODEL *pmod, SavedObjectFlags flags,
			   FILE *fp)
{
    int k = pmod->ncoeff;
    int m = k * (k + 1) / 2;
    int err = 0;

    fprintf(fp, "<gretl-model ID=\"%d\" name=\"%s\" saveflags=\"%d\" ", 
	    pmod->ID, (pmod->name == NULL)? "none" : pmod->name,
	    (int) flags);

    if (pmod->depvar != NULL) {
	fprintf(fp, "depvar=\"%s\" ", pmod->depvar);
    }

    fprintf(fp, "t1=\"%d\" t2=\"%d\" nobs=\"%d\" ",
	    pmod->t1, pmod->t2, pmod->nobs);
    fprintf(fp, "full_n=\"%d\" ncoeff=\"%d\" dfn=\"%d\" dfd=\"%d\" ", 
	    pmod->full_n, pmod->ncoeff, pmod->dfn, pmod->dfd);
    fprintf(fp, "ifc=\"%d\" ci=\"%s\" nwt=\"%d\" aux=\"%d\" ", 
	    pmod->ifc, gretl_command_word(pmod->ci), pmod->nwt, pmod->aux);

    gretl_push_c_numeric_locale();

    gretl_xml_put_double("ess", pmod->ess, fp);
    gretl_xml_put_double("tss", pmod->tss, fp);
    gretl_xml_put_double("sigma", pmod->sigma, fp);
    gretl_xml_put_double("rsq", pmod->rsq, fp);
    gretl_xml_put_double("adjrsq", pmod->adjrsq, fp);
    gretl_xml_put_double("fstt", pmod->fstt, fp);
    gretl_xml_put_double("lnL", pmod->lnL, fp);
    gretl_xml_put_double("ybar", pmod->ybar, fp);
    gretl_xml_put_double("sdy", pmod->sdy, fp);

    gretl_xml_put_double("crit0", pmod->criterion[0], fp);
    gretl_xml_put_double("crit1", pmod->criterion[1], fp);
    gretl_xml_put_double("crit2", pmod->criterion[2], fp);

    gretl_xml_put_double("dw", pmod->dw, fp);
    gretl_xml_put_double("rho", pmod->rho, fp);

    fputs(">\n", fp);

    fprintf(fp, "<sample t1=\"%d\" t2=\"%d\"/>\n",
	    pmod->smpl.t1, pmod->smpl.t2);

    gretl_push_c_numeric_locale();

    gretl_xml_put_double_array("coeff", pmod->coeff, k, fp);
    gretl_xml_put_double_array("sderr", pmod->sderr, k, fp);

    if (pmod->uhat != NULL) {
	gretl_xml_put_double_array("uhat", pmod->uhat, pmod->full_n, fp);
    }

    if (pmod->yhat != NULL) {
	gretl_xml_put_double_array("yhat", pmod->yhat, pmod->full_n, fp);
    }

    if (pmod->submask != NULL) {
	write_model_submask(pmod, fp);
    }

    if (pmod->missmask != NULL) {
	fputs("<missmask>", fp);
	fputs(pmod->missmask, fp);
	fputs("</missmask>\n", fp);
    }
	
    if (pmod->xpx != NULL) {
	gretl_xml_put_double_array("xpx", pmod->xpx, m, fp);
    }

    if (pmod->vcv != NULL) {
	gretl_xml_put_double_array("vcv", pmod->vcv, m, fp);
    }

    if (pmod->arinfo != NULL && pmod->arinfo->arlist != NULL) {
	int r = pmod->arinfo->arlist[0];

	fputs("<arinfo>\n", fp);
	gretl_xml_put_tagged_list("arlist", pmod->arinfo->arlist, fp);
	gretl_xml_put_double_array("rho", pmod->arinfo->rho, r, fp);
	gretl_xml_put_double_array("sderr", pmod->arinfo->sderr, r, fp);
	fputs("</arinfo>\n", fp);
    }

    if (pmod->ntests > 0 && pmod->tests != NULL) {
	serialize_model_tests(pmod, fp);
    }

    if (pmod->nparams > 0 && pmod->params != NULL) {
	gretl_xml_put_strings_array("params", 
				    (const char **) pmod->params, 
				    pmod->nparams, fp);
    } 

    if (pmod->list != NULL) {
	gretl_xml_put_tagged_list("list", pmod->list, fp);
    }

    if (pmod->n_data_items > 0) {
	serialize_model_data_items(pmod, fp);
    }

    fputs("</gretl-model>\n", fp);

    /* note: the DATASET element of pmod is not
       handled here */

    gretl_pop_c_numeric_locale();

    return err;
}

/* next block: functions for reconstituting a model from
   its XML representation */

static int model_submask_from_xml (xmlNodePtr node, xmlDocPtr doc,
				   MODEL *pmod)
{
    char *mask;
    int err;

    err = gretl_xml_get_submask(node, doc, &mask, NULL);
    if (!err) {
	pmod->submask = mask;
    }

    return err;
}

static int 
retrieve_model_coeff_separator (xmlNodePtr cur, MODEL *pmod)
{
    CoeffSep *cs = malloc(sizeof *cs);
    char *tmp = NULL;
    int err = 0;

    if (cs != NULL) {
	cs->str[0] = '\0';
	gretl_xml_get_prop_as_int(cur, "pos", &cs->pos);
	gretl_xml_get_prop_as_string(cur, "string", &tmp);
	if (tmp != NULL) {
	    strcpy(cs->str, tmp);
	    free(tmp);
	}
	err = gretl_model_set_data(pmod, "coeffsep", cs, 
				   MODEL_DATA_STRUCT, 
				   sizeof *cs);
    }

    return err;
}

static int model_data_items_from_xml (xmlNodePtr node, xmlDocPtr doc,
				      MODEL *pmod)
{
    xmlNodePtr cur;
    int n_items;
    int err = 0;

    if (!gretl_xml_get_prop_as_int(node, "count", &n_items)) {
	return 1;
    }

    cur = node->xmlChildrenNode;

    while (cur != NULL && !err) {
	char *key = NULL;
	char *typestr = NULL;
	int t, nelem = 0;

	if (!gretl_xml_get_prop_as_string(cur, "type", &typestr) ||
	    !gretl_xml_get_prop_as_string(cur, "key", &key)) {
	    err = 1;
	    break;
	}

	t = type_from_type_string(typestr);

	if (!strcmp(key, "coeffsep") || !strcmp(key, "10")) {
	    /* special, with backward compatibility */
	    err = retrieve_model_coeff_separator(cur, pmod);
	} else if (t == MODEL_DATA_INT) {
	    int ival;

	    if (!gretl_xml_node_get_int(cur, doc, &ival)) {
		err = 1;
	    } else {
		err = gretl_model_set_int(pmod, key, ival);
	    }
	} else if (t == MODEL_DATA_DOUBLE) {
	    double xval;

	    if (!gretl_xml_node_get_double(cur, doc, &xval)) {
		err = 1;
	    } else {
		err = gretl_model_set_double(pmod, key, xval);
	    }
	} else if (t == MODEL_DATA_LIST) {
	    if (!strcmp(key, "xlist")) {
		/* ad hoc (for forecasting): will be recreated if need be */
		;
	    } else {
		int *list = gretl_xml_node_get_list(cur, doc, &err);

		if (!err && list != NULL) {
		    err = gretl_model_set_list_as_data(pmod, key, list);
		} 
	    }
	} else if (t == MODEL_DATA_STRING) {
	    char *s;

	    if (!gretl_xml_node_get_string(cur, doc, &s)) {
		err = 1;
	    } else {
		err = gretl_model_set_string_as_data(pmod, key, s);
	    }
	} else if (t == MODEL_DATA_INT_ARRAY) {
	    int *ivals = gretl_xml_get_int_array(cur, doc, &nelem, &err);

	    if (nelem > 0) {
		err = gretl_model_set_data(pmod, key, ivals, t,
					   nelem * sizeof *ivals);
	    }
	} else if (t == MODEL_DATA_DOUBLE_ARRAY) {
	    double *xvals = gretl_xml_get_double_array(cur, doc, &nelem, &err);

	    if (nelem > 0) {
		err = gretl_model_set_data(pmod, key, xvals, t,
					   nelem * sizeof *xvals);
	    }
	} else if (t == MODEL_DATA_CMPLX_ARRAY) {
	    cmplx *cvals = gretl_xml_get_cmplx_array(cur, doc, &nelem, &err);

	    if (nelem > 0) {
		err = gretl_model_set_data(pmod, key, cvals, t,
					   nelem * sizeof *cvals);
	    }	    
	} 

	xmlFree(key);
	xmlFree(typestr);

	cur = cur->next;
    }

    return err;
}

static int arinfo_from_xml (xmlNodePtr node, xmlDocPtr doc,
			    MODEL *pmod)
{
    xmlNodePtr cur;
    int n, err = 0;

    if (gretl_model_add_arinfo(pmod, 0)) {
	return 1;
    }

    cur = node->xmlChildrenNode;

    while (cur != NULL && !err) {
	if (!xmlStrcmp(cur->name, (XUC) "arlist")) {
	    pmod->arinfo->arlist = gretl_xml_node_get_list(cur, doc, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "rho")) {
	    pmod->arinfo->rho = gretl_xml_get_double_array(cur, doc, &n, &err);
	} else if (!xmlStrcmp(cur->name, (XUC) "sderr")) {
	    pmod->arinfo->sderr = gretl_xml_get_double_array(cur, doc, &n, &err);
	}
	cur = cur->next;
    }

    return err;
}

/**
 * gretl_model_from_XML:
 * @node: XML node from which to read.
 * @doc: pointer to XML document.
 * @err: location to receive error code.
 *
 * Reads info on a gretl model from the given XML node
 * and doc, and reconstitutes the model in memory.
 *
 * Returns: allocated model, or %NULL on failure.
 */

MODEL *gretl_model_from_XML (xmlNodePtr node, xmlDocPtr doc, int *err)
{
    MODEL *pmod;
    char *buf = NULL;
    xmlNodePtr cur;
    int n, got = 0;

    pmod = gretl_model_new();
    if (pmod == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    got += gretl_xml_get_prop_as_int(node, "ID", &pmod->ID);
    got += gretl_xml_get_prop_as_int(node, "t1", &pmod->t1);
    got += gretl_xml_get_prop_as_int(node, "t2", &pmod->t2);
    got += gretl_xml_get_prop_as_int(node, "nobs", &pmod->nobs);
    got += gretl_xml_get_prop_as_int(node, "full_n", &pmod->full_n);
    got += gretl_xml_get_prop_as_int(node, "ncoeff", &pmod->ncoeff);
    got += gretl_xml_get_prop_as_int(node, "dfn", &pmod->dfn);
    got += gretl_xml_get_prop_as_int(node, "dfd", &pmod->dfd);
    got += gretl_xml_get_prop_as_int(node, "ifc", &pmod->ifc);
    got += gretl_xml_get_prop_as_int(node, "nwt", &pmod->nwt);
    got += gretl_xml_get_prop_as_int(node, "aux", &pmod->aux);

    got += gretl_xml_get_prop_as_string(node, "ci", &buf);

    if (got < 12) {
	fprintf(stderr, "gretl_model_from_XML: got(1) = %d (expected 12)\n", got);
	if (buf != NULL) {
	    free(buf);
	}
	*err = E_DATA;
	goto bailout;
    }

    if (isdigit(*buf)) {
	/* backward compatibility */
	pmod->ci = atoi(buf);
    } else {
	pmod->ci = gretl_command_number(buf);
    }
    free(buf);

    gretl_push_c_numeric_locale();

    got = 0;
    got += gretl_xml_get_prop_as_double(node, "ess", &pmod->ess);
    got += gretl_xml_get_prop_as_double(node, "tss", &pmod->tss);
    got += gretl_xml_get_prop_as_double(node, "sigma", &pmod->sigma);
    got += gretl_xml_get_prop_as_double(node, "rsq", &pmod->rsq);
    got += gretl_xml_get_prop_as_double(node, "adjrsq", &pmod->adjrsq);
    got += gretl_xml_get_prop_as_double(node, "fstt", &pmod->fstt);
    got += gretl_xml_get_prop_as_double(node, "lnL", &pmod->lnL);
    got += gretl_xml_get_prop_as_double(node, "ybar", &pmod->ybar);
    got += gretl_xml_get_prop_as_double(node, "sdy", &pmod->sdy);

    got += gretl_xml_get_prop_as_double(node, "crit0", &pmod->criterion[0]);
    got += gretl_xml_get_prop_as_double(node, "crit1", &pmod->criterion[1]);
    got += gretl_xml_get_prop_as_double(node, "crit2", &pmod->criterion[2]);

    got += gretl_xml_get_prop_as_double(node, "dw", &pmod->dw);
    got += gretl_xml_get_prop_as_double(node, "rho", &pmod->rho);

    if (got < 14) {
	fprintf(stderr, "gretl_model_from_XML: got(2) = %d (expected 14)\n", got);
	*err = E_DATA;
	gretl_pop_c_numeric_locale();
	goto bailout;
    }

    gretl_xml_get_prop_as_string(node, "name", &pmod->name);
    gretl_xml_get_prop_as_string(node, "depvar", &pmod->depvar);

    cur = node->xmlChildrenNode;

    while (cur != NULL && !*err) {
	if (!xmlStrcmp(cur->name, (XUC) "sample")) {
	    gretl_xml_get_prop_as_int(cur, "t1", &pmod->smpl.t1);
	    gretl_xml_get_prop_as_int(cur, "t2", &pmod->smpl.t2);
	} else if (!xmlStrcmp(cur->name, (XUC) "coeff")) {
	    pmod->coeff = gretl_xml_get_double_array(cur, doc, &n, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "sderr")) {
	    pmod->sderr = gretl_xml_get_double_array(cur, doc, &n, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "uhat")) {
	    pmod->uhat = gretl_xml_get_double_array(cur, doc, &n, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "yhat")) {
	    pmod->yhat = gretl_xml_get_double_array(cur, doc, &n, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "xpx")) {
	    pmod->xpx = gretl_xml_get_double_array(cur, doc, &n, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "vcv")) {
	    pmod->vcv = gretl_xml_get_double_array(cur, doc, &n, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "list")) {
	    pmod->list = gretl_xml_node_get_list(cur, doc, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "tests")) {
	    *err = attach_model_tests_from_xml(pmod, cur);
	} else if (!xmlStrcmp(cur->name, (XUC) "params")) {
	    *err = attach_model_params_from_xml(cur, doc, pmod);
	} else if (!xmlStrcmp(cur->name, (XUC) "submask")) {
	    *err = model_submask_from_xml(cur, doc, pmod);
	} else if (!xmlStrcmp(cur->name, (XUC) "missmask")) {
	    if (!gretl_xml_node_get_string(cur, doc, &pmod->missmask)) {
		*err = 1;
	    }
	} else if (!xmlStrcmp(cur->name, (XUC) "arinfo")) {
	    *err = arinfo_from_xml(cur, doc, pmod);
	} else if (!xmlStrcmp(cur->name, (XUC) "data-items")) {
	    *err = model_data_items_from_xml(cur, doc, pmod);
	}
	if (*err) {
	    fprintf(stderr, "gretl_model_from_XML: block 3: err = %d on %s\n", *err, cur->name);
	}
	cur = cur->next;
    }

    gretl_pop_c_numeric_locale();

    /* backward compatibility */
    if (!*err && pmod->nparams > pmod->ncoeff && pmod->depvar == NULL) {
	int i, np = pmod->nparams - 1;

	pmod->depvar = pmod->params[0];
	for (i=0; i<np; i++) {
	    pmod->params[i] = pmod->params[i+1];
	}
	pmod->params[np] = NULL;
	pmod->nparams = np;
    }

 bailout:

    if (*err) {
	if (pmod != NULL) {
	    gretl_model_free(pmod);
	    pmod = NULL;
	}
    }

    return pmod;
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
 * @targ: pointer to target #MODEL.
 * @src: pointer to source #MODEL.
 *
 * Swaps the content of the two model pointers.
 */

void swap_models (MODEL *targ, MODEL *src)
{
    MODEL tmp = *targ;

#if MDEBUG
    fprintf(stderr, "swap_models: %p <-> %p\n", targ, src);
#endif

    *targ = *src;
    *src = tmp;
}

int is_model_cmd (const char *s)
{
    int ret = 0;

    /* FIXME mle? */

    if (s == NULL || *s == '\0') {
	return 0;
    }

    if (!strcmp(s, "ols")  ||
	!strcmp(s, "mpols")  ||
	!strcmp(s, "corc") ||
	!strcmp(s, "hilu") ||
	!strcmp(s, "wls")  ||
	!strcmp(s, "pwe")  ||
	!strcmp(s, "hccm") ||
	!strcmp(s, "heckit")  ||
	!strcmp(s, "hsk")  ||
	!strcmp(s, "add")  ||
	!strcmp(s, "lad")  ||
	!strcmp(s, "omit") ||
	!strcmp(s, "tsls") ||
	!strcmp(s, "logit")  ||
	!strcmp(s, "probit") ||
	!strcmp(s, "tobit") ||
	!strcmp(s, "poisson") ||
	!strcmp(s, "panel") ||
	!strcmp(s, "pooled") ||
	!strcmp(s, "garch") ||
	!strcmp(s, "logistic") ||
	!strcmp(s, "endnls") ||
	!strcmp(s, "arma") ||
	!strcmp(s, "arima") ||
	!strcmp(s, "arbond") ||
	!strcmp(s, "arch") ||
	!strcmp(s, "ar")) {
	ret = 1;
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
 * @opt: option for command to be tested.
 * @model_ci: command index of a gretl model (for example,
 * %OLS, %WLS or %CORC).
 *
 * Returns: 1 if the model-related command in question is
 * meaningful and acceptable in the context of the specific
 * sort of model indentified by @model_ci, otherwise 0.
 */

int command_ok_for_model (int test_ci, gretlopt opt, int model_ci)
{
    int ok = 1;

    if (model_ci == MLE || model_ci == GMM) {
	return 0;
    }

    switch (test_ci) {

    case ADD:
    case ADDTO:
	if (model_ci == NLS || model_ci == ARMA || 
	    model_ci == GARCH || model_ci == HECKIT) {
	    ok = 0;
	}
	break;

    case OMIT:
    case OMITFROM:
	if (model_ci == NLS || model_ci == ARMA || 
	    model_ci == GARCH) {
	    ok = 0;
	}
	break;

    case VIF:
	if (model_ci == NLS || model_ci == TSLS ||
	    model_ci == ARMA || model_ci == GARCH ||
	    model_ci == PANEL || model_ci == ARBOND) {
	    ok = 0;
	}
	break;

    case EQNPRINT:
	if (model_ci == ARMA || model_ci == NLS ||
	    model_ci == ARBOND || model_ci == MLE ||
	    model_ci == GMM) {
	    ok = 0; 
	}
	break;

    case LMTEST:
	if (opt & OPT_H) {
	    ok = (model_ci != ARCH);
	} else if (model_ci != OLS) {
	    ok = 0; /* FIXME */
	}
	break;

    case CHOW:
    case CUSUM:
    case QLRTEST:
    case LEVERAGE:
    case RESET:
	if (model_ci != OLS) ok = 0;
	break;

    case HAUSMAN:
	if (model_ci != OLS) ok = 0;
	break;

    case RESTRICT:
	if (model_ci == LAD || model_ci == NLS || 
	    model_ci == MLE || model_ci == GMM) {
	    ok = 0;
	}
	break;

    case TESTUHAT:
	/* do we really need to exclude garch? */
	if (model_ci == TOBIT || model_ci == PROBIT ||
	    model_ci == LOGIT || model_ci == GARCH) {
	    ok = 0;
	}
	break;

    default:
	break;
    }

    return ok;
}

/**
 * model_test_ok:
 * @ci: index of a model test command.
 * @opt: option associated with test command, if any.
 * @pmod: the model to be tested.
 * @pdinfo: dataset information.
 *
 * A more rigorous version of command_ok_for_model().  Use
 * this function if the extra information is available.
 * 
 * Returns: 1 if the test command @ci (with possible option
 * @opt) is acceptable in the context of the model @pmod, and 
 * the dataset described by @pdinfo, otherwise 0.
 */

int model_test_ok (int ci, gretlopt opt, const MODEL *pmod, 
		   const DATAINFO *pdinfo)
{
    int ok = command_ok_for_model(ci, opt, pmod->ci);

    if (ok && pmod->missmask != NULL) {
	/* can't do these with embedded missing obs */
	if (ci == CUSUM || 
	    (ci == LMTEST && (opt & (OPT_A | OPT_H)))) {
	    ok = 0;
	}
    }

    if (ok && pmod->ncoeff == 1) {
	if (ci == OMIT || ci == OMITFROM || ci == COEFFSUM) {
	    ok = 0;
	} else if (pmod->ifc && ci == LMTEST) {
	    /* const only: rule out squares, logs, White's */
	    if (opt & (OPT_W | OPT_S | OPT_L)) {
		ok = 0;
	    }
	}
    }

    if (ok && !dataset_is_time_series(pdinfo)) {
	/* time-series-only tests */
	if (ci == CHOW || ci == CUSUM || ci == QLRTEST || 
	    (ci == LMTEST && (opt & OPT_H))) {
	    ok = 0;
	}
    }

    if (ok && !dataset_is_time_series(pdinfo) &&
	!dataset_is_panel(pdinfo)) {
	/* time-series or panel tests */
	if (ci == LMTEST && (opt & OPT_A)) {
	    ok = 0;
	}
    }

    if (ok && !dataset_is_panel(pdinfo)) {
	/* panel-only tests */
	if (ci == HAUSMAN || (ci == LMTEST && (opt & OPT_P))) {
	    ok = 0;
	}
    }

    if (ok && pmod->ncoeff - pmod->ifc <= 1 && ci == VIF) {
	/* needs at least two independent vars */
	ok = 0;
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

    if (pmod->ci == MLE || pmod->ci == GMM) {
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

double coeff_pval (int ci, double x, int df)
{
    double p = NADBL;

    if (!xna(x)) {
	if (ASYMPTOTIC_MODEL(ci)) {
	    p = normal_pvalue_2(x);
	} else {
	    p = student_pvalue_2(x, df);
	}
    }

    return p;
}

/**
 * gretl_model_allocate_params:
 * @pmod: pointer to target model.
 * @k: number of strings to allocate.
 * 
 * Allocate an array of @k strings to hold the names given to
 * the associated  coefficients, in a model where these strings 
 * are not simply given by the names of the independent variables.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_model_allocate_params (MODEL *pmod, int k)
{
    pmod->params = strings_array_new_with_length(k, VNAMELEN);
    if (pmod->params == NULL) {
	pmod->errcode = E_ALLOC;
    }

    if (!pmod->errcode) {
	pmod->nparams = k;
    }

    return pmod->errcode;
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
    int np = p + P + q + Q + r + pmod->ifc;
    int xstart;
    int i, j;

    pmod->depvar = gretl_strdup(pdinfo->varname[yno]);
    if (pmod->depvar == NULL) {
	pmod->errcode = E_ALLOC;
	return 1;
    }	

    pmod->params = strings_array_new_with_length(np, VNAMELEN);
    if (pmod->params == NULL) {
	free(pmod->depvar);
	pmod->depvar = NULL;
	pmod->errcode = E_ALLOC;
	return 1;
    }

    pmod->nparams = np;

    if (pmod->ifc) {
	strcpy(pmod->params[0], pdinfo->varname[0]);
	j = 1;
    } else {
	j = 0;
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

    xstart = arma_depvar_pos(pmod) + 1;

    for (i=0; i<r; i++) {
	strcpy(pmod->params[j++], pdinfo->varname[pmod->list[xstart+i]]); 
    }   

    return 0;
}

/**
 * gretl_model_add_panel_varnames:
 * @pmod: pointer to target model.
 * @pdinfo: dataset information.
 * @ulist: list of index numbers of cross-sectional
 * units included in the model.
 * 
 * Composes a set of names to be given to the regressors in an
 * panel model.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_model_add_panel_varnames (MODEL *pmod, const DATAINFO *pdinfo,
				    const int *ulist)
{
    int np = pmod->ncoeff;
    int i, j, v;

    pmod->depvar = gretl_strdup(pdinfo->varname[pmod->list[1]]);
    if (pmod->depvar == NULL) {
	pmod->errcode = E_ALLOC;
	return 1;
    }	

    pmod->params = strings_array_new_with_length(np, VNAMELEN);
    if (pmod->params == NULL) {
	pmod->errcode = E_ALLOC;
	return 1;
    }

    pmod->nparams = np;

    j = 1;
    for (i=0; i<np; i++) {
	v = pmod->list[i+2];
	if (v < pdinfo->v) {
	    strcpy(pmod->params[i], pdinfo->varname[v]);
	} else if (ulist != NULL) {
	    sprintf(pmod->params[i], "ahat_%d", ulist[j++]);
	} else {
	    sprintf(pmod->params[i], "ahat_%d", j++);
	}
    }

    return 0;
}

/**
 * gretl_model_add_allocated_varnames:
 * @pmod: pointer to target model.
 * @vnames: array of names of independent variables.
 * 
 * Attaches an allocated set of variable names to be used
 * when printing model results, for use in special cases
 * where we can't just reference names from the list of
 * regressors attached to the model.  The number of strings
 * must match the number of coefficients, given by the
 * %ncoeff member of @pmod.
 *
 * Note that @pmod "takes charge" of the array @vnames:
 * this will be freed when the model is destroyed.
 */

void gretl_model_add_allocated_varnames (MODEL *pmod, char **vnames)
{
    pmod->nparams = pmod->ncoeff;
    pmod->params = vnames;
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
		ret = gretl_strdup_printf("yformula: %g/(1.0+exp(-(%g%s%g*x)))",
					  lmax, pmod->coeff[0], 
					  (pmod->coeff[1] >= 0)? "+" : "",
					  pmod->coeff[1]);
	    }
	}
    } else if (!pmod->ifc && pmod->ncoeff == 1 && xvar == pmod->list[2]) {
	ret = gretl_strdup_printf("yformula: %g*x", pmod->coeff[0]);
    } else if (pmod->ifc && pmod->ncoeff == 2 && xvar == pmod->list[3]) {
	ret = gretl_strdup_printf("yformula: %g%s%g*x", pmod->coeff[0], 
				  (pmod->coeff[1] >= 0)? "+" : "",
				  pmod->coeff[1]);
    } else if (pmod->ifc && pmod->ncoeff == 3 && xvar == pmod->list[3]) {
	if (model_is_quadratic(pmod, mZ, mdinfo)) {
	    ret = gretl_strdup_printf("yformula: %g%s%g*x%s%g*x**2", pmod->coeff[0], 
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

/**
 * gretl_model_get_scalar:
 * @pmod: pointer to target model.
 * @idx: index for the scalar value that is wanted.
 * @err: location to receive error code (required).
 * 
 * Retrieves a specified scalar statistic from @pmod:
 * @idx must be less than %M_SCALAR_MAX.
 * 
 * Returns: the requested statistic, or #NADBL on failure,
 * in which case @err will contain a non-zero error code.
 */

double gretl_model_get_scalar (const MODEL *pmod, ModelDataIndex idx, 
			       int *err)
{
    double x = NADBL;

    if (pmod == NULL) {
	fprintf(stderr, "model get scalar: model is NULL\n");
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
	x = pmod->sigma;
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
	fprintf(stderr, "model get scalar: x is NA\n");
	*err = E_BADSTAT;
    }

    return x;
}

/**
 * gretl_model_get_series:
 * @pmod: pointer to target model.
 * @idx: index for the series that is wanted.
 * @err: location to receive error code (required).
 * 
 * Retrieves from @pmod a copy of a specified series (for
 * example, regression residuals); @idx must be greater than
 * %M_ELEM_MAX and less than %M_SERIES_MAX.
 * 
 * Returns: the allocated series, or %NULL on failure,
 * in which case @err will contain a non-zero error code.
 */

double *
gretl_model_get_series (const MODEL *pmod, const DATAINFO *pdinfo, 
			ModelDataIndex idx, int *err)
{
    double *x = NULL;
    double *mdata = NULL;
    int t;

    if (pmod->t2 - pmod->t1 + 1 > pdinfo->n || 
	model_sample_problem(pmod, pdinfo)) {
	strcpy(gretl_errmsg, 
	       (idx == M_UHAT)? 
	       _("Can't retrieve uhat: data set has changed") :
	       (idx == M_YHAT)?
	       _("Can't retrieve yhat: data set has changed") :
	       (idx == M_H)?
	       _("Can't retrieve ht: data set has changed") :
	       _("Can't retrieve series: data set has changed"));
	*err = 1;
	return NULL;
    }   

    if ((idx == M_UHAT && pmod->uhat == NULL) ||
	(idx == M_YHAT && pmod->yhat == NULL)) {
	*err = 1;
	return NULL;
    }

    if (idx == M_AHAT) {
	mdata = gretl_model_get_data(pmod, "ahat");
	if (mdata == NULL) {
	    strcpy(gretl_errmsg, _("Can't retrieve intercepts"));
	    *err = 1;
	    return NULL;
	}
    } else if (idx == M_H) {
	mdata = gretl_model_get_data(pmod, "garch_h");
	if (mdata == NULL) {
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
	    } else if (idx == M_AHAT || idx == M_H) {
		x[t] = mdata[t];
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

static gretl_matrix *
model_get_special_vec (const MODEL *pmod, ModelDataIndex idx, int *err)
{
    gretl_matrix *v = NULL;
    double *mdata = NULL;
    int t;

    if (idx == M_AHAT) {
	mdata = gretl_model_get_data(pmod, "ahat");
    } else if (idx == M_H) {
	mdata = gretl_model_get_data(pmod, "garch_h");
    }

    if (mdata == NULL) {
	fprintf(stderr, "model_get_special_vec: mdata is NULL\n");
	*err = E_BADSTAT;
    }

    v = gretl_column_vector_alloc(pmod->t2 - pmod->t1 + 1);
    if (v == NULL) {
	*err = E_ALLOC;
    } else {
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    /* FIXME: is indexation right? */
	    gretl_vector_set(v, t - pmod->t1, mdata[t]);
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

/**
 * gretl_model_get_matrix:
 * @pmod: pointer to target model.
 * @idx: index for the matrix that is wanted.
 * @err: location to receive error code (required).
 * 
 * Retrieves from @pmod a copy of a specified matrix (for
 * example, regression residuals); @idx must be greater than
 * %M_ELEM_MAX and less than %M_MAX.
 * 
 * Returns: the allocated matrix, or %NULL on failure,
 * in which case @err will contain a non-zero error code.
 */

gretl_matrix *gretl_model_get_matrix (MODEL *pmod, ModelDataIndex idx, 
				      int *err)
{
    gretl_matrix *M = NULL;

    if (pmod == NULL) {
	fprintf(stderr, "gretl_model_get_matrix: pmod is NULL\n");
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
    case M_AHAT:
	if (gretl_model_get_data(pmod, "ahat") == NULL) {
	    *err = E_BADSTAT;
	} else {
	    M = model_get_special_vec(pmod, M_AHAT, err);
	}
	break;
    case M_H:
	if (gretl_model_get_data(pmod, "garch_h") == NULL) {
	    *err = E_BADSTAT;
	} else {
	    M = model_get_special_vec(pmod, M_H, err);
	}
	break;
    case M_RHO:
	M = model_get_rhovec(pmod, err);
	break;
    default:
	fprintf(stderr, "gretl_model_get_matrix: got to default\n");
	*err = E_BADSTAT;
	break;
    }

    if (M == NULL && !*err) {
	*err = E_ALLOC;
    }

    return M;
}

static double get_vcv_element (MODEL *pmod, const char *s, 
			       const DATAINFO *pdinfo)
{
    char v1str[VNAMELEN], v2str[VNAMELEN];
    int v1 = 0, v2 = 0;
    int i, j, k, gotit;
    double ret = NADBL;

    if (pmod == NULL || pmod->list == NULL) {
	return NADBL;
    }

    if (sscanf(s, "%15[^,],%15s", v1str, v2str) != 2) {
	return NADBL;
    }

    v1 = gretl_model_get_param_number(pmod, pdinfo, v1str);
    v2 = gretl_model_get_param_number(pmod, pdinfo, v2str);

    if (v1 < 0 || v2 < 0) {
	return NADBL;
    }

    /* make model vcv matrix if need be */
    if (pmod->vcv == NULL && makevcv(pmod, pmod->sigma)) {
	return NADBL;
    }

    /* now find the right entry */
    if (v1 > v2) {
	k = v1;
	v1 = v2;
	v2 = k;
    }

    gotit = 0;
    k = 0;
    for (i=0; i<pmod->ncoeff && !gotit; i++) {
	for (j=0; j<pmod->ncoeff; j++) {
	    if (j < i) {
		continue;
	    }
	    if (i == v1 && j == v2) {
		ret = pmod->vcv[k];
		gotit = 1;
		break;
	    }
	    k++;
	}
    }

    return ret;
}

/* retrieve a specific element from one of the arrays of data
   on a model */

double 
gretl_model_get_data_element (MODEL *pmod, int idx, const char *s,
			      const DATAINFO *pdinfo, int *err)
{
    GretlObjType type;
    double x = NADBL;
    int vi = 0;

    if (pmod == NULL) {
	pmod = get_genr_model(&type);
	if (pmod == NULL || type != GRETL_OBJ_EQN) {
	    pmod = get_last_model(&type);
	    if (pmod == NULL || type != GRETL_OBJ_EQN) {
		*err = E_INVARG;
		return x;
	    }
	}
    }

    /* FIXME 0-based versus 1-based indexing */

    if (idx == M_RHO) {
	if (!(numeric_string(s))) {
	    *err = E_INVARG;
	} else if (dot_atof(s) == 1 && AR1_MODEL(pmod->ci)) {
	    x = gretl_model_get_double(pmod, "rho_in");
	} else if (pmod->ci != AR && dot_atof(s) == 1) {
	    x = pmod->rho;
	} else if (pmod->arinfo == NULL || 
		   pmod->arinfo->arlist == NULL || 
		   pmod->arinfo->rho == NULL) {
	    *err = E_INVARG;
	} else if (!(vi = gretl_list_position(atoi(s), pmod->arinfo->arlist))) {
	    *err = E_INVARG;
	} else {
	    x = pmod->arinfo->rho[vi-1];
	}
    } else if (idx == M_VCV) {
	x = get_vcv_element(pmod, s, pdinfo);
	if (na(x)) {
	    *err = E_INVARG;
	}
    } else if (idx == M_COEFF || idx == M_SE) {
	vi = gretl_model_get_param_number(pmod, pdinfo, s);
	if (vi < 0) {
	    *err = E_INVARG;
	} else {
	    if (idx == M_COEFF && pmod->coeff != NULL) { 
		x = pmod->coeff[vi];
	    } else if (idx == M_SE && pmod->sderr != NULL) {
		x = pmod->sderr[vi];
	    } else {
		*err = E_INVARG;
	    }
	}
    } 

    return x;
}

