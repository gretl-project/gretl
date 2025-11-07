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

#define FULL_XML_HEADERS

#include "libgretl.h"
#include "gretl_xml.h"
#include "matrix_extra.h"
#include "libset.h"
#include "gretl_func.h"
#include "gretl_cmatrix.h"
#include "gretl_panel.h"
#include "gretl_array.h"
#include "qr_estimate.h"

/**
 * SECTION:gretl_model
 * @short_description: handling of the MODEL struct
 * @title: Model structure
 * @include: libgretl.h
 *
 * Provides underlying functionality for gretl's MODEL datatype.
 */

#define MDEBUG 0

typedef void (*DESTFUNC) (void *);

struct model_data_item_ {
    char *key;
    void *ptr;
    int type;
    size_t size;
    DESTFUNC destructor;
};

struct ModelTest_ {
    int type;
    int order;
    char *param;
    unsigned char teststat;
    int dfn;
    double dfd;
    double value;
    double pvalue;
    double crit;
    double alpha;
    gretlopt opt;
};

typedef struct CoeffSep_ CoeffSep;

#define CSLEN 64

struct CoeffSep_ {
    char str[CSLEN];
    int pos;
};

#define PNAMELEN 16 /* for parameter names */

static gretl_bundle *bundlize_test (const ModelTest *src);
static int discard_model_data_item (MODEL *pmod, const char *key,
				    int free_data);

static void gretl_test_init (ModelTest *test, ModelTestType ttype)
{
    test->type = ttype;
    test->order = 0;
    test->param = NULL;
    test->teststat = GRETL_STAT_NONE;
    test->dfn = 0;
    test->dfd = 0;
    test->value = test->pvalue = NADBL;
    test->crit = test->alpha = NADBL;
    test->opt = OPT_NONE;
}

static const char *test_type_key (ModelTestType t)
{
    if (t == GRETL_TEST_ADD) {
	return "add_test";
    } else if (t == GRETL_TEST_ARCH) {
	return "arch_test";
    } else if (t == GRETL_TEST_AUTOCORR) {
	return "autocorr_test";
    } else if (t == GRETL_TEST_CHOW ||
	       t == GRETL_TEST_CHOWDUM) {
	return "chow_test";
    } else if (t == GRETL_TEST_CUSUM) {
	return "cusum_test";
    } else if (t == GRETL_TEST_QLR) {
	return "qlr_test";
    } else if (t == GRETL_TEST_GROUPWISE) {
	return "grpwise_test";
    } else if (t == GRETL_TEST_LOGS) {
	return "logs_test";
    } else if (t == GRETL_TEST_NORMAL) {
	return "normality_test";
    } else if (t == GRETL_TEST_OMIT) {
	return "omit_test";
    } else if (t == GRETL_TEST_RESET) {
	return "reset_test";
    } else if (t == GRETL_TEST_SQUARES) {
	return "squares_test";
    } else if (t == GRETL_TEST_WHITES) {
	return "whites_test";
    } else if (t == GRETL_TEST_SARGAN) {
	return "sargan_test";
    } else if (t == GRETL_TEST_IV_HAUSMAN ||
	       t == GRETL_TEST_PANEL_HAUSMAN) {
	return "hausman_test";
    } else if (t == GRETL_TEST_PANEL_F ||
	       t == GRETL_TEST_PANEL_WELCH) {
	return "fixed_effects_F";
    } else if (t == GRETL_TEST_PANEL_BP ||
	       t == GRETL_TEST_BP) {
	return "bp_test";
    } else if (t == GRETL_TEST_PANEL_TIMEDUM) {
	return "timedum_test";
    } else if (t == GRETL_TEST_HET_1) {
	return "het1_test";
    } else if (t == GRETL_TEST_COMFAC) {
	return "comfac_test";
    } else if (t == GRETL_TEST_INDEP) {
	return "independence_test";
    } else if (t == GRETL_TEST_RE) {
	return "rho_test";
    } else if (t == GRETL_TEST_WITHIN_F) {
	return "within_F";
    } else if (t == GRETL_TEST_RE_WALD) {
	return "re_wald_test";
    } else if (t == GRETL_TEST_XDEPEND) {
	return "cross_sectional_dependence_test";
    } else if (t == GRETL_TEST_PANEL_AR) {
	return "panel_ar_test";
    } else {
	fprintf(stderr, "test_type_key(): type %d has no key!\n", t);
	return NULL;
    }
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
					  GretlType type, size_t size,
					  DESTFUNC destructor)
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
    int err = 0;

    if (item != NULL) {
	item->key = gretl_strdup(orig->key);
	if (item->key == NULL) {
	    free(item);
	    item = NULL;
	}
    }

    if (item != NULL) {
	if (orig->type == GRETL_TYPE_MATRIX) {
	    item->ptr = gretl_matrix_copy(orig->ptr);
	} else if (orig->type == GRETL_TYPE_ARRAY) {
	    item->ptr = gretl_array_copy(orig->ptr, &err);
	} else {
	    item->ptr = malloc(orig->size);
	}
	if (item->ptr == NULL) {
	    free(item->key);
	    free(item);
	    item = NULL;
	}
    }

    if (item != NULL) {
	if (orig->type != GRETL_TYPE_MATRIX &&
	    orig->type != GRETL_TYPE_ARRAY) {
	    memcpy(item->ptr, orig->ptr, orig->size);
	}
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

static int gretl_model_set_data_with_destructor (MODEL *pmod,
                                                 const char *key,
                                                 void *ptr,
                                                 GretlType type,
                                                 size_t size,
                                                 DESTFUNC destructor)
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
 * copied to the model; simply, @ptr is recorded on the model,
 * meaning that @pmod takes ownership of the data. The data
 * pointer will be freed when @pmod is cleared with clear_model().
 * If the data item has structure that requires special treatment
 * on freeing, use a type-specific function instead, for example
 * gretl_model_set_matrix_as_data().
 *
 * The @size is needed just in case the model is copied with
 * copy_model(), in which case the target of the copying
 * operation receives a newly allocated copy of the data in
 * question.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_data (MODEL *pmod, const char *key, void *ptr,
			  GretlType type, size_t size)
{
    return gretl_model_set_data_with_destructor(pmod, key, ptr, type,
						size, NULL);
}

static void matrix_free_callback (void *p)
{
    gretl_matrix_free((gretl_matrix *) p);
}

/**
 * gretl_model_set_matrix_as_data:
 * @pmod: pointer to #MODEL.
 * @key: key string, used in retrieval.
 * @m: matrix to attach.
 *
 * Attaches @m to @pmod as data, recoverable via the key @key
 * using gretl_model_get_data().
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_matrix_as_data (MODEL *pmod, const char *key,
				    gretl_matrix *m)
{
    return gretl_model_set_data_with_destructor(pmod, key, (void *) m,
						GRETL_TYPE_MATRIX, 0,
						matrix_free_callback);
}

static void array_free_callback (void *p)
{
    gretl_array_destroy((gretl_array *) p);
}

/**
 * gretl_model_set_array_as_data:
 * @pmod: pointer to #MODEL.
 * @key: key string, used in retrieval.
 * @A: array to attach.
 *
 * Attaches @A to @pmod as data, recoverable via the key @key
 * using gretl_model_get_data().
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_array_as_data (MODEL *pmod, const char *key,
				   gretl_array *A)
{
    return gretl_model_set_data_with_destructor(pmod, key, (void *) A,
						GRETL_TYPE_ARRAY, 0,
						array_free_callback);
}

/**
 * gretl_model_set_list_as_data:
 * @pmod: pointer to #MODEL.
 * @key: key string, used in retrieval.
 * @list: list to attach.
 *
 * Attaches @list to @pmod as data, recoverable via the key @key
 * using gretl_model_get_list().  Note that the model takes
 * ownership of the supplied list.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_list_as_data (MODEL *pmod, const char *key, int *list)
{
    size_t size = (list[0] + 1) * sizeof *list;

    return gretl_model_set_data_with_destructor(pmod, key, (void *) list,
						GRETL_TYPE_LIST, size,
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
						GRETL_TYPE_STRING, size,
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
    int err = 0;

    /* if value is already set, reset it */
    valp = gretl_model_get_data(pmod, key);
    if (valp != NULL) {
	*valp = val;
	return 0;
    }

    /* otherwise start from scratch */
    valp = malloc(sizeof *valp);
    if (valp == NULL) {
        err = E_ALLOC;
    } else {
        *valp = val;
        err = gretl_model_set_data(pmod, key, valp, GRETL_TYPE_INT,
                                   sizeof(int));
        if (err) {
            free(valp);
        }
    }

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

    err = gretl_model_set_data(pmod, key, valp, GRETL_TYPE_DOUBLE,
			       sizeof(double));
    if (err) free(valp);

    return err;
}

static VCVInfo *vcv_info_new (void)
{
    VCVInfo *vi;

    vi = malloc(sizeof *vi);

    if (vi != NULL) {
	vi->vmaj = vi->vmin = 0;
	vi->order = vi->flags = 0;
	vi->bw = NADBL;
	vi->cv1 = NULL;
	vi->cv2 = NULL;
    }

    return vi;
}

static void vcv_info_free (void *data)
{
    VCVInfo *vi = data;

    if (vi != NULL) {
	free(vi->cv1);
	free(vi->cv2);
	free(vi);
    }
}

/**
 * gretl_model_set_full_vcv_info:
 *
 * Returns: 0 on success, 1 on failure.
 */

static int
gretl_model_set_full_vcv_info (MODEL *pmod, int vmaj, int vmin,
			       int order, int flags, double bw,
			       const char *cv1, const char *cv2)
{
    VCVInfo *vi;
    int prev = 0;
    int err = 0;

    vi = gretl_model_get_data(pmod, "vcv_info");

    if (vi == NULL) {
	vi = vcv_info_new();
	if (vi == NULL) {
	    return E_ALLOC;
	}
    } else {
	prev = 1;
	free(vi->cv1);
	free(vi->cv2);
	vi->cv1 = vi->cv2 = NULL;
    }

    vi->vmaj = vmaj;
    vi->vmin = vmin;
    vi->order = order;
    vi->flags = flags;
    vi->bw = bw;

    if (cv1 != NULL) {
	vi->cv1 = gretl_strdup(cv1);
    }
    if (cv2 != NULL) {
	vi->cv2 = gretl_strdup(cv2);
    }

    if (!prev) {
	err = gretl_model_set_data_with_destructor(pmod, "vcv_info", vi,
						   GRETL_TYPE_STRUCT,
						   sizeof *vi,
						   vcv_info_free);
    }

    return err;
}

/**
 * gretl_model_set_vcv_info:
 * @pmod: pointer to model.
 * @vmaj: top-level VCV type.
 * @vmin: variant under @vmaj, if applicable.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_vcv_info (MODEL *pmod, int vmaj, int vmin)
{
    return gretl_model_set_full_vcv_info(pmod, vmaj, vmin,
					 0, 0, 0, NULL, NULL);
}

/**
 * gretl_model_set_cluster_vcv_info:
 * @pmod: pointer to model.
 * @cv1: name of (first) cluster variable.
 * @cv2: name of second cluster variable or NULL.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_cluster_vcv_info (MODEL *pmod,
				      const char *cv1,
				      const char *cv2)
{
    return gretl_model_set_full_vcv_info(pmod, VCV_CLUSTER,
					 0, 0, 0, 0,
					 cv1, cv2);
}

/**
 * gretl_model_set_hac_vcv_info:
 * @pmod: pointer to model.
 * @kern: kernel type.
 * @order: lag order.
 * @flags: bitflags.
 * @bw: QS bandwidth, or #NADBL if not applicable.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_hac_vcv_info (MODEL *pmod, int kern,
				  int order, int flags,
				  double bw)
{
    return gretl_model_set_full_vcv_info(pmod, VCV_HAC, kern,
					 order, flags, bw,
					 NULL, NULL);
}

/**
 * gretl_model_set_hac_order:
 * @pmod: pointer to model.
 * @order: lag order.
 *
 * Returns: 0 on success, 1 on failure.
 */

int gretl_model_set_hac_order (MODEL *pmod, int order)
{
    VCVInfo *vi = gretl_model_get_data(pmod, "vcv_info");

    if (vi != NULL) {
	vi->order = order;
	return 0;
    } else {
	return E_DATA;
    }
}

/**
 * gretl_model_get_vcv_type:
 * @pmod: pointer to model.
 *
 * Returns: index of variance-covariance matrix type.
 */

int gretl_model_get_vcv_type (const MODEL *pmod)
{
    VCVInfo *vi = gretl_model_get_data(pmod, "vcv_info");

    if (vi != NULL) {
	return vi->vmaj;
    } else {
	return 0;
    }
}

/**
 * gretl_model_get_hc_version:
 * @pmod: pointer to model.
 *
 * Returns: the "HC" (HCCME) variant employed in @pmod,
 * or -1 if this does not apply.
 */

int gretl_model_get_hc_version (const MODEL *pmod)
{
    VCVInfo *vi = gretl_model_get_data(pmod, "vcv_info");

    if (vi != NULL && vi->vmaj == VCV_HC) {
	return vi->vmin;
    } else {
	return -1;
    }
}

/**
 * gretl_model_get_cluster_vname:
 * @pmod: pointer to model.
 *
 * Returns: the name of the (first) clustering variable used for
 * the variance-covariance matrix in @pmod, or NULL if there is no
 * such variable.
 */

const char *gretl_model_get_cluster_vname (const MODEL *pmod)
{
    VCVInfo *vi = gretl_model_get_data(pmod, "vcv_info");

    if (vi != NULL && vi->vmaj == VCV_CLUSTER) {
	return vi->cv1;
    } else {
	return NULL;
    }
}

/**
 * gretl_model_get_cluster_vname2:
 * @pmod: pointer to model.
 *
 * Returns: the name of the second clustering variable used for
 * the variance-covariance matrix in @pmod, or NULL if there is no
 * such variable.
 */

const char *gretl_model_get_cluster_vname2 (const MODEL *pmod)
{
    VCVInfo *vi = gretl_model_get_data(pmod, "vcv_info");

    if (vi != NULL && vi->vmaj == VCV_CLUSTER) {
	return vi->cv2;
    } else {
	return NULL;
    }
}

/**
 * gretl_model_get_data_full:
 * @pmod: pointer to model.
 * @key: key string.
 * @copied: location to receive flag indicating whether the
 * return value is an allocated copy of the original data.
 * @type: location to receive data type.
 * @sz: location to receive the size of the data.
 *
 * Returns: the data pointer identified by @key, or %NULL on failure.
 * If a non-zero value is written to @copied this indicates that the
 * return value is a copy of the original (and therefore it is the
 * caller's responsibility to free the data when it is no longer
 * required).
 */

void *gretl_model_get_data_full (const MODEL *pmod, const char *key,
				 GretlType *type, int *copied,
				 size_t *sz)
{
    void *ret = NULL;
    GretlType itype = 0;
    size_t isize = 0;
    int alloced = 0;
    int i, found = 0;

    if (pmod == NULL) {
	return NULL;
    }

    for (i=0; i<pmod->n_data_items; i++) {
	if (!strcmp(key, pmod->data_items[i]->key)) {
	    ret = pmod->data_items[i]->ptr;
	    itype = pmod->data_items[i]->type;
	    isize = pmod->data_items[i]->size;
	    found = 1;
	    break;
	}
    }

    if (!found && pmod->tests != NULL) {
	const ModelTest *test;
	const char *tkey;
	gretl_bundle *b;

	for (i=0; i<pmod->ntests; i++) {
	    test = &pmod->tests[i];
	    tkey = test_type_key(test->type);
	    if (tkey != NULL && !strcmp(key, tkey)) {
		b = bundlize_test(test);
		if (b != NULL) {
		    ret = b;
		    itype = GRETL_TYPE_BUNDLE;
		    alloced = 1;
		}
		break;
	    }
	}
    }

    if (ret != NULL) {
	if (type != NULL) {
	    *type = itype;
	}
	if (sz != NULL) {
	    *sz = isize;
	}
	if (copied != NULL) {
	    *copied = alloced;
	}
    }

    return ret;
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
    return gretl_model_get_data_full(pmod, key, NULL, NULL, NULL);
}

/**
 * gretl_model_steal_data:
 * @pmod: pointer to model.
 * @key: key string.
 *
 * Returns: the data pointer identified by @key, or %NULL on failure.
 * The caller takes ownership.
 */

void *gretl_model_steal_data (MODEL *pmod, const char *key)
{
    void *ret = NULL;
    int i;

    if (pmod == NULL || key == NULL) {
	return NULL;
    }

    for (i=0; i<pmod->n_data_items; i++) {
	if (!strcmp(key, pmod->data_items[i]->key)) {
	    ret = pmod->data_items[i]->ptr;
            discard_model_data_item(pmod, key, 0);
	    break;
	}
    }

    return ret;
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

    if (pmod == NULL) {
	return 0;
    }

    for (i=0; i<pmod->n_data_items; i++) {
	if (pmod->data_items[i]->type != GRETL_TYPE_INT) {
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

    if (pmod == NULL) {
	return NADBL;
    }

    for (i=0; i<pmod->n_data_items; i++) {
	if (pmod->data_items[i]->type != GRETL_TYPE_DOUBLE) {
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
 * gretl_model_get_double_default:
 * @pmod: pointer to model.
 * @key: key string.
 * @deflt: default value
 *
 * Returns: the double-precision value identified by @key, or
 * @deflt if there is no such value.
 */

double gretl_model_get_double_default (const MODEL *pmod,
				       const char *key,
				       double deflt)
{
    double *valp = NULL;
    double ret = deflt;
    int i;

    if (pmod == NULL) {
	return NADBL;
    }

    for (i=0; i<pmod->n_data_items; i++) {
	if (pmod->data_items[i]->type != GRETL_TYPE_DOUBLE) {
	    continue;
	}
	if (!strcmp(key, pmod->data_items[i]->key)) {
	    valp = (double *) pmod->data_items[i]->ptr;
	    ret = *valp;
	    break;
	}
    }

    return ret;
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

    if (pmod == NULL) {
	return NULL;
    }

    for (i=0; i<pmod->n_data_items; i++) {
	if (pmod->data_items[i]->type != GRETL_TYPE_LIST) {
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

    err = gretl_model_set_data(pmod, "coeffsep", cs, GRETL_TYPE_STRUCT,
			       sizeof *cs);
    if (err) {
	free(cs);
    }

    return err;
}

/**
 * gretl_model_get_coeff_separator:
 * @pmod: pointer to model.
 * @ps: location to receive string, if any, or NULL.
 * @ppos: location to receive position, if any, or NULL.
 *
 * Retrieves information that has been set on @pmod regarding the
 * insertion of an extra line when printing the coefficients, if any.
 *
 * Returns: 1 if such information is present, 0 otherwise.
 */

int gretl_model_get_coeff_separator (const MODEL *pmod, char **ps, int *ppos)
{
    CoeffSep *cs = gretl_model_get_data(pmod, "coeffsep");

    if (cs == NULL) {
	return 0;
    }

    if (ps != NULL) {
	*ps = gretl_strdup(cs->str);
    }
    if (ppos != NULL) {
	*ppos = cs->pos;
    }

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

static void plain_ar_coeff_name (char *targ, const MODEL *pmod, int i)
{
    int k = i - pmod->ncoeff;

    if (pmod->arinfo != NULL && pmod->arinfo->arlist != NULL &&
	k >= 0 && k < pmod->arinfo->arlist[0]) {
	sprintf(targ, "u_%d", pmod->arinfo->arlist[k+1]);
    } else {
	strcpy(targ, "unknown");
    }
}

/**
 * gretl_model_get_param_name:
 * @pmod: pointer to model.
 * @dset: dataset information.
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

char *gretl_model_get_param_name (const MODEL *pmod,
				  const DATASET *dset,
				  int i, char *targ)
{
    *targ = '\0';

    if (pmod != NULL) {
	int j = i + 2;
	int k = -1;

#if 0
	/* should we do this? */
	if (pmod->dataset != NULL) {
	    dset = pmod->dataset;
	}
#endif

	if (pmod->aux == AUX_ARCH) {
	    make_cname(dset->varname[pmod->list[j]], targ);
	} else if (pmod->ci == PANEL && (pmod->opt & OPT_H)) {
	    /* --unit-weights */
	    strcpy(targ, dset->varname[pmod->list[j]]);
	} else if (NONLIST_MODEL(pmod->ci) ||
		   pmod->ci == ARMA || pmod->ci == PANEL ||
		   pmod->ci == DPANEL || pmod->ci == GARCH ||
		   pmod->ci == BIPROBIT) {
	    k = i;
	} else if (pmod->ci == MPOLS && pmod->params != NULL) {
	    k = i;
	} else if ((pmod->ci == PROBIT || pmod->ci == LOGIT ||
		    pmod->ci == HECKIT) && pmod->params != NULL) {
	    k = i;
	} else if (pmod->ci == AR && i >= pmod->ncoeff) {
	    plain_ar_coeff_name(targ, pmod, i);
	} else if (pmod->ci == ARCH && i >= pmod->ncoeff) {
	    sprintf(targ, "alpha(%d)", i - pmod->ncoeff);
	} else if (pmod->list == NULL || j > pmod->list[0] ||
		   pmod->list[j] >= dset->v) {
	    k = i;
	} else {
	    strcpy(targ, dset->varname[pmod->list[j]]);
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

gretl_array *gretl_model_get_param_names (const MODEL *pmod,
					  const DATASET *dset,
					  int *err)
{
    gretl_array *names = NULL;

    if (pmod == NULL) {
	return NULL;
    }

    names = gretl_array_new(GRETL_TYPE_STRINGS,
			    pmod->ncoeff, err);

    if (names != NULL) {
	char targ[VNAMELEN];
	int i;

	for (i=0; i<pmod->ncoeff; i++) {
	    gretl_model_get_param_name(pmod, dset, i, targ);
	    gretl_array_set_string(names, i, targ, 1);
	}
    }

    return names;
}

/**
 * gretl_model_get_param_number:
 * @pmod: pointer to model.
 * @dset: dataset information.
 * @s: name of model parameter.
 *
 * Returns: the zero-based index of the coefficient in @pmod
 * corresponding to @s, or -1 if @s is not the name
 * of a parameter.
 */

int
gretl_model_get_param_number (const MODEL *pmod, const DATASET *dset,
			      const char *s)
{
    char pname[VNAMELEN], tmp[VNAMELEN];
    int i, idx = -1;

    if (pmod == NULL || s == NULL) {
	return -1;
    }

    if (!strcmp(s, "0")) {
	strcpy(pname, "const");
    } else {
	strcpy(pname, s);
    }

    for (i=0; i<pmod->ncoeff; i++) {
	gretl_model_get_param_name(pmod, dset, i, tmp);
	if (!strcmp(pname, tmp)) {
	    idx = i;
	    break;
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
 * reset_coeff_intervals:
 * @cf: pointer to confidence intervals struct.
 * @alpha: nominal non-coverage, as decimal.
 *
 * Recomputes the intervals in @cf using the given value
 * of @alpha.
 *
 * Returns: 0 on success, non-zero if @alpha is out of
 * bounds.
 */

int reset_coeff_intervals (CoeffIntervals *cf, double alpha)
{
    double se, tbak = cf->t;
    int i;

    if (alpha <= 0 || alpha >= 1) {
	return E_DATA;
    }

    if (cf->asy) {
	cf->t = normal_cdf_inverse(1 - alpha / 2);
    } else {
	cf->t = student_cdf_inverse(cf->df, 1 - alpha / 2);
    }

    for (i=0; i<cf->ncoeff; i++) {
	if (cf->maxerr[i] > 0) {
	    se = cf->maxerr[i] / tbak;
	    cf->maxerr[i] = se * cf->t;
	}
    }

    cf->alpha = alpha;

    return 0;
}

/**
 * gretl_model_get_coeff_intervals:
 * @pmod: pointer to gretl model.
 * @dset: dataset information.
 * @opt: TBA.
 *
 * Save the 95 percent confidence intervals for the parameter
 * estimates in @pmod.
 *
 * Returns: pointer to #CONFINT struct containing the results.
 */

CoeffIntervals *
gretl_model_get_coeff_intervals (const MODEL *pmod,
				 const DATASET *dset,
				 gretlopt opt)
{
    CoeffIntervals *cf;
    char pname[32];
    int nc = pmod->ncoeff;
    const double *b = pmod->coeff;
    const double *se = pmod->sderr;
    int offset = 0;
    int i, err = 0;

    cf = malloc(sizeof *cf);
    if (cf == NULL) {
	return NULL;
    }

    if ((opt & OPT_O) && pmod->ifc) {
	/* doing odds ratios: skip intercept */
	offset = 1;
	b++;
	se++;
	nc--;
    }

    cf->ncoeff = nc;
    cf->df = pmod->dfd;
    cf->opt = opt;
    cf->coeff = cf->maxerr = NULL;
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

    cf->alpha = .05;

    if (ASYMPTOTIC_MODEL(pmod->ci)) {
	cf->asy = 1;
	cf->t = normal_cdf_inverse(0.975);
    } else {
	cf->asy = 0;
	cf->t = tcrit95(pmod->dfd);
    }

    for (i=0; i<cf->ncoeff && !err; i++) {
	cf->coeff[i] = b[i];
	cf->maxerr[i] = (se[i] > 0)? cf->t * se[i] : 0.0;
	gretl_model_get_param_name(pmod, dset, i+offset, pname);
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

gretl_matrix *conf_intervals_matrix (CoeffIntervals *cf)
{
    gretl_matrix *ret;
    const double *b = cf->coeff;
    const double *me = cf->maxerr;
    double se, lo, hi, pc;
    char **S = NULL;
    char head[32];
    int odds = (cf->opt & OPT_O);
    int nrows = cf->ncoeff;
    int ncols = 3;
    int se_col = 0;
    int lo_col = 1;
    int hi_col = 2;
    int i, j;

    if (cf->opt & OPT_E) {
	se_col = 1;
	lo_col++;
	hi_col++;
	ncols++;
    }

    ret = gretl_matrix_alloc(nrows, ncols);
    gretl_matrix_set_rownames(ret, strings_array_dup(cf->names, nrows));
    S = strings_array_new(ncols);
    pc = 100 * (1 - cf->alpha);

    for (i=0; i<nrows; i++) {
	if (na(me[i])) {
	    lo = hi = NADBL;
	} else {
	    lo = b[i] - me[i];
	    hi = b[i] + me[i];
	}
	for (j=0; j<ncols; j++) {
	    if (j == 0) {
		/* coeff or odds ratio */
		gretl_matrix_set(ret, i, j, odds ? exp(b[i]) : b[i]);
                if (i == 0) {
                    S[j] = gretl_strdup(odds? "odds ratio" : "coefficient");
                }
	    } else if (j == se_col) {
		/* standard error */
		se = me[i] / cf->t;
		gretl_matrix_set(ret, i, j, odds ? exp(b[i]) * se : se);
                if (i == 0) {
                    S[j] = gretl_strdup("std. error");
                }
	    } else if (j == lo_col) {
		/* lower limit */
		gretl_matrix_set(ret, i, j, odds ? exp(lo) : lo);
                if (i == 0) {
                    sprintf(head, "low%g", pc);
                    S[j] = gretl_strdup(head);
                }
	    } else if (j == hi_col) {
		/* upper limit */
		gretl_matrix_set(ret, i, j, odds ? exp(hi) : hi);
                if (i == 0) {
                    sprintf(head, "high%g", pc);
                    S[j] = gretl_strdup(head);
                }
	    }
	}
    }

    gretl_matrix_set_colnames(ret, S);

    return ret;
}

int gretl_is_simple_OLS (const MODEL *pmod)
{
    if (pmod->ci != OLS || pmod->list == NULL) {
        return 0;
    } else if (pmod->list[0] == 3 && pmod->list[2] == 0) {
        return 1;
    } else if (pmod->list[0] == 2 && pmod->list[2] != 0) {
        return 1;
    } else {
	return 0;
    }
}

int gretl_is_arima_model (const MODEL *pmod)
{
    int d = gretl_model_get_int(pmod, "arima_d");
    int D = gretl_model_get_int(pmod, "arima_D");

    return (d > 0 || D > 0);
}

int gretl_is_between_model (const MODEL *pmod)
{
    if (pmod->ci == PANEL && (pmod->opt & OPT_B)) {
	return 1;
    } else {
	return 0;
    }
}

/* "regular" here means pooled OLS, fixed effects or
   random effects
*/

int gretl_is_regular_panel_model (const MODEL *pmod)
{
    if (pmod->ci == OLS && gretl_model_get_int(pmod, "pooled")) {
	return 1;
    } else if (pmod->ci != PANEL) {
	return 0;
    } else {
	gretlopt opt = pmod->opt & (OPT_P | OPT_F | OPT_U);

	return opt != 0;
    }
}

#define arma_included(m,i) (m == NULL || m[i] == '1')

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

static int arma_included_lags (int k, const char *mask)
{
    int i, nk = k;

    if (mask != NULL) {
	nk = 0;
	for (i=0; i<k; i++) {
	    if (mask[i] == '1') nk++;
	}
    }

    return nk;
}

static int arma_AR_lags (const MODEL *pmod)
{
    const char *pmask = gretl_model_get_data(pmod, "pmask");
    int p = arma_model_nonseasonal_AR_order(pmod);

    return arma_included_lags(p, pmask);
}

static int arma_MA_lags (const MODEL *pmod)
{
    const char *qmask = gretl_model_get_data(pmod, "qmask");
    int q = arma_model_nonseasonal_MA_order(pmod);

    return arma_included_lags(q, qmask);
}

/**
 * arma_model_AR_MA_coeffs:
 * @pmod: pointer to gretl model.
 * @phi_star: pointer to receive AR coeff vector.
 * @theta_star: pointer to receive MA coeff vector.
 * @opt: use OPT_I to get the "integrated" variant, which
 * handles differencing.
 *
 * Creates consolidated versions of the AR and MA coefficient vectors
 * from @pmod.  If @pmod includes seasonal ARMA terms, the vectors are
 * suitably expanded to include the interactions between seasonal and
 * non-seasonal terms.  If OPT_I is given and the dependent variable
 * has been differenced, the AR coefficients are integrated to account
 * for the differencing.  These are the \Phi^* and \Theta^* as used by
 * Box and Jenkins for forecasting.
 *
 * Returns: 0 on success, non-zero on error.
 */

int arma_model_AR_MA_coeffs (const MODEL *pmod,
			     gretl_vector **phi_star,
			     gretl_vector **theta_star,
			     gretlopt opt)
{
    gretl_vector *ac = NULL;
    gretl_vector *mc = NULL;
    int *ainfo;
    int err = 0;

    if (pmod->ci != ARMA) {
	err = E_DATA;
    } else if ((ainfo = gretl_model_get_data(pmod, "ainfo")) == NULL) {
	fprintf(stderr, "AR_MA_coeffs: no 'ainfo' available!\n");
	err = E_DATA;
    } else {
	const double *phi = NULL, *Phi = NULL;
	const double *theta = NULL, *Theta = NULL;
	const char *pmask = gretl_model_get_data(pmod, "pmask");
	const char *qmask = gretl_model_get_data(pmod, "qmask");
	int p  = ainfo[1];
	int q  = ainfo[2];
	int P  = ainfo[3];
	int Q  = ainfo[4];
	int np = ainfo[5];
	int nq = ainfo[6];
	int d  = ainfo[7];
	int D  = ainfo[8];
	int s  = ainfo[9];
	int pmax = p + s * P;
	int pstar = pmax + d + s * D;
	int qmax = q + s * Q;
	double x, y;
	int i, j, k, ii;

	if (pstar > 0) {
	    /* we have some AR terms */
	    ac = gretl_zero_matrix_new(pstar + 1, 1);
	    if (ac == NULL) {
		err = E_ALLOC;
	    }
	}

	if (!err && qmax > 0) {
	    /* we have some MA terms */
	    mc = gretl_zero_matrix_new(qmax + 1, 1);
	    if (mc == NULL) {
		gretl_vector_free(ac);
		ac = NULL;
		err = E_ALLOC;
	    }
	}

	if (!err) {
	    phi = pmod->coeff + pmod->ifc; /* non-seasonal AR coeffs */
	    Phi = phi + np;                /* seasonal AR coeffs */
	    theta = Phi + P;               /* non-seasonal MA coeffs */
	    Theta = theta + nq;            /* seasonal MA coeffs */
	}

	if (ac != NULL) {
	    for (j=0; j<=P; j++) {
		x = (j == 0)? -1 : Phi[j-1];
		k = 0;
		for (i=0; i<=p; i++) {
		    y = (i == 0)? -1 :
			(arma_included(pmask, i-1))? phi[k++] : 0;
		    ii = i + s * j;
		    ac->val[ii] -= x * y;
		}
	    }
	    if ((opt & OPT_I) && (D > 0 || d > 0)) {
		ar_coeff_integrate(ac->val, d, D, s, pmax);
	    }
	}

	if (mc != NULL) {
	    for (j=0; j<=Q; j++) {
		x = (j == 0)? 1 : Theta[j-1];
		k = 0;
		for (i=0; i<=q; i++) {
		    y = (i == 0)? 1 :
			(arma_included(qmask, i-1))? theta[k++] : 0;
		    ii = i + s * j;
		    mc->val[ii] += x * y;
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

/*
  arma_model_spectrum computes

   s2    C(exp(i*omega)) * C(exp(-i*omega))
  ---- * ----------------------------------
  2*pi   A(exp(i*omega)) * A(exp(-i*omega))

  for omega that goes from 0 to pi in T steps. The
  @phi matrix should contain the AR parameters and
  the @theta matrix the MA parameters.
*/

static gretl_matrix *arma_model_spectrum (gretl_matrix *phi,
					  gretl_matrix *theta,
					  double s2, int T,
					  int *err)
{
    gretl_vector *ret = NULL;
    double num, den, nre_t, nim_t, dre_t, dim_t, xt;
    double c, s;
    double ar, ma;
    double scale;
    int nar, nma;
    int i, n, t;

    ret = gretl_matrix_alloc(T, 2);

    if (ret == NULL) {
	*err = E_ALLOC;
	return ret;
    }

    nar = phi != NULL ? phi->rows : 0;
    nma = theta != NULL ? theta->rows : 0;

    n = MAX(nar, nma);
    scale = s2/M_2PI;

    for (t=0; t<T; t++) {
	xt = t * M_PI / (T-1);
	ar = nar > 0 ? gretl_vector_get(phi, 0) : 1.0;
	ma = nma > 0 ? gretl_vector_get(theta, 0) : 1.0;
	nre_t = ma * ma;
	dre_t = ar * ar;
	nim_t = dim_t = 0;

	for (i=1; i<n; i++) {
	    ar = i < nar ? gretl_vector_get(phi, i) : 0.0;
	    ma = i < nma ? gretl_vector_get(theta, i) : 0.0;
	    c = cos(i * xt);
	    s = sin(i * xt);
	    nre_t += c * ma;
	    nim_t += s * ma;
	    dre_t += c * ar;
	    dim_t += s * ar;
	}

	num = nre_t * nre_t + nim_t * nim_t;
	den = dre_t * dre_t + dim_t * dim_t;
	gretl_matrix_set(ret, t, 0, xt);
	gretl_matrix_set(ret, t, 1, scale * num/den);
    }

    return ret;
}

static double arima_diff (const double *x, int t,
			  int *delta, int k, int *err)
{
    double dxt = x[t];
    int i, p;

    if (na(dxt)) {
	*err = E_MISSDATA;
    }

    for (i=0; i<k && !*err; i++) {
	if (delta[i] != 0) {
	    p = t - i - 1;
	    if (p < 0 || na(x[p])) {
		dxt = NADBL;
		*err = E_MISSDATA;
	    } else {
		dxt -= delta[i] * x[p];
	    }
	}
    }

    return dxt;
}

/* get_arma_yvec: here we reconstruct the dependent variable from an
   ARMA model, differencing as we go in the ARIMA case.  If the model
   includes any regressors besides the constant we have to net out
   their effect in order to produce a periodogram that is comparable
   with the ARMA spectrum. A complication is that in the ARIMAX case
   the regressors may also have to be differenced: this is indicated
   by the presence of the "xdiff" flag on the model.
*/

static gretl_vector *get_arma_yvec (const MODEL *pmod,
				    const DATASET *dset,
				    int *err)
{
    gretl_vector *y = NULL;
    const double *beta = NULL;
    double yt, xti;
    int *xlist = NULL;
    int *delta = NULL;
    int yno = gretl_model_get_depvar(pmod);
    int T = pmod->nobs;
    int i, j, t, t1 = pmod->t1;
    int s, d, D, nx = 0;
    int k = 0, xdiff = 0;

    if (yno < 1 || yno >= dset->v) {
	*err = E_DATA;
	return NULL;
    }

    xlist = gretl_model_get_x_list(pmod);

    if (xlist != NULL) {
	for (i=1; i<=xlist[0]; i++) {
	    if (xlist[i] >= dset->v) {
		*err = E_DATA;
		break;
	    } else if (xlist[i] != 0) {
		nx++;
	    }
	}
	if (!*err && nx > 0) {
	    beta = arma_model_get_x_coeffs(pmod);
	}
    }

    if (*err) {
	free(xlist);
	return NULL;
    }

    d = gretl_model_get_int(pmod, "arima_d");
    D = gretl_model_get_int(pmod, "arima_D");

    if (d > 0 || D > 0) {
	delta = arima_delta_coeffs(d, D, dset->pd);
	if (delta == NULL) {
	    *err = E_ALLOC;
	} else {
	    k = d + dset->pd * D;
	    xdiff = gretl_model_get_int(pmod, "xdiff");
	}
    }

    if (!*err) {
	y = gretl_matrix_alloc(T, 1);
	if (y == NULL) {
	    *err = E_ALLOC;
	}
    }

    if (!*err) {
	for (t=0; t<T && !*err; t++) {
	    s = t + t1;
	    if (delta != NULL) {
		yt = arima_diff(dset->Z[yno], s, delta, k, err);
	    } else {
		yt = dset->Z[yno][s];
	    }
	    if (beta != NULL) {
		/* ARMAX: net out X_t\beta */
		for (i=1, j=0; i<=xlist[0]; i++) {
		    if (xlist[i] != 0) {
			if (xdiff) {
			    xti = arima_diff(dset->Z[xlist[i]], s,
					     delta, k, err);
			} else {
			    xti = dset->Z[xlist[i]][s];
			}
			yt -= beta[j++] * xti;
		    }
		}
	    }
	    gretl_vector_set(y, t, yt);
	}
    }

    free(xlist);
    free(delta);

    if (*err) {
	gretl_matrix_free(y);
	y = NULL;
    }

    return y;
}

gretl_matrix *arma_spectrum_plot_data (const MODEL *pmod,
				       const DATASET *dset,
				       int *err)
{
    gretl_matrix *pdata = NULL;
    gretl_matrix *phi = NULL;
    gretl_matrix *theta = NULL;
    gretl_matrix *spec = NULL;
    gretl_matrix *pergm = NULL;
    gretl_matrix *y = NULL;
    int i, T = pmod->nobs;
    int grid = (T-1)/2;
    double s2 = pmod->sigma * pmod->sigma;

    *err = arma_model_AR_MA_coeffs(pmod, &phi, &theta, OPT_NONE);
    if (*err) {
	return NULL;
    }

    /* different sign convention from forecasting */
    if (phi != NULL) {
	gretl_matrix_multiply_by_scalar(phi, -1.0);
    }

    spec = arma_model_spectrum(phi, theta, s2, grid, err);
    if (*err) {
	goto bailout;
    }

    y = get_arma_yvec(pmod, dset, err);

    if (!*err) {
	pergm = gretl_matrix_fft(y, err);
    }

    if (!*err) {
	pdata = gretl_matrix_alloc(grid, 4);
	if (pdata == NULL) {
	    *err = E_ALLOC;
	} else {
	    double complex z;

	    for (i=0; i<grid; i++) {
		gretl_matrix_set(pdata, i, 0, gretl_matrix_get(spec, i, 0));
		gretl_matrix_set(pdata, i, 1, gretl_matrix_get(spec, i, 1));
		z = gretl_cmatrix_get(pergm, i+1, 0);
		gretl_matrix_set(pdata, i, 2, creal(z));
		gretl_matrix_set(pdata, i, 3, cimag(z));
	    }
	}
    }

 bailout:

    gretl_matrix_free(phi);
    gretl_matrix_free(theta);
    gretl_matrix_free(spec);
    gretl_matrix_free(pergm);
    gretl_matrix_free(y);

    return pdata;
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
    const char *pmask = gretl_model_get_data(pmod, "pmask");
    int p = arma_model_nonseasonal_AR_order(pmod);
    int np = arma_included_lags(p, pmask);
    int P = gretl_model_get_int(pmod, "arma_P");
    int s = gretl_model_get_int(pmod, "arma_pd");
    const double *phi = NULL, *Phi = NULL;
    double *ac = NULL;
    double x, y;
    int i, j, k, ii, pmax;
    int err = 0;

    pmax = p + s * P;

#if 0
    fprintf(stderr, "regarma_model_AR_coeffs: pmax = %d\n", pmax);
#endif

    if (pmax == 0) {
	*pp = 0;
	return 0;
    }

    ac = malloc((pmax + 1) * sizeof *ac);
    if (ac == NULL) {
	return E_ALLOC;
    }

    phi = pmod->coeff + pmod->ifc;
    Phi = phi + np;

    for (i=0; i<=pmax; i++) {
	ac[i] = 0.0;
    }

    for (j=0; j<=P; j++) {
	x = (j == 0)? -1 : Phi[j-1];
	k = 0;
	for (i=0; i<=p; i++) {
	    y = (i == 0)? -1 :
		(arma_included(pmask, i-1))? phi[k++] : 0;
	    ii = i + s * j;
	    ac[ii] -= x * y;
	}
    }

    *phi0 = ac;
    *pp = pmax;

    return err;
}

/**
 * arima_delta_coeffs:
 * @d: order of non-seasonal differencing (<= 2)
 * @D: order of seasonal differencing (<= 2)
 * @s: seasonal periodicity
 *
 * Returns: array of d + s * D coefficients of the lag operator in
 * the expansion of (1-L)^d * (1-L^s)^D. These are given in the
 * negative; for example, if d = 1 then c[0] = 1.
 */

int *arima_delta_coeffs (int d, int D, int s)
{
    int i, k = d + s * D;
    int *c = malloc(k * sizeof *c);

    if (c == NULL) {
	return NULL;
    }

    for (i=0; i<k; i++) {
	c[i] = 0;
    }

    if (d == 1) {
	c[0] = 1;
    } else if (d == 2) {
	c[0] = 2;
	c[1] = -1;
    }

    if (D > 0) {
	c[s-1] += 1;
	if (d > 0) {
	    c[s] -= 1;
	}
	if (d == 2) {
	    c[s] -= 1;
	    c[s+1] += 1;
	}
	if (D == 2) {
	    c[s-1] += 1;
	    c[2*s-1] -= 1;
	    if (d > 0) {
		c[s] -= 1;
		c[2*s] += 1;
	    }
	    if (d == 2) {
		c[s] -= 1;
		c[2*s] += 1;
		c[s+1] += 1;
		c[2*s+1] -= 1;
	    }
	}
    }

    return c;
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
	xc += arma_AR_lags(pmod);
	xc += arma_MA_lags(pmod);
	xc += gretl_model_get_int(pmod, "arma_P");
	xc += gretl_model_get_int(pmod, "arma_Q");
    }

    return xc;
}

/**
 * arma_model_get_n_arma_coeffs:
 * @pmod: pointer to gretl model.
 *
 * Returns: the sum of the numbers of AR and MA coefficients in
 * @pmod.
 */

int arma_model_get_n_arma_coeffs (const MODEL *pmod)
{
    int npq = 0;

    if (pmod->ci == ARMA) {
	npq += arma_AR_lags(pmod);
	npq += arma_MA_lags(pmod);
	npq += gretl_model_get_int(pmod, "arma_P");
	npq += gretl_model_get_int(pmod, "arma_Q");
    }

    return npq;
}

/**
 * gretl_model_ahat_vec:
 * @pmod: pointer to gretl model.
 * @err: location to receive error code.
 *
 * If @pmod was estimated on panel data via fixed or random
 * effects, provides a column vector holding the "ahat" values
 * or estimated individual effects for the individuals
 * included in the active dataset when the model was estimated.
 * Some of these values may be NA if some units were not
 * actually included in estimation.
 *
 * Returns: allocated column vector on success or NULL on failure.
 */

gretl_matrix *gretl_model_ahat_vec (const MODEL *pmod,
				    int *err)
{
    gretl_vector *a_vec = NULL;
    double *ahat;
    int i, s, t, N, T;

    ahat = gretl_model_get_data(pmod, "ahat");
    T = gretl_model_get_int(pmod, "panel_T");

    if (ahat == NULL || T == 0) {
	*err = E_BADSTAT;
	return NULL;
    }

    /* number of units included in model dataset */
    N = pmod->full_n / T;

    a_vec = gretl_column_vector_alloc(N);
    if (a_vec == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    for (i=0; i<N; i++) {
	a_vec->val[i] = NADBL;
	for (t=0; t<T; t++) {
	    s = i * T + t;
	    if (!na(ahat[s])) {
		/* pick up first non-missing value */
		a_vec->val[i] = ahat[s];
		break;
	    }
	}
    }

    return a_vec;
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
    int sep1 = 0, sep2 = 0;
    int seasonal = 0;
    int arima = 0;
    int i, dvpos;

    for (i=1; i<pmod->list[0]; i++) {
	if (pmod->list[i] == LISTSEP) {
	    if (sep1 == 0) {
		sep1 = i;
	    } else {
		sep2 = i;
	    }
	}
    }

    if (sep2) {
	dvpos = sep2 + 1;
    } else if (sep1) {
	dvpos = sep1 + 1;
    } else {
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
	} else if (pmod->ci == DPANEL) {
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
 * @dset: dataset information.
 *
 * Returns: the name of the dependent variable in @pmod.
 */

const char *gretl_model_get_depvar_name (const MODEL *pmod,
					 const DATASET *dset)
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
	    } else if (pmod->ci == DPANEL) {
		dv = arbond_get_depvar(pmod);
	    } else {
		dv = pmod->list[1];
	    }
	}
    }

    return dset->varname[dv];
}

static int ordered_model (const MODEL *pmod)
{
    if (pmod->ci == PROBIT || pmod->ci == LOGIT) {
	if (gretl_model_get_int(pmod, "ordered")) {
	    return 1;
	}
    }

    return 0;
}

static int longlist_model (const MODEL *pmod)
{
    if (pmod->ci == LOGIT ||
	pmod->ci == NEGBIN ||
	pmod->ci == DURATION ||
	pmod->ci == PANEL) {
	return 1;
    } else if (pmod->ci == PROBIT && (pmod->opt & OPT_E)) {
	/* random-effects probit */
	return 1;
    } else {
	return 0;
    }
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

    if (pmod == NULL || pmod->errcode || pmod->list == NULL) {
	return NULL;
    }

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
    } else if (pmod->ci == DPANEL) {
	int ifc = gretl_model_get_int(pmod, "ifc");
	int vi, sep = 0;

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
		vi = pmod->list[i];
		if (vi > 0 || ifc) {
		    /* don't include const if it was dropped */
		    list = gretl_list_append_term(&list, vi);
		    if (list == NULL) {
			return NULL;
		    }
		}
	    }
	}
    } else if (pmod->ci == BIPROBIT) {
	/* not sure what we ought to do here, but for now we'll
	   return the list of regressors for the first equation
	*/
	nx = 0;
	for (i=3; i<=pmod->list[0]; i++) {
	    if (pmod->list[i] == LISTSEP) {
		nx = i - 3;
		break;
	    }
	}
	if (nx > 0) {
	    list = gretl_list_new(nx);
	    if (list != NULL) {
		for (i=1; i<=nx; i++) {
		    list[i] = pmod->list[i+2];
		}
	    }
	}
    } else if (pmod->ci == MIDASREG) {
	int *lfx = gretl_model_get_list(pmod, "lfxlist");

	if (lfx != NULL) {
	    list = gretl_list_copy(lfx);
	}
    } else if (ordered_model(pmod)) {
	nx = pmod->list[0] - 1;
	list = gretl_list_new(nx);
	if (list != NULL) {
	    for (i=1; i<=list[0]; i++) {
		list[i] = pmod->list[i + 1];
	    }
	}
    } else if (!NONLIST_MODEL(pmod->ci)) {
	if (pmod->ci == HECKIT) {
	    nx = gretl_model_get_int(pmod, "base-coeffs");
	} else if (longlist_model(pmod)) {
	    /* models in which the array of coefficients
	       is (or may be) longer than the list of
	       regressors
	    */
	    nx = pmod->list[0] - 1;
	} else {
	    nx = pmod->ncoeff;
	}
	if (nx > 0) {
	    if (nx > pmod->list[0] - 1) {
		/* don't read off the end of pmod->list! */
		nx = pmod->list[0] - 1;
	    }
	    list = gretl_list_new(nx);
	    if (list != NULL) {
		for (i=1; i<=list[0]; i++) {
		    list[i] = pmod->list[i + 1];
		}
	    }
	}
    }

    return list;
}

/**
 * gretl_model_get_secondary_list:
 * @pmod: model to examine.
 *
 * Retrieve an allocated copy of the secondary list from
 * @pmod: e.g. the list of instruments in the case of %IVREG, or
 * the selection equation list for %HECKIT.
 *
 * Returns: allocated list or %NULL on error.
 */

int *gretl_model_get_secondary_list (const MODEL *pmod)
{
    int *list = NULL;
    int pos;

    pos = gretl_list_separator_position(pmod->list);

    if (pos > 0) {
	list = gretl_list_copy_from_pos(pmod->list, pos + 1);
    }

    return list;
}

/**
 * gretl_model_get_y_list:
 * @pmod: model to examine.
 *
 * Retrieve an allocated copy of the list of dependent variables
 * for @pmod: in almost all cases this will have a single
 * element; an exception is biprobit.
 *
 * Returns: allocated list or %NULL on error.
 */

int *gretl_model_get_y_list (const MODEL *pmod)
{
    int *list = NULL;

    if (pmod->list == NULL) {
	return NULL;
    }

    if (pmod->ci == BIPROBIT) {
	list = gretl_list_new(2);
	if (list != NULL) {
	    list[1] = pmod->list[1];
	    list[2] = pmod->list[2];
	}
    } else {
	list = gretl_list_new(1);
	if (list != NULL) {
	    list[1] = gretl_model_get_depvar(pmod);
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
	int i;

	pmod->coeff = malloc(k * sizeof *pmod->coeff);
	if (pmod->coeff == NULL) {
	    return E_ALLOC;
	}
	pmod->sderr = malloc(k * sizeof *pmod->sderr);
	if (pmod->sderr == NULL) {
	    return E_ALLOC;
	}
	for (i=0; i<k; i++) {
	    pmod->coeff[i] = pmod->sderr[i] = NADBL;
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

static double vcv_get_se (double vii, int restricted)
{
    if (restricted) {
	vii = fabs(vii) < 1.0e-17 ? 0.0 : vii;
    }

    return (na(vii) || vii < 0)? NADBL : sqrt(vii);
}

/**
 * gretl_model_write_vcv:
 * @pmod: pointer to model.
 * @V: covariance matrix.
 *
 * Write the covariance matrix @V into the model @pmod, using the
 * special packed format that is required by the MODEL struct,
 * and set the standard errors to the square root of the diagonal
 * elements of this matrix.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_model_write_vcv (MODEL *pmod, const gretl_matrix *V)
{
    int i, j, k, n;
    double x, *tmp;
    int err = 0;

    if (gretl_is_null_matrix(V)) {
	return 0; /* no-op */
    }

    if (V->cols != V->rows) {
	return E_NONCONF;
    }

    k = V->rows;
    n = (k * k + k) / 2;

    /* reallocate model vcv in case it's wrongly sized */
    tmp = realloc(pmod->vcv, n * sizeof *tmp);
    if (tmp == NULL) {
	err = E_ALLOC;
    } else {
	pmod->vcv = tmp;
    }

    /* same for standard errors array */
    if (!err) {
	tmp = realloc(pmod->sderr, k * sizeof *tmp);
	if (tmp == NULL) {
	    err = E_ALLOC;
	} else {
	    pmod->sderr = tmp;
	}
    }

    if (!err) {
	int restricted, idx = 0;

	restricted = gretl_model_get_int(pmod, "restricted");

	for (i=0; i<k; i++) {
	    for (j=i; j<k; j++) {
		x = gretl_matrix_get(V, i, j);
		pmod->vcv[idx++] = x;
		if (i == j) {
		    pmod->sderr[i] = vcv_get_se(x, restricted);
		}
	    }
	}
    }

    return err;
}

static gretl_matrix *make_cluster_vals (MODEL *pmod, int cvar,
					const DATASET *dset,
					int *err)
{
    gretl_matrix *cvals = NULL;
    double ct, *cdata;
    int s, t;

    cdata = malloc(pmod->nobs * sizeof *cdata);

    if (cdata == NULL) {
	*err = E_ALLOC;
    } else {
	s = 0;
	for (t=pmod->t1; t<=pmod->t2 && !*err; t++) {
	    if (!model_missing(pmod, t)) {
		ct = dset->Z[cvar][t];
		if (na(ct)) {
		    *err = E_MISSDATA;
		} else {
		    cdata[s++] = ct;
		}
	    }
	}
    }

    if (!*err) {
	cvals = gretl_matrix_values(cdata, pmod->nobs, OPT_S, err);
	if (!*err && gretl_vector_get_length(cvals) < 2) {
	    gretl_matrix_free(cvals);
	    cvals = NULL;
	    *err = E_DATA;
	}
    }

    free(cdata);

    return cvals;
}

static int model_make_clustered_GG (MODEL *pmod, int ci,
				    const gretl_matrix *G,
				    gretl_matrix *GG,
				    const DATASET *dset,
				    int *pcvar,
				    int *pnc)
{
    gretl_matrix *cvals = NULL;
    gretl_matrix *Gi = NULL;
    const double *cZ = NULL;
    const char *cname;
    int cvar, n_c = 0;
    int k = G->cols;
    int i, j, t;
    int err = 0;

    cname = get_optval_string(ci, OPT_C);
    if (cname == NULL) {
	return E_PARSE;
    }

    cvar = current_series_index(dset, cname);

    if (cvar < 1 || cvar >= dset->v) {
	err = E_UNKVAR;
    } else {
	cvals = make_cluster_vals(pmod, cvar, dset, &err);
    }

    if (!err) {
	Gi = gretl_matrix_alloc(1, k);
	if (Gi == NULL) {
	    gretl_matrix_free(cvals);
	    err = E_ALLOC;
	}
    }

    if (err) {
	return err;
    }

    n_c = gretl_vector_get_length(cvals);
    cZ = dset->Z[cvar];

    for (i=0; i<n_c; i++) {
	double cvi = cvals->val[i];
	int s = 0;

	gretl_matrix_zero(Gi);

	for (t=pmod->t1; t<=pmod->t2; t++) {
	    if (!model_missing(pmod, t)) {
		if (cZ[t] == cvi) {
		    for (j=0; j<k; j++) {
			Gi->val[j] += gretl_matrix_get(G, s, j);
		    }
		}
		s++;
	    }
	}

	gretl_matrix_multiply_mod(Gi, GRETL_MOD_TRANSPOSE,
				  Gi, GRETL_MOD_NONE,
				  GG, GRETL_MOD_CUMULATE);
    }

    if (!err) {
	*pcvar = cvar;
	*pnc = n_c;
    }

    gretl_matrix_free(cvals);
    gretl_matrix_free(Gi);

    return err;
}

static int do_biprobit_adjustment (MODEL *pmod,
				    gretl_matrix *V)
{
    void (*biprobit_adj) (MODEL *, gretl_matrix *);

    biprobit_adj = get_plugin_function("biprobit_adjust_vcv");

    if (biprobit_adj == NULL) {
	return E_FOPEN;
    } else {
	(*biprobit_adj) (pmod, V);
	return 0;
    }
}

/**
 * gretl_model_add_QML_vcv:
 * @pmod: pointer to model.
 * @ci: command index for model.
 * @H: inverse of the (negative) Hessian, k x k.
 * @G: score matrix, T x k.
 * @dset: pointer to dataset (can be NULL if not doing
 * clustering).
 * @opt: may include OPT_C for cluster-robust variant.
 * @pV: location to receive the full covariance matrix, or NULL.
 *
 * Write a QML covariance matrix into the model @pmod, and set
 * the standard errors to the square root of the diagonal
 * elements of this matrix. The @ci argument, specifying the
 * estimator for @pmod, is required only if @opt includes
 * OPT_C; otherwise it is ignored.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_model_add_QML_vcv (MODEL *pmod, int ci,
			     const gretl_matrix *H,
			     const gretl_matrix *G,
			     const DATASET *dset,
			     gretlopt opt,
			     gretl_matrix **pV)
{
    gretl_matrix *GG = NULL;
    gretl_matrix *V = NULL;
    int cvar = 0, n_c = 0;
    int k = H->rows;
    int err = 0;

    V = gretl_matrix_alloc(k, k);

    if (V == NULL) {
	return E_ALLOC;
    }

    if (!err) {
	if (opt & OPT_C) {
	    /* clustered */
	    GG = gretl_zero_matrix_new(k, k);
	    if (GG == NULL) {
		err = E_ALLOC;
	    } else {
		err = model_make_clustered_GG(pmod, ci, G, GG,
					      dset, &cvar, &n_c);
	    }
	} else if (opt & OPT_N) {
	    /* Newey-West */
	    GG = newey_west_OPG(G, &err);
	} else {
	    /* regular QML using OPG */
	    GG = gretl_matrix_XTX_new(G);
	    if (GG == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    if (!err) {
	err = gretl_matrix_qform(H, GRETL_MOD_NONE, GG,
				 V, GRETL_MOD_NONE);
	if (!err && (opt & OPT_C)) {
	    /* clustering: use stata-style df adjustment */
	    double dfc = n_c / (n_c - 1.0);

	    gretl_matrix_multiply_by_scalar(V, dfc);
	}
    }

    if (!err) {
	if (ci == BIPROBIT) {
	    do_biprobit_adjustment(pmod, V);
	}
	err = gretl_model_write_vcv(pmod, V);
    }

    if (!err) {
	if (opt & OPT_C) {
	    const char *cname = dset->varname[cvar];

	    gretl_model_set_int(pmod, "n_clusters", n_c);
	    gretl_model_set_cluster_vcv_info(pmod, cname, NULL);
	    pmod->opt |= OPT_C;
	} else {
	    MLVCVType vt = (opt & OPT_N)? ML_HAC : ML_QML;

	    gretl_model_set_vcv_info(pmod, VCV_ML, vt);
	}
	pmod->opt |= OPT_R;
    }

    gretl_matrix_free(GG);

    if (pV != NULL) {
	*pV = V;
    } else {
	gretl_matrix_free(V);
    }

    return err;
}

/**
 * gretl_model_add_hessian_vcv:
 * @pmod: pointer to model.
 * @H: inverse of the (negative) Hessian.
 *
 * Write @H into the model @pmod as its covariance matrix, and
 * set the standard errors to the square roots of the diagonal
 * elements of this matrix.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_model_add_hessian_vcv (MODEL *pmod,
				 const gretl_matrix *H)
{
    int err = gretl_model_write_vcv(pmod, H);

    if (!err) {
	gretl_model_set_vcv_info(pmod, VCV_ML, ML_HESSIAN);
    }

    return err;
}

/**
 * gretl_model_add_OPG_vcv:
 * @pmod: pointer to model.
 * @G: T x k gradient matrix.
 * @pV: location to receive the full covariance matrix, or NULL.
 *
 * Compute (G'G)^{-1}, write this into @pmod as its covariance matrix,
 * and set the standard errors to the square roots of the diagonal
 * elements of this matrix.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_model_add_OPG_vcv (MODEL *pmod,
			     const gretl_matrix *G,
			     gretl_matrix **pV)
{
    gretl_matrix *GG = gretl_matrix_XTX_new(G);
    int err = 0;

    if (GG == NULL) {
	return E_ALLOC;
    }

    err = gretl_invert_symmetric_matrix(GG);

    if (!err) {
	err = gretl_model_write_vcv(pmod, GG);
	if (!err) {
	    gretl_model_set_vcv_info(pmod, VCV_ML, ML_OP);
	}
    } else {
	/* diagnostic: recreate and print G'G */
	gretl_matrix_multiply_mod(G, GRETL_MOD_TRANSPOSE,
				  G, GRETL_MOD_NONE,
				  GG, GRETL_MOD_NONE);
	gretl_matrix_print(GG, "non-p.d. G'G");
    }

    if (pV != NULL) {
	*pV = GG;
    } else {
	gretl_matrix_free(GG);
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
 * @dset: dataset information.
 *
 * Supplies the caller with a copy of the variance-covariance
 * matrix for the parameter estimates in @pmod, in a format
 * suitable for printing.  See also free_vcv().  To get the
 * covariance matrix as a gretl_matrix, see
 * gretl_vcv_matrix_from_model().
 *
 * Returns: #VMatrix struct or %NULL on error.
 */

VMatrix *gretl_model_get_vcv (MODEL *pmod, const DATASET *dset)
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
	gretl_model_get_param_name(pmod, dset, i, varname);
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
 * gretl_model_get_vcv_element:
 * @pmod: pointer to model.
 * @i: row (0-based).
 * @j: column(0-based).
 * @np: the number of parameters represented in @pmod's
 * covariance matrix.
 *
 * Returns: the (@i, @j) element of @pmod's covariance matrix,
 * or #NADBL on failure.
 */

double gretl_model_get_vcv_element (const MODEL *pmod,
				    int i, int j,
				    int np)
{
    if (pmod->vcv == NULL) {
	return NADBL;
    } else {
	int k = ijton(i, j, np);

	return pmod->vcv[k];
    }
}

/**
 * gretl_model_write_coeffs:
 * @pmod: pointer to model.
 * @b: array of coefficients.
 * @k: number of elements in @b.
 *
 * Write the coefficients @b into the model @pmod, whose
 * coefficient array is resized appropriately if need be.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_model_write_coeffs (MODEL *pmod, double *b, int k)
{
    size_t sz = k * sizeof(double);
    int err = 0;

    if (pmod->coeff == NULL || pmod->ncoeff != k) {
	double *tmp = realloc(pmod->coeff, sz);

	if (tmp == NULL) {
	    err = E_ALLOC;
	} else {
	    pmod->coeff = tmp;
	}
    }

    if (!err) {
	memcpy(pmod->coeff, b, sz);
	pmod->ncoeff = k;
    }

    return err;
}

/**
 * impose_model_smpl:
 * @pmod: pointer to model.
 * @dset: dataset information.
 *
 * Sets on @dset the sample range (starting and ending
 * observations) that was in effect when @pmod was estimated.
 * This is not always the same as the data range over which
 * @pmod was actually estimated (e.g. in case of
 * autoregressive models, where observations are dropped
 * to allow for lags).
 */

void impose_model_smpl (const MODEL *pmod, DATASET *dset)
{
    dset->t1 = pmod->smpl.t1;
    dset->t2 = pmod->smpl.t2;
#if MDEBUG
    fprintf(stderr, "impose_model_smpl: set t1=%d, t2=%d\n",
	    dset->t1, dset->t2);
#endif
}

#ifdef HAVE_X12A

static void maybe_delete_x12_file (const MODEL *pmod)
{
    char *fname = NULL;

    fname = gretl_model_get_data(pmod, "x12a_output");
    if (fname != NULL) {
	gretl_remove(fname);
    }
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
	int n_items = pmod->n_data_items - 1;

	if (n_items == 0) {
	    free(pmod->data_items);
	    pmod->data_items = NULL;
	} else {
	    model_data_item **items;

	    for (i=targ; i<n_items; i++) {
		pmod->data_items[i] = pmod->data_items[i+1];
	    }
	    items = realloc(pmod->data_items, n_items * sizeof *items);
	    if (items != NULL) {
		pmod->data_items = items;
	    }
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

static void model_stats_init (MODEL *pmod)
{
    int i;

    pmod->ess = pmod->tss = NADBL;
    pmod->sigma = NADBL;
    pmod->rsq = pmod->adjrsq = NADBL;
    pmod->fstt = pmod->chisq = NADBL;
    pmod->lnL = NADBL;
    pmod->ybar = pmod->sdy = NADBL;
    pmod->dw = pmod->rho = NADBL;

    for (i=0; i<C_MAX; i++) {
	pmod->criterion[i] = NADBL;
    }
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
 * @dset: pointer to dataset.
 *
 * Initializes a gretl #MODEL, including setting its pointer members
 * to %NULL. This initialization should be done if the caller has
 * declared a #MODEL struct directly, rather than obtaining a pointer to
 * #MODEL using gretl_model_new() (in which case the initialization is
 * done automatically).
 */

void gretl_model_init (MODEL *pmod, const DATASET *dset)
{
    if (pmod == NULL) return;

#if MDEBUG
    fprintf(stderr, "gretl_model_init: pmod at %p\n", (void *) pmod);
#endif

    pmod->ID = 0;
    pmod->refcount = 0;
    pmod->ci = 0;
    pmod->opt = OPT_NONE;
    pmod->full_n = 0;
    pmod->t1 = pmod->t2 = 0;
    pmod->nobs = 0;

    if (dset != NULL) {
	pmod->smpl.t1 = dset->t1;
	pmod->smpl.t2 = dset->t2;
	pmod->smpl.rseed = dset->rseed;
    } else {
	pmod->smpl.t1 = 0;
	pmod->smpl.t2 = 0;
	pmod->smpl.rseed = 0;
    }

    pmod->ncoeff = 0;
    pmod->ntests = 0;
    pmod->nparams = 0;
    pmod->errcode = 0;
    pmod->ifc = 0;
    pmod->nwt = 0;
    pmod->aux = AUX_NONE;
    pmod->esttime = 0;

    model_stats_init(pmod);

    gretl_model_init_pointers(pmod);
    pmod->n_data_items = 0;
}

/**
 * gretl_model_smpl_init:
 * @pmod: pointer to model.
 * @dset: dataset information.
 *
 * Records the start and end of the current sample range in
 * the model @pmod, which may be necessary for future reference
 * if a hypothesis test is to be performed.  Note that this
 * sample range may not be the same as the data range over
 * which the model is actually estimated (for example, in the
 * case of autoregressive models where observations have to
 * be dropped to allow for lags).
 */

void gretl_model_smpl_init (MODEL *pmod, const DATASET *dset)
{
    pmod->smpl.t1 = dset->t1;
    pmod->smpl.t2 = dset->t2;
    pmod->smpl.rseed = dset->rseed;
#if MDEBUG
    fprintf(stderr, "gretl_model_smpl_init: set t1=%d, t2=%d\n",
	    dset->t1, dset->t2);
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
	gretl_model_init(pmod, NULL);
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
 * allocate_working_model:
 *
 * Allocates memory for gretl #MODEL struct and initializes it.
 * The model is "protected" against deletion.
 *
 * Returns: pointer to model (or %NULL if allocation fails).
 */

MODEL *allocate_working_model (void)
{
    MODEL *model = gretl_model_new();

    if (model != NULL) {
	int err = gretl_model_protect(model);

	if (err) {
	    gretl_model_free(model);
	    model = NULL;
	}
    }

    return model;
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

void destroy_working_model (MODEL *model)
{
    if (model != NULL) {
	gretl_model_free_on_exit(model);
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
	    " pmod->depvar = %p\n"
	    " pmod->params = %p\n"
	    " pmod->arinfo = %p\n"
	    " pmod->tests = %p\n"
	    " pmod->data_items = %p\n", msg,
	    (void *) pmod, (void *) pmod->list,
	    (void *) pmod->submask, (void *) pmod->missmask,
	    (void *) pmod->coeff, (void *) pmod->sderr,
	    (void *) pmod->yhat, (void *) pmod->uhat,
	    (void *) pmod->xpx,
	    (void *) pmod->vcv, (void *) pmod->name,
	    (void *) pmod->depvar, (void *) pmod->params,
	    (void *) pmod->arinfo, (void *) pmod->tests,
	    (void *) pmod->data_items);
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
 * declared directly by the caller; in that case the caller should
 * pass the address of the #MODEL in question.
 */

void clear_model (MODEL *pmod)
{
    if (pmod != NULL) {
#if MDEBUG
	fprintf(stderr, "clear model: model at %p\n", (void *) pmod);
#endif

#if MDEBUG > 1
	debug_print_model_info(pmod, "Doing clear_model");
#endif
	if (pmod->list != NULL) free(pmod->list);
	if (pmod->missmask != NULL) free(pmod->missmask);
	if (pmod->coeff != NULL) free(pmod->coeff);
	if (pmod->sderr != NULL) free(pmod->sderr);
	if (pmod->yhat != NULL) free(pmod->yhat);
	if (pmod->uhat != NULL) free(pmod->uhat);
	if (pmod->xpx != NULL) free(pmod->xpx);
	if (pmod->vcv != NULL) free(pmod->vcv);
	if (pmod->name != NULL) free(pmod->name);
	if (pmod->depvar != NULL) free(pmod->depvar);

	if (pmod->submask != NULL) {
	    free_subsample_mask(pmod->submask);
	}

	if (pmod->arinfo != NULL) {
	    clear_ar_info(pmod);
	}
	if (pmod->params != NULL) {
	    strings_array_free(pmod->params, pmod->nparams);
	}

	destroy_dataset(pmod->dataset);
	gretl_model_destroy_tests(pmod);
	destroy_all_data_items(pmod);
    }

    /* this may be redundant */
#if MDEBUG
    fprintf(stderr, "clear_model, calling gretl_model_init\n");
#endif
#if 1
    gretl_model_init(pmod, NULL);
#endif
}

void clear_model_xpx (MODEL *pmod)
{
    if (pmod->xpx != NULL) {
	free(pmod->xpx);
	pmod->xpx = NULL;
    }
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
	remove_model_from_stack_on_exit(pmod);
	clear_model(pmod);
	free(pmod);
    }
}

static void copy_test (ModelTest *targ, const ModelTest *src)
{
    targ->type = src->type;

    if (src->param != NULL && *src->param != '\0') {
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
    targ->opt = src->opt;
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
    int slop = (pmod->ci == GARCH);
    int np = 0;
    int err = 0;

    S = gretl_xml_get_strings_array(node, doc, &np, slop, &err);

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

    while (cur != NULL && !err) {
	gretl_test_init(&test, 0);
	got = 0;
	got += gretl_xml_get_prop_as_int(cur, "type", &test.type);
	got += gretl_xml_get_prop_as_uchar(cur, "teststat", &test.teststat);
	got += gretl_xml_get_prop_as_int(cur, "dfn", &test.dfn);
	got += gretl_xml_get_prop_as_double(cur, "dfd", &test.dfd);
	got += gretl_xml_get_prop_as_int(cur, "order", &test.order);
	got += gretl_xml_get_prop_as_double(cur, "value", &test.value);
	got += gretl_xml_get_prop_as_double(cur, "pvalue", &test.pvalue);
	got += gretl_xml_get_prop_as_double(cur, "crit", &test.crit);
	got += gretl_xml_get_prop_as_double(cur, "alpha", &test.alpha);
	if (got < 7) {
	    err = E_DATA;
	} else {
	    gretl_xml_get_prop_as_opt(cur, "opt", &test.opt);
	    gretl_xml_get_prop_as_string(cur, "param", &test.param);
	    err = real_add_test_to_model(pmod, &test);
	    free(test.param); /* copied by real_add_test_to_model */
	}
	cur = cur->next;
    }

    return err;
}

static void serialize_test (const ModelTest *src, PRN *prn)
{
    pprintf(prn, "<test type=\"%d\" ", src->type);

    if (src->param != NULL && *src->param != '\0') {
	pprintf(prn, "param=\"%s\" ", src->param);
    }

    pprintf(prn, "teststat=\"%d\" ", (int) src->teststat);
    pprintf(prn, "dfn=\"%d\" ", src->dfn);
    pprintf(prn, "dfd=\"%g\" ", src->dfd);
    pprintf(prn, "order=\"%d\" ", src->order);
    pprintf(prn, "value=\"%.15g\" ", src->value);
    pprintf(prn, "pvalue=\"%.15g\" ", src->pvalue);

    if (!na(src->crit)) {
	pprintf(prn, "crit=\"%g\" ", src->crit);
	pprintf(prn, "alpha=\"%g\" ", src->alpha);
    }

    if (src->opt != OPT_NONE) {
	pprintf(prn, "opt=\"%u\" ", (unsigned) src->opt);
    }

    pputs(prn, "/>\n");
}

static gretl_bundle *bundlize_test (const ModelTest *src)
{
    gretl_bundle *b = gretl_bundle_new();

    if (b == NULL) {
	return NULL;
    }

    if (src->param != NULL && *src->param != '\0') {
	gretl_bundle_set_string(b, "param", src->param);
    }

    if (src->dfn > 0) {
	gretl_bundle_set_scalar(b, "dfn", (double) src->dfn);
    }
    if (src->dfd > 0) {
	gretl_bundle_set_scalar(b, "dfd", src->dfd);
    }
    if (src->order > 0) {
	gretl_bundle_set_scalar(b, "order", src->order);
    }
    if (!na(src->value)) {
	gretl_bundle_set_scalar(b, "test", src->value);
    }
    if (!na(src->pvalue)) {
	gretl_bundle_set_scalar(b, "pvalue", src->pvalue);
    }
    if (!na(src->crit)) {
	gretl_bundle_set_scalar(b, "crit", src->crit);
	gretl_bundle_set_scalar(b, "alpha", src->alpha);
    }

    return b;
}

static int serialize_model_tests (const MODEL *pmod, PRN *prn)
{
    int i, n = pmod->ntests;

    if (n <= 0 || pmod->tests == NULL) {
	return 0;
    }

    pprintf(prn, "<tests count=\"%d\">\n", pmod->ntests);

    for (i=0; i<n; i++) {
	serialize_test(&pmod->tests[i], prn);
    }

    pputs(prn, "</tests>\n");

    return 0;
}

static int testvals_differ (double x, double y)
{
    double eq_tol = 1.0e-10;
    double reldiff;

    if (x == 0.0) {
	reldiff = fabs(y);
    } else if (y == 0.0) {
	reldiff = fabs(x);
    } else if (x > y) {
	reldiff = fabs((x - y) / y);
    } else {
	reldiff = fabs((y - x) / x);
    }

    return reldiff > eq_tol;
}

static int test_params_differ (const char *p, const char *s)
{
    int ret = 0;

    if (s == NULL && p != NULL) {
	ret = 1;
    } else if (p == NULL && s != NULL) {
	ret = 1;
    } else if (p != NULL && s != NULL) {
	ret = strcmp(p, s);
    }

    return ret;
}

static int model_tests_differ (ModelTest *mt1, ModelTest *mt2)
{
    int ret = 0;

    if (mt1->type != mt2->type) {
	ret = 1;
    } else if (mt1->order != mt2->order) {
	ret = 1;
    } else if (mt1->teststat != mt2->teststat) {
	ret = 1;
    } else if (test_params_differ(mt1->param, mt2->param)) {
	ret = 1;
    } else if (testvals_differ(mt1->value, mt2->value)) {
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
    int i, done = 0, add = 0;

    if (test == NULL || test->teststat == GRETL_STAT_NONE) {
	return 0;
    }

    for (i=0; i<pmod->ntests; i++) {
	if (!model_tests_differ(test, &pmod->tests[i])) {
	    done = 1;
	}
    }

    if (!done) {
	int n = pmod->ntests + 1;
	ModelTest *tests;

	tests = realloc(pmod->tests, n * sizeof *tests);

	if (tests != NULL) {
	    pmod->tests = tests;
	    copy_test(&pmod->tests[n-1], test);
	    pmod->ntests += 1;
	    add = 1;
	}
    }

    free(test->param);
    free(test);

    return add;
}

int gretl_model_get_normality_test (const MODEL *pmod, PRN *prn)
{
    ModelTest *test = NULL;
    int i, err = 0;

    for (i=0; i<pmod->ntests; i++) {
	if (pmod->tests[i].type == GRETL_TEST_NORMAL) {
	    test = &pmod->tests[i];
	    break;
	}
    }

    if (test == NULL) {
	err = E_BADSTAT;
    } else {
	record_test_result(test->value, test->pvalue);
	gretl_model_test_print(pmod, i, prn);
    }

    return err;
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

void model_test_set_dfd (ModelTest *test, double df)
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

void model_test_set_opt (ModelTest *test, gretlopt opt)
{
    test->opt = opt;
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
    { GRETL_TEST_CHOWDUM,
      N_("Chow test for structural difference with respect to %s"),
      N_("no structural difference") },
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
    { GRETL_TEST_BP,
      N_("Breusch-Pagan test for heteroskedasticity"),
      N_("heteroskedasticity not present") },
    { GRETL_TEST_SARGAN,
      N_("Sargan over-identification test"),
      N_("all instruments are valid") },
    { GRETL_TEST_IV_HAUSMAN,
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
      N_("Wald joint test on time dummies"),
      N_("No time effects") },
    { GRETL_TEST_HET_1,
      N_("Pesaran-Taylor test for heteroskedasticity"),
      N_("heteroskedasticity not present") },
    { GRETL_TEST_COMFAC,
      N_("Test of common factor restriction"),
      N_("restriction is acceptable") },
    { GRETL_TEST_INDEP,
      N_("Test of independence"),
      N_("rho = 0") },
    { GRETL_TEST_RE,
      N_("LR test for rho = 0"),
      NULL },
    { GRETL_TEST_WITHIN_F,
      N_("Joint test on named regressors"),
      NULL },
    { GRETL_TEST_RE_WALD,
      N_("Joint test on named regressors"),
      NULL },
    { GRETL_TEST_PANEL_WELCH,
      N_("Robust test for differing group intercepts"),
      N_("The groups have a common intercept") },
    { GRETL_TEST_XDEPEND,
      N_("Pesaran CD test for cross-sectional dependence"),
      N_("No cross-sectional dependence") },
    { GRETL_TEST_PANEL_AR,
      N_("Wooldridge test for autocorrelation in panel data"),
      N_("No first-order autocorrelation (rho = -0.5)") },
    { GRETL_TEST_MAX, NULL, NULL }
};

static const char *test_get_descrip (const ModelTest *test)
{
    static const char *dfree =
	N_("Distribution free Wald test for heteroskedasticity");
    static const char *lr_indep =
	N_("Likelihood ratio test of independence");

    if (test->type == GRETL_TEST_GROUPWISE &&
	test->teststat == GRETL_STAT_WALD_CHISQ) {
	return dfree;
    } else if (test->type == GRETL_TEST_INDEP &&
	       test->teststat == GRETL_STAT_LR) {
	return lr_indep;
    } else {
	int i;

	for (i=0; tstrings[i].ID < GRETL_TEST_MAX; i++) {
	    if (test->type == tstrings[i].ID) {
		return tstrings[i].descrip;
	    }
	}
    }

    return NULL;
}

static const char *test_get_opt_string (const ModelTest *test)
{
    if (test->type == GRETL_TEST_RESET) {
	if (test->opt & OPT_U) {
	    return N_("squares only");
	} else if (test->opt & OPT_C) {
	    return N_("cubes only");
	}
    } else if (test->type == GRETL_TEST_WHITES) {
	if (test->opt & OPT_X) {
	    return N_("squares only");
	}
    } else if (test->type == GRETL_TEST_BP) {
	if (test->opt & OPT_R) {
	    return N_("robust variant");
	}
    }

    return NULL;
}

static void print_test_opt (const ModelTest *test, PRN *prn)
{
    const char *optstr = test_get_opt_string(test);

    if (optstr != NULL) {
	pprintf(prn, " (%s)", _(optstr));
    }
}

static int gretl_test_print_heading (const ModelTest *test, PRN *prn)
{
    const char *descrip, *param = NULL;
    char ordstr[16];

    descrip = test_get_descrip(test);
    if (descrip == NULL) {
	return 1;
    }

    if (test->order > 0) {
	sprintf(ordstr, "%d", test->order);
	param = ordstr;
    } else if (test->type == GRETL_TEST_CHOW ||
	       test->type == GRETL_TEST_CHOWDUM) {
	param = test->param;
    }

    if (param != NULL) {
	if (plain_format(prn)) {
	    pprintf(prn, _(descrip), param);
	} else {
	    pprintf(prn, _(descrip), param);
	}
    } else {
	if (plain_format(prn)) {
	    pputs(prn, _(descrip));
	} else {
	    pputs(prn, _(descrip));
	}
    }

    if (test->opt) {
	print_test_opt(test, prn);
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

const char *get_h0_string_for_test (ModelTestType ttype)
{
    const char *ret = NULL;
    int i;

    for (i=0; tstrings[i].ID < GRETL_TEST_MAX; i++) {
	if (ttype == tstrings[i].ID) {
	    ret = tstrings[i].H0;
	}
    }

    return ret;
}

static void gretl_test_print_h_0 (const ModelTest *test, int heading,
				  PRN *prn)
{
    const char *H0 = NULL;
    int i;

    if (test->type == GRETL_TEST_PANEL_AR &&
	test->teststat == GRETL_STAT_STUDENT) {
	H0 = N_("No first-order autocorrelation (rho = 0)");
    } else {
	for (i=0; tstrings[i].ID < GRETL_TEST_MAX; i++) {
	    if (test->type == tstrings[i].ID) {
		H0 = tstrings[i].H0;
	    }
	}
    }

    if (test->type == GRETL_TEST_INDEP &&
	test->teststat == GRETL_STAT_WALD_CHISQ &&
	H0 != NULL && !strcmp(H0, "rho = 0")) {
	/* adjust for QML biprobit */
	H0 = N_("atanh(rho) = 0");
    }

    if (H0 == NULL && test->type != GRETL_TEST_WITHIN_F &&
	test->type != GRETL_TEST_RE_WALD) {
	return;
    }

    if (heading) {
	pputs(prn, " -");
	if (tex_format(prn)) {
	    pputc(prn, '-');
	}
    }

    if (H0 == NULL) {
	return;
    }

    if (plain_format(prn)) {
	pprintf(prn, "\n  %s: ", _("Null hypothesis"));
	pputs(prn, _(H0));
    } else if (tex_format(prn)) {
	pprintf(prn, "\\\\\n\\quad %s: ", _("Null hypothesis"));
	if (!strcmp(H0, "rho = 0")) {
	    pputs(prn, "$\\rho = 0$");
	} else {
	    pputs(prn, _(H0));
	}
    } else if (rtf_format(prn)) {
	pprintf(prn, "\\par\n %s: ", _("Null hypothesis"));
	pputs(prn, _(H0));
    }

    if (test->type == GRETL_TEST_ADD || test->type == GRETL_TEST_OMIT) {
	print_add_omit_varnames(test->param, prn);
    }
}

static gchar *get_test_stat_string (const ModelTest *test, PRN *prn)
{
    int tex = tex_format(prn);
    gchar *ret = NULL;

    switch (test->teststat) {
    case GRETL_STAT_LM:
	ret = g_strdup_printf("LM = %g", test->value);
	break;
    case GRETL_STAT_F:
    case GRETL_STAT_RESET:
	if (tex) {
	    ret = g_strdup_printf("$F(%d, %g)$ = %g", test->dfn, test->dfd, test->value);
	} else {
	    ret = g_strdup_printf("F(%d, %g) = %g", test->dfn, test->dfd, test->value);
	}
	break;
    case GRETL_STAT_SUP_WALD:
	if (tex) {
	    ret = g_strdup_printf("max $\\chi^2(%d)$ = %g (%s)", test->dfn,
				  test->value, test->param);
	} else {
	    ret = g_strdup_printf(_("chi-square(%d) = %g at observation %s"),
				  test->dfn, test->value, test->param);
	}
	break;
    case GRETL_STAT_LMF:
	ret = g_strdup_printf("LMF = %g", test->value);
	break;
    case GRETL_STAT_WF:
	if (tex) {
	    ret = g_strdup_printf("Welch $F(%d, %.1f)$ = %g", test->dfn, test->dfd, test->value);
	} else {
	    ret = g_strdup_printf("Welch F(%d, %.1f) = %g", test->dfn, test->dfd, test->value);
	}
	break;
    case GRETL_STAT_HARVEY_COLLIER:
	if (tex) {
	    ret = g_strdup_printf("Harvey--Collier $t(%d)$ = %g", test->dfn, test->value);
	} else {
	    ret = g_strdup_printf("Harvey-Collier t(%d) = %g", test->dfn, test->value);
	}
	break;
    case GRETL_STAT_NORMAL_CHISQ:
	if (tex) {
	    ret = g_strdup_printf("$\\chi^2(2)$ = %g", test->value);
	} else {
	    ret = g_strdup_printf("%s(2) = %g", _("Chi-square"), test->value);
	}
	break;
    case GRETL_STAT_LR:
    case GRETL_STAT_WALD_CHISQ:
    case GRETL_STAT_LB_CHISQ:
	if (tex) {
	    if (na(test->value)) {
		ret = g_strdup_printf("$\\chi^2(%d)$ = NA (%s)", test->dfn, _("failed"));
	    } else {
		ret = g_strdup_printf("$\\chi^2(%d)$ = %g", test->dfn, test->value);
	    }
	} else {
	    if (na(test->value)) {
		ret = g_strdup_printf("%s(%d) = NA (%s)", _("Chi-square"), test->dfn,
				      _("failed"));
	    } else {
		ret = g_strdup_printf("%s(%d) = %g", _("Chi-square"), test->dfn, test->value);
	    }
	}
	break;
    case GRETL_STAT_Z:
	if (tex) {
	    ret = g_strdup_printf("$z$ = %g", test->value);
	} else {
	    ret = g_strdup_printf("z = %g", test->value);
	}
	break;
    case GRETL_STAT_STUDENT:
	if (tex) {
	    ret = g_strdup_printf("$t$(%d) = %g", test->dfn, test->value);
	} else {
	    ret = g_strdup_printf("t(%d) = %g", test->dfn, test->value);
	}
	break;
    default:
	break;
    }

    return ret;
}

static gchar *get_test_pval_string (const ModelTest *test, PRN *prn)
{
    int tex = tex_format(prn);
    gchar *ret = NULL;

    switch (test->teststat) {
    case GRETL_STAT_LM:
	if (tex) {
	    ret = g_strdup_printf("$P$($\\chi^2(%d) >$ %g) = %g",
				  test->dfn, test->value, test->pvalue);
	} else {
	    ret = g_strdup_printf("P(%s(%d) > %g) = %g", _("Chi-square"),
				  test->dfn, test->value, test->pvalue);
	}
	break;
    case GRETL_STAT_F:
    case GRETL_STAT_RESET:
	if (tex) {
	    ret = g_strdup_printf("$P$($F(%d, %g) >$ %g) = %g",
				  test->dfn, test->dfd, test->value, test->pvalue);
	} else {
	    ret = g_strdup_printf("P(F(%d, %g) > %g) = %g",
				  test->dfn, test->dfd, test->value, test->pvalue);
	}
	break;
    case GRETL_STAT_WF:
	if (tex) {
	    ret = g_strdup_printf("$P$($F(%d, %.1f) >$ %g) = %g",
				  test->dfn, test->dfd, test->value, test->pvalue);
	} else {
	    ret = g_strdup_printf("P(F(%d, %.1f) > %g) = %g",
				  test->dfn, test->dfd, test->value, test->pvalue);
	}
	break;
    case GRETL_STAT_LMF:
	if (tex) {
	    ret = g_strdup_printf("$P$($F(%d, %g) >$ %g) = %g",
				  test->dfn, test->dfd, test->value, test->pvalue);
	} else {
	    ret = g_strdup_printf("P(F(%d, %g) > %g) = %g",
				  test->dfn, test->dfd, test->value, test->pvalue);
	}
	break;
    case GRETL_STAT_HARVEY_COLLIER:
	if (tex) {
	    ret = g_strdup_printf("$P$($t_{%d} >$ %g) = %g",
				  test->dfn, test->value, test->pvalue);
	} else {
	    ret = g_strdup_printf("P(t(%d) > %g) = %g",
				  test->dfn, test->value, test->pvalue);
	}
	break;
    case GRETL_STAT_STUDENT:
	if (tex) {
	    ret = g_strdup_printf("$P$($|t| >$ %g) = %g",
				  fabs(test->value), test->pvalue);
	} else {
	    ret = g_strdup_printf("P(|t| > %g) = %g",
				  fabs(test->value), test->pvalue);
	}
	break;
    case GRETL_STAT_NORMAL_CHISQ:
    case GRETL_STAT_LR:
    case GRETL_STAT_WALD_CHISQ:
    case GRETL_STAT_SUP_WALD:
    case GRETL_STAT_Z:
	if (na(test->value)) {
	    ; /* leave @ret NULL */
	} else if (na(test->pvalue)) {
	    ret = g_strdup("NA");
	} else {
	    ret = g_strdup_printf("%g", test->pvalue);
	}
	break;
    default:
	break;
    }

    return ret;
}

static void csv_print_test (const ModelTest *src, PRN *prn)
{
    const char *descrip = NULL;
    char c = prn_delim(prn);

    descrip = test_get_descrip(src);

    if (descrip != NULL) {
	const char *optstr = test_get_opt_string(src);

	if (optstr != NULL) {
	    pprintf(prn, "\"%s, %s\"\n", descrip, optstr);
	} else {
	    pprintf(prn, "\"%s\"\n", descrip);
	}
    }

    if (src->param != NULL && *src->param != '\0') {
	pprintf(prn, "\"%s\"%c\"%s\"\n", _("parameter"), c, src->param);
    }

    if (src->dfn > 0 && src->dfd > 0) {
	pprintf(prn, "\"%s\"%c%d\n", _("dfn"), c, src->dfn);
	pprintf(prn, "\"%s\"%c%g\n", _("dfd"), c, src->dfd);
    } else if (src->dfn > 0) {
	pprintf(prn, "\"%s\"%c%d\n", _("df"), c, src->dfn);
    }
    if (src->order) {
	pprintf(prn, "\"%s\"%c%d\n", _("lag order"), c, src->order);
    }
    pprintf(prn, "\"%s\"%c%.15g\n", _("test statistic"), c, src->value);
    if (!na(src->pvalue)) {
	pprintf(prn, "\"%s\"%c%.15g\n", _("p-value"), c, src->pvalue);
    }

    if (!na(src->crit)) {
	double a = 100 * src->alpha;
	gchar *buf;

	buf = g_strdup_printf(_("%g percent critical value"), a);
	pprintf(prn, "\"%s\"%c%g\n", buf, c, src->crit);
	g_free(buf);
    }

    pputc(prn, '\n');
}

#define asy_test(t) (t == GRETL_STAT_WALD_CHISQ || \
		     t == GRETL_STAT_Z)

#define asy_pval(test) (test->teststat == GRETL_STAT_SUP_WALD || \
			(test->opt & OPT_A))

#define boot_pval(test) (test->opt & OPT_B)

void gretl_model_test_print_direct (const ModelTest *test, int heading, PRN *prn)
{
    const char *tstr;
    gchar *buf = NULL;

    if (rtf_format(prn)) {
	pputs(prn, "\\par \\ql ");
    }

    if (heading) {
	gretl_test_print_heading(test, prn);
    }

    gretl_test_print_h_0(test, heading, prn);

    tstr = asy_test(test->teststat)? N_("Asymptotic test statistic") :
	N_("Test statistic");

    buf = get_test_stat_string(test, prn);

    if (plain_format(prn)) {
	pprintf(prn, "\n  %s: %s\n", _(tstr), buf);
    } else if (tex_format(prn)) {
	pprintf(prn, "\\\\\n\\quad %s: %s\\\\\n", _(tstr), buf);
    } else if (rtf_format(prn)) {
	pprintf(prn, "\\par\n %s: %s\\par\n", _(tstr), buf);
    }

    g_free(buf);
    buf = get_test_pval_string(test, prn);

    if (buf != NULL) {
	const char *pvstr =
	    asy_pval(test) ? N_("with asymptotic p-value") :
	    boot_pval(test) ? N_("with bootstrap p-value") :
	    N_("with p-value");

	if (plain_format(prn)) {
	    pprintf(prn, "  %s = %s\n\n", _(pvstr), buf);
	} else if (tex_format(prn)) {
	    pprintf(prn, "\\quad %s = %s\\\\\n", _(pvstr), buf);
	} else if (rtf_format(prn)) {
	    pprintf(prn, " %s = %s\\par\n\n", _(pvstr), buf);
	}
    } else if (!na(test->crit) && !na(test->alpha)) {
	double a = test->alpha * 100.0;

	if (plain_format(prn)) {
	    buf = g_strdup_printf(_("%g percent critical value"), a);
	    pprintf(prn, "  (%s = %.2f)\n\n", buf, test->crit);
	} else if (tex_format(prn)) {
	    buf = g_strdup_printf(_("%g percent critical value"), a);
	    pprintf(prn, "\\quad (%s = %.2f)\\\\\n", buf, test->crit);
	} else if (rtf_format(prn)) {
	    buf = g_strdup_printf(_("%g percent critical value"), a);
	    pprintf(prn, " (%s = %.2f)\\par\n\n", buf, test->crit);
	}
    } else {
	pputc(prn, '\n');
    }

    g_free(buf);
}

void gretl_model_test_print (const MODEL *pmod, int i, PRN *prn)
{
    if (i >= 0 && i < pmod->ntests) {
	if (csv_format(prn)) {
	    csv_print_test(&pmod->tests[i], prn);
	} else {
	    gretl_model_test_print_direct(&pmod->tests[i], 1, prn);
	}
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
    int err = 0;

    if (src->nparams > 0) {
	targ->params = strings_array_dup(src->params, src->nparams);

	if (targ->params == NULL) {
	    err = E_ALLOC;
	} else {
	    targ->nparams = src->nparams;
	}
    }

    return err;
}

static void serialize_coeff_sep (model_data_item *item, PRN *prn)
{
    CoeffSep *cs = (CoeffSep *) item->ptr;

    pprintf(prn, " pos=\"%d\"", cs->pos);
    if (*cs->str != '\0') {
	pputs(prn, " string=\"");
	gretl_xml_put_string(cs->str, prn);
	pputc(prn, '"');
    }
    pputs(prn, "/>\n");
}

static void serialize_vcv_info (model_data_item *item, PRN *prn)
{
    VCVInfo *vi = (VCVInfo *) item->ptr;

    pprintf(prn, " vmaj=\"%d\"", vi->vmaj);
    pprintf(prn, " vmin=\"%d\"", vi->vmin);
    if (vi->order > 0) {
	pprintf(prn, " order=\"%d\"", vi->order);
    }
    if (vi->flags > 0) {
	pprintf(prn, " flags=\"%d\"", vi->flags);
    }
    if (vi->bw > 0 && !na(vi->bw)) {
	pprintf(prn, " bw=\"%.12g\"", vi->bw);
    }
    if (vi->cv1 != NULL) {
	pprintf(prn, " cv1=\"%s\"", vi->cv1);
    }
    if (vi->cv2 != NULL) {
	pprintf(prn, " cv2=\"%s\"", vi->cv2);
    }
    pputs(prn, "/>\n");
}

/* FIXME updating and placement of these functions */

struct type_mapper {
    int type;
    const char *name;
    const char *compat;
};

static struct type_mapper mapper[] = {
    { GRETL_TYPE_INT,          "int",         "1" },
    { GRETL_TYPE_LIST,         "list",        "2" },
    { GRETL_TYPE_DOUBLE,       "double",      "3" },
    { GRETL_TYPE_INT_ARRAY,    "intarray",    "4" },
    { GRETL_TYPE_DOUBLE_ARRAY, "doublearray", "5" },
    { GRETL_TYPE_STRING,       "string",      "6" },
    { GRETL_TYPE_CMPLX_ARRAY,  "cmplxarray",  "8" },
    { GRETL_TYPE_STRUCT,       "struct" ,     "9" },
    { GRETL_TYPE_MATRIX,       "matrix" ,    "10" },
    { GRETL_TYPE_NONE,         "none",        "0" },
};

static int type_from_type_string (const char *s)
{
    int i;

    if (isdigit(*s)) {
	/* backward compatibility */
	for (i=0; mapper[i].type != GRETL_TYPE_NONE; i++) {
	    if (*s == mapper[i].compat[0]) {
		return mapper[i].type;
	    }
	}
	if (*s == '7') {
	    /* note: was "chararray" */
	    return GRETL_TYPE_STRING;
	}
    } else {
	for (i=0; mapper[i].type != GRETL_TYPE_NONE; i++) {
	    if (!strcmp(s, mapper[i].name)) {
		return mapper[i].type;
	    }
	}
	if (!strcmp(s, "chararray")) {
	    return GRETL_TYPE_STRING;
	}
    }

    return GRETL_TYPE_NONE;
}

static const char *gretl_type_name (int t)
{
    int i;

    for (i=0; mapper[i].type != GRETL_TYPE_NONE; i++) {
	if (t == mapper[i].type) {
	    return mapper[i].name;
	}
    }

    return "unknown";
}

static void serialize_model_data_items (const MODEL *pmod, PRN *prn)
{
    model_data_item *item;
    int i, j, nelem;

    pprintf(prn, "<data-items count=\"%d\">\n", pmod->n_data_items);

    for (i=0; i<pmod->n_data_items; i++) {
	item = pmod->data_items[i];
	nelem = 0;

	pprintf(prn, "<data-item key=\"%s\" type=\"%s\"",
		item->key, gretl_type_name(item->type));

	if (!strcmp(item->key, "coeffsep")) {
	    serialize_coeff_sep(item, prn);
	    continue;
	} else if (!strcmp(item->key, "vcv_info")) {
	    serialize_vcv_info(item, prn);
	    continue;
	}

	if (item->type == GRETL_TYPE_INT_ARRAY) {
	    nelem = item->size / sizeof(int);
	} else if (item->type == GRETL_TYPE_DOUBLE_ARRAY) {
	    nelem = item->size / sizeof(double);
	} else if (item->type == GRETL_TYPE_CMPLX_ARRAY) {
	    nelem = item->size / sizeof(cmplx);
	}

	if (nelem > 0) {
	    pprintf(prn, " count=\"%d\">\n", nelem);
	} else if (item->type == GRETL_TYPE_STRING) {
	    pputc(prn, '>');
	} else {
	    pputs(prn, ">\n");
	}

	if (item->type == GRETL_TYPE_INT) {
	    pprintf(prn, "%d", *(int *) item->ptr);
	} else if (item->type == GRETL_TYPE_DOUBLE) {
	    double x = *(double *) item->ptr;

	    if (na(x)) {
		pputs(prn, "NA");
	    } else {
		pprintf(prn, "%.15g", x);
	    }
	} else if (item->type == GRETL_TYPE_INT_ARRAY) {
	    int *vals = (int *) item->ptr;

	    for (j=0; j<nelem; j++) {
		pprintf(prn, "%d ", vals[j]);
	    }
	} else if (item->type == GRETL_TYPE_DOUBLE_ARRAY) {
	    double *vals = (double *) item->ptr;

	    for (j=0; j<nelem; j++) {
		if (na(vals[j])) {
		    pputs(prn, "NA ");
		} else {
		    pprintf(prn, "%.15g ", vals[j]);
		}
	    }
	} else if (item->type == GRETL_TYPE_CMPLX_ARRAY) {
	    cmplx *vals = (cmplx *) item->ptr;
	    double x;

	    for (j=0; j<nelem; j++) {
		x = vals[j].r;
		if (na(x)) {
		    pputs(prn, "NA ");
		} else {
		    pprintf(prn, "%.15g ", x);
		}
		x = vals[j].i;
		if (na(x)) {
		    pputs(prn, "NA ");
		} else {
		    pprintf(prn, "%.15g ", x);
		}
	    }
	} else if (item->type == GRETL_TYPE_LIST) {
	    int *list = (int *) item->ptr;

	    for (j=0; j<=list[0]; j++) {
		pprintf(prn, "%d ", list[j]);
	    }
	} else if (item->type == GRETL_TYPE_STRING) {
	    gretl_xml_put_string((const char *) item->ptr, prn);
	} else if (item->type == GRETL_TYPE_MATRIX) {
	    gretl_matrix *m = (gretl_matrix *) item->ptr;

	    gretl_matrix_serialize(m, NULL, prn);
	} else if (item->type == GRETL_TYPE_ARRAY) {
	    gretl_array *A = (gretl_array *) item->ptr;

	    gretl_array_serialize(A, prn);
	} else {
	    ; /* no-op: not handled */
	}

	pputs(prn, "</data-item>\n");
    }

    pputs(prn, "</data-items>\n");
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

/* ensure we get a valid identifier for the bundle key */

static char *item_key_to_bundle_key (char *targ,
				     const char *src)
{
    strcpy(targ, src);
    return gretl_charsub(targ, '-', '_');
}

static gretl_matrix *my_matrix_from_list (model_data_item *item)
{
    int *list = item->ptr;
    gretl_matrix *m;
    int i;

    m = gretl_matrix_alloc(1, list[0]);

    if (m != NULL) {
	for (i=0; i<list[0]; i++) {
	    m->val[i] = list[i+1];
	}
    }

    return m;
}

static gretl_matrix *matrix_from_cmplx (model_data_item *item)
{
    int nr = item->size / sizeof(cmplx);
    cmplx *c = item->ptr;
    gretl_matrix *m;
    int i;

    m = gretl_matrix_alloc(nr, 2);

    if (m != NULL) {
	for (i=0; i<nr; i++) {
	    gretl_matrix_set(m, i, 0, c[i].r);
	    gretl_matrix_set(m, i, 1, c[i].i);
	}
    }

    return m;
}

int bundlize_model_data_items (const MODEL *pmod, gretl_bundle *b)
{
    model_data_item *item;
    char bkey[VNAMELEN];
    gretl_matrix *m;
    double xval;
    int i, ival;
    int err = 0;

    for (i=0; i<pmod->n_data_items && !err; i++) {
	item = pmod->data_items[i];
	item_key_to_bundle_key(bkey, item->key);
	if (gretl_bundle_has_key(b, bkey)) {
	    continue; /* item already present */
	} else if (item->type == GRETL_TYPE_INT) {
	    ival = *(int *) item->ptr;
	    err = gretl_bundle_set_scalar(b, bkey, ival);
	} else if (item->type == GRETL_TYPE_DOUBLE) {
	    xval = *(double *) item->ptr;
	    err = gretl_bundle_set_scalar(b, bkey, xval);
	} else if (item->type == GRETL_TYPE_MATRIX) {
	    err = gretl_bundle_set_matrix(b, bkey, item->ptr);
	} else if (item->type == GRETL_TYPE_LIST) {
	    if (pmod->ci == MIDASREG && !strcmp(bkey, "seplist")) {
		; /* internal use only, skip it */
	    } else if (pmod->ci == ARMA && !strcmp(bkey, "ainfo")) {
		m = my_matrix_from_list(item);
		if (m != NULL) {
		    err = gretl_bundle_donate_data(b, bkey, m,
						   GRETL_TYPE_MATRIX,
						   0);
		}
	    } else {
		err = gretl_bundle_set_list(b, bkey, item->ptr);
	    }
	} else if (item->type == GRETL_TYPE_ARRAY) {
	    gretl_bundle_set_data(b, bkey, item->ptr, item->type, 0);
	} else if (item->type == GRETL_TYPE_CMPLX_ARRAY) {
	    m = matrix_from_cmplx(item);
	    if (m != NULL) {
		err = gretl_bundle_donate_data(b, bkey, m,
					       GRETL_TYPE_MATRIX,
					       0);
	    }
	}
    }

    if (!err && pmod->tests != NULL) {
	const ModelTest *test;
	const char *tkey;
	gretl_bundle *bt;

	for (i=0; i<pmod->ntests && !err; i++) {
	    test = &pmod->tests[i];
	    tkey = test_type_key(test->type);
	    if (tkey != NULL && !gretl_bundle_has_key(b, tkey)) {
		bt = bundlize_test(test);
		if (bt != NULL) {
		    err = gretl_bundle_donate_data(b, tkey, bt,
						   GRETL_TYPE_BUNDLE, 0);
		}
	    }
	}
    }

    if (!err) {
	/* estimation start and end, 1-based */
	gretl_bundle_set_int(b, "t1", pmod->t1 + 1);
	gretl_bundle_set_int(b, "t2", pmod->t2 + 1);
	/* contemporaneous sample range, 1-based */
	gretl_bundle_set_int(b, "sample_t1", pmod->smpl.t1 + 1);
	gretl_bundle_set_int(b, "sample_t2", pmod->smpl.t2 + 1);
    }

    if (pmod->ci == PANEL && (pmod->opt & OPT_F) &&
	!gretl_bundle_has_key(b, "within_rsq")) {
	/* fixed effects: add within R-squared */
	gretl_bundle_set_scalar(b, "within_rsq", pmod->adjrsq);
    }

    if (!na(pmod->dw) && !gretl_bundle_has_key(b, "dw")) {
	/* add Durbin-Watson statistic if available */
	gretl_bundle_set_scalar(b, "dw", pmod->dw);
    }

    gretl_bundle_set_int(b, "asymptotic", ASYMPTOTIC_MODEL(pmod->ci));

    return err;
}

void display_model_data_items (const MODEL *pmod)
{
    model_data_item *item;
    int i, n = pmod->n_data_items;

    fprintf(stderr, "model has %d data items\n", n);

    for (i=0; i<n; i++) {
	item = pmod->data_items[i];
	fprintf(stderr, "%d '%s': ", i, item->key);
	if (item->type == GRETL_TYPE_INT) {
	    fprintf(stderr, "%d\n", *(int *) item->ptr);
	} else if (item->type == GRETL_TYPE_DOUBLE) {
	    fprintf(stderr, "%.15g\n", *(double *) item->ptr);
	} else {
	    fprintf(stderr, "%p\n", (void *) item->ptr);
	}
    }
}

static int copy_model (MODEL *targ, MODEL *src)
{
    int k = src->ncoeff;
    int m = k * (k + 1) / 2;
    int err = 0;

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
	(targ->submask = copy_subsample_mask(src->submask, &err)) == NULL) {
	return 1;
    }

    if (src->missmask != NULL &&
	(targ->missmask = gretl_strdup(src->missmask)) == NULL) {
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
	err = copy_model_params(targ, src);
	if (err) {
	    return err;
	}
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

    if (gretl_is_between_model(src) && src->dataset != NULL) {
	/* special for "between" panel model: transfer the
	   group-means dataset to @targ
	*/
	targ->dataset = src->dataset;
	src->dataset = NULL;
    }

    return 0;
}

int gretl_model_serialize (const MODEL *pmod, SavedObjectFlags flags,
			   PRN *prn)
{
    char *xmlname = NULL;
    int k = pmod->ncoeff;
    int m = k * (k + 1) / 2;
    int err = 0;

    if (pmod->name == NULL || *pmod->name == '\0') {
	xmlname = gretl_strdup("none");
    } else {
	xmlname = gretl_xml_encode(pmod->name);
    }

    pprintf(prn, "<gretl-model ID=\"%d\" name=\"%s\" saveflags=\"%d\" ",
	    pmod->ID, xmlname, (int) flags);
    free(xmlname);

    if (pmod->depvar != NULL) {
	/* In an nls model pmod->depvar may contain '<' or '>' */
	if (gretl_xml_validate(pmod->depvar)) {
	    pprintf(prn, "depvar=\"%s\" ", pmod->depvar);
	} else {
	    char *tmp = gretl_xml_encode(pmod->depvar);

	    pprintf(prn, "depvar=\"%s\" ", tmp);
	    free(tmp);
	}
    }

    pprintf(prn, "t1=\"%d\" t2=\"%d\" nobs=\"%d\" ",
	    pmod->t1, pmod->t2, pmod->nobs);
    pprintf(prn, "full_n=\"%d\" ncoeff=\"%d\" dfn=\"%d\" dfd=\"%d\" ",
	    pmod->full_n, pmod->ncoeff, pmod->dfn, pmod->dfd);
    pprintf(prn, "ifc=\"%d\" ci=\"%s\" opt=\"%u\" nwt=\"%d\" aux=\"%d\" ",
	    pmod->ifc, gretl_command_word(pmod->ci), pmod->opt,
	    pmod->nwt, pmod->aux);

    gretl_push_c_numeric_locale();

    gretl_xml_put_double("ess", pmod->ess, prn);
    gretl_xml_put_double("tss", pmod->tss, prn);
    gretl_xml_put_double("sigma", pmod->sigma, prn);
    gretl_xml_put_double("rsq", pmod->rsq, prn);
    gretl_xml_put_double("adjrsq", pmod->adjrsq, prn);
    gretl_xml_put_double("fstt", pmod->fstt, prn);
    gretl_xml_put_double("chi-square", pmod->chisq, prn);
    gretl_xml_put_double("lnL", pmod->lnL, prn);
    gretl_xml_put_double("ybar", pmod->ybar, prn);
    gretl_xml_put_double("sdy", pmod->sdy, prn);

    gretl_xml_put_double("crit0", pmod->criterion[0], prn);
    gretl_xml_put_double("crit1", pmod->criterion[1], prn);
    gretl_xml_put_double("crit2", pmod->criterion[2], prn);

    gretl_xml_put_double("dw", pmod->dw, prn);
    gretl_xml_put_double("rho", pmod->rho, prn);

    pputs(prn, ">\n");

    if (pmod->smpl.rseed > 0) {
	pprintf(prn, "<sample t1=\"%d\" t2=\"%d\" rseed=\"%u\"/>\n",
		pmod->smpl.t1, pmod->smpl.t2, pmod->smpl.rseed);
    } else {
	pprintf(prn, "<sample t1=\"%d\" t2=\"%d\"/>\n",
		pmod->smpl.t1, pmod->smpl.t2);
    }

    gretl_xml_put_double_array("coeff", pmod->coeff, k, prn);
    gretl_xml_put_double_array("sderr", pmod->sderr, k, prn);

    if (pmod->uhat != NULL) {
	gretl_xml_put_double_array("uhat", pmod->uhat, pmod->full_n, prn);
    }

    if (pmod->yhat != NULL) {
	gretl_xml_put_double_array("yhat", pmod->yhat, pmod->full_n, prn);
    }

    if (pmod->submask != NULL) {
	write_model_submask(pmod, prn);
    }

    if (pmod->missmask != NULL) {
	pputs(prn, "<missmask>");
	pputs(prn, pmod->missmask);
	pputs(prn, "</missmask>\n");
    }

    if (pmod->xpx != NULL) {
	gretl_xml_put_double_array("xpx", pmod->xpx, m, prn);
    }

    if (pmod->vcv != NULL) {
	gretl_xml_put_double_array("vcv", pmod->vcv, m, prn);
    }

    if (pmod->arinfo != NULL && pmod->arinfo->arlist != NULL) {
	int r = pmod->arinfo->arlist[0];

	pputs(prn, "<arinfo>\n");
	gretl_xml_put_tagged_list("arlist", pmod->arinfo->arlist, prn);
	gretl_xml_put_double_array("rho", pmod->arinfo->rho, r, prn);
	gretl_xml_put_double_array("sderr", pmod->arinfo->sderr, r, prn);
	pputs(prn, "</arinfo>\n");
    }

    if (pmod->ntests > 0 && pmod->tests != NULL) {
	serialize_model_tests(pmod, prn);
    }

    if (pmod->nparams > 0 && pmod->params != NULL) {
	gretl_xml_put_strings_array("params",
				    (const char **) pmod->params,
				    pmod->nparams, prn);
    }

    if (pmod->list != NULL) {
	gretl_xml_put_tagged_list("list", pmod->list, prn);
    }

    if (pmod->n_data_items > 0) {
	serialize_model_data_items(pmod, prn);
    }

    pputs(prn, "</gretl-model>\n");

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

    err = gretl_xml_get_submask(node, doc, &mask);
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
				   GRETL_TYPE_STRUCT,
				   sizeof *cs);
    }

    return err;
}

static int
retrieve_model_vcv_info (xmlNodePtr cur, MODEL *pmod)
{
    VCVInfo *vi = vcv_info_new();
    int err = 0;

    if (vi != NULL) {
	char *s = NULL;
	int ival;
	double x;

	gretl_xml_get_prop_as_int(cur, "vmaj", &vi->vmaj);
	gretl_xml_get_prop_as_int(cur, "vmin", &vi->vmin);
	if (gretl_xml_get_prop_as_int(cur, "order", &ival)) {
	    vi->order = ival;
	}
	if (gretl_xml_get_prop_as_int(cur, "flags", &ival)) {
	    vi->flags = ival;
	}
	if (gretl_xml_get_prop_as_double(cur, "bw", &x)) {
	    vi->bw = x;
	}
	if (gretl_xml_get_prop_as_string(cur, "cv1", &s)) {
	    vi->cv1 = s;
	}
	if (gretl_xml_get_prop_as_string(cur, "cv2", &s)) {
	    vi->cv2 = s;
	}
	err = gretl_model_set_data_with_destructor(pmod, "vcv_info", vi,
						   GRETL_TYPE_STRUCT,
						   sizeof *vi,
						   vcv_info_free);
    }

    return err;
}

/* backward compatibility: transcribe old-style "gretl_set_int"
   settings to bits in pmod->opt, where applicable
*/

static int gretl_model_set_int_compat (MODEL *pmod,
				       const char *key,
				       int ival)
{
    gretlopt opt = pmod->opt;
    int err = 0;

    if (ival == 0) {
	return 0;
    }

    if (pmod->ci == PANEL) {
	if (!strcmp(key, "between")) {
	    pmod->opt |= OPT_B;
	} else if (!strcmp(key, "fixed-effects")) {
	    pmod->opt |= OPT_F;
	} else if (!strcmp(key, "random-effects")) {
	    pmod->opt |= OPT_U;
	} else if (!strcmp(key, "unit-weights")) {
	    pmod->opt |= OPT_H;
	}
    } else if (pmod->ci == AR1) {
	if (!strcmp(key, "hilu")) {
	    pmod->opt |= OPT_H;
	} else if (!strcmp(key, "pwe")) {
	    pmod->opt |= OPT_P;
	}
    } else if (pmod->ci == GMM) {
	if (!strcmp(key, "iterated")) {
	    pmod->opt |= OPT_I;
	} else if (!strcmp(key, "two-step")) {
	    pmod->opt |= OPT_T;
	}
    } else if (pmod->ci == HECKIT) {
	if (!strcmp(key, "two-step")) {
	    pmod->opt |= OPT_T;
	}
    } else if (pmod->ci == LOGIT || pmod->ci == PROBIT) {
	if (!strcmp(key, "show-pvals")) {
	    pmod->opt |= OPT_P;
	}
    }

    if (!strcmp(key, "robust")) {
	pmod->opt |= OPT_R;
    } else if (!strcmp(key, "time-dummies")) {
	pmod->opt |= OPT_D;
    } else if (!strcmp(key, "no-df-corr")) {
	pmod->opt |= OPT_N;
    }

    if (pmod->opt == opt) {
	/* not handled above */
	err = gretl_model_set_int(pmod, key, ival);
    }

    return err;
}

/* backward compat: convert ad hoc settings to VCVInfo */

static void compat_compose_vcv_info (MODEL *pmod)
{
    VCVInfo vi = { 0, 0, 0, 0, NADBL, NULL, NULL };
    int ival;

    if (gretl_model_get_int(pmod, "hc")) {
	vi.vmaj = VCV_HC;
	vi.vmin = gretl_model_get_int(pmod, "hc_version");
    } else if ((ival = gretl_model_get_int(pmod, "ml_vcv"))) {
	vi.vmaj = VCV_ML;
	vi.vmin = ival;
    } else if (gretl_model_get_int(pmod, "panel_hac")) {
	vi.vmaj = VCV_PANEL;
	vi.vmin = PANEL_HAC;
    } else if (gretl_model_get_int(pmod, "panel_bk")) {
	vi.vmaj = VCV_PANEL;
	vi.vmin = PANEL_BK;
    } else if (gretl_model_get_int(pmod, "using_hac") ||
	       gretl_model_get_int(pmod, "hac_kernel") ||
	       gretl_model_get_int(pmod, "hac_lag")) {
	vi.vmaj = VCV_HAC;
	vi.vmin = gretl_model_get_int(pmod, "hac_kernel");
	vi.order = gretl_model_get_int(pmod, "hac_lag");
	vi.flags = gretl_model_get_int(pmod, "hac_prewhiten");
	if (vi.vmin == KERNEL_QS) {
	    vi.bw = gretl_model_get_double(pmod, "qs_bandwidth");
	}
    } else if (pmod->ci == LAD && gretl_model_get_int(pmod, "rq")) {
	double a = gretl_model_get_double(pmod, "rq_alpha");

	if (na(a)) {
	    /* doing VCV, not confidence intervals */
	    vi.vmaj = VCV_RQ;
	    vi.vmin = (gretl_model_get_int(pmod, "rq_nid"))?
		RQ_NID : RQ_ASY;
	}
    }

    if (vi.vmaj > 0) {
	VCVInfo *pvi = vcv_info_new();

	if (pvi != NULL) {
	    *pvi = vi;
	    gretl_model_set_data_with_destructor(pmod, "vcv_info", pvi,
						 GRETL_TYPE_STRUCT,
						 sizeof *pvi,
						 vcv_info_free);
	}
    }
}

/* backward compat: convert old list separator */

#define old_listsep(v) (v == 999 || v == 9999)

static void maybe_convert_listsep (int *list, const DATASET *dset)
{
    int i;

    for (i=1; i<=list[0]; i++) {
	if (old_listsep(list[i]) && list[i] >= dset->v) {
	    list[i] = LISTSEP;
	}
    }
}

static gretl_matrix *data_list_to_matrix (const int *list)
{
    gretl_matrix *y = gretl_column_vector_alloc(list[0]);
    int i;

    if (y != NULL) {
	for (i=0; i<y->rows; i++) {
	    y->val[i] = list[i+1];
	}
    }

    return y;
}

#define XDEBUG 0

static int model_data_items_from_xml (xmlNodePtr node, xmlDocPtr doc,
				      MODEL *pmod, int fixopt)
{
    xmlNodePtr cur;
    int n_items;
    int got_vcv = 0;
    int err = 0;

    if (!gretl_xml_get_prop_as_int(node, "count", &n_items)) {
	return 1;
    }

    cur = node->xmlChildrenNode;

    while (cur != NULL && !err) {
	char *key = NULL;
	char *typestr = NULL;
	int ignore_err = 0;
	int t, nelem = 0;

	if (!gretl_xml_get_prop_as_string(cur, "type", &typestr) ||
	    !gretl_xml_get_prop_as_string(cur, "key", &key)) {
	    err = 1;
	    break;
	}

#if XDEBUG
	fprintf(stderr, "data_item: type='%s', key='%s'\n", typestr, key);
#endif

	t = type_from_type_string(typestr);

	if (!strcmp(key, "roots")) {
	    /* could be a backward-compat problem */
	    ignore_err = 1;
	}

	if (!strcmp(key, "coeffsep") || !strcmp(key, "10")) {
	    /* special, with backward compatibility */
	    err = retrieve_model_coeff_separator(cur, pmod);
	} else if (!strcmp(key, "vcv_info")) {
	    err = retrieve_model_vcv_info(cur, pmod);
	    if (!err) {
		got_vcv = 1;
	    }
	} else if (t == GRETL_TYPE_INT) {
	    int ival;

	    if (!gretl_xml_node_get_int(cur, doc, &ival)) {
		err = 1;
	    } else if (fixopt) {
		err = gretl_model_set_int_compat(pmod, key, ival);
	    } else {
		err = gretl_model_set_int(pmod, key, ival);
	    }
	} else if (t == GRETL_TYPE_DOUBLE) {
	    double xval;

	    if (!gretl_xml_node_get_double(cur, doc, &xval)) {
		err = 1;
	    } else {
		err = gretl_model_set_double(pmod, key, xval);
	    }
	} else if (t == GRETL_TYPE_LIST) {
	    if (!strcmp(key, "xlist")) {
		/* ad hoc (for forecasting): will be recreated if need be */
		;
	    } else if (!strcmp(key, "yvals")) {
		/* backward compatibility for mnlogit models */
		int *list = gretl_xml_get_list(cur, doc, &err);

		if (!err && list != NULL) {
		    gretl_matrix *yv = data_list_to_matrix(list);

		    if (yv != NULL) {
			err = gretl_model_set_matrix_as_data(pmod, key, yv);
		    }
		}
		free(list);
	    } else {
		int *list = gretl_xml_get_list(cur, doc, &err);

		if (!err && list != NULL) {
		    err = gretl_model_set_list_as_data(pmod, key, list);
		}
	    }
	} else if (t == GRETL_TYPE_STRING) {
	    char *s;
	    int gots;

	    if (!strcmp(key, "pmask") || !strcmp(key, "qmask")) {
		/* these masks must not have leading white space */
		gots = gretl_xml_node_get_trimmed_string(cur, doc, &s);
	    } else {
		gots = gretl_xml_node_get_string(cur, doc, &s);
	    }
	    if (!gots) {
		err = 1;
	    } else {
		err = gretl_model_set_string_as_data(pmod, key, s);
	    }
	} else if (t == GRETL_TYPE_INT_ARRAY) {
	    int *ivals = gretl_xml_get_int_array(cur, doc, &nelem, &err);

	    if (nelem > 0) {
		err = gretl_model_set_data(pmod, key, ivals, t,
					   nelem * sizeof *ivals);
	    }
	} else if (t == GRETL_TYPE_DOUBLE_ARRAY) {
	    double *xvals = gretl_xml_get_double_array(cur, doc, &nelem, &err);

	    if (nelem > 0) {
		err = gretl_model_set_data(pmod, key, xvals, t,
					   nelem * sizeof *xvals);
	    }
	} else if (t == GRETL_TYPE_CMPLX_ARRAY) {
	    cmplx *cvals = gretl_xml_get_cmplx_array(cur, doc, &nelem, &err);

	    if (nelem > 0) {
		err = gretl_model_set_data(pmod, key, cvals, t,
					   nelem * sizeof *cvals);
	    }
	} else if (t == GRETL_TYPE_MATRIX) {
	    xmlNodePtr child = cur->xmlChildrenNode;
	    gretl_matrix *m;

	    if (child == NULL) {
		err = E_DATA;
	    } else {
		m = gretl_xml_get_matrix(child, doc, &err);
		if (!err) {
		    err = gretl_model_set_matrix_as_data(pmod, key, m);
		}
	    }
	}

	if (err && ignore_err) {
	    fprintf(stderr, "ignoring error reloading %s\n", key);
	    err = 0;
	}

	xmlFree(key);
	xmlFree(typestr);

	cur = cur->next;
    }

    if (!err && !got_vcv) {
	compat_compose_vcv_info(pmod);
    }

#if XDEBUG
    fprintf(stderr, "data items from XML, returning %d\n", err);
#endif

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
	    pmod->arinfo->arlist = gretl_xml_get_list(cur, doc, &err);
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
 * @dset: dataset information.
 * @err: location to receive error code.
 *
 * Reads info on a gretl model from the given XML @node
 * and @doc, and reconstitutes the model in memory.
 *
 * Returns: allocated model, or %NULL on failure.
 */

MODEL *gretl_model_from_XML (xmlNodePtr node, xmlDocPtr doc,
			     const DATASET *dset,
			     int *err)
{
    MODEL *pmod;
    char *buf = NULL;
    char *modname = NULL;
    xmlNodePtr cur;
    int fixopt = 0;
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

    got += gretl_xml_get_prop_as_string(node, "ci", &buf);

    if (got < 10) {
	fprintf(stderr, "gretl_model_from_XML: got(1) = %d (expected 12)\n", got);
	if (buf != NULL) {
	    free(buf);
	}
	*err = E_DATA;
	goto bailout;
    }

    gretl_xml_get_prop_as_int(node, "nwt", &pmod->nwt);
    gretl_xml_get_prop_as_int(node, "aux", &pmod->aux);

    if (gretl_xml_get_prop_as_int(node, "opt", &n)) {
	pmod->opt = (unsigned) n;
    } else {
	/* backward compat may be needed */
	fixopt = 1;
    }

    if (isdigit(*buf)) {
	/* backward compatibility: the model command is given
	   as a number, not a string */
	pmod->ci = atoi(buf);
    } else {
	pmod->ci = gretl_command_number(buf);
    }

    free(buf);

    gretl_push_c_numeric_locale();

    gretl_xml_get_prop_as_double(node, "ess", &pmod->ess);
    gretl_xml_get_prop_as_double(node, "tss", &pmod->tss);
    gretl_xml_get_prop_as_double(node, "sigma", &pmod->sigma);
    gretl_xml_get_prop_as_double(node, "rsq", &pmod->rsq);
    gretl_xml_get_prop_as_double(node, "adjrsq", &pmod->adjrsq);
    gretl_xml_get_prop_as_double(node, "fstt", &pmod->fstt);
    gretl_xml_get_prop_as_double(node, "chi-square", &pmod->chisq);
    gretl_xml_get_prop_as_double(node, "lnL", &pmod->lnL);
    gretl_xml_get_prop_as_double(node, "ybar", &pmod->ybar);
    gretl_xml_get_prop_as_double(node, "sdy", &pmod->sdy);

    gretl_xml_get_prop_as_double(node, "crit0", &pmod->criterion[0]);
    gretl_xml_get_prop_as_double(node, "crit1", &pmod->criterion[1]);
    gretl_xml_get_prop_as_double(node, "crit2", &pmod->criterion[2]);

    gretl_xml_get_prop_as_double(node, "dw", &pmod->dw);
    gretl_xml_get_prop_as_double(node, "rho", &pmod->rho);

    gretl_xml_get_prop_as_string(node, "name", &modname);
    gretl_xml_get_prop_as_string(node, "depvar", &pmod->depvar);

    cur = node->xmlChildrenNode;

#if XDEBUG
    fprintf(stderr, "reading model child nodes...\n");
#endif

    while (cur != NULL) {
	if (!xmlStrcmp(cur->name, (XUC) "sample")) {
	    gretl_xml_get_prop_as_int(cur, "t1", &pmod->smpl.t1);
	    gretl_xml_get_prop_as_int(cur, "t2", &pmod->smpl.t2);
	    gretl_xml_get_prop_as_unsigned_int(cur, "rseed", &pmod->smpl.rseed);
	} else if (!xmlStrcmp(cur->name, (XUC) "coeff")) {
	    pmod->coeff = gretl_xml_get_double_array(cur, doc, &n, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "sderr")) {
	    pmod->sderr = gretl_xml_get_double_array(cur, doc, &n, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "uhat")) {
	    pmod->uhat = gretl_xml_get_double_array(cur, doc, &n, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "yhat")) {
	    pmod->yhat = gretl_xml_get_double_array(cur, doc, &n, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "xpx")) {
	    /* note: allow this one to fail silently */
	    pmod->xpx = gretl_xml_get_double_array(cur, doc, &n, NULL);
	} else if (!xmlStrcmp(cur->name, (XUC) "vcv")) {
	    pmod->vcv = gretl_xml_get_double_array(cur, doc, &n, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "list")) {
	    pmod->list = gretl_xml_get_list(cur, doc, err);
	} else if (!xmlStrcmp(cur->name, (XUC) "tests")) {
	    *err = attach_model_tests_from_xml(pmod, cur);
	} else if (!xmlStrcmp(cur->name, (XUC) "params")) {
	    *err = attach_model_params_from_xml(cur, doc, pmod);
	} else if (!xmlStrcmp(cur->name, (XUC) "submask")) {
	    fprintf(stderr, "getting submask...\n");
	    *err = model_submask_from_xml(cur, doc, pmod);
	    fprintf(stderr, "got it\n");
	} else if (!xmlStrcmp(cur->name, (XUC) "missmask")) {
	    if (!gretl_xml_node_get_string(cur, doc, &pmod->missmask)) {
		*err = 1;
	    }
	} else if (!xmlStrcmp(cur->name, (XUC) "arinfo")) {
	    *err = arinfo_from_xml(cur, doc, pmod);
	} else if (!xmlStrcmp(cur->name, (XUC) "data-items")) {
	    *err = model_data_items_from_xml(cur, doc, pmod, fixopt);
	}

	if (*err) {
	    fprintf(stderr, "gretl_model_from_XML: block 3: err = %d reading '%s'\n",
		    *err, cur->name);
	    break;
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

    if (!*err && pmod->opt != 0) {
	/* compat for some options whose flags have been changed */
	if (pmod->opt & OPT_E) {
	    if (pmod->ci == GARCH) {
		/* standardized residuals */
		pmod->opt &= ~OPT_E;
		pmod->opt |= OPT_Z;
	    }
	}
	if (pmod->opt & OPT_W) {
	    if (pmod->ci == PANEL || pmod->ci == IVREG) {
		/* weights */
		pmod->opt &= ~OPT_W;
		pmod->opt |= OPT_H;
	    } else if (pmod->ci == DURATION) {
		/* Weibull */
		pmod->opt &= ~OPT_W;
		pmod->opt |= OPT_B;
	    }
	}
    }

    if (!*err && pmod->list != NULL) {
	maybe_convert_listsep(pmod->list, dset);
    }

    if (modname != NULL) {
	gretl_model_set_name(pmod, modname);
	free(modname);
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

MODEL *gretl_model_copy (MODEL *pmod)
{
    MODEL *new = malloc(sizeof *new);

#if MDEBUG
    fprintf(stderr, "gretl_model_copy: copying %p, allocated at %p\n",
	    pmod, new);
#endif

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

#define normality_test(ci, opt) (ci == MODTEST && (opt & OPT_N))

/**
 * command_ok_for_model:
 * @test_ci: index of command to be tested.
 * @opt: option for command to be tested.
 * @pmod: the gretl model in question.
 *
 * Check to see if the model-related command in question is
 * meaningful and acceptable in the context of the estimator
 * associated with @pmod. Note, though, that this function may
 * give a "false positive": to be quite sure, we may need to
 * know more about the model (e.g. specific options used).
 * See also model_test_ok().
 *
 * Returns: 1 if the command seems OK, otherwise 0.
 */

int command_ok_for_model (int test_ci, gretlopt opt,
			  const MODEL *pmod)
{
    int between = 0;
    int regular_panel = 0;
    int mci = pmod->ci;
    int ok = 1;

    if (mci == MIDASREG) {
	/* treat as a case of NLS */
	mci = NLS;
    }

    if (mci == NLS && test_ci == FCAST) {
	return 1;
    }

    if (mci == BIPROBIT) {
	/* can we support anything else? */
	return test_ci == RESTRICT;
    }

    if (test_ci == BKW) {
	/* most models should be OK */
	return pmod->ncoeff > 1 && (pmod->vcv != NULL || pmod->xpx != NULL);
    }

    if (NONLIST_MODEL(mci)) {
	return (test_ci == RESTRICT ||
		test_ci == TABPRINT ||
		(mci != MLE && normality_test(test_ci, opt)));
    }

    between = gretl_is_between_model(pmod);
    regular_panel = gretl_is_regular_panel_model(pmod);

    switch (test_ci) {

    case ADD:
	if (mci == ARMA || mci == GARCH ||
	    mci == HECKIT || mci == INTREG) {
	    ok = 0;
	} else if (mci == PANEL && (pmod->opt & OPT_B)) {
	    ok = 0;
	} else if (opt & (OPT_L | OPT_A)) {
	    /* --lm and --auto variants: OLS only */
	    ok = (mci == OLS);
        }
	break;

    case OMIT:
	if (mci == ARMA || mci == GARCH || mci == INTREG ||
	    mci == DPANEL) {
	    ok = 0;
	} else if (between) {
	    ok = 0;
	}
	break;

    case VIF:
	if (mci == IVREG || mci == ARMA || mci == GARCH ||
	    mci == PANEL || mci == DPANEL) {
	    ok = 0;
	}
	break;

    case EQNPRINT:
	if (mci == ARMA || mci == DPANEL ||
	    mci == HECKIT || mci == INTREG) {
	    ok = 0;
	}
	break;

    case MODTEST:
	if (opt & OPT_H) {
	    /* ARCH */
	    ok = (mci != ARCH && mci != GARCH);
	} else if (opt & OPT_C) {
	    /* common factor restriction */
	    ok = (mci == AR1);
	} else if (opt & OPT_D) {
	    /* x-sectional dependence */
	    ok = !between;
	} else if (opt & OPT_N) {
	    /* normality */
	    if (mci == LOGIT || mci == HECKIT || mci == DURATION) {
		/* POISSON, NEGBIN? */
		ok = 0;
	    }
	} else if (mci != OLS) {
	    if (mci == IVREG && (opt & (OPT_A | OPT_W))) {
		/* Autocorr. and H'sked. supported for IVREG */
		ok = 1;
	    } else if (mci == ARMA && (opt & OPT_A)) {
		ok = 1;
	    } else if (mci == PANEL && (opt & OPT_P)) {
		ok = 1;
	    } else if (regular_panel && (opt & OPT_A)) {
		ok = 1;
	    } else {
		ok = 0;
	    }
	}
	break;

    case CHOW:
    case CUSUM:
    case QLRTEST:
    case LEVERAGE:
    case RESET:
    case PANSPEC:
	/* OLS-only tests */
	ok = (mci == OLS);
	break;

    case RESTRICT:
	if (mci == LAD || mci == QUANTREG) {
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
 * @dset: dataset information.
 *
 * A more rigorous version of command_ok_for_model().  Use
 * this function if the extra information is available.
 *
 * Returns: 1 if the test command @ci (with possible option
 * @opt) is acceptable in the context of the model @pmod, and
 * the dataset described by @dset, otherwise 0.
 */

int model_test_ok (int ci, gretlopt opt, const MODEL *pmod,
		   const DATASET *dset)
{
    int ok = command_ok_for_model(ci, opt, pmod);

    /* for now we'll treat MIDASREG as a case of NLS */
    if (ci == MIDASREG) {
	ci = NLS;
    }

    if (ok && pmod->missmask != NULL) {
	/* can't do these with embedded missing obs? */
	if (gretl_is_regular_panel_model(pmod)) {
	    ; /* OK? */
	} else if (ci == CUSUM || ci == BDS ||
	    (ci == MODTEST && (opt & (OPT_A | OPT_H)))) {
	    ok = 0;
	}
    }

    if (ok && pmod->ncoeff == 1) {
	if (ci == COEFFSUM) {
	    ok = 0;
	} else if (pmod->ifc && ci == MODTEST) {
	    /* const only: rule out squares, logs, h'sked */
	    if (opt & (OPT_W | OPT_B | OPT_S | OPT_L)) {
		ok = 0;
	    }
	}
    }

    if (ci == MODTEST && (opt & OPT_A) &&
	gretl_is_regular_panel_model(pmod)) {
	return 1;
    }

    if (ok && !dataset_is_time_series(dset)) {
	/* time-series-only tests */
	if (ci == CUSUM || ci == QLRTEST || ci == BDS ||
	    (ci == MODTEST && (opt & (OPT_H | OPT_A)))) {
	    ok = 0;
	}
    }

    if (ok && !dataset_is_panel(dset)) {
	/* panel-only tests */
	if (ci == PANSPEC) {
	    ok = 0;
	} else if (ci == MODTEST && (opt & (OPT_P | OPT_D))) {
	    ok = 0;
	}
    }

    if (ok && pmod->ncoeff - pmod->ifc <= 1 && ci == VIF) {
	/* needs at least two independent vars */
	ok = 0;
    }

    if (ok && ci == MODTEST && (opt & OPT_C)) {
	/* common factor test */
	if (pmod->opt & OPT_P) {
	    /* ?? check what this means! */
	    ok = 0;
	}
    }

    if (ok && ci == RESTRICT) {
	if (gretl_model_get_int(pmod, "null-model")) {
	    ok = 0;
	}
    }

    return ok;
}

static int gretl_model_count = 0;

int get_model_count (void)
{
    return gretl_model_count;
}

void set_model_count (int c)
{
    gretl_model_count = c;
}

int model_count_plus (void)
{
    return ++gretl_model_count;
}

void model_count_minus (MODEL *pmod)
{
    if (pmod == NULL) {
	--gretl_model_count;
    } else if (pmod->ID > 0) {
	gretl_model_count = pmod->ID - 1;
	pmod->ID = 0;
    }
}

void set_model_id (MODEL *pmod, gretlopt opt)
{
    /* always record model estimation time */
     pmod->esttime = gretl_monotonic_time();

    /* Experimental: limit setting of sequential model
       number to "visible" models estimated in main
       script or session (not in functions).
    */
    if (opt & OPT_A) {
	/* An auxiliary model? Likely, but OPT_A has
	   special meaning for a few estimators
	*/
	if (pmod->ci != DPANEL && pmod->ci != GARCH &&
	    pmod->ci != ARMA) {
	    return;
	}
    }

    if ((opt & OPT_Q) && !(opt & OPT_W)) {
	/* --quiet and not --window, so "invisible" */
	return;
    } else if (gretl_function_depth() > 0) {
	/* model was estimated inside a function */
	return;
    }

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
				   const DATASET *dset)
{
    const char *vname;
    int i, v, vmax = 0;
    int gotsep = 0;

    if (pmod->ci == MLE || pmod->ci == GMM || pmod->list == NULL) {
	return 0;
    }

    for (i=1; i<=pmod->list[0]; i++) {
	v = pmod->list[i];
	if (v == LISTSEP) {
	    gotsep = 1;
	    continue;
	}
	if (v >= dset->v) {
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

    /* clustered standard errors? */
    vname = gretl_model_get_cluster_vname(pmod);
    if (vname != NULL) {
	v = current_series_index(dset, vname);
	if (v > vmax) vmax = v;
    }
    vname = gretl_model_get_cluster_vname2(pmod);
    if (vname != NULL) {
	v = current_series_index(dset, vname);
	if (v > vmax) vmax = v;
    }

    /* auxiliary variables for some model types */

    if (pmod->ci == WLS) {
	if (pmod->nwt > vmax) vmax = pmod->nwt;
    } else if (pmod->ci == INTREG) {
	v = gretl_model_get_int(pmod, "lovar");
	if (v > vmax) vmax = v;
	v = gretl_model_get_int(pmod, "hivar");
	if (v > vmax) vmax = v;
    } else if (COUNT_MODEL(pmod->ci)) {
	v = gretl_model_get_int(pmod, "offset_var");
	if (v > vmax) vmax = v;
    } else if (pmod->ci == DURATION) {
	v = gretl_model_get_int(pmod, "cens_var");
	if (v > vmax) vmax = v;
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
	int n = pmod->nobs;

	pmod->criterion[C_AIC] = -2.0 * pmod->lnL + 2.0 * k;
	pmod->criterion[C_BIC] = -2.0 * pmod->lnL + k * log(n);
	pmod->criterion[C_HQC] = -2.0 * pmod->lnL + 2 * k * log(log(n));
    }

    return err;
}

int model_use_zscore (const MODEL *pmod)
{
    if (pmod == NULL) {
	/* modprint */
	return 1;
    } else if (gretl_model_get_int(pmod, "dfcorr")) {
	/* override ASYMPTOTIC_MODEL if need be */
	return 0;
    } else if (pmod->ci == OLS && (pmod->opt & OPT_N)) {
	/* OLS with --no-df-corr */
	return 1;
    } else if (ASYMPTOTIC_MODEL(pmod->ci)) {
	return 1;
    } else if (pmod->ci == PANEL && (pmod->opt & (OPT_U | OPT_N))) {
	return 1;
    } else if ((pmod->opt & OPT_R) && libset_get_bool(ROBUST_Z)) {
	return 1;
    } else if (gretl_model_get_int(pmod, "asy")) {
        return 1;
    } else {
	return 0;
    }
}

double coeff_pval (int ci, double x, int df)
{
    double p = NADBL;

    if (!na(x)) {
	if (df == 0 || ASYMPTOTIC_MODEL(ci)) {
	    p = normal_pvalue_2(x);
	} else {
	    p = student_pvalue_2(df, x);
	}
    }

    return p;
}

double model_coeff_pval (const MODEL *pmod, double x)
{
    double p = NADBL;

    if (!na(x)) {
	if (model_use_zscore(pmod)) {
	    p = normal_pvalue_2(x);
	} else if (pmod->ci == AR) {
	    /* backward compat (?) */
	    int k = pmod->arinfo->arlist[0];
	    int dfd = pmod->dfd;

	    if (k > 1) {
		dfd += pmod->ncoeff - k;
	    }
	    p = student_pvalue_2(dfd, x);
	} else {
	    p = student_pvalue_2(pmod->dfd, x);
	}
    }

    return p;
}

/**
 * gretl_model_allocate_param_names:
 * @pmod: pointer to target model.
 * @k: number of strings to allocate.
 *
 * Allocate an array of @k strings to hold the names given to
 * the associated  coefficients, in a model where these strings
 * are not simply given by the names of the independent variables.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_model_allocate_param_names (MODEL *pmod, int k)
{
    pmod->params = strings_array_new_with_length(k, PNAMELEN);

    if (pmod->params == NULL) {
	pmod->errcode = E_ALLOC;
    }

    if (!pmod->errcode) {
	pmod->nparams = k;
    }

    return pmod->errcode;
}

/**
 * gretl_model_set_param_name:
 * @pmod: pointer to target model.
 * @i: 0-based index of value to set.
 * @name: string value to set.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_model_set_param_name (MODEL *pmod, int i, const char *name)
{
    if (pmod->params == NULL || i < 0 || i >= pmod->nparams ||
	name == NULL) {
	return E_DATA;
    } else {
	pmod->params[i][0] = '\0';
	if (strlen(name) >= PNAMELEN) {
	    strncat(pmod->params[i], name, PNAMELEN - 2);
	    strcat(pmod->params[i], "~");
	} else {
	    strncat(pmod->params[i], name, PNAMELEN - 1);
	}
	return 0;
    }
}

static int count_coeffs (int p, const char *pmask,
			 int q, const char *qmask,
			 int P, int Q, int r, int ifc)
{
    int i, n = P + Q + r + ifc;

    for (i=0; i<p; i++) {
	if (arma_included(pmask, i)) {
	    n++;
	}
    }

    for (i=0; i<q; i++) {
	if (arma_included(qmask, i)) {
	    n++;
	}
    }

    return n;
}

/**
 * gretl_model_add_arma_varnames:
 * @pmod: pointer to target model.
 * @dset: dataset information.
 * @yno: ID number of dependent variable.
 * @p: non-seasonal (max) AR order.
 * @q: non-seasonal (max) MA order.
 * @pmask: for masking out specific AR lags.
 * @qmask: for masking out specific MA lags.
 * @P: seasonal AR order.
 * @Q: seasonal MA order.
 * @r: number of exogenous regressors (excluding the constant).
 *
 * Composes a set of names to be given to the regressors in an
 * ARMA model.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_model_add_arma_varnames (MODEL *pmod, const DATASET *dset,
				   int yno, int p, int q,
				   const char *pmask, const char *qmask,
				   int P, int Q,
				   int r)
{
    int nullmod = 0;
    int nc, xstart;
    int i, j;

    nc = count_coeffs(p, pmask, q, qmask, P, Q, r, pmod->ifc);

    if (pmod->depvar != NULL) {
	free(pmod->depvar);
    }

    pmod->depvar = gretl_strdup(dset->varname[yno]);
    if (pmod->depvar == NULL) {
	pmod->errcode = E_ALLOC;
	return 1;
    }

    if (pmod->nparams > 0 && pmod->params != NULL) {
	for (i=0; i<pmod->nparams; i++) {
	    free(pmod->params[i]);
	}
	free(pmod->params);
    }

    nullmod = gretl_model_get_int(pmod, "null-model");

    if (nc == 0 && nullmod) {
	/* special case of null model */
	nc = 1;
    }

    pmod->params = strings_array_new_with_length(nc, VNAMELEN);
    if (pmod->params == NULL) {
	free(pmod->depvar);
	pmod->depvar = NULL;
	pmod->errcode = E_ALLOC;
	return 1;
    }

    pmod->nparams = nc;

    if (pmod->ifc || nullmod) {
	strcpy(pmod->params[0], dset->varname[0]);
	j = 1;
    } else {
	j = 0;
    }

    for (i=0; i<p; i++) {
	if (arma_included(pmask, i)) {
	    sprintf(pmod->params[j++], "phi_%d", i + 1);
	}
    }

    for (i=0; i<P; i++) {
	sprintf(pmod->params[j++], "Phi_%d", i + 1);
    }

    for (i=0; i<q; i++) {
	if (arma_included(qmask, i)) {
	    sprintf(pmod->params[j++], "theta_%d", i + 1);
	}
    }

    for (i=0; i<Q; i++) {
	sprintf(pmod->params[j++], "Theta_%d", i + 1);
    }

    xstart = arma_depvar_pos(pmod) + 1;

    for (i=0; i<r; i++) {
	strcpy(pmod->params[j++], dset->varname[pmod->list[xstart+i]]);
    }

    return 0;
}

/**
 * gretl_model_add_panel_varnames:
 * @pmod: pointer to target model.
 * @dset: dataset information.
 * @ulist: list of index numbers of cross-sectional
 * units included in the model.
 *
 * Composes a set of names to be given to the regressors in an
 * panel model.
 *
 * Returns: 0 on success, non-zero on error.
 */

int gretl_model_add_panel_varnames (MODEL *pmod, const DATASET *dset,
				    const int *ulist)
{
    int np = pmod->ncoeff;
    int i, j, v;

    pmod->depvar = gretl_strdup(dset->varname[pmod->list[1]]);
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
	if (v < dset->v) {
	    strcpy(pmod->params[i], dset->varname[v]);
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

/**
 * gretl_model_add_y_median:
 * @pmod: pointer to target model.
 * @y: array containing the dependent variable.
 *
 * Calculates the median of @y using the valid observations
 * with the model's sample range and attaches the median
 * to the model as data under the key %ymedian.
 *
 * Returns: 0 on success or error code on error.
 */

int gretl_model_add_y_median (MODEL *pmod, const double *y)
{
    int T = pmod->t2 - pmod->t1 + 1;
    double *sy, m;
    int t, n, ok, n2p;

    sy = malloc(T * sizeof *sy);

    if (sy == NULL) {
	return E_ALLOC;
    }

    n = 0;
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (pmod->uhat != NULL) {
	    ok = !na(pmod->uhat[t]);
	} else {
	    ok = !model_missing(pmod, t);
	}
	if (ok) {
	    sy[n++] = y[t];
	}
    }

    if (n == 0) {
	free(sy);
	return E_DATA;
    }

    qsort(sy, n, sizeof *sy, gretl_compare_doubles);

    n2p = (T = n / 2) + 1;
    m = (n % 2)? sy[n2p - 1] : 0.5 * (sy[T - 1] + sy[n2p - 1]);

    gretl_model_set_double(pmod, "ymedian", m);

    free(sy);

    return 0;
}

int gretl_model_add_normality_test (MODEL *pmod, double X2)
{
    ModelTest *test = model_test_new(GRETL_TEST_NORMAL);
    int err = 0;

    if (test != NULL) {
        model_test_set_teststat(test, GRETL_STAT_NORMAL_CHISQ);
        model_test_set_dfn(test, 2);
        model_test_set_value(test, X2);
        model_test_set_pvalue(test, chisq_cdf_comp(2, X2));
        maybe_add_test_to_model(pmod, test);
    } else {
        err = E_ALLOC;
    }

    return err;
}

ModelTest *gretl_model_get_test (MODEL *pmod, ModelTestType ttype)
{
    int i;

    for (i=0; i<pmod->ntests; i++) {
	if (pmod->tests[i].type == ttype) {
	    return &pmod->tests[i];
	}
    }

    return NULL;
}

static int xneq (double x, double y)
{
    double reldiff;

    if (x == 0.0) {
	reldiff = fabs(y);
    } else if (y == 0.0) {
	reldiff = fabs(x);
    } else if (x > y) {
	reldiff = fabs((x - y) / y);
    } else {
	reldiff = fabs((y - x) / x);
    }

    return reldiff > 1.5e-12;
}

/* try to tell if an OLS model with two independent variables is
   actually a quadratic model (x_2 = x_1^2). */

static int model_is_quadratic (const MODEL *pmod, const int *xlist,
			       const DATASET *dset)
{
    const double *x1 = dset->Z[xlist[2]];
    const double *x2 = dset->Z[xlist[3]];
    int t;

    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (!na(x1[t]) && xneq(x2[t], x1[t] * x1[t])) {
	    return 0;
	}
    }

    return 1;
}

/* we'll be a bit conservative here so as not to make
   fools of ourselves */

#define fitted_formula_ok(c) (c == OLS || c == WLS || \
                              c == HSK || c == IVREG || \
                              c == LAD || c == LOGISTIC)

/**
 * gretl_model_get_fitted_formula:
 * @pmod: pointer to target model.
 * @xvar: ID number of variable that _may_ be "x" in the model.
 * @dset: dataset information.
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
				      const DATASET *dset)
{
    const DATASET *mset;
    int *xlist = NULL;
    char *ret = NULL;

    if (xvar == 0 || pmod->ncoeff > 3 || !fitted_formula_ok(pmod->ci)) {
	return NULL;
    }

    xlist = gretl_model_get_x_list(pmod);
    if (xlist == NULL) {
	return NULL;
    }

    if (pmod->dataset != NULL) {
	mset = pmod->dataset;
    } else {
	mset = dset;
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
    } else if (!pmod->ifc && pmod->ncoeff == 1 && xvar == xlist[1]) {
	ret = gretl_strdup_printf("yformula: %g*x", pmod->coeff[0]);
    } else if (pmod->ifc && pmod->ncoeff == 2 && xvar == xlist[2]) {
	ret = gretl_strdup_printf("yformula: %g%s%g*x", pmod->coeff[0],
				  (pmod->coeff[1] >= 0)? "+" : "",
				  pmod->coeff[1]);
    } else if (pmod->ifc && pmod->ncoeff == 3 && xvar == xlist[2]) {
	if (model_is_quadratic(pmod, xlist, mset)) {
	    ret = gretl_strdup_printf("yformula: %g%s%g*x%s%g*x**2", pmod->coeff[0],
				      (pmod->coeff[1] >= 0)? "+" : "",
				      pmod->coeff[1],
				      (pmod->coeff[2] >= 0)? "+" : "",
				      pmod->coeff[2]);
	}
    }

    gretl_pop_c_numeric_locale();

    free(xlist);

    return ret;
}

/**
 * gretl_model_set_name:
 * @pmod: pointer to target model.
 * @name: the name to give the model.
 *
 * Sets the name of the given model; this is used in
 * printing the model and in displaying it in the
 * gretl GUI. Note that a model's name must be no more
 * than #MAXSAVENAME bytes in length, including the
 * terminating NUL byte; @name is truncated if it is
 * too long.
 */

void gretl_model_set_name (MODEL *pmod, const char *name)
{
    if (name == pmod->name) {
	return;
    }

    if (pmod->name == NULL) {
	pmod->name = calloc(1, MAXSAVENAME);
    }

    if (pmod->name != NULL) {
	*pmod->name = '\0';
	if (strlen(name) >= MAXSAVENAME) {
	    char *tmp = g_strdup(name);

	    strcpy(pmod->name, gretl_utf8_truncate(tmp, MAXSAVENAME-1));
	    g_free(tmp);
	} else {
	    strcpy(pmod->name, name);
	}
    }
}

/**
 * gretl_model_get_name:
 * @pmod: pointer to gretl model.
 *
 * Returns: the name that has been set for @pmod, if any.
 * Note that the value returned may be %NULL.
 */

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
 * @dset: dataset struct.
 * @err: location to receive error code (required).
 *
 * Retrieves a specified scalar statistic from @pmod:
 * @idx must be less than %M_SCALAR_MAX.
 *
 * Returns: the requested statistic, or #NADBL on failure,
 * in which case @err will contain a non-zero error code.
 */

double gretl_model_get_scalar (MODEL *pmod, ModelDataIndex idx,
			       DATASET *dset, int *err)
{
    double x = NADBL;

    if (pmod == NULL) {
	*err = E_BADSTAT;
	return x;
    }

    if (idx == M_GMMCRIT && pmod->ci != GMM) {
	*err = E_BADSTAT;
	return x;
    }

    switch (idx) {
    case M_ESS:
    case M_GMMCRIT:
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
    case M_DW:
	x = pmod->dw;
	break;
    case M_DWPVAL:
	x = get_DW_pvalue_for_model(pmod, dset, err);
	break;
    case M_FSTT:
	x = pmod->fstt;
	break;
    case M_CHISQ:
	x = pmod->chisq;
	break;
    default:
	break;
    }

    if (idx != M_DWPVAL && na(x)) {
	*err = E_BADSTAT;
    }

    return x;
}

/**
 * gretl_model_get_series:
 * @x: series to fill (must be of length dset->n).
 * @pmod: pointer to target model.
 * @dset: dataset information.
 * @idx: index for the series that is wanted.
 *
 * Retrieves from @pmod a copy of a specified series (for
 * example, regression residuals); @idx must be greater than
 * %M_SCALAR_MAX and less than %M_SERIES_MAX.
 *
 * Returns: 0 on success, non-zero code on error.
 */

int gretl_model_get_series (double *x, MODEL *pmod,
			    const DATASET *dset,
			    ModelDataIndex idx)
{
    const double *src = NULL;
    int t;

    if (pmod->t2 - pmod->t1 + 1 > dset->n ||
	model_sample_problem(pmod, dset)) {
	gretl_errmsg_set(
	       (idx == M_UHAT)?
	       _("Can't retrieve uhat: data set has changed") :
	       (idx == M_YHAT)?
	       _("Can't retrieve yhat: data set has changed") :
	       (idx == M_H)?
	       _("Can't retrieve ht: data set has changed") :
	       _("Can't retrieve series: data set has changed"));
	return E_BADSTAT;
    }

    if (pmod->ci == BIPROBIT && (idx == M_UHAT || idx == M_YHAT)) {
	return E_BADSTAT;
    }

    if (idx == M_UHAT) {
	src = pmod->uhat;
    } else if (idx == M_YHAT) {
	src = pmod->yhat;
    } else if (idx == M_LLT) {
	src = gretl_model_get_data(pmod, "llt");
    } else if (idx == M_AHAT) {
	src = gretl_model_get_data(pmod, "ahat");
    } else if (idx == M_H) {
	src = gretl_model_get_data(pmod, "garch_h");
    }

    if (src == NULL && idx != M_SAMPLE) {
	return E_BADSTAT;
    }

    /* allow for internal "just testing" usage */
    if (x == NULL) {
	return 0;
    }

    if (idx == M_SAMPLE) {
	for (t=0; t<dset->n; t++) {
	    if (t < pmod->t1 || t > pmod->t2) {
		x[t] = 0.0;
	    } else {
		x[t] = model_missing(pmod, t)? 0 : 1;
	    }
	}
    } else {
	for (t=0; t<dset->n; t++) {
	    if (t < pmod->t1 || t > pmod->t2) {
		x[t] = NADBL;
	    } else {
		x[t] = src[t];
	    }
	}
    }

    return 0;
}

static gretl_matrix *
model_get_estvec (const MODEL *pmod, int idx, int *err)
{
    const double *src;
    gretl_matrix *v = NULL;
    int i;

    if (gretl_model_get_data(pmod, "rq_tauvec")) {
	*err = E_BADSTAT;
	return NULL;
    }

    src = (idx == M_COEFF)? pmod->coeff : pmod->sderr;

    if (src == NULL) {
	*err = E_BADSTAT;
	return NULL;
    }

    v = gretl_column_vector_alloc(pmod->ncoeff);
    if (v == NULL) {
	*err = E_ALLOC;
    } else {
	for (i=0; i<pmod->ncoeff; i++) {
	    gretl_vector_set(v, i, src[i]);
	}
    }

    return v;
}

static gretl_matrix *
model_get_hatvec (const MODEL *pmod, int idx, int *err)
{
    gretl_matrix *v = NULL;
    const double *src = NULL;
    int t;

    if (idx == M_UHAT) {
	src = pmod->uhat;
    } else if (idx == M_YHAT) {
	src = pmod->yhat;
    } else if (idx == M_LLT) {
	src = gretl_model_get_data(pmod, "llt");
    }

    if (src == NULL) {
	*err = E_BADSTAT;
	return NULL;
    }

    /* FIXME: we reject NAs when creating a gretl_matrix -- but should
       we just substitute NaNs?
    */
    for (t=pmod->t1; t<=pmod->t2; t++) {
	if (na(src[t])) {
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
		gretl_vector_set(v, t - pmod->t1, src[t]);
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
	gretl_matrix_set_t1(v, pmod->t1);
	gretl_matrix_set_t2(v, pmod->t2);
	for (t=pmod->t1; t<=pmod->t2; t++) {
	    /* FIXME: is this indexation always right? */
	    gretl_vector_set(v, t - pmod->t1, mdata[t]);
	}
    }

    return v;
}

static gretl_matrix *model_get_rhovec (const MODEL *pmod, int *err)
{
    gretl_matrix *r = NULL;

    if (pmod->ci == AR1) {
	double x = gretl_model_get_double(pmod, "rho_gls");

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

static gretl_matrix *
model_get_intervals_matrix (const MODEL *pmod, int *err)
{
    gretl_matrix *m, *ci;

    ci = gretl_model_get_data(pmod, "coeff_intervals");
    if (ci == NULL) {
	*err = E_BADSTAT;
	return NULL;
    }

    m = gretl_matrix_copy(ci);

    if (m == NULL) {
	*err = E_ALLOC;
    }

    return m;
}

static gretl_matrix *model_get_special_test (const MODEL *pmod,
					     int type, int *err)
{
    gretl_matrix *r = NULL;
    int i, found = 0;

    for (i=0; i<pmod->ntests; i++) {
	found = (pmod->tests[i].type == type);
	if (found) {
	    r = gretl_vector_alloc(3);
	    if (r == NULL) {
		*err = E_ALLOC;
		return NULL;
	    }
	    r->val[0] = pmod->tests[i].value;
	    r->val[1] = pmod->tests[i].dfn;
	    r->val[2] = pmod->tests[i].pvalue;
	    break;
	}
    }

    if (!found) {
	*err = E_BADSTAT;
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
 * %M_SERIES_MAX and less than %M_MATRIX_MAX.
 *
 * Returns: the allocated matrix, or %NULL on failure,
 * in which case @err will contain a non-zero error code.
 */

gretl_matrix *gretl_model_get_matrix (MODEL *pmod, ModelDataIndex idx,
				      int *err)
{
    gretl_matrix *M = NULL;

    if (pmod == NULL) {
	*err = E_BADSTAT;
	return NULL;
    }

    if ((idx == M_UHAT || idx == M_YHAT) &&
	(pmod->ci == BIPROBIT || gretl_is_between_model(pmod))) {
	/* special: uhat and yhat are matrices */
	const char *utag = pmod->ci == BIPROBIT ? "bp_uhat" : "uhat";
	const char *ytag = pmod->ci == BIPROBIT ? "bp_yhat" : "yhat";
	gretl_matrix *P;

	if (idx == M_UHAT) {
	    P = gretl_model_get_data(pmod, utag);
	} else {
	    P = gretl_model_get_data(pmod, ytag);
	}

	if (P == NULL) {
	    *err = E_BADSTAT;
	} else {
	    M = gretl_matrix_copy(P);
	    if (M == NULL) {
		*err = E_ALLOC;
	    }
	}
	return M;
    }

    switch (idx) {
    case M_UHAT:
    case M_YHAT:
    case M_LLT:
	M = model_get_hatvec(pmod, idx, err);
	break;
    case M_COEFF:
    case M_SE:
	M = model_get_estvec(pmod, idx, err);
	break;
    case M_COEFF_CI:
	M = model_get_intervals_matrix(pmod, err);
	break;
    case M_VCV:
	M = gretl_vcv_matrix_from_model(pmod, NULL, err);
	break;
    case M_AHAT:
	if (gretl_model_get_data(pmod, "ahat") == NULL) {
	    *err = E_BADSTAT;
	} else {
	    M = model_get_special_vec(pmod, M_AHAT, err);
	}
	break;
    case M_H:
    case M_SIGMA:
	if (gretl_model_get_data(pmod, "garch_h") == NULL) {
	    *err = E_BADSTAT;
	} else {
	    M = model_get_special_vec(pmod, M_H, err);
	}
	break;
    case M_RHO:
	M = model_get_rhovec(pmod, err);
	break;
    case M_HAUSMAN:
	if (pmod->ci == IVREG) {
	    M = model_get_special_test(pmod, GRETL_TEST_IV_HAUSMAN, err);
	} else if (pmod->ci == PANEL) {
	    M = model_get_special_test(pmod, GRETL_TEST_PANEL_HAUSMAN, err);
	} else {
	    *err = E_BADSTAT;
	}
	break;
    case M_SARGAN:
	M = model_get_special_test(pmod, GRETL_TEST_SARGAN, err);
	break;
    case M_EHAT:
	M = gretl_model_get_data(pmod, "ehat");
	if (M == NULL) {
	    *err = E_BADSTAT;
	} else {
	    M = gretl_matrix_copy(M);
	    if (M == NULL) {
		*err = E_ALLOC;
	    }
	}
	break;
    case M_ODDSRATIOS:
	M = gretl_model_get_data(pmod, "oddsratios");
	if (M == NULL) {
	    *err = E_BADSTAT;
	} else {
	    M = gretl_matrix_copy(M);
	    if (M == NULL) {
		*err = E_ALLOC;
	    }
	}
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

static double get_vcv_element (MODEL *pmod, const char *s,
			       const DATASET *dset)
{
    char v1str[VNAMELEN], v2str[VNAMELEN];
    char fmt[16];
    int v1 = 0, v2 = 0;
    int i, j, k, gotit;
    double ret = NADBL;

    if (pmod == NULL || pmod->list == NULL) {
	return NADBL;
    }

    sprintf(fmt, "%%%d[^,],%%%ds", VNAMELEN-1, VNAMELEN-1);

    if (sscanf(s, fmt, v1str, v2str) != 2) {
	return NADBL;
    }

    v1 = gretl_model_get_param_number(pmod, dset, v1str);
    v2 = gretl_model_get_param_number(pmod, dset, v2str);

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
			      const DATASET *dset, int *err)
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

    if (idx == M_RHO) {
	int k = atoi(s);

	if (k == 1 && pmod->ci == AR1) {
	    x = gretl_model_get_double(pmod, "rho_gls");
	} else if (k == 1 && pmod->ci != AR) {
	    x = pmod->rho;
	} else if (pmod->arinfo == NULL ||
		   pmod->arinfo->arlist == NULL ||
		   pmod->arinfo->rho == NULL) {
	    *err = E_INVARG;
	} else {
	    vi = in_gretl_list(pmod->arinfo->arlist, k);
	    if (vi > 0) {
		x = pmod->arinfo->rho[vi-1];
	    } else {
		*err = E_INVARG;
	    }
	}
    } else if (idx == M_VCV) {
	x = get_vcv_element(pmod, s, dset);
	if (na(x)) {
	    *err = E_INVARG;
	}
    } else if (idx == M_COEFF || idx == M_SE) {
	vi = gretl_model_get_param_number(pmod, dset, s);
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

int exact_fit_check (const MODEL *pmod, PRN *prn)
{
    if (pmod->rsq == 1.0) {
	pputs(prn, _("The model exhibits an exact linear fit"));
	pputc(prn, '\n');
	return 1;
    }

    return 0;
}

void maybe_suppress_time_dummies (MODEL *pmod, int ndum)
{
    const char *s = get_optval_string(pmod->ci, OPT_D);

    if (s != NULL && !strcmp(s, "noprint")) {
	gretl_model_set_int(pmod, "skipdums", ndum);
    }
}
