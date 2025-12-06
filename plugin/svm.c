/*
 *  gretl -- Gnu Regression, Econometrics and Time-series Library
 *  Copyright (C) 2017 Allin Cottrell and Riccardo "Jack" Lucchetti
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

/* gretl interface to libsvm for machine learning */

#include "libgretl.h"
#include "libset.h"
#include "gretl_mt.h"
#include "version.h"

#ifdef HAVE_MPI
# include "gretl_mpi.h"
# include "gretl_foreign.h"
# include "matrix_extra.h"
#endif

#ifdef _OPENMP
# include <omp.h>
#endif

#include "svmlib.h"

#define DATA_DEBUG 0

typedef struct svm_problem sv_data;
typedef struct svm_node sv_cell;
typedef struct svm_model sv_model;
typedef struct svm_parameter sv_parm;
typedef struct sv_wrapper_ sv_wrapper;
typedef struct sv_grid_ sv_grid;
typedef struct grid_row_ grid_row;

/* the following need to be in sync with svmlib.h */

static const char *svm_type_names[] = {
    "C-SVC", "nu-SVC", "one-class", "eps-SVR", "nu-SVR",
    "C-rnk", NULL
};

static const char *kernel_type_names[] = {
    "linear", "polynomial", "RBF", "sigmoid",
    "stump", "perc", "laplace", "expo", NULL
};

enum {
    W_SAVEMOD = 1 << 0, /* saving a model? */
    W_LOADMOD = 1 << 1, /* loading a model? */
    W_QUIET   = 1 << 2, /* quiet operation? */
    W_SEARCH  = 1 << 3, /* doing parameter search? */
    W_STDFMT  = 1 << 4, /* use plain libsvm format for ranges? */
    W_SVPARM  = 1 << 5, /* saving tuned params to bundle? */
    W_FOLDVAR = 1 << 6, /* caller-supplied folds variable? */
    W_YSCALE  = 1 << 7, /* scaling the dependent var? */
    W_CONSEC  = 1 << 8, /* using consecutive folds? */
    W_REFOLD  = 1 << 9, /* folds differ across x-validation calls? */
    W_INTDEP  = 1 << 10, /* dependent variable is integer-valued */
    W_NOTRAIN = 1 << 11, /* not doing any training */
    W_MPI     = 1 << 12  /* use MPI for cross validation */
};

enum {
    REG_MSE = 1,
    REG_MAD,
    REG_ROUND_MAD,
    REG_ROUND_MISS
};

struct sv_wrapper_ {
    int auto_type;
    int flags;
    int scaling;
    int t1;
    int t2;
    int t2_train;
    int k;
    int nfold;    /* number of folds for x-validation */
    int predict;  /* are we doing (how much) prediction? */
    int nproc;    /* number of processes (MPI only) */
    int rank;     /* "this" MPI rank, if applicable */
    int regcrit;  /* regression optimality criterion */
    int do_probs; /* are we doing probability estimation? */
    double ymin;  /* min value of dependent variable */
    double ymax;  /* max value of dependent variable */
    gretl_matrix *ranges;
    char *ranges_outfile;
    char *ranges_infile;
    char *model_outfile;
    char *model_infile;
    char *data_outfile;
    char foldname[VNAMELEN];
    sv_grid *grid;
    char *plot;
    gretl_matrix *xdata;
    gretl_matrix *Ptrain;
    gretl_matrix *Ptest;
    double svr_sigma;
    int *flist;  /* array of folds IDs, if applicable */
    int *fsize;  /* array of fold sizes, if applicable */
    unsigned seed; /* FIXME */
};

struct sv_parm_info {
    const char *key;
    GretlType type;
};

enum {
    G_C,
    G_g,
    G_p
};

struct grid_row_ {
    double start;
    double stop;
    double step;
};

struct sv_grid_ {
    grid_row row[3];
    int null[3];
    int n[3];
    int linear[3];
};

#define N_PARMS 15

#define doing_regression(p) (p->svm_type == EPSILON_SVR || p->svm_type == NU_SVR)

static int gui_mode;

#ifdef HAVE_MPI
static PRN *worker_prn;
#endif

static void sv_wrapper_init (sv_wrapper *w, const DATASET *dset)
{
    w->auto_type = EPSILON_SVR;
    w->flags = 0;
    w->scaling = 1;
    w->t1 = dset->t1;
    w->t2 = dset->t2;
    w->t2_train = 0;
    w->k = 0;
    w->nfold = 0;
    w->predict = 2;
    w->nproc = 0;
    w->rank = -1;
    w->regcrit = 0;
    w->do_probs = 0;
    w->ymin = 0.0;
    w->ymax = 0.0;
    w->ranges = NULL;
    w->ranges_outfile = NULL;
    w->ranges_infile = NULL;
    w->model_outfile = NULL;
    w->model_infile = NULL;
    w->data_outfile = NULL;
    w->foldname[0] = '\0';
    w->grid = NULL;
    w->plot = NULL;
    w->xdata = NULL;
    w->Ptrain = NULL;
    w->Ptest = NULL;
    w->svr_sigma = NADBL;
    w->flist = NULL;
    w->fsize = NULL;
    w->seed = time(NULL);
}

static void sv_wrapper_free (sv_wrapper *w)
{
    gretl_matrix_free(w->ranges);
    free(w->ranges_outfile);
    free(w->ranges_infile);
    free(w->model_outfile);
    free(w->model_infile);
    free(w->data_outfile);
    free(w->grid);
    free(w->plot);
    gretl_matrix_free(w->xdata);
    gretl_matrix_free(w->Ptrain);
    gretl_matrix_free(w->Ptest);
    free(w->flist);
    free(w->fsize);
}

static void maybe_set_svm_seed (const sv_wrapper *w)
{
    if (w->flags & W_REFOLD) {
	/* set the seed just once */
	static int set;

	if (!set) {
	    gretl_alt_rand_set_seed(w->seed);
	    set = 1;
	}
    } else {
	/* set the seed on each x-validation run */
	gretl_alt_rand_set_seed(w->seed);
    }
}

static int loading_model (const sv_wrapper *w)
{
    return (w->flags & W_LOADMOD) || w->model_infile != NULL;
}

static int saving_model (const sv_wrapper *w)
{
    return (w->flags & W_SAVEMOD) || w->model_outfile != NULL;
}

static int uses_epsilon (const sv_parm *parm)
{
    return parm->svm_type == EPSILON_SVR;
}

static int uses_gamma (const sv_parm *parm)
{
    return parm->kernel_type == POLY ||
	parm->kernel_type == RBF ||
	parm->kernel_type == SIGMOID ||
	parm->kernel_type == LAPLACE ||
	parm->kernel_type == EXPO;
}

static int uses_coef0 (const sv_parm *parm)
{
    return parm->kernel_type == POLY ||
	parm->kernel_type == SIGMOID ||
	parm->kernel_type == STUMP ||
	parm->kernel_type == PERC;
}

static int uses_nu (const sv_parm *parm)
{
    return parm->svm_type == NU_SVC ||
	parm->svm_type == NU_SVR || parm->svm_type == ONE_CLASS;
}

static void set_sv_parm_defaults (sv_parm *parm)
{
    parm->svm_type = -1; /* mark as unknown for now */
    parm->kernel_type = RBF;
    parm->degree = 3;   /* for polynomial */
    parm->gamma = 0;    /* poly/RBF/sigmoid: default 1.0 / num_features */
    parm->coef0 = 0;    /* for use in kernel function */

    /* training-only variables */
    parm->cache_size = 1024;   /* cache size in MB */
    parm->eps = 0.001;         /* stopping criterion */
    parm->C = 1;               /* cost: for C_SVC, EPSILON_SVR and NU_SVR */
    parm->nr_weight = 0;       /* for C_SVC */
    parm->weight_label = NULL; /* for C_SVC */
    parm->weight = NULL;       /* for C_SVC */
    parm->nu = 0.5;            /* for NU_SVC, ONE_CLASS, and NU_SVR */
    parm->p = 0.1;             /* for EPSILON_SVR */
    parm->shrinking = 1;       /* use shrinking heuristics */
    parm->probability = 0;     /* do probability estimates */
}

static struct sv_parm_info pinfo[N_PARMS] = {
    { "svm_type",     GRETL_TYPE_INT },
    { "kernel_type",  GRETL_TYPE_INT },
    { "degree",       GRETL_TYPE_INT },
    { "gamma",        GRETL_TYPE_DOUBLE },
    { "coef0",        GRETL_TYPE_DOUBLE },
    { "cachesize",    GRETL_TYPE_DOUBLE },
    { "toler",        GRETL_TYPE_DOUBLE },
    { "C",            GRETL_TYPE_DOUBLE },
    { "nr_weight",    GRETL_TYPE_INT },
    { "weight_label", GRETL_TYPE_SERIES },
    { "weight",       GRETL_TYPE_SERIES },
    { "nu",           GRETL_TYPE_DOUBLE },
    { "epsilon",      GRETL_TYPE_DOUBLE },
    { "shrinking",    GRETL_TYPE_BOOL },
    { "probability",  GRETL_TYPE_BOOL }
};

static int is_sv_parm (const char *s)
{
    int i;

    for (i=0; i<N_PARMS; i++) {
	if (!strcmp(s, pinfo[i].key)) {
	    return 1;
	}
    }

    return 0;
}

static int type_value_from_string (const char *s, int i,
				   PRN *prn, int *err)
{
    int j;

    if (i == 0) {
	/* SVM type */
	for (j=0; svm_type_names[j] != NULL; j++) {
	    if (!g_ascii_strcasecmp(s, svm_type_names[j])) {
		return j;
	    }
	}
	pprintf(prn, _("%s: unrecognized SVM type\n"), s);
	*err = E_INVARG;
    } else {
	/* kernel type */
	for (j=0; kernel_type_names[j] != NULL; j++) {
	    if (!g_ascii_strcasecmp(s, kernel_type_names[j])) {
		return j;
	    }
	}
	pprintf(prn, _("%s: unrecognized kernel type\n"), s);
	*err = E_INVARG;
    }

    return -1;
}

static int set_or_store_sv_parm (sv_parm *parm, gretl_bundle *b,
				 int store, PRN *prn)
{
    void *elem[N_PARMS] = {
	&parm->svm_type,
	&parm->kernel_type,
	&parm->degree,
	&parm->gamma,
	&parm->coef0,
	&parm->cache_size,
	&parm->eps,
	&parm->C,
	&parm->nr_weight,
	&parm->weight_label,
	&parm->weight,
	&parm->nu,
	&parm->p,
	&parm->shrinking,
	&parm->probability
    };
    int got_kspec = 0;
    int i, ival;
    double xval;
    int err = 0;

    if (store) {
	/* we're supposed to store @parm in @b */
	for (i=0; i<N_PARMS && !err; i++) {
	    if (pinfo[i].type == GRETL_TYPE_DOUBLE) {
		xval = *(double *) elem[i];
		gretl_bundle_set_scalar(b, pinfo[i].key, xval);
	    } else {
		ival = *(int *) elem[i];
		gretl_bundle_set_int(b, pinfo[i].key, ival);
	    }
	}
	/* done */
	return 0;
    }

    set_sv_parm_defaults(parm);

    for (i=0; i<N_PARMS && !err; i++) {
	if (gretl_bundle_has_key(b, pinfo[i].key)) {
	    if (i < 2) {
		/* types of svm and kernel: int or string */
		GretlType type = 0;
		void *ptr = gretl_bundle_get_data(b, pinfo[i].key,
						  &type, NULL, &err);
		if (type == GRETL_TYPE_INT) {
		    *(int *) elem[i] = *(int *) ptr;
		} else if (type == GRETL_TYPE_DOUBLE) {
		    *(int *) elem[i] = *(double *) ptr;
		} else if (type == GRETL_TYPE_STRING) {
		    *(int *) elem[i] = type_value_from_string(ptr, i, prn, &err);
		} else {
		    fprintf(stderr, "parameter %d, bad option type\n", i);
		    err = E_TYPES;
		}
		if (i == 1) {
		    got_kspec = 1;
		}
	    } else if (i >= 8 && i <= 10) {
		pputs(prn, _("Sorry, svm weighting not handled yet\n"));
		err = E_INVARG;
	    } else if (pinfo[i].type == GRETL_TYPE_DOUBLE) {
		xval = gretl_bundle_get_scalar(b, pinfo[i].key, &err);
		if (!err) {
		    *(double *) elem[i] = xval;
		}
	    } else if (pinfo[i].type == GRETL_TYPE_INT ||
		       pinfo[i].type == GRETL_TYPE_BOOL) {
		ival = gretl_bundle_get_int(b, pinfo[i].key, &err);
		if (!err) {
		    if (pinfo[i].type == GRETL_TYPE_BOOL) {
			*(int *) elem[i] = (ival != 0);
		    } else {
			*(int *) elem[i] = ival;
		    }
		}
	    }
	}
    }

    if (!err && !got_kspec) {
	/* the user didn't specify a kernel type: we default to
	   Gaussian RBF above, but PERC seems to work better for
	   ranking
	*/
	if (parm->svm_type == C_RNK) {
	    parm->kernel_type = PERC;
	}
    }

    return err;
}

static int set_svm_parm (sv_parm *parm, gretl_bundle *b,
			 PRN *prn)
{
    return set_or_store_sv_parm(parm, b, 0, prn);
}

/* printing apparatus */

static PRN *svm_prn;

static void gretl_libsvm_print (const char *s)
{
    if (svm_prn != NULL) {
	pputs(svm_prn, s);
	gretl_flush(svm_prn);
    } else {
	fputs(s, stdout);
	fflush(stdout);
    }
}

static void gretl_libsvm_noprint (const char *s)
{
    return;
}

/* end printing apparatus */

static void gretl_destroy_svm_model (sv_model *model)
{
    if (model != NULL) {
	if (model->l > 0 && model->SV != NULL &&
	    model->SV[0] != NULL) {
	    free(model->SV[0]);
	}
	if (model->sv_coef != NULL) {
	    doubles_array_free(model->sv_coef, model->nr_class-1);
	}
	free(model->SV);
	free(model->rho);
	free(model->label);
	free(model->probA);
	free(model->probB);
	free(model->sv_indices);
	free(model->nSV);
	free(model);
    }
}

static int bundle_as_matrix (gretl_bundle *b, const char *key,
			     double *xvals, int n)
{
    gretl_matrix *m = gretl_matrix_alloc(n, 1);

    if (m == NULL) {
	return E_ALLOC;
    } else {
	memcpy(m->val, xvals, n * sizeof *xvals);
	gretl_bundle_donate_data(b, key, m, GRETL_TYPE_MATRIX, 0);
	return 0;
    }
}

static int bundle_as_list (gretl_bundle *b, const char *key,
			   int *ivals, int n)
{
    int *list = gretl_list_new(n);

    if (list == NULL) {
	return E_ALLOC;
    } else {
	memcpy(list + 1, ivals, n * sizeof *ivals);
	gretl_bundle_donate_data(b, key, list, GRETL_TYPE_LIST, 0);
	return 0;
    }
}

static int svm_model_save_to_bundle (const sv_model *model,
				     gretl_bundle *b)
{
    const sv_parm *parm = &model->param;
    gretl_matrix *m = NULL;
    int i, j, nc, l, ntr;
    int err = 0;

    gretl_bundle_void_content(b);

    gretl_bundle_set_int(b, "svm_type", parm->svm_type);
    gretl_bundle_set_int(b, "kernel_type", parm->kernel_type);

    if (parm->kernel_type == POLY) {
	gretl_bundle_set_int(b, "degree", parm->degree);
    }

    if (uses_gamma(parm)) {
	gretl_bundle_set_scalar(b, "gamma", parm->gamma);
    }

    if (uses_coef0(parm)) {
	gretl_bundle_set_scalar(b, "coef0", parm->coef0);
    }

    nc = model->nr_class;
    l = model->l;

    gretl_bundle_set_int(b, "nr_class", nc);
    gretl_bundle_set_int(b, "l", l);

    /* number of triangular elements */
    ntr = nc * (nc - 1) / 2;

    bundle_as_matrix(b, "rho", model->rho, ntr);

    if (model->label != NULL) {
	bundle_as_list(b, "label", model->label, nc);
    }
    if (model->probA != NULL) {
	bundle_as_matrix(b, "probA", model->probA, ntr);
    }
    if (model->probB != NULL) {
	bundle_as_matrix(b, "probB", model->probB, ntr);
    }
    if (model->nSV != NULL) {
	bundle_as_list(b, "nr_sv", model->nSV, nc);
    }

    /* store the SVs */

    m = gretl_matrix_alloc(l, nc - 1);
    if (m != NULL) {
	for (j=0; j<nc-1; j++) {
	    /* j = class index */
	    for (i=0; i<l; i++) {
		/* i = row index */
		gretl_matrix_set(m, i, j, model->sv_coef[j][i]);
	    }
	}
	gretl_bundle_donate_data(b, "sv_coef", m,
				 GRETL_TYPE_MATRIX, 0);
    }

    {
	gretl_array *aidx, *avec = NULL;
	gretl_matrix *vec;
	int *idx;
	int n_elements = 0;

	aidx = gretl_array_new(GRETL_TYPE_LISTS, l, &err);
	if (!err) {
	    avec = gretl_array_new(GRETL_TYPE_MATRICES, l, &err);
	}

	for (i=0; i<l && !err; i++) {
	    const sv_cell *p = model->SV[i];
	    int k, ni = 0;

	    /* count the nodes on this row */
	    while (p->index != -1) {
		ni++;
		p++;
	    }
	    idx = gretl_list_new(ni);
	    vec = gretl_matrix_alloc(1, ni);
	    if (idx == NULL || vec == NULL) {
		err = E_ALLOC;
		break;
	    }
	    p = model->SV[i]; /* reset */
	    for (k=0; k<ni; k++) {
		idx[k+1] = p[k].index;
		vec->val[k] = p[k].value;
	    }
	    gretl_array_set_list(aidx, i, idx, 0);
	    gretl_array_set_matrix(avec, i, vec, 0);
	    n_elements += ni + 1;
	}

	if (err) {
	    gretl_array_destroy(aidx);
	    gretl_array_destroy(avec);
	} else {
	    gretl_bundle_set_int(b, "n_elements", n_elements);
	    gretl_bundle_donate_data(b, "SV_indices", aidx,
				     GRETL_TYPE_ARRAY, 0);
	    gretl_bundle_donate_data(b, "SV_vecs", avec,
				     GRETL_TYPE_ARRAY, 0);
	}
    }

    if (err) {
	gretl_bundle_void_content(b);
    }

    return err;
}

static void save_results_to_bundle (const sv_parm *parm,
				    sv_wrapper *w,
				    gretl_bundle *b)
{
    if (w->xdata != NULL) {
	int ns = w->xdata->cols;
	char **S = strings_array_new(ns);

	if (S != NULL) {
	    S[0] = gretl_strdup("C");
	    S[1] = gretl_strdup("gamma");
	    if (ns == 4) {
		const char *s = uses_epsilon(parm) ? "epsilon" : "nu";

		S[2] = gretl_strdup(s);
		S[3] = gretl_strdup("score");
	    } else {
		S[2] = gretl_strdup("score");
	    }
	    gretl_matrix_set_colnames(w->xdata, S);
	}
	gretl_bundle_donate_data(b, "xvalid_results", w->xdata,
				 GRETL_TYPE_MATRIX, 0);
	w->xdata = NULL;
    }
}

/* If we did probability estimation we should have results
   in w->Ptrain and/or w->Ptest (in the classification
   case) or a single Laplace scale value in w->svr_sigma
   (in the regression case). We now stuff this info into
   the courier bundle supplied by the caller.
*/

static void save_probs_to_bundle (sv_wrapper *w,
				  gretl_bundle *b)
{
    if (w->Ptrain != NULL) {
	gretl_bundle_donate_data(b, "Ptrain", w->Ptrain,
				 GRETL_TYPE_MATRIX, 0);
	w->Ptrain = NULL;
    }
    if (w->Ptest != NULL) {
	gretl_bundle_donate_data(b, "Ptest", w->Ptest,
				 GRETL_TYPE_MATRIX, 0);
	w->Ptest = NULL;
    }
    if (!na(w->svr_sigma)) {
	gretl_bundle_set_scalar(b, "svr_sigma", w->svr_sigma);
    }
}

static double *array_from_bundled_matrix (gretl_bundle *b,
					  const char *key,
					  int required,
					  int *err)
{
    double *ret = NULL;

    if (*err) return NULL;

    if (gretl_bundle_has_key(b, key)) {
	gretl_matrix *m = gretl_bundle_get_matrix(b, key, err);

	if (m != NULL) {
	    int n = m->rows * m->cols;

	    ret = malloc(n * sizeof *ret);
	    if (ret == NULL) {
		*err = E_ALLOC;
	    } else {
		memcpy(ret, m->val, n * sizeof *ret);
	    }
	}
    } else if (required) {
	gretl_errmsg_sprintf(_("svm model: required matrix %s was not found"), key);
	*err = E_DATA;
    }

    return ret;
}

static int *array_from_bundled_list (gretl_bundle *b,
				     const char *key,
				     int required,
				     int *err)
{
    int *ret = NULL;

    if (*err) return NULL;

    if (gretl_bundle_has_key(b, key)) {
	int *list = gretl_bundle_get_list(b, key, err);

	if (list != NULL) {
	    int n = list[0];

	    ret = malloc(n * sizeof *ret);
	    if (ret == NULL) {
		*err = E_ALLOC;
	    } else {
		memcpy(ret, list + 1, n * sizeof *ret);
	    }
	}
    } else if (required) {
	gretl_errmsg_sprintf(_("svm model: required list %s was not found"), key);
	*err = E_DATA;
    }

    return ret;
}

static sv_model *svm_model_from_bundle (gretl_bundle *b,
					int *err)
{
    struct svm_parameter *parm;
    sv_model *model;
    gretl_matrix *m;
    int n_elements = 0;
    int nc = 0, l = 0;
    int i, j;

    model = malloc(sizeof *model);
    if (model == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    memset(model, 0, sizeof *model);
    parm = &model->param;
    *err = set_svm_parm(parm, b, NULL);

    if (!*err) {
	nc = model->nr_class = gretl_bundle_get_int(b, "nr_class", err);
	l = model->l = gretl_bundle_get_int(b, "l", err);
	n_elements = gretl_bundle_get_int(b, "n_elements", err);
	if (nc <= 0 || l <= 0 || n_elements <= 0) {
	    *err = E_DATA;
	}
    }

    model->rho = array_from_bundled_matrix(b, "rho", 1, err);

    /* not always present */
    model->label = array_from_bundled_list(b, "label", 0, err);
    model->probA = array_from_bundled_matrix(b, "probA", 0, err);
    model->probB = array_from_bundled_matrix(b, "probB", 0, err);
    model->nSV = array_from_bundled_list(b, "nr_sv", 0, err);

    /* load the SVs */

    if (!*err) {
	m = gretl_bundle_get_matrix(b, "sv_coef", err);
	if (m == NULL) {
	    *err = E_DATA;
	} else {
	    model->sv_coef = doubles_array_new(nc-1, l);
	    if (model->sv_coef == NULL) {
		*err = E_ALLOC;
	    } else {
		double *val = m->val;

		for (j=0; j<nc-1; j++) {
		    memcpy(model->sv_coef[j], val, l * sizeof *val);
		    val += l;
		}
	    }
	}
    }

    if (!*err) {
	sv_cell *p = NULL, *x_space = NULL;
	gretl_array *aidx = NULL;
	gretl_array *avec = NULL;
	gretl_matrix *vec;
	int *idx;

	model->SV = malloc(l * sizeof *model->SV);
	if (model->SV == NULL) {
	    *err = E_ALLOC;
	} else {
	    x_space = malloc(n_elements * sizeof *x_space);
	    if (x_space == NULL) {
		*err = E_ALLOC;
	    } else {
		model->SV[0] = p = x_space;
	    }
	}

	if (!*err) {
	    aidx = gretl_bundle_get_array(b, "SV_indices", NULL);
	    avec = gretl_bundle_get_array(b, "SV_vecs", NULL);
	    if (gretl_array_get_type(aidx) != GRETL_TYPE_LISTS ||
		gretl_array_get_type(avec) != GRETL_TYPE_MATRICES) {
		*err = E_DATA;
	    }
	}

	for (i=0; i<l && !*err; i++) {
	    int ni;

	    model->SV[i] = p;
	    idx = gretl_array_get_element(aidx, i, NULL, err);
	    vec = gretl_array_get_element(avec, i, NULL, err);
	    ni = idx[0];
	    for (j=0; j<ni; j++) {
		p[j].index = idx[j+1];
		p[j].value = vec->val[j];
	    }
	    /* add -1 sentinel */
	    p[j].index = -1;
	    p += ni + 1;
	}
    }

    if (*err) {
	gretl_destroy_svm_model(model);
	model = NULL;
    }

    return model;
}

/* can use for testing against svm-scale */

static int write_ranges (sv_wrapper *w)
{
    const char *fname;
    int libsvm_format = 0;
    double lo, hi;
    int i, idx, vi;
    FILE *fp;

    fname = gretl_maybe_switch_dir(w->ranges_outfile);
    fp = gretl_fopen(fname, "wb");
    if (fp == NULL) {
	return E_FOPEN;
    }

    if (w->flags & W_STDFMT) {
	libsvm_format = 1;
    }

    gretl_push_c_numeric_locale();

    if (libsvm_format) {
	fprintf(fp, "x\n%d %d\n",
		(int) gretl_matrix_get(w->ranges, 0, 0),
		(int) gretl_matrix_get(w->ranges, 0, 1));
    } else {
	fprintf(fp, "x\n%d %d %d\n",
		(int) gretl_matrix_get(w->ranges, 0, 0),
		(int) gretl_matrix_get(w->ranges, 0, 1),
		(int) gretl_matrix_get(w->ranges, 0, 2));
    }

    for (i=1; i<w->ranges->rows; i++) {
	idx = gretl_matrix_get(w->ranges, i, 0);
	lo  = gretl_matrix_get(w->ranges, i, 1);
	hi  = gretl_matrix_get(w->ranges, i, 2);
	if (libsvm_format) {
	    fprintf(fp, "%d %.16g %.16g\n", idx, lo, hi);
	} else {
	    vi = gretl_matrix_get(w->ranges, i, 3);
	    fprintf(fp, "%d %.16g %.16g %d\n", idx, lo, hi, vi);
	}
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return 0;
}

static int read_ranges (sv_wrapper *w)
{
    const char *fname;
    FILE *fp;
    char line[512];
    double lo, hi, j;
    int read_lims = 0;
    int ncols = 4;
    int i, vi, idx, n = 0;
    int err = 0;

    fname = gretl_maybe_switch_dir(w->ranges_infile);
    fp = gretl_fopen(fname, "rb");
    if (fp == NULL) {
	return E_FOPEN;
    }

    gretl_push_c_numeric_locale();

    while (fgets(line, sizeof line, fp) && !err) {
	if (*line == 'x') {
	    read_lims = 1;
	    continue;
	}
	if (read_lims) {
	    if (ncols == 3) {
		n = sscanf(line, "%lf %lf\n", &lo, &hi);
	    } else {
		n = sscanf(line, "%lf %lf %lf\n", &lo, &hi, &j);
	    }
	    if (n != ncols - 1) {
		err = E_DATA;
	    }
	    read_lims = 0;
	} else if (!string_is_blank(line)) {
	    n++;
	}
    }

    w->ranges = gretl_matrix_alloc(n+1, ncols);
    if (w->ranges == NULL) {
	err = E_ALLOC;
    } else {
	gretl_matrix_set(w->ranges, 0, 0, lo);
	gretl_matrix_set(w->ranges, 0, 1, hi);
	gretl_matrix_set(w->ranges, 0, 2, j);
	if (ncols == 4) {
	    gretl_matrix_set(w->ranges, 0, 3, 0);
	}
	rewind(fp);
	i = 1;
    }

    while (fgets(line, sizeof line, fp) && !err) {
	if (*line == 'x') {
	    if (fgets(line, sizeof line, fp) == NULL) {
		err = E_DATA;
		break;
	    } else {
		continue;
	    }
	}
	if (ncols == 3) {
	    n = sscanf(line, "%d %lf %lf\n", &idx, &lo, &hi);
	} else {
	    n = sscanf(line, "%d %lf %lf %d\n", &idx, &lo, &hi, &vi);
	}
	if (n != ncols) {
	    err = E_DATA;
	} else {
	    gretl_matrix_set(w->ranges, i, 0, idx);
	    gretl_matrix_set(w->ranges, i, 1, lo);
	    gretl_matrix_set(w->ranges, i, 2, hi);
	    if (ncols == 4) {
		gretl_matrix_set(w->ranges, i, 3, vi);
	    }
	    i++;
	}
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return err;
}

/* can use for testing against svm-scale */

static int write_problem (sv_data *p, sv_wrapper *w)
{
    const char *fname;
    FILE *fp;
    int i, t, idx;
    double val;

    fname = gretl_maybe_switch_dir(w->data_outfile);
    fp = gretl_fopen(fname, "wb");
    if (fp == NULL) {
	return E_FOPEN;
    }

    gretl_push_c_numeric_locale();

    for (t=0; t<p->l; t++) {
	fprintf(fp, "%g ", p->y[t]);
	for (i=0; i<w->k; i++) {
	    idx = p->x[t][i].index;
	    val = p->x[t][i].value;
	    if (val != 0) {
		fprintf(fp, "%d:%g ", idx, val);
	    }
	}
	fputc('\n', fp);
    }

    gretl_pop_c_numeric_locale();

    fclose(fp);

    return 0;
}

static void gretl_sv_data_destroy (sv_data *p, sv_cell *x_space)
{
    if (p != NULL) {
	free(p->y);
	free(p->x);
	free(p);
    }
    if (x_space != NULL) {
	free(x_space);
    }
}

static sv_data *gretl_sv_data_alloc (int T, int k,
				     sv_cell **px_space,
				     int *err)
{
    sv_data *p = malloc(sizeof *p);

    if (p != NULL) {
	p->l = T;
	p->y = malloc(T * sizeof *p->y);
	p->x = malloc(T * sizeof *p->x);
	if (p->y == NULL || p->x == NULL) {
	    *err = E_ALLOC;
	} else {
	    /* we need an extra cell on each row to hold a
	       sentinel index value of -1
	    */
	    *px_space = malloc(T * (k+1) * sizeof(sv_cell));
	    if (*px_space == NULL) {
		*err = E_ALLOC;
	    }
	}
	if (*err) {
	    gretl_sv_data_destroy(p, NULL);
	    p = NULL;
	}
    } else {
	*err = E_ALLOC;
    }

    return p;
}

/* initial discovery of data ranges using the training data */

static int get_data_ranges (const int *list,
			    const DATASET *dset,
			    sv_wrapper *w)
{
    const double *x;
    double xmin, xmax;
    int lmax = list[0];
    int i, j, vi;
    int err = 0;

    if (w->flags & W_FOLDVAR) {
	/* the last series is not a regressor */
	lmax--;
    }

    w->ranges = gretl_matrix_alloc(lmax, 4);
    if (w->ranges == NULL) {
	return E_ALLOC;
    }

    /* scaling limits */
    xmin = w->scaling == 2 ? 0 : -1;
    gretl_matrix_set(w->ranges, 0, 0, xmin); /* lower */
    gretl_matrix_set(w->ranges, 0, 1, 1);    /* upper */

    /* padding */
    gretl_matrix_set(w->ranges, 0, 2, 0);
    gretl_matrix_set(w->ranges, 0, 3, 0);

    j = 0;
    for (i=2; i<=lmax; i++) {
	vi = list[i];
	x = dset->Z[vi];
	gretl_minmax(dset->t1, dset->t2, x, &xmin, &xmax);
	if (xmin != xmax) {
	    j++;
	    gretl_matrix_set(w->ranges, j, 0, j);
	    gretl_matrix_set(w->ranges, j, 1, xmin);
	    gretl_matrix_set(w->ranges, j, 2, xmax);
	    gretl_matrix_set(w->ranges, j, 3, vi);
	} else {
	    fprintf(stderr, "training data: dropping var %d (%s)\n",
		    vi, dset->varname[vi]);
	}
    }

    if (w->flags & W_YSCALE) {
	/* we'll arrange to scale y onto [-1, +1] */
	gretl_minmax(dset->t1, dset->t2, dset->Z[list[1]],
		     &w->ymin, &w->ymax);
    }

    /* record number of rows actually occupied, which
       could be less than the number allocated
    */
    gretl_matrix_set(w->ranges, 0, 2, j + 1);

    return err;
}

/* The following requires a gretl-special "ranges" matrix
   with 4 columns, including the dataset IDs of the series.
*/

static int check_test_data (const DATASET *dset,
			    gretl_matrix *ranges,
			    int k)
{
    double xmin, xmax;
    int i, n, vi;
    int err = 0;

    n = 0;
    for (i=1; i<=k; i++) {
	vi = gretl_matrix_get(ranges, i, 3);
	gretl_minmax(dset->t1, dset->t2, dset->Z[vi], &xmin, &xmax);
	if (xmin != xmax) {
	    n++;
	} else {
	    fprintf(stderr, "test data: dropping var %d (%s)\n",
		    vi, dset->varname[vi]);
	    /* arrange to exclude this variable by setting the
	       record of its series ID to zero */
	    gretl_matrix_set(ranges, i, 3, 0);
	}
    }

    if (n != k) {
	fprintf(stderr, "test data: number of usable variables (%d) "
		"differs from training data (%d)\n", n, k);
    }

    return err;
}

/* apply scaling as per the svm-scale binary */

static double scale_x (double val, double lo, double hi,
		       double scalemin, double scalemax)
{
    if (val == lo) {
	val = scalemin;
    } else if (val == hi) {
	val = scalemax;
    } else {
	val = scalemin + (scalemax - scalemin) *
	    (val - lo) / (hi - lo);
    }

    return val;
}

static double scale_y (double y, sv_wrapper *w)
{
    return -1 + 2 * (y - w->ymin) / (w->ymax - w->ymin);
}

static double unscale_y (double y, sv_wrapper *w)
{
    return w->ymin + (w->ymax - w->ymin) * (y + 1) / 2.0;
}

static int sv_data_fill (sv_data *prob,
			 sv_cell *x_space,
			 sv_wrapper *w,
			 const int *list,
			 const DATASET *dset,
			 int pass)
{
    double scalemin, scalemax;
    double yt, xit, xmin, xmax;
    int i, j, s, t, idx;
    int vf = 0, pos = 0;
    int vi = list[1];
    int all_ints = 0;
    int k = w->k;

#if DATA_DEBUG
    fprintf(stderr, "sv_data_fill, pass %d\n", pass);
#endif

    /* deal with the LHS variable */
    if (pass == 1) {
	if (gretl_isdummy(dset->t1, dset->t2, dset->Z[vi]) ||
	    series_is_coded(dset, vi)) {
	    /* classification, not regression */
	    w->auto_type = C_SVC;
	}
	all_ints = 1;
    }

    for (i=0, t=dset->t1; t<=dset->t2; t++, i++) {
	yt = dset->Z[vi][t];
        if (w->flags & W_YSCALE) {
	    prob->y[i] = scale_y(yt, w);
	} else {
	    prob->y[i] = yt;
	}
	if (all_ints && prob->y[i] != floor(prob->y[i])) {
	    all_ints = 0;
	}
    }

    if (pass == 1 && (w->flags & W_FOLDVAR)) {
	w->flist = gretl_list_new(prob->l);
	if (w->flist != NULL) {
	    /* record the ID of the folds series */
	    vf = list[list[0]];
	}
    }

    if (pass == 1 && all_ints) {
	w->flags |= W_INTDEP;
    }

    /* retrieve the global x-scaling limits */
    scalemin = gretl_matrix_get(w->ranges, 0, 0);
    scalemax = gretl_matrix_get(w->ranges, 0, 1);

#if DATA_DEBUG
    fprintf(stderr, "scalemin = %g, scalemax = %g\n", scalemin, scalemax);
#endif

    /* write the scaled x-data into the problem struct */
    for (s=0, t=dset->t1; t<=dset->t2; t++, s++) {
	if (vf > 0) {
	    /* record the specified fold for this obs */
	    w->flist[s+1] = dset->Z[vf][t];
	}
	prob->x[s] = &x_space[pos];
	j = 0;
	for (i=1; i<=k; i++) {
	    if (w->ranges->cols == 4) {
		vi = (int) gretl_matrix_get(w->ranges, i, 3);
		if (vi <= 0) {
		    /* may happen when we get to the test data */
		    continue;
		}
	    }
	    idx = (int) gretl_matrix_get(w->ranges, i, 0);
	    xmin = gretl_matrix_get(w->ranges, i, 1);
	    xmax = gretl_matrix_get(w->ranges, i, 2);
	    xit = dset->Z[vi][t];
	    if (na(xit)) {
		/* catch missing values */
		fprintf(stderr, "skipping NA for var %d, obs %d\n", vi, t+1);
		xit = 0;
	    } else if (w->scaling != 0) {
		xit = scale_x(xit, xmin, xmax, scalemin, scalemax);
	    }
	    if (xit == 0) {
		/* fprintf(stderr, "skipping a 0 data value (var %d)\n", vi); */
		continue;
	    }
	    prob->x[s][j].index = idx;
	    prob->x[s][j].value = xit;
	    pos++;
	    j++;
	}
	/* end-of-row sentinel */
	prob->x[s][j].index = -1;
	prob->x[s][j].value = 0;
	pos++;
    }

    return 0;
}

/* apparatus for sorting labels into ascending order */

struct lsort {
    int val;
    int pos;
};

static int ls_compare (const void *a, const void *b)
{
    const struct lsort *pa = a;
    const struct lsort *pb = b;

    return (pa->val > pb->val) - (pa->val < pb->val);
}

/* Allocate a matrix (either w->Ptrain or w->Ptest)
   into which we'll write the per-outcome probabilities.
   Set the sorted labels as column names for this
   matrix, and return in location @pls the sorting
   info that will enable us to write the probabilities
   into the correct corresponding columns.
*/

static gretl_matrix *get_probs_matrix (sv_wrapper *w,
				       sv_model *model,
				       struct lsort **pls,
				       int training)
{
    gretl_matrix **targ;
    int nc = model->nr_class;
    int T, pt1 = 0, pt2 = 0;

    if (training) {
	targ = &w->Ptrain;
	pt1 = w->t1;
	pt2 = w->t2_train;
    } else {
	targ = &w->Ptest;
	pt1 = w->t2_train + 1;
	pt2 = w->t2;
    }

    T = pt2 - pt1 + 1;
    *targ = gretl_matrix_alloc(T, nc);
    if (*targ == NULL) {
	return NULL;
    }

    gretl_matrix_set_t1(*targ, pt1);
    gretl_matrix_set_t2(*targ, pt2);

    if (model->label != NULL) {
	char **S = strings_array_new(nc);
	struct lsort *ls = malloc(nc * sizeof *ls);

	if (S != NULL && ls != NULL) {
	    char labnum[32];
	    int i;

	    for (i=0; i<nc; i++) {
		ls[i].val = model->label[i];
		ls[i].pos = i;
	    }
	    qsort(ls, nc, sizeof *ls, ls_compare);
	    for (i=0; i<nc; i++) {
		sprintf(labnum, "%d", ls[i].val);
		S[i] = gretl_strdup(labnum);
	    }
	    gretl_matrix_set_colnames(*targ, S);
	    *pls = ls;
	}
    }

    return *targ;
}

static int real_svm_predict (double *yhat,
			     sv_data *prob,
			     sv_wrapper *w,
			     sv_model *model,
			     int training,
			     const DATASET *dset,
			     PRN *prn)
{
    const char *datastr;
    gretl_matrix *P = NULL;
    struct lsort *ls = NULL;
    int n_correct = 0;
    int regression = 0;
    int nr_class = 0;
    int get_sigma = 0;
    double ymean = 0.0;
    double TSS = 0.0;
    double SSR = 0.0;
    double MAD = 0.0;
    int misses = 0;
    double dev, yhi, yi;
    double *pi = NULL;
    sv_cell *x;
    int i, j, err = 0;

    if (model->param.svm_type == EPSILON_SVR ||
	model->param.svm_type == NU_SVR) {
	regression = 1;
	if (w->flags & W_YSCALE) {
	    for (i=0; i<prob->l; i++) {
		ymean += unscale_y(prob->y[i], w);
	    }
	    ymean /= prob->l;
	} else {
	    ymean = gretl_mean(0, prob->l - 1, prob->y);
	}
    }

    if (model->param.probability) {
	if (model->probA == NULL) {
	    fprintf(stderr, "probability requested but no probA!\n");
	    w->do_probs = 0;
	} else if (regression) {
	    get_sigma = training;
	} else {
	    /* classification */
	    nr_class = svm_get_nr_class(model);
	    P = get_probs_matrix(w, model, &ls, training);
	    if (P == NULL) {
		err = E_ALLOC;
	    } else {
		pi = malloc(nr_class * sizeof *pi);
		w->do_probs = 1;
	    }
	}
    }

    if (err) {
	return err;
    }

    if (w->do_probs) {
	maybe_set_svm_seed(w);
    }

    pprintf(prn, _("Calling prediction function (this may take a while)\n"));
    gretl_flush(prn);
    for (i=0; i<prob->l; i++) {
	x = prob->x[i];
	if (w->do_probs) {
	    yhi = svm_predict_probability(model, x, pi);
	    for (j=0; j<nr_class; j++) {
		/* transcribe probability estimates */
		if (ls != NULL) {
		    /* re-order the columns */
		    gretl_matrix_set(P, i, j, pi[ls[j].pos]);
		} else {
		    gretl_matrix_set(P, i, j, pi[j]);
		}
	    }
	} else {
	    yhi = svm_predict(model, x);
	}
	yi = prob->y[i];
	if (!regression) {
	    n_correct += (yhi == yi);
	}
	if (w->flags & W_YSCALE) {
	    yhi = unscale_y(yhi, w);
	    yi = unscale_y(yi, w);
	}
	yhat[dset->t1 + i] = yhi;
	if (regression) {
	    /* deviation from mean */
	    dev = yi - ymean;
	    TSS += dev * dev;
	    /* actual minus predicted */
	    dev = yi - yhi;
	    SSR += dev * dev;
	    if (w->regcrit == REG_ROUND_MISS) {
		misses += yi != round(yhi);
	    } else if (w->regcrit == REG_ROUND_MAD) {
		MAD += fabs(yi - round(yhi));
	    } else {
		MAD += fabs(dev);
	    }
	}
    }

    if (pi != NULL) {
	free(pi);
	free(ls);
    } else if (get_sigma) {
	/* retrieve estimate of Laplace scale */
	w->svr_sigma = svm_get_svr_probability(model);
    }

    datastr = training ? N_("Training data") : N_("Test data");

    if (regression) {
	MAD /= prob->l;
	if (w->regcrit == REG_ROUND_MISS) {
	    pprintf(prn, _("%s: miss ratio = %g (MSE = %g, R^2 = %g)\n"), _(datastr),
		    misses / (double) prob->l, SSR / prob->l, 1.0 - SSR / TSS);
	} else if (w->regcrit > REG_MSE) {
	    pprintf(prn, "%s: MAD = %g (MSE = %g, R^2 = %g)\n", _(datastr),
		    MAD, SSR / prob->l, 1.0 - SSR / TSS);
	} else {
	    pprintf(prn, "%s: MSE = %g, R^2 = %g (MAD = %g)\n", _(datastr),
		    SSR / prob->l, 1.0 - SSR / TSS, MAD);
	}
    } else {
	pprintf(prn, "%s: %s = %d (%.1f %s)\n", _(datastr),
		_("correct predictions"), n_correct, 100 * n_correct / (double) prob->l, 
		_("percent"));
    }

    return 0;
}

static void print_xvalid_iter (sv_parm *parm,
			       sv_wrapper *w,
			       double val,
			       const char *label,
			       int iter,
			       PRN *prn)
{
    if (iter >= 0) {
	pprintf(prn, "[%d] ", iter + 1);
    } else {
	pprintf(prn, "\n%s:\n ", _("Cross validation"));
    }
    pprintf(prn, "C = %g", parm->C);
    if (uses_gamma(parm)) {
	pprintf(prn, ", gamma = %g", parm->gamma);
    }
    if (uses_epsilon(parm)) {
	pprintf(prn, ", epsilon = %g", parm->p);
    } else if (uses_nu(parm)) {
	pprintf(prn, ", nu = %g", parm->nu);
    }
    pprintf(prn, ": %s = %#.8g\n", label, val);
    gretl_flush(prn);
}

static int *get_fold_sizes (const sv_data *data,
			    sv_wrapper *w,
			    const DATASET *dset)
{
    int *ret = gretl_list_new(w->nfold);
    int i, t;

    if (w->flags & W_CONSEC) {
	/* use equally sized consecutive blocks */
	int t2 = w->t2_train > 0 ? w->t2_train : dset->t2;
	int n = t2 - dset->t1 + 1;
	int ni = n / w->nfold;

	for (i=1; i<=w->nfold; i++) {
	    ret[i] = ni;
	}
	ret[w->nfold] += n % w->nfold;
    } else {
	/* use the supplied folds series */
	for (i=1; i<=w->nfold; i++) {
	    ret[i] = 0;
	    for (t=0; t<data->l; t++) {
		if (w->flist[t+1] == i) {
		    ret[i] += 1;
		}
	    }
	}
    }

    return ret;
}

/* Carry out cross validation in the case where the user has provided
   a series to specify the "folds", or has specified a given number of
   consecutive blocks, as opposed to the default random subsetting.
*/

static void custom_xvalidate (const sv_data *prob,
			      const sv_parm *parm,
			      const sv_wrapper *w,
			      double *targ)
{
    int i, vi, ni;

    for (i=0; i<w->nfold; i++) {
	struct svm_problem subprob;
	struct svm_model *submodel;
	int jmin = 0, jmax = 0;
	int j, k, useobs;

	vi = i + 1;
	ni = w->fsize[vi];
	subprob.l = prob->l - ni;
	subprob.x = malloc(subprob.l * sizeof *subprob.x);
	subprob.y = malloc(subprob.l * sizeof *subprob.y);

	if (w->flags & W_CONSEC) {
	    /* find start and end points for fold */
	    jmin = i * w->fsize[1];
	    jmax = jmin + ni;
	}

	/* set the training subsample, excluding fold i */
	k = 0;
	for (j=0; j<prob->l; j++) {
	    if (w->flags & W_CONSEC) {
		useobs = j < jmin || j >= jmax;
	    } else {
		useobs = w->flist[j+1] != vi;
	    }
	    if (useobs) {
		subprob.x[k] = prob->x[j];
		subprob.y[k] = prob->y[j];
		k++;
	    }
	}

	/* train on the given subsample */
	submodel = svm_train(&subprob, parm);

	/* predict on the complementary subsample (fold i only) */
	if (w->flags & W_CONSEC) {
	    /* we don't have to scan the whole prob->x array */
	    for (j=jmin; j<jmax; j++) {
		targ[j] = svm_predict(submodel, prob->x[j]);
	    }
	} else {
	    /* the values we want may be interspersed */
	    for (j=0; j<prob->l; j++) {
		if (w->flist[j+1] == vi) {
		    targ[j] = svm_predict(submodel, prob->x[j]);
		}
	    }
	}

	svm_free_and_destroy_model(&submodel);
	free(subprob.x);
	free(subprob.y);
    }
}

/* implement a single cross validation pass */

static int xvalidate_once (sv_data *prob,
			   sv_parm *parm,
			   sv_wrapper *w,
			   double *targ,
			   double *crit,
			   int iter,
			   PRN *prn)
{
    int i, n = prob->l;

    if (w->fsize != NULL) {
	custom_xvalidate(prob, parm, w, targ);
    } else {
	maybe_set_svm_seed(w);
	svm_cross_validation(prob, parm, w->nfold, targ);
    }

    if (doing_regression(parm)) {
	double yi, yhi, dev, minimand = 0;

	for (i=0; i<prob->l; i++) {
	    yi = prob->y[i];
	    yhi = targ[i];
	    if (w->flags & W_YSCALE) {
		yi = unscale_y(yi, w);
		yhi = unscale_y(yhi, w);
	    }
	    dev = yi - yhi;
	    if (w->regcrit == REG_ROUND_MISS) {
		minimand += (yi != round(yhi));
	    } else if (w->regcrit == REG_ROUND_MAD) {
		minimand += fabs(yi - round(yhi));
	    } else if (w->regcrit == REG_MAD) {
		minimand += fabs(dev);
	    } else {
		minimand += dev * dev;
	    }
	}
	minimand /= n;
	if (prn != NULL) {
	    const char *s = (w->regcrit == REG_MSE)? "MSE" :
		(w->regcrit == REG_ROUND_MISS)? "miss ratio" : "MAD";

	    print_xvalid_iter(parm, w, minimand, s, iter, prn);
	}
	*crit = -minimand;
    } else {
	/* classification */
	double pc_correct = 0;
	int n_correct = 0;

	for (i=0; i<n; i++) {
	    if (targ[i] == prob->y[i]) {
		n_correct++;
	    }
	}
	pc_correct = 100.0 * n_correct / (double) n;
	if (prn != NULL) {
	    print_xvalid_iter(parm, w, pc_correct,
			      "percent correct", iter, prn);
	}
	*crit = pc_correct;
    }

    return 0;
}

/* For comparison, from libsvm's grid.py:

  -log2c {begin,end,step | "null"} : set the range of c (default -5,15,2)
    begin,end,step -- c_range = 2^{begin,...,begin+k*step,...,end}
    "null"         -- do not grid with c
  -log2g {begin,end,step | "null"} : set the range of g (default 3,-15,-2)
    begin,end,step -- g_range = 2^{begin,...,begin+k*step,...,end}
    "null"         -- do not with g
*/

static int grid_set_dimensions (sv_grid *g, const gretl_matrix *m)
{
    int i, linvals = 0;
    double x;

    if (m != NULL && m->cols == 4) {
	linvals = 1;
    }

    for (i=0; i<=G_p; i++) {
	if ((g->row[i].stop < g->row[i].start && g->row[i].step >= 0) ||
	    (g->row[i].stop > g->row[i].start && g->row[i].step <= 0)) {
	    return E_INVARG;
	}
	g->null[i] = g->n[i] = g->linear[i] = 0;
	if (g->row[i].start == 0 && g->row[i].stop == 0 && g->row[i].step == 0) {
	    /* flag clamping of this row */
	    g->null[i] = g->n[i] = 1;
	} else if (g->row[i].stop >= g->row[i].start) {
	    for (x=g->row[i].start; x<=g->row[i].stop; x+=g->row[i].step) {
		g->n[i] += 1;
	    }
	} else {
	    for (x=g->row[i].start; x>=g->row[i].stop; x+=g->row[i].step) {
		g->n[i] += 1;
	    }
	}
	if (linvals && i < m->rows && gretl_matrix_get(m, i, 3) != 0) {
	    g->linear[i] = 1;
	}
    }

    return 0;
}

static void sv_grid_default (sv_grid *g, sv_parm *parm)
{
    int gsearch = parm == NULL || uses_gamma(parm);
    int altgrid = 0;

    if (altgrid) {
	g->row[G_C].start = -2;
	g->row[G_C].stop = 6;
	g->row[G_C].step = 1;
    } else {
	g->row[G_C].start = -5;
	g->row[G_C].stop = 8; /* grid.py has 15 (too big?) */
	g->row[G_C].step = 2;
    }

    if (gsearch) {
	g->row[G_g].start = 3;
	g->row[G_g].stop = -15;
	g->row[G_g].step = -2;
    } else {
	g->row[G_g].start = 0;
	g->row[G_g].stop = 0;
	g->row[G_g].step = 0;
    }

    g->row[G_p].start = 0;
    g->row[G_p].stop = 0;
    g->row[G_p].step = 0;

    grid_set_dimensions(g, NULL);
}

static void print_grid (sv_grid *g, sv_parm *parm, PRN *prn)
{
    const char *labels[] = {
	"C",
	"gamma",
	""
    };
    int i, imax = 2;

    if (!g->null[G_p]) {
	imax = 3;
	labels[2] = uses_nu(parm) ? "nu" : "eps";
    }

    pputs(prn, _("parameter search grid (start, stop, step):\n"));

    for (i=0; i<imax; i++) {
	if (!g->null[i]) {
	    pprintf(prn, " %-8s %g, %g, %g", labels[i],
		    g->row[i].start, g->row[i].stop,
		    g->row[i].step);
	    if (g->n[i] > 1) {
		pprintf(prn, " (%d %s, %s)\n", g->n[i], _("values"),
			g->linear[i] ? _("linear") : _("log2-based"));
	    } else {
		pputc(prn, '\n');
	    }
	}
    }
    pputc(prn, '\n');
}

static double grid_get_C (sv_grid *g, int i)
{
    double Ci = g->row[G_C].start + i * g->row[G_C].step;

    return g->linear[G_C] ? Ci : pow(2.0, Ci);
}

static double grid_get_g (sv_grid *g, int j)
{
    double gj = g->row[G_g].start + j * g->row[G_g].step;

    return g->linear[G_g] ? gj : pow(2.0, gj);
}

static double grid_get_p (sv_grid *g, int k)
{
    double pk = g->row[G_p].start + k * g->row[G_p].step;

    return g->linear[G_p] ? pk : pow(2.0, pk);
}

static int sv_wrapper_add_grid (sv_wrapper *w,
				sv_parm *parm,
				const gretl_matrix *m)
{
    sv_grid *g;
    int i, err = 0;

    g = calloc(1, sizeof *g);

    if (g == NULL) {
	err = E_ALLOC;
    } else if (m != NULL) {
	if (m->rows < 1 || m->cols < 3) {
	    err = E_INVARG;
	} else {
	    for (i=0; i<m->rows; i++) {
		g->row[i].start = gretl_matrix_get(m, i, 0);
		g->row[i].stop  = gretl_matrix_get(m, i, 1);
		g->row[i].step  = gretl_matrix_get(m, i, 2);
	    }
	    err = grid_set_dimensions(g, m);
	}
    } else {
	sv_grid_default(g, parm);
    }

    if (err) {
	if (g != NULL) {
	    free(g);
	}
    } else {
	w->grid = g;
    }

    return err;
}

static void sv_wrapper_remove_grid (sv_wrapper *w)
{
    if (w != NULL && w->grid != NULL) {
	free(w->grid);
	w->grid = NULL;
    }
}

static int write_plot_file (sv_wrapper *w,
			    sv_parm *parm,
			    double cmax)
{
    gretl_matrix *m = w->xdata;
    const char *zlabel = "MSE";
    double xprev, x[3], best[3] = {0};
    int critcol = m->cols - 1;
    int display, iact = 0;
    int i, j, err = 0;
    FILE *fp;

    set_optval_string(GNUPLOT, OPT_U, w->plot);

    if (!strcmp(w->plot, "display")) {
	display = iact = 1;
    }

    if (gui_mode) {
	fp = open_3d_plot_input_file(&iact);
    } else {
	fp = open_plot_input_file(PLOT_USER, 0, &err);
    }
    if (fp == NULL) {
	err = E_FOPEN;
	return err;
    }

    if (!doing_regression(parm)) {
	/* must be classification */
	zlabel = "correct";
    }

    gretl_push_c_numeric_locale();

    if (w->grid->linear[G_C]) {
	fputs("set xlabel 'C'\n", fp);
    } else {
	fputs("set xlabel 'log2(C)'\n", fp);
    }
    if (w->grid->linear[G_g]) {
	fputs("set ylabel 'gamma'\n", fp);
    } else {
	fputs("set ylabel 'log2(gamma)'\n", fp);
    }
    fprintf(fp, "set zlabel '%s'\n", zlabel);
    fputs("splot '-' using 1:2:3 title '' w pm3d ,\\\n", fp);
    fputs(" '-' using 1:2:3 title 'best' w p lt 1 pt 8\n", fp);
    xprev = 0;
    for (i=0; i<m->rows; i++) {
	x[0] = gretl_matrix_get(m, i, 0);
	x[1] = gretl_matrix_get(m, i, 1);
	x[2] = gretl_matrix_get(m, i, critcol);
	if (x[2] == cmax) {
	    for (j=0; j<3; j++) {
		best[j] = x[j];
	    }
	}
	if (i > 0 && x[0] != xprev) {
	    fputc('\n', fp);
	}
	xprev = x[0];
	fprintf(fp, "%g %g %g\n",
		w->grid->linear[0] ? x[0] : log2(x[0]),
		w->grid->linear[1] ? x[1] : log2(x[1]),
		x[2]);
    }
    fputs("e\n", fp);
    fprintf(fp, "%g %g %g\n",
	    w->grid->linear[0] ? best[0] : log2(best[0]),
	    w->grid->linear[1] ? best[1] : log2(best[1]),
	    best[2]);
    fputs("e\n", fp);

    gretl_pop_c_numeric_locale();

    if (gui_mode && display && iact) {
	fputs("pause mouse close\n", fp);
	err = finalize_3d_plot_input_file(fp);
	if (!err && gui_mode) {
	    manufacture_gui_callback(GP_ASYNC);
	}
    } else {
	err = finalize_plot_input_file(fp);
	if (!err && gui_mode) {
	    /* wanted? */
	    manufacture_gui_callback(GNUPLOT);
	}
    }

    return err;
}

/* get ready to do parameter search */

static int do_search_prep (sv_data *data,
			   sv_parm *parm,
			   sv_wrapper *w,
			   const DATASET *dset,
			   PRN *prn)
{
    int err = 0;

    if (0 && parm->kernel_type != RBF) {
	pputs(prn, _("Non-RBF kernel, don't know how to tune parameters\n"));
	sv_wrapper_remove_grid(w);
	err = E_INVARG;
    } else if (w->grid == NULL) {
	sv_wrapper_add_grid(w, parm, NULL);
    }

    if (!err && (w->flags & (W_FOLDVAR | W_CONSEC))) {
	w->fsize = get_fold_sizes(data, w, dset);
    }

    return err;
}

static void maybe_hush (sv_wrapper *w)
{
    if (!(w->flags & W_QUIET)) {
	/* cut out excessive verbosity */
	svm_set_print_string_function(gretl_libsvm_noprint);
    }
}

static void maybe_resume_printing (sv_wrapper *w)
{
    if (!(w->flags & W_QUIET)) {
	svm_set_print_string_function(gretl_libsvm_print);
    }
}

static int can_write_plot (sv_wrapper *w)
{
    /* for now we handle only the case of a 2D grid
       with C and gamma, but this should be generalized
       at some point
    */
    if (w->xdata == NULL) {
	return 0;
    } else if (w->grid->null[G_C] || w->grid->null[G_g]) {
	return 0;
    } else if (!w->grid->null[G_p]) {
	return 0;
    } else {
	return 1;
    }
}

static void param_search_finalize (sv_parm *parm,
				   sv_wrapper *w,
				   double cmax,
				   PRN *prn)
{
    sv_grid *grid = w->grid;

    if (w->plot != NULL && can_write_plot(w)) {
	write_plot_file(w, parm, fabs(cmax));
    }

    pprintf(prn, _("*** Criterion optimized at %g: C=%g"), fabs(cmax), parm->C);
    if (uses_gamma(parm)) {
	pprintf(prn, ", gamma=%g", parm->gamma);
    }
    if (grid->null[G_p]) {
	; /* neither epsilon nor nu used in x-validation */
    } else if (uses_epsilon(parm)) {
	pprintf(prn, ", epsilon=%g", parm->p);
    } else if (uses_nu(parm)) {
	pprintf(prn, ", nu=%g", parm->nu);
    }
    pputs(prn, " ***\n");
}

#ifdef HAVE_MPI

static int cross_validate_worker_task (sv_data *data,
				       sv_parm *parm,
				       sv_wrapper *w,
				       double *targ)
{
    gretl_matrix *task = NULL;
    double crit;
    int i, j, seq;
    int err = 0;

    /* get our sub-task matrix */
    task = gretl_matrix_mpi_receive(0, &err);

    /* do the actual cross validation */
    for (i=0; i<task->rows && !err; i++) {
	j = 0;
	parm->C     = gretl_matrix_get(task, i, j++);
	parm->gamma = gretl_matrix_get(task, i, j++);
	if (parm->svm_type == EPSILON_SVR) {
	    parm->p = gretl_matrix_get(task, i, j++);
	} else if (parm->svm_type == NU_SVR) {
	    parm->nu = gretl_matrix_get(task, i, j++);
	}
	seq = gretl_matrix_get(task, i, task->cols - 1);
	err = xvalidate_once(data, parm, w, targ, &crit, seq, worker_prn);
	if (!err) {
	    gretl_matrix_set(task, i, j, crit);
	}
    }

    /* send results back to root */
    if (!err) {
	err = gretl_matrix_mpi_send(task, 0);
    }

    /* local clean up */
    gretl_matrix_free(task);

    return err;
}

static gretl_matrix **allocate_submatrices (int n,
					    int cols,
					    int r1,
					    int nr1)
{
    gretl_matrix **M = malloc(n * sizeof *M);
    int i, err = 0;

    if (M == NULL) {
	err = E_ALLOC;
    } else {
	for (i=0; i<n && !err; i++) {
	    if (i < nr1) {
		M[i] = gretl_matrix_alloc(r1, cols);
	    } else {
		M[i] = gretl_matrix_alloc(r1 + 1, cols);
	    }
	    if (M[i] == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    if (err && M != NULL) {
	free(M);
	M = NULL;
    }

    return M;
}

static int process_results (sv_parm *parm,
			    sv_wrapper *w,
			    gretl_matrix *m,
			    double *pcmax)
{
    double crit, cmax = -DBL_MAX;
    int i, imax = 0;

    for (i=0; i<m->rows; i++) {
	crit = gretl_matrix_get(m, i, m->cols - 1);
	if (crit > cmax) {
	    cmax = crit;
	    imax = i;
	}
	gretl_matrix_set(m, i, m->cols - 1, fabs(crit));
    }

    parm->C = gretl_matrix_get(m, imax, 0);
    parm->gamma = gretl_matrix_get(m, imax, 1);
    if (!w->grid->null[G_p]) {
	if (parm->svm_type == EPSILON_SVR) {
	    parm->p = gretl_matrix_get(m, imax, 2);
	} else {
	    parm->nu = gretl_matrix_get(m, imax, 2);
	}
    }

    /* wrapper takes ownership of @m */
    w->xdata = m;

    *pcmax = cmax;

    return 0;
}

/* The following is called only by the root process,
   when doing cross validation via MPI */

static int carve_up_xvalidation (sv_data *data,
				 sv_parm *parm,
				 sv_wrapper *w,
				 double *targ,
				 double *pcmax,
				 PRN *prn)
{
    sv_grid *grid = w->grid;
    gretl_matrix **M = NULL;
    gretl_matrix *m = NULL;
    int *row = NULL;
    int nC = grid->n[G_C];
    int ng = grid->n[G_g];
    int np = grid->n[G_p];
    double C, g, p, crit;
    double *p3 = NULL;
    int ncom, nproc = w->nproc;
    int i, j, k, ii;
    int r1, nr1, seq;
    int rem, ncols = 4;
    int err = 0;

    if (uses_epsilon(parm)) {
	p = parm->p;
	ncols++;
	if (!grid->null[G_p]) {
	    p3 = &parm->p;
	}
    } else if (uses_nu(parm)) {
	p = parm->nu;
	ncols++;
	if (!grid->null[G_p]) {
	    p3 = &parm->nu;
	}
    } else {
	p = 0;
    }

    ncom = nC * ng * np;  /* number of param combinations */
    r1 = ncom / nproc;    /* initial rows per matrix */
    rem = ncom % nproc;   /* remainder, if any */
    nr1 = nproc - rem;    /* number of matrices with r1 rows */

    if (prn != NULL) {
	pprintf(prn, _("MPI: dividing %d parameter specifications "
		"between %d processes\n"), ncom, nproc);
    }

    /* matrices to store per-worker parameter sets */
    M = allocate_submatrices(nproc, ncols, r1, nr1);

    if (M == NULL) {
	err = E_ALLOC;
    } else {
	/* array to keep track of which row we're on for
	   each of the submatrices */
	row = malloc(nproc * sizeof *row);
	if (row == NULL) {
	    err = E_ALLOC;
	} else {
	    for (i=0; i<nproc; i++) {
		row[i] = 0;
	    }
	}
    }

    if (err) {
	goto bailout;
    }

    m = M[nproc-1]; /* convenience pointer to root's submatrix */
    C = parm->C;
    g = parm->gamma;

    /* We'll dole out the full parameter matrix row by row: this
       is an attempt to even up the tasks, given that doing the
       cross validation is much more expensive for large values
       of C. We add an extra column holding "seq" so that we can
       easily put the reassembled matrix into "standard order".
    */

    ii = seq = 0;
    for (i=0; i<nC; i++) {
	if (!grid->null[G_C]) {
	    C = grid_get_C(grid, i);
	}
	for (j=0; j<ng; j++) {
	    if (!grid->null[G_g]) {
		g = grid_get_g(grid, j);
	    }
	    for (k=0; k<np; k++) {
		if (!grid->null[G_p]) {
		    p = grid_get_p(grid, k);
		}
		if (row[ii] == M[ii]->rows) {
		    /* skip to the larger matrices */
		    ii = nr1;
		}
		gretl_matrix_set(M[ii], row[ii], 0, C);
		gretl_matrix_set(M[ii], row[ii], 1, g);
		if (ncols == 5) {
		    gretl_matrix_set(M[ii], row[ii], 2, p);
		}
		gretl_matrix_set(M[ii], row[ii], ncols - 2, 0);
		gretl_matrix_set(M[ii], row[ii], ncols - 1, seq);
		seq++;
		row[ii] += 1;
		if (ii == nproc - 1) {
		    ii = 0;
		} else {
		    ii++;
		}
	    }
	}
    }

    /* send matrices to workers, keeping the last for root */
    for (i=0; i<nproc-1; i++) {
	gretl_matrix_mpi_send(M[i], i+1);
	gretl_matrix_free(M[i]);
	M[i] = NULL;
    }

    maybe_hush(w);

    /* do root's share of the cross validation */
    for (i=0; i<m->rows && !err; i++) {
	parm->C = gretl_matrix_get(m, i, 0);
	parm->gamma = gretl_matrix_get(m, i, 1);
	if (!grid->null[G_p]) {
	    *p3 = gretl_matrix_get(m, i, 2);
	}
	seq = gretl_matrix_get(m, i, ncols-1);
	err = xvalidate_once(data, parm, w, targ, &crit, seq, prn);
	if (!err) {
	    gretl_matrix_set(m, i, ncols-2, crit);
	}
    }

    maybe_resume_printing(w);

    /* get results back from workers and process */
    if (!err) {
	gretl_matrix *Tmp, *Res = NULL;
	gretl_matrix *mi;
	int row = 0;

	Tmp = gretl_matrix_alloc(ncom, 5);
	if (Tmp == NULL) {
	    err = E_ALLOC;
	}
     	for (i=1; i<nproc && !err; i++) {
	    mi = gretl_matrix_mpi_receive(i, &err);
	    if (mi != NULL) {
		gretl_matrix_inscribe_matrix(Tmp, mi, row, 0, GRETL_MOD_NONE);
		row += mi->rows;
		gretl_matrix_free(mi);
	    }
	}
	if (!err) {
	    gretl_matrix_inscribe_matrix(Tmp, m, row, 0, GRETL_MOD_NONE);
	    Res = gretl_matrix_sort_by_column(Tmp, ncols - 1, &err);
	}
	gretl_matrix_free(Tmp);
	if (Res != NULL) {
	    /* hide the sequence column */
	    gretl_matrix_reuse(Res, -1, ncols - 1);
	    err = process_results(parm, w, Res, pcmax);
	}
    }

 bailout:

    free(M);
    free(row);

    return err;
}

static int call_mpi_cross_validation (sv_data *data,
				      sv_parm *parm,
				      sv_wrapper *w,
				      PRN *prn)
{
    double cmax = -DBL_MAX;
    double *targ;
    int err = 0;

    targ = malloc(data->l * sizeof *targ);
    if (targ == NULL) {
	return E_ALLOC;
    }

    /* arrange to get rand() in sync across the processes */
    gretl_mpi_bcast(&w->seed, GRETL_TYPE_UINT64, 0);

    if (w->rank == 0) {
	if (prn != NULL) {
	    print_grid(w->grid, parm, prn);
	}
	err = carve_up_xvalidation(data, parm, w, targ, &cmax, prn);
    } else {
	err = cross_validate_worker_task(data, parm, w, targ);
    }

    if (w->rank == 0) {
	param_search_finalize(parm, w, cmax, prn);
    }

    free(targ);

    return err;
}

#endif /* HAVE_MPI */

static int call_cross_validation (sv_data *data,
				  sv_parm *parm,
				  sv_wrapper *w,
				  double *targ,
				  PRN *prn)
{
    double crit;
    int err = 0;

    pprintf(prn, "%s ", _("Cross-validation"));
    if (w->flags & W_CONSEC) {
	pprintf(prn, _("using %d consecutive blocks"), w->nfold);
    } else if (w->flags & W_FOLDVAR) {
	pprintf(prn, _("using %d values of %s"), w->nfold, w->foldname);
    } else {
	pprintf(prn, _("using %d random folds"), w->nfold);
    }
    pprintf(prn, " %s\n", _("(may take a while)"));
    gretl_flush(prn);

    if (w->grid != NULL) {
	sv_grid *grid = w->grid;
	double cmax = -DBL_MAX;
	double *p3 = NULL;
	int nC = grid->n[G_C];
	int ng = grid->n[G_g];
	int np = grid->n[G_p];
	int ibest = 0, jbest = 0, kbest = 0;
	int iter = 0;
	int ncols = 3;
	int i, j, k;

	if (prn != NULL) {
	    print_grid(grid, parm, prn);
	}

	if (uses_epsilon(parm)) {
	    ncols++;
	    p3 = &parm->p;
	} else if (uses_nu(parm)) {
	    ncols++;
	    p3 = &parm->nu;
	}

	if ((w->plot != NULL || (w->flags & W_SVPARM)) && (nC > 1 || ng > 1)) {
	    w->xdata = gretl_matrix_alloc(nC * ng * np, ncols);
	}

	maybe_hush(w);

	for (i=0; i<nC; i++) {
	    if (!grid->null[G_C]) {
		parm->C = grid_get_C(grid, i);
	    }
	    for (j=0; j<ng; j++) {
		if (!grid->null[G_g]) {
		    parm->gamma = grid_get_g(grid, j);
		}
		for (k=0; k<np; k++) {
		    if (!grid->null[G_p]) {
			*p3 = grid_get_p(grid, k);
		    }
		    xvalidate_once(data, parm, w, targ, &crit, iter, prn);
		    if (crit > cmax) {
			cmax = crit;
			ibest = i;
			jbest = j;
			kbest = k;
		    }
		    if (w->xdata != NULL) {
			gretl_matrix_set(w->xdata, iter, 0, parm->C);
			gretl_matrix_set(w->xdata, iter, 1, parm->gamma);
			if (p3 != NULL) {
			    gretl_matrix_set(w->xdata, iter, 2, *p3);
			}
			gretl_matrix_set(w->xdata, iter, ncols - 1, fabs(crit));
		    }
		    iter++;
		}
	    }
	}

	maybe_resume_printing(w);

	if (!grid->null[G_C]) {
	    parm->C = grid_get_C(grid, ibest);
	}
	if (!grid->null[G_g]) {
	    parm->gamma = grid_get_g(grid, jbest);
	}
	if (!grid->null[G_p]) {
	    *p3 = grid_get_p(grid, kbest);
	}
	param_search_finalize(parm, w, cmax, prn);
    } else {
	/* no search: just a single cross validation pass
	   FIXME save the criterion to optional bundle arg?
	*/
	err = xvalidate_once(data, parm, w, targ, &crit, -1, prn);
    }

    return err;
}

/* If the caller provided a "folds" series, check it for
   validity: the values must be consecutive integers
   starting at 1, and the number of folds must be within
   bounds.
*/

static int check_folds_series (const int *list,
			       const DATASET *dset,
			       sv_wrapper *w,
			       PRN *prn)
{
    gretl_matrix *fvec;
    const double *x;
    int v = list[list[0]];
    int i, t2, n;
    int err = 0;

    x = dset->Z[v] + dset->t1;
    t2 = w->t2_train > 0 ? w->t2_train : dset->t2;
    n = t2 - dset->t1 + 1;

    fvec = gretl_matrix_values(x, n, OPT_S, &err);
    if (!err) {
	n = gretl_vector_get_length(fvec);
	x = fvec->val;

	for (i=0; i<n && !err; i++) {
	    if (i > 0) {
		if (x[i] != x[i-1] + 1) {
		    fputs("foldvar must hold consecutive integers\n", stderr);
		    err = E_DATA;
		}
	    } else if (x[i] != 1.0) {
		fputs("foldvar values must start at 1\n", stderr);
		err = E_DATA;
	    }
	}

	if (!err) {
	    int nf = x[n-1];

	    pprintf(prn, _("%s: found %d folds\n"), dset->varname[v], nf);
	    if (nf < 2 || (w->nfold > 0 && nf != w->nfold)) {
		fprintf(stderr, "invalid number of folds %d\n", nf);
		err = E_DATA;
	    } else {
		strcpy(w->foldname, dset->varname[v]);
		w->nfold = nf;
	    }
	}
    }

    gretl_matrix_free(fvec);

    return err;
}

/* It's OK if there's no @key value in bundle @b, but if it
   is present it should hold an integer value. Return 1
   if key found and integer retrieved OK, else 0.
*/

static int get_optional_int (gretl_bundle *b, const char *key,
			     int *ival, int *err)
{
    GretlType type = 0;
    void *ptr;

    ptr = gretl_bundle_get_data(b, key, &type, NULL, NULL);

    if (ptr != NULL) {
	if (type == GRETL_TYPE_INT) {
	    *ival = *(int *) ptr;
	    return 1;
	} else if (type == GRETL_TYPE_DOUBLE) {
	    *ival = *(double *) ptr;
	    return 1;
	} else if (type == GRETL_TYPE_UINT32) {
	    *ival = *(guint32 *) ptr;
	    return 1;
	} else if (err != NULL) {
	    *err = E_TYPES;
	}
    }

    return 0;
}

static int get_optional_unsigned (gretl_bundle *b, const char *key,
				  unsigned int *uval, int *err)
{
    GretlType type = 0;
    void *ptr;

    ptr = gretl_bundle_get_data(b, key, &type, NULL, NULL);

    if (ptr != NULL) {
	if (type == GRETL_TYPE_INT) {
	    *uval = *(int *) ptr;
	    return 1;
	} else if (type == GRETL_TYPE_DOUBLE) {
	    *uval = *(double *) ptr;
	    return 1;
	} else if (type == GRETL_TYPE_UINT32) {
	    *uval = *(unsigned *) ptr;
	    return 1;
	} else {
	    *err = E_TYPES;
	}
    }

    return 0;
}

/* Determine if @s is a recognized parameter key: we
   do this so we can flag anything that may be a
   mistyped key.
*/

static int is_w_parm (const char *s)
{
    const char *wparms[] = {
	"loadmod", "scaling", "predict", "n_train",
	"folds", "seed", "quiet", "search",
	"foldvar", "consecutive", "yscale",
	"search_only", "grid", "ranges_outfile",
	"data_outfile", "ranges_infile", "model_outfile",
	"model_infile", "plot", "range_format", "refold",
	"autoseed", "regcrit", "use_mpi", NULL
    };
    int i;

    for (i=0; wparms[i] != NULL; i++) {
	if (!strcmp(s, wparms[i])) {
	    return 1;
	}
    }

    return 0;
}

static int check_user_params (gretl_array *A)
{
    char **pstrs;
    int i, np, err = 0;

    pstrs = gretl_array_get_strings(A, &np);
    if (pstrs == NULL) {
	return 0; /* ? */
    }

    for (i=0; i<np; i++) {
	if (!is_w_parm(pstrs[i]) && !is_sv_parm(pstrs[i])) {
	    gretl_errmsg_sprintf(_("svm: unrecognized parameter '%s'"),
				 pstrs[i]);
	    err = E_BADOPT;
	    break;
	}
    }

    return err;
}

/* process @bparm, which contains the parameters provided by
   the caller */

static int read_params_bundle (gretl_bundle *bparm,
			       gretl_bundle *bmod,
			       sv_wrapper *wrap,
			       sv_parm *parm,
			       const int *list,
			       const DATASET *dset,
			       PRN *prn)
{
    gretl_array *A;
    int no_savemod = 0;
    unsigned int uval;
    int ival, err = 0;

    /* got any bad keys in @bparm? */
    A = gretl_bundle_get_keys(bparm, NULL);
    if (A != NULL) {
	err = check_user_params(A);
	gretl_array_destroy(A);
	if (err) {
	    return err;
	}
    }

    /* start by reading some info that's not included in
       the libsvm @parm struct
    */

    if (get_optional_int(bparm, "loadmod", &ival, &err)) {
	if (ival != 0 && bmod == NULL) {
	    gretl_errmsg_sprintf(_("svm: invalid 'loadmod' value %d\n"), ival);
	    err = E_INVARG;
	} else if (ival != 0) {
	    wrap->flags |= W_LOADMOD;
	    no_savemod = 1;
	}
    }

    if (get_optional_int(bparm, "scaling", &ival, &err)) {
	if (ival < 0 || ival > 2) {
	    gretl_errmsg_sprintf(_("svm: invalid 'scaling' value %d\n"), ival);
	    err = E_INVARG;
	} else {
	    wrap->scaling = ival;
	}
    }

    if (get_optional_int(bparm, "predict", &ival, &err)) {
	if (ival < 0 || ival > 2) {
	    gretl_errmsg_sprintf(_("svm: invalid 'predict' value %d\n"), ival);
	    err = E_INVARG;
	} else {
	    wrap->predict = ival;
	}
    }

    if (get_optional_int(bparm, "quiet", &ival, &err) && ival != 0) {
	wrap->flags |= W_QUIET;
    }

    if (get_optional_int(bparm, "n_train", &ival, &err)) {
	/* number of training observations: this sets a range
	   starting at the first complete observation in the
	   incoming sample, which was recorded as wrap->t1
	   after trimming any missing values
	*/
	if (ival == 0) {
	    wrap->flags |= W_NOTRAIN;
	    pprintf(prn, "n_train = 0, %s\n", _("no training will be done"));
	} else {
	    int nmax = wrap->t2 - wrap->t1 + 1;

	    if (ival < list[0]) {
		gretl_errmsg_sprintf(_("svm: n_train must be at least %d"), list[0]);
		err = E_INVARG;
	    } else if (ival > nmax) {
		gretl_errmsg_sprintf(_("svm: n_train cannot exceed the number of "
				     "complete observations, %d"), nmax);
		err = E_INVARG;
	    } else {
		char obs1[OBSLEN], obs2[OBSLEN];

		wrap->t2_train = wrap->t1 + ival - 1;
		get_obs_string(obs1, wrap->t1, dset);
		get_obs_string(obs2, wrap->t2_train, dset);
                if (!(wrap->flags & W_QUIET)) {
                    pprintf(prn, _("n_train = %d; use observations %s to %s for training\n"),
                            ival, obs1, obs2);
                }
	    }
	}
    }

    if (get_optional_int(bparm, "folds", &ival, &err)) {
	if (ival < 2) {
	    gretl_errmsg_sprintf(_("svm: invalid 'folds' value %d\n"), ival);
	    err = E_INVARG;
	} else {
	    wrap->nfold = ival;
	}
    }

    if (get_optional_unsigned(bparm, "seed", &uval, &err)) {
	wrap->seed = uval;
    }

    if (get_optional_int(bparm, "refold", &ival, &err) && ival != 0) {
	wrap->flags |= W_REFOLD;
    }

    if (get_optional_int(bparm, "search", &ival, &err) && ival != 0) {
	wrap->flags |= W_SEARCH;
    }

    if (get_optional_int(bparm, "use_mpi", &ival, &err) && ival != 0) {
	wrap->flags |= W_MPI;
    }

    if (get_optional_int(bparm, "foldvar", &ival, &err) && ival != 0) {
	wrap->flags |= W_FOLDVAR;
    }

    if (get_optional_int(bparm, "consecutive", &ival, &err) && ival != 0) {
	wrap->flags |= W_CONSEC;
    }

    if (get_optional_int(bparm, "yscale", &ival, &err) && ival != 0) {
	wrap->flags |= W_YSCALE;
    }

    if (get_optional_int(bparm, "search_only", &ival, &err) && ival != 0) {
	wrap->flags |= W_SEARCH;
	if (bmod != NULL) {
	    wrap->flags |= W_SVPARM;
	    no_savemod = 1;
	}
    }

    if (get_optional_int(bparm, "regcrit", &ival, &err) && ival != 0) {
	if (ival < REG_MSE || ival > REG_ROUND_MISS) {
	    gretl_errmsg_sprintf(_("svm: invalid regcrit value %d"), ival);
	    err = E_INVARG;
	} else {
	    wrap->regcrit = ival;
	}
    }

    if (!err) {
	const gretl_matrix *m;

	m = gretl_bundle_get_matrix(bparm, "grid", NULL);
	if (m != NULL) {
	    err = sv_wrapper_add_grid(wrap, NULL, m);
	    if (!err) {
		wrap->flags |= W_SEARCH;
	    }
	}
    }

    if (!err) {
	const char *strval;

	strval = gretl_bundle_get_string(bparm, "ranges_outfile", NULL);
	if (strval != NULL && *strval != '\0') {
	    wrap->ranges_outfile = gretl_strdup(strval);
	}
	strval = gretl_bundle_get_string(bparm, "data_outfile", NULL);
	if (strval != NULL && *strval != '\0') {
	    wrap->data_outfile = gretl_strdup(strval);
	}
	strval = gretl_bundle_get_string(bparm, "ranges_infile", NULL);
	if (strval != NULL && *strval != '\0') {
	    wrap->ranges_infile = gretl_strdup(strval);
	}
	strval = gretl_bundle_get_string(bparm, "model_outfile", NULL);
	if (strval != NULL && *strval != '\0') {
	    wrap->model_outfile = gretl_strdup(strval);
	}
	strval = gretl_bundle_get_string(bparm, "model_infile", NULL);
	if (strval != NULL && *strval != '\0') {
	    wrap->model_infile = gretl_strdup(strval);
	}
	strval = gretl_bundle_get_string(bparm, "plot", NULL);
	if (strval != NULL && *strval != '\0') {
	    wrap->plot = gretl_strdup(strval);
	}
	strval = gretl_bundle_get_string(bparm, "range_format", NULL);
	if (strval != NULL && !strcmp(strval, "libsvm")) {
	    wrap->flags |= W_STDFMT;
	}
    }

    if (!err) {
	if (bmod != NULL && !no_savemod) {
	    /* implicitly, a model should be saved to @bmod */
	    wrap->flags |= W_SAVEMOD;
	}
	if ((wrap->flags & W_SEARCH) && !(wrap->flags & W_FOLDVAR) &&
	    wrap->nfold == 0) {
	    /* default to 5 folds */
	    wrap->nfold = 5;
	}
    }

    if (!err && (wrap->flags & W_FOLDVAR)) {
	if (wrap->flags & W_CONSEC) {
	    gretl_errmsg_set(_("The foldvar and consecutive options "
			     "are incompatible"));
	    err = E_BADOPT;
	} else {
	    err = check_folds_series(list, dset, wrap, prn);
	}
    }

    if (!err) {
	/* if we're still OK, fill out the libsvm @parm struct */
	err = set_svm_parm(parm, bparm, prn);
    }

#if 0 /* FIXME? */
    if (!err) {
	/* param search: at present we only support the RBF kernel */
	if ((wrap->flags & W_SEARCH) && parm->kernel_type != RBF) {
	    gretl_errmsg_set("parameter search is only supported for "
			     "the RBF kernel at present");
	    err = E_INVARG;
	}
    }
#endif

    return err;
}

/* check in advance if the incoming list variable is supposed to
   hold a fold variable in last place */

static int peek_foldvar (gretl_bundle *bparm, int *fvar)
{
    int ival, err = 0;

    if (get_optional_int(bparm, "foldvar", &ival, &err) && ival != 0) {
	*fvar = 1;
    }

    return err;
}

#ifdef HAVE_MPI

/* Here we're trying to tell if we should invoke automated
   local-machine MPI -- having already ascertained that it's
   feasible in general terms.
*/

static int peek_use_mpi (gretl_bundle *bparm,
			 gretl_bundle *bmodel,
			 gretl_bundle *bprob,
			 PRN *prn)
{
    int ival = 0, ret = 0;

    if (bmodel != NULL || bprob != NULL) {
	/* these cases are not handled yet */
	return 0;
    }

    if (get_optional_int(bparm, "use_mpi", &ival, NULL)) {
	/* if the user specified @use_mpi, respect it provisionally */
	ret = ival;
    } else {
	/* otherwise default to using MPI subject to the check below? */
	ret = 1;
    }

    if (ret > 0) {
	/* MPI is actually used internally only in parameter search,
	   so let's peek at that.
	*/
	int search = 0;

	ival = 0;
	if (get_optional_int(bparm, "search", &ival, NULL) && ival > 0) {
	    search = 1;
	} else if (gretl_bundle_get_matrix(bparm, "grid", NULL) != NULL) {
	    search = 1;
	}
	if (!search && ival > 0) {
	    pputs(prn, _("svm: 'use_mpi' option ignored in the absence of search"));
	    pputc(prn, '\n');
	    ret = 0;
	}
    }

    return ret;
}

#endif

/* wrap a couple of libsvm I/0 functions so they respect
   @workdir
*/

static sv_model *svm_load_model_wrapper (const char *fname,
					 int *err)
{
    sv_model *model = NULL;

    fname = gretl_maybe_switch_dir(fname);
    model = svm_load_model(fname);

    if (model == NULL) {
	*err = E_FOPEN;
    }

    return model;
}

static int svm_save_model_wrapper (const char *fname,
				   const sv_model *model)
{
    int err = 0;

    fname = gretl_maybe_switch_dir(fname);
    err = svm_save_model(fname, model);

    return err ? E_FOPEN : 0;
}

static void report_result (int err, PRN *prn)
{
    if (prn != NULL) {
	if (err) {
	    pprintf(prn, "err = %d\n", err);
	} else {
	    pputs(prn, "OK\n");
	}
    }
}

static int get_svm_ranges (const int *list,
			   DATASET *dset,
			   sv_wrapper *wrap,
			   PRN *prn)
{
    int err = 0;

    if (wrap->ranges_infile != NULL) {
	pprintf(prn, _("Getting data ranges from %s... "),
		wrap->ranges_infile);
	err = read_ranges(wrap);
	report_result(err, prn);
	gretl_flush(prn);
    } else {
	pprintf(prn, _("Getting data ranges (sample = %d to %d)... "),
		dset->t1 + 1, dset->t2 + 1);
	err = get_data_ranges(list, dset, wrap);
	report_result(err, prn);
	if (!err && wrap->ranges_outfile != NULL) {
	    err = write_ranges(wrap);
	    report_result(err, prn);
	}
	gretl_flush(prn);
    }

    return err;
}

/* load an svm model, either from a bundle in memory
   or a text file in libsvm format */

static sv_model *do_load_model (sv_wrapper *w,
				gretl_bundle *b,
				PRN *prn,
				int *err)
{
    sv_model *model = NULL;

    if (w->flags & W_LOADMOD) {
	pputs(prn, _("Loading svm model from bundle... "));
	model = svm_model_from_bundle(b, err);
    } else {
	pprintf(prn, _("Loading svm model from %s... "), w->model_infile);
	model = svm_load_model_wrapper(w->model_infile, err);
    }

    return model;
}

/* save an svm model, either to a bundle in memory
   or a text file in libsvm format */

static int do_save_model (sv_model *model, sv_wrapper *w,
			  gretl_bundle *b, PRN *prn)
{
    int err = 0;

    if (w->flags & W_SAVEMOD) {
	pputs(prn, _("Saving svm model to bundle... "));
	err = svm_model_save_to_bundle(model, b);
    }
    if (w->model_outfile != NULL) {
	pprintf(prn, _("Saving svm model as %s... "), w->model_outfile);
	err = svm_save_model_wrapper(w->model_outfile, model);
    }

    return err;
}

/* here we call the libsvm parameter-checking
   function; if all is well and we're not in quiet
   mode, we print some information on the primary
   parameters
*/

static int check_svm_params (sv_data *data,
			     sv_parm *parm,
			     sv_wrapper *w,
			     PRN *prn)
{
    const char *msg = svm_check_parameter(data, parm);
    int prob_ok;
    int err = 0;

    prob_ok = parm->svm_type == C_SVC ||
	parm->svm_type == NU_SVC ||
	parm->svm_type == EPSILON_SVR ||
	parm->svm_type == NU_SVR;

    pputs(prn, _("Checking parameter values... "));
    if (msg != NULL) {
	pputs(prn, "problem\n");
	gretl_errmsg_sprintf("svm: %s", msg);
	err = E_INVARG;
    } else if (parm->probability && !prob_ok) {
	pputs(prn, "problem\n");
	gretl_errmsg_set(_("svm: probability estimates not supported "
			 "for this specification"));
	err = E_INVARG;
    } else if (!doing_regression(parm) && w->regcrit > 0) {
	pputs(prn, "problem\n");
	gretl_errmsg_sprintf(_("svm: cross validation criterion %d not applicable"),
			     w->regcrit);
	err = E_INVARG;
    } else if (w->regcrit > 2 && !(w->flags & W_INTDEP)) {
	pputs(prn, "problem\n");
	gretl_errmsg_sprintf(_("svm: cross validation criterion %d not applicable"),
			     w->regcrit);
	err = E_INVARG;
    } else if (prn != NULL) {
	pputs(prn, "OK\n");
	pprintf(prn, "svm_type %s, kernel_type %s\n",
		svm_type_names[parm->svm_type],
		kernel_type_names[parm->kernel_type]);
	pprintf(prn, _("initial params: C = %g"), parm->C);
	if (uses_gamma(parm)) {
	    pprintf(prn, ", gamma = %g", parm->gamma);
	}
	if (uses_epsilon(parm)) {
	    pprintf(prn, ", epsilon = %g", parm->p);
	}
	if (uses_nu(parm)) {
	    pprintf(prn, ", nu = %g", parm->nu);
	}
	pputs(prn, "\n\n");
    }

    if (!err && doing_regression(parm) && w->regcrit == 0) {
	w->regcrit = REG_MSE;
    }

    return err;
}

/* note: if we're in auto-MPI mode, only the rank 0 process
   executes this function
*/

static int svm_predict_main (const int *list,
			     gretl_bundle *bmodel,
			     double *yhat,
			     int *yhat_written,
			     DATASET *dset,
			     PRN *prn,
			     sv_parm *parm,
			     sv_wrapper *wrap,
			     sv_data *prob1)
{
    sv_data *prob2 = NULL;
    sv_cell *x_space2 = NULL;
    sv_model *model = NULL;
    int do_training = 1;
    int err = 0;

    gui_mode = gretl_in_gui_mode();

    if (!err && wrap->data_outfile != NULL && wrap->rank <= 0) {
	err = write_problem(prob1, wrap);
	report_result(err, prn);
    }

    if (wrap->flags & W_NOTRAIN) {
	do_training = 0;
    }

    if (!err && loading_model(wrap)) {
	do_training = 0;
	model = do_load_model(wrap, bmodel, prn, &err);
	report_result(err, prn);
    }

    if (!err && wrap->nfold > 0 && model == NULL) {
#ifdef HAVE_MPI
	if (wrap->nproc > 1) {
	    err = call_mpi_cross_validation(prob1, parm, wrap, prn);
	} else {
	    err = call_cross_validation(prob1, parm, wrap, yhat + dset->t1,
					prn);
	    if (!err) {
		*yhat_written = 1;
	    }
	}
#else
	err = call_cross_validation(prob1, parm, wrap, yhat + dset->t1, prn);
	if (!err) {
	    *yhat_written = 1;
	}
#endif
	if (wrap->flags & W_SEARCH) {
	    /* continue to train using tuned params? */
	    if (wrap->flags & W_SVPARM) {
		/* no, just save the tuned parms */
		save_results_to_bundle(parm, wrap, bmodel);
		wrap->predict = do_training = 0;
	    }
	} else {
	    /* not searching, we're done */
	    wrap->predict = do_training = 0;
	}
    }

#if defined(HAVE_MPI) && defined(_OPENMP)
    if (wrap->nproc > 1) {
	/* we're running more than one process, but from this point
	   only rank 0 is working, so we can assign more than one OMP
	   thread to it without worrying about contention?
	*/
	int nc = gretl_n_physical_cores();

	omp_set_num_threads(nc);
    }
#endif

    if (!err && do_training) {
	pputs(prn, _("Calling training function (this may take a while)\n"));
	gretl_flush(prn);
	if (wrap->do_probs) {
	    maybe_set_svm_seed(wrap);
	}
	model = svm_train(prob1, parm);
	if (model == NULL) {
	    err = E_DATA;
	}
	pprintf(prn, _("Training done, err = %d\n"), err);
	gretl_flush(prn);
    }

    if (model != NULL && saving_model(wrap)) {
	err = do_save_model(model, wrap, bmodel, prn);
	report_result(err, prn);
    }

    if (!err && wrap->predict > 0) {
	/* Now do some prediction */
	int training = (wrap->t2_train > 0);
	int T_os = -1;

	err = real_svm_predict(yhat, prob1, wrap, model, training, dset, prn);
	if (err) {
	    goto bailout;
	}
	*yhat_written = 1;
	dset->t2 = wrap->t2;
	if (training && wrap->predict > 1) {
	    T_os = dset->t2 - wrap->t2_train;
	}
	if (T_os >= 10) {
	    /* If we have some out-of-sample data, go
	       ahead and predict out of sample.
	    */
	    int T;

	    dset->t1 = wrap->t2_train + 1;
	    T = sample_size(dset);
	    pprintf(prn, _("Found %d testing observations\n"), T);
	    err = check_test_data(dset, wrap->ranges, wrap->k);
	    if (!err) {
		prob2 = gretl_sv_data_alloc(T, wrap->k, &x_space2, &err);
	    }
	    if (!err) {
		sv_data_fill(prob2, x_space2, wrap, list, dset, 2);
		err = real_svm_predict(yhat, prob2, wrap, model, 0, dset, prn);
	    }
	}
    }

 bailout:

    gretl_sv_data_destroy(prob2, x_space2);
    if (wrap->flags & W_LOADMOD) {
	gretl_destroy_svm_model(model);
    } else {
	svm_free_and_destroy_model(&model);
    }
    svm_prn = NULL;

    return err;
}

static int sv_trim_missing (int *list, int fvar, DATASET *dset)
{
    int t1 = dset->t1;
    int t2 = dset->t2;
    int T, nmiss = 0;
    int err = 0;

    if (fvar) {
	/* Temporary adjustment, since we don't require that
	   the folds variable (in last place) has valid values
	   at _all_ observations, but only in the training set;
	   and we'll check for that later.
	*/
	list[0] -= 1;
    }

    list_adjust_sample(list, &t1, &t2, dset, &nmiss);
    T = t2 - t1 + 1 - nmiss;

    if (fvar) {
	/* reset to the original value */
	list[0] += 1;
    }

    if (T < sample_size(dset)) {
	/* we can't handle NAs for the dependent variable within the
	   (possibly reduced) sample range
	*/
	int t, yno = list[1];

	for (t=t1; t<=t2 && !err; t++) {
	    if (na(dset->Z[yno][t])) {
		gretl_errmsg_sprintf(_("Dependent variable is NA at obs %d"), t+1);
		err = E_MISSDATA;
	    }
	}
	if (!err) {
	    fprintf(stderr, "sv_trim_missing: t1=%d, t2=%d, excluding %d observations\n",
		    t1+1, t2+1, sample_size(dset) - T);
	}
    }

    if (!err) {
	if (T > list[0]) {
	    dset->t1 = t1;
	    dset->t2 = t2;
	} else {
	    err = E_MISSDATA;
	}
    }

    return err;
}

/* If w->seed was automatic, and was actually used, save it
   to the incoming bundle in case it may be of interest.
*/

static void maybe_save_auto_seed (sv_wrapper *w, gretl_bundle *b)
{
    if (gretl_bundle_has_key(b, "seed")) {
	return; /* no, an incoming seed was specified */
    } else if (gretl_bundle_has_key(b, "foldvar") ||
	       gretl_bundle_has_key(b, "consecutive")) {
	return; /* no, randomization of folds not employed */
    } else if (gretl_bundle_has_key(b, "folds") ||
	       gretl_bundle_has_key(b, "search") ||
	       gretl_bundle_has_key(b, "grid")) {
	/* OK, randomized folds employed */
	gretl_bundle_set_uint32(b, "autoseed", w->seed);
    }
}

#ifdef HAVE_MPI

static int forward_to_gretlmpi (const int *list,
				gretl_bundle *bparams,
				DATASET *dset,
				int np,
				double *yhat,
				int *got_yhat,
				PRN *prn)
{
    gchar *fname;
    int err = 0;

    /* write data file for gretlmpi */
    fname = gretl_make_dotpath("svmtmp.gdt");
    err = write_data(fname, (int *) list, dset, OPT_NONE, NULL);
    g_free(fname);

    if (!err) {
	/* write parameter bundle for gretlmpi */
	err = gretl_bundle_write_to_file(bparams, "svmtmp.xml", 1);
    }

    if (!err) {
	/* compose and execute MPI script */
	if (np > 1) {
	    set_optval_int(MPI, OPT_N, np);
	    err = foreign_start(MPI, NULL, OPT_N, prn);
	} else {
	    err = foreign_start(MPI, NULL, OPT_NONE, prn);
	}
	if (!err) {
	    foreign_append("open @dotdir/svmtmp.gdt -q", MPI);
	    foreign_append("bundle parms = bread(\"@dotdir/svmtmp.xml\")", MPI);
	    foreign_append("list All = dataset", MPI);
	    foreign_append("yhat_svm = svm(All, parms)", MPI);
	    foreign_append("if $mpirank < 1", MPI);
	    foreign_append("  mwrite({yhat_svm}, \"@dotdir/svmtmp.mat\")", MPI);
	    foreign_append("endif", MPI);
	    err = foreign_execute(dset, OPT_L | OPT_Q, prn); /* note: --local */
	}
    }

    if (!err) {
	gretl_matrix *m;
	int s, t;

	m = gretl_matrix_read_from_file("svmtmp.mat", 1, &err);
	if (m == NULL) {
	    err = E_DATA;
	} else {
	    if (m->rows == sample_size(dset)) {
		for (t=dset->t1, s=0; t<=dset->t2; t++, s++) {
		    yhat[t] = m->val[s];
		}
		*got_yhat = 1;
	    } else {
		err = E_DATA;
	    }
	    gretl_matrix_free(m);
	}
    }

    return err;
}

#endif

int gretl_svm_driver (const int *list,
		      gretl_bundle *bparams,
		      gretl_bundle *bmodel,
		      gretl_bundle *bprob,
		      double *yhat,
		      int *yhat_written,
		      DATASET *dset,
		      PRN *inprn)
{
    sv_parm parm;
    sv_wrapper wrap;
    sv_data *prob = NULL;
    sv_cell *x_space = NULL;
    PRN *prn = inprn;
    int save_t1 = dset->t1;
    int save_t2 = dset->t2;
    int fvar = 0, lmin = 2;
    int T, err = 0;

#ifdef HAVE_MPI
    if (gretl_mpi_rank() > 0) {
	/* cut out chatter from ranks other than 0 */
	prn = NULL;
    } else if (auto_mpi_ok()) {
	int np = peek_use_mpi(bparams, bmodel, bprob, prn);

	if (np > 0) {
	    return forward_to_gretlmpi(list, bparams, dset, np,
				       yhat, yhat_written, prn);
	}
    }
#endif

    /* look ahead to see if we're supposed to be getting
       a "folds" variable at the end of @list */
    err = peek_foldvar(bparams, &fvar);
    if (err) {
	return err;
    } else if (fvar) {
	lmin = 3;
    }

    /* general checking and initialization */
    if (list == NULL || list[0] < lmin) {
	gretl_errmsg_set(_("svm: invalid list argument"));
	err = E_INVARG;
    } else {
	/* adjust the sample range for any NAs fore or aft */
	err = sv_trim_missing((int *) list, fvar, dset);
	if (!err) {
	    sv_wrapper_init(&wrap, dset);
	    err = read_params_bundle(bparams, bmodel, &wrap, &parm,
				     list, dset, prn);
	}
    }

    if (err) {
	/* restore incoming sample and get out */
	dset->t1 = save_t1;
	dset->t2 = save_t2;
	return err;
    }

#ifdef HAVE_MPI
    wrap.nproc = gretl_mpi_n_processes();
    if (wrap.nproc > 1) {
	/* we're actually in MPI mode */
	wrap.rank = gretl_mpi_rank();
	worker_prn = NULL;
	if (wrap.rank > 0) {
	    /* let rank 0 do (almost) all the talking */
	    if (!(wrap.flags & W_QUIET)) {
		worker_prn = inprn;
		wrap.flags |= W_QUIET;
	    }
	}
    }
#endif

    if (wrap.flags & W_QUIET) {
	svm_set_print_string_function(gretl_libsvm_noprint);
	prn = NULL;
    } else {
	svm_prn = prn;
	svm_set_print_string_function(gretl_libsvm_print);
    }

    if (wrap.t2_train > 0) {
	/* adjust sample temporarily */
	dset->t2 = wrap.t2_train;
    }
    T = sample_size(dset);

    err = get_svm_ranges(list, dset, &wrap, prn);

    if (!err) {
	wrap.k = (int) gretl_matrix_get(wrap.ranges, 0, 2) - 1;
	pputs(prn, _("Allocating problem space... "));
	prob = gretl_sv_data_alloc(T, wrap.k, &x_space, &err);
	report_result(err, prn);
    }
    gretl_flush(prn);

    if (!err) {
	/* fill out the "problem" data */
	pputs(prn, _("Scaling and transcribing data... "));
	sv_data_fill(prob, x_space, &wrap, list, dset, 1);
	if (parm.svm_type < 0) {
	    parm.svm_type = wrap.auto_type;
	}
	if (parm.gamma == 0) {
	    parm.gamma = 1.0 / wrap.k;
	}
	pputs(prn, "OK\n");
	/* we're now ready to run a full check on @parm */
	err = check_svm_params(prob, &parm, &wrap, prn);
	if (!err && parm.probability) {
	    wrap.do_probs = 1;
	}
    }

    if (!err && (wrap.flags & W_SEARCH)) {
	err = do_search_prep(prob, &parm, &wrap, dset, prn);
    }

    if (err) {
	goto getout;
    }

    if (wrap.nproc < 2 || !(wrap.flags & W_SEARCH)) {
	/* not doing parameter search via MPI */
	err = svm_predict_main(list, bmodel,
			       yhat, yhat_written,
			       dset, prn, &parm,
			       &wrap, prob);
	goto getout;
    }

#ifdef HAVE_MPI
    if (wrap.rank == 0) {
	err = svm_predict_main(list, bmodel,
			       yhat, yhat_written,
			       dset, prn, &parm,
			       &wrap, prob);
    } else {
	err = call_mpi_cross_validation(prob, &parm, &wrap, NULL);
    }
#endif

 getout:

    if (!err && bprob != NULL) {
	save_probs_to_bundle(&wrap, bprob);
    }

    if (!err && bparams != NULL) {
	maybe_save_auto_seed(&wrap, bparams);
    }

    gretl_sv_data_destroy(prob, x_space);
    svm_destroy_param(&parm);
    sv_wrapper_free(&wrap);

    dset->t1 = save_t1;
    dset->t2 = save_t2;

    return err;
}
