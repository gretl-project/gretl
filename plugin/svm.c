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

/* interface to libsvm for support vector machines */

#include "libgretl.h"
#include "version.h"

#include <svm.h>

typedef struct svm_problem sv_data;
typedef struct svm_node sv_cell;
typedef struct svm_model sv_model;

/*
   key enumerations from svm.h

   svm_type:
   C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR

   kernel_type:
   LINEAR, POLY, RBF, SIGMOID, PRECOMPUTED
*/

static void set_svm_param_defaults (struct svm_parameter *parm,
				    int mtype)
{
    parm->svm_type = mtype;
    parm->kernel_type = RBF;
    parm->degree = 3; /* for polynomial */
    parm->gamma = 0;  /* poly/rbf/sigmoid: default 1.0 / num_features */
    parm->coef0 = 0;  /* for use in kernel function */

    /* training-only variables */
    parm->cache_size = 1000;   /* cache size in MB */
    parm->eps = 0.001;         /* stopping criterion */
    parm->C = 1;               /* cost: for C_SVC, EPSILON_SVR and NU_SVR */
    parm->nr_weight = 0;       /* for C_SVC */
    parm->weight_label = NULL; /* for C_SVC */
    parm->weight = NULL;       /* for C_SVC */
    parm->nu = 0.5;            /* for NU_SVC, ONE_CLASS, and NU_SVR */
    parm->p = 0.1;             /* for EPSILON_SVR */
    parm->shrinking = 1;       /* use the shrinking heuristics */
    parm->probability = 0;     /* do probability estimates */
}

/* for testing against svm-scale */

static void print_ranges (const gretl_matrix *ranges, PRN *prn)
{
    double lo, hi;
    int i, idx;

    pprintf(prn, "x\n%d %d\n", (int) gretl_matrix_get(ranges, 0, 0),
	    (int) gretl_matrix_get(ranges, 0, 1));

    for (i=1; i<ranges->rows; i++) {
	idx = gretl_matrix_get(ranges, i, 0);
	lo = gretl_matrix_get(ranges, i, 1);
	hi = gretl_matrix_get(ranges, i, 2);
	pprintf(prn, "%d %.16g %.16g\n", idx, lo, hi);
    }
}

/* also for testing against svm-scale */

static void print_problem (sv_data *p, int k, PRN *prn)
{
    int i, t, idx;
    double val;

    for (t=0; t<p->l; t++) {
	pprintf(prn, "%g ", p->y[t]);
	for (i=0; i<k; i++) {
	    idx = p->x[t][i].index;
	    val = p->x[t][i].value;
	    if (val != 0) {
		pprintf(prn, "%d:%g ", idx, val);
	    }
	}
	pputc(prn, '\n');
    }
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
				     sv_cell **px_space)
{
    sv_data *p = malloc(sizeof *p);
    int err = 0;

    if (p != NULL) {
	p->l = T;
	p->y = malloc(T * sizeof *p->y);
	p->x = malloc(T * sizeof *p->x);
	if (p->y == NULL || p->x == NULL) {
	    err = E_ALLOC;
	} else {
	    /* we need an extra cell on each row to hold a
	       sentinel index value of -1
	    */
	    *px_space = malloc(T * (k+1) * sizeof(sv_cell));
	    if (*px_space == NULL) {
		err = E_ALLOC;
	    }
	}
	if (err) {
	    gretl_sv_data_destroy(p, NULL);
	    p = NULL;
	}
    }

    return p;
}

/* initial discovery of ranges of the RHS data using the
   training data */

static gretl_matrix *get_data_ranges (const int *list,
				      int scaling,
				      const DATASET *dset,
				      int *err)
{
    gretl_matrix *ranges;
    const double *x;
    double xmin, xmax;
    int k = list[0] - 1;
    int i, j, vi;

    ranges = gretl_matrix_alloc(k+1, 4);
    if (ranges == NULL) {
	*err = E_ALLOC;
	return NULL;
    }

    /* scaling limits */
    xmin = scaling == 2 ? 0 : -1;
    gretl_matrix_set(ranges, 0, 0, xmin); /* lower */
    gretl_matrix_set(ranges, 0, 1, 1);    /* upper */

    /* padding */
    gretl_matrix_set(ranges, 0, 2, 0);
    gretl_matrix_set(ranges, 0, 3, 0);

    /* note: we make no provision for scaling y at present */

    j = 0;
    for (i=2; i<=list[0]; i++) {
	vi = list[i];
	x = dset->Z[vi];
	gretl_minmax(dset->t1, dset->t2, x, &xmin, &xmax);
	if (xmin != xmax) {
	    j++;
	    gretl_matrix_set(ranges, j, 0, j); /* or... ? */
	    gretl_matrix_set(ranges, j, 1, xmin);
	    gretl_matrix_set(ranges, j, 2, xmax);
	    gretl_matrix_set(ranges, j, 3, vi);
	} else {
	    fprintf(stderr, "training data: dropping var %d (%s)\n",
		    vi, dset->varname[vi]);
	}
    }

    /* FIXME check for pathologies! (NAs, no non-constant
       series or whatever) */

    /* record number of rows actually occupied, which
       could be less than the number allocated
    */
    gretl_matrix_set(ranges, 0, 2, j + 1);

    return ranges;
}

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
	       record of its series ID to zero
	    */
	    gretl_matrix_set(ranges, i, 3, 0);
	}
    }

    if (n != k) {
	fprintf(stderr, "test data: number of usable variables (%d) "
		"differs from training data (%d)\n", n, k);
    } else {
	fprintf(stderr, "test data: number of usable variables "
		"agrees with training data\n");
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

#if 0
    if (value != 0) {
	new_num_nonzeros++;
    }
#endif

    return val;
}

static int sv_data_fill (sv_data *prob,
			 sv_cell *x_space, int k,
			 const gretl_matrix *ranges,
			 int scaling,
			 int *svm_type,
			 const int *list,
			 const DATASET *dset)
{
    double scalemin, scalemax;
    double xit, xmin, xmax;
    int i, j, s, t, vi, idx;
    int pos = 0;

    /* deal with the LHS variable */
    vi = list[1];
    if (svm_type != NULL &&
	(gretl_isdummy(dset->t1, dset->t2, dset->Z[vi]) ||
	 series_is_coded(dset, vi))) {
	/* classification, not regression */
	*svm_type = C_SVC;
    }
    for (i=0, t=dset->t1; t<=dset->t2; t++, i++) {
	prob->y[i] = dset->Z[vi][t];
    }

    /* retrieve the global x-scaling limits */
    scalemin = gretl_matrix_get(ranges, 0, 0);
    scalemax = gretl_matrix_get(ranges, 0, 1);

    /* write the scaled x-data into the problem struct */
    for (s=0, t=dset->t1; t<=dset->t2; t++, s++) {
	prob->x[s] = &x_space[pos];
	j = 0;
	for (i=1; i<=k; i++) {
	    vi = (int) gretl_matrix_get(ranges, i, 3);
	    if (vi <= 0) {
		/* may happen when we get to the test data */
		continue;
	    }
	    idx = (int) gretl_matrix_get(ranges, i, 0);
	    xmin = gretl_matrix_get(ranges, i, 1);
	    xmax = gretl_matrix_get(ranges, i, 2);
	    xit = dset->Z[vi][t];
	    if (scaling != 0) {
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

static int real_svm_predict (double *yhat,
			     sv_data *prob,
			     sv_model *model,
			     int training,
			     const DATASET *dset,
			     PRN *prn)
{
    const char *label;
    int n_correct = 0;
    int regression = 0;
    double ymean = 0;
    double TSS = 0.0;
    double SSR = 0.0;
    double dev, yhi;
    sv_cell *x;
    int i;

    if (model->param.svm_type == EPSILON_SVR ||
	model->param.svm_type == NU_SVR) {
	regression = 1;
	ymean = gretl_mean(0, prob->l - 1, prob->y);
    }

    pprintf(prn, "Calling prediction function (this may take a while)\n");
    gretl_gui_flush();
    gretl_print_flush_stream(prn);
    for (i=0; i<prob->l; i++) {
	x = prob->x[i];
	yhi = svm_predict(model, x);
	yhat[dset->t1 + i] = yhi;
	if (regression) {
	    dev = prob->y[i] - ymean;
	    TSS += dev * dev;
	    dev = prob->y[i] - yhi;
	    SSR += dev * dev;
	} else {
	    n_correct += (yhi == prob->y[i]);
	}
    }

    label = training ? "Training data" : "Test data";

    if (regression) {
	double r;
	
	r = gretl_corr(0, prob->l, prob->y, yhat + dset->t1, NULL);
	pprintf(prn, "%s: MSE = %g, R^2 = %g, squared corr = %g\n", label,
		SSR / prob->l, 1.0 - SSR / TSS, r * r);
    } else {
	pprintf(prn, "%s: correct predictions = %d (%.1f percent)\n", label,
		n_correct, 100 * n_correct / (double) prob->l);
    }	

    return 0;
}

static int parse_params_bundle (gretl_bundle *b,
				int *scaling,
				int *t2_train,
				const int *list,
				const DATASET *dset)
{
    int ival, err = 0;

    if (gretl_bundle_has_key(b, "scaling")) {
	ival = gretl_bundle_get_int(b, "scaling", &err);
	if (!err && (ival < 0 || ival > 2)) {
	    fprintf(stderr, "invalid 'scaling' arg %d\n", ival);
	    err = E_INVARG;
	}
	if (!err) {
	    *scaling = ival;
	}
    }

    if (!err) {
	ival = gretl_bundle_get_int(b, "t2_train", &err);
	if (!err && (ival < list[0] || ival > dset->n)) {
	    fprintf(stderr, "invalid 't2_train' arg %d\n", ival);
	    err = E_INVARG;
	}
	if (!err) {
	    *t2_train = ival - 1; /* zero-based */
	}
    }

    /* FIXME: add a lot more here, sooner or later! */

    return err;
}

static PRN *svm_prn;

static void gretl_libsvm_print (const char *s)
{
    if (svm_prn != NULL) {
	pputs(svm_prn, s);
	gretl_gui_flush();
    } else {
	fputs(s, stdout);
    }
}

int gretl_svm_predict (const int *list,
		       gretl_bundle *bparams,
		       gretl_bundle **bmodel,
		       double *yhat,
		       DATASET *dset,
		       PRN *prn)
{
    struct svm_parameter parm;
    gretl_matrix *ranges;
    sv_data *prob1 = NULL;
    sv_data *prob2 = NULL;
    sv_cell *x_space1 = NULL;
    sv_cell *x_space2 = NULL;
    sv_model *model = NULL;
    int svm_type = EPSILON_SVR;
    int save_t2 = dset->t2;
    int scaling = 1;
    int t2_train = 0;
    int T, k = 0;
    int err = 0;

    if (list == NULL || list[0] < 2) {
	fprintf(stderr, "svm: invalid list argument\n");
	err = E_INVARG;
    } else {
	err = parse_params_bundle(bparams, &scaling, &t2_train,
				  list, dset);
    }

    if (err) {
	return err;
    }

    svm_prn = prn;
    svm_set_print_string_function(gretl_libsvm_print);

    dset->t2 = t2_train;
    T = sample_size(dset);

    pprintf(prn, "Getting data ranges (sample = %d to %d)... ",
	    dset->t1 + 1, dset->t2 + 1);
    ranges = get_data_ranges(list, scaling, dset, &err);
    if (err) {
	pprintf(prn, "done, err = %d\n", err);
    } else {
	pputs(prn, "OK\n");
    }
    gretl_gui_flush();

    if (0 && !err) {
	/* just for testing */
	print_ranges(ranges, prn);
    }

    if (!err) {
	k = (int) gretl_matrix_get(ranges, 0, 2) - 1;
	pputs(prn, "Allocating problem space... ");
	prob1 = gretl_sv_data_alloc(T, k, &x_space1);
	if (prob1 == NULL) {
	    pprintf(prn, "err = %d\n", err);
	    err = E_ALLOC;
	} else {
	    pputs(prn, "OK\n");
	}
    }
    gretl_gui_flush();

    if (!err) {
	/* fill out the "problem" data */
	pputs(prn, "Scaling and transcribing data... ");
	sv_data_fill(prob1, x_space1, k, ranges, scaling,
		     &svm_type, list, dset);
	pputs(prn, "OK\n");
    }

    if (0 && !err) {
	/* just for testing */
	print_problem(prob1, k, prn);
    }

    if (!err) {
	set_svm_param_defaults(&parm, svm_type);
	if (parm.gamma == 0) {
	    parm.gamma = 1.0 / k;
	}
	pprintf(prn, "Calling training function (this may take a while)\n");
	gretl_gui_flush();
	model = svm_train(prob1, &parm);
	if (model == NULL) {
	    err = E_DATA;
	}
	gretl_gui_flush();
	pprintf(prn, "Training done, err = %d\n", err);
	gretl_gui_flush();
    }

    if (!err) {
	int T_os;

	real_svm_predict(yhat, prob1, model, 1, dset, prn);
	dset->t2 = save_t2;
	T_os = dset->t2 - t2_train;
	if (T_os >= t2_train) {
	    /* If we have enough out-of-sample data, go
	       ahead and predict out of sample.
	    */
	    dset->t1 = t2_train + 1;
	    T = sample_size(dset);
	    pprintf(prn, "Found %d testing observations\n", T);
	    err = check_test_data(dset, ranges, k);
	    if (!err) {
		prob2 = gretl_sv_data_alloc(T, k, &x_space2);
		if (prob2 == NULL) {
		    err = E_ALLOC;
		}
	    }
	    if (!err) {
		sv_data_fill(prob2, x_space2, k, ranges, scaling,
			     NULL, list, dset);
		real_svm_predict(yhat, prob2, model, 0, dset, prn);
	    }
	}
    }

    dset->t2 = save_t2;

    gretl_matrix_free(ranges);
    gretl_sv_data_destroy(prob1, x_space1);
    gretl_sv_data_destroy(prob2, x_space2);
    svm_free_and_destroy_model(&model);
    svm_destroy_param(&parm);
    svm_prn = NULL;

    return err;
}
