/* bootstrapped confidence intervals for impulse response functions */

#include "matrix_extra.h"

#if BDEBUG
# define BOOT_ITERS 5
#else
# define BOOT_ITERS 999
#endif

typedef struct irfboot_ irfboot;

struct irfboot_ {
    int n;              /* number of observations */
    int neqns;          /* number of equations */
    int order;          /* VAR lag order */
    int ecm;            /* 0 for an ordinary VAR, 1 for VECM */
    int rank;           /* rank, in case of VECM */
    gretlopt opt;       /* VECM options flags */
    int ncoeff;         /* number of coefficients per equation */
    int ifc;            /* equations include an intercept? */
    int t1;             /* starting observation */
    int t2;             /* ending observation */
    int horizon;        /* horizon for impulse responses */
    double **Z;         /* artificial dataset */
    DATAINFO *dinfo;    /* datainfo for artificial dataset */
    int **lists;        /* regression lists for models */
    gretl_matrix *A;    /* augmented coefficient matrix */
    gretl_matrix *C;    /* error covariance matrix */
    gretl_matrix *E;    /* matrix of residuals */
    gretl_matrix *S;    /* covariance matrix of residuals */
    gretl_matrix *rE;   /* matrix of resampled original residuals */
    gretl_matrix *rtmp; /* temporary storage */
    gretl_matrix *ctmp; /* temporary storage */
    gretl_matrix *resp; /* impulse response matrix */
    gretl_matrix *C0;   /* initial coefficient estimates (VECM only) */
    int *sample;        /* resampling array */
};

static void irf_boot_free (irfboot *boot)
{
    int i, nl = (boot->ecm)? 2 : boot->neqns;

    destroy_dataset(boot->Z, boot->dinfo);

    if (boot->lists != NULL) {
	for (i=0; i<nl; i++) {
	    free(boot->lists[i]);
	}
	free(boot->lists);
    }

    gretl_matrix_free(boot->A);
    gretl_matrix_free(boot->C);
    gretl_matrix_free(boot->E);
    gretl_matrix_free(boot->S);
    gretl_matrix_free(boot->rE);
    gretl_matrix_free(boot->rtmp);
    gretl_matrix_free(boot->ctmp);
    gretl_matrix_free(boot->resp);
    gretl_matrix_free(boot->C0);

    free(boot->sample);
}

static int boot_A_C_init (irfboot *boot)
{
    int n = boot->neqns * boot->order;
    int i, j, err = 0;
    
    boot->A = gretl_matrix_alloc(n, n);

    if (boot->A == NULL) {
	err = E_ALLOC;
    } else {
	for (i=boot->neqns; i<n; i++) {
	    for (j=0; j<n; j++) {
		gretl_matrix_set(boot->A, i, j, (j == i - boot->neqns)? 1.0 : 0.0);
	    }
	}
    }

    if (!err) {
	boot->C = gretl_matrix_alloc(n, boot->neqns);
	if (boot->C == NULL) {
	    err = E_ALLOC;
	    gretl_matrix_free(boot->A);
	    boot->A = NULL;
	} else {
	    gretl_matrix_zero(boot->C);
	}
    }

    return err;
}

static int boot_tmp_init (irfboot *boot)
{
    int n = boot->neqns * (boot->order + boot->ecm);
    int err = 0;

    boot->rtmp = gretl_matrix_alloc(n, boot->neqns);
    if (boot->rtmp == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	boot->ctmp = gretl_matrix_alloc(n, boot->neqns);
	if (boot->ctmp == NULL) {
	    gretl_matrix_free(boot->rtmp);
	    boot->rtmp = NULL;
	    err = E_ALLOC;
	}
    }

    return err;
}

static gretlopt boot_get_opt (const GRETL_VAR *var)
{
    gretlopt opt = OPT_NONE;

    if (var->jinfo != NULL) {
	if (var->jinfo->code == J_NO_CONST) {
	    opt = OPT_N;
	} else if (var->jinfo->code == J_REST_CONST) {
	    opt = OPT_R;
	} else if (var->jinfo->code == J_REST_TREND) {
	    opt = OPT_A;
	} else if (var->jinfo->code == J_UNREST_TREND) {
	    opt = OPT_T;
	}

	if (var->jinfo->seasonals) {
	    opt |= OPT_D;
	}

	/* OPT_S: save (don't delete) the variables added
	   to the dataset in the Johansen process.
	   OPT_B: flags the fact that we're doing impulse
	   response bootstrap.
	*/
	opt |= (OPT_S | OPT_B);
    }

    return opt;
}

static int irf_boot_init (irfboot *boot, const GRETL_VAR *var,
			  int periods)
{
    int err = 0;

    boot->Z = NULL;
    boot->dinfo = NULL;
    boot->lists = NULL;

    boot->C = NULL;
    boot->A = NULL;
    boot->E = NULL;
    boot->S = NULL;
    boot->rE = NULL;
    boot->rtmp = NULL;
    boot->ctmp = NULL;
    boot->resp = NULL;
    boot->C0 = NULL;

    boot->sample = NULL;

    boot->n = var->T;
    boot->neqns = var->neqns;
    boot->order = var->order;
    boot->ecm = var->ecm;
    boot->ifc = var->ifc;

    boot->t1 = var->t1;
    boot->t2 = var->t2;

    boot->horizon = periods;

    if (var->jinfo != NULL) {
	boot->rank = var->jinfo->rank;
	boot->ncoeff = var->ncoeff + (boot->neqns - boot->rank) + 
	    restricted(var);
    } else {
	boot->rank = 0;
	boot->ncoeff = var->ncoeff;
    }

    boot->opt = boot_get_opt(var);

#if BDEBUG
    fprintf(stderr, "boot: t1=%d, t2=%d, ncoeff=%d, horizon=%d\n",
	    boot->t1, boot->t2, boot->ncoeff, boot->horizon);
    fprintf(stderr, " n=%d, neqns=%d, order=%d, ecm=%d, ifc=%d\n",
	    boot->n, boot->neqns, boot->order, boot->ecm, boot->ifc);
#endif

    err = boot_tmp_init(boot);
    if (err) {
	goto bailout;
    }

    if (!boot->ecm) {
	/* we don't need these matrices for VECMs */
	err = boot_A_C_init(boot);
	if (err) {
	    goto bailout;
	}
	boot->E = gretl_matrix_alloc(boot->n, boot->neqns);
	boot->S = gretl_matrix_alloc(boot->neqns, boot->neqns);
	if (boot->E == NULL || boot->S == NULL) {
	    err = E_ALLOC;
	    goto bailout;
	}	    
    }

    boot->rE = gretl_matrix_alloc(boot->n, boot->neqns);
    boot->resp = gretl_matrix_alloc(boot->horizon, BOOT_ITERS);
    boot->sample = malloc(boot->n * sizeof *boot->sample);
    if (boot->rE == NULL || boot->resp == NULL || 
	boot->sample == NULL) {
	err = E_ALLOC;
    }

 bailout:
    if (err) {
	irf_boot_free(boot);
    }

    return err;
}

static int 
recalculate_impulse_responses (irfboot *boot, const gretl_matrix *A, 
			       const gretl_matrix *C,
			       int targ, int shock, int iter)
{
    double rt;
    int t, err = 0;

    for (t=0; t<boot->horizon && !err; t++) {
	if (t == 0) {
	    /* initial estimated responses */
	    err = gretl_matrix_copy_values(boot->rtmp, C);
	} else {
	    /* calculate further estimated responses */
	    err = gretl_matrix_multiply(A, boot->rtmp, boot->ctmp);
	    gretl_matrix_copy_values(boot->rtmp, boot->ctmp);
	}

	if (!err) {
	    rt = gretl_matrix_get(boot->rtmp, targ, shock);
	    gretl_matrix_set(boot->resp, t, iter, rt);
	}
    }

    return err;
}

/* transcribe residuals and coefficients from re-estimated
   VAR model to the boot->E and boot->A matrices */

static void irf_record_VAR_data (irfboot *boot, const MODEL *pmod, int k)
{
    int v = 0, lag = 0;
    int start = pmod->ifc;
    int rowmax = start + boot->neqns * boot->order;
    int i, j;

    /* record residuals */
    for (i=0; i<boot->n; i++) {
	gretl_matrix_set(boot->E, i, k, pmod->uhat[pmod->t1 + i]);
    }

    /* record coefficients */
    for (i=start; i<rowmax; i++) {
	if ((i - start) % boot->order == 0) {
	    v++;
	    lag = 1;
	} else {
	    lag++;
	}
	j = (lag - 1) * boot->neqns + v - 1;
	gretl_matrix_set(boot->A, k, j, pmod->coeff[i]);
    }
}

/* In re-estimation of VAR or VECM we'll tolerate at most 4 cases of
   near-perfect collinearity (which can arise by chance): maybe this
   should be more flexible? 
*/

#define MAXSING 5
#define VAR_FATAL(e,i,s) (e && (e != E_SINGULAR || i == 0 || s >= MAXSING))

static int 
re_estimate_VECM (irfboot *boot, GRETL_VAR *jvar, int targ, int shock, 
		  int iter, int scount)
{
    static int (*jbr) (GRETL_VAR *, double **, DATAINFO *, int) = NULL;
    static void *handle = NULL;
    int err = 0;

    gretl_error_clear();

    if (iter == 0) {
	/* first round: open the Johansen plugin */
	jbr = get_plugin_function("johansen_boot_round", &handle);
	if (jbr == NULL) {
	    err = 1;
	}
    }

    if (!err) {
	/* FIXME restricted vecm */
	err = johansen_driver(jvar, NULL, (const double **) boot->Z, 
			      boot->dinfo, boot->opt, NULL);
    }

    if (!err) {
	/* call the plugin function */
	err = jbr(jvar, (const double **) boot->Z, boot->dinfo, iter);
    }

    if (!err) {   
#if BDEBUG
	gretl_matrix_print(jvar->S, "jvar->S (Omega)");
#endif
	err = gretl_VAR_do_error_decomp(jvar->S, jvar->C);
    }

    if (!err) {
	recalculate_impulse_responses(boot, jvar->A, jvar->C, targ, shock, iter);
    }

    if (iter == BOOT_ITERS - 1 || VAR_FATAL(err, iter, scount)) {
	if (handle != NULL) {
	    close_plugin(handle);
	    handle = NULL;
	    jbr = NULL;
	}
    }

    return err;
}

static int re_estimate_VAR (irfboot *boot, int targ, int shock, int iter)
{
    MODEL var_model;
    gretlopt lsqopt = OPT_A | OPT_Z;
    int i, err = 0;

    for (i=0; i<boot->neqns && !err; i++) {
	var_model = lsq(boot->lists[i], &boot->Z, boot->dinfo, VAR, lsqopt);
	err = var_model.errcode;
	if (!err) {
	    irf_record_VAR_data(boot, &var_model, i);
	}
	clear_model(&var_model);
    }

    if (!err) {    
	gretl_matrix_multiply_mod(boot->E, GRETL_MOD_TRANSPOSE,
				  boot->E, GRETL_MOD_NONE,
				  boot->S, GRETL_MOD_NONE);
	gretl_matrix_divide_by_scalar(boot->S, boot->n);
	err = gretl_VAR_do_error_decomp(boot->S, boot->C);
    }

    if (!err) {
	recalculate_impulse_responses(boot, boot->A, boot->C, targ, shock, iter);
    }

    return err;
}

/* Allocate storage for the regression lists that will be used for the
   bootstrap VAR models; or in case of a VECM, create the "master" lists
   of endogenous and exogenous variables.
*/

static int allocate_bootstrap_lists (irfboot *boot, const GRETL_VAR *var)
{
    int nl = (var->ecm)? 2 : boot->neqns;
    int i, err = 0;

    boot->lists = malloc(nl * sizeof *boot->lists);

    if (boot->lists == NULL) {
	return E_ALLOC;
    }

    for (i=0; i<nl; i++) {
	boot->lists[i] = NULL;
    }

    if (boot->ecm) {
#if BDEBUG
	printlist(var->jinfo->list, "var->jinfo->list");
	printlist(var->jinfo->exolist, "var->jinfo->exolist");
#endif
	/* these lists will be adjusted shortly */
	boot->lists[0] = gretl_list_copy(var->jinfo->list);
	boot->lists[1] = gretl_list_copy(var->jinfo->exolist);
	if (boot->lists[0] == NULL || 
	    (var->jinfo->exolist != NULL && boot->lists[1] == NULL)) {
	    err = 1;
	}
    } else {
	for (i=0; i<nl && !err; i++) {
	    boot->lists[i] = gretl_list_new(boot->ncoeff + 1);
	    if (boot->lists[i] == NULL) {
		err = E_ALLOC;
	    }
	}
    }

    if (err) {
	for (i=0; i<nl; i++) {
	    free(boot->lists[i]);
	}
	free(boot->lists);
	boot->lists = NULL;
    }

    return err;
}

/* Initial construction of temporary dataset for use with the
   bootstrap procedure, by copying variables from the main dataset.
   Build the appropriate lists while we're at it.
*/

static int make_bs_dataset_and_lists (irfboot *boot,
				      const GRETL_VAR *var,
				      const double **Z,
				      const DATAINFO *pdinfo)
{
    double **bZ = NULL;
    DATAINFO *binfo = NULL;
    int ns, nd, nv;
    int i, j, k;
    int v, t;
    int err = 0;

    err = allocate_bootstrap_lists(boot, var);
    if (err) {
	return err;
    }

    if (var->ecm) {
	ns = var->neqns;
	nd = (boot->lists[1] != NULL)? boot->lists[1][0] : 0;
	nv = ns + nd + 1;
    } else {
	ns = var->order * var->neqns;     /* stochastic regressors per equation */
	nd = boot->ncoeff - ns;           /* number of deterministic regressors */
	nv = boot->ncoeff + boot->neqns;  /* total variables required */
	if (!boot->ifc) {
	    /* a gretl dataset includes a constant, regardless */
	    nv++;
	}
    }

#if BDEBUG
    fprintf(stderr, "make_bs_dataset_and_lists: ns=%d, nd=%d, nv=%d\n", ns, nd, nv);
#endif

    binfo = create_new_dataset(&bZ, nv, pdinfo->n, 0);
    if (binfo == NULL) {
	return E_ALLOC;
    }

    copy_dataset_obs_info(binfo, pdinfo);
    binfo->t1 = pdinfo->t1;
    binfo->t2 = pdinfo->t2;
    boot->Z = bZ;
    boot->dinfo = binfo;

    if (var->ecm) {
	j = 1;
	for (i=0; i<var->neqns; i++) {
	    v = boot->lists[0][i+1];
	    for (t=0; t<pdinfo->n; t++) {
		bZ[j][t] = Z[v][t];
	    }
	    boot->lists[0][i+1] = j;
	    j++;
	}

	for (i=0; i<nd; i++) {
	    v = boot->lists[1][i+1];
	    for (t=0; t<pdinfo->n; t++) {
		bZ[j][t] = Z[v][t];
	    }
	    boot->lists[1][i+1] = j;
	    j++;
	}	    
    } else {
	int lv = var->neqns + 1;
	int dv = lv + ns;

	for (i=0; i<var->neqns; i++) {
	    /* copy stochastic vars into positions 1 to var->neqns */
	    v = var->ylist[i+1];
	    for (t=0; t<pdinfo->n; t++) {
		bZ[i+1][t] = Z[v][t];
	    }

	    /* create lags */
	    for (j=1; j<=var->order; j++) {
		for (t=0; t<pdinfo->n; t++) {
		    bZ[lv][t] = (t - j >= 0)? Z[v][t-j] : NADBL;
		}
		lv++;
	    }
	}

	/* exogenous vars */
	if (var->xlist != NULL) {
	    for (j=1; j<=var->xlist[0]; j++) {
		v = var->xlist[j];
		for (t=0; t<pdinfo->n; t++) {
		    bZ[dv][t] = Z[v][t];
		}
		dv++;
	    }
	}

	/* seasonals? */
	if (var->detflags & DET_SEAS) {
	    for (t=0; t<pdinfo->n; t++) {
		for (j=0; j<pdinfo->pd - 1; j++) {
		    /* FIXME */
		    bZ[dv+j][t] = (1)? 1 : 0;
		}
	    }
	    dv += pdinfo->pd - 1;
	}

	/* trend? */
	if (var->detflags & DET_TREND) {
	    for (t=0; t<pdinfo->n; t++) {
		bZ[dv][t] = (double) (t + 1);
	    }
	}		

	/* compose lists */
	for (i=0; i<var->neqns; i++) {
	    lv = var->neqns + 1;
	    dv = lv + ns;
	    boot->lists[i][1] = i+1;
	    k = 2;
	    if (boot->ifc) {
		boot->lists[i][k++] = 0;
	    }
	    for (j=0; j<ns; j++) {
		boot->lists[i][k++] = lv++;
	    }
	    for (j=0; j<nd-boot->ifc; j++) {
		boot->lists[i][k++] = dv++;
	    }
	}
    }

    for (i=1; i<binfo->v; i++) {
	sprintf(binfo->varname[i], "bZ%d", i);
    }

#if BDEBUG
    if (1) {
	int nl = (boot->ecm)? 2 : boot->neqns;

	for (i=0; i<nl; i++) {
	    printlist(boot->lists[i], "IRF boot list");
	}
    }	
#endif

    return err;
}

/* VECM: make a vector, rbeta, containing the coefficients on the
   "restricted" constant or trend, if this is required.
*/

static gretl_matrix *make_restricted_coeff_vector (const GRETL_VAR *var)
{
    gretl_matrix *b = gretl_column_vector_alloc(var->jinfo->rank);
    gretl_matrix *rbeta;
    double x;
    int j;

    if (b == NULL) {
	return NULL;
    }

    rbeta = gretl_column_vector_alloc(var->neqns);
    if (rbeta == NULL) {
	gretl_matrix_free(b);
	return NULL;
    }
    
    for (j=0; j<var->jinfo->rank; j++) {
	/* extract the last row of \beta, transposed */
	x = gretl_matrix_get(var->jinfo->Beta, var->neqns, j);
	gretl_vector_set(b, j, x);
    }

    gretl_matrix_multiply(var->jinfo->Alpha, b, rbeta);

#if BDEBUG > 1
    gretl_matrix_print(b, "restricted var row of beta'");
    gretl_matrix_print(rbeta, "coeffs for restricted term");
#endif

    gretl_matrix_free(b);

    return rbeta;
}

/* VECM: make a record of all the original coefficients in the VAR
   representation of the VECM, so we have these on hand in convenient
   form for recomputing the dataset after resampling the residuals,
   or for forecasting.  

   We get these coefficients from three sources: (1) the A or
   "companion" matrix; (2) the VECM models (unrestricted constant
   and/or trend, seasonals if applicable); and (3) the implied
   coefficient on a restricted constant or trend.

   Note that the construction here depends on the order in which
   variables are added to the dataset, and to the regression lists, in
   johansen_VAR_prepare(): that order can't be changed without
   breaking stuff here.
*/

static gretl_matrix *VAR_coeff_matrix_from_VECM (const GRETL_VAR *var)
{
    gretl_matrix *C0 = NULL;
    gretl_matrix *rbeta = NULL;
    int order = var->order + 1;
    int nexo = (var->jinfo->exolist != NULL)? var->jinfo->exolist[0] : 0;
    int ndelta = var->order * var->neqns;
    int nseas = var->jinfo->seasonals;
    int ncoeff;
    int X0, S0, T0;
    double aij;
    int i, j, k;

    /* total coeffs in VAR representation */
    ncoeff = var->ncoeff + (var->neqns - var->jinfo->rank) + 
	restricted(var);

    if (restricted(var)) {
	rbeta = make_restricted_coeff_vector(var);
	if (rbeta == NULL) {
	    return NULL;
	}
    }

    C0 = gretl_matrix_alloc(var->neqns, ncoeff);
    if (C0 == NULL) {
	gretl_matrix_free(rbeta);
	return NULL;
    }

    /* position of first exog var coeff in VECM model */
    X0 = var->ifc + ndelta;
    /* position of first seasonal coeff in VECM model */
    S0 = X0 + nexo;
    /* position of trend coeff in VECM model */
    T0 = S0 + nseas;

    for (i=0; i<var->neqns; i++) {
	const MODEL *pmod = var->models[i];
	int col = 0;

	/* constant, if present */
	if (var->ifc) {
	    gretl_matrix_set(C0, i, col++, pmod->coeff[0]);
	}

	/* endogenous vars: use companion matrix */
	for (j=0; j<var->neqns; j++) {
	    for (k=0; k<order; k++) {
		aij = gretl_matrix_get(var->A, i, k * var->neqns + j);
		gretl_matrix_set(C0, i, col++, aij);
	    }
	}

	/* exogenous vars */
	for (j=0; j<nexo; j++) {
	    gretl_matrix_set(C0, i, col++, pmod->coeff[X0+j]);
	}

	/* seasonals, if present */
	for (j=0; j<nseas; j++) {
	    gretl_matrix_set(C0, i, col++, pmod->coeff[S0+j]);
	}

	/* unrestricted trend, if present */
	if (jcode(var) == J_UNREST_TREND) {
	    gretl_matrix_set(C0, i, col++, pmod->coeff[T0]);
	}

	/* restricted term (const or trend), if present */
	if (rbeta != NULL) {
	    aij = gretl_vector_get(rbeta, i);
	    gretl_matrix_set(C0, i, col, aij);
	} 	
    }

    if (rbeta != NULL) {
	gretl_matrix_free(rbeta);
    }

#if BDEBUG > 1
    gretl_matrix_print(var->A, "var->A");
    gretl_matrix_print(C0, "C0");
#endif

    return C0;
}

#if BDEBUG > 1
static void print_boot_dataset (const irfboot *boot)
{
    PRN *prn = gretl_print_new(GRETL_PRINT_STDERR);
    int t1 = boot->dinfo->t1;
    int *list;
    int i;

    list = gretl_list_new(boot->dinfo->v - 1);
    for (i=1; i<boot->dinfo->v; i++) {
	list[i] = i;
    }
    fprintf(stderr, "boot->dinfo->t1 = %d, boot->dinfo->t2 = %d\n",
	    boot->dinfo->t1, boot->dinfo->t2);
    boot->dinfo->t1 = 0;
    printdata(list, (const double **) boot->Z, boot->dinfo, OPT_O, prn);
    boot->dinfo->t1 = t1;
    gretl_print_destroy(prn);
    free(list);
}
#endif

static void 
compute_VECM_dataset (irfboot *boot, const GRETL_VAR *var, int iter)
{
    int order = var->order + 1;
    int nexo = (boot->lists[1] != NULL)? boot->lists[1][0] : 0;
    int nseas = var->jinfo->seasonals;
    int i, j, k, vj, t;

#if BDEBUG
    fprintf(stderr, "compute_VECM_dataset: order=%d, nexo=%d, nseas=%d, t1=%d\n",
	    order, nexo, nseas, boot->t1);
#endif

    for (t=boot->t1; t<=boot->t2; t++) {
	for (i=0; i<boot->neqns; i++) {
	    double cij, eti, bti = 0.0;
	    int col = 0;

	    /* unrestricted constant, if present */
	    if (var->ifc) {
		cij = gretl_matrix_get(boot->C0, i, col++);
		bti += cij;
	    }

	    /* lags of endogenous vars */
	    for (j=0; j<boot->neqns; j++) {
		vj = boot->lists[0][j+1];
		for (k=0; k<order; k++) {
		    cij = gretl_matrix_get(boot->C0, i, col++);
		    bti += cij * boot->Z[vj][t-k-1];
		}
	    }

	    /* exogenous vars, if present */
	    for (j=0; j<nexo; j++) {
		vj = boot->lists[1][j+1];
		cij = gretl_matrix_get(boot->C0, i, col++);
		bti += cij * boot->Z[vj][t];
	    }

	    /* seasonals, if present */
	    for (j=0; j<nseas; j++) {
		vj = boot->neqns + nexo + 1 + j;
		cij = gretl_matrix_get(boot->C0, i, col++);
		bti += cij * boot->Z[vj][t];
	    }

	    if (jcode(var) == J_UNREST_TREND) {
		/* unrestricted trend */
		cij = gretl_matrix_get(boot->C0, i, col);
		bti += cij * (t + 1);
	    } else if (jcode(var) == J_REST_CONST) {
		/* restricted constant */
		bti += gretl_matrix_get(boot->C0, i, col);
	    } else if (jcode(var) == J_REST_TREND) {
		/* restricted trend */
		cij = gretl_matrix_get(boot->C0, i, col);
		bti += cij * t;
	    }

	    /* set value of dependent var to fitted + re-sampled error */
	    eti = gretl_matrix_get(boot->rE, t - boot->t1, i);
	    boot->Z[i+1][t] = bti + eti;
	}
    }

    if (iter > 0) {
	int vl = 1 + var->neqns + nexo + nseas;
	int vd = vl + var->neqns;

	if (nseas > 0) {
	    /* allow for the unused seasonal */
	    vl++;
	    vd++;
	}

	/* recompute first lags and first differences */
	for (i=0; i<boot->neqns; i++) {
	    for (t=1; t<boot->dinfo->n; t++) {
		boot->Z[vl][t] = boot->Z[i+1][t-1];
		boot->Z[vd][t] = boot->Z[i+1][t] - boot->Z[i+1][t-1];
	    }
	    vl++;
	    vd++;
	}
    }

#if BDEBUG > 1
    print_boot_dataset(boot);
#endif
}

/* (Re-)fill the bootstrap dataset with artificial data, based on the
   re-sampled residuals from the original VAR (simple VAR, not VECM).
*/

static void compute_VAR_dataset (irfboot *boot, const GRETL_VAR *var)
{
    const MODEL *pmod;
    int ns = boot->order * boot->neqns;
    double xti, bti, eti;
    int i, j, vj, t;

    for (t=boot->t1; t<=boot->t2; t++) {
	int lv = boot->neqns + 1;

	for (i=0; i<boot->neqns; i++) {
	    int m, lag = 1;

	    pmod = var->models[i];
	    bti = 0.0;

	    for (j=0; j<boot->ncoeff; j++) {
		vj = boot->lists[i][j+2];
		if (j < ns + boot->ifc && vj > 0) {
		    /* stochastic variable */
		    m = (j - boot->ifc) / boot->order;
		    vj = boot->lists[m][1];
		    xti = boot->Z[vj][t-lag];
		    lag++;
		    if (lag > boot->order) {
			lag = 1;
		    }
		} else {
		    /* exogenous variable */
		    xti = boot->Z[vj][t];
		}
		bti += pmod->coeff[j] * xti;
	    }

	    /* set value of dependent var to forecast + re-sampled error */
	    eti = gretl_matrix_get(boot->rE, t - boot->t1, i);
	    boot->Z[i+1][t] = bti + eti;

	    /* and recreate the lags */
	    for (j=1; j<=boot->order; j++) {
		if (t - j >= boot->t1) {
		    boot->Z[lv][t] = boot->Z[i+1][t-j];
		}
		lv++;
	    }
	}
    }

#if BDEBUG > 1
    print_boot_dataset(boot);
#endif
}

/* Resample the original VAR or VECM residuals.  Note the option to
   "resample" _without_ actually changing the order, if BDEBUG > 1.
   This is useful for checking that the IRF bootstrap rounds are
   idempotent: we should then get exactly the same set of responses as
   in the original estimation of the VAR/VECM.
*/

static void resample_resids (irfboot *boot, const GRETL_VAR *var)
{
    double e;
    int i, s;

    /* construct sampling array */
    for (s=0; s<boot->n; s++) {
#if BDEBUG > 1
	boot->sample[s] = s;
#else
	boot->sample[s] = gretl_rand_int_max(boot->n);
#endif
#if BDEBUG
	fprintf(stderr, "boot->sample[%d] = %d\n", s, boot->sample[s]);
#endif
    }

    /* draw sample from the original residuals */
    for (s=0; s<boot->n; s++) {
	for (i=0; i<boot->neqns; i++) {
	    e = gretl_matrix_get(var->E, boot->sample[s], i);
	    gretl_matrix_set(boot->rE, s, i, e);
	}
    }
}

static int irf_boot_quantiles (irfboot *boot, gretl_matrix *R)
{
    double alpha = 0.05;
    double *rk;
    int k, ilo, ihi;

    rk = malloc(BOOT_ITERS * sizeof *rk);
    if (rk == NULL) {
	return E_ALLOC;
    }

#if BDEBUG
    ilo = 1;
    ihi = BOOT_ITERS;
    fprintf(stderr, "IRF bootstrap (%d iters), min and max values\n", BOOT_ITERS);
#else
    ilo = (BOOT_ITERS + 1) * alpha / 2.0;
    ihi = (BOOT_ITERS + 1) * (1.0 - alpha / 2.0);
#endif

    for (k=0; k<boot->horizon; k++) {
	gretl_matrix_row_to_array(boot->resp, k, rk);
	qsort(rk, BOOT_ITERS, sizeof *rk, gretl_compare_doubles);
	gretl_matrix_set(R, k, 1, rk[ilo-1]);
	gretl_matrix_set(R, k, 2, rk[ihi-1]);
    }

    free(rk);

    return 0;
}

static gretl_matrix *irf_bootstrap (const GRETL_VAR *var, 
				    int targ, int shock, int periods,
				    const double **Z, 
				    const DATAINFO *pdinfo)
{
    GRETL_VAR *jvar = NULL;
    gretl_matrix *R;
    int *list = NULL; /* FIXME */
    irfboot boot;
    int scount = 0;
    int iter, err = 0;

#if BDEBUG
    fprintf(stderr, "\n*** irf_bootstrap() called\n");
#endif

    R = gretl_matrix_alloc(periods, 3);
    if (R == NULL) {
	return NULL;
    }

    err = irf_boot_init(&boot, var, periods);
    if (err) {
	goto bailout;
    }

    err = make_bs_dataset_and_lists(&boot, var, Z, pdinfo);

    if (var->ecm && !err) {
	boot.C0 = VAR_coeff_matrix_from_VECM(var);
	if (boot.C0 == NULL) {
	    err = E_ALLOC;
	}
	if (!err) {
	    jvar = johansen_VAR_new(boot.rank, boot.order + 1, list, 
				    (const double **) boot.Z, 
				    boot.dinfo, boot.opt, &err);
	    if (jvar != NULL) {
		err = jvar->err;
	    }
	}
    }

    for (iter=0; iter<BOOT_ITERS && !err; iter++) {
#if BDEBUG
	fprintf(stderr, "starting iteration %d\n", iter);
#endif
	resample_resids(&boot, var);
	if (var->ecm) {
	    compute_VECM_dataset(&boot, var, iter);
	    err = re_estimate_VECM(&boot, jvar, targ, shock, iter, scount);
	} else {
	    compute_VAR_dataset(&boot, var);
	    err = re_estimate_VAR(&boot, targ, shock, iter);
	}
	if (err && !VAR_FATAL(err, iter, scount)) {
	    /* excessive collinearity: try again, unless this is becoming a habit */
	    scount++;
	    iter--;
	    err = 0;
	}
    }

    if (err && scount == MAXSING) {
	strcpy(gretl_errmsg, "Excessive collinearity in resampled datasets");
    }

    if (!err) {
	err = irf_boot_quantiles(&boot, R);
    }

    irf_boot_free(&boot);

    if (jvar != NULL) {
	gretl_VAR_free(jvar);
    }

 bailout:

    if (err) {
	gretl_matrix_free(R);
	R = NULL;
    }

    return R;
}


