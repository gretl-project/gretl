/* bootstrapped confidence intervals for impulse response functions */

#define BDEBUG 1

#define BOOT_ITERS 999

typedef struct irfboot_ irfboot;

struct irfboot_ {
    int n;              /* number of observations */
    int neqns;          /* number of equations */
    int order;          /* VAR lag order */
    int nc;             /* number of coefficients per equation */
    int ifc;            /* equations include intercept? */
    int t1;             /* starting observation */
    int t2;             /* ending observation */
    int horizon;        /* horizon for impulse responses */
    double **Z;         /* artificial dataset */
    DATAINFO *dinfo;    /* datainfo for artificial dataset */
    int **lists;        /* regression lists for VAR models */
    gretl_matrix *C;    /* error covariance matrix */
    gretl_matrix *A;    /* augmented coefficient matrix */
    gretl_matrix *E;    /* matrix of residuals */
    gretl_matrix *S;    /* covariance matrix of residuals */
    gretl_matrix *rE;   /* matrix of resampled original residuals */
    gretl_matrix *rtmp; /* temporary storage */
    gretl_matrix *ctmp; /* temporary storage */
    gretl_matrix *resp; /* impulse response matrix */
    int *sample;        /* resampling array */
};

static void irf_boot_free (irfboot *boot)
{
    int i;

    destroy_dataset(boot->Z, boot->dinfo);

    if (boot->lists != NULL) {
	for (i=0; i<boot->neqns; i++) {
	    free(boot->lists[i]);
	}
	free(boot->lists);
    }

    gretl_matrix_free(boot->C);
    gretl_matrix_free(boot->A);
    gretl_matrix_free(boot->E);
    gretl_matrix_free(boot->S);
    gretl_matrix_free(boot->rE);
    gretl_matrix_free(boot->rtmp);
    gretl_matrix_free(boot->ctmp);
    gretl_matrix_free(boot->resp);

    free(boot->sample);
}

static int boot_A_C_init (irfboot *boot)
{
    int n = boot->neqns * boot->order; /* FIXME vecm */
    int err = 0;

    boot->A = gretl_matrix_alloc(n, n);
    if (boot->A == NULL) {
	err = 1;
    } else {
	int i, j;

	for (i=boot->neqns; i<n; i++) {
	    for (j=0; j<n; j++) {
		gretl_matrix_set(boot->A, i, j, (j == i - boot->neqns)? 1.0 : 0.0);
	    }
	}
    }

    if (!err) {
	boot->C = gretl_matrix_alloc(n, boot->neqns);
	if (boot->C == NULL) {
	    err = 1;
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
    int rows = boot->neqns * boot->order;
    int err = 0;

    boot->rtmp = gretl_matrix_alloc(rows, boot->neqns);
    if (boot->rtmp == NULL) {
	err = E_ALLOC;
    }

    if (!err) {
	boot->ctmp = gretl_matrix_alloc(rows, boot->neqns);
	if (boot->ctmp == NULL) {
	    gretl_matrix_free(boot->rtmp);
	    boot->rtmp = NULL;
	    err = E_ALLOC;
	}
    }

    return err;
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

    boot->sample = NULL;

    boot->n = var->T;
    boot->neqns = var->neqns;
    boot->order = var->order;
    boot->ifc = var->ifc;

    boot->t1 = var->models[0]->t1;
    boot->t2 = var->models[0]->t2;
    boot->nc = var->models[0]->list[0] - 1;

    boot->horizon = periods;

#if BDEBUG
    fprintf(stderr, "boot: t1=%d, t2=%d, nc=%d, horizon=%d\n",
	    boot->t1, boot->t2, boot->nc, boot->horizon);
    fprintf(stderr, " n=%d, neqns=%d, order=%d, ifc=%d\n",
	    boot->n, boot->neqns, boot->order, boot->ifc);
#endif

    err = boot_A_C_init(boot);
    if (err) {
	goto bailout;
    }

    err = boot_tmp_init(boot);
    if (err) {
	goto bailout;
    }

    boot->rE = gretl_matrix_alloc(boot->n, boot->neqns);
    if (boot->rE == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    boot->E = gretl_matrix_alloc(boot->n, boot->neqns);
    if (boot->E == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    boot->S = gretl_matrix_alloc(boot->neqns, boot->neqns);
    if (boot->S == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    boot->resp = gretl_matrix_alloc(boot->horizon, BOOT_ITERS);
    if (boot->resp == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    boot->sample = malloc(boot->n * sizeof *boot->sample);
    if (boot->sample == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

 bailout:
    if (err) {
	irf_boot_free(boot);
    }

    return err;
}

static int
recalculate_impulse_responses (irfboot *boot, int targ, int shock, int iter)
{
    double rt;
    int t, err = 0;

#if BDEBUG > 1
    fprintf(stderr, "\nrecalculate_impulse_responses: iteration %d\n", iter);
#endif

    for (t=0; t<boot->horizon && !err; t++) {
	if (t == 0) {
	    /* calculate initial estimated responses */
	    err = gretl_matrix_copy_values(boot->rtmp, boot->C);
	} else {
	    /* calculate further estimated responses */
	    err = gretl_matrix_multiply(boot->A, boot->rtmp, boot->ctmp);
	    gretl_matrix_copy_values(boot->rtmp, boot->ctmp);
	}

	if (!err) {
	    rt = gretl_matrix_get(boot->rtmp, targ, shock);
	    gretl_matrix_set(boot->resp, t, iter, rt);
	}

#if BDEBUG > 1
	fprintf(stderr, "resp[%d] = %g\n", t, gretl_matrix_get(boot->resp, t, iter));
#endif
    }

    return err;
}

static void irf_record_model_data (irfboot *boot, const MODEL *pmod, int k)
{
    int v = 0, lag = 0;
    int start = pmod->ifc;
    int rowmax = boot->neqns * boot->order + start;
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

static int re_estimate_VAR (irfboot *boot, int targ, int shock, int iter)
{
    MODEL var_model;
    int i, err = 0;

    /* changes needed here for VECM */

    for (i=0; i<boot->neqns && !err; i++) {
	var_model = lsq(boot->lists[i], &boot->Z, boot->dinfo, VAR, OPT_A, 0.0);
	err = var_model.errcode;
	if (!err) {
	    irf_record_model_data(boot, &var_model, i);
	}
	clear_model(&var_model);
    }

    if (!err) {    
	gretl_matrix_multiply_mod(boot->E, GRETL_MOD_TRANSPOSE,
				  boot->E, GRETL_MOD_NONE,
				  boot->S);
	gretl_matrix_divide_by_scalar(boot->S, boot->n);
	err = gretl_VAR_do_error_decomp(boot->neqns, boot->S, boot->C);
    }

    if (!err) {
	recalculate_impulse_responses(boot, targ, shock, iter);
    }

#if BDEBUG > 1
    fprintf(stderr, "\ntarg=%d, shock = %d:\n", targ, shock);
    for (i=0; i<boot->horizon && !err; i++) {
	fprintf(stderr, "resp[%d] = %g\n", i, boot->resp[i]);
    }
#endif

    return err;
}

/* Allocate storage for the regression lists that will be
   used for the bootstrap VAR models */

static int allocate_bootstrap_lists (irfboot *boot)
{
    int i, err = 0;

    boot->lists = malloc(boot->neqns * sizeof *boot->lists);
    if (boot->lists == NULL) {
	err = E_ALLOC;
    }

    for (i=0; i<boot->neqns && !err; i++) {
	boot->lists[i] = gretl_list_new(boot->nc + 1);
	if (boot->lists[i] == NULL) {
	    int j;

	    for (j=0; j<i; j++) {
		free(boot->lists[j]);
	    }
	    free(boot->lists);
	    boot->lists = NULL;
	    err = E_ALLOC;
	}
    }

    return err;
}

/* Construct a temporary dataset for use with the bootstrap
   procedure.  Also build the VAR regression lists.
*/

static int make_bs_dataset_and_lists (irfboot *boot,
				      const GRETL_VAR *var,
				      const double **Z,
				      const DATAINFO *pdinfo)
{
    const MODEL *pmod;
    double **bZ = NULL;
    DATAINFO *binfo = NULL;
    int ns, nd, nv;
    int i, j, k;
    int v, t;
    int err = 0;

    ns = var->order * var->neqns; /* stochastic regressors per equation */
    nd = boot->nc - ns;           /* number of deterministic regressors */
    nv = boot->nc + boot->neqns;  /* total variables required */

    if (!boot->ifc) {
	/* a gretl dataset includes a constant, regardless */
	nv++;
    }

    err = allocate_bootstrap_lists(boot);
    
    if (!err) {
	binfo = create_new_dataset(&bZ, nv, pdinfo->n, 0);
	if (binfo == NULL) {
	    err = E_ALLOC;
	}
    }

    if (!err) {
	copy_dataset_obs_info(binfo, pdinfo);
	binfo->t1 = pdinfo->t1;
	binfo->t2 = pdinfo->t2;
	boot->Z = bZ;
	boot->dinfo = binfo;
    }

    if (!err) {
	int lv = var->neqns + 1;
	int dv = lv + ns;

	fprintf(stderr, "dv = %d\n", dv);

	for (i=0; i<var->neqns; i++) {
	    /* copy stochastic vars into positions 1 to var->neqns */
	    pmod = var->models[i];
	    v = pmod->list[1];
#if BDEBUG
	    fprintf(stderr, "Copying %s as var %d\n", pdinfo->varname[v], i+1);
#endif
	    for (t=0; t<pdinfo->n; t++) {
		bZ[i+1][t] = Z[v][t];
	    }
	    /* create lags */
	    for (j=1; j<=var->order; j++) {
#if BDEBUG
		fprintf(stderr, "Creating lag %d of %s as var %d\n", 
			j, pdinfo->varname[v], lv);
#endif
		for (t=0; t<pdinfo->n; t++) {
		    bZ[lv][t] = (t - j >= 0)? Z[v][t-j] : NADBL;
		}
		lv++;
	    }
	    /* deterministic vars: copy once */
	    if (i == 0) {
		for (j=ns+2+pmod->ifc; j<=pmod->list[0]; j++) {
		    v = pmod->list[j];
#if BDEBUG
		    fprintf(stderr, "Copying %s as var %d\n", 
			    pdinfo->varname[v], dv);
#endif
		    for (t=0; t<pdinfo->n; t++) {
			bZ[dv][t] = Z[v][t];
		    }
		    dv++;
		}
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
#if BDEBUG
	    printlist(boot->lists[i], "IRF boot list");
#endif
	}
    }

    for (i=1; i<binfo->v; i++) {
	sprintf(binfo->varname[i], "bZ%d", i);
    }

    return err;
}

/* (Re-)fill the bootstrap dataset with artificial data, based on the
   re-sampled residuals from the original VAR.
*/

static void compute_bootstrap_dataset (irfboot *boot, const GRETL_VAR *var)
{
    const MODEL *pmod;
    int ns = boot->order * boot->neqns;
    double xti, bti, eti;
    int i, j, vj, t;

    for (t=boot->t1; t<=boot->t2; t++) {
	int lv = boot->neqns + 1;

	for (i=0; i<boot->neqns; i++) {
	    int m, k = 0, lag = 1;

	    pmod = var->models[i];
	    bti = 0.0;

	    for (j=0; j<boot->nc; j++) {
		vj = boot->lists[i][j+2];
		if (j < ns + boot->ifc && vj > 0) {
		    /* stochastic variable */
		    m = (j - boot->ifc) / boot->order;
		    vj = boot->lists[m][1];
		    xti = boot->Z[vj][t-lag];
		    lag++;
		    if (lag > boot->order) {
			lag = 1;
			k++;
		    }
		} else {
		    /* deterministic variable */
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
    if (1) {
	PRN *prn;
	int list[5];

	list[0] = 4;
	for (i=1; i<5; i++) {
	    list[i] = i;
	}
	prn = gretl_print_new(GRETL_PRINT_STDERR);
	printdata(list, (const double **) boot->Z, boot->dinfo, OPT_O, prn);
	gretl_print_destroy(prn);
    }
#endif
}

static void resample_resids (irfboot *boot, const GRETL_VAR *var)
{
    double e;
    int i, s;

    /* construct sampling array */
    for (s=0; s<boot->n; s++) {
	boot->sample[s] = gretl_rand_int_max(boot->n);
#if BDEBUG > 1
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

    ilo = (BOOT_ITERS + 1) * alpha / 2.0;
    ihi = (BOOT_ITERS + 1) * (1.0 - alpha / 2.0);

#if 1
    fprintf(stderr, "IRF bootstrap, 0.025 and 0.975 quantiles\n");
    fprintf(stderr, " based on %d iterations, ilo = %d, ihi = %d\n", 
	    BOOT_ITERS, ilo, ihi);
#endif

    for (k=0; k<boot->horizon; k++) {
	gretl_matrix_row_to_array(boot->resp, k, rk);
	qsort(rk, BOOT_ITERS, sizeof *rk, gretl_compare_doubles);
#if 1
	fprintf(stderr, "%10g, %10g\n", rk[ilo-1], rk[ihi-1]);
#endif
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
    gretl_matrix *R;
    irfboot boot;
    int i, err = 0;

    R = gretl_matrix_alloc(periods, 3);
    if (R == NULL) {
	return NULL;
    }

    err = irf_boot_init(&boot, var, periods);
    if (err) {
	goto bailout;
    }

    err = make_bs_dataset_and_lists(&boot, var, Z, pdinfo);

    for (i=0; i<BOOT_ITERS && !err; i++) {
	resample_resids(&boot, var);
	compute_bootstrap_dataset(&boot, var);
	err = re_estimate_VAR(&boot, targ, shock, i);
    }

    if (!err) {
	irf_boot_quantiles(&boot, R);
    }

    irf_boot_free(&boot);

 bailout:

    if (err) {
	gretl_matrix_free(R);
	R = NULL;
    }

    return R;
}


