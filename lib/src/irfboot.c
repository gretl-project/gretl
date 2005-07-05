
/* bootstrapped confidence intervals for impulse response functions */

typedef struct irfboot_ irfboot;

struct irfboot_ {
    int n;
    int neqns;
    int order;
    int ifc;
    double **Z;
    DATAINFO *dinfo;
    int **lists;
    gretl_matrix *C;
    gretl_matrix *A;
    gretl_matrix *E;
    gretl_matrix *rE;
    gretl_matrix *rtmp;
    gretl_matrix *ctmp;
    double *resp;
    int *sample;
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
    gretl_matrix_free(boot->rE);
    gretl_matrix_free(boot->rtmp);
    gretl_matrix_free(boot->ctmp);

    free(boot->resp);
    free(boot->sample);
}

static int boot_A_C_init (irfboot *boot)
{
    int rows = boot->neqns * boot->order;
    int err = 0;

    boot->A = gretl_matrix_alloc(rows, rows);
    if (boot->A == NULL) {
	err = 1;
    } else {
	pad_var_coeff_matrix(boot->A, boot->neqns, boot->order);
    }

    if (!err) {
	boot->C = gretl_matrix_alloc(rows, boot->neqns);
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

static int irf_boot_init (irfboot *boot, const GRETL_VAR *var)
{
    int err = 0;

    boot->Z = NULL;
    boot->dinfo = NULL;
    boot->lists = NULL;

    boot->C = NULL;
    boot->A = NULL;
    boot->E = NULL;
    boot->rE = NULL;
    boot->rtmp = NULL;
    boot->ctmp = NULL;

    boot->resp = NULL;
    boot->sample = NULL;

    boot->n = var->n;
    boot->neqns = var->neqns;
    boot->order = var->order;
    boot->ifc = var->ifc;

    err = boot_A_C_init(boot);
    if (err) {
	goto bailout;
    }

    err = boot_tmp_init(boot);
    if (err) {
	goto bailout;
    }

    boot->rE = gretl_matrix_alloc(boot->n, boot->neqns);
    if (boot->E == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    boot->E = gretl_matrix_alloc(boot->n, boot->neqns);
    if (boot->E == NULL) {
	err = E_ALLOC;
	goto bailout;
    }

    boot->sample = malloc(var->n * sizeof *boot->sample);
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

int
recalculate_impulse_responses (int targ, int shock, irfboot *boot)
{
    int horizon = default_VAR_horizon(boot->dinfo);
    int t, err = 0;

    for (t=0; t<horizon && !err; t++) {
	if (t == 0) {
	    /* calculate initial estimated responses */
	    err = gretl_matrix_copy_values(boot->rtmp, boot->C);
	} else {
	    /* calculate further estimated responses */
	    err = gretl_matrix_multiply(boot->A, boot->rtmp, boot->ctmp);
	    gretl_matrix_copy_values(boot->rtmp, boot->ctmp);
	}

	if (!err) {
	    boot->resp[t] = gretl_matrix_get(boot->rtmp, targ, shock);
	}
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

static int re_estimate_VAR (irfboot *boot)
{
    MODEL var_model;
    int i, err = 0;

    for (i=0; i<boot->neqns && !err; i++) {
	var_model = lsq(boot->lists[i], &boot->Z, boot->dinfo, VAR, OPT_A, 0.0);
	err = var_model.errcode;

	if (!err) {
	    irf_record_model_data(boot, &var_model, i);
	}

	clear_model(&var_model);
    }

    if (!err) {
	err = gretl_var_do_error_decomp(boot->n, boot->neqns, boot->E, boot->C);
    }

#if 0
    if (!err) {
	for (i=0; i<boot->neqns; i++) {
	    err = real_var_get_impulse_responses(targ, shock, periods,
						 boot->A, boot->C,
						 boot->rtmp, boot->ctmp,
						 resp);
	}
    }
#endif

    return err;
}

static int allocate_bootstrap_lists (irfboot *boot, int nl)
{
    int i, err = 0;

    boot->lists = malloc(boot->neqns * sizeof *boot->lists);
    if (boot->lists == NULL) {
	err = E_ALLOC;
    }

    for (i=0; i<boot->neqns && !err; i++) {
	boot->lists[i] = gretl_list_new(nl);
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

static int make_bs_dataset_and_lists (const GRETL_VAR *var,
				      const double **Z,
				      const DATAINFO *pdinfo,
				      irfboot *boot)
{
    const MODEL *pmod;

    double **bZ = NULL;
    DATAINFO *binfo = NULL;
    int ns, nd, nv;
    int i, j, k, p;
    int m, v, t;
    int err = 0;

    pmod = var->models[0]; 

    m = pmod->list[0] - 1;        /* total regressors per equation */
    ns = var->order * var->neqns; /* stochastic regressors per equation */
    nd = m - ns;                  /* number of deterministic regressors */
    nv = m + var->neqns;          /* total variables required */

    err = allocate_bootstrap_lists(boot, m + 1);
    
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
	for (i=0; i<var->neqns; i++) {
	    p = 1;
	    /* copy stochastic vars into positions 1 to var->neqns */
	    pmod = var->models[i];
	    v = pmod->list[1];
	    k = i + 1;
	    for (t=0; t<pdinfo->n; t++) {
		bZ[k][t] = Z[v][t];
	    }
	    boot->lists[i][p++] = k;
	    /* create lags */
	    for (j=1; j<=var->order; j++) {
		k = var->neqns + (i + 1) * j;
		for (t=0; t<pdinfo->n; t++) {
		    bZ[k][t] = (t - j >= 0)? Z[v][t-j] : NADBL;
		}
		boot->lists[i][p++] = k;
	    }
	    /* deal with deterministic vars */
	    k = var->neqns * (var->order + 1);
	    for (j=ns+1; j<=pmod->list[0]; j++) {
		v = pmod->list[j];
		if (v == 0) {
		    boot->lists[i][p++] = 0;
		    continue;
		}
		k++;
		boot->lists[i][p++] = k;
		if (i == 0) {
		    for (t=0; t<pdinfo->n; t++) {
			bZ[k][t] = Z[v][t];
		    }
		}
	    }
	}
    }

    return err;
}

static int fill_bootstrap_dataset (irfboot *boot, const GRETL_VAR *var)
{
    const MODEL *pmod;

    double bti, xti;
    int i, j, k, s, t;
    int nc, t1, t2;
    int ns, lag, vj, m;
    int err = 0;

    ns = boot->order * boot->neqns;

    pmod = var->models[0];
    nc = pmod->list[0] - 1;

    t1 = pmod->t1;
    t2 = pmod->t2;

    for (t=t1; t<=t2; t++) {
	int miss = 0;

	s = t - t1;
	for (i=0; i<boot->neqns; i++) {
	    pmod = var->models[i];
	    bti = 0.0;
	    lag = 1;
	    k = 0;
	    for (j=0; j<nc; j++) {
		vj = boot->lists[i][j+2];
		if (j < ns + pmod->ifc && vj > 0) {
		    /* stochastic var */
		    if (s - lag < 0) { /* FIXME? */
			/* pre-simulation value */
			m = (j - boot->ifc) / boot->order;
			vj = boot->lists[m][1];
			if (t - lag < 0) {
			    xti = NADBL;
			} else {
			    xti = boot->Z[vj][t-lag];
			}
			if (na(xti)) {
			    miss = 1;
			}
		    } else {
			/* prior simulated value */
			xti = boot->Z[k+1][t-lag]; /* ?? */
			/* xti = gretl_matrix_get(F, s - lag, k); */
		    }
		    lag++;
		    if (lag > boot->order) {
			lag = 1;
			k++;
		    }
		} else {
		    /* deterministic var */
		    xti = boot->Z[vj][t];
		    if (na(xti)) {
			miss = 1;
		    }
		}
		if (miss) {
		    bti = NADBL;
		} else {
		    bti += pmod->coeff[j] * xti;
		}
	    }
	    /* FIXME need to fill out lags too */
	    boot->Z[i+1][t] = bti + gretl_matrix_get(boot->rE, s, i);
	}
    }

    return err;
}

static void resample_resids (irfboot *boot, const GRETL_VAR *var)
{
    double e;
    int i, s;

    for (s=0; s<boot->n; s++) {
	boot->sample[s] = gretl_rand_int_max(boot->n);
    }

    for (s=0; s<boot->n; s++) {
	for (i=0; i<boot->neqns; i++) {
	    e = gretl_matrix_get(var->E, boot->sample[s], i);
	    gretl_matrix_set(boot->rE, s, i, e);
	}
    }
}

#define ITERS 499

int irf_bootstrap_driver (const GRETL_VAR *var, 
			  const double **Z, 
			  DATAINFO *pdinfo)
{
    irfboot boot;
    int i, err = 0;

    err = irf_boot_init(&boot, var);

    if (!err) {
	err = make_bs_dataset_and_lists(var, Z, pdinfo, &boot);
    }

    for (i=0; i<ITERS && !err; i++) {
	resample_resids(&boot, var);
	fill_bootstrap_dataset(&boot, var);
	err = re_estimate_VAR(&boot);
    }

    irf_boot_free(&boot);

    return err;
}


