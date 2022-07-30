/* Write the square root of diagonal of square matrix @src
   into row @t of @targ, starting at column offset @j
*/

static void record_to_diag (gretl_matrix *targ,
                            const gretl_matrix *src,
                            int t, int j)
{
    int i, n = gretl_vector_get_length(src);
    double x;

    for (i=0; i<n; i++) {
        x = gretl_matrix_get(src, i, i);
        if (x <= 0.0) {
            gretl_matrix_set(targ, t, i+j, 0.0);
        } else {
            gretl_matrix_set(targ, t, i+j, sqrt(x));
        }
    }
}

enum {
    SM_TYPE_NONE, /* not doing smoothing */
    SM_STATE_STD, /* regular state smoother */
    SM_DIST_BKWD, /* disturbance smoother, backward pass */
    SM_DIST_FRWD, /* disturbance smoother, forward pass */
};

static int ksmooth_refresh_matrices (kalman *K, PRN *prn)
{
    gretl_matrix **mptr[] = {
        &K->VS, &K->VY
    };
    int idx[] = {
        K_VS, K_VY
    };
    int cross_update = 0;
    int i, j, err = 0;

    if (kalman_xcorr(K)) {
        mptr[0] = &K->H;
        mptr[1] = &K->G;
    }

    if (K->matcall != NULL) {
        err = kalman_update_matrices(K, prn);
    }

    for (i=0; i<2 && !err; i++) {
        j = idx[i];
        if (matrix_is_varying(K, j)) {
            if (kalman_xcorr(K) && (j == K_VS || j == K_VY)) {
                /* handle revised H and/or G */
                cross_update = 1;
            } else {
                err = check_matrix_dims(K, *mptr[i], j);
            }
            if (err) {
                fprintf(stderr, "ksmooth_refresh_matrices: err = %d at t = %d\n",
                        err, K->t);
            }
        }
    }

    if (!err && cross_update) {
        /* cross-correlated case */
        err = kalman_update_crossinfo(K, UPDATE_STEP);
    }

    return err;
}

/* For use with smoothing: load what we need for step @t, from the
   record that was kept on the prior forecasting pass. What we need
   from the forward pass depends on what exactly we're computing,
   which is conveyed by @smtype. If @pnt is non-NULL, use it to pass
   back the effective size of the observables vector.

   We need to handle some complications here. In particular we must
   update anything that's time-varying, and (unless we're doing
   disturbance smoothing) deal with any incomplete observations, which
   have implications for the dimensions of various matrices.
*/

static int old_load_filter_data (kalman *K, int smtype)
{
    int err = 0;

    /* load the forecast error */
    load_from_row(K->vt, K->V, K->t);

    if (smtype == SM_DIST_BKWD) {
        /* disturbances */
        if (filter_is_varying(K)) {
            ksmooth_refresh_matrices(K, NULL);
        }
    } else {
        /* state: get T_t and/or Z_t if need be */
        if (matrix_is_varying(K, K_T)) {
            err = retrieve_Tt(K);
        }
        if (!err && matrix_is_varying(K, K_ZT)) {
            err = retrieve_Zt(K);
        }
        if (err) {
            return err;
        }
    }

    if (smtype < SM_DIST_BKWD) {
        /* load the state and its MSE */
        load_from_row(K->a0, K->A, K->t);
        load_from_vech(K->P0, K->P, K->r, K->t, GRETL_MOD_NONE);
    }

    if (smtype == SM_STATE_STD || smtype == SM_DIST_BKWD) {
        /* load the gain and F^{-1} */
        load_from_vec(K->Kt, K->K, K->t);
        load_from_vech(K->iFt, K->F, K->n, K->t, GRETL_MOD_NONE);
    }

    return err;
}

/* For disturbance smoothing: ensure we have on hand matrices that are
   correctly sized to hold estimates of the variance of the
   disturbance(s) in the state and (if applicable) observation
   equations.
*/

static int maybe_resize_dist_mse (kalman *K,
                                  gretl_matrix **Vwt,
                                  gretl_matrix **Vut)
{
    int n = K->VY == NULL ? 0 : K->n;
    int k, err = 0;

    /* combined results: how many columns do we need? */
    k = K->r + n;

    if (K->Vsd == NULL) {
        K->Vsd = gretl_matrix_alloc(K->N, k);
        if (K->Vsd == NULL) {
            err = E_ALLOC;
        }
    } else if (K->Vsd->rows != K->N || K->Vsd->cols != k) {
        err = gretl_matrix_realloc(K->Vsd, K->N, k);
    }

    if (!err) {
        /* step-t square state matrix */
        *Vwt = gretl_matrix_alloc(K->r, K->r);
        if (*Vwt == NULL) {
            err = E_ALLOC;
        }
    }

    if (!err && n > 0) {
        /* step-t square obs matrix */
        *Vut = gretl_matrix_alloc(K->n, K->n);
        if (*Vut == NULL) {
            err = E_ALLOC;
        }
    }

    return err;
}

/* Calculate the variance of the smoothed disturbances for the
   cross-correlated case. See Koopman, Shephard and Doornik (1998),
   page 19, var(\varepsilon_t|Y_n).
*/

static int combined_dist_variance (kalman *K,
                                   gretl_matrix *D,
                                   gretl_matrix *Nt,
                                   gretl_matrix *Vwt,
                                   gretl_matrix *Vut,
                                   gretl_matrix_block *BX,
                                   int DKstyle)
{
    gretl_matrix *DG, *KN, *Veps, *NH, *NK;

    DG   = gretl_matrix_block_get_matrix(BX, 0);
    KN   = gretl_matrix_block_get_matrix(BX, 1);
    Veps = gretl_matrix_block_get_matrix(BX, 2);
    NH   = gretl_matrix_block_get_matrix(BX, 3);

    /* First chunk of Veps in Koopman's notation:
       G_t'(D_t*G_t - K_t'*N_t*H_t)
    */
    KN = gretl_matrix_reuse(KN, K->n, K->r);
    gretl_matrix_multiply(D, K->G, DG);
    gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_TRANSPOSE,
                              Nt, GRETL_MOD_NONE,
                              KN, GRETL_MOD_NONE);
    gretl_matrix_multiply_mod(KN, GRETL_MOD_NONE,
                              K->H, GRETL_MOD_NONE,
                              DG, GRETL_MOD_DECREMENT);
    gretl_matrix_multiply_mod(K->G, GRETL_MOD_TRANSPOSE,
                              DG, GRETL_MOD_NONE,
                              Veps, GRETL_MOD_NONE);

    /* Second chunk of Veps, to be added to the above
       H_t'(N_t*H_t - N_t*K_t*G_t)
    */
    NK = gretl_matrix_reuse(KN, K->r, K->n);
    gretl_matrix_multiply(Nt, K->H, NH);
    gretl_matrix_multiply_mod(Nt, GRETL_MOD_TRANSPOSE,
                              K->Kt, GRETL_MOD_NONE,
                              NK, GRETL_MOD_NONE);
    gretl_matrix_multiply_mod(NK, GRETL_MOD_NONE,
                              K->G, GRETL_MOD_NONE,
                              NH, GRETL_MOD_DECREMENT);
    gretl_matrix_multiply_mod(K->H, GRETL_MOD_TRANSPOSE,
                              NH, GRETL_MOD_NONE,
                              Veps, GRETL_MOD_CUMULATE);

    if (DKstyle) {
        /* Veps = I_p - Veps */
        double vii;
        int i;

        gretl_matrix_multiply_by_scalar(Veps, -1.0);
        for (i=0; i<K->p; i++) {
            vii = gretl_matrix_get(Veps, i, i);
            gretl_matrix_set(Veps, i, i, 1.0 + vii);
        }
    }

    /* Veps (p x p) holds the variance of \epsilon_t
       conditional on Y_n: now form the per-equation
       disturbance variance matrices, @Vwt and @Vut,
       for this time-step.
    */
    gretl_matrix_qform(K->H, GRETL_MOD_NONE, Veps,
                       Vwt, GRETL_MOD_NONE);
    gretl_matrix_qform(K->G, GRETL_MOD_NONE, Veps,
                       Vut, GRETL_MOD_NONE);

    return 0;
}

static int dist_variance (kalman *K,
			  gretl_matrix *D,
			  gretl_matrix *Vwt,
			  gretl_matrix *Vut,
			  gretl_matrix *Nt,
			  gretl_matrix_block *BX,
			  int DKstyle)
{
    int err = 0;

    if (D != NULL) {
	/* needed only in presence of obs disturbance */
	/* D_t = F_t^{-1} + K_t' N_t K_t */
	fast_copy_values(D, K->iFt);
	if (K->t < K->N - 1) {
	    gretl_matrix_qform(K->Kt, GRETL_MOD_TRANSPOSE,
			       Nt, D, GRETL_MOD_CUMULATE);
	}
    }

    if (K->p == 0) {
	/* variance of state disturbance */
	if (DKstyle) {
	    /* HH' - HH' N_t HH' */
	    fast_copy_values(Vwt, K->VS);
	    gretl_matrix_qform(K->VS, GRETL_MOD_TRANSPOSE,
			       Nt, Vwt, GRETL_MOD_DECREMENT);
	} else {
	    /* HH' N_t HH' */
	    gretl_matrix_qform(K->VS, GRETL_MOD_TRANSPOSE,
			       Nt, Vwt, GRETL_MOD_NONE);
	}
	record_to_diag(K->Vsd, Vwt, K->t, 0);

        /* variance of obs disturbance */
        if (DKstyle) {
            /* GG' - GG D_t GG' */
            fast_copy_values(Vut, K->VY);
            gretl_matrix_qform(K->VY, GRETL_MOD_TRANSPOSE,
                               D, Vut, GRETL_MOD_DECREMENT);
        } else {
            /* GG' D_t GG' */
            gretl_matrix_qform(K->VY, GRETL_MOD_TRANSPOSE,
                               D, Vut, GRETL_MOD_NONE);
        }
        record_to_diag(K->Vsd, Vut, K->t, K->r);
    } else {
        /* cross-correlated disturbance variance */
        err = combined_dist_variance(K, D, Nt, Vwt, Vut, BX,
                                     DKstyle);
        if (!err) {
            record_to_diag(K->Vsd, Vwt, K->t, 0);
            record_to_diag(K->Vsd, Vut, K->t, K->r);
        }
    }

    return err;
}

/* This iteration is in common between the state smoother
   (Anderson-Moore) and the disturbance smoother (Koopman).
*/

static void LrN_iteration (kalman *K,
			   gretl_matrix *L,
			   gretl_matrix *n1,
			   gretl_matrix *r0,
			   gretl_matrix *r1,
			   gretl_matrix *N0,
			   gretl_matrix *N1)
{
    if (K->t < K->N - 1) {
	/* L_t = T_t - K_t Z_t */
	fast_copy_values(L, K->T);
	gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
				  K->ZT, GRETL_MOD_TRANSPOSE,
				  L, GRETL_MOD_DECREMENT);
    }

    /* r_{t-1} = Z_t' F_t^{-1} v_t + L_t' r_t */
    gretl_matrix_multiply(K->iFt, K->vt, n1);
    if (K->t == K->N - 1) {
	gretl_matrix_multiply(K->ZT, n1, r0);
    } else {
	gretl_matrix_multiply(K->ZT, n1, r1);
	gretl_matrix_multiply_mod(L, GRETL_MOD_TRANSPOSE,
				  r0, GRETL_MOD_NONE,
				  r1, GRETL_MOD_CUMULATE);
	fast_copy_values(r0, r1);
    }

    /* N_{t-1} = Z_t' F_t^{-1} Z_t + L_t' N_t L_t */
    if (K->t == K->N - 1) {
	gretl_matrix_qform(K->ZT, GRETL_MOD_NONE,
			   K->iFt, N0, GRETL_MOD_NONE);
    } else {
	gretl_matrix_qform(K->ZT, GRETL_MOD_NONE,
			   K->iFt, N1, GRETL_MOD_NONE);
	gretl_matrix_qform(L, GRETL_MOD_TRANSPOSE,
			   N0, N1, GRETL_MOD_CUMULATE);
	fast_copy_values(N0, N1);
    }
}

/* Initial smoothed state: a + P*r0 */

static void koopman_calc_a0 (kalman *K, gretl_matrix *r0)
{
    if (K->Pini != NULL) {
        gretl_matrix_multiply(K->Pini, r0, K->a0);
    } else {
        set_initial_statevar(K);
        gretl_matrix_multiply(K->P0, r0, K->a0);
    }
    if (K->aini != NULL) {
        gretl_matrix_add_to(K->a0, K->aini);
    }
    record_to_row(K->A, K->a0, 0);
}

/* Disturbance smoothing -- see Koopman, Shephard and Doornik
   (SsfPack doc), section 4.4.
*/

static int koopman_smooth (kalman *K, int DKstyle)
{
    gretl_matrix_block *B, *BX = NULL;
    gretl_matrix *u, *L, *R;
    gretl_matrix *r0, *r1, *N0, *N1, *n1, *tr;
    gretl_matrix *D = NULL;
    gretl_matrix *Vwt = NULL;
    gretl_matrix *Vut = NULL;
    gretl_matrix *DG = NULL;
    gretl_matrix *KN = NULL;
    gretl_matrix *RZS = NULL;
    gretl_matrix *NH = NULL;
    gretl_matrix *Ut = NULL;
    int ft_min = 0;
    int t, err = 0;

    if (trace) {
        printf("koopman_smooth(), legacy, DKstyle = %d\n", DKstyle);
    }

    B = gretl_matrix_block_new(&u,  K->n, 1,
                               &L,  K->r, K->r,
                               &R,  K->N, K->r,
                               &r0, K->r, 1,
                               &r1, K->r, 1,
                               &N0, K->r, K->r,
                               &N1, K->r, K->r,
                               &n1, K->n, 1,
                               &tr, K->r, 1,
                               NULL);

    if (B == NULL) {
        return E_ALLOC;
    }

    /* for variance of smoothed disturbances */
    err = maybe_resize_dist_mse(K, &Vwt, &Vut);

    if (K->VY != NULL) {
	/* for variance of observable */
	D = gretl_matrix_alloc(K->n, K->n);
    }

    if (K->p > 0) {
	/* cross-correlated disturbances */
        BX = gretl_matrix_block_new(&DG,  K->n, K->p,
                                    &KN,  K->n, K->r,
                                    &RZS, K->p, K->p,
                                    &NH,  K->r, K->p,
                                    &Ut,  K->p, 1,
                                    NULL);
        if (BX == NULL) {
            err = E_ALLOC;
        }
    }

    if (err) {
        gretl_matrix_block_destroy(B);
        gretl_matrix_block_destroy(BX);
        gretl_matrix_free(Vwt);
        gretl_matrix_free(Vut);
	gretl_matrix_free(D);
        return err;
    }

    gretl_matrix_zero(r0);
    gretl_matrix_zero(N0);

    /* The backward recursion */

    for (t=K->N-1; t>=0 && !err; t--) {
	K->t = t;
	err = old_load_filter_data(K, SM_DIST_BKWD);
        if (err) {
            break;
        }

        /* u_t = F_t^{-1} v_t - K_t' r_t */
        gretl_matrix_multiply(K->iFt, K->vt, u);
        if (t < K->N - 1) {
            gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_TRANSPOSE,
                                      r0, GRETL_MOD_NONE,
                                      u, GRETL_MOD_DECREMENT);
        }
        /* Store u_t values in K->V: these are needed in
           the forward pass to compute the smoothed
           disturbances. Also save r_t in R.
        */
        record_to_row(K->V, u, t);
	record_to_row(R, r0, t);

	/* compute variance of disturbances */
	err = dist_variance(K, D, Vwt, Vut, N0, BX, DKstyle);
	if (err) {
	    break;
	}

	/* compute r_{t-1}, N_{t-1} */
	LrN_iteration(K, L, n1, r0, r1, N0, N1);

	if (t == 0) {
	    /* compute initial smoothed state */
	    koopman_calc_a0(K, r0);
        }
    }

    /* Forward iteration for smoothed disturbances, all time steps,
       plus smoothed state from t = 1 onward.
    */
    for (t=ft_min; t<K->N; t++) {
	K->t = t;
        err = old_load_filter_data(K, SM_DIST_FRWD);
        if (err) {
            break;
        }

	/* state disturbance */
        load_from_row(r0, R, t);
        if (K->p > 0) {
            gretl_matrix_multiply_mod(K->H, GRETL_MOD_TRANSPOSE,
                                      r0, GRETL_MOD_NONE,
                                      Ut, GRETL_MOD_NONE);
            load_from_row(K->vt, K->V, t);
            gretl_matrix_multiply_mod(K->G, GRETL_MOD_TRANSPOSE,
                                      K->vt, GRETL_MOD_NONE,
                                      Ut, GRETL_MOD_CUMULATE);
            gretl_matrix_multiply(K->H, Ut, r1);
        } else {
            gretl_matrix_multiply(K->VS, r0, r1);
        }
        record_to_row(R, r1, t);
	record_to_row_offset(K->U, r1, t, 0);

	/* observation disturbance */
	if (K->VY != NULL) {
	    if (K->p > 0) {
		gretl_matrix_multiply(K->G, Ut, n1);
	    } else {
		gretl_matrix_multiply(K->VY, K->vt, n1);
	    }
	    record_to_row_offset(K->U, n1, t, K->r);
	}

	/* state: a_{t+1} = T a_t + w_t (or + H*eps_t) */
	load_from_row(K->a0, K->A, t-1);
	gretl_matrix_multiply(K->T, K->a0, K->a1);
	vector_from_row_mod(K->a1, R, t-1, GRETL_MOD_CUMULATE);
	if (K->mu != NULL) {
	    gretl_matrix_add_to(K->a1, K->mu);
	}
	record_to_row(K->A, K->a1, t);
    }

    gretl_matrix_block_destroy(B);
    gretl_matrix_block_destroy(BX);
    gretl_matrix_free(Vwt);
    gretl_matrix_free(Vut);
    gretl_matrix_free(D);

    return err;
}

/* Anderson-Moore state smoothing */

static int anderson_moore_smooth (kalman *K)
{
    gretl_matrix_block *B;
    gretl_matrix *r0, *r1, *N0, *N1, *n1, *L;
    int t, err = 0;

    if (trace) {
        printf("anderson_moore_smooth(), legacy\n");
    }

    B = gretl_matrix_block_new(&r0,  K->r, 1,
                               &r1,  K->r, 1,
                               &N0,  K->r, K->r,
                               &N1,  K->r, K->r,
                               &n1,  K->n, 1,
                               &L,   K->r, K->r,
                               NULL);
    if (B == NULL) {
        return E_ALLOC;
    }

    gretl_matrix_zero(r0);
    gretl_matrix_zero(N0);

    for (t=K->N-1; t>=0 && !err; t--) {
	K->t = t;
        err = old_load_filter_data(K, SM_STATE_STD);
        if (err) {
            break;
        }

	/* compute r_{t-1}, N_{t-1} */
	LrN_iteration(K, L, n1, r0, r1, N0, N1);

        /* a_{t|T} = a_{t|t-1} + P_{t|t-1} r_{t-1} */
        fast_copy_values(K->a1, K->a0);
        gretl_matrix_multiply_mod(K->P0, GRETL_MOD_NONE,
                                  r0, GRETL_MOD_NONE,
                                  K->a1, GRETL_MOD_CUMULATE);
        record_to_row(K->A, K->a1, t);

        /* P_{t|T} = P_{t|t-1} - P_{t|t-1} N_{t-1} P_{t|t-1} */
        fast_copy_values(K->P1, K->P0);
        gretl_matrix_qform(K->P0, GRETL_MOD_NONE, N0,
                           K->P1, GRETL_MOD_DECREMENT);
        record_to_vech(K->P, K->P1, K->r, t);
    }

    gretl_matrix_block_destroy(B);

    return err;
}
