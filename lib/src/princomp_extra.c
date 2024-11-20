/*
  Note 2024-11-18: This file is (hopefully) a temporary hack, an
  intermediate step on the way to a proper renovation of gretl's
  princomp() function.
*/

static int matrix_is_incomplete (const gretl_matrix *m)
{
    int i, n = m->rows * m->cols;

    for (i=0; i<n; i++) {
        if (na(m->val[i])) {
            return 1;
        }
    }

    return 0;
}

static gint8 *makemask (const gretl_matrix *X, int *nfull, int *ndrop)
{
    gint8 *ret;
    double xij;
    int i, j;
    int nna;

    ret = malloc(X->rows * sizeof *ret);

    for (i=0; i<X->rows; i++) {
        nna = 0;
        for (j=0; j<X->cols; j++) {
            xij = gretl_matrix_get(X, i, j);
            if (na(xij)) {
                nna++;
            }
        }
        if (nna == 0) {
            ret[i] = 1;
            *nfull += 1;
        } else if (nna == X->cols) {
            ret[i] = -1;
            *ndrop += 1;
          } else {
            ret[i] = 0;
        }
    }

    return ret;
}

static gretl_matrix *select_rows_from (const gretl_matrix *M,
                                       int nr,
                                       const gint8 *mask)
{
    gretl_matrix *X = gretl_matrix_alloc(nr, M->cols);
    int i, j, k = 0;
    double mij;

    for (i=0; i<M->rows; i++) {
        if (mask[i] > 0) {
            for (j=0; j<M->cols; j++) {
                mij = gretl_matrix_get(M, i, j);
                gretl_matrix_set(X, k, j, mij);
            }
            k++;
        }
    }

    return X;
}

static void F_insert_P (gretl_matrix *F,
                        const gretl_matrix *P,
                        const gint8 *mask)
{
    double pkj;
    int i, j, k = 0;

    for (i=0; i<F->rows; i++) {
        if (mask[i] > 0) {
            for (j=0; j<F->cols; j++) {
                pkj = gretl_matrix_get(P, k, j);
                gretl_matrix_set(F, i, j, pkj);
            }
            k++;
        }
    }
}

/* Copy from @a to @targ unless an @a element is missing,
   in which case copy that element from @b instead.
*/

static void alt_transcribe_values (gretl_matrix *targ,
                                   const gretl_matrix *a,
                                   const gretl_matrix *b)
{
    int i, n = targ->rows * targ->cols;

    for (i=0; i<n; i++) {
        targ->val[i] = na(a->val[i]) ? b->val[i] : a->val[i];
    }
}

static double max_abs_value (const gretl_matrix *E)
{
    int i, n = E->rows * E->cols;
    double ai, ret = 0;

    for (i=0; i<n; i++) {
        ai = fabs(E->val[i]);
        if (ai > ret) {
            ret = ai;
        }
    }

    return ret;
}

static int demean_with_NAs (gretl_matrix *X, int stdize)
{
    double *Xj = X->val;
    double cmean, csd;
    int i, j, Tj;
    int err = 0;

    for (j=0; j<X->cols; j++) {
        cmean = csd = 0.0;
        Tj = 0;
        for (i=0; i<X->rows; i++) {
            if (!na(Xj[i])) {
                cmean += Xj[i];
                Tj++;
            }
        }
        if (Tj == 0) {
            gretl_errmsg_sprintf("Column %d is entirely missing", j+1);
            err = E_MISSDATA;
            break;
        }
        cmean /= Tj;
        for (i=0; i<X->rows; i++) {
            if (!na(Xj[i])) {
                Xj[i] -= cmean;
                if (stdize) {
                    csd += Xj[i] * Xj[i];
                }
            }
        }
        if (stdize) {
            csd = sqrt(csd / Tj);
            for (i=0; i<X->rows; i++) {
                if (!na(Xj[i])) {
                    Xj[i] /= csd;
                }
            }
        }
        Xj += X->rows;
    }

    return err;
}

static gretl_matrix *princomp_with_NAs (const gretl_matrix *M,
                                        int p, gretlopt opt,
                                        int *err)
{
    gretl_matrix *X = NULL;
    gretl_matrix *Tmp = NULL;
    gretl_matrix *F = NULL;
    gretl_matrix *L = NULL;
    gretl_matrix *FL = NULL;
    gretl_matrix *F0 = NULL;
    gretl_matrix *E = NULL;
    gretlopt pca_opt = 0;
    double tol = libset_get_user_tolerance(BHHH_TOLER);
    double crit;
    gint8 *mask = NULL;
    int maxiter = 2048;
    int stdize = 1;
    int conv = 0;
    int nfull = 0;
    int ndrop = 0;
    int iter = 0;

    if (opt & OPT_C) {
        stdize = 0;
        pca_opt |= OPT_V; /* FIXME? */
    }

    mask = makemask(M, &nfull, &ndrop);

    if (ndrop > 0) {
        *err = E_MISSDATA;
        gretl_errmsg_sprintf("%d rows were entirely missing", ndrop);
    } else if (nfull < p) {
        *err = E_TOOFEW;
        gretl_errmsg_sprintf("%d complete rows are needed but only %d were found",
                             p, nfull);
    } else {
        X = gretl_matrix_copy(M);
        if (X == NULL) {
            *err = E_ALLOC;
        }
    }

    if (*err) {
        free(mask);
        return NULL;
    }

    /* demean, and possibly standardize */
    *err = demean_with_NAs(X, stdize);

    if (!*err) {
        /* initialize loadings via PCA on complete observations */
        gretl_matrix *P;

        Tmp = select_rows_from(X, nfull, mask);
        F = gretl_zero_matrix_new(X->rows, p);
        P = real_gretl_matrix_pca(Tmp, p, pca_opt, err);
        F_insert_P(F, P, mask);
        L = gretl_zero_matrix_new(P->cols, Tmp->cols);
        *err = gretl_matrix_multi_ols(Tmp, P, L, NULL, NULL);
        gretl_matrix_free(P);
        gretl_matrix_free(Tmp);
        Tmp = NULL;
    }

    if (!*err) {
        /* additional matrix allocation */
        Tmp = gretl_matrix_alloc(X->rows, X->cols);
        F0 = gretl_matrix_alloc(F->rows, F->cols);
        FL = gretl_matrix_alloc(X->rows, X->cols);
        E = gretl_matrix_alloc(X->rows, F->cols);
    }

    /* EM iteration */
    while (!conv && !*err && iter < maxiter) {
        gretl_matrix_copy_values(F0, F);
        gretl_matrix_multiply(F, L, FL);
        alt_transcribe_values(Tmp, X, FL);
        gretl_matrix_free(F);
        F = real_gretl_matrix_pca(Tmp, p, pca_opt, err);
        if (!*err) {
            *err = gretl_matrix_multi_ols(F, F0, NULL, E, NULL);
        }
        if (*err) {
            break;
        }
        crit = max_abs_value(E);
#if 0
	  fprintf(stderr, "iter %04d: crit = %g, tol = %g\n", iter, crit, tol);
#endif
        if (crit < tol) {
            conv = 1;
        }
        *err = gretl_matrix_multi_ols(Tmp, F, L, NULL, NULL);
        iter++;
    }

    if (!*err && !conv) {
	fprintf(stderr, "PC EM algo didn't converge: crit = %g\n", crit);
        *err = E_NOCONV;
    }
    if (*err && F != NULL) {
        gretl_matrix_free(F);
        F = NULL;
    }

    gretl_matrix_free(X);
    gretl_matrix_free(Tmp);
    gretl_matrix_free(FL);
    gretl_matrix_free(F0);
    gretl_matrix_free(E);
    free(mask);

    return F;
}
