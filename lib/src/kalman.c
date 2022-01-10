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

#define FULL_XML_HEADERS 1

#include "libgretl.h"
#include "uservar.h"
#include "gretl_func.h"
#include "gretl_xml.h"
#include "matrix_extra.h"
#include "libset.h"
#include "kalman.h"

/**
 * SECTION:kalman
 * @short_description: The Kalman filter
 * @title: Kalman
 * @include: gretl/libgretl.h, gretl/kalman.h
 *
 */

#define KDEBUG 0

/* exact initial smoothing: work in progress */
#define EXACT_SM 0

/*
   State: a_{t+1} = T_t*a_t + u_t,     E(u_t*u_t') = Q

   Obs:   y_t = B*x_t + Z*a_t + w_t,   E(w_t*w_t') = R
*/

#define K_TINY 1.0e-15

typedef struct stepinfo_ stepinfo;

struct stepinfo_ {
    gretl_matrix *T;   /* N x (r * r) */
    gretl_matrix *ZT;  /* N x (r * n) */
};

struct kalman_ {
    int flags;  /* for recording any options */
    int exact;  /* exact initial iterations? */

    int r;   /* rows of a = number of elements in state */
    int n;   /* columns of y = number of observables */
    int k;   /* columns of B = number of exogenous vars in obs eqn */
    int p;   /* length of combined disturbance vector */
    int N;   /* rows of y = number of observations */
    int okN; /* N - number of missing observations */
    int t;   /* current time step, when filtering */
    int d;   /* time-step at which standard iterations start */

    int ifc; /* boolean: obs equation includes an implicit constant? */

    double SSRw;    /* \sum_{t=1}^N v_t^{\prime} F_t^{-1} v_t */
    double loglik;  /* log-likelihood */
    double s2;      /* = SSRw / k */

    /* continuously updated matrices */
    gretl_matrix *a0; /* r x 1: state vector, before updating */
    gretl_matrix *a1; /* r x 1: state vector, after updating */
    gretl_matrix *P0; /* r x r: MSE matrix, before updating */
    gretl_matrix *P1; /* r x r: MSE matrix, after updating */
    gretl_matrix *v;  /* n x 1: one-step forecast error(s), time t */

    /* for the exact diffuse case */
    gretl_matrix *Pk0;
    gretl_matrix *Pk1;
    gretl_matrix *Fk;
    gretl_matrix *Ck;
    gretl_matrix *PK; /* for smoothing, not actually used yet */

    /* input data matrices: note that the order matters for various
       functions, including matrix_is_varying()
    */
    gretl_matrix *T;  /* r x r: state transition matrix */
    gretl_matrix *BT; /* k x n: B', coeffs on exogenous vars, obs eqn */
    gretl_matrix *ZT; /* r x n: Z', coeffs on state variables, obs eqn */
    gretl_matrix *HH; /* r x r: contemp covariance matrix, state eqn */
    gretl_matrix *GG; /* n x n: contemp covariance matrix, obs eqn */
    gretl_matrix *mu; /* r x 1: constant term in state transition */
    gretl_matrix *y;  /* N x n: dependent variable vector (or matrix) */
    gretl_matrix *x;  /* N x k: independent variables matrix */
    gretl_matrix *aini; /* r x 1: a_0 */
    gretl_matrix *Pini; /* r x r: P_0 */

    /* user inputs for cross-correlated disturbances */
    gretl_matrix *H;  /* r x p: HH' ("Q") */
    gretl_matrix *G;  /* n x p: GG' ("R") */
    /* and the cross-matrix itself */
    gretl_matrix *HG; /* r x n */

    /* apparatus for registering time-variation of matrices */
    char *matcall;
    char *varying;

    /* optional matrices for recording extra info */
    gretl_matrix *LL;  /* N x 1: loglikelihood, all time-steps */

    /* optional run-time export matrices */
    gretl_matrix *V;   /* N x n: forecast errors, all time-steps */
    gretl_matrix *F;   /* N x nn: MSE for observables, all time-steps */
    gretl_matrix *A;   /* N x r: state vector, all time-steps */
    gretl_matrix *P;   /* N x nr: MSE for state, all time-steps */
    gretl_matrix *K;   /* N x rn: gain matrix, all time-steps */
    gretl_matrix *U;   /* N x ??: smoothed disturbances */
    gretl_matrix *Vsd; /* Variance of smoothed disturbance */

    /* structure needed only when smoothing in the time-varying case */
    stepinfo *step;

    /* workspace matrices */
    gretl_matrix_block *Blk; /* holder for the following */
    gretl_matrix *PZ;
    gretl_matrix *Ft;
    gretl_matrix *iFt;
    gretl_matrix *Bx;
    gretl_matrix *Kt;
    gretl_matrix *Mt;
    gretl_matrix *Ct;

    gretl_bundle *b; /* the bundle of which this struct is a member */
    void *data;      /* handle for attaching additional info */
    PRN *prn;        /* verbose printer */
};

/* max number of time-varying matrices: T, B, Z, Q, R, mu */
#define K_N_MATCALLS 6

#define set_kalman_running(K) (K->flags |= KALMAN_FORWARD)
#define set_kalman_stopped(K) (K->flags &= ~KALMAN_FORWARD)
#define kalman_is_running(K)  (K->flags & KALMAN_FORWARD)
#define kalman_simulating(K)  (K->flags & KALMAN_SIM)
#define kalman_checking(K)    (K->flags & KALMAN_CHECK)
#define kalman_xcorr(K)       (K->flags & KALMAN_CROSS)
#define kalman_ssfsim(K)      (K->flags & KALMAN_SSFSIM)
#define kalman_diffuse(K)     (K->flags & KALMAN_DIFFUSE)
#define kalman_smoothing(K)   (K->flags & KALMAN_SMOOTH)

#define filter_is_varying(K) (K->matcall != NULL)

static const char *kalman_matrix_name (int sym);
static int kalman_revise_variance (kalman *K);
static int check_for_matrix_updates (kalman *K, ufunc *uf);

/* symbolic identifiers for input matrices: note that potentially
   time-varying matrices must appear first in the enumeration, and
   the order must match the order of the "input data matrices" in
   the Kalman struct (above).
*/

enum {
    K_T = 0,
    K_BT,
    K_ZT,
    K_Q,
    K_R,
    K_m,
    K_y,
    K_x,
    K_a,
    K_P,
    K_MMAX /* sentinel */
};

void free_stepinfo (kalman *K)
{
    if (K->step != NULL) {
        gretl_matrix_free(K->step->T);
        gretl_matrix_free(K->step->ZT);
        free(K->step);
        K->step = NULL;
    }
}

void kalman_free (kalman *K)
{
    if (K == NULL) {
        return;
    }

    /* internally allocated matrices */
    gretl_matrix_free(K->a0);
    gretl_matrix_free(K->a1);
    gretl_matrix_free(K->P0);
    gretl_matrix_free(K->P1);
    gretl_matrix_free(K->v);
    gretl_matrix_free(K->LL);
    gretl_matrix_free(K->Pk0);
    gretl_matrix_free(K->Pk1);
    gretl_matrix_free(K->Fk);
    gretl_matrix_free(K->Ck);

    /* internally allocated workspace */
    gretl_matrix_block_destroy(K->Blk);

    if (K->flags & KALMAN_BUNDLE) {
        gretl_matrix **mptr[] = {
            &K->T, &K->BT, &K->ZT, &K->HH, &K->GG,
            &K->mu, &K->y, &K->x, &K->aini, &K->Pini
        };
        int i;

        for (i=0; i<K_MMAX; i++) {
            gretl_matrix_free(*mptr[i]);
        }

        /* @K also owns these "export" matrices */
        gretl_matrix_free(K->V);
        gretl_matrix_free(K->F);
        gretl_matrix_free(K->A);
        gretl_matrix_free(K->P);
        gretl_matrix_free(K->K);
        gretl_matrix_free(K->U);
        gretl_matrix_free(K->Vsd);
    }

    /* time-variation info */
    free(K->matcall);
    free(K->varying);

    if (kalman_xcorr(K)) {
        /* correlated errors info */
        gretl_matrix_free(K->H);
        gretl_matrix_free(K->G);
        gretl_matrix_free(K->HG);
    }

    if (K->step != NULL) {
        free_stepinfo(K);
    }

    free(K);
}

static kalman *kalman_new_empty (int flags)
{
    kalman *K = malloc(sizeof *K);

    if (K != NULL) {
	K->exact = 0;
        K->aini = K->Pini = NULL;
        K->a0 = K->a1 = NULL;
        K->P0 = K->P1 = NULL;
        K->Pk0 = K->Pk1 = NULL;
	K->PK = NULL;
        K->Fk = K->Ck = NULL;
        K->LL = NULL;
        K->v = NULL;
        K->Blk = NULL;
        K->T = K->BT = K->ZT = NULL;
        K->HH = K->GG = K->HG = NULL;
        K->H = K->G = NULL;
        K->V = K->F = K->A = K->K = K->P = NULL;
        K->y = K->x = NULL;
        K->mu = NULL;
        K->U = NULL;
        K->Vsd = NULL;
        K->matcall = NULL;
        K->varying = NULL;
        K->step = NULL;
        K->flags = flags;
        K->t = 0;
        K->prn = NULL;
        K->data = NULL;
        K->b = NULL;
        K->d = 0;
    }

    return K;
}

#define kappa 1.0e7 /* 1.0e6? */

static int diffuse_Pini (kalman *K)
{
    gretl_matrix_zero(K->P0);

    if (K->exact) {
	if (K->Pk0 == NULL) {
	    K->Pk0 = gretl_identity_matrix_new(K->r);
	    K->Pk1 = gretl_matrix_alloc(K->r, K->r);
	    K->Fk  = gretl_matrix_alloc(K->n, K->n);
	    K->Ck  = gretl_matrix_alloc(K->r, K->r);
	    if (K->Pk0 == NULL || K->Pk1 == NULL ||
		K->Fk == NULL || K->Ck == NULL) {
		return E_ALLOC;
	    }
	} else {
	    gretl_matrix_inscribe_I(K->Pk0, 0, 0, K->r);
	}
    } else {
	/* old-style */
	int i;

	for (i=0; i<K->r; i++) {
	    gretl_matrix_set(K->P0, i, i, kappa);
	}
    }

    return 0;
}

static int statemat_out_of_bounds (kalman *K)
{
    gretl_matrix *evals;
    double r, c, x;
    int i, err = 0;

    evals = gretl_general_matrix_eigenvals(K->T, &err);

    for (i=0; i<evals->rows && !err; i++) {
        r = gretl_matrix_get(evals, i, 0);
        c = gretl_matrix_get(evals, i, 1);
        x = sqrt(r*r + c*c);
        if (x >= 1.0) {
            fprintf(stderr, "T: modulus of eigenvalue %d = %g\n", i+1, x);
            err = E_SINGULAR;
        }
    }

    gretl_matrix_free(evals);

    return err;
}

/* If the user has not given an initial value for P_{1|0}, compute
   this automatically as per Hamilton, ch 13, p. 378.  This works only
   if the eigenvalues of K->T lie inside the unit circle.  Failing
   that, or if the --diffuse option is given for the user Kalman
   filter, we apply a diffuse initialization.
*/

static int construct_Pini (kalman *K)
{
    gretl_matrix *Svar;
    gretl_matrix *vQ;
    int r2, err = 0;

    if (K->flags & KALMAN_DIFFUSE) {
        return diffuse_Pini(K);
    }

    r2 = K->r * K->r;

    Svar = gretl_matrix_alloc(r2, r2);
    vQ = gretl_column_vector_alloc(r2);

    if (Svar == NULL || vQ == NULL) {
        gretl_matrix_free(Svar);
        gretl_matrix_free(vQ);
        return E_ALLOC;
    }

    gretl_matrix_kronecker_product(K->T, K->T, Svar);
    gretl_matrix_I_minus(Svar);
    gretl_matrix_vectorize(vQ, K->HH);

    err = gretl_LU_solve(Svar, vQ);
    if (err) {
        /* failed: are some of the eigenvalues out of bounds? */
        err = statemat_out_of_bounds(K);
        if (err == E_SINGULAR) {
            err = diffuse_Pini(K);
            K->flags |= KALMAN_DIFFUSE;
        }
    } else {
        gretl_matrix_unvectorize(K->P0, vQ);
    }

    gretl_matrix_free(Svar);
    gretl_matrix_free(vQ);

    return err;
}

static int check_matrix_dims (kalman *K, const gretl_matrix *m, int i)
{
    int r = 0, c = 0, symm = (i == K_Q || i == K_R);
    int err = 0;

    if (i == K_T || i == K_Q || i == K_P) {
        r = c = K->r;
    } else if (i == K_BT) {
        r = K->k;
        c = K->n;
    } else if (i == K_ZT) {
        r = K->r;
        c = K->n;
    } else if (i == K_R)  {
        r = c = K->n;
    } else if (i == K_a || i == K_m) {
        r = K->r;
        c = 1;
    }

    if (m->rows != r || m->cols != c) {
        gretl_errmsg_sprintf("kalman: %s is %d x %d, should be %d x %d\n",
                             kalman_matrix_name(i), m->rows, m->cols, r, c);
        err = E_NONCONF;
    } else if (symm && !gretl_matrix_is_symmetric(m)) {
        gretl_errmsg_sprintf("kalman: %s is not symmetric\n",
                kalman_matrix_name(i));
        err = E_NONCONF;
    }

    return err;
}

enum {
    K_V,
    K_F,
    K_A,
    K_BIG_P,
    K_K,
    K_LL,
};

static int maybe_resize_export_matrix (kalman *K, gretl_matrix *m, int i)
{
    int rows = K->N, cols = 0;
    int err = 0;

    if (i == K_V) {
        cols = K->n;
    } else if (i == K_F) {
        cols = (K->n * K->n + K->n) / 2;
    } else if (i == K_A) {
        cols = K->r;
    } else if (i == K_BIG_P) {
        cols = (K->r * K->r + K->r) / 2;
    } else if (i == K_LL) {
        cols = 1;
    } else if (i == K_K) {
        cols = K->r * K->n;
    } else {
        err = E_DATA;
    }

    if (!err && (m->rows != rows || m->cols != cols)) {
        err = gretl_matrix_realloc(m, rows, cols);
        if (!err) {
            gretl_matrix_zero(m);
        }
    }

    return err;
}

static int missing_matrix_error (const char *name)
{
    if (name == NULL) {
        gretl_errmsg_set(_("kalman: a required matrix is missing"));
    } else {
        gretl_errmsg_sprintf(_("kalman: required matrix %s is missing"),
                             name);
    }
    return E_DATA;
}

static int kalman_check_dimensions (kalman *K)
{
    int err = 0;

    if (K->r < 1 || K->n < 1 || K->N < 2) {
        /* the state and observation vectors must have at least one
           element, and there must be at least two observations
        */
        err = E_DATA;
    }

    /* T is mandatory, should be r x r */
    if (!err) {
        err = check_matrix_dims(K, K->T, K_T);
    }

    /* Z' is mandatory, should be r x n */
    if (!err) {
        err = check_matrix_dims(K, K->ZT, K_ZT);
    }

    if (K->H != NULL) {
        /* cross-correlated disturbances */
        if (K->G == NULL) {
            err = missing_matrix_error("obsymat");
        } else if (K->H->rows != K->r || K->H->cols != K->p ||
                   K->G->rows != K->n || K->G->cols != K->p) {
            err = E_NONCONF;
        }
    } else {
        /* "Q" = HH' is mandatory, should be r x r and symmetric */
        if (!err) {
            err = check_matrix_dims(K, K->HH, K_Q);
        }

        /* "R" = GG' should be n x n and symmetric, if present */
        if (!err && K->GG != NULL) {
            err = check_matrix_dims(K, K->GG, K_R);
        }
    }

    /* initial a should be r x 1, if present */
    if (!err && K->aini != NULL) {
        err = check_matrix_dims(K, K->aini, K_a);
    }

    /* initial P should be r x r, if present */
    if (!err && K->Pini != NULL) {
        err = check_matrix_dims(K, K->Pini, K_P);
    }

    /* B' should be k x n, if present (BT->rows defines k) */
    if (!err) {
        if (K->BT != NULL) {
            err = check_matrix_dims(K, K->BT, K_BT);
        } else if (K->x != NULL) {
            /* B is NULL => can't have a non-NULL x */
            err = E_NONCONF;
        }
    }

    /* mu should be r x 1, if present */
    if (!err && K->mu != NULL) {
        err = check_matrix_dims(K, K->mu, K_m);
    }

    if (err) {
        goto bailout;
    }

    K->ifc = 0;

    /* x should have T rows to match y; and it should have either k or k - 1
       columns (the latter case indicating an implicit const) */
    if (K->x != NULL) {
        if (K->x->rows < K->N) {
            fprintf(stderr, "kalman: %s has %d rows, should have %d\n",
                    kalman_matrix_name(K_x), K->x->rows, K->N);
            return E_NONCONF;
        } else if (K->x->cols != K->k && K->x->cols != K->k - 1) {
            fprintf(stderr, "kalman: %s has %d columns, should have %d or %d\n",
                    kalman_matrix_name(K_x), K->x->cols, K->k, K->k - 1);
            return E_NONCONF;
        } else if (K->x->cols == K->k - 1) {
            /* register the implicit const */
            K->ifc = 1;
        }
    } else if (K->k == 1) {
        /* B has one row but there's no x => implicit const */
        K->ifc = 1;
    } else if (K->k > 1) {
        /* B has more than one row but there's no x => error */
        return missing_matrix_error("obsxmat");
    }

    /* Below we have the optional "export" matrices for shipping out
       results. If these are present but not sized correctly we'll
       try to fix them up -- but note that they are not used in a
       simulation run.
    */

    if (kalman_simulating(K)) {
        return err;
    }

    /* big V should be T x n */
    if (K->V != NULL) {
        err = maybe_resize_export_matrix(K, K->V, K_V);
    }

    /* big F should be T x nn */
    if (!err && K->F != NULL) {
        err = maybe_resize_export_matrix(K, K->F, K_F);
    }

    /* big A should be T x r */
    if (!err && K->A != NULL) {
        err = maybe_resize_export_matrix(K, K->A, K_A);
    }

    /* big P should be T x nr */
    if (!err && K->P != NULL) {
        err = maybe_resize_export_matrix(K, K->P, K_BIG_P);
    }

    /* LL should be T x 1 */
    if (!err && K->LL != NULL) {
        err = maybe_resize_export_matrix(K, K->LL, K_LL);
    }

    /* K (gain) should be T x (r * n) */
    if (!err && K->K != NULL) {
        err = maybe_resize_export_matrix(K, K->K, K_K);
    }

 bailout:

    if (err) {
        fprintf(stderr, "kalman_check_dimensions: err = %d\n", err);
    }

    return err;
}

/* variant of gretl_matrix_copy_values for us when we already
   know that the matrices are non-NULL and conformable
*/

static inline void fast_copy_values (gretl_matrix *B, const gretl_matrix *A)
{
    memcpy(B->val, A->val, B->rows * B->cols * sizeof(double));
}

/* Write the vech of @src into row @t of @targ */

static void load_to_vech (gretl_matrix *targ,
                          const gretl_matrix *src,
                          int n, int t)
{
    int i, j, m = 0;
    double x;

    for (i=0; i<n; i++) {
        for (j=i; j<n; j++) {
            x = gretl_matrix_get(src, i, j);
            gretl_matrix_set(targ, t, m++, x);
        }
    }
}

/* Write the vec of @src into row @t of @targ */

static void load_to_vec (gretl_matrix *targ,
                         const gretl_matrix *src,
                         int t)
{
    int i;

    for (i=0; i<targ->cols; i++) {
        gretl_matrix_set(targ, t, i, src->val[i]);
    }
}

/* Write the square root of diagonal of square matrix @src
   into row @t of @targ, starting at column offset @j
*/

static void load_to_diag (gretl_matrix *targ,
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

/* copy from vector @src into row @t of @targ */

static void load_to_row (gretl_matrix *targ,
                         const gretl_vector *src,
                         int t)
{
    int i;

    for (i=0; i<targ->cols; i++) {
        gretl_matrix_set(targ, t, i, src->val[i]);
    }
}

/* copy from vector @src into row @t of @targ,
   starting at column offset @j in @targ */

static void load_to_row_offset (gretl_matrix *targ,
                                const gretl_vector *src,
                                int t, int j)
{
    int i, n = gretl_vector_get_length(src);

    for (i=0; i<n; i++) {
        gretl_matrix_set(targ, t, i+j, src->val[i]);
    }
}

static void set_row_to_value (gretl_matrix *targ, int t, double x)
{
    int i;

    for (i=0; i<targ->cols; i++) {
        gretl_matrix_set(targ, t, i, x);
    }
}

static int kalman_init (kalman *K)
{
    int err = 0;

    K->SSRw = NADBL;
    K->loglik = NADBL;
    K->s2 = NADBL;

    clear_gretl_matrix_err();

    if (K->aini != NULL) {
        K->a0 = gretl_matrix_copy(K->aini);
        K->a1 = gretl_matrix_copy(K->aini);
    } else {
        K->a0 = gretl_zero_matrix_new(K->r, 1);
        K->a1 = gretl_zero_matrix_new(K->r, 1);
    }

    if (K->Pini != NULL) {
        K->P0 = gretl_matrix_copy(K->Pini);
        K->P1 = gretl_matrix_copy(K->Pini);
    } else {
        K->P0 = gretl_zero_matrix_new(K->r, K->r);
        K->P1 = gretl_zero_matrix_new(K->r, K->r);
    }

    /* forecast error vector, per observation */
    K->v = gretl_matrix_alloc(K->n, 1);

    err = get_gretl_matrix_err();
    if (err) {
        return err;
    }

    K->Blk = gretl_matrix_block_new(&K->PZ,  K->r, K->n, /* P*Z */
                                    &K->Ft,  K->n, K->n, /* (Z'*P*Z + R)^{-1} */
                                    &K->iFt, K->n, K->n, /* Ft-inverse */
                                    &K->Bx,  K->n, 1,    /* B'*x at obs t */
                                    &K->Kt,  K->r, K->n, /* gain at t */
                                    &K->Mt,  K->r, K->n, /* intermediate term */
                                    &K->Ct,  K->r, K->r, /* intermediate term */
                                    NULL);

    if (K->Blk == NULL) {
        err = E_ALLOC;
    }

    if (!err && K->Pini == NULL && !(K->flags & KALMAN_USER)) {
        /* in the "user" case we do this later */
        err = construct_Pini(K);
    }

    return err;
}

static void kalman_set_dimensions (kalman *K)
{
    K->r = gretl_matrix_rows(K->T);  /* T->rows defines r */
    K->k = gretl_matrix_rows(K->BT); /* BT->rows defines k */
    K->n = gretl_matrix_cols(K->y);  /* y->cols defines n */

    if (!kalman_simulating(K)) {
        /* y->rows defines N, except when simulating */
        K->N = gretl_matrix_rows(K->y);
    }

    K->okN = K->N;

    /* K->p is non-zero only under cross-correlation; in that case the
       matrix given as 'Q' in Kalman set-up in fact represents H (as
       in v_t = H \varepsilon_t) and it must be r x p, where p is the
       number of elements in the "combined" disturbance vector
       \varepsilon_t.
    */
    K->p = (K->H != NULL)? gretl_matrix_cols(K->H): 0;
}

/* supports hansl function for creating a named Kalman bundle */

kalman *kalman_new_minimal (gretl_matrix *M[], int copy[],
                            int nmat, int *err)
{
    gretl_matrix **targ[5];
    kalman *K;
    int i;

    *err = 0;

    if (M[0] == NULL || M[1] == NULL || M[2] == NULL || M[3] == NULL) {
        fprintf(stderr, "kalman_new_minimal: nmat=%d, y=%p, Z=%p, T=%p, Q=%p\n",
                nmat, (void *) M[0], (void *) M[1], (void *) M[2], (void *) M[3]);
        *err = missing_matrix_error(NULL);
        return NULL;
    }

    K = kalman_new_empty(KALMAN_USER | KALMAN_BUNDLE);
    if (K == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    targ[0] = &K->y;
    targ[1] = &K->ZT;
    targ[2] = &K->T;

    if (nmat == 5) {
        K->flags |= KALMAN_CROSS;
        targ[3] = &K->H;
        targ[4] = &K->G;
    } else {
        targ[3] = &K->HH; /* "Q" */
    }

    for (i=0; i<nmat; i++) {
        if (copy[i]) {
            *targ[i] = gretl_matrix_copy(M[i]);
        } else {
            *targ[i] = M[i];
        }
    }

    kalman_set_dimensions(K);

    if (K->p > 0) {
        *err = kalman_revise_variance(K);
    }
    if (!*err) {
        *err = kalman_check_dimensions(K);
    }
    if (!*err) {
        *err = kalman_init(K);
    }

    if (*err) {
        kalman_free(K);
        K = NULL;
    } else {
        gretl_matrix_zero(K->v);
    }

    return K;
}

static int matrix_is_varying (kalman *K, int i)
{
    if (K->matcall != NULL) {
        if (K->varying == NULL) {
            check_for_matrix_updates(K, NULL);
        }
        if (K->varying != NULL) {
            return K->varying[i];
        }
    }

    return 0;
}

enum {
    UPDATE_INIT, /* initialization of matrices */
    UPDATE_STEP  /* refreshing matrices per time-step */
};

/* After reading 'Q' = H and 'R' = G from the user, either at
   (re-)initialization or at a given time-step in the case where either
   of these matrices is time-varying, record the user input in K->H
   and K->G and form the 'real' Q and R.  But note that in the
   time-step case it may be that only one of Q, R needs to be treated
   in this way (if only one is time-varying, only one will have
   been redefined via a function call).
*/

static int kalman_update_crossinfo (kalman *K, int mode)
{
    int err = 0;

    /* Note that H and G may be needed as such for simulation */

    if (mode == UPDATE_INIT || matrix_is_varying(K, K_Q)) {
        /* recreate HH' using modified H */
        err = gretl_matrix_multiply_mod(K->H, GRETL_MOD_NONE,
                                        K->H, GRETL_MOD_TRANSPOSE,
                                        K->HH, GRETL_MOD_NONE);
    }

    if (!err && (mode == UPDATE_INIT || matrix_is_varying(K, K_R))) {
        /* recreate GG' using modified G */
        err = gretl_matrix_multiply_mod(K->G, GRETL_MOD_NONE,
                                        K->G, GRETL_MOD_TRANSPOSE,
                                        K->GG, GRETL_MOD_NONE);
    }

    if (!err && (mode == UPDATE_INIT || matrix_is_varying(K, K_Q) ||
                 matrix_is_varying(K, K_R))) {
        /* recreate HG' using modified H and/or G */
        err = gretl_matrix_multiply_mod(K->H, GRETL_MOD_NONE,
                                        K->G, GRETL_MOD_TRANSPOSE,
                                        K->HG, GRETL_MOD_NONE);
    }

    return err;
}

/* kalman_revise_variance: the user has actually given H in place of Q
   and G in place of R; so we have to form Q = HH', R = GG', and HG'
   (cross-correlated disturbances).

   This function is called in the course of initial set-up of a
   filter, and also when the Kalman matrices are being re-checked at
   the start of filtering, smoothing or simulation.
*/

static int kalman_revise_variance (kalman *K)
{
    int err = 0;

    if (K->H == NULL || K->G == NULL) {
        return missing_matrix_error("'statevar' or 'obsvar'");
    }

    err = kalman_update_crossinfo(K, UPDATE_INIT);

    if (err) {
        fprintf(stderr, "kalman_revise_variance: err = %d\n", err);
    }

    return err;
}

#if KDEBUG > 1
static void kalman_print_state (kalman *K)
{
   int j;

    /* if (t > 5) return; */

    fprintf(stderr, "t = %d:\n", K->t);

    for (j=0; j<K->n; j++) {
        fprintf(stderr, "y[%d] = %.10g, err[%d] = %.10g\n", j,
                gretl_matrix_get(K->y, K->t, j),
                j, gretl_vector_get(K->v, j));
    }

    gretl_matrix_print(K->a0, "K->a0");
    gretl_matrix_print(K->P0, "K->P0");
}
#endif

static int kalman_record_state (kalman *K)
{
    int err = 0;

    if (K->A != NULL) {
        load_to_row(K->A, K->a0, K->t);
    }
    if (K->P != NULL) {
	load_to_vech(K->P, K->P0, K->r, K->t);
    }

#if EXACT_SM
    if (K->PK != NULL && K->t < K->r) {
	fprintf(stderr, "HERE record_state for smoother: t=%d, d=%d\n",
		K->t, K->d);
	gretl_matrix_print(K->P0,  "P0");
	gretl_matrix_print(K->Pk0, "Pk0");
	gretl_matrix_print(K->ZT,  "Zt'");
	gretl_matrix_print(K->T,   "Tt");
	gretl_matrix_print(K->Ft,  "Ft");
	gretl_matrix_print(K->Fk,  "Fkt");
	gretl_matrix_print(K->v,   "vt");
	gretl_matrix_print(K->HH,  "HH");
	gretl_matrix_print(K->GG,  "GG");
	//load_to_vech(K->PK, K->Pk0, K->r, K->t);
    }
#endif

    return err;
}

/* Read from the appropriate row of x (N x k) and multiply by B' to
   form B'x_t.  Note: the flag K->ifc is used to indicate that the
   observation equation has an implicit constant, with an entry in
   the B matrix (the first) but no explicit entry in the x matrix.

   The case where x is NULL and B an n-vector (implicit constant)
   is also handled.

   Return 1 if missing values encountered, otherwise 0.
*/

static int kalman_set_Bx (kalman *K)
{
    double xjt, bxi;
    int i, j, missobs = 0;

    for (i=0; i<K->n && !missobs; i++) {
        if (K->x == NULL) {
            /* the implicit constant case */
            bxi = K->BT->val[i];
        } else {
            bxi = 0;
            for (j=0; j<K->k && !missobs; j++) {
                if (K->ifc) {
                    xjt = (j == 0)? 1.0 : gretl_matrix_get(K->x, K->t, j-1);
                } else {
                    xjt = gretl_matrix_get(K->x, K->t, j);
                }
                if (na(xjt)) {
                    missobs = 1;
                } else {
                    bxi += xjt * gretl_matrix_get(K->BT, j, i);
                }
            }
        }
        gretl_vector_set(K->Bx, i, bxi);
    }

    return missobs;
}

/* Compute the one-step ahead forecast error:

   v_t = y_t - Bx - Za

   Return 1 if missing values found, otherwise 0.
*/

static int compute_forecast_error (kalman *K)
{
    int i, missobs = 0;
    double yti;

    /* initialize v_t to y_t */
    for (i=0; i<K->n && !missobs; i++) {
        yti = gretl_matrix_get(K->y, K->t, i);
        K->v->val[i] = yti;
        if (na(yti)) {
            missobs = 1;
        }
    }

    if (K->BT != NULL && !missobs) {
        /* subtract effect of exogenous terms, if any */
        missobs = kalman_set_Bx(K);
        if (!missobs) {
            gretl_matrix_subtract_from(K->v, K->Bx);
        }
    }

    if (!missobs) {
        /* subtract contribution from state */
        gretl_matrix_multiply_mod(K->ZT,  GRETL_MOD_TRANSPOSE,
                                  K->a0, GRETL_MOD_NONE,
                                  K->v, GRETL_MOD_DECREMENT);
    }

    return missobs;
}

/* Given a unified function to update one or more of the
   potentially time-varying matrices, try to figure out
   which matrix or matrices are actually modified by
   this function.
*/

static int check_for_matrix_updates (kalman *K, ufunc *uf)
{
    char **lines;
    int nlines = 0;

    if (K->varying != NULL) {
        free(K->varying);
        K->varying = NULL;
    }

    if (uf == NULL) {
        uf = get_user_function_by_name(K->matcall);
        if (uf == NULL) {
            gretl_errmsg_sprintf("Couldn't find function '%s'", K->matcall);
            return E_DATA;
        }
    }

    K->varying = calloc(K_N_MATCALLS, 1);

    lines = gretl_function_retrieve_code(uf, &nlines);

    if (lines != NULL) {
        const char *bname = fn_param_name(uf, 0);
        char test[VNAMELEN+1];
        const char *s;
        int n = strlen(bname) + 1;
        int i, j;

        sprintf(test, "%s.", bname);
        for (i=0; i<nlines; i++) {
            if (!strncmp(lines[i], test, n)) {
                for (j=K_T; j<=K_m; j++) {
                    s = kalman_matrix_name(j);
                    if (!strncmp(lines[i] + n, s, strlen(s))) {
                        fprintf(stderr, "matrix %s is varying\n", s);
                        K->varying[j] = 1;
                        break;
                    }
                }
            }
        }
        free(lines);
    }

    return 0;
}

/* Function to update any time-varying matrices, for use
   with a kalman bundle. Bypasses the regular "genr" apparatus,
   passing the attached bundle directly to the given user
   function after is has been found by name.
*/

static int kalman_update_matrices (kalman *K, PRN *prn)
{
    ufunc *uf;
    fncall *fc;
    int err = 0;

    uf = get_user_function_by_name(K->matcall);
    if (uf == NULL) {
        gretl_errmsg_sprintf("Couldn't find function '%s'", K->matcall);
        return E_DATA;
    }

    if (K->varying == NULL) {
        check_for_matrix_updates(K, uf);
    }

    fc = fncall_new(uf, 0);
    err = push_anon_function_arg(fc, GRETL_TYPE_BUNDLE_REF, K->b);

    if (!err) {
        err = gretl_function_exec(fc, GRETL_TYPE_NONE, NULL, NULL,
                                  NULL, prn);
    }
    if (err) {
        fprintf(stderr, "kalman_update_matrices: call='%s', err=%d\n",
                K->matcall, err);
    }

    return err;
}

/* If we have any time-varying coefficient matrices, refresh these for
   the current time step. This is called on a forward filtering pass.
*/

static int kalman_refresh_matrices (kalman *K, PRN *prn)
{
    gretl_matrix **mptr[] = {
        &K->T, &K->BT, &K->ZT, &K->HH, &K->GG, &K->mu
    };
    int cross_update = 0;
    int i, err = 0;

    if (kalman_xcorr(K)) {
        mptr[3] = &K->H;
        mptr[4] = &K->G;
    }

    if (K->matcall != NULL) {
        err = kalman_update_matrices(K, prn);
    }

    for (i=0; i<K_N_MATCALLS && !err; i++) {
        if (matrix_is_varying(K, i)) {
            if (kalman_xcorr(K) && (i == K_Q || i == K_R)) {
                /* handle revised H and/or G */
                cross_update = 1;
            } else {
                err = check_matrix_dims(K, *mptr[i], i);
            }
            if (err) {
                fprintf(stderr, "kalman_refresh_matrices: err = %d at t = %d\n",
                        err, K->t);
            }
        }
    }

    if (!err && K->step != NULL) {
        /* keep a record of T and/or Z' at the given time step */
        if (K->step->T != NULL) {
            load_to_vec(K->step->T, K->T, K->t);
        }
        if (K->step->ZT != NULL) {
            load_to_vec(K->step->ZT, K->ZT, K->t);
        }
    }

    if (!err && cross_update) {
        /* cross-correlated case */
        err = kalman_update_crossinfo(K, UPDATE_STEP);
    }

    return err;
}

/* Variant of the above for use when Koopman-smoothing */

static int ksmooth_refresh_matrices (kalman *K, PRN *prn)
{
    gretl_matrix **mptr[] = {
        &K->HH, &K->GG
    };
    int idx[] = {
        K_Q, K_R
    };
    int cross_update = 0;
    int i, ii, err = 0;

    if (kalman_xcorr(K)) {
        mptr[0] = &K->H;
        mptr[1] = &K->G;
    }

    if (K->matcall != NULL) {
        err = kalman_update_matrices(K, prn);
    }

    for (i=0; i<2 && !err; i++) {
        ii = idx[i];
        if (matrix_is_varying(K, ii)) {
            if (kalman_xcorr(K) && (ii == K_Q || ii == K_R)) {
                /* handle revised H and/or G */
                cross_update = 1;
            } else {
                err = check_matrix_dims(K, *mptr[i], ii);
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

/* exact initial iteration for multivariate y_t */

static int koopman_exact_general (kalman *K,
                                  double *ldet,
                                  double *qt)
{
    gretl_matrix_block *B;
    gretl_matrix *Mk;
    gretl_matrix *Fmk;
    gretl_matrix *Fmt;
    gretl_matrix *MFk;
    gretl_matrix *PkZ;
    gretl_matrix *V, *J;
    gretl_matrix *l = NULL;
    double xij, rlj;
    int i, j, ns;
    size_t sz;
    int err = 0;

    B = gretl_matrix_block_new(&Mk,  K->r, K->n,
                               &Fmk, K->n, K->n,
                               &Fmt, K->n, K->n,
                               &MFk, K->r, K->n,
                               &PkZ, K->r, K->n,
                               &V,   K->n, K->n,
                               &J,   K->n, K->n,
                               NULL);

    if (B == NULL) {
        err = E_ALLOC;
    } else {
        /* Mk = T*Pk*Z' */
        gretl_matrix_multiply(K->Pk0, K->ZT, PkZ);
        gretl_matrix_multiply(K->T, PkZ, Mk);
        /* eigenanalysis */
        l = gretl_gensymm_eigenvals(K->Fk, K->Ft, V, &err);
    }

    if (err) {
        goto bailout;
    }

    /* @ns counts the "zero" eigenvalues */
    ns = 0;
    for (i=0; i<K->n; i++) {
        if (l->val[i] < 1.0e-12) {
            ns++;
        } else {
            break;
        }
    }

    /* J = V[,1:ns] */
    gretl_matrix_reuse(J, K->n, ns);
    sz = K->n * ns * sizeof(double);
    memcpy(J->val, V->val, sz);
    /* Fmt = J*J' */
    gretl_matrix_multiply_mod(J, GRETL_MOD_NONE,
                              J, GRETL_MOD_TRANSPOSE,
                              Fmt, GRETL_MOD_NONE);

    /* J = V[,ns+1:] ./ sqrt(l[ns+1:]') */
    gretl_matrix_reuse(J, K->n, K->n-ns);
    sz = K->n * J->cols * sizeof(double);
    memcpy(J->val, V->val + K->n*ns, sz);
    for (j=0; j<J->cols; j++) {
        rlj = sqrt(l->val[ns+j]);
        for (i=0; i<K->n; i++) {
            xij = gretl_matrix_get(J, i, j);
            gretl_matrix_set(J, i, j, xij / rlj);
        }
    }
    /* Fmk = J*J' */
    gretl_matrix_multiply_mod(J, GRETL_MOD_NONE,
                              J, GRETL_MOD_TRANSPOSE,
                              Fmk, GRETL_MOD_NONE);
    /* copy to iFt for recording */
    fast_copy_values(K->iFt, Fmk);

    /* MFk = Mk*Fmk */
    gretl_matrix_multiply(Mk, Fmk, MFk);

    /* Kt = Mt*Fmt + MFk */
    fast_copy_values(K->Kt, MFk);
    gretl_matrix_multiply_mod(K->Mt, GRETL_MOD_NONE,
                              Fmt, GRETL_MOD_NONE,
                              K->Kt, GRETL_MOD_CUMULATE);

    /* Ct = Mt*Kt' */
    gretl_matrix_multiply_mod(K->Mt, GRETL_MOD_NONE,
                              K->Kt, GRETL_MOD_TRANSPOSE,
                              K->Ct, GRETL_MOD_NONE);

    /* Mt <- Mt - MFk*Fmt */
    gretl_matrix_multiply_mod(MFk, GRETL_MOD_NONE,
                              Fmt, GRETL_MOD_NONE,
                              K->Mt, GRETL_MOD_DECREMENT);

    /* Ct += MFk*(Mt - MFk*Fmt)' */
    gretl_matrix_multiply_mod(MFk, GRETL_MOD_NONE,
                              K->Mt, GRETL_MOD_TRANSPOSE,
                              K->Ct, GRETL_MOD_CUMULATE);

    /* Ck = MFk*Mk' */
    gretl_matrix_multiply_mod(MFk, GRETL_MOD_NONE,
                              Mk,  GRETL_MOD_TRANSPOSE,
                              K->Ck, GRETL_MOD_NONE);

    gretl_matrix_add_to(Fmk, Fmt);
    *ldet = gretl_matrix_log_determinant(Fmk, &err);
    *qt = gretl_scalar_qform(K->v, Fmt, &err);

 bailout:

    gretl_matrix_block_destroy(B);

    return err;
}

/* exact initial iteration for univariate y_t,
   when F^\infty (Fk) > 0
*/

static int koopman_exact_nonsingular (kalman *K)
{
    gretl_matrix *Mk;
    gretl_matrix *PkZ;
    int err = 0;

    Mk  = gretl_matrix_alloc(K->r, K->n);
    PkZ = gretl_matrix_alloc(K->r, K->n);

    /* Mk = T*Pk*Z' */
    gretl_matrix_multiply(K->Pk0, K->ZT, PkZ);
    gretl_matrix_multiply(K->T, PkZ, Mk);

    /* Kt = Mk ./ Fk[1] */
    fast_copy_values(K->Kt, Mk);
    gretl_matrix_divide_by_scalar(K->Kt, K->Fk->val[0]);

    /* Ct = Mt*Kt' */
    gretl_matrix_multiply_mod(K->Mt, GRETL_MOD_NONE,
                              K->Kt, GRETL_MOD_TRANSPOSE,
                              K->Ct, GRETL_MOD_NONE);

    /* Mt <- Mt - Kt*Ft */
    gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
                              K->Ft, GRETL_MOD_NONE,
                              K->Mt, GRETL_MOD_DECREMENT);

    /* Ct += Kt*(Mt - Kt*Ft)' */
    gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
                              K->Mt, GRETL_MOD_TRANSPOSE,
                              K->Ct, GRETL_MOD_CUMULATE);

    /* Ck = Kt * Mk' */
    gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
                              Mk, GRETL_MOD_TRANSPOSE,
                              K->Ck, GRETL_MOD_NONE);

    gretl_matrix_free(Mk);
    gretl_matrix_free(PkZ);

    return err;
}

static double max_val (const gretl_matrix *m)
{
    int i, n = m->rows * m->cols;
    double ret = 0;

    for (i=0; i<n; i++) {
        if (m->val[i] > ret) {
            ret = m->val[i];
        }
    }

    return ret;
}

static void P_infty_update (kalman *K, int C_zero)
{
    gretl_matrix_qform(K->T, GRETL_MOD_NONE, K->Pk0,
                       K->Pk1, GRETL_MOD_NONE);
    if (!C_zero) {
	gretl_matrix_subtract_from(K->Pk1, K->Ck);
    }
    fast_copy_values(K->Pk0, K->Pk1);

    if (K->d == 0 && max_val(K->Pk0) < K_TINY) {
	K->d = K->t + 1;
    }
}

/* Handling of missing y_t or x_t, for the case of
   univariate y_t (or all y_t values missing).
*/

static void handle_missing_obs (kalman *K)
{
    /* state update: a1 = T*a0 */
    gretl_matrix_multiply(K->T, K->a0, K->a1);
    /* handle stconst if present */
    if (K->mu != NULL) {
        gretl_matrix_add_to(K->a1, K->mu);
    }
    fast_copy_values(K->a0, K->a1);

    /* var(state) update: P1 = T*P0*T' + HH (C = 0) */
    fast_copy_values(K->P1, K->HH);
    gretl_matrix_qform(K->T, GRETL_MOD_NONE, K->P0,
		       K->P1, GRETL_MOD_CUMULATE);
    fast_copy_values(K->P0, K->P1);

    if (K->exact && max_val(K->Pk0) > K_TINY) {
	/* update P∞ (what about K->d?) */
	gretl_matrix_qform(K->T, GRETL_MOD_NONE, K->Pk0,
			   K->Pk1, GRETL_MOD_NONE);
	fast_copy_values(K->Pk0, K->Pk1);
    }

    /* record stuff if wanted */
    if (K->F != NULL) {
        if (kalman_smoothing(K)) {
            set_row_to_value(K->F, K->t, 0.0);
        } else {
            set_row_to_value(K->F, K->t, NADBL); /* ? */
        }
    }
    if (K->LL != NULL) {
        gretl_vector_set(K->LL, K->t, 0.0); /* ? */
    }
    if (K->K != NULL) {
        set_row_to_value(K->K, K->t, 0.0);
    }
    if (K->V != NULL) {
        if (kalman_smoothing(K)) {
            set_row_to_value(K->V, K->t, 0.0);
        } else {
            set_row_to_value(K->V, K->t, NADBL); /* ? */
        }
    }
}

/**
 * kalman_forecast:
 * @K: pointer to Kalman struct.
 * @prn: printing apparatus (or NULL).
 *
 * Generates a series of one-step ahead forecasts for y, based on
 * information in the kalman struct @K.
 *
 * Returns: 0 on success, non-zero on error.
 */

int kalman_forecast (kalman *K, PRN *prn)
{
    double ll0 = K->n * LN_2_PI;
    int err = 0;

#if KDEBUG
    fprintf(stderr, "kalman_forecast: N = %d\n", K->N);
#endif

    if (kalman_diffuse(K) && !K->exact) {
	K->d = K->r;
    } else {
	K->d = 0;
    }
    K->SSRw = K->loglik = 0.0;
    K->s2 = NADBL;
    K->okN = K->N;
    set_kalman_running(K);

    for (K->t = 0; K->t < K->N && !err; K->t += 1) {
        int Kt_done = 0;
        int missobs = 0;
        double llt = NADBL;
        double ldet = 0, qt = 0;

#if KDEBUG > 1
        kalman_print_state(K);
#endif

        if (K->A != NULL || K->P != NULL) {
            kalman_record_state(K);
        }

        if (filter_is_varying(K)) {
            /* we have time-varying coefficients */
            err = kalman_refresh_matrices(K, prn);
            if (err) {
                K->loglik = NADBL;
                break;
            }
        }

        /* calculate v_t */
        missobs = compute_forecast_error(K);
        if (missobs) {
            K->okN -= 1;
            /* FIXME the case of multivariate y_t with not all
               elements missing is more complicated than this.
            */
            handle_missing_obs(K);
            continue;
        }

        /* calculate F_t = ZPZ' [+ GG'] */
	if (K->GG != NULL) {
	    fast_copy_values(K->Ft, K->GG);
	    gretl_matrix_qform(K->ZT, GRETL_MOD_TRANSPOSE, K->P0,
			       K->Ft, GRETL_MOD_CUMULATE);
	} else {
	    gretl_matrix_qform(K->ZT, GRETL_MOD_TRANSPOSE, K->P0,
			       K->Ft, GRETL_MOD_NONE);
	}

        /* calculate M_t = TPZ' [+ HG'] */
	gretl_matrix_multiply(K->P0, K->ZT, K->PZ);
        gretl_matrix_multiply(K->T, K->PZ, K->Mt);
        if (K->HG != NULL) {
            gretl_matrix_add_to(K->Mt, K->HG);
        }

        if (K->exact && max_val(K->Pk0) > K_TINY) {
            /* handle initial exact iterations */
	    int C_zero = 0;

            gretl_matrix_qform(K->ZT, GRETL_MOD_TRANSPOSE,
                               K->Pk0, K->Fk, GRETL_MOD_NONE);
            if (K->n == 1 && K->Fk->val[0] < K_TINY) {
		/* univariate F∞ = 0 */
                ldet = log(K->Ft->val[0]);
                K->iFt->val[0] = 1.0 / K->Ft->val[0];
                qt = K->v->val[0] * K->v->val[0] * K->iFt->val[0];
		C_zero = 1;
            } else if (K->n == 1) {
		/* univariate F∞ nonsingular */
                ldet = log(K->Fk->val[0]);
		K->iFt->val[0] = 0;
                qt = 0;
                err = koopman_exact_nonsingular(K);
                Kt_done = 1;
            } else {
		/* multivariate case */
                err = koopman_exact_general(K, &ldet, &qt);
                Kt_done = 1;
            }
	    P_infty_update(K, C_zero);
        } else {
            /* standard Kalman procedure */
            fast_copy_values(K->iFt, K->Ft);
            err = gretl_invert_symmetric_matrix2(K->iFt, &ldet);
            if (err) {
                fprintf(stderr, "kalman_forecast: failed to invert F\n");
            } else {
                qt = gretl_scalar_qform(K->v, K->iFt, &err);
            }
        }

        if (K->F != NULL) {
            /* we're recording F_t for all t */
            if (kalman_smoothing(K)) {
                /* record inverse */
                load_to_vech(K->F, K->iFt, K->n, K->t);
            } else {
                /* record F_t itself */
                load_to_vech(K->F, K->Ft, K->n, K->t);
            }
        }

        /* determine and record loglikelihood */
        if (err) {
            K->loglik = NADBL;
            break;
        } else {
            llt = -0.5 * (ll0 + ldet + qt);
            if (na(llt)) {
#if 0
                fprintf(stderr, "kfilter: t=%d, ldet %g, qt %g\n", K->t, ldet, qt);
#endif
                K->loglik = NADBL;
                break;
            }
            K->SSRw += qt;
            K->loglik += llt;
        }
        if (K->LL != NULL) {
            gretl_vector_set(K->LL, K->t, llt);
        }

        if (!Kt_done) {
            /* Calculate gain K_t = M_t F_t^{-1}, and C matrix.
               Note that this will have been done correctly already
               in two of the "exact initial" cases above.
            */
            gretl_matrix_multiply(K->Mt, K->iFt, K->Kt);
            gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
                                      K->Mt, GRETL_MOD_TRANSPOSE,
                                      K->Ct, GRETL_MOD_NONE);
        }

        if (!err && K->K != NULL) {
            /* record the gain */
            load_to_vec(K->K, K->Kt, K->t);
        }

        /* update state: a1 = T a0 + K_t v_t [+ mu] */
        gretl_matrix_multiply(K->T, K->a0, K->a1);
        if (K->mu != NULL) {
            gretl_matrix_add_to(K->a1, K->mu);
        }
        if (!missobs) {
            gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
                                      K->v,  GRETL_MOD_NONE,
                                      K->a1, GRETL_MOD_CUMULATE);
        }
        fast_copy_values(K->a0, K->a1);

        /* update var(state): P1 = TPT' + HH' - C */
        fast_copy_values(K->P1, K->HH);
        gretl_matrix_qform(K->T, GRETL_MOD_NONE,
                           K->P0, K->P1, GRETL_MOD_CUMULATE);
        gretl_matrix_subtract_from(K->P1, K->Ct);
        fast_copy_values(K->P0, K->P1);

        /* record forecast errors if wanted */
        if (!err && K->V != NULL) {
            load_to_row(K->V, K->v, K->t);
        }
    }

    set_kalman_stopped(K);

    if (na(K->loglik)) {
	err = E_NAN;
    } else {
        K->s2 = K->SSRw / (K->n * K->okN - K->d);
    }

#if KDEBUG
    fprintf(stderr, "kalman_forecast: err=%d, ll=%#.12g, d=%d\n",
            err, K->loglik, K->d);
#endif

    return err;
}

struct K_input_mat {
    int sym;
    const char *name;
};

/* mapping to names used in setting elements of kalman bundle */

struct K_input_mat K_input_mats[] = {
    { K_y,  "obsy" },
    { K_ZT, "obsymat" },
    { K_x,  "obsx" },
    { K_BT,  "obsxmat" },
    { K_R,  "obsvar" },
    { K_T,  "statemat" },
    { K_Q,  "statevar" },
    { K_m,  "stconst" },
    { K_a,  "inistate" },
    { K_P,  "inivar" }
};

int extra_mats[] = {
    K_BT,
    K_R,
    K_m,
    K_x,
    K_a,
    K_P,
};

static int n_extra_mats = G_N_ELEMENTS(extra_mats);

static int obsy_check (kalman *K)
{
    if (K->y == NULL) {
        return missing_matrix_error("obsy");
    } else if (K->y->rows != K->N || K->y->cols != K->n) {
        fprintf(stderr, "obsy_check: K->y should be %d x %d, is %d x %d\n",
                K->N, K->n, gretl_matrix_rows(K->y),
                gretl_matrix_cols(K->y));
        return E_NONCONF;
    } else {
        return 0;
    }
}

static int kalman_bundle_recheck_matrices (kalman *K, PRN *prn)
{
    int err = 0;

    K->flags |= KALMAN_CHECK;

    if (filter_is_varying(K)) {
        err = kalman_update_matrices(K, prn);
    }

    K->flags ^= KALMAN_CHECK;

    if (!err && (K->ZT == NULL || K->T == NULL || K->HH == NULL)) {
        fprintf(stderr, "kalman_bundle_kalman_recheck_matrices: Z=%p, T=%p, Q=%p\n",
                K->ZT, K->T, K->HH);
        err = missing_matrix_error(NULL);
    }

    if (err) {
        return err;
    }

    /* redundant? */
    kalman_set_dimensions(K);

    if (gretl_matrix_rows(K->T) != K->r ||
        gretl_matrix_rows(K->BT) != K->k) {
        err = E_NONCONF;
    } else if (!kalman_simulating(K)) {
        err = obsy_check(K);
    }

    if (!err && K->p > 0) {
        err = kalman_revise_variance(K);
    }

    if (!err) {
        err = kalman_check_dimensions(K);
    }

    if (!err) {
        if (K->aini != NULL) {
            gretl_matrix_copy_values(K->a0, K->aini);
        } else {
            gretl_matrix_zero(K->a0);
        }
        if (K->Pini != NULL) {
            gretl_matrix_copy_values(K->P0, K->Pini);
        } else {
            err = construct_Pini(K);
        }
    }

    return err;
}

static const char *kalman_matrix_name (int sym)
{
    int i;

    for (i=0; i<K_MMAX; i++) {
        if (K_input_mats[i].sym == sym) {
            return K_input_mats[i].name;
        }
    }

    /* failed */
    return "matrix";
}

static int kalman_ensure_output_matrices (kalman *K)
{
    int err = 0;

    if (K->V == NULL) {
        K->V = gretl_null_matrix_new();
    }
    if (K->F == NULL) {
        K->F = gretl_null_matrix_new();
    }
    if (K->A == NULL) {
        K->A = gretl_null_matrix_new();
    }
    if (K->P == NULL) {
        K->P = gretl_null_matrix_new();
    }
    if (K->K == NULL) {
        K->K = gretl_null_matrix_new();
    }

    if (K->V == NULL || K->F == NULL || K->A == NULL ||
        K->P == NULL || K->K == NULL) {
        err = E_ALLOC;
    }

    return err;
}

/* Implements the user-space kfilter() function */

int kalman_bundle_run (gretl_bundle *b, PRN *prn, int *errp)
{
    kalman *K = gretl_bundle_get_private_data(b);
    int err;

    K->b = b; /* attach bundle pointer */
    err = kalman_ensure_output_matrices(K);

    if (!err) {
        gretl_matrix_zero(K->v);
        err = kalman_bundle_recheck_matrices(K, prn);
    }

    if (!err && K->LL == NULL) {
        K->LL = gretl_matrix_alloc(K->N, 1);
        if (K->LL == NULL) {
            err = E_ALLOC;
        }
    }

    if (!err) {
        err = kalman_forecast(K, prn);
    }

    if (err != E_NAN) {
        *errp = err;
    } else {
        /* we'll flag E_NAN with a return value of 1 but
           won't count it as a 'true' error */
        *errp = 0;
    }

    return err;
}

/* Copy row @t from @src into @targ; or add row @t of @src to
   @targ; or subtract row @t of @src from @targ.  We allow the
   possibility that the length of vector @targ is less than
   the number of columns in @src, but not the converse.
*/

static int load_from_row (gretl_vector *targ,
                          const gretl_matrix *src,
                          int t, GretlMatrixMod mod)
{
    int i, n = gretl_vector_get_length(targ);
    double x;

    if (n > src->cols) {
        fprintf(stderr, "load_from_row: targ length = %d, but src "
                "has %d columns\n", n, src->cols);
        return E_NONCONF;
    }

    for (i=0; i<n; i++) {
        x = gretl_matrix_get(src, t, i);
        if (mod == GRETL_MOD_CUMULATE) {
            targ->val[i] += x;
        } else if (mod == GRETL_MOD_DECREMENT) {
            targ->val[i] -= x;
        } else {
            targ->val[i] = x;
        }
    }

    return 0;
}

/* As load_from_row(), except that a column offset, @j, is
   supported for the reading of a row from @src, and we
   don't support the @mod option.
*/

static int load_from_row_offset (gretl_vector *targ,
                                 const gretl_matrix *src,
                                 int t, int j)
{
    int i, n = gretl_vector_get_length(targ);
    double x;

    if (n > src->cols - j) {
        return E_NONCONF;
    }

    for (i=0; i<n; i++) {
        x = gretl_matrix_get(src, t, i + j);
        targ->val[i] = x;
    }

    return 0;
}

/* Row @t of @src represents the vech of an n x n matrix: extract the
   row and apply the inverse operation of vech to reconstitute the
   matrix in @targ -- or subtract the newly reconstituted matrix
   from @targ.
*/

static void load_from_vech (gretl_matrix *targ, const gretl_matrix *src,
                            int n, int t, int mod)
{
    int i, j, m = 0;
    double x;

    for (i=0; i<n; i++) {
        for (j=i; j<n; j++) {
            x = gretl_matrix_get(src, t, m++);
            if (mod == GRETL_MOD_DECREMENT) {
                x = gretl_matrix_get(targ, i, j) - x;
            }
            gretl_matrix_set(targ, i, j, x);
            if (i != j) {
                gretl_matrix_set(targ, j, i, x);
            }
        }
    }
}

/* Row @t of @src represents the vec of a certain matrix: extract the
   row and reconstitute the matrix in @targ.
*/

static int load_from_vec (gretl_matrix *targ,
                          const gretl_matrix *src,
                          int t)
{
    int i, k = targ->rows * targ->cols;

    for (i=0; i<k; i++) {
        targ->val[i] = gretl_matrix_get(src, t, i);
    }

    return 0;
}

/* For disturbance smoothing: ensure we have on hand
   matrices that are correctly sized to hold estimates
   of the variance of the disturbance(s) in the state
   and (if applicable) observation equations.
*/

static int maybe_resize_dist_mse (kalman *K,
                                  gretl_matrix **vvt,
                                  gretl_matrix **vwt)
{
    int n = K->GG == NULL ? 0 : K->n;
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
        *vvt = gretl_matrix_alloc(K->r, K->r);
        if (*vvt == NULL) {
            err = E_ALLOC;
        }
    }

    if (!err && n > 0) {
        /* step-t square obs matrix */
        *vwt = gretl_matrix_alloc(K->n, K->n);
        if (*vwt == NULL) {
            err = E_ALLOC;
        }
    }

    return err;
}

static int retrieve_Tt (kalman *K, int t)
{
    if (K->step == NULL || K->step->T == NULL) {
        return E_DATA;
    } else {
        return load_from_vec((gretl_matrix *) K->T, K->step->T, t);
    }
}

static int retrieve_Zt (kalman *K, int t)
{
    if (K->step == NULL || K->step->ZT == NULL) {
        return E_DATA;
    } else {
        return load_from_vec((gretl_matrix *) K->ZT, K->step->ZT, t);
    }
}

/* Calculate the variance of the smoothed disturbances
   for the cross-correlated case. See Koopman, Shephard
   and Doornik (1998), page 19, var(\varepsilon_t|Y_n).
*/

static int combined_dist_variance (kalman *K,
                                   gretl_matrix *D,
                                   gretl_matrix *N,
                                   gretl_matrix *vv,
                                   gretl_matrix *vw,
                                   gretl_matrix_block *BX,
                                   int dkstyle)
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
                              N, GRETL_MOD_NONE,
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
    gretl_matrix_multiply(N, K->H, NH);
    gretl_matrix_multiply_mod(N, GRETL_MOD_TRANSPOSE,
                              K->Kt, GRETL_MOD_NONE,
                              NK, GRETL_MOD_NONE);
    gretl_matrix_multiply_mod(NK, GRETL_MOD_NONE,
                              K->G, GRETL_MOD_NONE,
                              NH, GRETL_MOD_DECREMENT);
    gretl_matrix_multiply_mod(K->H, GRETL_MOD_TRANSPOSE,
                              NH, GRETL_MOD_NONE,
                              Veps, GRETL_MOD_CUMULATE);

    if (dkstyle) {
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
       disturbance variance matrices, @vv and @vw, for
       this time-step.
    */
    gretl_matrix_qform(K->H, GRETL_MOD_NONE, Veps,
                       vv, GRETL_MOD_NONE);
    gretl_matrix_qform(K->G, GRETL_MOD_NONE, Veps,
                       vw, GRETL_MOD_NONE);

    return 0;
}

/* See Koopman, Shephard and Doornik, section 4.4 */

static int koopman_smooth (kalman *K, int dkstyle)
{
    gretl_matrix_block *B, *BX = NULL;
    gretl_matrix *u, *D, *L, *R;
    gretl_matrix *r0, *r1, *r2, *N1, *N2, *n1, *tr;
    gretl_matrix *Vvt = NULL;
    gretl_matrix *Vwt = NULL;
    gretl_matrix *DG = NULL;
    gretl_matrix *KN = NULL;
    gretl_matrix *RZS = NULL;
    gretl_matrix *NH = NULL;
    gretl_matrix *Ut = NULL;
    double x;
    int i, t, err = 0;

    B = gretl_matrix_block_new(&u,  K->n, 1,
                               &D,  K->n, K->n,
                               &L,  K->r, K->r,
                               &R,  K->N, K->r,
                               &r0, K->r, 1,
                               &r1, K->r, 1,
                               &r2, K->r, 1,
                               &N1, K->r, K->r,
                               &N2, K->r, K->r,
                               &n1, K->n, 1,
                               &tr, K->r, 1,
                               NULL);

    if (B == NULL) {
        return E_ALLOC;
    }

    if (K->b != NULL) {
        /* for variance of smoothed disturbances */
        err = maybe_resize_dist_mse(K, &Vvt, &Vwt);
    }

    if (K->b != NULL && K->p > 0) {
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
        gretl_matrix_free(Vvt);
        gretl_matrix_free(Vwt);
        return err;
    }

    gretl_matrix_zero(r1);
    gretl_matrix_zero(N1);

    /* The backward recursion */

    for (t=K->N-1; t>=0 && !err; t--) {
        /* load vt, Ft and Kt for time t */
        load_from_row(K->v, K->V, t, GRETL_MOD_NONE);
        load_from_vech(K->Ft, K->F, K->n, t, GRETL_MOD_NONE);
        load_from_vec(K->Kt, K->K, t);

        /* get T_t and/or Z_t if need be */
        if (matrix_is_varying(K, K_T)) {
            err = retrieve_Tt(K, t);
        }
        if (!err && matrix_is_varying(K, K_ZT)) {
            err = retrieve_Zt(K, t);
        }
        if (err) {
            break;
        }

        if (filter_is_varying(K)) {
            /* K->HH and/or K->GG may be time-varying */
            K->t = t;
            ksmooth_refresh_matrices(K, NULL);
        }

        /* u_t = F_t^{-1} v_t - K_t' r_t */
        gretl_matrix_multiply(K->Ft, K->v, u);
        if (t < K->N - 1) {
            gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_TRANSPOSE,
                                      r1, GRETL_MOD_NONE,
                                      u, GRETL_MOD_DECREMENT);
        }

        /* save u_t values in V */
        load_to_row(K->V, u, t);

        if (K->Vsd != NULL && K->p == 0) {
            /* variance of state disturbance */
            if (dkstyle) {
                /* Q_t - Q_t N_t Q_t */
                fast_copy_values(Vvt, K->HH);
                gretl_matrix_qform(K->HH, GRETL_MOD_TRANSPOSE,
                                   N1, Vvt, GRETL_MOD_DECREMENT);

            } else {
                /* Q_t N_t Q_t */
                gretl_matrix_qform(K->HH, GRETL_MOD_TRANSPOSE,
                                   N1, Vvt, GRETL_MOD_NONE);
            }
            load_to_diag(K->Vsd, Vvt, t, 0);
        }

        /* D_t = F_t^{-1} + K_t' N_t K_t */
        fast_copy_values(D, K->Ft);
        if (t < K->N - 1) {
            gretl_matrix_qform(K->Kt, GRETL_MOD_TRANSPOSE,
                               N1, D, GRETL_MOD_CUMULATE);
        }

        if (K->GG != NULL && K->Vsd != NULL && K->p == 0) {
            /* variance of obs disturbance */
            if (dkstyle) {
                /* R_t - R_t D_t R_t */
                fast_copy_values(Vwt, K->GG);
                gretl_matrix_qform(K->GG, GRETL_MOD_TRANSPOSE,
                                   D, Vwt, GRETL_MOD_DECREMENT);

            } else {
                /* R_t D_t R_t */
                gretl_matrix_qform(K->GG, GRETL_MOD_TRANSPOSE,
                                   D, Vwt, GRETL_MOD_NONE);
            }
            load_to_diag(K->Vsd, Vwt, t, K->r);
        }

        if (K->Vsd != NULL && K->p > 0) {
            /* variance of combined disturbance */
            err = combined_dist_variance(K, D, N1, Vvt, Vwt, BX,
                                         dkstyle);
            if (err) {
                break;
            } else {
                load_to_diag(K->Vsd, Vvt, t, 0);
                load_to_diag(K->Vsd, Vwt, t, K->r);
            }
        }

        /* L_t = T - KZ' */
        fast_copy_values(L, K->T);
        gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
                                  K->ZT, GRETL_MOD_TRANSPOSE,
                                  L, GRETL_MOD_DECREMENT);

        /* save r_t values in R */
        load_to_row(R, r1, t);

        /* r_{t-1} = Z F_t^{-1}v_t + L_t' r_t */
        gretl_matrix_multiply(K->ZT, K->Ft, tr);
        gretl_matrix_multiply(tr, K->v, r2);
        if (t < K->N - 1) {
            gretl_matrix_multiply_mod(L, GRETL_MOD_TRANSPOSE,
                                      r1, GRETL_MOD_NONE,
                                      r2, GRETL_MOD_CUMULATE);
        }
        /* transcribe for next step */
        fast_copy_values(r1, r2);

        /* preserve r_0 for smoothing of state */
        if (t == 0) {
            fast_copy_values(r0, r2);
        }

        /* N_{t-1} = Z F_t^{-1}Z' + L' N L */
        gretl_matrix_qform(K->ZT, GRETL_MOD_NONE,
                           K->Ft, N2, GRETL_MOD_NONE);
        if (t < K->N - 1) {
            gretl_matrix_qform(L, GRETL_MOD_TRANSPOSE,
                               N1, N2, GRETL_MOD_CUMULATE);
        }
        /* transcribe for next step */
        fast_copy_values(N1, N2);
    }

#if 0
    gretl_matrix_print(K->V, "u_t, all t");
    gretl_matrix_print(R, "r_t, all t");
#endif

    /* Smoothed disturbances, all time steps */

    if (K->p > 0) {
        /* eps = B' r_t + G' e_t */
        for (t=0; t<K->N; t++) {
            if (filter_is_varying(K)) {
                K->t = t;
                ksmooth_refresh_matrices(K, NULL);
            }
            load_from_row(r1, R, t, GRETL_MOD_NONE);
            gretl_matrix_multiply_mod(K->H, GRETL_MOD_TRANSPOSE,
                                      r1, GRETL_MOD_NONE,
                                      Ut, GRETL_MOD_NONE);
            load_from_row(K->v, K->V, t, GRETL_MOD_NONE);
            gretl_matrix_multiply_mod(K->G, GRETL_MOD_TRANSPOSE,
                                      K->v, GRETL_MOD_NONE,
                                      Ut, GRETL_MOD_CUMULATE);
            gretl_matrix_multiply(K->H, Ut, r2);
            load_to_row(R, r2, t);
            gretl_matrix_multiply(K->G, Ut, n1);
            for (i=0; i<K->r; i++) {
                x = gretl_vector_get(r2, i);
                gretl_matrix_set(K->U, t, i, x);
            }
            for (i=0; i<K->n; i++) {
                x = gretl_vector_get(n1, i);
                gretl_matrix_set(K->U, t, K->r + i, x);
            }
        }
    } else {
        /* the independent case */
        for (t=0; t<K->N; t++) {
            if (filter_is_varying(K)) {
                K->t = t;
                ksmooth_refresh_matrices(K, NULL);
            }
            load_from_row(r1, R, t, GRETL_MOD_NONE);
            gretl_matrix_multiply(K->HH, r1, r2);
            load_to_row(R, r2, t);
            if (K->U != NULL) {
                for (i=0; i<K->r; i++) {
                    x = gretl_vector_get(r2, i);
                    gretl_matrix_set(K->U, t, i, x);
                }
            }
            if (n1 != NULL) {
                load_from_row(K->v, K->V, t, GRETL_MOD_NONE);
                gretl_matrix_multiply(K->GG, K->v, n1);
                for (i=0; i<K->n; i++) {
                    x = gretl_vector_get(n1, i);
                    gretl_matrix_set(K->U, t, K->r + i, x);
                }
            }
        }
    }

    /* Write initial smoothed state into first row of S */

    if (K->Pini != NULL) {
        gretl_matrix_multiply(K->Pini, r0, K->a0);
    } else {
        construct_Pini(K);
        gretl_matrix_multiply(K->P0, r0, K->a0);
    }
    if (K->aini != NULL) {
        gretl_matrix_add_to(K->a0, K->aini);
    }
    load_to_row(K->A, K->a0, 0);

    /* Smoothed state, remaining time steps */

    for (t=1; t<K->N; t++) {
        /* a_{t+1} = Ta_t + v_t (or + H*eps_t) */
        load_from_row(K->a0, K->A, t-1, GRETL_MOD_NONE);
        gretl_matrix_multiply(K->T, K->a0, K->a1);
        load_from_row(K->a1, R, t-1, GRETL_MOD_CUMULATE);
        load_to_row(K->A, K->a1, t);
    }

    gretl_matrix_block_destroy(B);
    gretl_matrix_block_destroy(BX);
    gretl_matrix_free(Vvt);
    gretl_matrix_free(Vwt);

    return err;
}

/* Anderson-Moore Kalman smoothing: see Iskander Karibzhanov's
   exposition at http://karibzhanov.com/help/kalcvs.htm
   This is much the clearest account I have seen (AC 2009-04-14,
   URL updated 2016-03-24).

   This method uses a_{t|t-1} and P_{t|t-1} for all t, but we can
   overwrite these with the smoothed values as we go. We also need
   stored values for the prediction error, its MSE, and the gain at
   each time step.  Note that r_t and N_t are set to zero for
   t = N - 1.
*/

static int anderson_moore_smooth (kalman *K)
{
    gretl_matrix_block *B;
    gretl_matrix *r0, *r1, *N0, *N1, *iFv, *L;
    gretl_matrix *atT, *PtT;
    int t, err = 0;

    B = gretl_matrix_block_new(&atT, K->r, 1,
                               &PtT, K->r, K->r,
                               &r0,  K->r, 1,
                               &r1,  K->r, 1,
                               &N0,  K->r, K->r,
                               &N1,  K->r, K->r,
                               &iFv, K->n, 1,
                               &L,   K->r, K->r,
                               NULL);
    if (B == NULL) {
        return E_ALLOC;
    }

    gretl_matrix_zero(r0);
    gretl_matrix_zero(N0);

    for (t=K->N-1; t>=0 && !err; t--) {
        /* get T_t and/or Z_t if need be */
        if (matrix_is_varying(K, K_T)) {
            err = retrieve_Tt(K, t);
        }
        if (!err && matrix_is_varying(K, K_ZT)) {
            err = retrieve_Zt(K, t);
        }
        if (err) {
            break;
        }

#if EXACT_SM
	if (t <= K->d + 1) {
	    fprintf(stderr, "HERE smoothing, t=%d\n", t);
	    gretl_matrix_print(r0, "r0");
	    gretl_matrix_print(N0, "N0");
	}
#endif

        /* L_t = T_t - K_t Z_t */
        fast_copy_values(L, K->T);
        load_from_vec(K->Kt, K->K, t);
        gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
                                  K->ZT, GRETL_MOD_TRANSPOSE,
                                  L, GRETL_MOD_DECREMENT);

        /* r_{t-1} = Z_t' F^{-1}_t v_t + L_t' r_t */
        load_from_vech(K->iFt, K->F, K->n, t, GRETL_MOD_NONE);
        load_from_row(K->v, K->V, t, GRETL_MOD_NONE);
        gretl_matrix_multiply(K->iFt, K->v, iFv);
        gretl_matrix_multiply(K->ZT, iFv, r1);
        if (t == K->N - 1) {
            gretl_matrix_multiply(K->ZT, iFv, r0);
        } else {
            gretl_matrix_multiply_mod(L, GRETL_MOD_TRANSPOSE,
                                      r0, GRETL_MOD_NONE,
                                      r1, GRETL_MOD_CUMULATE);
            fast_copy_values(r0, r1);
        }

        /* N_{t-1} = Z_t' F^{-1}_t Z_t + L_t' N_t L_t */
        if (t == K->N - 1) {
            gretl_matrix_qform(K->ZT, GRETL_MOD_NONE,
                               K->iFt, N0, GRETL_MOD_NONE);
        } else {
            gretl_matrix_qform(K->ZT, GRETL_MOD_NONE,
                               K->iFt, N1, GRETL_MOD_NONE);
            gretl_matrix_qform(L, GRETL_MOD_TRANSPOSE,
                               N0, N1, GRETL_MOD_CUMULATE);
            fast_copy_values(N0, N1);
        }

	/* FIXME: the following calculations have to be modified
	   for the exact diffuse case and t <= d. See Durbin and
	   Koopman (2012), chapter 5, section 3.
	*/

        /* a_{t|T} = a_{t|t-1} + P_{t|t-1} r_{t-1} */
        load_from_row(atT, K->A, t, GRETL_MOD_NONE);
	load_from_vech(K->P0, K->P, K->r, t, GRETL_MOD_NONE);
        gretl_matrix_multiply_mod(K->P0, GRETL_MOD_NONE,
                                  r0, GRETL_MOD_NONE,
                                  atT, GRETL_MOD_CUMULATE);
#if 0 // EXACT_SM
	if (t <= K->d + 1) {
	    fprintf(stderr, "HERE smoothing, t=%d\n", t);
	    gretl_matrix_print(K->P0, "K->P0");
	    gretl_matrix_print(r0, "r0");
	    gretl_matrix_print(atT, "atT");
	}
#endif
        load_to_row(K->A, atT, t);

        /* P_{t|T} = P_{t|t-1} - P_{t|t-1} N_{t-1} P_{t|t-1} */
        fast_copy_values(PtT, K->P0);
        gretl_matrix_qform(K->P0, GRETL_MOD_NONE,
                           N0, PtT, GRETL_MOD_DECREMENT);
        load_to_vech(K->P, PtT, K->r, t);
    }

    gretl_matrix_block_destroy(B);

    return err;
}

/* If we're doing smoothing for a system that has time-varying
   coefficients in K->T or K->ZT we'll record the vec of the
   coefficient matrices for each time-step on the forward pass.
   Here we allocate the required storage.
*/

static int kalman_add_stepinfo (kalman *K)
{
    int err = 0;

    K->step = malloc(sizeof *K->step);

    if (K->step == NULL) {
        return E_ALLOC;
    }

    K->step->T = K->step->ZT = NULL;

    if (matrix_is_varying(K, K_T)) {
        K->step->T = gretl_matrix_alloc(K->N, K->r * K->r);
        if (K->step->T == NULL) {
            err = E_ALLOC;
        }
    }

    if (!err && matrix_is_varying(K, K_ZT)) {
        K->step->ZT = gretl_matrix_alloc(K->N, K->r * K->n);
        if (K->step->ZT == NULL) {
            err = E_ALLOC;
        }
    }

    if (err) {
        free_stepinfo(K);
    }

    return err;
}

/* optional matrix to hold smoothed disturbances a la Koopman */

static int ensure_U_matrix (kalman *K)
{
    int Ucols = K->r;
    int Urows = K->N;
    int err = 0;

    if (K->GG != NULL) {
        Ucols += K->n;
    }

    if (K->U == NULL) {
        K->U = gretl_matrix_alloc(Urows, Ucols);
        if (K->U == NULL) {
            err = E_ALLOC;
        }
    } else if (K->U->rows != Urows || K->U->cols != Ucols) {
        err = gretl_matrix_realloc(K->U, Urows, Ucols);
    }

    return err;
}

/**
 * kalman_smooth:
 * @K: pointer to kalman struct.
 * @pP: pointer to matrix in which to retrieve the MSE of the
 * smoothed state (or NULL if this is not required).
 * @pU: pointer to matrix in which to retrieve the smoothed
 * disturbances (or NULL if this is not required).
 * @err: location to receive error code.
 *
 * Runs a filtering pass followed by a backward, smoothing pass.
 * At present the @pU argument is experimental and a bodge: it will
 * not actually do anything unless @pP is left NULL.
 *
 * Returns: matrix containing the smoothed estimate of the
 * state, or NULL on error.
 */

gretl_matrix *kalman_smooth (kalman *K,
                             gretl_matrix **pP,
                             gretl_matrix **pU,
                             int *err)
{
    gretl_matrix *V, *S;
    gretl_matrix *G, *F;
    gretl_matrix *P = NULL;
    int nr, nn;

    if (pP == NULL && pU != NULL) {
        /* optional accessor for smoothed disturbances a la Koopman:
           experimental, and available only if @pP is not given
        */
        *err = ensure_U_matrix(K);
        if (*err) {
            return NULL;
        }
    }

    nr = (K->r * K->r + K->r) / 2;
    nn = (K->n * K->n + K->n) / 2;

    /* Set up the matrices we need to store computed results from all
       time steps on the forward pass: prediction error, (inverse)
       error variance, gain, state a_{t|t-1} and MSE of state,
       P_{t|t-1}.
    */
    V = gretl_matrix_alloc(K->N, K->n);
    F = gretl_matrix_alloc(K->N, nn);
    G = gretl_matrix_alloc(K->N, K->r * K->n);
    S = gretl_matrix_alloc(K->N, K->r);
    P = gretl_matrix_alloc(K->N, nr);

    if (V == NULL || F == NULL || G == NULL ||
        S == NULL || P == NULL) {
        *err = E_ALLOC;
        goto bailout;
    }

    if (matrix_is_varying(K, K_T) || matrix_is_varying(K, K_ZT)) {
        /* add recorder for T_t and/or Z_t */
        *err = kalman_add_stepinfo(K);
        if (*err) {
            goto bailout;
        }
    }

#if EXACT_SM
    if (K->exact) {
	fprintf(stderr, "smoothing, K->exact, add PK\n");
	K->PK = gretl_zero_matrix_new(K->r, nr);
    }
#endif

    /* attach all export matrices to Kalman */
    K->V = V;
    K->F = F;
    K->K = G;
    K->A = S;
    K->P = P;

#if 0
    /* and recheck dimensions */
    *err = user_kalman_recheck_matrices(K, prn);
#endif

    if (!*err) {
        /* forward pass */
        K->flags |= KALMAN_SMOOTH;
        *err = kalman_forecast(K, NULL);
        K->flags &= ~KALMAN_SMOOTH;
    }

    K->t = 0;

    if (!*err) {
        /* bodge */
        if (K->U != NULL) {
            *err = koopman_smooth(K, 0);
        } else {
            *err = anderson_moore_smooth(K);
        }
    }

    /* detach matrices */
    K->V = NULL;
    K->F = NULL;
    K->K = NULL;
    K->A = NULL;
    K->P = NULL;

    /* and trash the "stepinfo" storage */
    free_stepinfo(K);

 bailout:

    gretl_matrix_free(V);
    gretl_matrix_free(F);
    gretl_matrix_free(G);

    if (K->PK != NULL) {
	gretl_matrix_free(K->PK);
	K->PK = NULL;
    }

    if (!*err && pP != NULL) {
        *pP = P;
    } else {
        gretl_matrix_free(P);
    }

    if (!*err && pU != NULL) {
        *pU = K->U;
    } else {
        gretl_matrix_free(K->U);
    }

    K->U = NULL;

    if (*err) {
        gretl_matrix_free(S);
        S = NULL;
    }

    return S;
}

int kalman_bundle_smooth (gretl_bundle *b, int dist, PRN *prn)
{
    kalman *K = gretl_bundle_get_private_data(b);
    int err;

    if (K == NULL) {
        fprintf(stderr, "kalman_bundle_smooth: K is NULL\n");
        return E_DATA;
    }

    K->b = b; /* attach bundle pointer */

    err = kalman_ensure_output_matrices(K);

    if (!err && dist) {
        err = ensure_U_matrix(K);
    }

    if (err) {
        return err;
    }

    if (matrix_is_varying(K, K_T) || matrix_is_varying(K, K_ZT)) {
        /* add recorder for T_t and/or Z_t */
        err = kalman_add_stepinfo(K);
        if (err) {
            goto bailout;
        }
    }

#if EXACT_SM
    if (K->exact) {
	int rr = (K->r * K->r + K->r) / 2;

	K->PK = gretl_zero_matrix_new(K->r, rr);
    }
#endif

    if (!err) {
        err = kalman_bundle_recheck_matrices(K, prn);
    }

    if (!err) {
        /* forward pass */
        K->flags |= KALMAN_SMOOTH;
        err = kalman_forecast(K, NULL);
        K->flags &= ~KALMAN_SMOOTH;
    }

    K->t = 0;

    if (!err) {
        if (dist > 1) {
            err = koopman_smooth(K, 1);
        } else if (dist == 1) {
            err = koopman_smooth(K, 0);
        } else {
            err = anderson_moore_smooth(K);
        }
    }

#if EXACT_SM
    if (K->PK != NULL) {
	gretl_matrix_free(K->PK);
	K->PK = NULL;
    }
#endif

 bailout:

    /* trash the "stepinfo" storage */
    free_stepinfo(K);

    return err;
}

static gretl_matrix *extract_Q (kalman *K,
                                const gretl_matrix *Sim0)
{
    gretl_matrix *Q;
    double x;
    int i, j;

    Q = gretl_matrix_alloc(K->r, K->r);

    if (Q != NULL) {
        for (i=0; i<K->r; i++) {
            for (j=0; j<K->r; j++) {
                x = gretl_matrix_get(Sim0, i, j);
                gretl_matrix_set(Q, i, j, x);
            }
        }
    }

    return Q;
}

/* See the account in Koopman, Shephard and Doornik, Econometrics
   Journal, 1999 (volume 2, pp. 113-166), section 4.2, regarding
   the initialization of the state under simulation.
*/

static int sim_state_0 (kalman *K, const gretl_matrix *U,
                        const gretl_matrix *Sim0)
{
    gretl_matrix *Q, *v0 = NULL, *bv = NULL;
    int getroot = 1;
    int i, err = 0;

    if (!kalman_ssfsim(K)) {
        if (Sim0 != NULL) {
            /* Sim0 contains the state for t = 1 */
            err = gretl_matrix_copy_values(K->a0, Sim0);
        }
        /* error or not, we're done */
        return err;
    }

    /* now we're in the "ssfsim" case, emulating ssfpack */

    if (Sim0 != NULL) {
        /* Sim0 contains state variance factor
           plus the state for t = 0
        */
        Q = extract_Q(K, Sim0);
        getroot = 0;
    } else {
        Q = gretl_matrix_copy(K->P0);
    }

    if (Q == NULL) {
        err = E_ALLOC;
    } else if (getroot) {
        err = gretl_matrix_psd_root(Q, 0);
    }

    if (!err) {
        int vlen = K->p > 0 ? K->p : K->r;

        v0 = gretl_matrix_alloc(vlen, 1);
        if (v0 == NULL) {
            err = E_ALLOC;
        }
    }

    if (K->p > 0) {
        bv = gretl_matrix_alloc(K->r, 1);
        if (bv == NULL) {
            err = E_ALLOC;
        }
    }

    if (!err && Sim0 != NULL) {
        /* set a0 from last row of Sim0 */
        for (i=0; i<K->r; i++) {
            K->a0->val[i] = gretl_matrix_get(Sim0, K->r, i);
        }
    }

    if (!err) {
        /* handle the t = 0 disturbance */
        load_from_row(v0, U, 0, GRETL_MOD_NONE);
        if (K->p > 0) {
            /* cross-correlated */
            gretl_matrix_multiply(K->H, v0, bv);
            gretl_matrix_multiply_mod(Q, GRETL_MOD_NONE,
                                      bv, GRETL_MOD_NONE,
                                      K->a0, GRETL_MOD_CUMULATE);
        } else {
            gretl_matrix_multiply_mod(Q, GRETL_MOD_NONE,
                                      v0, GRETL_MOD_NONE,
                                      K->a0, GRETL_MOD_CUMULATE);
        }
    }

    gretl_matrix_free(Q);
    gretl_matrix_free(v0);
    gretl_matrix_free(bv);

    return err;
}

/* note: it's OK for @S to be NULL (if the simulated
   state is not wanted), so watch out for that!
*/

static int kalman_simulate (kalman *K,
                            const gretl_matrix *U,
                            const gretl_matrix *Sim0,
                            gretl_matrix *Y,
                            gretl_matrix *S,
                            PRN *prn)
{
    gretl_matrix *yt, *et = NULL;
    int obs_offset = 0;
    int obsdist = 0;
    int tmin = 0;
    int err = 0;

    yt = gretl_zero_matrix_new(K->n, 1);
    if (yt == NULL) {
        return E_ALLOC;
    }

    if (K->p > 0) {
        et = gretl_matrix_alloc(K->p, 1);
        if (et == NULL) {
            gretl_matrix_free(yt);
            return E_ALLOC;
        }
    }

    if (Y->cols == K->r + K->n) {
        /* combined (state, obs) in @Y */
        S = Y;
        obs_offset = K->r;
    }

    err = sim_state_0(K, U, Sim0);

    if (!err && kalman_ssfsim(K)) {
        if (S != NULL) {
            load_to_row_offset(S, K->a0, 0, 0);
        }
        load_to_row_offset(Y, yt, 0, obs_offset);
        /* the first row of output is handled */
        tmin = 1;
    }

    if (K->p == 0 && K->GG != NULL) {
        /* we want to read observation disturbances */
        obsdist = 1;
    }

    for (K->t = tmin; K->t < K->N && !err; K->t += 1) {
        int missobs = 0;

        if (filter_is_varying(K)) {
            err = kalman_refresh_matrices(K, prn);
            if (err) {
                break;
            }
        }

        /* y_t = B'*x_t + Z'*a_t + w_t */
        gretl_matrix_multiply_mod(K->ZT, GRETL_MOD_TRANSPOSE,
                                  K->a0, GRETL_MOD_NONE,
                                  yt, GRETL_MOD_NONE);
        if (K->BT != NULL) {
            missobs = kalman_set_Bx(K);
            if (!missobs) {
                gretl_matrix_add_to(yt, K->Bx);
            }
        }
        if (K->p > 0) {
            /* G \varepsilon_t */
            load_from_row(et, U, K->t, GRETL_MOD_NONE);
            gretl_matrix_multiply(K->G, et, K->v);
        } else if (obsdist) {
            load_from_row_offset(K->v, U, K->t, K->r);
        }
        gretl_matrix_add_to(yt, K->v);

        /* record the t-dated observables */
        load_to_row_offset(Y, yt, K->t, obs_offset);

        /* record the t-dated state? */
        if (S != NULL && tmin == 0) {
            load_to_row_offset(S, K->a0, K->t, 0);
        }

        /* a_{t+1} = T*a_t + v_t */
        gretl_matrix_multiply(K->T, K->a0, K->a1);
        if (K->p > 0) {
            /* H \varepsilon_t */
            gretl_matrix_multiply_mod(K->H, GRETL_MOD_NONE,
                                      et, GRETL_MOD_NONE,
                                      K->a1, GRETL_MOD_CUMULATE);
        } else {
            load_from_row(K->a1, U, K->t, GRETL_MOD_CUMULATE);
        }

        if (K->mu != NULL) {
            gretl_matrix_add_to(K->a1, K->mu);
        }

        /* record the (t+1)-dated state? */
        if (S != NULL && tmin == 1) {
            load_to_row_offset(S, K->a1, K->t, 0);
        }

        fast_copy_values(K->a0, K->a1);
    }

    gretl_matrix_free(yt);
    gretl_matrix_free(et);

    return err;
}

static int check_simul_inputs (kalman *K,
                               const gretl_matrix *U,
                               const gretl_matrix *Sim0,
                               const gretl_matrix *SimX,
                               int ssfsim,
                               PRN *prn)
{
    int err = 0;

    if (U == NULL) {
        err = missing_matrix_error("U");
    } else {
        int ncols;

        if (K->p > 0) {
            /* cross-correlated */
            ncols = K->p;
        } else {
            ncols = K->GG == NULL ? K->r : K->r + K->n;
        }

        if (U->cols != ncols) {
            pprintf(prn, "U should have %d columns but has %d\n",
                    ncols, U->cols);
            err = E_NONCONF;
        }
    }

    if (!err && Sim0 != NULL) {
        int r = ssfsim ? K->r + 1 : K->r;
        int c = ssfsim ? K->r : 1;

        if (Sim0->rows != r || Sim0->cols != c) {
            pprintf(prn, "simstart should be %d x %d, is %d x %d\n",
                    r, c, Sim0->rows, Sim0->cols);
        }
    }

    if (!err && K->x != NULL) {
        /* do we have enough "obsx" data? */
        const gretl_matrix *X = SimX != NULL ? SimX : K->x;

        if (X->rows < U->rows) {
            pprintf(prn, "obsx should have %d rows but has %d\n",
                    U->rows, X->rows);
            err = E_NONCONF;
        }
    }

    return err;
}

gretl_matrix *kalman_bundle_simulate (gretl_bundle *b,
                                      const gretl_matrix *U,
                                      int get_state,
                                      PRN *prn, int *err)
{
    kalman *K = gretl_bundle_get_private_data(b);
    const gretl_matrix *Sim0 = NULL;
    const gretl_matrix *SimX = NULL;
    gretl_matrix *Ret = NULL;
    gretl_matrix *savex = NULL;
    double ssfx;
    int ssfsim = 0;
    int saveN;

    if (K == NULL) {
        *err = E_DATA;
        return NULL;
    }

    /* try accessing auxiliary info from the bundle */
    Sim0 = gretl_bundle_get_matrix(b, "simstart", NULL);
    ssfx = gretl_bundle_get_scalar(b, "ssfsim", NULL);
    if (K->x != NULL) {
        SimX = gretl_bundle_get_matrix(b, "simx", NULL);
    }

    ssfsim = !na(ssfx) && ssfx != 0;

    *err = check_simul_inputs(K, U, Sim0, SimX, ssfsim, prn);
    if (*err) {
        return NULL;
    }

    K->b = b; /* attach bundle pointer */

    saveN = K->N;
    savex = K->x;

    /* we let U temporarily define the sample length */
    K->N = U->rows;

    /* and we allow temporary replacement of K->x */
    if (SimX != NULL) {
        K->x = (gretl_matrix *) SimX;
    }

    /* set state */
    if (ssfsim) {
        K->flags |= (KALMAN_SIM | KALMAN_SSFSIM);
    } else {
        K->flags |= KALMAN_SIM;
    }

    /* now, are the other needed matrices in place? */
    *err = kalman_bundle_recheck_matrices(K, prn);

    /* matrix to hold simulated observables, and state
       if wanted */
    if (!*err) {
        int ncols = get_state ? (K->r + K->n) : K->n;

        Ret = gretl_matrix_alloc(K->N, ncols);
        if (Ret == NULL) {
            *err = E_ALLOC;
        }
    }

    if (!*err) {
        *err = kalman_simulate(K, U, Sim0, Ret, NULL, prn);
    }

    if (*err) {
        gretl_matrix_free(Ret);
        Ret = NULL;
    }

    /* restore state */
    K->flags &= ~KALMAN_SIM;
    K->flags &= ~KALMAN_SSFSIM;
    K->N = saveN;
    K->x = savex;

    return Ret;
}

static int matrix_is_diagonal (const gretl_matrix *m)
{
    double x;
    int i, j;

    for (j=0; j<m->cols; j++) {
        for (i=0; i<m->rows; i++) {
            x = gretl_matrix_get(m, i, j);
            if (i != j && x != 0.0) return 0;
        }
    }

    return 1;
}

static int simdata_refresh_QR (kalman *K, PRN *prn)
{
    int err = 0;

    if (matrix_is_varying(K, K_Q) || matrix_is_varying(K, K_R)) {
        err = kalman_update_matrices(K, prn);
    }

    return err;
}

/* Return a matrix in which the standard normal variates
   in @U are scaled according to K->HH, and K->GG if present.
*/

gretl_matrix *kalman_bundle_simdata (gretl_bundle *b,
                                     const gretl_matrix *U,
                                     PRN *prn, int *err)
{
    kalman *K = gretl_bundle_get_private_data(b);
    gretl_matrix *E = NULL;

    if (K == NULL || U == NULL) {
        *err = E_DATA;
        return NULL;
    }

    if (K->p > 0) {
        if (U->cols != K->p) {
            *err = E_DATA;
            return NULL;
        } else {
            /* nothing to be done, really */
            E = gretl_matrix_copy(U);
        }
    } else {
        int n = K->GG == NULL ? 0 : K->n;
        int t, j, rn = K->r + n;
        int T = U->rows;
        int varying = 0;
        double vjj, utj;

        if (U->cols != rn) {
            *err = E_DATA;
            return NULL;
        }

        E = gretl_matrix_alloc(U->rows, rn);
        if (E == NULL) {
            *err = E_ALLOC;
            return NULL;
        }

        if (matrix_is_varying(K, K_Q) || matrix_is_varying(K, K_R)) {
            varying = 1;
        }

        K->b = b;
        set_kalman_running(K);

        if (matrix_is_diagonal(K->HH) &&
            (K->GG == NULL || matrix_is_diagonal(K->GG))) {
            for (t=0; t<T && !*err; t++) {
                if (varying) {
                    K->t = t;
                    *err = simdata_refresh_QR(K, prn);
                }
                for (j=0; j<rn && !*err; j++) {
                    if (j < K->r) {
                        vjj = gretl_matrix_get(K->HH, j, j);
                    } else {
                        vjj = gretl_matrix_get(K->GG, j-K->r, j-K->r);
                    }
                    utj = gretl_matrix_get(U, t, j);
                    gretl_matrix_set(E, t, j, sqrt(vjj) * utj);
                }
            }
        } else {
            gretl_matrix *V = gretl_zero_matrix_new(rn, rn);
            gretl_matrix *Ut = NULL;
            gretl_matrix *Et = NULL;

            if (V == NULL) {
                *err = E_ALLOC;
                goto bailout;
            }

            if (varying) {
                Ut = gretl_matrix_alloc(1, rn);
                Et = gretl_matrix_alloc(1, rn);
                if (Ut == NULL || Et == NULL) {
                    gretl_matrix_free(V);
                    *err = E_ALLOC;
                    goto bailout;
                }
            }

            if (varying) {
                for (t=0; t<T && !*err; t++) {
                    K->t = t;
                    *err = simdata_refresh_QR(K, prn);
                    if (!*err) {
                        gretl_matrix_inscribe_matrix(V, K->HH, 0, 0,
                                                     GRETL_MOD_NONE);
                        if (n > 0) {
                            gretl_matrix_inscribe_matrix(V, K->GG, K->r, K->r,
                                                         GRETL_MOD_NONE);
                        }
                        *err = gretl_matrix_psd_root(V, 0);
                        if (*err) {
                            gretl_errmsg_set("Failed to compute factor of Omega_t");
                        } else {
                            load_from_row(Ut, U, t, GRETL_MOD_NONE);
                            gretl_matrix_multiply_mod(Ut, GRETL_MOD_NONE,
                                                      V, GRETL_MOD_TRANSPOSE,
                                                      Et, GRETL_MOD_NONE);
                            load_to_row(E, Et, t);
                        }
                    }
                }

                gretl_matrix_free(Ut);
                gretl_matrix_free(Et);
            } else {
                gretl_matrix_inscribe_matrix(V, K->HH, 0, 0,
                                             GRETL_MOD_NONE);
                if (n > 0) {
                    gretl_matrix_inscribe_matrix(V, K->GG, K->r, K->r,
                                                 GRETL_MOD_NONE);
                }
                *err = gretl_matrix_psd_root(V, 0);
                if (*err) {
                    gretl_errmsg_set("Failed to compute factor of Omega");
                } else {
                    gretl_matrix_multiply_mod(U, GRETL_MOD_NONE,
                                              V, GRETL_MOD_TRANSPOSE,
                                              E, GRETL_MOD_NONE);
                }
            }

            gretl_matrix_free(V);
        }
    }

 bailout:

    K->t = 0;
    set_kalman_stopped(K);

    if (E == NULL && !*err) {
        *err = E_ALLOC;
    } else if (E != NULL && *err) {
        gretl_matrix_free(E);
        E = NULL;
    }

    return E;
}

/*
   below: functions to support the "wrapping" of a Kalman
   struct in a gretl bundle
*/

static int check_replacement_dims (const gretl_matrix *orig,
                                   const gretl_matrix *repl,
                                   int id)
{
    int err = 0;

    if (id == K_ZT || id == K_R || id == K_T || id == K_Q ||
        id == K_BT || id == K_m || id == K_a || id == K_P) {
        if (repl->rows != orig->rows || repl->cols != orig->cols) {
            err = E_DATA;
        }
    } else if (id == K_y || id == K_x) {
        if (repl->cols != orig->cols) {
            err = E_DATA;
        }
    }

    if (err) {
        gretl_errmsg_set("You cannot resize a state-space "
                         "system matrix");
    }

    return err;
}

/* On input @targ is a pointer to the matrix to be replaced, if
   there's already a matrix corresponding to @id in place; @src is
   the new matrix, to be copied in if @copy is non-zero or
   otherwise just attached. On output @targ points to the new
   matrix.
*/

static int add_or_replace_k_matrix (kalman *K,
                                    gretl_matrix **targ,
                                    gretl_matrix *src,
                                    int id, int copy)
{
    int err = 0;

    if (*targ != src) {
        if (*targ != NULL) {
            err = check_replacement_dims(*targ, src, id);
            if (err) {
                return err;
            }
            /* destroy old Kalman-owned matrix */
            gretl_matrix_free(*targ);
        }
        if (copy) {
            *targ = gretl_matrix_copy(src);
            if (*targ == NULL) {
                err = E_ALLOC;
            }
        } else {
            *targ = src;
        }
    }

    return err;
}

static gretl_matrix **
get_input_matrix_target_by_id (kalman *K, int i)
{
    gretl_matrix **targ = NULL;

    if (i == K_T) {
        targ = &K->T;
    } else if (i == K_BT) {
        targ = &K->BT;
    } else if (i == K_ZT) {
        targ = &K->ZT;
    } else if (i == K_Q) {
        if (kalman_xcorr(K)) {
            targ = &K->H;
        } else {
            targ = &K->HH;
        }
    } else if (i == K_R) {
        if (kalman_xcorr(K)) {
            targ = &K->G;
        } else {
            targ = &K->GG;
        }
    } else if (i == K_m) {
        targ = &K->mu;
    } else if (i == K_y) {
        targ = &K->y;
    } else if (i == K_x) {
        targ = &K->x;
    } else if (i == K_a) {
        targ = &K->aini;
    } else if (i == K_P) {
        targ = &K->Pini;
    }

    return targ;
}

/* Try attaching a matrix to a Kalman bundle: similar
   to kalman_bundle_set_matrix() below, but allowing
   for the possibility that the @data input of type
   @vtype has to be converted first.
*/

static int
kalman_bundle_try_set_matrix (kalman *K, void *data,
                              GretlType vtype, int id,
                              int copy)
{
    gretl_matrix **targ;
    int err = 0;

    /* determine the location for this matrix */
    targ = get_input_matrix_target_by_id(K, id);

    if (targ == NULL) {
        err = E_DATA;
    } else {
        gretl_matrix *m;

        if (vtype == GRETL_TYPE_MATRIX) {
            m = data;
            err = add_or_replace_k_matrix(K, targ, m, id, copy);
        } else if (vtype == GRETL_TYPE_DOUBLE) {
            m = gretl_matrix_from_scalar(*(double *) data);
            if (m == NULL) {
                err = E_ALLOC;
            } else {
                err = add_or_replace_k_matrix(K, targ, m, id, 0);
            }
        } else {
            err = E_TYPES;
        }
    }

    return err;
}

/* Called by kalman_deserialize() when reconstructing a
   kalman bundle from XML, and also by kalman_bundle_copy()
   when duplicating such a bundle. In these contexts we
   know we have a gretl_matrix on input.
*/

static int
kalman_bundle_set_matrix (kalman *K, gretl_matrix *m,
                          int i, int copy)
{
    gretl_matrix **targ;
    int err = 0;

    targ = get_input_matrix_target_by_id(K, i);

    if (targ == NULL) {
        err = E_DATA;
    } else {
        err = add_or_replace_k_matrix(K, targ, m, i, copy);
    }

    return err;
}

static gretl_matrix **kalman_output_matrix (kalman *K,
                                            const char *key)
{
    gretl_matrix **pm = NULL;

    if (!strcmp(key, "prederr")) {
        pm = &K->V;
    } else if (!strcmp(key, "pevar")) {
        pm = &K->F;
    } else if (!strcmp(key, "state")) {
        pm = &K->A;
    } else if (!strcmp(key, "stvar")) {
        pm = &K->P;
    } else if (!strcmp(key, "gain")) {
        pm = &K->K;
    } else if (!strcmp(key, "llt")) {
        pm = &K->LL;
    } else if (!strcmp(key, "smdist")) {
        pm = &K->U;
    } else if (!strcmp(key, "smdisterr")) {
        pm = &K->Vsd;
    } else if (!strcmp(key, "uhat")) {
        pm = &K->v;
    }

    return pm;
}

#define K_N_OUTPUTS 9

static const char *kalman_output_matrix_names[K_N_OUTPUTS] = {
    "prederr",
    "pevar",
    "state",
    "stvar",
    "gain",
    "llt",
    "smdist",
    "smdisterr",
    "uhat"
};

static int output_matrix_slot (const char *s)
{
    int i;

    for (i=0; i<K_N_OUTPUTS; i++) {
        if (!strcmp(s, kalman_output_matrix_names[i])) {
            return i;
        }
    }

    return -1;
}

#define K_N_SCALARS 10

enum {
    Ks_t = 0,
    Ks_DIFFUSE,
    Ks_CROSS,
    Ks_S2,
    Ks_LNL,
    Ks_r,
    Ks_n,
    Ks_N,
    Ks_p,
    Ks_d
};

static const char *kalman_output_scalar_names[K_N_SCALARS] = {
    "t",
    "diffuse",
    "cross",
    "s2",
    "lnl",
    "r",
    "n",
    "N",
    "p",
    "d"
};

static double *kalman_output_scalar (kalman *K,
                                     const char *key)
{
    /* static storage for on-the-fly scalars */
    static double retval[K_N_SCALARS];
    int i, idx = -1;

    for (i=0; i<K_N_SCALARS; i++) {
        if (!strcmp(key, kalman_output_scalar_names[i])) {
            idx = i;
            break;
        }
    }

    if (idx < 0 && !strcmp(key, "T")) {
        /* backward compatibility */
        idx = Ks_N;
    }

    if (idx < 0) {
        return NULL;
    }

    switch (idx) {
    case Ks_t:
        if (kalman_is_running(K)) {
            retval[idx] = K->t + 1;
        } else {
            retval[idx] = kalman_checking(K) ? 1 : 0;
        }
        break;
    case Ks_DIFFUSE:
        retval[idx] = (K->flags & KALMAN_DIFFUSE)? 1 : 0;
        break;
    case Ks_CROSS:
        retval[idx] = (K->flags & KALMAN_CROSS)? 1 : 0;
        break;
    case Ks_S2:
        retval[idx] = K->s2;
        break;
    case Ks_LNL:
        retval[idx] = K->loglik;
        break;
    case Ks_r:
        retval[idx] = K->r;
        break;
    case Ks_n:
        retval[idx] = K->n;
        break;
    case Ks_N:
        retval[idx] = K->N;
        break;
    case Ks_p:
        retval[idx] = K->p;
        break;
    case Ks_d:
        retval[idx] = K->d;
        break;
    default:
        break;
    }

    return &retval[idx];
}

static const gretl_matrix *k_input_matrix_by_id (kalman *K, int i)
{
    const gretl_matrix *m = NULL;

    if (i == K_T) {
        m = K->T;
    } else if (i == K_BT) {
        m = K->BT;
    } else if (i == K_ZT) {
        m = K->ZT;
    } else if (i == K_Q) {
        if (kalman_xcorr(K)) {
            m = K->H;
        } else {
            m = K->HH;
        }
    } else if (i == K_R) {
        if (kalman_xcorr(K)) {
            m = K->G;
        } else {
            m = K->GG;
        }
    } else if (i == K_m) {
        m = K->mu;
    } else if (i == K_y) {
        m = K->y;
    } else if (i == K_x) {
        m = K->x;
    } else if (i == K_a) {
        m = K->aini;
    } else if (i == K_P) {
        m = K->Pini;
    }

    return m;
}

static int input_matrix_slot (const char *s)
{
    int i;

    for (i=0; i<K_MMAX; i++) {
        if (!strcmp(s, K_input_mats[i].name)) {
            return K_input_mats[i].sym;
        }
    }

    return -1;
};

static GretlType kalman_extra_type (const char *key)
{
    if (!strcmp(key, "ssfsim")) {
        return GRETL_TYPE_DOUBLE;
    } else if (!strcmp(key, "simstart") ||
               !strcmp(key, "simx")) {
        return GRETL_TYPE_MATRIX;
    } else {
        return GRETL_TYPE_NONE;
    }
}

/* respond to "diffuse" setting via bundle apparatus */

static int kalman_set_diffuse (kalman *K, int d)
{
    if (d) {
	K->exact = (d == 2);
	K->flags |= KALMAN_DIFFUSE;
	return diffuse_Pini(K);
    } else {
	K->exact = 0;
	K->flags &= ~KALMAN_DIFFUSE;
	if (K->Pk0 != NULL) {
	    gretl_matrix_free(K->Pk0);
	    gretl_matrix_free(K->Pk1);
	    gretl_matrix_free(K->Fk);
	    K->Pk0 = K->Pk1 = K->Fk = NULL;
	}
	return 0;
    }
}

/* Called by real_bundle_set_data() in gretl_bundle.c.
   The return value indicates whether the putative
   setting was handled (1) or not (0). Not being
   handled here is not necessarily an error.

   The @copy flag here is inherited from the specific
   caller of real_bundle_set_data(): @copy = 1 if the
   caller was gretl_bundle_set_data(), 0 if it was
   gretl_bundle_donate_data(). Either way the kalman
   struct takes ownership.
*/

int maybe_set_kalman_element (void *kptr,
			      const char *key,
			      void *vptr,
			      GretlType vtype,
			      int copy,
			      int *err)
{
    GretlType targtype;
    kalman *K = kptr;
    int fncall = 0;
    int i, id = -1;
    int done = 0;

    if (K == NULL) {
	*err = E_DATA;
	return 0;
    }

    /* Check for optional "extra" kalman items that
       live outside of the kalman struct itself.
    */
    targtype = kalman_extra_type(key);
    if (targtype != GRETL_TYPE_NONE) {
	if (vtype != targtype) {
	    *err = E_TYPES;
	}
	return 0;
    }

    if (!strcmp(key, "diffuse")) {
	if (vtype == GRETL_TYPE_DOUBLE) {
	    double d = *(double *) vptr;

	    if (d != 0 && d != 1 && d != 2) {
		*err = E_INVARG;
	    } else {
		*err = kalman_set_diffuse(K, (int) d);
	    }
	} else {
	    *err = E_TYPES;
	}
	return 1; /* done, error or no */
    }

    if (!strcmp(key, "timevar_call")) {
        /* try for a function call specifier (string) */
        if (vtype == GRETL_TYPE_STRING) {
            fncall = 1;
        } else {
            *err = E_TYPES;
        }
    } else {
        /* try for a matrix specifier */
        i = input_matrix_slot(key);
        if (i >= 0) {
            if (vtype == GRETL_TYPE_MATRIX ||
                vtype == GRETL_TYPE_DOUBLE) {
                id = i;
            } else {
                *err = E_TYPES;
            }
        }
    }

    if (*err) {
        return 0;
    } else if (fncall) {
        if (copy) {
            K->matcall = gretl_strdup((char *) vptr);
        } else {
            K->matcall = (char *) vptr;
        }
        /* re-evaluate what's actually varying */
        *err = check_for_matrix_updates(K, NULL);
        if (!*err) {
            done = 1;
        }
    } else if (id < 0) {
        if (kalman_output_matrix(K, key) != NULL ||
            kalman_output_scalar(K, key) != NULL) {
            *err = E_DATA;
            gretl_errmsg_sprintf("The member %s is read-only", key);
        }
    } else {
        *err = kalman_bundle_try_set_matrix(K, vptr, vtype, id, copy);
        done = (*err == 0);
    }

    return done;
}

int maybe_delete_kalman_element (void *kptr,
                                 const char *key,
                                 int *err)
{
    kalman *K = kptr;
    gretl_matrix **pm;
    int done = 0;

    if (K == NULL) {
        return 0;
    }

    if (kalman_output_scalar(K, key) != NULL ||
        input_matrix_slot(key) >= 0 || !strcmp(key, "uhat")) {
        /* note: the matrix under the key "uhat" is part of
           the internal kalman apparatus */
        gretl_errmsg_sprintf("%s: cannot be deleted", key);
        *err = E_DATA;
    } else if ((pm = kalman_output_matrix(K, key)) != NULL) {
        /* OK to delete a user-output matrix */
        gretl_matrix_free(*pm);
        *pm = NULL;
    } else if (!strcmp(key, "timevar_call")) {
        /* OK to delete time-variation call */
        if (K->matcall != NULL) {
            free(K->matcall);
            K->matcall = NULL;
            free(K->varying);
            K->varying = NULL;
            done = 1;
        } else {
            *err = E_DATA;
        }
    }

    return done;
}

void *maybe_retrieve_kalman_element (void *kptr,
                                     const char *key,
                                     GretlType *type,
                                     int *reserved,
                                     int *err)
{
    kalman *K = kptr;
    void *ret = NULL;
    int i, id = -1;

    *type = GRETL_TYPE_NONE;

    if (K == NULL) {
        *err = E_DATA;
        return NULL;
    }

    if (!strcmp(key, "timevar_call")) {
        /* function call specifier? */
        *reserved = 1;
        if (K->matcall != NULL) {
            ret = K->matcall;
            *type = GRETL_TYPE_STRING;
        }
    } else {
        /* try for an input matrix specifier */
        for (i=0; i<K_MMAX; i++) {
            if (!strcmp(key, K_input_mats[i].name)) {
                id = K_input_mats[i].sym;
                ret = (gretl_matrix *) k_input_matrix_by_id(K, id);
                if (ret != NULL) {
                    *type = GRETL_TYPE_MATRIX;
                }
                break;
            }
        }
        if (id < 0) {
            /* try for an output matrix */
            gretl_matrix **pm = kalman_output_matrix(K, key);

            if (pm != NULL) {
                *reserved = 1;
                if (*pm != NULL) {
                    ret = *pm;
                    *type = GRETL_TYPE_MATRIX;
                }
            }
        }
        if (id < 0 && *reserved == 0) {
            /* nothing matched yet: try scalar member */
            ret = kalman_output_scalar(K, key);
            if (ret != NULL) {
                *type = GRETL_TYPE_DOUBLE;
            }
        }
    }

    if (id >= 0 && *reserved == 0) {
        /* flag the fact that @key was a kalman-reserved
           identifier */
        *reserved = 1;
    }

    if (*reserved && ret == NULL) {
        gretl_errmsg_sprintf("\"%s\": %s", key, _("no such item"));
        *err = E_DATA;
    }

    return ret;
}

static int output_matrix_count (kalman *K)
{
    gretl_matrix **pm;
    int i, n = 0;

    for (i=0; i<K_N_OUTPUTS; i++) {
        pm = kalman_output_matrix(K, kalman_output_matrix_names[i]);
        n += (pm != NULL && *pm != NULL);
    }

    return n;
}

int print_kalman_bundle_info (void *kptr, PRN *prn)
{
    kalman *K = kptr;
    int err = 0;

    if (K == NULL) {
        pputs(prn, "Kalman struct: empty\n");
        err = E_DATA;
    } else {
        const gretl_matrix *m;
        gretl_matrix **pm;
        double *px;
        const char *name;
        int i, id;

        pputs(prn, "\nKalman input matrices\n");

        for (i=0; i<K_MMAX; i++) {
            id = K_input_mats[i].sym;
            m = k_input_matrix_by_id(K, id);
            if (m != NULL) {
                pprintf(prn, " %s: ", K_input_mats[i].name);
                pprintf(prn, "%d x %d\n", m->rows, m->cols);
            }
        }

        if (output_matrix_count(K) > 0) {
            pputs(prn, "\nKalman output matrices\n");
            for (i=0; i<K_N_OUTPUTS; i++) {
                name = kalman_output_matrix_names[i];
                pm = kalman_output_matrix(K, name);
                if (pm != NULL && *pm != NULL) {
                    m = *pm;
                    pprintf(prn, " %s: ", name);
                    pprintf(prn, "%d x %d\n", m->rows, m->cols);
                }
            }
        }

        pputs(prn, "\nKalman scalars\n");

        for (i=0; i<K_N_SCALARS; i++) {
            name = kalman_output_scalar_names[i];
            pprintf(prn, " %s: ", name);
            px = kalman_output_scalar(K, name);
            if (px == NULL || na(*px)) {
                pputs(prn, "NA\n");
            } else {
                pprintf(prn, "%g\n", *px);
            }
        }

        if (K->matcall != NULL) {
            pputs(prn, "\nKalman strings\n");
            pprintf(prn, " timevar_call: %s\n", K->matcall);
        }
    }

    return err;
}

/* for use in context of a kalman bundle: serialize the
   information in the kalman struct to XML
*/

int kalman_serialize (void *kptr, PRN *prn)
{
    kalman *K = kptr;
    const gretl_matrix *m;
    gretl_matrix **pm;
    double *px;
    const char *name;
    int i, err = 0;

    if (K == NULL) {
        fputs("kalman_serialize: got NULL\n", stderr);
        return E_DATA;
    }

    pputs(prn, "<gretl-kalman>\n");

    for (i=0; i<K_MMAX; i++) {
        m = k_input_matrix_by_id(K, K_input_mats[i].sym);
        if (m != NULL) {
            gretl_matrix_serialize(m, K_input_mats[i].name, prn);
        }
    }

    for (i=0; i<K_N_OUTPUTS; i++) {
        name = kalman_output_matrix_names[i];
        pm = kalman_output_matrix(K, name);
        if (pm != NULL && *pm != NULL) {
            gretl_matrix_serialize(*pm, name, prn);
        }
    }

    for (i=0; i<=Ks_LNL; i++) {
        name = kalman_output_scalar_names[i];
        px = kalman_output_scalar(K, name);
        if (px != NULL && !na(*px)) {
            gretl_finite_scalar_serialize(*px, name, prn);
        }
    }

    if (K->matcall != NULL) {
        gretl_string_serialize(K->matcall, "timevar_call", prn);
    }

    pputs(prn, "</gretl-kalman>\n");

    return err;
}

static int required_matrix_slot (const char *s)
{
    if (!strcmp(s, "obsy"))     return 0;
    if (!strcmp(s, "obsymat"))  return 1;
    if (!strcmp(s, "statemat")) return 2;
    if (!strcmp(s, "statevar")) return 3;
    if (!strcmp(s, "obsvar"))   return 4;
    return -1;
};

/* for use in context of a kalman bundle: deserialize the
   kalman struct from XML
*/

gretl_bundle *kalman_deserialize (void *p1, void *p2, int *err)
{
    xmlNodePtr cur, node = p1;
    xmlDocPtr doc = p2;
    gretl_matrix *Mreq[5] = {NULL};
    gretl_matrix *Mopt[K_MMAX] = {NULL};
    gretl_matrix *Mout[K_N_OUTPUTS] = {NULL};
    char *tvcall = NULL;
    double s2 = NADBL;
    double lnl = NADBL;
    int copy[5] = {0};
    int i, nmats = 0;
    int Kflags = 0;
    gretl_matrix *m;
    double x;
    char *key, *strv;
    gretl_bundle *b = NULL;

    while (node != NULL && !*err) {
        if (!xmlStrcmp(node->name, (XUC) "gretl-kalman")) {
            cur = node->xmlChildrenNode;
            while (cur != NULL && !*err) {
                key = (char *) xmlGetProp(cur, (XUC) "name");
                if (!xmlStrcmp(cur->name, (XUC) "gretl-matrix")) {
                    /* pick up kalman matrices */
                    m = gretl_xml_get_matrix(cur, doc, err);
                    if ((i = required_matrix_slot(key)) >= 0) {
                        nmats++;
                        Mreq[i] = m;
                    } else if ((i = input_matrix_slot(key)) >= 0) {
                        Mopt[i] = m;
                    } else if ((i = output_matrix_slot(key)) >= 0) {
                        Mout[i] = m;
                    }
                } else if (!xmlStrcmp(cur->name, (XUC) "scalar")) {
                    /* pick up kalman scalars */
                    if (gretl_xml_get_prop_as_double(cur, "value", &x)) {
                        if (!strcmp(key, "diffuse") && x > 0) {
                            Kflags |= KALMAN_DIFFUSE;
                        } else if (!strcmp(key, "cross") && x > 0) {
                            Kflags |= KALMAN_CROSS;
                        } else if (!strcmp(key, "s2")) {
                            s2 = x;
                        } else if (!strcmp(key, "lnl")) {
                            lnl = x;
                        }
                    }
                } else if (!xmlStrcmp(cur->name, (XUC) "string")) {
                    /* pick up kalman strings */
                    if (!strcmp(key, "timevar_call") &&
                        gretl_xml_get_prop_as_string(cur, "value", &strv)) {
                            tvcall = strv;
                    }
                }
                free(key);
                cur = cur->next;
            }
            break;
        }
        node = node->next;
    }

    if (nmats == 5 && !(Kflags & KALMAN_CROSS)) {
        /* drop obsvar from initialization */
        Mopt[K_R] = Mreq[4];
        Mreq[4] = NULL;
        nmats--;
    }

    if (((Kflags & KALMAN_CROSS) && nmats != 5) ||
        (!(Kflags & KALMAN_CROSS) && nmats != 4)) {
        *err = E_DATA;
    } else {
        b = kalman_bundle_new(Mreq, copy, nmats, err);
        if (b != NULL) {
            kalman *K = gretl_bundle_get_private_data(b);
            gretl_matrix **pm;
            const char *name;

            K->flags = Kflags;
            K->s2 = s2;
            K->loglik = lnl;

            for (i=0; i<K_MMAX; i++) {
                if (Mopt[i] != NULL) {
                    kalman_bundle_set_matrix(K, Mopt[i], i, 0);
                }
            }
            for (i=0; i<K_N_OUTPUTS; i++) {
                if (Mout[i] != NULL) {
                    name = kalman_output_matrix_names[i];
                    pm = kalman_output_matrix(K, name);
                    *pm = Mout[i];
                }
            }
            K->matcall = tvcall;
        }
    }

    if (*err) {
        /* clean up */
        for (i=0; i<5; i++) {
            gretl_matrix_free(Mreq[i]);
        }
        for (i=0; i<K_MMAX; i++) {
            gretl_matrix_free(Mopt[i]);
        }
        for (i=0; i<K_N_OUTPUTS; i++) {
            gretl_matrix_free(Mout[i]);
        }
        free(tvcall);
    }

    return b;
}

/* Called from gretl_bundle.c to meet the case where the user calls
   for a kalman bundle to be copied: here we create a new kalman
   struct and copy across the required elements (since they are
   not regular bundle members).
*/

gretl_bundle *kalman_bundle_copy (const gretl_bundle *src, int *err)
{
    kalman *K, *Knew;
    gretl_bundle *b = NULL;
    gretl_matrix *M[5] = {NULL};
    int copy[5] = {1, 1, 1, 1, 1};
    gretl_matrix *m, **pm, **pm1;
    const char *name;
    int i, id, k = 4;

    K = gretl_bundle_get_private_data((gretl_bundle *) src);

    if (K == NULL) {
        *err = E_DATA;
        return NULL;
    }

    /* set pointers to required Kalman matrices */
    M[0] = K->y;
    M[1] = K->ZT;
    M[2] = K->T;

    /* variants dependent on presence/absence of cross-correlation */
    if (kalman_xcorr(K)) {
        M[3] = K->H;
        M[4] = K->G;
        k = 5;
    } else {
        M[3] = K->HH;
    }

    b = kalman_bundle_new(M, copy, k, err);

    if (*err) {
        return b;
    }

    Knew = gretl_bundle_get_private_data(b);
    Knew->flags = K->flags;

    /* add any "extra" matrices, beyond the required ones */
    for (i=0; i<n_extra_mats && !*err; i++) {
        id = extra_mats[i];
        m = (gretl_matrix *) k_input_matrix_by_id(K, id);
        if (m != NULL) {
            *err = kalman_bundle_set_matrix(Knew, m, id, 1);
        }
    }

    for (i=0; i<K_N_OUTPUTS && !*err; i++) {
        name = kalman_output_matrix_names[i];
        pm = kalman_output_matrix(K, name);
        if (pm != NULL && *pm != NULL) {
            pm1 = kalman_output_matrix(Knew, name);
            *pm1 = gretl_matrix_copy(*pm);
        }
    }

    Knew->ifc = K->ifc;
    Knew->s2 = K->s2;
    Knew->loglik = K->loglik;

    if (K->flags & KALMAN_CROSS) {
        Knew->flags |= KALMAN_CROSS;
    }

    if (K->flags & KALMAN_DIFFUSE) {
        Knew->flags |= KALMAN_DIFFUSE;
    }

    if (K->matcall != NULL) {
        Knew->matcall = gretl_strdup(K->matcall);
    }

    return b;
}

/* for use in constructing GUI bundle save menu */

char **kalman_bundle_get_matrix_names (kalman *K, int *ns)
{
    char **S = NULL;
    gretl_matrix **pm;
    const char *name;
    int i, id, err = 0;

    *ns = 0;

    for (i=0; i<K_MMAX && !err; i++) {
        id = K_input_mats[i].sym;
        if (k_input_matrix_by_id(K, id) != NULL) {
            err = strings_array_add(&S, ns, K_input_mats[i].name);
        }
    }

    for (i=0; i<K_N_OUTPUTS && !err; i++) {
        name = kalman_output_matrix_names[i];
        pm = kalman_output_matrix(K, name);
        if (pm != NULL && *pm != NULL) {
            err = strings_array_add(&S, ns, name);
        }
    }

    return S;
}

/* also for use in constructing GUI bundle save menu */

char **kalman_bundle_get_scalar_names (kalman *K, int *ns)
{
    char **S;

    *ns = K_N_SCALARS -1 - na(K->s2) - na(K->loglik);
     S = strings_array_new(*ns);

    if (S != NULL) {
        int i = 0;

        /* flags */
        S[i++] = gretl_strdup("cross");
        S[i++] = gretl_strdup("diffuse");

        /* actual numerical outputs */
        if (!na(K->s2)) {
            S[i++] = gretl_strdup("s2");
        }
        if (!na(K->loglik)) {
            S[i++] = gretl_strdup("lnl");
        }

        /* system dimensions */
        S[i++] = gretl_strdup("r");
        S[i++] = gretl_strdup("n");
        S[i++] = gretl_strdup("N");
        S[i++] = gretl_strdup("p");
    }

    return S;
}

/* to support the nelem() function for kalman bundles */

int kalman_bundle_n_members (gretl_bundle *b)
{
    kalman *K = gretl_bundle_get_private_data(b);
    int n = 0;

    if (K != NULL) {
        int i, id;

        n = K_N_SCALARS;

        for (i=0; i<K_MMAX; i++) {
            id = K_input_mats[i].sym;
            n += (k_input_matrix_by_id(K, id) != NULL);
        }

        n += output_matrix_count(K);
    }

    return n;
}
