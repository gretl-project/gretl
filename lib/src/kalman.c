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
#include "gretl_bfgs.h"
#include "gretl_typemap.h"
#include "kalman.h"

/**
 * SECTION:kalman
 * @short_description: The Kalman filter
 * @title: Kalman
 * @include: gretl/libgretl.h, gretl/kalman.h
 *
 */

#define K_TINY 1.0e-12
#define P0_ZERO 1
#define SUPPORT_LEGACY 1

static int kdebug;
static int trace;
static int identify;

/* temporary testing thing */
static int exact_set;

enum {
    KALMAN_BUNDLE  = 1 << 0,  /* user-defined filter, inside bundle */
    KALMAN_DIFFUSE = 1 << 1,  /* using diffuse P_{1|0} */
    KALMAN_FORWARD = 1 << 2,  /* running forward filtering pass */
    KALMAN_SIM     = 1 << 3,  /* running simulation */
    KALMAN_CHECK   = 1 << 4,  /* checking user-defined matrices */
    KALMAN_SSFSIM  = 1 << 5,  /* on simulation, emulate SsfPack */
    KALMAN_ARMA_LL = 1 << 6,  /* filtering for ARMA estimation */
    KALMAN_EXTRA   = 1 << 7   /* recording extra info (att, Ptt) */
};

/* How exactly are the disturbance variances represented? */
enum {
    STD_VAR, /* standard, basic variance representation: VS, VY */
    DJ_VAR,  /* as per DeJong, supporting cross correlation */
    DK_VAR   /* as per Durbin and Koopman, with selection */
};

/* On filtering pass, what are we preparing for? */
enum {
    SM_NONE = 0, /* not preparing for smoothing */
    SM_STATE,    /* preparing for state smoothing */
    SM_DIST      /* preparing for disturbance smoothing */
};

/* Which filtering and smoothing code are we using? */
enum {
    K_LEGACY, /* legacy Kalman code */
    K_UNIVAR, /* univariate approach, a la KFAS */
    K_DEJONG  /* de Jong and Lin apporach */
};

typedef struct stepinfo_ stepinfo;

/* recorder for time-varying T and/or Z */
struct stepinfo_ {
    gretl_matrix *T;   /* (r * r) x N */
    gretl_matrix *ZT;  /* (r * n) x N */
};

typedef struct univar_info_ univar_info;

/* holder for data specific to univariate representation */
struct univar_info_ {
    gretl_matrix *g;    /* obs disturbance variance vector */
    gretl_matrix *Z;    /* untransposed Z matrix */
    gretl_matrix *Zi;   /* row i of Z */
    gretl_matrix *Kti;  /* column i of Kt */
    gretl_matrix *Pki;  /* P at step i */
    gretl_matrix *Kki;  /* K∞ at i */
    gretl_matrix *m1;   /* workspace */
    gretl_matrix *mm;   /* workspace */
    gretl_matrix *Linv; /* for diagonalization */
    gretl_matrix *Finf; /* recorder for smoothing */
    gretl_matrix *Kinf; /* ditto */
    gretl_matrix *y0;   /* original data (when diagonalizing) */
    gretl_matrix *UF;   /* univariate counterpart of F */
    int free_g;
    int free_Z;
    int Nd;
};

typedef struct dejong_info_ dejong_info;

/* holder for data specific to de Jong (exact) representation */
struct dejong_info_ {
    gretl_matrix *Vt;  /* for augmented recursion */
    gretl_matrix *At;  /* ditto */
    gretl_matrix *Vv;  /* ditto */
    gretl_matrix *Qt;  /* ditto */
    gretl_matrix *S;   /* north-west block of Qt */
    gretl_matrix *s;   /* vector embedded in Qt */
    gretl_matrix *A;   /* recorder for smoother */
    gretl_matrix *V;   /* ditto */
    gretl_matrix *Am;  /* ditto (redundant?) */
    double qm;         /* scalar embedded in Qt */
    double ldS;        /* log determinant of S-inverse */
    int Nd;            /* size of diffuse recorders */
};

struct kalman_ {
    int flags;   /* for recording any options */
    int code;    /* code collection used: K_LEGACY, K_UNIVAR or K_DEJONG */
    int exact;   /* exact initial iterations? */
    int vartype; /* STD_VAR, DJ_VAR or DK_VAR */

    int r;   /* rows of a = number of elements in state */
    int n;   /* columns of y = number of observables */
    int k;   /* columns of B = number of exogenous vars in obs eqn */
    int p;   /* length of combined disturbance vector */
    int q;   /* length of state disturbance vector */
    int N;   /* rows of y = number of observations */
    int okN; /* N - number of missing observations */
    int t;   /* current time step, when filtering */
    int d;   /* (diffuse) time-step at which standard iterations start */
    int j;   /* (univariate) observable at which standard iters start */
    int TI;  /* flag for K->T is identity matrix */

    int ifc; /* boolean: obs equation includes an implicit constant? */
    int dj_initted; /* flag for basic de Jong initialization done */
    int smo_prep;   /* preparing for smoothing (of type) */

    double SSRw;    /* \sum_{t=1}^N v_t^{\prime} F_t^{-1} v_t */
    double loglik;  /* log-likelihood */
    double s2;      /* = SSRw / k */

    /* continuously updated matrices */
    gretl_matrix *a0; /* r x 1: state vector, before updating */
    gretl_matrix *a1; /* r x 1: state vector, after updating */
    gretl_matrix *P0; /* r x r: MSE matrix, before updating */
    gretl_matrix *P1; /* r x r: MSE matrix, after updating */
    gretl_matrix *vt; /* n x 1: one-step forecast error(s), time t */

    /* for the exact diffuse univariate case */
    gretl_matrix *Pk0; /* P∞,t */
    gretl_matrix *PK;  /* P∞, all t */

    /* for the de Jong variant */
    gretl_matrix *Lt; /* r x r */
    gretl_matrix *Jt; /* r x p */
    gretl_matrix *rr; /* r x r */
    gretl_matrix *L;  /* recorder */
    gretl_matrix *J;  /* recorder */
    gretl_matrix *VS0; /* baseline for detecting change */
    gretl_matrix *VY0; /* ditto */

    /* input data matrices: note that the order matters for various
       functions, including matrix_is_varying()
    */
    gretl_matrix *T;  /* r x r: state transition matrix */
    gretl_matrix *BT; /* k x n: B', coeffs on exogenous vars, obs eqn */
    gretl_matrix *ZT; /* r x n: Z', coeffs on state variables, obs eqn */
    gretl_matrix *VS; /* r x r: contemp covariance matrix, state eqn */
    gretl_matrix *VY; /* n x n: contemp covariance matrix, obs eqn */
    gretl_matrix *mu; /* r x 1: constant term in state transition */
    gretl_matrix *y;  /* N x n: dependent variable vector (or matrix) */
    gretl_matrix *x;  /* N x k: independent variables matrix */
    gretl_matrix *aini; /* r x 1: a_0 */
    gretl_matrix *Pini; /* r x r: P_0 */

    /* user inputs for cross-correlated disturbances */
    gretl_matrix *H;  /* r x p: state var factor */
    gretl_matrix *G;  /* n x p: obs var factor */
    /* and the cross-matrix itself */
    gretl_matrix *HG; /* r x n */

    /* OR: Durbin-Koopman R and Q (as in statevar = RQR') */
    gretl_matrix *Q;   /* q x q, symmetric */
    gretl_matrix *R;   /* r x q, selection matrix */
    gretl_matrix *QRT; /* Q*R' */

    /* apparatus for registering time-variation of matrices */
    char *matcall;
    char *varying;

    /* optional matrices for recording extra info */
    gretl_matrix *LL;  /* N x 1: loglikelihood, all time-steps */

    /* optional run-time export matrices */
    gretl_matrix *V;   /* N x n: forecast errors, all time-steps */
    gretl_matrix *F;   /* N x n(n+1)/2: MSE for observables, all time-steps */
    gretl_matrix *A;   /* N x r: state vector, all time-steps */
    gretl_matrix *P;   /* N x r(r+1)/2: MSE for state, all time-steps */
    gretl_matrix *K;   /* N x rn: gain matrix, all time-steps */
    gretl_matrix *U;   /* N x r+n: smoothed disturbances */
    gretl_matrix *Vsd; /* Variance of smoothed disturbance */

    /* struct needed only when smoothing in the time-varying case */
    stepinfo *step;

    /* struct needed for univariate approach */
    univar_info *uinfo;

    /* struct needed for de Jong exact approach */
    dejong_info *djinfo;

    /* workspace matrices */
    gretl_matrix_block *Blk; /* holder for the following */
    gretl_matrix *PZ;
    gretl_matrix *Ft;
    gretl_matrix *iFt;
    gretl_matrix *Kt;
    gretl_matrix *Mt;
    gretl_matrix *Ct;

    gretl_bundle *b; /* the bundle of which this struct is a member */
    void *data;      /* handle for attaching additional info */
    PRN *prn;        /* verbose printer */
};

enum {
    K_V,
    K_F,
    K_A,
    K_BIG_P,
    K_K,
    K_LL,
};

/* max number of time-varying matrices: T, B, Z, VS, VY, mu */
#define K_N_MATCALLS 6

#define set_kalman_running(K) (K->flags |= KALMAN_FORWARD)
#define kalman_is_running(K)  (K->flags & KALMAN_FORWARD)
#define kalman_simulating(K)  (K->flags & KALMAN_SIM)
#define kalman_checking(K)    (K->flags & KALMAN_CHECK)
#define kalman_ssfsim(K)      (K->flags & KALMAN_SSFSIM)
#define kalman_diffuse(K)     ((K->flags & KALMAN_DIFFUSE)? 1 : 0)
#define kalman_arma_ll(K)     (K->flags & KALMAN_ARMA_LL)
#define kalman_extra(K)       (K->flags & KALMAN_EXTRA)
#define kalman_is_bundle(K)   ((K->flags & KALMAN_BUNDLE)? 1 : 0)

#define kalman_univariate(K)  (K->code == K_UNIVAR)
#define kalman_dejong(K)      (K->code == K_DEJONG)

#define kalman_djvar(K)       (K->vartype == DJ_VAR)
#define kalman_dkvar(K)       (K->vartype == DK_VAR)

#define filter_is_varying(K) (K->matcall != NULL)

static const char *kalman_matrix_name (int sym);
static int kalman_revise_variance (kalman *K);
static int check_for_matrix_updates (kalman *K, ufunc *uf);
static int set_initial_statevar (kalman *K);
static int kalman_set_diffuse (kalman *K, int d);
static int kfilter_univariate (kalman *K, PRN *prn);
static int kfilter_dejong (kalman *K, PRN *prn);
static int ksmooth_univariate (kalman *K, int dist);
static int state_smooth_dejong (kalman *K);
static int dist_smooth_dejong (kalman *K, int DKstyle);
static int kalman_add_univar_info (kalman *K);
static int kalman_add_dejong_info (kalman *K);
static int bundle_add_matrix (gretl_bundle *b,
                              const char *key,
                              gretl_matrix *m);
#if SUPPORT_LEGACY
static int anderson_moore_smooth (kalman *K);
static int koopman_smooth (kalman *K, int DKstyle);
#endif

int is_kalman_bundle (gretl_bundle *b)
{
    kalman *K = gretl_bundle_get_private_data(b);

    return K != NULL;
}

/* symbolic identifiers for input matrices: note that potentially
   time-varying matrices must appear first in the enumeration, and
   the order must match the order of the "input data matrices" in
   the Kalman struct (above).
*/

enum {
    K_T = 0,
    K_BT,
    K_ZT,
    K_VS,
    K_VY,
    K_m,
    K_y,
    K_x,
    K_a,
    K_P,
    K_R,
    K_MMAX /* sentinel */
};

/* variant of gretl_matrix_copy_values for us when we already
   know that the matrices are non-NULL and conformable
*/

static inline void fast_copy_values (gretl_matrix *B, const gretl_matrix *A)
{
    memcpy(B->val, A->val, B->rows * B->cols * sizeof(double));
}

static void fast_write_I (gretl_matrix *A)
{
    int i, j, k = 0;

    for (j=0; j<A->cols; j++) {
        for (i=0; i<A->rows; i++) {
            A->val[k++] = (i == j)? 1 : 0;
        }
    }
}

static int matrix_changed (const gretl_matrix *A, const gretl_matrix *B)
{
    int i, n = A->rows * A->cols;

    for (i=0; i<n; i++) {
	if (A->val[i] != B->val[i]) {
	    return 1;
	}
    }

    return 0;
}

static void set_kalman_stopped (kalman *K)
{
    K->flags &= ~KALMAN_FORWARD;
    K->t = 0;
}

static void kalman_check_env (kalman *K)
{
    char *s1 = getenv("KALMAN_DEBUG");
    char *s2 = getenv("KALMAN_TRACE");

    if (s1 != NULL) kdebug = atoi(s1);
    if (s2 != NULL) trace  = atoi(s2);

    if (getenv("KALMAN_EXACT")) {
	exact_set = 1;
    }

    s1 = getenv("KALMAN_UNIVAR");
    if (s1 != NULL) {
        if (K->vartype == DJ_VAR) {
            fprintf(stderr, "KALMAN_UNIVAR not applicable!\n");
        } else {
            K->code = K_UNIVAR;
        }
    } else {
        s2 = getenv("KALMAN_DEJONG");
        if (s2 != NULL) {
            K->code = K_DEJONG;
        }
    }
}

static void free_stepinfo (kalman *K)
{
    if (K->step != NULL) {
        gretl_matrix_free(K->step->T);
        gretl_matrix_free(K->step->ZT);
        free(K->step);
        K->step = NULL;
    }
}

static void free_univar_info (kalman *K)
{
    if (K->uinfo != NULL) {
        if (K->uinfo->free_g) {
            gretl_matrix_free(K->uinfo->g);
        }
        if (K->uinfo->free_Z) {
            gretl_matrix_free(K->uinfo->Z);
        }
        gretl_matrix_free(K->uinfo->Zi);
        gretl_matrix_free(K->uinfo->Kti);
        gretl_matrix_free(K->uinfo->Pki);
        gretl_matrix_free(K->uinfo->Kki);
        gretl_matrix_free(K->uinfo->m1);
        gretl_matrix_free(K->uinfo->mm);
        gretl_matrix_free(K->uinfo->Linv);
        gretl_matrix_free(K->uinfo->Finf);
        gretl_matrix_free(K->uinfo->Kinf);
        gretl_matrix_free(K->uinfo->y0);
        gretl_matrix_free(K->uinfo->UF);

        free(K->uinfo);
        K->uinfo = NULL;
    }
}

static void free_dejong_info (kalman *K)
{
    if (K->djinfo != NULL) {
        gretl_matrix_free(K->djinfo->Vt);
        gretl_matrix_free(K->djinfo->At);
        gretl_matrix_free(K->djinfo->Vv);
        gretl_matrix_free(K->djinfo->Qt);
        gretl_matrix_free(K->djinfo->S);
        gretl_matrix_free(K->djinfo->s);

        gretl_matrix_free(K->djinfo->A);
        gretl_matrix_free(K->djinfo->V);
        gretl_matrix_free(K->djinfo->Am);

        free(K->djinfo);
        K->djinfo = NULL;
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
    gretl_matrix_free(K->vt);
    gretl_matrix_free(K->LL);
    gretl_matrix_free(K->Pk0);
    gretl_matrix_free(K->L);
    gretl_matrix_free(K->J);

    /* internally allocated workspace */
    gretl_matrix_block_destroy(K->Blk);

    if (kalman_is_bundle(K)) {
        gretl_matrix **mptr[] = {
            &K->T, &K->BT, &K->ZT, &K->VS, &K->VY,
            &K->mu, &K->y, &K->x, &K->aini, &K->Pini,
            &K->R
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

    if (!kalman_univariate(K)) {
        gretl_matrix_free(K->H);
        gretl_matrix_free(K->G);
        gretl_matrix_free(K->HG);
        gretl_matrix_free(K->Jt);
	gretl_matrix_free(K->VS0);
	gretl_matrix_free(K->VY0);
    }

    if (kalman_dkvar(K)) {
        /* Durbin-Koopman errors info */
        gretl_matrix_free(K->Q);
        gretl_matrix_free(K->QRT);
    }

    if (K->step != NULL) {
        free_stepinfo(K);
    }
    if (K->uinfo != NULL) {
        free_univar_info(K);
    } else if (K->djinfo != NULL) {
        free_dejong_info(K);
    }

    free(K);
}

static kalman *kalman_new_empty (int flags)
{
    kalman *K = malloc(sizeof *K);

    if (K != NULL) {
        K->exact = 0;
	K->code = K_LEGACY;
        K->vartype = STD_VAR;
        K->aini = K->Pini = NULL;
        K->a0 = K->a1 = NULL;
        K->P0 = K->P1 = NULL;
        K->Pk0 = NULL;
        K->PK = NULL;
        K->LL = NULL;
        K->vt = NULL;
        K->Blk = NULL;
        K->T = K->BT = K->ZT = NULL;
        K->VS = K->VY = NULL;
        K->H = K->G = K->HG = K->Jt = NULL;  /* de Jong variance */
        K->Q = K->R = K->QRT = NULL; /* Durbin-Koopman variance */
        K->V = K->F = K->A = K->K = K->P = NULL;
	K->VS0 = K->VY0 = NULL;
        K->L = K->J = NULL;
        K->y = K->x = NULL;
        K->mu = NULL;
        K->U = NULL;
        K->Vsd = NULL;
        K->matcall = NULL;
        K->varying = NULL;
        K->step = NULL;
        K->uinfo = NULL;
        K->djinfo = NULL;
        K->flags = flags;
        K->t = 0;
        K->prn = NULL;
        K->data = NULL;
        K->b = NULL;
        K->d = 0;
        K->j = 0;
	K->smo_prep = 0;
	K->dj_initted = 0;
    }

    return K;
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

    /* K->q is non-zero only when using the Durbin-Koopman
       representation of the variance of the state disturbance as
       R*Q*R'
    */
    K->q = (K->Q != NULL)? gretl_matrix_rows(K->Q): 0;
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

static int check_matrix_dims (kalman *K, const gretl_matrix *m, int i)
{
    int r = 0, c = 0, symm = (i == K_VS || i == K_VY);
    int err = 0;

    if (i == K_T || i == K_VS || i == K_P) {
        r = c = K->r;
    } else if (i == K_BT) {
        r = K->k;
        c = K->n;
    } else if (i == K_ZT) {
        r = K->r;
        c = K->n;
    } else if (i == K_VY)  {
        r = c = K->n;
    } else if (i == K_a || i == K_m) {
        r = K->r;
        c = 1;
    }

    if (m->rows != r || m->cols != c) {
        gretl_errmsg_sprintf(_("kalman: %s is %d x %d, should be %d x %d\n"),
                             kalman_matrix_name(i), m->rows, m->cols, r, c);
        err = E_NONCONF;
    } else if (symm && !gretl_matrix_is_symmetric(m)) {
        gretl_errmsg_sprintf(_("kalman: %s is not symmetric\n"),
                kalman_matrix_name(i));
        err = E_NONCONF;
    }

    return err;
}

static int maybe_resize_export_matrix (kalman *K, gretl_matrix *m, int i)
{
    int rows = K->N, cols = 0;
    int err = 0;

    if (i == K_V) {
        cols = K->n;
    } else if (i == K_F) {
        /* vech(F) per row */
        cols = (K->n * K->n + K->n) / 2;
    } else if (i == K_A) {
        cols = K->r;
    } else if (i == K_BIG_P) {
	/* vech(P) per row */
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
        /* de Jong variance representation */
        if (K->H->rows != K->r || K->H->cols != K->p) {
            fprintf(stderr, "H is %d x %d, p=%d, r=%d, n=%d\n",
                    K->H->rows, K->H->cols, K->p, K->r, K->n);
            err = E_NONCONF;
        }
        if (!err && K->G != NULL) {
            if (K->G->rows != K->n || K->G->cols != K->p) {
                fprintf(stderr, "G is %d x %d, p=%d, r=%d, n=%d\n",
                        K->G->rows, K->G->cols, K->p, K->r, K->n);
                err = E_NONCONF;
            }
        }
    } else if (K->R != NULL) {
        /* state disturbances as per Durbin-Koopman */
        if (K->Q == NULL) {
            err = missing_matrix_error("statevar");
        } else if (K->R->rows != K->r || K->R->cols != K->Q->rows) {
            fprintf(stderr, "R is %d x %d, Q is %d x %d, r=%d\n",
                    K->R->rows, K->R->cols, K->Q->rows, K->Q->cols,
                    K->r);
            err = E_NONCONF;
        }
    } else {
        /* VS is mandatory, should be r x r and symmetric */
        if (!err) {
            err = check_matrix_dims(K, K->VS, K_VS);
        }
        /* VY should be n x n and symmetric, if present */
        if (!err && K->VY != NULL) {
            err = check_matrix_dims(K, K->VY, K_VY);
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

    /* x should have N rows to match y; and it should have either k or k - 1
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

    /* big V should be N x n */
    if (K->V != NULL) {
        err = maybe_resize_export_matrix(K, K->V, K_V);
    }

    /* big F should be N x n*n */
    if (!err && K->F != NULL) {
        err = maybe_resize_export_matrix(K, K->F, K_F);
    }

    /* big A should be N x r */
    if (!err && K->A != NULL) {
        err = maybe_resize_export_matrix(K, K->A, K_A);
    }

    /* big P should be N x nr, or N x r*r */
    if (!err && K->P != NULL) {
        err = maybe_resize_export_matrix(K, K->P, K_BIG_P);
    }

    /* LL should be N x 1 */
    if (!err && K->LL != NULL) {
        err = maybe_resize_export_matrix(K, K->LL, K_LL);
    }

    /* K (gain) should be N x (r * n) */
    if (!err && K->K != NULL) {
        err = maybe_resize_export_matrix(K, K->K, K_K);
    }

 bailout:

    if (err) {
        fprintf(stderr, "kalman_check_dimensions: err = %d\n", err);
    }

    return err;
}

static int kalman_init (kalman *K)
{
    int err = 0;

    if (trace) {
        printf("kalman_init(): aini %p, Pini %p\n",
               (void *) K->aini, (void *) K->Pini);
    }

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
    K->vt = gretl_matrix_alloc(K->n, 1);

    err = get_gretl_matrix_err();
    if (err) {
        return err;
    }

    K->Blk = gretl_matrix_block_new(&K->PZ,  K->r, K->n, /* P*Z */
                                    &K->Ft,  K->n, K->n, /* (Z'*P*Z + R)^{-1} */
                                    &K->iFt, K->n, K->n, /* Ft-inverse */
                                    &K->Kt,  K->r, K->n, /* gain at t */
                                    &K->Mt,  K->r, K->n, /* intermediate term */
                                    &K->Ct,  K->r, K->r, /* intermediate term */
                                    &K->Lt,  K->r, K->r, /* ditto */
                                    &K->rr,  K->r, K->r, /* ditto */
                                    NULL);

    if (K->Blk == NULL) {
        err = E_ALLOC;
    }

    if (!err && !kalman_is_bundle(K)) {
        /* in the "user-bundle" case we do this later */
        err = set_initial_statevar(K);
    }

    return err;
}

/* The following section includes functions that support the "plain C"
   Kalman API, as opposed to the bundle-based interface that subserves
   userland state-space functionality. We use this API in gretl's arma
   plugin in the special case of ARIMA with non-zero order of
   integration along with missing values. It's also possible that
   third-party users of libgretl might wish to use it.
*/

/**
 * kalman_new:
 * @a: r x 1 initial state vector.
 * @P: r x r initial precision matrix.
 * @T: r x r state transition matrix.
 * @BT: n x k matrix of coefficients on exogenous variables in the
 * observation equation, transposed.
 * @ZT: n x r matrix of coefficients on the state variables in the
 * observation equation, transposed.
 * @VS: r x r contemporaneous covariance matrix for the errors in the
 * state equation.
 * @VY: n x n contemporaneous covariance matrix for the errors in the
 * observation equation (or NULL if this is not applicable).
 * @y: T x n matrix of observable variable(s).
 * @x: T x k matrix of exogenous variable(s).  May be NULL if there
 * are no exogenous variables, or if there's only a constant.
 * @mu: r x 1 vector of constants in the state transition, or NULL.
 * @V: T x n matrix in which to record forecast errors (or NULL if
 * this is not required).
 * @err: location to receive error code.
 *
 * Allocates and initializes a Kalman struct, which can subsequently
 * be used for forecasting with kalman_forecast().
 *
 * Returns: pointer to allocated struct, or NULL on failure, in
 * which case @err will receive a non-zero code.
 */

kalman *kalman_new (gretl_matrix *a, gretl_matrix *P,
                    gretl_matrix *T, gretl_matrix *BT,
                    gretl_matrix *ZT, gretl_matrix *VS,
                    gretl_matrix *VY, gretl_matrix *y,
                    gretl_matrix *x, gretl_matrix *mu,
                    gretl_matrix *V, int *err)
{
    kalman *K;

    *err = 0;

    if (y == NULL || T == NULL || ZT == NULL || VS == NULL) {
        fprintf(stderr, "kalman_new: y=%p, T=%p, ZT=%p, VS=%p\n",
                (void *) y, (void *) T, (void *) ZT, (void *) VS);
        *err = missing_matrix_error(NULL);
        return NULL;
    }

    K = kalman_new_empty(0);
    if (K == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    /* use pointers for input matrices, don't copy */
    K->T = T;
    K->BT = BT;
    K->ZT = ZT;
    K->VS = VS;
    K->VY = VY;
    K->y = y;
    K->x = x;
    K->aini = a;
    K->Pini = P;
    K->mu = mu;

    /* output, but again use external pointer */
    K->V = V;

    kalman_set_dimensions(K);

    *err = kalman_check_dimensions(K);
    if (*err) {
        free(K);
        return NULL;
    }

    *err = kalman_init(K);

    if (*err) {
        kalman_free(K);
        K = NULL;
    } else {
        gretl_matrix_zero(K->vt);
    }

    return K;
}

/**
 * kalman_get_loglik:
 * @K: pointer to Kalman struct.
 *
 * Retrieves the log-likelhood calculated via a run of
 * kalman_forecast().
 *
 * Returns: ll value, or #NADBL on failure.
 */

double kalman_get_loglik (const kalman *K)
{
    return K->loglik;
}

/**
 * kalman_get_arma_variance:
 * @K: pointer to Kalman struct.
 *
 * Retrieves the estimated variance for an ARMA model
 * estimated using the Kalman filter.
 *
 * Returns: sigma-squared value, or #NADBL on failure.
 */

double kalman_get_arma_variance (const kalman *K)
{
    if (na(K->SSRw)) {
        return NADBL;
    } else {
        return K->SSRw / K->okN;
    }
}

/**
 * kalman_set_initial_state_vector:
 * @K: pointer to Kalman struct.
 * @a: vector of values to set.
 *
 * Resets the initial value of the state vector in a Kalman
 * struct, using the values from @a.  See also kalman_new().
 *
 * Returns: 0 on success, non-zero on error.
 */

int kalman_set_initial_state_vector (kalman *K, const gretl_vector *a)
{
    return gretl_matrix_copy_values(K->a0, a);
}

/**
 * kalman_set_initial_MSE_matrix:
 * @K: pointer to Kalman struct.
 * @P: matrix of values to set.
 *
 * Resets the initial value of the MSE matrix in a Kalman
 * struct, using the values from @P.  See also kalman_new().
 *
 * Returns: 0 on success, non-zero on error.
 */

int kalman_set_initial_MSE_matrix (kalman *K, const gretl_matrix *P)
{
    return gretl_matrix_copy_values(K->P0, P);
}

void kalman_attach_data (kalman *K, void *data)
{
    if (K != NULL) {
        K->data = data;
    }
}

void *kalman_get_data (const kalman *K)
{
    return (K != NULL)? K->data : NULL;
}

void kalman_attach_printer (kalman *K, PRN *prn)
{
    if (K != NULL) {
        K->prn = prn;
    }
}

PRN *kalman_get_printer (const kalman *K)
{
    return (K != NULL)? K->prn : NULL;
}

void kalman_set_arma_ll (kalman *K)
{
    K->flags |= KALMAN_ARMA_LL;
}

/* end of functions dedicated to non-bundle C API */

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

#define kappa 1.0e7 /* 1.0e6? */

static int diffuse_Pini (kalman *K)
{
    if (trace) {
        printf("diffuse_Pini(), exact = %d\n", K->exact);
    }

    if (K->exact) {
        if (kalman_dejong(K)) {
            fast_copy_values(K->P0, K->VS); /* H*H' */
        } else {
            /* univariate case */
#if P0_ZERO
            gretl_matrix_zero(K->P0);
#else
            fast_copy_values(K->P0, K->VS);
#endif
            /* initialize P∞ */
            if (K->Pk0 == NULL) {
                K->Pk0 = gretl_identity_matrix_new(K->r);
                if (K->Pk0 == NULL) {
                    return E_ALLOC;
                }
            } else {
                fast_write_I(K->Pk0);
            }
        }
    } else {
        /* old-style, using big kappa */
        int i;

        gretl_matrix_zero(K->P0);
        for (i=0; i<K->r; i++) {
            gretl_matrix_set(K->P0, i, i, kappa);
        }
    }

    return 0;
}

static int hamilton_Pini (kalman *K)
{
    gretl_matrix *Svar;
    gretl_matrix *vQ;
    int r2 = K->r * K->r;
    int err = 0;

    if (trace) {
        printf("hamilton_Pini()\n");
    }

    Svar = gretl_matrix_alloc(r2, r2);
    vQ = gretl_column_vector_alloc(r2);

    if (Svar == NULL || vQ == NULL) {
        gretl_matrix_free(Svar);
        gretl_matrix_free(vQ);
        return E_ALLOC;
    }

    gretl_matrix_kronecker_product(K->T, K->T, Svar);
    gretl_matrix_I_minus(Svar);
    gretl_matrix_vectorize(vQ, K->VS);

    err = gretl_LU_solve(Svar, vQ);
    if (err) {
        /* failed: are some of the eigenvalues out of bounds? */
        err = statemat_out_of_bounds(K);
        if (err == E_SINGULAR) {
            err = diffuse_Pini(K);
            K->flags |= KALMAN_DIFFUSE;
        }
	if (kdebug) {
	    printf("Hamilton-style P0 failed, set KALMAN_DIFFUSE\n");
	}
    } else {
	if (kdebug) {
	    printf("Hamilton-style P0 succeeded\n");
	}
        gretl_matrix_unvectorize(K->P0, vQ);
    }

    gretl_matrix_free(Svar);
    gretl_matrix_free(vQ);

    return err;
}

/* If the user has not given an initial value for P_{1|0}, compute
   this automatically as per Hamilton, ch 13, p. 378.  This works only
   if the eigenvalues of K->T lie inside the unit circle.  Failing
   that, or if the --diffuse option is given for the user Kalman
   filter, we apply a diffuse initialization.
*/

static int set_initial_statevar (kalman *K)
{
    if (K->Pini != NULL) {
        fast_copy_values(K->P0, K->Pini);
        return 0;
    } else if (K->flags & KALMAN_DIFFUSE) {
        return diffuse_Pini(K);
    } else {
        return hamilton_Pini(K);
    }
}

/* Write the vech of @src into row @t of @targ */

static void record_to_vech (gretl_matrix *targ,
			    const gretl_matrix *src,
			    int n, int t)
{
    int i, j, k = 0;
    double x;

    for (j=0; j<n; j++) {
        for (i=j; i<n; i++) {
            x = gretl_matrix_get(src, i, j);
            gretl_matrix_set(targ, t, k++, x);
        }
    }
}

/* Write the vec of @src into row @t of @targ */

static void record_to_vec (gretl_matrix *targ,
                           const gretl_matrix *src,
                           int t)
{
    int i;

    for (i=0; i<targ->cols; i++) {
        gretl_matrix_set(targ, t, i, src->val[i]);
    }
}

/* Write the vec of @src into column @t of @targ */

static void record_to_col (gretl_matrix *targ,
                           const gretl_matrix *src,
                           int t)
{
    size_t sz = targ->rows * sizeof(double);

    memcpy(targ->val + t * targ->rows, src->val, sz);
}

/* copy from vector @src into row @t of @targ */

static void record_to_row (gretl_matrix *targ,
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

static void record_to_row_offset (gretl_matrix *targ,
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

/* supports hansl function for creating a named Kalman bundle */

kalman *kalman_new_minimal (gretl_matrix *M[], int copy[],
                            int nmat, int dkvar, int *err)
{
    gretl_matrix **targ[5];
    kalman *K;
    int i;

    if (trace) {
        printf("kalman_new_minimal (dkvar = %d)\n", dkvar);
    }

    *err = 0;

    if (M[0] == NULL || M[1] == NULL || M[2] == NULL || M[3] == NULL) {
        fprintf(stderr, "kalman_new_minimal: nmat=%d, y=%p, Z=%p, T=%p, Q=%p\n",
                nmat, (void *) M[0], (void *) M[1], (void *) M[2], (void *) M[3]);
        *err = missing_matrix_error(NULL);
        return NULL;
    }

    K = kalman_new_empty(KALMAN_BUNDLE);
    if (K == NULL) {
        *err = E_ALLOC;
        return NULL;
    }

    kalman_check_env(K);

    targ[0] = &K->y;
    targ[1] = &K->ZT;
    targ[2] = &K->T;

    if (nmat == 5) {
        if (dkvar) {
            K->vartype = DK_VAR;
            targ[3] = &K->Q;
            targ[4] = &K->R;
        } else {
            K->vartype = DJ_VAR;
            targ[3] = &K->H;
            targ[4] = &K->G;
        }
    } else {
        targ[3] = &K->VS; /* statevar */
    }

    for (i=0; i<nmat; i++) {
        if (copy[i]) {
            *targ[i] = gretl_matrix_copy(M[i]);
        } else {
            *targ[i] = M[i];
        }
    }

    kalman_set_dimensions(K);

    if (K->vartype != STD_VAR) {
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
        gretl_matrix_zero(K->vt);
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

/* Ensure we have suitable matrices into which to write
   VS, VY and HG'
*/

static int ensure_cross_covariance_matrices (kalman *K)
{
    int r[] = {K->H->rows, K->G->rows, K->H->rows};
    int c[] = {K->H->rows, K->G->rows, K->G->rows};
    gretl_matrix **V[] = {&K->VS, &K->VY, &K->HG};
    gretl_matrix **targ;
    int i, err = 0;

    for (i=0; i<3 && !err; i++) {
        targ = V[i];
        if (*targ == NULL) {
            *targ = gretl_matrix_alloc(r[i], c[i]);
        } else {
            gretl_matrix *Vi = *targ;

            if (Vi->rows != r[i] || Vi->cols != c[i]) {
                gretl_matrix_realloc(Vi, r[i], c[i]);
            }
        }
        if (*targ == NULL) {
            err = E_ALLOC;
        }
    }

    return err;
}

static int ensure_DK_covariance_matrices (kalman *K)
{
    int r[] = {K->r, K->q};
    int c[] = {K->r, K->r};
    gretl_matrix **V[] = {&K->VS, &K->QRT};
    gretl_matrix **targ;
    int i, err = 0;

    for (i=0; i<2 && !err; i++) {
        targ = V[i];
        if (*targ == NULL) {
            *targ = gretl_matrix_alloc(r[i], c[i]);
        } else {
            gretl_matrix *Vi = *targ;

            if (Vi->rows != r[i] || Vi->cols != c[i]) {
                gretl_matrix_realloc(Vi, r[i], c[i]);
            }
        }
        if (*targ == NULL) {
            err = E_ALLOC;
        }
    }

    return err;
}

/* After reading H and G from the user, either at (re-)initialization
   or at a given time-step in the case where either of these matrices
   is time-varying, record the user input in K->H and K->G and form
   the 'real' VS and VY.  But note that in the time-step case it may be
   that only one of VS, VY needs to be treated in this way (if only one
   is time-varying, only one will have been redefined via a function
   call).
*/

static int kalman_update_crossinfo (kalman *K, int mode)
{
    int err = 0;

    /* Note that H and G may be needed as such for simulation */

    if (mode == UPDATE_INIT) {
        err = ensure_cross_covariance_matrices(K);
        if (err) {
            return err;
        }
    }

    if (mode == UPDATE_INIT || matrix_is_varying(K, K_VS)) {
        /* (re)create VS using modified H */
        err = gretl_matrix_multiply_mod(K->H, GRETL_MOD_NONE,
                                        K->H, GRETL_MOD_TRANSPOSE,
                                        K->VS, GRETL_MOD_NONE);
    }

    if (!err && (mode == UPDATE_INIT || matrix_is_varying(K, K_VY))) {
        /* (re)create VY using modified G */
        err = gretl_matrix_multiply_mod(K->G, GRETL_MOD_NONE,
                                        K->G, GRETL_MOD_TRANSPOSE,
                                        K->VY, GRETL_MOD_NONE);
    }

    if (!err && (mode == UPDATE_INIT || matrix_is_varying(K, K_VS) ||
                 matrix_is_varying(K, K_VY))) {
        /* (re)create HG' using modified H and/or G */
        err = gretl_matrix_multiply_mod(K->H, GRETL_MOD_NONE,
                                        K->G, GRETL_MOD_TRANSPOSE,
                                        K->HG, GRETL_MOD_NONE);
    }

    return err;
}

static int kalman_update_dkvar (kalman *K, int mode)
{
    int err = 0;

    if (mode == UPDATE_INIT) {
        err = ensure_DK_covariance_matrices(K);
        if (err) {
            return err;
        }
    }

    if (mode == UPDATE_INIT || matrix_is_varying(K, K_VS)) {
        /* (re)create VS and QR' using modified Q */
        //fprintf(stderr, "*** kalman_update_dkvar, mode %d\n", mode);
        gretl_matrix_qform(K->R, GRETL_MOD_NONE, K->Q,
                           K->VS, GRETL_MOD_NONE);
        err = gretl_matrix_multiply_mod(K->Q, GRETL_MOD_NONE,
                                        K->R, GRETL_MOD_TRANSPOSE,
                                        K->QRT, GRETL_MOD_NONE);
     }

    return err;
}

/* kalman_revise_DJ_variance: the user has actually given H in place
   of VS and G in place of VY; so we have to form VS, VY, and HG'
   (cross-correlated disturbances).

   This function is called in the course of initial set-up of a
   filter, and also when the Kalman matrices are being re-checked at
   the start of filtering, smoothing or simulation.
*/

static int kalman_revise_DJ_variance (kalman *K)
{
    int err = 0;

    if (K->H == NULL || K->G == NULL) {
        fprintf(stderr, "K->H %p, K->G %p\n", (void *) K->H, (void *) K->G);
        return missing_matrix_error("'statevar' or 'obsvar'");
    }

    err = kalman_update_crossinfo(K, UPDATE_INIT);

    if (err) {
        fprintf(stderr, "kalman_revise_cross_variance: err = %d\n", err);
    }

    return err;
}

static int kalman_revise_DK_variance (kalman *K)
{
    int err = 0;

    if (K->Q == NULL || K->R == NULL) {
        return missing_matrix_error("'statevar' or 'statsel'");
    }

    err = kalman_update_dkvar(K, UPDATE_INIT);

    if (err) {
        fprintf(stderr, "kalman_revise_DK_variance: err = %d\n", err);
    }

    return err;
}

static int kalman_revise_variance (kalman *K)
{
    if (K->vartype == DJ_VAR) {
        return kalman_revise_DJ_variance(K);
    } else {
        return kalman_revise_DK_variance(K);
    }
}

static void kalman_print_state (kalman *K)
{
   int j;

    if (K->t > 4) return;

    fprintf(stderr, "t = %d:\n", K->t);

    for (j=0; j<K->n; j++) {
        fprintf(stderr, "y[%d] = %.10g, err[%d] = %.10g\n", j,
                gretl_matrix_get(K->y, K->t, j),
                j, gretl_vector_get(K->vt, j));
    }

    gretl_matrix_print(K->Kt, "K->Kt");
    gretl_matrix_print(K->a0, "K->a0");
    gretl_matrix_print(K->P0, "K->P0");
}

/* On filtering: record the state and/or its variance, as
   wanted.
*/

static int kalman_record_state (kalman *K)
{
    int err = 0;

    if (K->A != NULL) {
        record_to_row(K->A, K->a0, K->t);
    }
    if (K->P != NULL) {
        record_to_vech(K->P, K->P0, K->r, K->t);
    }
    if (K->exact && K->djinfo != NULL && K->djinfo->A != NULL) {
	/* record "augmented" quantities */
        record_to_col(K->djinfo->A, K->djinfo->At, K->t);
        record_to_col(K->djinfo->V, K->djinfo->Vt, K->t);
    }

    return err;
}

static double get_xjt (kalman *K, int j)
{
    if (K->ifc) {
        return (j == 0)? 1.0 : gretl_matrix_get(K->x, K->t, j-1);
    } else {
        return gretl_matrix_get(K->x, K->t, j);
    }
}

/* Read from the appropriate row of x (N x k) and multiply by B' to
   form B'x_t.  Note: the flag K->ifc is used to indicate that the
   observation equation has an implicit constant, with an entry in
   the B matrix (the first) but no explicit entry in the x matrix.

   The case where x is NULL and B an n-vector (implicit constant)
   is also handled.

   There's no need to store B'x_t: we either want to add this to,
   or subtract it from, the n-vector @targ (the operator being
   indicated by the @mod argument).

   Returns: the number of missing values encountered.
*/

static int kalman_do_Bx (kalman *K, gretl_matrix *targ,
                         GretlMatrixMod mod)
{
    double xjt, bji, bxi;
    int i, j, missvals = 0;

    for (i=0; i<K->n; i++) {
        if (K->x == NULL) {
            /* the implicit constant case */
            bxi = K->BT->val[i];
        } else {
            bxi = 0;
            for (j=0; j<K->k; j++) {
                xjt = get_xjt(K, j);
                bji = gretl_matrix_get(K->BT, j, i);
                if (bji != 0) {
                    /* here we'll implicitly take 0 * NA as 0 */
                    if (na(xjt)) {
                        missvals++;
                    } else {
                        bxi += xjt * bji;
                    }
                }
            }
        }
        if (mod == GRETL_MOD_DECREMENT) {
            targ->val[i] -= bxi;
        } else {
            targ->val[i] += bxi;
        }
    }

    return missvals;
}

/* version of the above for use with the univariate approach */

static double kalman_Bx_uni (kalman *K, int i)
{
    double xjt, bji, bxi;
    int j;

    if (K->x == NULL) {
        bxi = K->BT->val[i];
    } else {
        bxi = 0;
        for (j=0; j<K->k; j++) {
            xjt = get_xjt(K, j);
            bji = gretl_matrix_get(K->BT, j, i);
            if (bji != 0 && na(xjt)) {
                return NADBL;
            } else {
                bxi += xjt * bji;
            }
        }
    }

    return bxi;
}

/* Compute the one-step ahead forecast error:

   v_t = y_t - B_t*x_t - Z_t*a_t

   Returns the number of valid observables at time step t,
   or 0 if the forecast error cannot be computed due to
   missing values.
*/

static int compute_forecast_error (kalman *K)
{
    int i;

    /* initialize v_t to y_t */
    for (i=0; i<K->n; i++) {
        K->vt->val[i] = gretl_matrix_get(K->y, K->t, i);
    }

    if (K->BT != NULL) {
        /* subtract effect of exogenous terms, if any */
        kalman_do_Bx(K, K->vt, GRETL_MOD_DECREMENT);
    }

    /* check for any missing values in v_t */
    for (i=0; i<K->n; i++) {
        if (na(K->vt->val[i])) {
            return 0;
        }
    }

    /* subtract contribution from state */
    gretl_matrix_multiply_mod(K->ZT,  GRETL_MOD_TRANSPOSE,
                              K->a0, GRETL_MOD_NONE,
                              K->vt, GRETL_MOD_DECREMENT);

    return K->n;
}

/* Given a unified function to update one or more of the potentially
   time-varying matrices, try to figure out which matrix or matrices
   are actually modified by this function.
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
            gretl_errmsg_sprintf(_("Couldn't find function '%s'"), K->matcall);
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
			if (kdebug) {
			    fprintf(stderr, "matrix %s is varying\n", s);
			}
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

/* Function to update any time-varying matrices, for use with a kalman
   bundle. Bypasses the regular "genr" apparatus, passing the attached
   bundle directly to the given user function after it has been found
   by name.
*/

static int kalman_update_matrices (kalman *K, PRN *prn)
{
    ufunc *uf;
    fncall *fc;
    int err = 0;

    uf = get_user_function_by_name(K->matcall);
    if (uf == NULL) {
        gretl_errmsg_sprintf(_("Couldn't find function '%s'"), K->matcall);
        return E_DATA;
    }

    if (K->varying == NULL) {
        check_for_matrix_updates(K, uf);
    }

    fc = fncall_new(uf, 0);
    err = push_anon_function_arg(fc, GRETL_TYPE_BUNDLE_REF, K->b);

    if (!err) {
        err = gretl_function_exec(fc, GRETL_TYPE_NONE, NULL, NULL, prn);
    }
    if (err) {
        fprintf(stderr, "kalman_update_matrices: call='%s', err=%d\n",
                K->matcall, err);
    }

    return err;
}

/* On getting a new K->ZT in case of time-variation, redo any
   transformations that may be needed.
*/

static int adjust_univariate_Z (kalman *K)
{
    int err = 0;

    if (K->r > 1 || K->n > 1) {
        if (kdebug > 1) {
            fprintf(stderr, "doing adjust_univariate_Z\n");
        }
        /* transposition needed */
        if (K->uinfo->Z == NULL) {
            fprintf(stderr, "K->uinfo->Z not allocated!\n");
            err = E_DATA;
        } else {
            err = gretl_matrix_transpose(K->uinfo->Z, K->ZT);
        }
        if (!err && K->uinfo->Linv != NULL) {
            /* diagonalization needed (FIXME conditionality) */
            gretl_matrix_multiply_mod(K->uinfo->Linv, GRETL_MOD_NONE,
                                      K->ZT, GRETL_MOD_TRANSPOSE,
                                      K->uinfo->Z, GRETL_MOD_NONE);
        }
    }

    return err;
}

/* If we have any time-varying coefficient matrices, refresh these for
   the current time step. This is called on a forward filtering pass.
*/

static int kalman_refresh_matrices (kalman *K, PRN *prn)
{
    gretl_matrix **mptr[] = {
        &K->T, &K->BT, &K->ZT, &K->VS, &K->VY, &K->mu
    };
    int cross_update = 0;
    int dkvar_update = 0;
    int i, err = 0;

    if (kalman_djvar(K)) {
        mptr[3] = &K->H;
        mptr[4] = &K->G;
    } else if (kalman_dkvar(K)) {
        mptr[3] = &K->Q;
        mptr[4] = &K->R;
    }

    if (K->matcall != NULL) {
        err = kalman_update_matrices(K, prn);
    }

    for (i=0; i<K_N_MATCALLS && !err; i++) {
        if (matrix_is_varying(K, i)) {
            if (kalman_djvar(K) && (i == K_VS || i == K_VY)) {
                /* handle revised H and/or G */
                cross_update = 1;
            } else if (kalman_dkvar(K) && i == K_VS) {
                /* handle revised Q and/or R */
                dkvar_update = 1;
            } else {
                err = check_matrix_dims(K, *mptr[i], i);
            }
            if (err) {
                fprintf(stderr, "kalman_refresh_matrices: err = %d at t = %d\n",
                        err, K->t);
            } else if (i == K_ZT && kalman_univariate(K)) {
                err = adjust_univariate_Z(K);
            } else if (i == K_T) {
                K->TI = gretl_is_identity_matrix(K->T);
            }
        }
    }

    if (!err && K->step != NULL) {
        /* keep a record of T and/or Z' at the given time step */
        if (K->step->T != NULL) {
            record_to_col(K->step->T, K->T, K->t);
        }
        if (K->step->ZT != NULL) {
            if (K->uinfo != NULL) {
                record_to_col(K->step->ZT, K->uinfo->Z, K->t);
            } else {
                record_to_col(K->step->ZT, K->ZT, K->t);
            }
        }
    }

    if (!err) {
        if (cross_update) {
            /* cross-correlated case */
            err = kalman_update_crossinfo(K, UPDATE_STEP);
        } else if (dkvar_update) {
            /* Durbin-Koopman case */
            err = kalman_update_dkvar(K, UPDATE_STEP);
        }
    }

    return err;
}

/* Handling of missing y_t or x_t, for the case of univariate y_t
   (or all y_t values missing).
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

    /* var(state) update: P1 = T*P0*T' + VS (C = 0) */
    fast_copy_values(K->P1, K->VS);
    gretl_matrix_qform(K->T, GRETL_MOD_NONE, K->P0,
                       K->P1, GRETL_MOD_CUMULATE);
    fast_copy_values(K->P0, K->P1);

    /* record stuff if wanted */
    if (K->F != NULL) {
        if (K->smo_prep) {
            set_row_to_value(K->F, K->t, 0.0);
        } else {
            set_row_to_value(K->F, K->t, NADBL); /* 0 ? */
        }
    }
    if (K->LL != NULL) {
        gretl_vector_set(K->LL, K->t, 0.0); /* ? */
    }
    if (K->K != NULL) {
        set_row_to_value(K->K, K->t, 0.0);
    }
    if (K->V != NULL) {
        if (K->smo_prep) {
            set_row_to_value(K->V, K->t, 0.0);
        } else {
            set_row_to_value(K->V, K->t, NADBL); /* 0 ? */
        }
    }
    if (K->L != NULL) {
	record_to_col(K->L, K->T, K->t);
    }
    if (K->J != NULL) {
	record_to_col(K->J, K->H, K->t);
    }
    /* more to be done if we're in the de Jong diffuse phase? */
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
    double sumldet = 0;
    int nt, err = 0;

    if (trace) {
        printf("kalman_forecast (legacy)\n");
    }
    if (kdebug > 1) {
        fprintf(stderr, "\n*** kalman_forecast: N=%d, n=%d ***\n",
                K->N, K->n);
    }

    K->SSRw = K->loglik = 0.0;
    K->s2 = NADBL;
    K->okN = K->N;
    set_kalman_running(K);

    for (K->t = 0; K->t < K->N && !err; K->t += 1) {
        double llt = NADBL;
        double ldet = 0, qt = 0;

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

        if (kdebug > 1) {
            kalman_print_state(K);
        }

        /* calculate v_t, checking for missing values */
        nt = compute_forecast_error(K);
        if (nt == 0) {
            /* skip this observation */
            K->okN -= 1;
            handle_missing_obs(K);
            continue;
        }

        /* calculate F_t = ZPZ' [+ VY] */
        if (K->VY != NULL) {
            fast_copy_values(K->Ft, K->VY);
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

        /* standard Kalman procedure */
        fast_copy_values(K->iFt, K->Ft);
        err = gretl_invert_symmetric_matrix2(K->iFt, &ldet);
        if (err) {
            fprintf(stderr, "kalman_forecast: failed to invert Ft\n");
        } else {
            qt = gretl_scalar_qform(K->vt, K->iFt, &err);
        }

        if (K->F != NULL) {
            /* we're recording F_t for all t */
            if (K->smo_prep) {
                /* record inverse */
                record_to_vech(K->F, K->iFt, nt, K->t);
            } else {
                /* record F_t itself */
                record_to_vech(K->F, K->Ft, nt, K->t);
            }
        }

        /* determine and record loglikelihood */
        if (err) {
            K->loglik = NADBL;
            break;
        } else {
            llt = -0.5 * (ll0 + ldet + qt);
            if (na(llt)) {
                K->loglik = NADBL;
                break;
            }
            K->loglik += llt;
            K->SSRw += qt;
            sumldet += ldet;
        }
        if (K->LL != NULL) {
            gretl_vector_set(K->LL, K->t, llt);
        }

        /* Calculate gain K_t = M_t F_t^{-1}, and C matrix */
        gretl_matrix_multiply(K->Mt, K->iFt, K->Kt);
        gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
                                  K->Mt, GRETL_MOD_TRANSPOSE,
                                  K->Ct, GRETL_MOD_NONE);
        if (!err && K->K != NULL) {
            /* record the gain */
            record_to_vec(K->K, K->Kt, K->t);
        }

        /* update state: a1 = T a0 + K_t v_t [+ mu] */
        gretl_matrix_multiply(K->T, K->a0, K->a1);
        if (K->mu != NULL) {
            gretl_matrix_add_to(K->a1, K->mu);
        }
        gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
                                  K->vt, GRETL_MOD_NONE,
                                  K->a1, GRETL_MOD_CUMULATE);
        fast_copy_values(K->a0, K->a1);

        /* update var(state): P1 = TPT' + VS - C */
        fast_copy_values(K->P1, K->VS);
        gretl_matrix_qform(K->T, GRETL_MOD_NONE,
                           K->P0, K->P1, GRETL_MOD_CUMULATE);
        gretl_matrix_subtract_from(K->P1, K->Ct);
        fast_copy_values(K->P0, K->P1);

        /* record forecast errors if wanted */
        if (!err && K->V != NULL) {
            record_to_row(K->V, K->vt, K->t);
        }
    }

    set_kalman_stopped(K);

    if (na(K->loglik)) {
        err = E_NAN;
    } else if (kalman_arma_ll(K)) {
        double ll1 = 1.0 + LN_2_PI + log(K->SSRw / K->okN);

        K->loglik = -0.5 * (K->okN * ll1 + sumldet);
    } else {
	int d = kalman_diffuse(K) ? K->r : 0;

        K->s2 = K->SSRw / (K->n * K->okN - d);
    }

    if (kdebug) {
        fprintf(stderr, "kalman_forecast: err=%d, ll=%#.8g\n",
		err, K->loglik);
    }

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
    { K_BT, "obsxmat" },
    { K_VY, "obsvar" },
    { K_T,  "statemat" },
    { K_VS, "statevar" },
    { K_m,  "stconst" },
    { K_a,  "inistate" },
    { K_P,  "inivar" },
    { K_R,  "statesel" }
};

int extra_mats[] = {
    K_BT,
    K_VY,
    K_m,
    K_x,
    K_a,
    K_P,
    K_R
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

    if (!err && (K->ZT == NULL || K->T == NULL || K->VS == NULL)) {
        fprintf(stderr, "kalman_bundle_kalman_recheck_matrices: Z=%p, T=%p, Q=%p\n",
                K->ZT, K->T, K->VS);
        err = missing_matrix_error(NULL);
    }

    if (err) {
        return err;
    }

    /* maybe redundant? */
    kalman_set_dimensions(K);

    if (gretl_matrix_rows(K->T) != K->r ||
        gretl_matrix_rows(K->BT) != K->k) {
        err = E_NONCONF;
    } else if (!kalman_simulating(K)) {
        err = obsy_check(K);
    }

    if (!err && K->vartype != STD_VAR) {
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
        err = set_initial_statevar(K);
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

static int ensure_univariate_info (kalman *K)
{
    if (K->uinfo == NULL) {
        return kalman_add_univar_info(K);
    } else {
        return 0;
    }
}

static int ensure_dejong_info (kalman *K)
{
    if (K->djinfo == NULL) {
        return kalman_add_dejong_info(K);
    } else {
        return 0;
    }
}

int kalman_run (kalman *K, PRN *prn, int *errp)
{
    int err = kalman_ensure_output_matrices(K);

    if (trace) {
        printf("kalman_run()\n");
    }

    if (!err) {
        gretl_matrix_zero(K->vt);
        /* includes setting of K->P0 */
        err = kalman_bundle_recheck_matrices(K, prn);
    }

    if (!err && K->LL == NULL) {
        K->LL = gretl_matrix_alloc(K->N, 1);
        if (K->LL == NULL) {
            err = E_ALLOC;
        }
    }

    if (!err) {
        if (kalman_univariate(K)) {
            err = ensure_univariate_info(K);
        } else if (K->exact) {
            /* only if diffuse? */
            err = ensure_dejong_info(K);
        }
    }

    if (K->exact && kalman_univariate(K)) {
        /* in case we're re-running the filter */
        fast_write_I(K->Pk0);
    }

    if (!err) {
        if (kalman_univariate(K)) {
            err = kfilter_univariate(K, prn);
        } else if (kalman_dejong(K)) {
	    err = kfilter_dejong(K, prn);
	} else {
	    err = kalman_forecast(K, prn);
	}
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

/* Implements the userland kfilter() function */

int kalman_bundle_filter (gretl_bundle *b, PRN *prn, int *errp)
{
    kalman *K = gretl_bundle_get_private_data(b);

    if (trace) {
        printf("kalman_bundle_filter()\n");
    }

    K->b = b; /* attach bundle pointer */

    return kalman_run(K, prn, errp);
}

/* write row @t of matrix @src into matrix @targ */

static void load_from_row (gretl_matrix *targ,
                           const gretl_matrix *src,
                           int t)
{
    int j;

    for (j=0; j<src->cols; j++) {
        targ->val[j] = gretl_matrix_get(src, t, j);
    }
}

/* Copy row @t from @src into @targ; or add row @t of @src to @targ;
   or subtract row @t of @src from @targ.  We allow the possibility
   that the length of vector @targ is less than the number of columns
   in @src, but not the converse.
*/

static int vector_from_row_mod (gretl_vector *targ,
                                const gretl_matrix *src,
                                int t, GretlMatrixMod mod)
{
    int i, n = gretl_vector_get_length(targ);
    double x;

    if (n > src->cols) {
        fprintf(stderr, "vector_from_row: targ length = %d, but src "
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

/* Similar to vector_from_row_offset() except that a column offset,
   @j, is supported for the reading of a row from @src, and we don't
   support the @mod option.
*/

static int vector_from_row_offset (gretl_vector *targ,
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
   matrix in @targ.
*/

static void load_from_vech (gretl_matrix *targ, const gretl_matrix *src,
                            int n, int t)
{
    int i, j, k = 0;
    double x;

    for (j=0; j<n; j++) {
        for (i=j; i<n; i++) {
            x = gretl_matrix_get(src, t, k++);
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

/* Column @t of @src represents the vec of a certain matrix: extract the
   column and reconstitute the matrix in @targ.
*/

static void load_from_col (gretl_matrix *targ,
                           const gretl_matrix *src,
                           int t)
{
    size_t sz = src->rows * sizeof(double);

    memcpy(targ->val, src->val + t * src->rows, sz);
}

static int retrieve_Tt (kalman *K)
{
    if (K->step == NULL || K->step->T == NULL) {
        return E_DATA;
    } else {
        load_from_col((gretl_matrix *) K->T, K->step->T, K->t);
        return 0;
    }
}

static int retrieve_Zt (kalman *K)
{
    if (K->step == NULL || K->step->ZT == NULL) {
        return E_DATA;
    } else {
        if (K->uinfo != NULL) {
            load_from_col((gretl_matrix *) K->uinfo->Z, K->step->ZT, K->t);
        } else {
            load_from_col((gretl_matrix *) K->ZT, K->step->ZT, K->t);
        }
        return 0;
    }
}

static int load_filter_data (kalman *K, int no_collapse, int *err)
{
    int missing = 0;

    /* FIXME maybe we should have a byte array of length K->N
       to identify missing observations (constructed on the
       filtering step if needed). This could speed things
       up below,
    */

    /* load the forecast error into K->vt */
    load_from_row(K->vt, K->V, K->t);

    /* state: get T_t and/or Z_t if need be */
    if (matrix_is_varying(K, K_T)) {
        *err = retrieve_Tt(K);
    }
    if (!*err && matrix_is_varying(K, K_ZT)) {
        *err = retrieve_Zt(K);
    }
    if (*err) {
        return 0;
    }

    /* load the state and its MSE */
    load_from_row(K->a0, K->A, K->t);
    if (no_collapse) {
        /* FIXME? looks weird but seems to help */
        if (K->t == 0) {
            gretl_matrix_zero(K->P0);
        } else {
            load_from_vech(K->P0, K->P, K->r, K->t-1);
        }
    } else {
        load_from_vech(K->P0, K->P, K->r, K->t);
    }

    /* load the gain and F^{-1} */
    load_from_vec(K->Kt, K->K, K->t);
    load_from_vech(K->iFt, K->F, K->n, K->t);

    /* heuristic for missing observable (see note above) */
    missing = gretl_is_zero_matrix(K->Kt);

    return missing ? 0 : K->n;
}

/* kalman_ensure_stepinfo() is called in two cases:

  (a) If we're doing univariate filtering, K->ZT is time-varying, and
   K->n > 1, we need a record of K->ZT per time step for building a
   proper F matrix (pevar).

  (b) If we're doing smoothing for a system that has time-varying
   coefficients in K->T or K->ZT we record the vec of these matrices
   for each time-step on the forward pass.
*/

static int kalman_ensure_stepinfo (kalman *K, int smoothing)
{
    int alloc_T = 0;
    int alloc_Z = 0;
    int err = 0;

    if (K->step != NULL) {
        /* already started */
        if (smoothing && matrix_is_varying(K, K_T)) {
            if (K->step->T != NULL &&
                K->step->T->rows == K->r * K->r &&
                K->step->T->cols == K->N) {
                ; /* T already OK */
            } else {
                gretl_matrix_free(K->step->T);
                alloc_T = 1;
            }
        }
        if (matrix_is_varying(K, K_ZT)) {
            if (K->step->ZT != NULL &&
                K->step->ZT->rows == K->r * K->n &&
                K->step->ZT->cols == K->N) {
                ; /* ZT already OK */
            } else {
                gretl_matrix_free(K->step->ZT);
                alloc_Z = 1;
            }
        }
    } else {
        /* starting from scratch */
        K->step = malloc(sizeof *K->step);
        if (K->step == NULL) {
            return E_ALLOC;
        }
        K->step->T = K->step->ZT = NULL;
        alloc_T = smoothing && matrix_is_varying(K, K_T);
        alloc_Z = matrix_is_varying(K, K_ZT);
    }

    if (alloc_T) {
        K->step->T = gretl_matrix_alloc(K->r * K->r, K->N);
        if (K->step->T == NULL) {
            err = E_ALLOC;
        }
    }
    if (!err && alloc_Z) {
        K->step->ZT = gretl_matrix_alloc(K->r * K->n, K->N);
        if (K->step->ZT == NULL) {
            err = E_ALLOC;
        }
    }

    if (err) {
        free_stepinfo(K);
    }

    return err;
}

static int kalman_add_univar_info (kalman *K)
{
    int err = 0;

    K->uinfo = calloc(1, sizeof *K->uinfo);
    if (K->uinfo == NULL) {
        return E_ALLOC;
    }

    K->uinfo->Zi  = gretl_matrix_alloc(1, K->r);
    K->uinfo->Kti = gretl_matrix_alloc(K->r, 1);
    K->uinfo->Pki = gretl_matrix_alloc(K->r, K->r);
    K->uinfo->Kki = gretl_matrix_alloc(K->r, 1);
    K->uinfo->m1  = gretl_matrix_alloc(K->r, 1);
    K->uinfo->mm  = gretl_matrix_alloc(K->r, K->r);

    if (K->n > 1) {
        K->uinfo->UF = gretl_matrix_alloc(K->N, K->n);
    } else {
        K->uinfo->UF = NULL;
    }

    K->uinfo->free_Z = 0;
    K->uinfo->free_g = 0;
    K->uinfo->Linv = NULL;
    K->uinfo->y0 = NULL;

    if (K->r > 1 || K->n > 1) {
        /* we'll need untransposed Z */
        K->uinfo->Z = gretl_zero_matrix_new(K->n, K->r);
        K->uinfo->free_Z = 1;
    } else {
        K->uinfo->Z = K->ZT;
    }

    if (K->VY != NULL) {
        if (K->n > 1) {
            K->uinfo->g = gretl_zero_matrix_new(K->n, 1);
            K->uinfo->free_g = 1;
            /* we may need this too */
            K->uinfo->Linv = gretl_matrix_copy(K->VY);
        } else {
            K->uinfo->g = K->VY;
        }
    } else {
        K->uinfo->g = NULL;
    }

    K->uinfo->Finf = NULL;
    K->uinfo->Kinf = NULL;
    K->uinfo->Nd = 0;

    /* FIXME proper error checking above */

    return err;
}

static int kalman_add_dejong_info (kalman *K)
{
    int err = 0;

    K->djinfo = calloc(1, sizeof *K->djinfo);
    if (K->djinfo == NULL) {
        return E_ALLOC;
    }

    K->djinfo->Vt = gretl_zero_matrix_new(K->n, K->r);
    K->djinfo->At = gretl_zero_matrix_new(K->r, K->r);
    K->djinfo->Vv = gretl_matrix_alloc(K->n, K->r + 1);
    K->djinfo->Qt = gretl_zero_matrix_new(K->r + 1, K->r + 1);
    K->djinfo->S  = gretl_matrix_alloc(K->r, K->r);
    K->djinfo->s  = gretl_matrix_alloc(K->r, 1);
    K->djinfo->A  = NULL;
    K->djinfo->V  = NULL;
    K->djinfo->Am = NULL;

    return err;
}

/* optional matrix (old-style) to hold smoothed disturbances */

static int ensure_U_matrix (kalman *K)
{
    int Ucols = K->r;
    int Urows = K->N;
    int err = 0;

    if (K->VY != NULL) {
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

static int real_kalman_smooth (kalman *K, int dist, PRN *prn)
{
    int err = kalman_ensure_output_matrices(K);

    if (trace) {
        printf("real_kalman_smooth()\n");
    }

    if (!err && dist && K->code == K_LEGACY) {
        err = ensure_U_matrix(K);
    }

    if (err) {
        return err;
    }

    if (matrix_is_varying(K, K_T) || matrix_is_varying(K, K_ZT)) {
        /* add recorder for T_t and/or Z_t */
        err = kalman_ensure_stepinfo(K, 1);
        if (err) {
            goto bailout;
        }
    }

    if (K->exact && kalman_univariate(K)) {
        int Nd = matrix_is_varying(K, K_ZT) ? K->N : 4 * K->r;
        int rows = K->r * K->r;

        K->PK = gretl_zero_matrix_new(rows, Nd);
    }

    if (!err) {
        err = kalman_bundle_recheck_matrices(K, prn);
    }

    if (!err) {
        if (kalman_univariate(K)) {
            err = ensure_univariate_info(K);
        } else if (K->exact) {
            err = ensure_dejong_info(K);
        }
    }

    if (!err) {
        /* prior forward pass */
	K->smo_prep = dist ? SM_DIST : SM_STATE;
        if (kalman_univariate(K)) {
            err = kfilter_univariate(K, NULL);
        } else if (kalman_dejong(K)) {
            err = kfilter_dejong(K, NULL);
        } else {
	    err = kalman_forecast(K, NULL);
	}
        K->smo_prep = SM_NONE;
    }

    K->t = 0;

    if (!err) {
        if (kalman_univariate(K)) {
            err = ksmooth_univariate(K, dist);
	} else if (kalman_dejong(K)) {
            if (dist) {
                err = dist_smooth_dejong(K, dist > 1);
	    } else {
		err = state_smooth_dejong(K);
	    }
	} else {
            if (dist) {
		err = koopman_smooth(K, dist > 1);
	    } else {
		err = anderson_moore_smooth(K);
	    }
        }
    }

    /* free "special case" storage */
    if (K->PK != NULL) {
        gretl_matrix_free(K->PK);
        K->PK = NULL;
    }

 bailout:

    /* trash the "stepinfo" storage */
    free_stepinfo(K);

    return err;
}

/* For use with userland bundle-based API */

int kalman_bundle_smooth (gretl_bundle *b, int dist, PRN *prn)
{
    kalman *K = gretl_bundle_get_private_data(b);

    if (trace) {
        printf("kalman_bundle_smooth(), dist = %d\n", dist);
    }

    K->b = b; /* attach bundle pointer */

    return real_kalman_smooth(K, dist, prn);
}

/* For use with "plain C" API. TODO: support disturbance
   smoothing via @opt -- right now only state smoothing
   is offered.
*/

gretl_matrix *kalman_smooth (kalman *K, gretlopt opt,
                             PRN *prn, int *err)
{
    gretl_matrix *S = NULL;

    *err = real_kalman_smooth(K, 0, prn);

    if (K->A != NULL && K->P != NULL) {
        int r = K->A->rows;
        int c = K->A->cols;
        size_t Asize = r * c * sizeof(double);
        size_t Psize = 0;

        if (opt & OPT_M) {
            /* include MSE of state */
            c += K->P->cols;
            Psize = r * K->P->cols * sizeof(double);
        }
        if (r > 0 && c > 0) {
            S = gretl_matrix_alloc(r, c);
            if (S == NULL) {
                *err = E_ALLOC;
            }
        } else {
            *err = E_DATA;
        }
        if (S != NULL) {
            memcpy(S->val, K->A->val, Asize);
            if (opt & OPT_M) {
                memcpy(S->val + r * K->A->cols, K->P->val, Psize);
            }
        }
    }

    return S;
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
   Journal, 1999 (volume 2, pp. 113-166), section 4.2, regarding the
   initialization of the state under simulation.
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
        load_from_row(v0, U, 0);
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

/* note: it's OK for @S to be NULL (if the simulated state is not
   wanted), so watch out for that!
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
            record_to_row_offset(S, K->a0, 0, 0);
        }
        record_to_row_offset(Y, yt, 0, obs_offset);
        /* the first row of output is handled */
        tmin = 1;
    }

    if (K->p == 0 && K->VY != NULL) {
        /* we want to read observation disturbances */
        obsdist = 1;
    }

    for (K->t = tmin; K->t < K->N && !err; K->t += 1) {
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
            /* handle missing values? */
            kalman_do_Bx(K, yt, GRETL_MOD_CUMULATE);
        }
        if (K->p > 0) {
            /* G \varepsilon_t */
            load_from_row(et, U, K->t);
            gretl_matrix_multiply(K->G, et, K->vt);
        } else if (obsdist) {
            vector_from_row_offset(K->vt, U, K->t, K->r);
        }
        gretl_matrix_add_to(yt, K->vt);

        /* record the t-dated observables */
        record_to_row_offset(Y, yt, K->t, obs_offset);

        /* record the t-dated state? */
        if (S != NULL && tmin == 0) {
            record_to_row_offset(S, K->a0, K->t, 0);
        }

        /* a_{t+1} = T*a_t + v_t */
        gretl_matrix_multiply(K->T, K->a0, K->a1);
        if (K->p > 0) {
            /* H \varepsilon_t */
            gretl_matrix_multiply_mod(K->H, GRETL_MOD_NONE,
                                      et, GRETL_MOD_NONE,
                                      K->a1, GRETL_MOD_CUMULATE);
        } else {
            vector_from_row_mod(K->a1, U, K->t, GRETL_MOD_CUMULATE);
        }

        if (K->mu != NULL) {
            gretl_matrix_add_to(K->a1, K->mu);
        }

        /* record the (t+1)-dated state? */
        if (S != NULL && tmin == 1) {
            record_to_row_offset(S, K->a1, K->t, 0);
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
            ncols = K->VY == NULL ? K->r : K->r + K->n;
        }

        if (U->cols != ncols) {
            pprintf(prn, _("U should have %d columns but has %d\n"),
                    ncols, U->cols);
            err = E_NONCONF;
        }
    }

    if (!err && Sim0 != NULL) {
        int r = ssfsim ? K->r + 1 : K->r;
        int c = ssfsim ? K->r : 1;

        if (Sim0->rows != r || Sim0->cols != c) {
            pprintf(prn, _("simstart should be %d x %d, is %d x %d\n"),
                    r, c, Sim0->rows, Sim0->cols);
        }
    }

    if (!err && K->x != NULL) {
        /* do we have enough "obsx" data? */
        const gretl_matrix *X = SimX != NULL ? SimX : K->x;

        if (X->rows < U->rows) {
            pprintf(prn, _("obsx should have %d rows but has %d\n"),
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
    int i, j;

    for (j=0; j<m->cols; j++) {
        for (i=0; i<m->rows; i++) {
            if (i != j && gretl_matrix_get(m, i, j) != 0.0) {
                return 0;
            }
        }
    }

    return 1;
}

static int simdata_refresh_QR (kalman *K, PRN *prn)
{
    int err = 0;

    if (matrix_is_varying(K, K_VS) || matrix_is_varying(K, K_VY)) {
        err = kalman_update_matrices(K, prn);
    }

    return err;
}

/* Return a matrix in which the standard normal variates
   in @U are scaled according to K->VS, and K->VY if present.
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
        int n = K->VY == NULL ? 0 : K->n;
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

        if (matrix_is_varying(K, K_VS) || matrix_is_varying(K, K_VY)) {
            varying = 1;
        }

        K->b = b;
        set_kalman_running(K);

        if (matrix_is_diagonal(K->VS) &&
            (K->VY == NULL || matrix_is_diagonal(K->VY))) {
            for (t=0; t<T && !*err; t++) {
                if (varying) {
                    K->t = t;
                    *err = simdata_refresh_QR(K, prn);
                }
                for (j=0; j<rn && !*err; j++) {
                    if (j < K->r) {
                        vjj = gretl_matrix_get(K->VS, j, j);
                    } else {
                        vjj = gretl_matrix_get(K->VY, j-K->r, j-K->r);
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
                        gretl_matrix_inscribe_matrix(V, K->VS, 0, 0,
                                                     GRETL_MOD_NONE);
                        if (n > 0) {
                            gretl_matrix_inscribe_matrix(V, K->VY, K->r, K->r,
                                                         GRETL_MOD_NONE);
                        }
                        *err = gretl_matrix_psd_root(V, 0);
                        if (*err) {
                            gretl_errmsg_set(_("Failed to compute factor of Omega_t"));
                        } else {
                            load_from_row(Ut, U, t);
                            gretl_matrix_multiply_mod(Ut, GRETL_MOD_NONE,
                                                      V, GRETL_MOD_TRANSPOSE,
                                                      Et, GRETL_MOD_NONE);
                            record_to_row(E, Et, t);
                        }
                    }
                }

                gretl_matrix_free(Ut);
                gretl_matrix_free(Et);
            } else {
                gretl_matrix_inscribe_matrix(V, K->VS, 0, 0,
                                             GRETL_MOD_NONE);
                if (n > 0) {
                    gretl_matrix_inscribe_matrix(V, K->VY, K->r, K->r,
                                                 GRETL_MOD_NONE);
                }
                *err = gretl_matrix_psd_root(V, 0);
                if (*err) {
                    gretl_errmsg_set(_("Failed to compute factor of Omega"));
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

    if (id == K_ZT || id == K_VY || id == K_T || id == K_VS ||
        id == K_BT || id == K_m || id == K_a || id == K_P || id == K_R) {
        if (repl->rows != orig->rows || repl->cols != orig->cols) {
            err = E_DATA;
        }
    } else if (id == K_y || id == K_x) {
        if (repl->cols != orig->cols) {
            err = E_DATA;
        }
    }

    if (err) {
        gretl_errmsg_set(_("You cannot resize a state-space system matrix"));
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
                fprintf(stderr, " bad kalman matrix replacement!\n");
                return err;
            }
            fast_copy_values(*targ, src);
            if (!copy) {
                gretl_matrix_free(src);
            }
        } else if (copy) {
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

static gretl_matrix **get_input_matrix_target_by_id (kalman *K, int i)
{
    gretl_matrix **targ = NULL;

    if (i == K_T) {
        targ = &K->T;
    } else if (i == K_BT) {
        targ = &K->BT;
    } else if (i == K_ZT) {
        targ = &K->ZT;
    } else if (i == K_VS) {
        /* variance of state */
        if (kalman_djvar(K)) {
            targ = &K->H;
        } else if (kalman_dkvar(K)) {
            targ = &K->Q;
        } else {
            targ = &K->VS;
        }
    } else if (i == K_VY) {
        /* variance of observable */
        if (kalman_djvar(K)) {
            targ = &K->G;
        } else {
            targ = &K->VY;
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
    } else if (i == K_R) {
        targ = &K->R;
    }

    return targ;
}

/* Try attaching a matrix to a Kalman bundle: similar to
   kalman_bundle_set_matrix() below, but allowing for the possibility
   that the @data input of type @vtype has to be converted first.
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

/* Called by kalman_deserialize() when reconstructing a kalman bundle
   from XML, and also by kalman_bundle_copy() when duplicating such a
   bundle. In these contexts we know we have a gretl_matrix on input.
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
        pm = &K->vt;
    }

    return pm;
}

static const char *kalman_output_matrix_names[] = {
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

#define K_N_OUTPUTS G_N_ELEMENTS(kalman_output_matrix_names)

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

#define K_N_SCALARS 15

enum {
    Ks_t = 0,
    Ks_DIFFUSE,
    Ks_EXACT,
    Ks_CROSS,
    Ks_DKVAR,
    Ks_UNI,
    Ks_S2,
    Ks_LNL,
    Ks_r,
    Ks_n,
    Ks_N,
    Ks_p,
    Ks_d,
    Ks_j,
    Ks_extra
};

static const char *kalman_output_scalar_names[K_N_SCALARS] = {
    "t",
    "diffuse",
    "exact",
    "cross",
    "dkvar",
    "univariate",
    "s2",
    "lnl",
    "r",
    "n",
    "N",
    "p",
    "d",
    "j",
    "extra"
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
    case Ks_EXACT:
        retval[idx] = K->exact;
        break;
    case Ks_CROSS:
        retval[idx] = (K->vartype == DJ_VAR);
        break;
    case Ks_DKVAR:
        retval[idx] = (K->vartype == DK_VAR);
        break;
    case Ks_UNI:
        retval[idx] = (K->code == K_UNIVAR)? 1 : 0;
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
    case Ks_j:
        retval[idx] = K->j;
        break;
    case Ks_extra:
        retval[idx] = (K->flags & KALMAN_EXTRA)? 1 : 0;
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
    } else if (i == K_VS) {
        if (kalman_djvar(K)) {
            m = K->H;
        } else if (kalman_dkvar(K)) {
            m = K->Q;
        } else {
            m = K->VS;
        }
    } else if (i == K_VY) {
        if (kalman_djvar(K)) {
            m = K->G;
        } else {
            m = K->VY;
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
    } else if (i == K_R) {
        m = K->R;
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

/* Respond to "diffuse" setting: 0 = off, 1 = traditional, 2 = exact */

static int kalman_set_diffuse (kalman *K, int d)
{
    /* temporary testing thing */
    if (exact_set && d > 0) {
	d = 2;
    }

    if (d != 0 && d != 1 && d != 2) {
        return E_INVARG;
    } else if (d > 0) {
	K->flags |= KALMAN_DIFFUSE;
	if (d == 2) {
	    /* exact initial diffuse mode */
	    K->exact = 1;
	    if (K->code == K_LEGACY) {
		/* the legacy code doesn't handle this */
		K->code = kalman_djvar(K) ? K_DEJONG : K_UNIVAR;
	    }
	} else {
	    /* traditional "kappa" mode */
	    K->exact = 0;
	}
        return diffuse_Pini(K);
    } else {
        K->exact = 0;
        K->flags &= ~KALMAN_DIFFUSE;
        if (K->Pk0 != NULL) {
            gretl_matrix_free(K->Pk0);
            K->Pk0 = NULL;
        }
        return 0;
    }
}

/* respond to "univariate" or "dejong" setting */

static int kalman_set_code (kalman *K, int code, int s)
{
    if (s && code == K_UNIVAR && K->vartype == DJ_VAR) {
        gretl_errmsg_set(_("kalman: the 'univariate' setting is not compatible with\n"
                         "cross-correlated disturbances"));
        return E_INVARG;
    } else if (s && (code != K->code)) {
        K->code = code;
#if 0 /* this should not be required */
	if (code == K_UNIVAR && kalman_diffuse(K) && !K->exact) {
	    K->exact = 1;
	}
#endif
	if (code == K_UNIVAR && kalman_diffuse(K) && K->Pk0 == NULL) {
	    K->Pk0 = gretl_identity_matrix_new(K->r);
	    if (K->Pk0 == NULL) {
		return E_ALLOC;
	    }
	}
    }

    return 0;
}

/* Called by real_bundle_set_data() in gretl_bundle.c.  The return
   value indicates whether the putative setting was handled (1) or not
   (0). Not being handled here is not necessarily an error.

   The @copy flag here is inherited from the specific caller of
   real_bundle_set_data(): @copy = 1 if the caller was
   gretl_bundle_set_data(), 0 if it was gretl_bundle_donate_data().
   Either way the kalman struct takes ownership.
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

    if (kdebug) {
	fprintf(stderr, "maybe_set_kalman_element: key '%s', type %s\n",
		key, gretl_type_get_name(vtype));
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

    if (!strcmp(key, "diffuse") || !strcmp(key, "univariate") ||
        !strcmp(key, "dejong") || !strcmp(key, "extra")) {
        /* scalar config settings */
        if (vtype == GRETL_TYPE_DOUBLE) {
            double v = *(double *) vptr;

            if (!strcmp(key, "diffuse")) {
                *err = kalman_set_diffuse(K, (int) v);
            } else if (!strcmp(key, "univariate")) {
                *err = kalman_set_code(K, K_UNIVAR, (int) v);
	    } else if (!strcmp(key, "dejong")) {
		*err = kalman_set_code(K, K_DEJONG, (int) v);
            } else {
                K->flags |= KALMAN_EXTRA;
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
		gretl_errmsg_sprintf(_("%s: expected matrix but got %s"),
				     key, gretl_type_get_name(vtype));
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
            gretl_errmsg_sprintf(_("The member %s is read-only"), key);
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
        gretl_errmsg_sprintf(_("%s: cannot be deleted"), key);
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

/* Take the vec representations of the smoothed disturbances
   K->b.etahat and K->b.epshat (if present) and write them into
   K->U in the format documented in the Guide.
*/

static gretl_matrix *fill_smdist (kalman *K, int *ownit)
{
    gretl_matrix *U1 = gretl_bundle_get_matrix(K->b, "etahat", NULL);
    gretl_matrix *U2 = gretl_bundle_get_matrix(K->b, "epshat", NULL);
    gretl_matrix *U = K->U;

    if (U1 != NULL) {
        int j, t, q = U1->cols;
        int n = (U2 != NULL)? U2->cols : 0;
        double utj;

        if (U == NULL || U->rows != K->N || U->cols != q + n) {
            /* we need new storage */
            U = gretl_matrix_alloc(K->N, q + n);
            if (U == NULL) {
                return NULL;
            } else {
                *ownit = 1;
            }
        }
        for (t=0; t<K->N; t++) {
            for (j=0; j<q; j++) {
                utj = gretl_matrix_get(U1, t, j);
                gretl_matrix_set(U, t, j, utj);
            }
            for (j=0; j<n; j++) {
                utj = gretl_matrix_get(U2, t, j);
                gretl_matrix_set(U, t, j+q, utj);
            }
        }
        return U;
    }

    return NULL;
}

static gretl_matrix *fill_smdisterr (kalman *K, int *ownit)
{
    gretl_matrix *V1 = gretl_bundle_get_matrix(K->b, "veta", NULL);
    gretl_matrix *V2 = gretl_bundle_get_matrix(K->b, "veps", NULL);
    gretl_matrix *V = K->Vsd;

    if (V1 == NULL) {
        return NULL;
    } else {
        int q = kalman_dkvar(K) ? K->q : K->r;
        int n = (V2 != NULL)? K->n : 0;
        int j, k, t;
        double vtj;

        if (V == NULL || V->rows != K->N || V->cols != q + n) {
            /* we need new storage */
            V = gretl_matrix_alloc(K->N, q + n);
            if (V == NULL) {
                return NULL;
            } else {
                *ownit = 1;
            }
        }
        for (t=0; t<K->N; t++) {
            for (j=0, k=0; j<q; j++) {
                vtj = gretl_matrix_get(V1, t, k);
                gretl_matrix_set(V, t, j, sqrt(vtj));
                k += q + 1;
            }
            for (j=0; j<n; j++) {
                vtj = gretl_matrix_get(V2, t, j);
                gretl_matrix_set(V, t, j+q, sqrt(vtj));
            }
        }
        return V;
    }
}

static int mv_transform (kalman *K)
{
    gretl_matrix *Vt = gretl_zero_matrix_new(K->n, 1);
    gretl_matrix *PM = gretl_zero_matrix_new(K->n, K->r);
    gretl_matrix *Z  = K->uinfo->Z;
    gretl_matrix *At = K->a0;
    gretl_matrix *Ft = K->Ft;
    gretl_matrix *Pt = K->P0;
    int t;

    if (Vt == NULL || PM == NULL) {
        return E_ALLOC;
    }

    for (t=0; t<K->N; t++) {
        if (matrix_is_varying(K, K_ZT)) {
            K->t = t;
            retrieve_Zt(K);
        }
        load_from_row(At, K->A, t);
        load_from_row(Vt, K->y, t);
        gretl_matrix_multiply_mod(Z,  GRETL_MOD_NONE,
                                  At, GRETL_MOD_NONE,
                                  Vt, GRETL_MOD_DECREMENT);
        record_to_row(K->V, Vt, t);
        load_from_vech(Pt, K->P, K->r, t);
        gretl_matrix_multiply(Z, Pt, PM);
        gretl_matrix_copy_values(Ft, K->VY);
        gretl_matrix_multiply_mod(PM, GRETL_MOD_NONE,
                                  Z,  GRETL_MOD_TRANSPOSE,
                                  Ft, GRETL_MOD_CUMULATE);
        record_to_vech(K->F, Ft, K->n, t);
    }

    gretl_matrix_free(Vt);
    gretl_matrix_free(PM);

    return 0;
}

static gretl_matrix *construct_kalman_matrix (kalman *K,
					      const char *key,
					      int *ownit)
{
    gretl_matrix *m = NULL;

    if (m == NULL && K->code != K_LEGACY) {
	if (!strcmp(key, "smdist")) {
	    m = fill_smdist(K, ownit);
	} else if (!strcmp(key, "smdisterr")) {
	    m = fill_smdisterr(K, ownit);
	}
    }

    return m;
}

void *maybe_retrieve_kalman_element (void *kptr,
                                     const char *key,
                                     GretlType *type,
                                     int *reserved,
                                     int *ownit,
                                     int *err)
{
    kalman *K = kptr;
    void *ret = NULL;
    int i, id = -1;

    if (K == NULL) {
        *err = E_DATA;
        return NULL;
    }

    *type = GRETL_TYPE_NONE;
    if (ownit != NULL) {
        *ownit = 0;
    }

    /* 2022-05-06: we'll want to set *ownit to 1 if the element in
       question is newly allocated, and will be 'owned' by the caller,
       rather than just being a pointer to something inside the kalman
       struct.  If @ownit is NULL we'll take it that this call is on
       behalf of gretl_bundle_get_member_type(), in which case there's
       no need to construct a 'constructable' element.
    */

    if (kdebug > 1) {
	fprintf(stderr, "maybe_retrieve_kalman_element: key '%s'\n", key);
	fprintf(stderr, " ownit %p, K->b %p, univariate %d\n",
                (void *) ownit, (void *) K->b, kalman_univariate(K));
    }

#if 0
    if (!strcmp(key, "pevar") && kalman_univariate(K) && K->n > 1) {
        /* special case: the univariate F-matrix is a vector */
        *reserved = 1;
        *type = GRETL_TYPE_MATRIX;
        if (ownit != NULL) {
            ret = fill_pevar(K, ownit);
        }
        goto finalize;
    }
#endif

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
	    gretl_matrix *m = NULL;

            if (pm != NULL && *pm != NULL) {
                m = *pm;
            } else if (pm != NULL && ownit != NULL) {
                m = construct_kalman_matrix(K, key, ownit);
	    }
	    if (m != NULL) {
		*reserved = 1;
		ret = m;
		*type = GRETL_TYPE_MATRIX;
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

    // finalize:

    if (*reserved && ret == NULL && ownit != NULL) {
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
        pputs(prn, _("Kalman struct: empty\n"));
        err = E_DATA;
    } else {
        const gretl_matrix *m;
        gretl_matrix **pm;
        double *px;
        const char *name;
        int i, id;

        pputs(prn, _("\nKalman input matrices\n"));

        for (i=0; i<K_MMAX; i++) {
            id = K_input_mats[i].sym;
            m = k_input_matrix_by_id(K, id);
            if (m != NULL) {
                pprintf(prn, " %s: ", K_input_mats[i].name);
                pprintf(prn, "%d x %d\n", m->rows, m->cols);
            }
        }

        if (output_matrix_count(K) > 0) {
            pputs(prn, _("\nKalman output matrices\n"));
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

        pputs(prn, _("\nKalman scalars\n"));

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

/* For use in context of a kalman bundle: serialize the information in
   the kalman struct to XML
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

/* For use in context of a kalman bundle: deserialize the kalman
   struct from XML
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
    int vtype = 0;
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
                            vtype = DJ_VAR;
                        } else if (!strcmp(key, "dkvar") && x > 0) {
                            vtype = DK_VAR;
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

    if (vtype == DK_VAR) {
        if (Mreq[K_VS] == NULL || Mopt[K_R] == NULL) {
            *err = E_DATA;
            goto bailout;
        } else {
            Mopt[K_VY] = Mreq[4];
            Mreq[4] = Mopt[K_R];
            Mopt[K_R] = NULL;
        }
    } else if (vtype == STD_VAR && nmats == 5) {
        /* drop obsvar from initialization */
        Mopt[K_VY] = Mreq[4];
        Mreq[4] = NULL;
        nmats--;
    }

    if ((vtype > 0 && nmats != 5) || (vtype == 0 && nmats != 4)) {
        *err = E_DATA;
    } else {
        b = kalman_bundle_new(Mreq, copy, nmats, vtype == DK_VAR, err);
        if (b != NULL) {
            kalman *K = gretl_bundle_get_private_data(b);
            gretl_matrix **pm;
            const char *name;

            K->flags = Kflags;
            K->vartype = vtype;
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

 bailout:

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
    int dkopt = 0;
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
    if (kalman_djvar(K)) {
        M[3] = K->H;
        M[4] = K->G;
        k = 5;
    } else if (kalman_dkvar(K)) {
        M[3] = K->Q;
        M[4] = K->R;
        k = 5;
        dkopt = 1;
    } else {
        M[3] = K->VS;
    }

    b = kalman_bundle_new(M, copy, k, dkopt, err);

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
    Knew->vartype = K->vartype;

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

    *ns = K_N_SCALARS - na(K->s2) - na(K->loglik);
    S = strings_array_new(*ns);

    if (S != NULL) {
        int i, j = 0;

        for (i=0; i<K_N_SCALARS; i++) {
            if (i == Ks_S2 && na(K->s2)) {
                continue;
            } else if (i == Ks_LNL && na(K->loglik)) {
                continue;
            } else {
                S[j++] = gretl_strdup(kalman_output_scalar_names[i]);
            }
        }
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

/* Start of "univariate" filtering and smoothing code. Based on Durbin
   and Koopman, Time Series Analysis by State Space Methods, 2e
   (Oxford, 2012), section 6.4 in particular. The C implementation
   here owes a lot to Jouni Helske's fortran code in his KFAS package
   for R (see kfilter2.f90). Some helper functions come first.
*/

static void dj_from_Finf (const gretl_matrix *Finf, int *d, int *j,
			  double tiny)
{
    int cmax = Finf->cols - 1;
    int rmax = Finf->rows - 1;
    int i, k;

    *d = *j = 0;

    for (k=cmax; k >= 0 && *d == 0; k--) {
        for (i=rmax; i>=0; i--) {
            if (gretl_matrix_get(Finf, i, k) > tiny) {
                *d = k + 1;
                *j = i + 1;
                break;
            }
        }
    }
}

static int bundle_add_matrix (gretl_bundle *b,
                              const char *key,
                              gretl_matrix *m)
{
    return gretl_bundle_donate_data(b, key, m, GRETL_TYPE_MATRIX, 0);
}

/* simple implementation assuming no error check required */

static double dotprod (const gretl_vector *a, const gretl_vector *b)
{
    int i, n = gretl_vector_get_length(a);
    double dp = 0;

    for (i=0; i<n; i++) {
        dp += a->val[i] * b->val[i];
    }

    return dp;
}

static void increment_state (gretl_matrix *a,
                             const gretl_matrix *b,
                             gretl_matrix *tmp,
                             double x)
{
    fast_copy_values(tmp, b);
    gretl_matrix_multiply_by_scalar(tmp, x);
    gretl_matrix_add_to(a, tmp);
}

static void decrement_state_var (gretl_matrix *P,
                                 const gretl_matrix *K,
                                 gretl_matrix *tmp,
                                 double x)
{
    gretl_matrix_multiply_mod(K, GRETL_MOD_NONE,
                              K, GRETL_MOD_TRANSPOSE,
                              tmp, GRETL_MOD_NONE);
    gretl_matrix_multiply_by_scalar(tmp, x);
    gretl_matrix_subtract_from(P, tmp);
}

static void state_cross_update (gretl_matrix *Pti,
                                const gretl_matrix *Kti,
                                const gretl_matrix *Kki,
                                gretl_matrix *tmp,
                                double Fti, double Fkinv)
{
    gretl_matrix_multiply_mod(Kki, GRETL_MOD_NONE,
                              Kki, GRETL_MOD_TRANSPOSE,
                              tmp, GRETL_MOD_NONE);
    gretl_matrix_multiply_by_scalar(tmp, Fkinv * Fti * Fkinv);
    gretl_matrix_add_to(Pti, tmp);
    gretl_matrix_multiply_mod(Kti, GRETL_MOD_NONE,
                              Kki, GRETL_MOD_TRANSPOSE,
                              tmp, GRETL_MOD_NONE);
    gretl_matrix_multiply_mod(Kki, GRETL_MOD_NONE,
                              Kti, GRETL_MOD_TRANSPOSE,
                              tmp, GRETL_MOD_CUMULATE);
    gretl_matrix_multiply_by_scalar(tmp, Fkinv);
    gretl_matrix_subtract_from(Pti, tmp);
}

static int ensure_DK_smo_recorders (kalman *K, int Nd)
{
    if (K->uinfo->Finf == NULL) {
        K->uinfo->Finf = gretl_zero_matrix_new(K->n, Nd);
        K->uinfo->Kinf = gretl_zero_matrix_new(K->n * K->r, Nd);
        K->uinfo->Nd = Nd;
    }

    return 0;
}

static int kalman_ldl (kalman *K)
{
    gretl_matrix *L = K->uinfo->Linv;
    gretl_matrix *g = K->uinfo->g;
    double gj, lij;
    int i, j, err;

    fast_copy_values(L, K->VY);
    err = gretl_matrix_cholesky_decomp(L);

    if (!err) {
        for (j=0; j<K->n; j++) {
            gj = gretl_matrix_get(L, j, j);
            g->val[j] = gj * gj;
            for (i=j; i<K->n; i++) {
                lij = gretl_matrix_get(L, i, j);
                gretl_matrix_set(L, i, j, lij / gj);
            }
        }
        /* record L so we can undo the transformation */
        gretl_bundle_set_matrix(K->b, "L", L);
        err = gretl_invert_triangular_matrix(L, 'L');
    }

    return err;
}

static int kalman_diagonalize (kalman *K)
{
    gretl_matrix *tmp = gretl_matrix_alloc(K->n, 1);
    gretl_matrix *yt = gretl_matrix_alloc(K->n, 1);
    int t, err;

    err = kalman_ldl(K);
    if (err) {
        gretl_matrix_free(tmp);
        gretl_matrix_free(yt);
        return err;
    }

    /* adjust Z */
    gretl_matrix_multiply_mod(K->uinfo->Linv, GRETL_MOD_NONE,
                              K->ZT, GRETL_MOD_TRANSPOSE,
                              K->uinfo->Z, GRETL_MOD_NONE);

    if (K->uinfo->y0 == NULL) {
        K->uinfo->y0 = gretl_matrix_copy(K->y);
        if (K->uinfo->y0 == NULL) {
            err = E_ALLOC;
        }
    }

    if (!err) {
        /* adjust the observable */
        for (t=0; t<K->N; t++) {
            load_from_row(tmp, K->uinfo->y0, t);
            gretl_matrix_multiply(K->uinfo->Linv, tmp, yt);
            record_to_vec(K->y, yt, t);
        }
    }

    gretl_matrix_free(tmp);
    gretl_matrix_free(yt);

    return err;
}

static int handle_univariate_transforms (kalman *K)
{
    int err = 0;

    err = gretl_matrix_transpose(K->uinfo->Z, K->ZT);

    if (K->n > 1 && K->VY != NULL) {
        /* obsvar needs transformation */
        int i;

        if (matrix_is_diagonal(K->VY)) {
            for (i=0; i<K->n; i++) {
                K->uinfo->g->val[i] = gretl_matrix_get(K->VY, i, i);
            }
        } else {
            if (kdebug) {
                fprintf(stderr, "kfilter: diagonalizing\n");
            }
            err = kalman_diagonalize(K);
        }
    }

    return err;
}

static double uni_get_vti (kalman *K, const gretl_matrix *Zi,
                           const gretl_matrix *ati, int i)
{
    double vti = gretl_matrix_get(K->y, K->t, i);

    if (!na(vti)) {
        if (K->BT != NULL) {
            vti -= kalman_Bx_uni(K, i);
        }
        if (!na(vti)) {
            vti -= dotprod(Zi, ati);
        }
    }

    return vti;
}

static double tiny_value (const gretl_matrix *Z)
{
    int i, n = Z->rows * Z->cols;
    double x, zmin = 1.0e+30;
    double macheps = 2.20e-16;

    for (i=0; i<n; i++) {
        x = fabs(Z->val[i]);
        if (x > 0 && x < zmin) {
            zmin = x;
        }
    }

    x = sqrt(macheps) * zmin * zmin;
#if 0
    fprintf(stderr, "Computed tiny value: %g (zmin = %g)\n", x, zmin);
    return K_TINY;
#else
    return x;
#endif
}

static int kfilter_univariate (kalman *K, PRN *prn)
{
    univar_info *ui = K->uinfo;
    gretl_matrix *Kkt = K->Mt;
    gretl_matrix *at = K->a0;
    gretl_matrix *Pt = K->P0;
    gretl_matrix *ati = K->a1;
    gretl_matrix *Pti = K->P1;
    gretl_matrix *Pk = K->Pk0;
    gretl_matrix *UF = K->F;

    gretl_matrix *Z   = ui->Z;
    gretl_matrix *g   = ui->g;
    gretl_matrix *Zi  = ui->Zi;
    gretl_matrix *Kti = ui->Kti;
    gretl_matrix *m1  = ui->m1;
    gretl_matrix *mm  = ui->mm;
    gretl_matrix *Pki = NULL;
    gretl_matrix *Kki = NULL;

    gretl_matrix *att = NULL;
    gretl_matrix *Ptt = NULL;

    double vti, Fti, llct;
    double Ftinv, Fkinv, Fki = 0;
    double qt, ldt;
    double k_tiny = 0;

    int exact_diffuse = kalman_diffuse(K) && K->exact;
    int rankPk = exact_diffuse ? K->r : 0;
    int d = K->exact ? 0 : -1;
    int j = K->exact ? 0 : -1;
    int Nd = 0;
    int t, i;
    int err = 0;

    if (identify || trace) {
        printf("kfilter_univariate(), diffuse = %d, exact = %d\n",
	       kalman_diffuse(K), K->exact);
    }

    if (K->exact) {
        Pki = ui->Pki;
        Kki = ui->Kki;
    }

    if (K->n > 1 && ) {
        err = kalman_ensure_stepinfo(K, 0);
        if (err) {
            return err;
        }
        UF = ui->UF;
    }

    if (K->n > 1 || K->r > 1) {
        handle_univariate_transforms(K);
    }

    if (K->smo_prep && K->exact) {
        Nd = matrix_is_varying(K, K_ZT) ? K->N : 4*K->r;
        ensure_DK_smo_recorders(K, Nd);
    }

    if (kalman_extra(K)) {
        att = gretl_matrix_alloc(K->N, K->r);
        Ptt = gretl_matrix_alloc(K->N, K->r * K->r);
    }

    K->TI = gretl_is_identity_matrix(K->T);
    K->okN = K->n * K->N;
    K->SSRw = 0;
    K->loglik = 0;
    set_kalman_running(K);

    if (trace) {
        printf("kfilter_univariate(), exact = %d\n", K->exact);
    }
    if (kdebug) {
        fprintf(stderr, "\n*** kfilter_univariate: r=%d, n=%d, exact=%d ***\n",
                K->r, K->n, K->exact);
    }

    if (K->exact) {
        k_tiny = tiny_value(Z);
    }

    for (t=0; t<K->N; t++) {
        K->t = t;
        record_to_vec(K->A, at, t);
        record_to_vech(K->P, Pt, K->r, t);
        fast_copy_values(ati, at);
        fast_copy_values(Pti, Pt);
        if (Pki != NULL) {
            fast_copy_values(Pki, Pk);
        }
        if (filter_is_varying(K)) {
            /* we have time-varying matrices */
            err = kalman_refresh_matrices(K, prn);
            if (err) {
                K->loglik = NADBL;
                break;
            }
        }
        qt = 0;
        ldt = 0;
        llct = 0;
        if (K->smo_prep && d == 0 && t < Nd) {
            record_to_col(K->PK, Pk, t);
        }
        for (i=0; i<K->n; i++) {
            if (kdebug && t < 2 && K->exact) {
                fprintf(stderr, "t,i = %d,%d, rankPk = %d\n", t, i, rankPk);
            }
            load_from_row(Zi, Z, i);
            gretl_matrix_multiply_mod(Pti, GRETL_MOD_NONE,
                                      Zi, GRETL_MOD_TRANSPOSE,
                                      Kti, GRETL_MOD_NONE);
            Fti = dotprod(Zi, Kti);
            if (g != NULL) {
                Fti += g->val[i];
            }
            gretl_matrix_set(UF, t, i, Fti);
            record_to_col(K->Kt, Kti, i);
            if (d == 0) {
                /* still initial */
                gretl_matrix_multiply_mod(Pki, GRETL_MOD_NONE,
                                          Zi, GRETL_MOD_TRANSPOSE,
                                          Kki, GRETL_MOD_NONE);
                Fki = dotprod(Zi, Kki);
                if (K->smo_prep && t < Nd) {
                    gretl_matrix_set(ui->Finf, i, t, Fki);
                    record_to_col(Kkt, Kki, i);
                }
            }

            vti = uni_get_vti(K, Zi, ati, i);
            gretl_matrix_set(K->V, t, i, vti);
            if (na(vti)) {
                K->okN -= 1;
                continue;
            }

            if (rankPk > 0) {
                if (Fki > k_tiny) {
                    K->okN -= 1;
                    Fkinv = 1.0 / Fki;
                    ldt += log(Fki);
                    increment_state(ati, Kki, m1, Fkinv * vti);
                    state_cross_update(Pti, Kti, Kki, mm, Fti, Fkinv);
                    decrement_state_var(Pki, Kki, mm, Fkinv);
                    --rankPk;
                } else if (Fti > 0) {
                    Ftinv = 1.0 / Fti;
                    qt += vti * vti * Ftinv;
                    ldt += log(Fti);
                    llct += LN_2_PI;
                    increment_state(ati, Kti, m1, Ftinv * vti);
                    decrement_state_var(Pti, Kti, mm, Ftinv);
                }
            } else {
                /* rankPk = 0 */
                if (j == 0) {
                    j = (K->n == 1)? 1 : i;
                    if (kdebug) {
                        fprintf(stderr, "*** filter: got j = %d (i=%d)\n", j, i);
                    }
                }
                Ftinv = 1/Fti;
                qt += vti * vti * Ftinv;
                ldt += log(Fti);
                llct += LN_2_PI;
                increment_state(ati, Kti, m1, Ftinv * vti);
                decrement_state_var(Pti, Kti, mm, Ftinv);
            }
        } /* observables at time t */

	if (K->K != NULL) {
	    /* record the gain */
	    record_to_vec(K->K, K->Kt, t);
	}
        if (K->smo_prep && d == 0 && t < Nd) {
	    /* record K∞ */
	    record_to_col(ui->Kinf, Kkt, t);
        }

        if (j > 0 && d == 0) {
            d = t + 1; /* note: 1-based */
            if (kdebug) {
                fprintf(stderr, "*** filter: got d = %d\n", d);
            }
        }

        /* update for t+1 */
        if (K->TI) {
            /* shortcut in case T = I_m */
            fast_copy_values(at, ati);
            fast_copy_values(Pt, Pti);
            gretl_matrix_add_to(Pt, K->VS);
            if (Pki != NULL) {
                fast_copy_values(Pk, Pki);
            }
        } else {
            gretl_matrix_multiply(K->T, ati, at);
            fast_copy_values(Pt, K->VS);
            gretl_matrix_qform(K->T, GRETL_MOD_NONE, Pti,
                               Pt, GRETL_MOD_CUMULATE);
            if (Pki != NULL) {
                gretl_matrix_qform(K->T, GRETL_MOD_NONE, Pki,
                                   Pk, GRETL_MOD_NONE);
            }
        }
        if (K->mu != NULL) {
            gretl_matrix_add_to(at, K->mu);
        }

        if (K->LL != NULL) {
            K->LL->val[t] = -0.5 * (llct + ldt + qt);
            K->loglik += K->LL->val[t];
        } else {
            K->loglik -= 0.5 * (llct + ldt + qt);
        }
        K->SSRw += qt;
        if (att != NULL) {
            record_to_vec(att, ati, t);
        }
        if (Ptt != NULL) {
            record_to_vec(Ptt, Pti, t);
        }
    } /* end of time loop */

    if (0) { /* Hmm, conditionality? */
        K->d = K->N;
        K->j = K->n;
    } else if (K->smo_prep && K->exact) {
        dj_from_Finf(ui->Finf, &d, &j, k_tiny);
#if 0 /* not sure about this, causes trouble? */
        if (d > 0 && d < ui->Finf->cols) {
            gretl_matrix_realloc(ui->Finf, ui->Finf->rows, d);
            gretl_matrix_realloc(ui->Kinf, ui->Kinf->rows, d);
        }
#endif
    }

    if (kalman_diffuse(K) && K->exact) {
        K->d = d > 0 ? d : 0;
        K->j = j > 0 ? j : 0;
    }

    set_kalman_stopped(K);

    K->s2 = K->SSRw / K->okN;

    if (kdebug) {
        fprintf(stderr, "*** after filtering: d=%d, j=%d ***\n", d, j);
        fprintf(stderr, "univariate: loglik = %#.8g\n", K->loglik);
    }

    if (att != NULL) {
        bundle_add_matrix(K->b, "att", att);
    }
    if (Ptt != NULL) {
        bundle_add_matrix(K->b, "Ptt", Ptt);
    }

    if (K->n > 1) {
        mv_transform(K);
    }

    return 0;
}

/* univariate smoothing follows */

/* cumulant matrices */

struct cumulants {
    gretl_matrix *r0;
    gretl_matrix *N0;
    gretl_matrix *L0;
    gretl_matrix *r1;
    gretl_matrix *N1;
    gretl_matrix *N2;
    gretl_matrix *Nt;
    gretl_matrix *mm;
    gretl_matrix *m1;
};

static int allocate_cumulants (struct cumulants *c,
                               kalman *K)
{
    int r = K->r;

    c->r0 = gretl_zero_matrix_new(r, 1);
    c->N0 = gretl_zero_matrix_new(r, r);
    c->L0 = gretl_matrix_alloc(r, r);
    c->Nt = gretl_zero_matrix_new(r, r);
    c->mm = gretl_matrix_alloc(r, r);
    c->m1 = gretl_matrix_alloc(r, 1);
    if (K->exact) {
        c->r1 = gretl_zero_matrix_new(r, 1);
        c->N1 = gretl_zero_matrix_new(r, r);
        c->N2 = gretl_zero_matrix_new(r, r);
    }

    return 0;
}

static void clear_cumulants (struct cumulants *c,
                             kalman *K)
{
    gretl_matrix_free(c->r0);
    gretl_matrix_free(c->N0);
    gretl_matrix_free(c->L0);
    gretl_matrix_free(c->Nt);
    gretl_matrix_free(c->mm);
    gretl_matrix_free(c->m1);
    if (K->exact) {
        gretl_matrix_free(c->r1);
        gretl_matrix_free(c->N1);
        gretl_matrix_free(c->N2);
    }
}

/* complex smoothing iteration when F∞ > 0 */

static void fkpos (double fkinv,
                   double Fti,
                   double vti,
                   const gretl_matrix *Kki,
                   const gretl_matrix *Kti,
                   const gretl_matrix *Zti,
                   struct cumulants *c)
{
    gretl_matrix *Linf, *Lmid;
    gretl_matrix *mm2;
    int m = Zti->cols;

    mm2 = gretl_matrix_alloc(m, m);

    /* Linf = I(m) - Kki * Zti * fkinv */
    Linf = gretl_identity_matrix_new(m);
    gretl_matrix_multiply(Kki, Zti, c->mm);
    gretl_matrix_multiply_by_scalar(c->mm, fkinv);
    gretl_matrix_subtract_from(Linf, c->mm);

    /* Lmid = Kki * Fti * fkinv - Kti */
    Lmid = gretl_matrix_copy(Kki);
    gretl_matrix_multiply_by_scalar(Lmid, Fti * fkinv);
    gretl_matrix_subtract_from(Lmid, Kti);

    /* L0 = Lmid * Zti * fkinv */
    gretl_matrix_multiply(Lmid, Zti, c->L0);
    gretl_matrix_multiply_by_scalar(c->L0, fkinv);

    /* r1 = Linf' * r1 + L0' * r0 */
    fast_copy_values(c->m1, c->r1);
    gretl_matrix_multiply_mod(Linf, GRETL_MOD_TRANSPOSE,
                              c->m1, GRETL_MOD_NONE,
                              c->r1, GRETL_MOD_NONE);
    gretl_matrix_multiply_mod(c->L0, GRETL_MOD_TRANSPOSE,
                              c->r0, GRETL_MOD_NONE,
                              c->r1, GRETL_MOD_CUMULATE);

    /* r1 += vti * fkinv * Zti' */
    gretl_matrix_copy_values_shaped(c->m1, Zti);
    gretl_matrix_multiply_by_scalar(c->m1, vti * fkinv);
    gretl_matrix_add_to(c->r1, c->m1);

    /* r0 = Linf' * r0 */
    gretl_matrix_multiply_mod(Linf, GRETL_MOD_TRANSPOSE,
                              c->r0, GRETL_MOD_NONE,
                              c->m1, GRETL_MOD_NONE);
    fast_copy_values(c->r0, c->m1);

    /* N2 = Linf' * N2 * Linf */
    gretl_matrix_multiply_mod(Linf, GRETL_MOD_TRANSPOSE,
                              c->N2, GRETL_MOD_NONE,
                              c->mm, GRETL_MOD_NONE);
    gretl_matrix_multiply(c->mm, Linf, c->N2);

    /* N2 -= Fti * fkinv * fkinv * Zti' * Zti */
    gretl_matrix_multiply_mod(Zti, GRETL_MOD_TRANSPOSE,
                              Zti, GRETL_MOD_NONE,
                              c->mm, GRETL_MOD_NONE);
    gretl_matrix_multiply_by_scalar(c->mm, Fti * fkinv * fkinv);
    gretl_matrix_subtract_from(c->N2, c->mm);

    /* N2 += L0' * N0 * L0 */
    gretl_matrix_multiply_mod(c->L0, GRETL_MOD_TRANSPOSE,
                              c->N0, GRETL_MOD_NONE,
                              c->mm, GRETL_MOD_NONE);
    gretl_matrix_multiply(c->mm, c->L0, mm2);
    gretl_matrix_add_to(c->N2, mm2);

    /* mm = Linf' * N1 * L0 */
    gretl_matrix_multiply_mod(Linf, GRETL_MOD_TRANSPOSE,
                              c->N1, GRETL_MOD_NONE,
                              mm2, GRETL_MOD_NONE);
    gretl_matrix_multiply(mm2, c->L0, c->mm);

    /* N2 += mm + mm' */
    gretl_matrix_add_to(c->N2, c->mm);
    gretl_matrix_add_transpose_to(c->N2, c->mm);

    /* N1 = Linf' * N1 * Linf */
    gretl_matrix_multiply(mm2, Linf, c->N1);

    /* N1 += fkinv * Zti' * Zti */
    gretl_matrix_multiply_mod(Zti, GRETL_MOD_TRANSPOSE,
                              Zti, GRETL_MOD_NONE,
                              c->mm, GRETL_MOD_NONE);
    gretl_matrix_multiply_by_scalar(c->mm, fkinv);
    gretl_matrix_add_to(c->N1, c->mm);

    /* N1 += L0' * N0 * Linf */
    gretl_matrix_multiply(c->N0, Linf, c->mm);
    gretl_matrix_multiply_mod(c->L0, GRETL_MOD_TRANSPOSE,
                              c->mm, GRETL_MOD_NONE,
                              c->N1, GRETL_MOD_CUMULATE);

    /* N0 = Linf' * N0 * Linf */
    gretl_matrix_multiply_mod(Linf, GRETL_MOD_TRANSPOSE,
                              c->N0, GRETL_MOD_NONE,
                              c->mm, GRETL_MOD_NONE);
    gretl_matrix_multiply(c->mm, Linf, c->N0);

    gretl_matrix_free(mm2);
    gretl_matrix_free(Linf);
    gretl_matrix_free(Lmid);
}

/* relatively simple smoothing iteration when F∞ = 0 */

static void fkzero (double ftinv,
                    double vti,
                    const gretl_matrix *Kti,
                    const gretl_matrix *Zti,
                    struct cumulants *c)
{
    /* L0 = I(m) - ftinv * Kti * Zti */
    fast_write_I(c->L0);
    gretl_matrix_multiply(Kti, Zti, c->mm);
    gretl_matrix_multiply_by_scalar(c->mm, ftinv);
    gretl_matrix_subtract_from(c->L0, c->mm);

    /* r0 = L0' * r0 + vti * ftinv * Zti' */
    gretl_matrix_multiply_mod(c->L0, GRETL_MOD_TRANSPOSE,
                              c->r0, GRETL_MOD_NONE,
                              c->m1, GRETL_MOD_NONE);
    fast_copy_values(c->r0, c->m1);
    gretl_matrix_copy_values_shaped(c->m1, Zti);
    gretl_matrix_multiply_by_scalar(c->m1, vti * ftinv);
    gretl_matrix_add_to(c->r0, c->m1);

    if (c->r1 != NULL) {
        /* r1 = L0' * r1 */
        fast_copy_values(c->m1, c->r1);
        gretl_matrix_multiply_mod(c->L0, GRETL_MOD_TRANSPOSE,
                                  c->m1, GRETL_MOD_NONE,
                                  c->r1, GRETL_MOD_NONE);
    }

    /* N0 = L0' * N0 * L0 + ftinv * Z[i,]' * Z[i,] */
    gretl_matrix_multiply_mod(c->L0, GRETL_MOD_TRANSPOSE,
                              c->N0, GRETL_MOD_NONE,
                              c->mm, GRETL_MOD_NONE);
    gretl_matrix_multiply(c->mm, c->L0, c->N0);
    gretl_matrix_multiply_mod(Zti, GRETL_MOD_TRANSPOSE,
                              Zti, GRETL_MOD_NONE,
                              c->mm, GRETL_MOD_NONE);
    gretl_matrix_multiply_by_scalar(c->mm, ftinv);
    gretl_matrix_add_to(c->N0, c->mm);

    if (c->N1 != NULL && c->N2 != NULL) {
        /* N1 = N1 * L0 */
        gretl_matrix_multiply(c->N1, c->L0, c->mm);
        fast_copy_values(c->N1, c->mm);
        /* N2 = N2 * L0 */
        gretl_matrix_multiply(c->N2, c->L0, c->mm);
        fast_copy_values(c->N2, c->mm);
    }
}

/* smoothing update of state and its variance */

static int dagger_calc (struct cumulants *c,
                        gretl_matrix *Pt,
                        gretl_matrix *Pk,
			kalman *K)
{
    gretl_matrix_block *B;
    gretl_matrix *rdag, *Ndag, *Pdag;
    gretl_matrix *at = c->m1;
    gretl_matrix *tmp;
    int m = c->N0->rows;
    int m2 = m * 2;
    int mm = m * m;
    size_t rsize = m * sizeof(double);
    size_t Nsize = mm * sizeof(double);

    /* workspace */
    B = gretl_matrix_block_new(&tmp, m, m2, &rdag, m2, 1,
                               &Ndag, m2, m2, &Pdag, m, m2,
                               NULL);
    if (B == NULL) {
        return E_ALLOC;
    }

    /* rdag = r0 | r1 */
    memcpy(rdag->val, c->r0->val, rsize);
    memcpy(rdag->val + m, c->r1->val, rsize);

    /* Ndag = (N0 ~ N1') | (N1 ~ N2) */
    gretl_matrix_inscribe_matrix(Ndag, c->N0, 0, 0, GRETL_MOD_NONE);
    gretl_matrix_inscribe_matrix(Ndag, c->N1, 0, m, GRETL_MOD_TRANSPOSE);
    gretl_matrix_inscribe_matrix(Ndag, c->N1, m, 0, GRETL_MOD_NONE);
    gretl_matrix_inscribe_matrix(Ndag, c->N2, m, m, GRETL_MOD_NONE);

    load_from_row(at, K->A, K->t);
    load_from_vech(Pt, K->P, K->r, K->t);
    load_from_col(Pk, K->PK, K->t);

    /* Pdag = Pt ~ Pk */
    memcpy(Pdag->val, Pt->val, Nsize);
    memcpy(Pdag->val + mm, Pk->val, Nsize);

#if 0
    fprintf(stderr, "dagger calc, t=%d\n", K->t);
    gretl_matrix_print(Pt, "Pt");
    gretl_matrix_print(Pdag, "Pdag");
#endif

    /* ahat[t,] = (at + Pdag * rdag)' */
    gretl_matrix_multiply_mod(Pdag, GRETL_MOD_NONE,
                              rdag, GRETL_MOD_NONE,
                              at, GRETL_MOD_CUMULATE);
    record_to_vec(K->A, at, K->t);

    /* Vhat[t,] = vec(Pt - Pdag * Ndag * Pdag')' */
    gretl_matrix_multiply(Pdag, Ndag, tmp);
    gretl_matrix_multiply_mod(tmp, GRETL_MOD_NONE,
                              Pdag, GRETL_MOD_TRANSPOSE,
                              Pt, GRETL_MOD_DECREMENT);
    record_to_vech(K->P, Pt, K->r, K->t);

    gretl_matrix_block_destroy(B);

    return 0;
}

static void eps_smooth_Ft (double ftinv, double gi,
                           double vti,
                           const gretl_matrix *Kti,
                           struct cumulants *c,
                           gretl_matrix *epshat,
                           gretl_matrix *veps,
                           int t, int i, int dist)
{
    double x = dotprod(Kti, c->r0);

    gretl_matrix_set(epshat, t, i, gi * ftinv * (vti - x));

    gretl_matrix_multiply(c->N0, Kti, c->m1);
    gretl_matrix_multiply_by_scalar(c->m1, ftinv * ftinv);
    if (dist == 2) {
        x = gi - gi*gi * (ftinv + dotprod(Kti, c->m1));
    } else {
        x = gi*gi * (ftinv + dotprod(Kti, c->m1));
    }
    gretl_matrix_set(veps, t, i, x);
}

static void eps_smooth_Fk (double fkinv, double gi,
                           const gretl_matrix *Kki,
                           struct cumulants *c,
                           gretl_matrix *epshat,
                           gretl_matrix *veps,
                           int t, int i, int dist)
{
    double x = dotprod(Kki, c->r0);

    gretl_matrix_set(epshat, t, i, -gi * x * fkinv);

    gretl_matrix_multiply(c->N0, Kki, c->m1);
    if (dist == 2) {
        x = gi - gi*gi * fkinv * fkinv * dotprod(Kki, c->m1);
    } else {
        x = gi*gi * fkinv * fkinv * dotprod(Kki, c->m1);
    }
    gretl_matrix_set(veps, t, i, x);
}

static void eta_smooth (const gretl_matrix *Q,
                        const gretl_matrix *QRT,
                        struct cumulants *c,
                        gretl_matrix *etahat,
                        gretl_matrix *veta,
                        int t, int dist)
{
    gretl_matrix *tmp;
    int q = Q->rows;
    int m = QRT->cols;

    tmp = gretl_matrix_reuse(c->mm, q, q);

    if (dist == 2) {
        /* Q - QRT * Nt * (QRT)' */
        fast_copy_values(tmp, Q);
        gretl_matrix_qform(QRT, GRETL_MOD_NONE, c->Nt,
                           tmp, GRETL_MOD_DECREMENT);
    } else {
        /* QRT * Nt * (QRT)' */
        gretl_matrix_qform(QRT, GRETL_MOD_NONE, c->Nt,
                           tmp, GRETL_MOD_NONE);
    }

    record_to_vec(veta, tmp, t);

    if (t > 0) {
        /* QRT * r */
        gretl_matrix_reuse(tmp, -1, 1);
        gretl_matrix_multiply(QRT, c->r0, tmp);
        record_to_vec(etahat, tmp, t-1);
    }

    gretl_matrix_reuse(c->mm, m, m);
}

static void state_smooth (const gretl_matrix *at,
                          const gretl_matrix *Pt,
                          struct cumulants *c,
			  kalman *K)
{
    /* Pt - Pt * N0 * Pt' */
    fast_copy_values(c->mm, Pt);
    gretl_matrix_qform(Pt, GRETL_MOD_NONE, c->N0,
                       c->mm, GRETL_MOD_DECREMENT);
    record_to_vech(K->P, c->mm, K->r, K->t);

    /* at + Pt * r0 */
    gretl_matrix_multiply(Pt, c->r0, c->m1);
    gretl_matrix_add_to(c->m1, at);
    record_to_vec(K->A, c->m1, K->t);
}

static void regular_backdate (struct cumulants *c,
                              gretl_matrix *T)
{
    /* N0 = T' * N0 * T */
    gretl_matrix_multiply_mod(T, GRETL_MOD_TRANSPOSE,
                              c->N0, GRETL_MOD_NONE,
                              c->mm, GRETL_MOD_NONE);
    gretl_matrix_multiply(c->mm, T, c->N0);

    /* r0 = T' * r0 */
    gretl_matrix_multiply_mod(T, GRETL_MOD_TRANSPOSE,
                              c->r0, GRETL_MOD_NONE,
                              c->m1, GRETL_MOD_NONE);
    fast_copy_values(c->r0, c->m1);
}

static void diffuse_backdate (struct cumulants *c,
                              gretl_matrix *T)
{
    /* N0 = T' * N0 * T */
    gretl_matrix_multiply_mod(T, GRETL_MOD_TRANSPOSE,
                              c->N0, GRETL_MOD_NONE,
                              c->mm, GRETL_MOD_NONE);
    gretl_matrix_multiply(c->mm, T, c->N0);

    /* N1 = T' * N1 * T */
    gretl_matrix_multiply_mod(T, GRETL_MOD_TRANSPOSE,
                              c->N1, GRETL_MOD_NONE,
                              c->mm, GRETL_MOD_NONE);
    gretl_matrix_multiply(c->mm, T, c->N1);

    /* N2 = T' * N2 * T */
    gretl_matrix_multiply_mod(T, GRETL_MOD_TRANSPOSE,
                              c->N2, GRETL_MOD_NONE,
                              c->mm, GRETL_MOD_NONE);
    gretl_matrix_multiply(c->mm, T, c->N2);

    /* r0 = T' * r0 */
    gretl_matrix_multiply_mod(T, GRETL_MOD_TRANSPOSE,
                              c->r0, GRETL_MOD_NONE,
                              c->m1, GRETL_MOD_NONE);
    fast_copy_values(c->r0, c->m1);
    /* r1 = T' * r1 */
    gretl_matrix_multiply_mod(T, GRETL_MOD_TRANSPOSE,
                              c->r1, GRETL_MOD_NONE,
                              c->m1, GRETL_MOD_NONE);
    fast_copy_values(c->r1, c->m1);
}

static void kalman_bundle_add_epshat (kalman *K, gretl_matrix *epshat)
{
    gretl_matrix *L;

    /* add epshat "as is" */
    bundle_add_matrix(K->b, "epshat", epshat);

    L = gretl_bundle_get_matrix(K->b, "L", NULL);

    if (L != NULL) {
        /* The triangular matrix L will be present only when the
           variance matrix of eps is not diagonal. In that case
           we'll try adding estimates of the untransformed
           disturbances.
        */
        gretl_matrix *real_epshat;
        int T = epshat->rows;
        int k = epshat->cols;

        real_epshat = gretl_matrix_alloc(T, k);

        if (real_epshat != NULL) {
            double eti, lij, etj;
            int i, j, jmax, t;

            /* real_epshat[t,] = (L * epshat[t,]')' */
            for (t=0; t<T; t++) {
                jmax = 1;
                for (i=0; i<k; i++) {
                    eti = 0.0;
                    for (j=0; j<jmax; j++) {
                        lij = gretl_matrix_get(L, i, j);
                        etj = gretl_matrix_get(epshat, t, j);
                        eti += lij * etj;
                    }
                    gretl_matrix_set(real_epshat, t, i, eti);
                    jmax++;
                }
            }
            bundle_add_matrix(K->b, "real_epshat", real_epshat);
        }
    }
}

static int state_dist_setup (kalman *K,
                             gretl_matrix **Q,
                             gretl_matrix **QRT,
                             gretl_matrix **etahat,
                             gretl_matrix **veta)
{
    int q = K->r;

    if (K->Q != NULL) {
        q = K->Q->rows;
        *Q = K->Q;
        *QRT = K->QRT;
    } else {
        *Q = K->VS;
        *QRT = K->VS;
    }

    *etahat = gretl_zero_matrix_new(K->N, q);
    *veta = gretl_zero_matrix_new(K->N, q*q);

    return 1;
}

static int ksmooth_univariate (kalman *K, int dist)
{
    gretl_matrix_block *B = NULL;
    univar_info *ui = K->uinfo;
    struct cumulants c = {0};
    gretl_matrix *Q = NULL;
    gretl_matrix *QRT = NULL;
    gretl_matrix *UF = NULL;
    gretl_matrix *g = NULL;

    int i, t, N = K->N;
    int d = K->d;
    int j = K->j;
    int eps_smo = 0;
    int eta_smo = 0;
    int err = 0;

    /* matrices to be returned in bundle */
    gretl_matrix *epshat = NULL;
    gretl_matrix *veps = NULL;
    gretl_matrix *etahat = NULL;
    gretl_matrix *veta = NULL;

    /* block-workspace matrices */
    gretl_matrix *Pt, *Pk, *vt;
    gretl_matrix *Ft, *Kt, *at;
    gretl_matrix *Fk, *Kk;
    gretl_matrix *Zti, *Kti, *Kki;

    /* workspace scalars */
    double vti, Fti, Fki;
    double ftinv, fkinv;
    double k_tiny = 0;

    /* convenience pointers */
    UF = K->n > 1 ? ui->UF : K->F;
    g  = ui->g;

    if (trace) {
        printf("ksmooth_univariate(), dist = %d\n", dist);
    }
    if (kdebug) {
        fprintf(stderr, "ksmooth_univariate: dist=%d, d=%d, j=%d\n",
                dist, K->d, K->j);
    }

    if (K->exact) {
        B = gretl_matrix_block_new(&Kt, K->r, K->n,
				   &Pt, K->r, K->r,
                                   &Ft, K->n, 1,
				   &vt, K->n, 1,
                                   &at, K->r, 1,
				   &Fk, K->n, 1,
                                   &Kk, K->r, K->n,
				   &Pk, K->r, K->r,
                                   &Zti, 1, K->r,
				   &Kti, K->r, 1,
                                   &Kki, K->r, 1, NULL);
    } else {
        B = gretl_matrix_block_new(&Kt, K->r, K->n,
				   &Pt, K->r, K->r,
                                   &Ft, K->n, 1,
				   &vt, K->n, 1,
                                   &at, K->r, 1,
				   &Zti, 1, K->r,
                                   &Kti, K->r, 1, NULL);
    }

    allocate_cumulants(&c, K);

    if (dist) {
        if (g != NULL) {
            /* smoothing of obs disturbances wanted */
            epshat = gretl_zero_matrix_new(N, K->n);
            veps = gretl_zero_matrix_new(N, K->n);
            eps_smo = 1;
        }
        /* smoothing of state disturbances wanted */
        eta_smo = state_dist_setup(K, &Q, &QRT, &etahat, &veta);
    }

    if (K->exact) {
	k_tiny = tiny_value(ui->Z); /* ? */
    }

    if (kdebug) {
        fprintf(stderr, "\n*** univariate smoother: d=%d, j=%d, dist=%d ***\n",
                d, j, dist);
    }

    for (t=K->N-1; t>=d && !err; t--) {
        K->t = t;
        if (matrix_is_varying(K, K_T)) {
            err = retrieve_Tt(K);
        }
        if (!err && matrix_is_varying(K, K_ZT)) {
            err = retrieve_Zt(K);
        }
        load_from_row(vt, K->V, t);
        load_from_row(Ft, UF, t);
        load_from_row(Kt, K->K, t);
        load_from_row(at, K->A, t);
        load_from_vech(Pt, K->P, K->r, K->t);

        /* loop across observables */
        for (i=K->n-1; i>=0; i--) {
            vti = vt->val[i];
            Fti = Ft->val[i];
            if (na(vti) || Fti == 0) {
                if (eps_smo) {
                    gretl_matrix_set(veps, t, i, g->val[i]);
                }
                continue;
            }
            ftinv = 1.0 / Fti;
            load_from_row(Zti, ui->Z, i);
            load_from_col(Kti, Kt, i);
            if (eps_smo) {
                eps_smooth_Ft(ftinv, g->val[i], vti, Kti, &c,
                              epshat, veps, t, i, dist);
            }
            fkzero(ftinv, vti, Kti, Zti, &c);
        }

        if (eta_smo) {
            eta_smooth(Q, QRT, &c, etahat, veta, t, dist);
        }

        /* smoothed state and its variance */
        state_smooth(at, Pt, &c, K);

        /* regular backdate for t-1 */
        if (t > 0) {
            fast_copy_values(c.Nt, c.N0);
            if (!K->TI) {
                regular_backdate(&c, K->T);
            }
        }
    } /* end first time loop */

    if (d > 0) {
        K->t = t = d - 1; /* first time-step in diffuse phase */
        if (matrix_is_varying(K, K_T)) {
            err = retrieve_Tt(K);
        }
        if (!err && matrix_is_varying(K, K_ZT)) {
            err = retrieve_Zt(K);
        }
        load_from_row(Kt, K->K, t);
        load_from_row(Ft, UF, t);
        load_from_row(vt, K->V, t);

        /* countdown to observable @j */
        for (i=K->n-1; i>=j; i--) {
            vti = vt->val[i];
            Fti = Ft->val[i];
            if (na(vti) || Fti == 0) {
                if (eps_smo) {
                    gretl_matrix_set(veps, t, i, g->val[i]);
                }
                continue;
            }
            ftinv = 1.0 / Fti;
            load_from_row(Zti, ui->Z, i);
            load_from_col(Kti, Kt, i);
            if (eps_smo) {
                eps_smooth_Ft(ftinv, g->val[i], vti, Kti, &c,
                              epshat, veps, t, i, dist);
            }
            fkzero(ftinv, vti, Kti, Zti, &c);
        }

        /* retrieve the diffuse values */
        load_from_col(Fk, ui->Finf, t);
        load_from_col(Kk, ui->Kinf, t);
        load_from_col(Pk, K->PK, t);

        /* further countdown to first observable */
        for (i=j-1; i>=0; i--) {
            Fki = Fk->val[i];
            Fti = Ft->val[i];
            vti = vt->val[i];
            if (kdebug) {
                fprintf(stderr, "t,i = %d,%d, i-countdown (2) j to 0, "
                        "Fki = %g\n", t, i, Fki);
            }
            if (na(vti)) {
                if (eps_smo) {
                    gretl_matrix_set(veps, t, i, g->val[i]);
                }
                continue;
            }
            if (Fki > k_tiny) { /* was K_TINY */
                fkinv = 1.0 / Fki;
                load_from_row(Zti, ui->Z, i);
                load_from_col(Kki, Kk, i);
                load_from_col(Kti, Kt, i);
                if (eps_smo) {
                    eps_smooth_Fk(fkinv, g->val[i], Kki, &c,
                                  epshat, veps, t, i, dist);
                }
                fkpos(fkinv, Fti, vti, Kki, Kti, Zti, &c);
            } else if (Fti > 0) {
                ftinv = 1.0 / Fti;
                load_from_row(Zti, ui->Z, i);
                load_from_col(Kti, Kt, i);
                if (eps_smo) {
                    eps_smooth_Ft(ftinv, g->val[i], vti, Kti, &c,
                                  epshat, veps, t, i, dist);
                }
                fkzero(ftinv, vti, Kti, Zti, &c);
            }
        }

        if (eta_smo) {
            eta_smooth(Q, QRT, &c, etahat, veta, t, dist);
        }
        if (kdebug) {
            fprintf(stderr, "dagger calc, call 1 (t=%d)\n", t);
        }
        dagger_calc(&c, Pt, Pk, K);

        if (t > 0) {
            if (kdebug) {
                fprintf(stderr, "call diffuse_backdate\n");
            }
            fast_copy_values(c.Nt, c.N0);
            if (!K->TI) {
                diffuse_backdate(&c, K->T);
            }
        }

        for (t=d-2; t>=0 && !err; t--) {
            K->t = t;
            // printf("t,i = %d,%d, t-countdown d-2 to 0\n", t, i);
            if (matrix_is_varying(K, K_T)) {
                err = retrieve_Tt(K);
            }
            if (!err && matrix_is_varying(K, K_ZT)) {
                err = retrieve_Zt(K);
            }
            load_from_row(vt, K->V, t);
            load_from_row(Ft, UF, t);
            load_from_row(Kt, K->K, t);
            load_from_col(Fk, ui->Finf, t);
            load_from_col(Kk, ui->Kinf, t);

            for (i=K->n-1; i>=0; i--) {
                if (kdebug) {
                    fprintf(stderr, "t,i = %d,%d, i-countdown p to 0\n", t, i);
                }
                vti = vt->val[i];
                Fti = Ft->val[i];
                Fki = Fk->val[i];
                if (na(vti)) {
                    if (eps_smo) {
                        gretl_matrix_set(veps, t, i, g->val[i]);
                    }
                    continue;
                }
                if (Fki > k_tiny) {
                    fkinv = 1.0 / Fki;
                    load_from_row(Zti, ui->Z, i);
                    load_from_col(Kki, Kk, i);
                    load_from_col(Kti, Kt, i);
                    if (eps_smo) {
                        eps_smooth_Fk(fkinv, g->val[i], Kki, &c,
                                      epshat, veps, t, i, dist);
                    }
                    fkpos(fkinv, Fti, vti, Kki, Kti, Zti, &c);
                } else if (Fti > 0) {
                    ftinv = 1.0 / Fti;
                    load_from_row(Zti, ui->Z, i);
                    load_from_col(Kti, Kt, i);
                    if (eps_smo) {
                        eps_smooth_Ft(ftinv, g->val[i], vti, Kti, &c,
                                      epshat, veps, t, i, dist);
                    }
                    fkzero(ftinv, vti, Kti, Zti, &c);
                }
            }

            if (kdebug) {
                fprintf(stderr, "dagger calc, call 2 (t=%d)\n", t);
            }
            dagger_calc(&c, Pt, Pk, K);

            if (eta_smo) {
                eta_smooth(Q, QRT, &c, etahat, veta, t, dist);
            }

            if (t > 0) {
                fast_copy_values(c.Nt, c.N0);
                if (!K->TI) {
                    diffuse_backdate(&c, K->T);
                }
            }
        } /* end t=d-2..0 */
    }

    /* add smoothed quantities to @b */
    if (eps_smo) {
        kalman_bundle_add_epshat(K, epshat);
        bundle_add_matrix(K->b, "veps",  veps);
    }
    if (eta_smo) {
        bundle_add_matrix(K->b, "etahat", etahat);
        bundle_add_matrix(K->b, "veta", veta);
    }

    gretl_matrix_block_destroy(B);
    clear_cumulants(&c, K);

    return 0;
}

/* end of univariate smoother functions */

/* Start de Jong filtering and smoothing functions.  These are based
   on the 2003 paper by de Jong and Lin, "Smoothing with an unknown
   initial condition" (Journal of Time Series Analysis 24:2). See also
   de Jong's original paper on this topic, "The diffuse Kalman filter"
   (Annals of Statistics, June 1991).
*/

static int ensure_dj_smo_recorders (kalman *K, int Nd)
{
    if (K->exact && K->djinfo->A == NULL) {
        K->djinfo->A = gretl_zero_matrix_new(K->r * K->r, Nd);
        K->djinfo->V = gretl_zero_matrix_new(K->n * K->r, Nd);
    }

    K->djinfo->Nd = Nd;

    return 0;
}

static int ensure_basic_smo_recorder (kalman *K)
{
    if (K->L == NULL) {
        K->L = gretl_zero_matrix_new(K->r * K->r, K->N);
    }
    if (K->J == NULL && K->smo_prep == SM_DIST) {
        K->J = gretl_zero_matrix_new(K->r * K->r, K->N);
    }

    return 0;
}

static int ensure_dj_variance_matrices (kalman *K)
{
    gretl_matrix *tmp;
    int err = 0;

    if (K->p == 0) {
        K->p = K->r;
        if (K->VY != NULL) {
            K->p += K->n;
        }
    }

    if (K->H == NULL) {
        /* create H from K->VS */
        tmp = gretl_matrix_copy(K->VS);
        if (tmp == NULL) {
            err = E_ALLOC;
        } else {
            err = gretl_matrix_psd_root(tmp, 0);
        }
        if (!err) {
	    int offset = K->VY != NULL ? K->n : 0;

            K->H = gretl_zero_matrix_new(K->r, K->p);
            err = gretl_matrix_inscribe_matrix(K->H, tmp, 0, offset, GRETL_MOD_NONE);
        }
	gretl_matrix_free(tmp);
    }

    if (!err && K->G == NULL && K->VY != NULL) {
        /* create G from K->VY */
        tmp = gretl_matrix_copy(K->VY);
        if (tmp == NULL) {
            err = E_ALLOC;
        } else {
            err = gretl_matrix_psd_root(tmp, 0);
        }
        if (!err) {
            K->G = gretl_zero_matrix_new(K->n, K->p);
            err = gretl_matrix_inscribe_matrix(K->G, tmp, 0, 0, GRETL_MOD_NONE);
        }
	gretl_matrix_free(tmp);
    }

    if (kalman_djvar(K) && !err && K->HG == NULL && K->H != NULL && K->G != NULL) {
        tmp = gretl_matrix_alloc(K->r, K->n);
        if (tmp == NULL) {
            err = E_ALLOC;
        } else {
            gretl_matrix_multiply_mod(K->H, GRETL_MOD_NONE,
                                      K->G, GRETL_MOD_TRANSPOSE,
                                      tmp, GRETL_MOD_NONE);
            if (gretl_is_zero_matrix(tmp)) {
                gretl_matrix_free(tmp);
            } else {
                K->HG = tmp;
            }
        }
    }

    if (!err && K->Jt == NULL) {
        K->Jt = gretl_zero_matrix_new(K->r, K->p);
    }

    if (!err && K->VS0 == NULL && !kalman_djvar(K)) {
	/* establish baseline for detecting changes */
	K->VS0 = gretl_matrix_copy(K->VS);
	if (K->VY != NULL) {
	    K->VY0 = gretl_matrix_copy(K->VY);
	}
    }

    return err;
}

/* The following function assumes that K->H (and K->G, if present) are
   already allocated at appropriate sizes. We come here only if the
   caller has not used the de Jong representation of the disturbance
   variances, in which case BFGS will be tweaking K->VS rather than
   K->H, so that K->H must be updated accordingly. And similarly for
   K->VY and K->G if applicable.
*/

static int update_dj_variance_matrices (kalman *K)
{
    int err = 0;

    if (matrix_changed(K->VS, K->VS0)) {
	/* update H from new K->VS */
	fast_copy_values(K->VS0, K->VS);
	err = gretl_matrix_psd_root(K->VS0, 0);
	if (!err) {
	    int offset = K->VY != NULL ? K->n : 0;

	    err = gretl_matrix_inscribe_matrix(K->H, K->VS0, 0, offset, GRETL_MOD_NONE);
	}
	fast_copy_values(K->VS0, K->VS);
    }

    if (err || K->VY == NULL) {
	return err;
    }

    if (matrix_changed(K->VY, K->VY0)) {
        /* update G from new K->VY */
	fast_copy_values(K->VY0, K->VY);
	err = gretl_matrix_psd_root(K->VY0, 0);
        if (!err) {
            err = gretl_matrix_inscribe_matrix(K->G, K->VY0, 0, 0, GRETL_MOD_NONE);
        }
	fast_copy_values(K->VY0, K->VY);
    }

    return err;
}

static int dejong_diffuse_filter_step (kalman *K)
{
    dejong_info *dj = K->djinfo;
    gretl_matrix *lam;
    double lmin = 0;
    int err = 0;

    /* Vt = -Zt * At */
    gretl_matrix_zero(dj->Vt);
    gretl_matrix_multiply_mod(K->ZT, GRETL_MOD_TRANSPOSE,
                              dj->At, GRETL_MOD_NONE,
                              dj->Vt, GRETL_MOD_DECREMENT);
    if (dj->V != NULL) {
        record_to_col(dj->V, dj->Vt, K->t);
    }

    /* Vv = Vt ~ vt */
    gretl_matrix_inscribe_matrix(dj->Vv, dj->Vt,
                                 0, 0, GRETL_MOD_NONE);
    gretl_matrix_inscribe_matrix(dj->Vv, K->vt,
                                 0, K->r, GRETL_MOD_NONE);

    /* Q[t+1] : Qt += Vv' * iFt * Vv */
    gretl_matrix_qform(dj->Vv, GRETL_MOD_TRANSPOSE,
                       K->iFt, dj->Qt, GRETL_MOD_CUMULATE);
    if (kdebug > 1) {
	gretl_matrix_print(dj->Qt, "Qt");
    }

    /* A[t+1]: At = T*At + Kt*Vt */
    gretl_matrix_multiply(K->Kt, dj->Vt, K->rr);
    if (K->TI) {
        gretl_matrix_add_to(dj->At, K->rr);
    } else {
        gretl_matrix_multiply_mod(K->T, GRETL_MOD_NONE,
                                  dj->At, GRETL_MOD_NONE,
                                  K->rr, GRETL_MOD_CUMULATE);
        fast_copy_values(dj->At, K->rr);
    }

    /* Sinv = Qt[1:r,1:r] */
    gretl_matrix_extract_matrix(dj->S, dj->Qt, 0, 0, GRETL_MOD_NONE);
    if (K->r == 1) {
        lmin = dj->S->val[0];
    } else {
        lam = gretl_symmetric_matrix_eigenvals(dj->S, 0, &err);
        lmin = lam->val[0];
        gretl_matrix_free(lam);
    }
    if (kdebug > 1) {
	fprintf(stderr, "dejong_diffuse_filter_step: t=%d, lmin=%g\n", K->t, lmin);
    }

    if (na(lmin)) {
	err = E_NAN;
    } else if (lmin > 1.0e-7) {
        /* handle the collapse */
        K->d = K->t + 1;
	if (kdebug) {
	    fprintf(stderr, "*** determined m = %d ***\n", K->d);
	}

        /* Sinv = Qt[1:r,1:r] */
        gretl_matrix_extract_matrix(dj->S, dj->Qt, 0, 0, GRETL_MOD_NONE);

        /* -s = Qt[1:r,r+1] */
        gretl_matrix_extract_matrix(dj->s, dj->Qt, 0, K->r, GRETL_MOD_NONE);
        gretl_matrix_multiply_by_scalar(dj->s, -1.0);

        /* S = invpd(Sinv) */
        gretl_invert_symmetric_matrix2(dj->S, &dj->ldS);

        /* at += At*S*s */
        gretl_matrix_multiply(dj->At, dj->S, K->rr);
        gretl_matrix_multiply_mod(K->rr, GRETL_MOD_NONE,
                                  dj->s, GRETL_MOD_NONE,
                                  K->a0, GRETL_MOD_CUMULATE);

        /* Pt += At*S*At' */
        gretl_matrix_qform(dj->At, GRETL_MOD_NONE,
                           dj->S, K->P0, GRETL_MOD_CUMULATE);

        /* qm = Qt[r+1,r+1] - s'*S*s */
        dj->qm = gretl_matrix_get(dj->Qt, K->r, K->r);
        dj->qm -= gretl_scalar_qform(dj->s, dj->S, &err); /* FIXME? */

        /* record */
	if (dj->Am == NULL) {
	    dj->Am = gretl_matrix_copy(dj->At);
	} else {
	    fast_copy_values(dj->Am, dj->At);
	}
    }

    return err;
}

/* Carry out some calculations that are required when the
   diffuse initial recursion never collapsed to the regular
   Kalman filter. @sum_preldet is the sum of log determinants
   of F_t, t = 1 to N, all of which observations were in the
   diffuse phase.
*/

static int handle_no_collapse_case (kalman *K, double nl2pi,
                                    double sum_preldet)
{
    dejong_info *dj = K->djinfo;
    double cterm = K->okN * nl2pi;
    double qm, qadj;
    int err = 0;

    qm = gretl_matrix_get(dj->Qt, K->r, K->r);
    gretl_matrix_extract_matrix(dj->s, dj->Qt, 0, K->r, GRETL_MOD_NONE);
    gretl_matrix_multiply_by_scalar(dj->s, -1.0);
    err = gretl_matrix_moore_penrose(dj->S, 1.0e-7);
    if (!err) {
        qadj = gretl_scalar_qform(dj->s, dj->S, &err);
    }
    if (!err) {
        K->loglik = -0.5 * (sum_preldet + qm - qadj + cterm);
    }
    if (dj->Am != NULL) {
        gretl_matrix_free(dj->Am);
	dj->Am = NULL;
    }

    return err;
}

static int kfilter_dejong (kalman *K, PRN *prn)
{
    double nl2pi = K->n * LN_2_PI;
    double sumldet = 0;
    double sum_preldet = 0;
    double llt, llct, ldt, qt, pldt;
    int exact_diffuse;
    int loglik_1991 = 0;
    int nt, Nd = K->N;
    int err = 0;

    exact_diffuse = kalman_diffuse(K) && K->exact;

    if (identify || trace) {
	printf("kfilter_dejong() exact_diffuse = %d\n", exact_diffuse);
    }

    if (K->dj_initted == 0) {
	err = ensure_dj_variance_matrices(K);
	if (err) {
	    return err;
	}
	K->dj_initted = 1;
    }

    if (K->smo_prep) {
        ensure_basic_smo_recorder(K);
        if (K->exact) {
            Nd = K->N; // ? matrix_is_varying(K, K_ZT) ? K->N : 4*K->r;
            ensure_dj_smo_recorders(K, Nd);
        }
    }

    if (exact_diffuse) {
	/* (re-)initialize djinfo */
        fast_write_I(K->djinfo->At);
	gretl_matrix_zero(K->djinfo->Qt);
	gretl_matrix_zero(K->djinfo->Vt);
    }

    if (kdebug) {
        fprintf(stderr, "\n*** kfilter_dejong: diffuse %d, exact %d, smo %d, djvar %d\n",
                kalman_diffuse(K) ? 1 : 0, K->exact, K->smo_prep ? 1 : 0, kalman_djvar(K));
    }

    if (exact_diffuse) {
        K->d = K->N;
    } else {
        K->d = -1;
    }
    K->SSRw = K->loglik = 0.0;
    K->s2 = NADBL;
    K->TI = gretl_is_identity_matrix(K->T);
    K->okN = K->N;

    if (!kalman_djvar(K) && (gretl_iteration_depth() > 0 || numgrad_in_progress())) {
	/* guard against "back door" changes in dist. variance parameters under mle */
        err = update_dj_variance_matrices(K);
    }

    set_kalman_running(K);

    for (K->t = 0; K->t < K->N && !err; K->t += 1) {
        llct = ldt = qt = 0;
        llt = NADBL;

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

        if (kdebug > 1) {
            kalman_print_state(K);
        }

        /* calculate v_t, checking for missing values */
        nt = compute_forecast_error(K);
        if (nt == 0) {
            /* skip this observation */
            K->okN -= 1;
            handle_missing_obs(K);
	    if (kdebug) {
		fprintf(stderr, "dejong filter: NA at t=%d\n", K->t);
	    }
            continue;
        }
        /* and record forecast errors if wanted */
        if (K->V != NULL) {
            record_to_row(K->V, K->vt, K->t);
        }

        /* calculate F_t = ZPZ' [+ VY] (where VY = GG') */
        if (K->VY != NULL) {
            fast_copy_values(K->Ft, K->VY);
            gretl_matrix_qform(K->ZT, GRETL_MOD_TRANSPOSE, K->P0,
                               K->Ft, GRETL_MOD_CUMULATE);
        } else {
            gretl_matrix_qform(K->ZT, GRETL_MOD_TRANSPOSE, K->P0,
                               K->Ft, GRETL_MOD_NONE);
        }

        /* calculate M_t = TPZ' [+ HG'] */
        gretl_matrix_multiply(K->P0, K->ZT, K->PZ);
        if (K->TI) {
            fast_copy_values(K->Mt, K->PZ);
        } else {
            gretl_matrix_multiply(K->T, K->PZ, K->Mt);
        }
        if (K->HG != NULL) {
            gretl_matrix_add_to(K->Mt, K->HG);
        }

        /* calculate F_t and its inverse */
        fast_copy_values(K->iFt, K->Ft);
	if (K->t < K->d && !loglik_1991) {
	    /* conditionality regarding the log-determinant? */
	    err = gretl_invert_symmetric_matrix2(K->iFt, &pldt);
            sum_preldet += pldt;
	} else {
	    err = gretl_invert_symmetric_matrix2(K->iFt, &ldt);
	}
        if (err) {
            fprintf(stderr, "kfilter_dejong: failed to invert Ft\n");
	    break;
        }

        if (K->F != NULL) {
            /* we're recording F_t for all t */
            if (K->smo_prep) {
                /* record inverse */
                record_to_vech(K->F, K->iFt, nt, K->t);
            } else {
                /* record F_t itself */
                record_to_vech(K->F, K->Ft, nt, K->t);
            }
        }

        /* gain: K_t = M_t F_t^{-1} */
        gretl_matrix_multiply(K->Mt, K->iFt, K->Kt);
        /* L_t = T_t - K_t Z_t */
        fast_copy_values(K->Lt, K->T);
        gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
                                  K->ZT, GRETL_MOD_TRANSPOSE,
                                  K->Lt, GRETL_MOD_DECREMENT);

        /* J_t = H_t - K_t * G_t */
        fast_copy_values(K->Jt, K->H);
        if (K->G != NULL) {
            gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
                                      K->G,  GRETL_MOD_NONE,
                                      K->Jt, GRETL_MOD_DECREMENT);
        }

        /* record matrices for the smoother */
        if (K->L != NULL) {
            record_to_col(K->L, K->Lt, K->t);
        }
        if (K->J != NULL) {
            record_to_col(K->J, K->Jt, K->t);
        }
        if (K->K != NULL) {
            record_to_vec(K->K, K->Kt, K->t);
        }

        /* update state: a1 = T a0 + K_t v_t [+ mu] */
        gretl_matrix_multiply(K->T, K->a0, K->a1);
        if (K->mu != NULL) {
            gretl_matrix_add_to(K->a1, K->mu);
        }
        gretl_matrix_multiply_mod(K->Kt, GRETL_MOD_NONE,
                                  K->vt, GRETL_MOD_NONE,
                                  K->a1, GRETL_MOD_CUMULATE);
        fast_copy_values(K->a0, K->a1);

        /* update var(state): P1 = TPL' + HJ' */
        gretl_matrix_multiply_mod(K->H, GRETL_MOD_NONE,
                                  K->Jt, GRETL_MOD_TRANSPOSE,
                                  K->P1, GRETL_MOD_NONE);
        gretl_matrix_multiply_mod(K->P0, GRETL_MOD_NONE,
                                  K->Lt, GRETL_MOD_TRANSPOSE,
                                  K->rr, GRETL_MOD_NONE);
        if (K->TI) {
            gretl_matrix_add_to(K->P1, K->rr);
        } else {
            gretl_matrix_multiply_mod(K->T,  GRETL_MOD_NONE,
                                      K->rr, GRETL_MOD_NONE,
                                      K->P1, GRETL_MOD_CUMULATE);
        }
        fast_copy_values(K->P0, K->P1);

        if (exact_diffuse && K->t < K->d) {
	    /* still in diffuse phase */
            err = dejong_diffuse_filter_step(K);
        } else {
            /* standard Kalman phase */
            qt = gretl_scalar_qform(K->vt, K->iFt, &err);
            llct = nl2pi;
        }

        /* determine and record loglikelihood */
        if (err) {
            K->loglik = NADBL;
            break;
        } else {
            llt = -0.5 * (llct + ldt + qt);
            if (na(llt)) {
                K->loglik = NADBL;
                break;
            }
            K->loglik += llt;
            K->SSRw += qt;
            sumldet += ldt;
        }
        if (K->LL != NULL) {
            gretl_vector_set(K->LL, K->t, llt);
        }
    }

    set_kalman_stopped(K);

    if (na(K->loglik)) {
        err = E_NAN;
    } else if (exact_diffuse) {
        if (K->d == K->N) {
            err = handle_no_collapse_case(K, nl2pi, sum_preldet);
        } else if (loglik_1991) {
            /* deal with extra terms set out in De Jong (1991) */
            double ll_adj = K->djinfo->qm + K->djinfo->ldS;

            K->loglik -= 0.5 * ll_adj;
        }
    } else if (kalman_arma_ll(K)) {
        double ll1 = 1.0 + LN_2_PI + log(K->SSRw / K->okN);

        K->loglik = -0.5 * (K->okN * ll1 + sumldet);
    } else {
        int d = (kalman_diffuse(K) && !K->exact)? K->r : 0;

        /* as per old kalman_forecast() code:
           do we want this at all?
        */
        K->s2 = K->SSRw / (K->n * K->okN - d);
    }

    if (kdebug) {
        fprintf(stderr, "kfilter_dejong: err=%d, ll=%#.8g, d=%d\n",
                err, K->loglik, K->d);
    }

    return err;
}

/* r_{t-1} = Z_t' F_t^{-1} v_t + L_t' r_t */

static void rt_recursion (kalman *K,
                          gretl_matrix *n1,
                          gretl_matrix *r0,
                          gretl_matrix *r1)
{
    gretl_matrix_multiply(K->iFt, K->vt, n1);
    if (K->t == K->N - 1) {
        gretl_matrix_multiply(K->ZT, n1, r0);
    } else {
        gretl_matrix_multiply(K->ZT, n1, r1);
        gretl_matrix_multiply_mod(K->Lt, GRETL_MOD_TRANSPOSE,
                                  r0, GRETL_MOD_NONE,
                                  r1, GRETL_MOD_CUMULATE);
        fast_copy_values(r0, r1);
    }
}

/* N_{t-1} = Z_t' F_t^{-1} Z_t + L_t' N_t L_t */

static void Nt_recursion (kalman *K,
                          gretl_matrix *N0,
                          gretl_matrix *N1)
{
    if (K->t == K->N - 1) {
        gretl_matrix_qform(K->ZT, GRETL_MOD_NONE,
                           K->iFt, N0, GRETL_MOD_NONE);
    } else {
        gretl_matrix_qform(K->ZT, GRETL_MOD_NONE,
                           K->iFt, N1, GRETL_MOD_NONE);
        load_from_col(K->Lt, K->L, K->t);
        gretl_matrix_qform(K->Lt, GRETL_MOD_TRANSPOSE,
                           N0, N1, GRETL_MOD_CUMULATE);
        fast_copy_values(N0, N1);
    }
}

static int record_dhat_and_Psi (kalman *K,
                                gretl_matrix *dhat,
                                gretl_matrix *Psi,
                                gretl_matrix *rt,
                                gretl_matrix *Nt)
{
    dejong_info *dj = K->djinfo;
    gretl_matrix *tmp;
    int err = 0;

    /* dhat = S * (s + Am' * rt) */
    tmp = gretl_matrix_copy(dj->s);
    gretl_matrix_multiply_mod(dj->Am, GRETL_MOD_TRANSPOSE,
                              rt, GRETL_MOD_NONE,
                              tmp, GRETL_MOD_CUMULATE);
    gretl_matrix_multiply(dj->S, tmp, dhat);
    gretl_matrix_free(tmp);

    /* Psi = S - S * Am' * Nt * Am * S */
    fast_copy_values(Psi, dj->S);
    gretl_matrix_qform(dj->Am, GRETL_MOD_TRANSPOSE,
                       Nt, K->rr, GRETL_MOD_NONE);
    gretl_matrix_qform(dj->S, GRETL_MOD_NONE, K->rr,
                       Psi, GRETL_MOD_DECREMENT);

    return err;
}

static int dejong_diffuse_smoother_step (kalman *K,
                                         gretl_matrix *rt,
                                         gretl_matrix *Rt,
                                         gretl_matrix *Bt,
                                         gretl_matrix *Nt,
                                         gretl_matrix *dhat,
                                         gretl_matrix *Psi,
                                         int no_collapse)
{
    dejong_info *dj = K->djinfo;
    gretl_matrix *tmp1 = NULL;
    gretl_matrix *tmp2 = NULL;
    gretl_matrix *Rr = NULL;
    int err = 0;

    load_from_col(dj->Vt, dj->V, K->t);
    load_from_col(dj->At, dj->A, K->t);
    load_from_col(K->Lt, K->L, K->t);

    tmp1 = gretl_matrix_alloc(K->r, K->r + 1);
    tmp2 = gretl_matrix_alloc(K->n, K->r + 1);

    /* Rr = Rt ~ rt */
    Rr = gretl_matrix_alloc(K->r, K->r + 1);
    gretl_matrix_inscribe_matrix(Rr, Rt, 0, 0, GRETL_MOD_NONE);
    gretl_matrix_inscribe_matrix(Rr, rt, 0, K->r, GRETL_MOD_NONE);

    /* Vv = Vt ~ vt */
    gretl_matrix_inscribe_matrix(dj->Vv, dj->Vt, 0, 0, GRETL_MOD_NONE);
    gretl_matrix_inscribe_matrix(dj->Vv, K->vt, 0, K->r, GRETL_MOD_NONE);

    /* Rr = Z'*iFt*Vv + Lt'*Rr */
    gretl_matrix_multiply_mod(K->Lt, GRETL_MOD_TRANSPOSE,
                              Rr, GRETL_MOD_NONE,
                              tmp1, GRETL_MOD_NONE);
    fast_copy_values(Rr, tmp1);
    gretl_matrix_multiply(K->iFt, dj->Vv, tmp2);
    gretl_matrix_multiply_mod(K->ZT, GRETL_MOD_NONE,
                              tmp2, GRETL_MOD_NONE,
                              Rr, GRETL_MOD_CUMULATE);

    /* Rt = Rr[,1:r] ; rt = Rr[,r+1] */
    gretl_matrix_extract_matrix(Rt, Rr, 0, 0, GRETL_MOD_NONE);
    gretl_matrix_extract_matrix(rt, Rr, 0, K->r, GRETL_MOD_NONE);

    /* Bt = At + Pt * Rt */
    fast_copy_values(Bt, dj->At);
    gretl_matrix_multiply_mod(K->P0, GRETL_MOD_NONE,
                              Rt, GRETL_MOD_NONE,
                              Bt, GRETL_MOD_CUMULATE);

    /* aht = at + Pt * rt */
    gretl_matrix_multiply_mod(K->P0, GRETL_MOD_NONE,
                              rt, GRETL_MOD_NONE,
                              K->a0, GRETL_MOD_CUMULATE);

    /* aht += Bt * dhat */
    gretl_matrix_multiply_mod(Bt, GRETL_MOD_NONE,
                              dhat, GRETL_MOD_NONE,
                              K->a0, GRETL_MOD_CUMULATE);

    /* vht = Pt - Pt * Nt * Pt */
    fast_copy_values(K->P1, K->P0);
    gretl_matrix_qform(K->P0, GRETL_MOD_TRANSPOSE,
                       Nt, K->P1, GRETL_MOD_DECREMENT);

    if (no_collapse) {
        /* vht += Bt * Psi * Bt' */
        gretl_matrix_qform(Bt, GRETL_MOD_NONE,
                           Psi, K->P1, GRETL_MOD_CUMULATE);
    } else {
        gretl_matrix *BSC = gretl_matrix_reuse(tmp1, K->r, K->r);

        /* Ct = Pt * (Rt + Nt * At) */
        gretl_matrix_multiply(Nt, dj->At, K->rr);
        gretl_matrix_add_to(K->rr, Rt);
        gretl_matrix_multiply(K->P0, K->rr, K->Ct);

        /* BSC = Bt * S * Ct' */
        gretl_matrix_multiply_mod(dj->S, GRETL_MOD_NONE,
                                  K->Ct, GRETL_MOD_TRANSPOSE,
                                  K->rr, GRETL_MOD_NONE);
        gretl_matrix_multiply(Bt, K->rr, BSC);

        /* vht += Bt * Psi * Bt' - BSC - BSC' */
        gretl_matrix_qform(Bt, GRETL_MOD_NONE,
                           Psi, K->rr, GRETL_MOD_NONE);
        gretl_matrix_subtract_from(K->rr, BSC);
        gretl_square_matrix_transpose(BSC);
        gretl_matrix_subtract_from(K->rr, BSC);
        gretl_matrix_add_to(K->P1, K->rr);
    }

    record_to_row(K->A, K->a0, K->t);
    record_to_vech(K->P, K->P1, K->r, K->t);

    gretl_matrix_free(tmp1);
    gretl_matrix_free(tmp2);
    gretl_matrix_free(Rr); /* FIXME make reusable? */

    return err;
}

static int state_smooth_dejong (kalman *K)
{
    gretl_matrix_block *B, *B2 = NULL;
    gretl_matrix *r0, *r1, *N0, *N1, *n1;
    gretl_matrix *Bt, *Rt;
    gretl_matrix *dhat, *Psi;
    int no_collapse = 0;
    int nt = K->n;
    int t, err = 0;

    no_collapse = (K->d == K->N);

    if (trace) {
        printf("state_smooth_dejong()\n");
    }
    if (kdebug) {
        fprintf(stderr, "state_smooth_dejong: no_collapse = %d\n", no_collapse);
    }

    B = gretl_matrix_block_new(&r0,  K->r, 1,
                               &r1,  K->r, 1,
                               &N0,  K->r, K->r,
                               &N1,  K->r, K->r,
                               &n1,  K->n, 1,
                               NULL);
    if (B == NULL) {
        err = E_ALLOC;
    }

    if (!err && K->exact) {
        B2 = gretl_matrix_block_new(&Bt,   K->r, K->r,
                                    &Rt,   K->r, K->r,
                                    &dhat, K->r, 1,
                                    &Psi,  K->r, K->r,
                                    NULL);
        if (B2 == NULL) {
            err = E_ALLOC;
        }
    }

    if (err) {
        return err;
    }

    gretl_matrix_zero(r0);
    gretl_matrix_zero(N0);
    if (B2 != NULL) {
        gretl_matrix_zero(Rt);
    }

    if (no_collapse) {
        gretl_matrix_multiply(K->djinfo->S, K->djinfo->s, dhat);
        gretl_matrix_copy_values(Psi, K->djinfo->S);
    }

    for (t=K->N-1; t>=0 && !err; t--) {
        K->t = t;
	nt = load_filter_data(K, no_collapse, &err);
        if (err) {
            break;
        } else if (kdebug && nt == 0) {
            fprintf(stderr, "state_smooth_dejong: nt=0 at t=%d\n", K->t);
        }

        /* compute N_{t-1} */
	Nt_recursion(K, N0, N1);

        if (t >= K->d) {
            /* regular recursion: compute r_{t-1} */
	    rt_recursion(K, n1, r0, r1);
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
            if (t == K->d) {
                record_dhat_and_Psi(K, dhat, Psi, r1, N1);
            }
        } else {
            /* t < K->d: diffuse phase (FIXME nt == 0?) */
            dejong_diffuse_smoother_step(K, r1, Rt, Bt, N1,
                                         dhat, Psi, no_collapse);
        }
    }

    gretl_matrix_block_destroy(B2);
    gretl_matrix_block_destroy(B);

    return err;
}

static int diffuse_dist_smooth_step (kalman *K,
                                     gretl_matrix *Rt,
                                     gretl_matrix *Bt,
                                     gretl_matrix *Ct,
                                     gretl_matrix *Nt,
                                     gretl_matrix *dhat,
                                     gretl_matrix *Psi,
                                     gretl_matrix *nr,
                                     gretl_matrix *nu_t,
                                     gretl_matrix *Vnu_t,
                                     int DKstyle,
                                     int no_collapse)
{
    dejong_info *dj = K->djinfo;
    gretl_matrix *BSC = NULL;
    gretl_matrix *pp = NULL;
    gretl_matrix *rp = NULL;
    int err = 0;

    if (!no_collapse) {
        BSC = gretl_matrix_alloc(K->p, K->p);
        rp = gretl_matrix_alloc(K->r, K->p);
    }
    pp = gretl_matrix_alloc(K->p, K->p);

    /* Bt = G' * iFt * Vt + Jt' * Rt */
    gretl_matrix_multiply_mod(K->Jt, GRETL_MOD_TRANSPOSE,
                              Rt, GRETL_MOD_NONE,
                              Bt, GRETL_MOD_NONE);
    if (K->G != NULL) {
        gretl_matrix_multiply(K->iFt, dj->Vt, nr);
        gretl_matrix_multiply_mod(K->G, GRETL_MOD_TRANSPOSE,
                                  nr, GRETL_MOD_NONE,
                                  Bt, GRETL_MOD_CUMULATE);
    }

    if (!no_collapse) {
        /* Ct = Jt' * (Rt + Nt * A_{t+1}) */
        gretl_matrix_multiply(Nt, dj->At, K->rr);
        gretl_matrix_add_to(K->rr, Rt);
        gretl_matrix_multiply_mod(K->Jt, GRETL_MOD_TRANSPOSE,
                                  K->rr, GRETL_MOD_NONE,
                                  Ct, GRETL_MOD_NONE);
    }

    /* nu_t += Bt * dhat */
    gretl_matrix_multiply_mod(Bt, GRETL_MOD_NONE,
                              dhat, GRETL_MOD_NONE,
                              nu_t, GRETL_MOD_CUMULATE);

    if (!no_collapse) {
        /* BSC = Bt * S * Ct' */
        gretl_matrix_multiply_mod(dj->S, GRETL_MOD_NONE,
                                  Ct, GRETL_MOD_TRANSPOSE,
                                  rp, GRETL_MOD_NONE);
        gretl_matrix_multiply(Bt, rp, BSC);
    }

    /* Vnu_adj = Bt * Psi * Bt' - BSC - BSC' */
    gretl_matrix_qform(Bt, GRETL_MOD_NONE, Psi,
                       pp, GRETL_MOD_NONE);
    if (!no_collapse) {
        gretl_matrix_subtract_from(pp, BSC);
        gretl_square_matrix_transpose(BSC);
        gretl_matrix_subtract_from(pp, BSC);
    }

    if (DKstyle) {
        gretl_matrix_add_to(Vnu_t, pp);
    } else {
        gretl_matrix_subtract_from(Vnu_t, pp);
    }

    if (!no_collapse) {
        gretl_matrix_free(BSC);
        gretl_matrix_free(rp);
    }
    gretl_matrix_free(pp);

    return err;
}

static int dist_smooth_dejong (kalman *K, int DKstyle)
{
    gretl_matrix_block *B, *B2 = NULL;
    gretl_matrix *nu_t, *vnu_t;
    gretl_matrix *r0, *r1, *N0, *N1, *n1;
    gretl_matrix *Bt, *Ct, *Rt, *nr, *nn;
    gretl_matrix *dhat, *Psi;
    gretl_matrix *nu = NULL, *vnu = NULL;
    gretl_matrix *H_nu, *G_nu = NULL;
    gretl_matrix *etahat, *veta;
    gretl_matrix *epshat = NULL;
    gretl_matrix *veps = NULL;
    dejong_info *dj = K->djinfo;
    int no_collapse = 0;
    int nt = K->n;
    int t, err = 0;

    if (trace) {
        printf("dist_smooth_dejong(), DKstyle = %d\n", DKstyle);
    }

    no_collapse = (K->d == K->N);

    if (1 /* K->HG != NULL ? */) {
        nu  = gretl_zero_matrix_new(K->p, K->N);
        vnu = gretl_zero_matrix_new(K->p * K->p, K->N);
    }

    etahat = gretl_matrix_alloc(K->N, K->r);
    veta = gretl_matrix_alloc(K->N, K->r * K->r);
    if (K->G != NULL) {
        epshat = gretl_matrix_alloc(K->N, K->n);
        veps = gretl_matrix_alloc(K->N, K->n * K->n);
        G_nu = gretl_matrix_alloc(K->n, 1);
    }

    B = gretl_matrix_block_new(&r0,  K->r, 1,
                               &r1,  K->r, 1,
                               &N0,  K->r, K->r,
                               &N1,  K->r, K->r,
                               &n1,  K->n, 1,
                               &nn,  K->n, K->n,
                               &nr,  K->n, K->r,
                               &nu_t,  K->p, 1,
                               &vnu_t, K->p, K->p,
                               &H_nu,  K->r, 1,
                               NULL);
    if (B == NULL) {
        err = E_ALLOC;
    }

    if (!err && K->exact) {
        B2 = gretl_matrix_block_new(&Bt,   K->p, K->r,
                                    &Ct,   K->p, K->r,
                                    &Rt,   K->r, K->r,
                                    &dhat, K->r, 1,
                                    &Psi,  K->r, K->r,
                                    NULL);
        if (B2 == NULL) {
            err = E_ALLOC;
        }
    }

    if (err) {
        return err;
    }

    gretl_matrix_zero(r0);
    gretl_matrix_zero(N0);
    if (B2 != NULL) {
        gretl_matrix_zero(Rt);
    }

    for (t=K->N-1; t>=0 && !err; t--) {
        K->t = t;
        nt = load_filter_data(K, 0, &err);
        if (err) {
            break;
        } else if (kdebug && nt == 0) {
	    fprintf(stderr, "dist_smooth_dejong: nt=0 at t=%d\n", K->t);
        }

	load_from_col(K->Jt, K->J, K->t);
	load_from_col(K->Lt, K->L, K->t);

	/* nu_t = G' * iFt * vt + Jt' * rt */
	gretl_matrix_multiply_mod(K->Jt, GRETL_MOD_TRANSPOSE,
				  r0, GRETL_MOD_NONE,
				  nu_t, GRETL_MOD_NONE);
	if (K->G != NULL) {
	    gretl_matrix_multiply(K->iFt, K->vt, n1);
	    gretl_matrix_multiply_mod(K->G, GRETL_MOD_TRANSPOSE,
				      n1, GRETL_MOD_NONE,
				      nu_t, GRETL_MOD_CUMULATE);
	}

        if (DKstyle) {
            /* vnu = I(r) - (G' * iFt * G + Jt' * Nt * Jt) */
            gretl_matrix_inscribe_I(vnu_t, 0, 0, K->p);
            gretl_matrix_qform(K->Jt, GRETL_MOD_TRANSPOSE,
                               N0, vnu_t, GRETL_MOD_DECREMENT); /* N1? */
            if (K->G != NULL) {
                gretl_matrix_qform(K->G, GRETL_MOD_TRANSPOSE,
                                   K->iFt, vnu_t, GRETL_MOD_DECREMENT);
            }
        } else {
	    /* vnu = G' * iFt * G + Jt' * Nt * Jt */
	    gretl_matrix_qform(K->Jt, GRETL_MOD_TRANSPOSE, N0,
			       vnu_t, GRETL_MOD_NONE); /* N1? */
	    if (K->G != NULL) {
		gretl_matrix_qform(K->G, GRETL_MOD_TRANSPOSE, K->iFt,
				   vnu_t, GRETL_MOD_CUMULATE);
	    }
        }

        if (t == K->d) {
            /* rm1 = Z' * iFt * vt + Lt' * rt */
            rt_recursion(K, n1, r0, r1);
            /* Nm1 = Z' * iFt * Z + Lt' * Nt * Lt */
            Nt_recursion(K, N0, N1);
            record_dhat_and_Psi(K, dhat, Psi, r1, N1); /* r0? N0? */
        } else if (t < K->d) {
            load_from_col(dj->Vt, dj->V, t);
            if (t == K->d - 1) {
                fast_copy_values(dj->At, dj->Am);
            } else {
                load_from_col(dj->At, dj->A, t);
            }
            diffuse_dist_smooth_step(K, Rt, Bt, Ct, N1, dhat, Psi,
                                     nr, nu_t, vnu_t, DKstyle,
                                     no_collapse);
        }

        /* record t-dated results */
        if (nu != NULL) {
            record_to_col(nu, nu_t, t);
            record_to_col(vnu, vnu_t, t);
        }
        gretl_matrix_multiply(K->H, nu_t, H_nu);
        record_to_row(etahat, H_nu, t);
        gretl_matrix_qform(K->H, GRETL_MOD_NONE, vnu_t,
                           K->rr, GRETL_MOD_NONE);
        record_to_vec(veta, K->rr, t);
        if (K->G != NULL) {
            gretl_matrix_multiply(K->G, nu_t, G_nu);
            record_to_row(epshat, G_nu, t);
	    if (nt == 0) {
		record_to_vec(veps, K->VY, t);
	    } else {
		gretl_matrix_qform(K->G, GRETL_MOD_NONE, vnu_t,
				   nn, GRETL_MOD_NONE);
		record_to_vec(veps, nn, t);
	    }
        }

        if (t > 1) {
            /* t-1 dated quantities */
            if (t == K->d) {
                fast_copy_values(r0, r1); /* or vice versa? */
                fast_copy_values(N0, N1); /* ditto */
            } else {
                rt_recursion(K, n1, r0, r1);
                Nt_recursion(K, N0, N1);
            }
            if (t < K->d) {
		/* Rt = Z' * iFt * Vt + Lt' * Rt */
		gretl_matrix_multiply_mod(K->Lt, GRETL_MOD_TRANSPOSE,
					  Rt, GRETL_MOD_NONE,
					  K->rr, GRETL_MOD_NONE);
		gretl_matrix_multiply(K->iFt, dj->Vt, nr);
		gretl_matrix_multiply_mod(K->ZT, GRETL_MOD_NONE,
					  nr, GRETL_MOD_NONE,
					  K->rr, GRETL_MOD_CUMULATE);
                fast_copy_values(Rt, K->rr);
            }
        }
    }

    /* add results to K->b */
    bundle_add_matrix(K->b, "etahat", etahat);
    bundle_add_matrix(K->b, "veta", veta);
    if (epshat != NULL) {
        bundle_add_matrix(K->b, "epshat", epshat);
        bundle_add_matrix(K->b, "veps",  veps);
    }
    if (nu != NULL) {
        bundle_add_matrix(K->b, "nuhat", nu);
        bundle_add_matrix(K->b, "vnu", vnu);
    }

    gretl_matrix_free(G_nu);
    gretl_matrix_block_destroy(B);
    gretl_matrix_block_destroy(B2);

    return err;
}

#if SUPPORT_LEGACY

/* Below: append legacy smoothing functions. These can be removed
   if/when we're confident that the new code is fully OK.
*/

#include "kalman_legacy.c"

#endif
