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

#include "libgretl.h"
#include "version.h"
#include "pvalues.h"
#include "gretl_matrix.h"
#include "matrix_extra.h"
#include "var.h"
#include "johansen.h"
#include "vartest.h"
#include "varprint.h"
#include "libset.h"
#include "jprivate.h"

#define JDEBUG 0

/* coefficient matrices for the trace test */

const double trace_m_coef[5][6] = {
    /* N^2  N       1   N==1   N==2  N^1/2 */
    { 2, -1.00,  0.07,  0.07,     0,     0 },
    { 2,  2.01,  0.00,  0.06,  0.05,     0 },
    { 2,  1.05, -1.55, -0.50, -0.23,     0 },
    { 2,  4.05,  0.50, -0.23, -0.07,     0 },
    { 2,  2.85, -5.10, -0.10, -0.06,  1.35 }
};

const double trace_v_coef[5][6] = {
    { 3, -0.33, -0.55,  0.0,  0.00,  0 },
    { 3,  3.60,  0.75, -0.4, -0.30,  0 },
    { 3,  1.80,  0.00, -2.8, -1.10,  0 },
    { 3,  5.70,  3.20, -1.3, -0.50,  0 },
    { 3,  4.00,  0.80, -5.8, -2.66,  0 }
};

/* trace-test sample size correction: log(m)-log(m-hat) */

const double trace_m_corr[5][7] = {
    /* sqrt(n)/T n/T  n^2/T^2  n==1/T     n==1     n==2      n==3 */
    { -0.101,  0.499,  0.896, -0.562,  0.00229, 0.00662,        0 },
    {      0,  0.465,  0.984, -0.273, -0.00244,       0,        0 },
    {  0.134,  0.422,   1.02,   2.17, -0.00182,       0, -0.00321 },
    { 0.0252,  0.448,   1.09, -0.353,        0,       0,        0 },
    { -0.819,  0.615,  0.896,   2.43,  0.00149,       0,        0 }
};

/* trace-test sample size correction: log(v)-log(v-hat) */

const double trace_v_corr[5][7] = {
    { -0.204,  0.980,  3.11,  -2.14,   0.0499,  -0.0103, -0.00902 },
    {  0.224,  0.863,  3.38, -0.807,        0,        0,  -0.0091 },
    {  0.422,  0.734,  3.76,   4.32, -0.00606,        0, -0.00718 },
    {      0,  0.836,  3.99,  -1.33, -0.00298, -0.00139, -0.00268 },
    {  -1.29,   1.01,  3.92,   4.67,  0.00484, -0.00127,  -0.0199 }
};

/* coefficient matrices for the lambdamax test */

const double maxev_m_coef[5][5] = {
    /*   N        1        N==1       N==2     N^1/2 */
    { 6.0019, -2.75580,  0.67185,  0.114900, -2.77640 },
    { 5.9498,  0.43402,  0.04836,  0.018198, -2.36690 },
    { 5.8271, -1.64870, -1.61180, -0.259490, -1.56660 },
    { 5.8658,  2.55950, -0.34443, -0.077991, -1.75520 },
    { 5.6364, -0.90531, -3.51660, -0.479660, -0.21447 }
};

const double maxev_v_coef[5][5] = {
    { 1.8806, -15.499,  1.11360,  0.070508, 14.714 },
    { 2.2231, -7.9064,  0.58592, -0.034324, 12.058 },
    { 2.0785, -9.7846, -3.36800, -0.245280, 13.074 },
    { 1.9955, -5.5428,  1.24250,  0.419490, 12.841 },
    { 2.0899, -5.3303, -7.15230, -0.252600, 12.393 }
};

static void fill_x_asy_array (double *x, int n)
{
    x[0] = n * n;
    x[1] = n;
    x[2] = 1.0;
    x[3] = (n == 1)? 1.0 : 0.0;
    x[4] = (n == 2)? 1.0 : 0.0;
    x[5] = sqrt((double) n);
}

static void fill_x_corr_array (double *x, int n, int T)
{
    x[0] = sqrt((double) n) / T;
    x[1] = n / (double) T;
    x[2] = (n * n) / (double) (T * T);
    x[3] = (n == 1)? 1.0 / T : 0.0;
    x[4] = (n == 1)? 1.0 : 0.0;
    x[5] = (n == 2)? 1.0 : 0.0;
    x[6] = (n == 3)? 1.0 : 0.0;;
}

static gretl_matrix *johansen_nullspace (const gretl_matrix *R,
					 int *err);

/*
   Asymptotic p-values for Johansen's likelihood ratio tests
   computed using J. Doornik's gamma approximation --
   see "Approximations to the Asymptotic Distributions of
   Cointegration Tests", Journal of Economic Surveys,
   12/5, 1998, pp. 573-593.

   trace, lmax: trace and lambdamax statistics

   det: index of setup of deterministic regressors
        J_NO_CONST     = no constant
        J_REST_CONST   = restricted constant
        J_UNREST_CONST = unrestricted constant
        J_REST_TREND   = restricted trend
        J_UNREST_TREND = unrestricted trend

   n: (>= 1), the number of potentially cointegrated
   variables minus the cointegration rank under the null

   pval: on output, array of p-values for the two tests
*/

static int
gamma_LR_asy_pvals (double trace, double lmax, JohansenCode det,
                    int n, double *pval)
{
    const double *tracem = trace_m_coef[det];
    const double *tracev = trace_v_coef[det];
    const double *maxevm = maxev_m_coef[det];
    const double *maxevv = maxev_v_coef[det];
    double mt = 0, vt = 0;
    double ml = 0, vl = 0;
    double x[6];
    int i;

    fill_x_asy_array(x, n);

    for (i=0; i<6; i++) {
        mt += x[i] * tracem[i];
        vt += x[i] * tracev[i];
        if (i > 0) {
            ml += x[i] * maxevm[i-1];
            vl += x[i] * maxevv[i-1];
        }
    }

    pval[0] = gamma_cdf_comp(mt, vt, trace, 2);
    pval[1] = gamma_cdf_comp(ml, vl, lmax, 2);

    return 0;
}

/* Variant of gamma_LR_asy_pvals() with sample-size correction.
   Only the trace test is handled.
*/

static double
gamma_LR_T_pval (double trace, JohansenCode det, int n, int T)
{
    const double *tracem = trace_m_coef[det];
    const double *tracev = trace_v_coef[det];
    const double *mcorr = trace_m_corr[det];
    const double *vcorr = trace_v_corr[det];
    double mt = 0, vt = 0;
    double mtcorr = 0, vtcorr = 0;
    double x[7];
    int i;

    fill_x_asy_array(x, n);

    for (i=0; i<6; i++) {
        mt += x[i] * tracem[i];
        vt += x[i] * tracev[i];
    }

    fill_x_corr_array(x, n, T);

    for (i=0; i<7; i++) {
        mtcorr += x[i] * mcorr[i];
        vtcorr += x[i] * vcorr[i];
    }

    mt = exp(log(mt) + mtcorr);
    vt = exp(log(vt) + vtcorr);

    return gamma_cdf_comp(mt, vt, trace, 2);
}

/*
  @n1: number of potentially cointegrated variables
  @n2: number of variables that are conditioned upon
  @r:  cointegrating rank under H0.
  @T:  sample size minus number of parameters, or
       0 for the asymptotic result

  The Doornik gamma-approximation approach with correction for the
  case of a "partial system" as discussed in Harbo, Johansen, Nielsen
  and Rahbek, "Asymptotic Inference on Cointegrating Rank in Partial
  Systems" (Journal of Business and Economic Statistics 16/4, October
  1998). The critical values were re-simulated in gretl with 50,000
  replications.
*/

static double
gamma_harbo_trace_pval (double trace, JohansenCode det,
                        int n1, int n2, int r, int T)
{
    const double *tracem = trace_m_coef[det];
    const double *tracev = trace_v_coef[det];
    double mt = 0, vt = 0;
    double cov, x[7];
    int n = n1 + n2;
    int i, m = n - r;

    fill_x_asy_array(x, m);

    for (i=0; i<6; i++) {
        mt += x[i] * tracem[i];
        vt += x[i] * tracev[i];
    }

    if (T > 0) {
        const double *mcorr = trace_m_corr[det];
        const double *vcorr = trace_v_corr[det];
        double mtcorr = 0, vtcorr = 0;

        fill_x_corr_array(x, m, T);

        for (i=0; i<7; i++) {
            mtcorr += x[i] * mcorr[i];
            vtcorr += x[i] * vcorr[i];
        }

        mt = exp(log(mt) + mtcorr);
        vt = exp(log(vt) + vtcorr);
    }

    if (det == J_REST_TREND) {
        cov = -1.35;
    } else if (det == J_REST_CONST) {
        cov = -1.066;
    } else {
        cov = -1.270;
    }

    mt *= (n1 - r) / (double) (n - r);
    vt *= (n1 - r) / (double) (n - r);
    vt -= (n - n1) * (n1 - r) * cov;

    return gamma_cdf_comp(mt, vt, trace, 2);
}

/* public function accessible via gretl's plugin apparatus

   In a script, do: pv = pvalue(J, n, det, T, trace)

   Use T = 0 for the plain asymptotic result.
*/

double trace_pvalue (double trace, int n, int det, int T)
{
    if (det < J_NO_CONST || det > J_UNREST_TREND || n < 1) {
        return NADBL;
    } else {
        const double *tracem = trace_m_coef[det];
        const double *tracev = trace_v_coef[det];
        double mt = 0, vt = 0;
        double x[7];
        int i;

        fill_x_asy_array(x, n);

        for (i=0; i<6; i++) {
            mt += x[i] * tracem[i];
            vt += x[i] * tracev[i];
        }

        if (T > 0 && T < 10000) {
            const double *mcorr = trace_m_corr[det];
            const double *vcorr = trace_v_corr[det];
            double mtcorr = 0, vtcorr = 0;

            fill_x_corr_array(x, n, T);

            for (i=0; i<7; i++) {
                mtcorr += x[i] * mcorr[i];
                vtcorr += x[i] * vcorr[i];
            }

            mt *= exp(mtcorr);
            vt *= exp(vtcorr);
        }

        return gamma_cdf_comp(mt, vt, trace, 2);
    }
}

/* Remove a possible excess zero from the end of a floating point
   number printed to the given precision p (working around a bug
   in the C library).
*/

static void fix_xstr (char *s, int p)
{
    int n = strlen(s);

    if (n > p && strspn(s + n - p, "0") == p) {
        s[n-1] = 0;
    }
}

#define ABMIN 1.0e-15

/* for cointegration test: print cointegrating vectors or adjustments,
   either "raw" or re-scaled */

static void print_beta_or_alpha (const GRETL_VAR *jvar, int k,
                                 const DATASET *dset, PRN *prn,
                                 int job, int rescale)
{
    JohansenInfo *jv = jvar->jinfo;
    gretl_matrix *c = (job == V_BETA)? jv->Beta : jv->Alpha;
    int rows = gretl_matrix_rows(c);
    int vnorm = libset_get_int(VECM_NORM);
    char xstr[32], tmp[NAMETRUNC];
    int n, namelen = 8;
    int i, j, row;
    double x, y;

    if (vnorm == NORM_NONE && rescale) {
        return;
    }

    if (rescale) {
        pprintf(prn, "\n%s\n", (job == V_BETA)?
                _("renormalized beta") :
                _("renormalized alpha"));
    } else {
        pprintf(prn, "\n%s\n", (job == V_BETA)?
                _("beta (cointegrating vectors)") :
                _("alpha (adjustment vectors)"));
    }

    for (i=0; i<rows; i++) {
        vecm_beta_varname(tmp, jvar, dset, i);
        n = strlen(tmp);
        if (n > namelen) {
            namelen = n;
        }
    }

    for (i=0; i<rows; i++) {
        vecm_beta_varname(tmp, jvar, dset, i);
        pprintf(prn, "%-*s", namelen + 2, tmp);
        for (j=0; j<k; j++) {
            x = gretl_matrix_get(c, i, j);
            if (rescale) {
                row = (vnorm == NORM_FIRST)? 0 : j;
                y = gretl_matrix_get(jv->Beta, row, j);
                if (job == V_BETA) {
                    x /= y;
                } else {
                    x *= y;
                }
            }
            if (x == -0.0 || fabs(x) < ABMIN) {
                x = 0.0;
            }
            sprintf(xstr, "%#.5g", x);
            fix_xstr(xstr, 5);
            pprintf(prn, "%12s ", xstr);
        }
        pputc(prn, '\n');
    }
}

/* Calculate \alpha (adjustments) matrix as per Johansen, 1991, eqn
   2.8, p. 1554.  Not needed when doing a VECM (without restrictions
   on alpha), in which case we get \alpha via VECM_estimate_full()
   below.  We use this to support verbose output when testing a
   restriction on \beta.
*/

static int compute_alpha (JohansenInfo *jv)
{
    const gretl_matrix *B = jv->Beta;
    gretl_matrix *alpha = NULL;
    gretl_matrix *BSB = NULL;
    gretl_matrix *Tmp = NULL;
    int p = jv->S01->rows;
    int r = B->cols;
    int err = 0;

    BSB = gretl_matrix_alloc(r, r);
    Tmp = gretl_matrix_alloc(B->rows, r);
    alpha = gretl_matrix_alloc(p, r);

    if (BSB == NULL || Tmp == NULL || alpha == NULL) {
        err = E_ALLOC;
    }

    if (!err) {
        err = gretl_matrix_qform(B, GRETL_MOD_TRANSPOSE, jv->S11,
                                 BSB, GRETL_MOD_NONE);
    }

    if (!err) {
        err = gretl_invert_symmetric_matrix(BSB);
    }

    if (!err) {
        gretl_matrix_multiply(B, BSB, Tmp);
        gretl_matrix_multiply(jv->S01, Tmp, alpha);
    }

    gretl_matrix_free(BSB);
    gretl_matrix_free(Tmp);

    if (!err) {
        gretl_matrix_replace(&jv->Alpha, alpha);
    } else {
        gretl_matrix_free(alpha);
    }

    return err;
}

/* print the long-run matrix, \alpha \beta' */

static int print_long_run_matrix (const GRETL_VAR *jvar,
                                  const DATASET *dset,
                                  PRN *prn)
{
    JohansenInfo *jv = jvar->jinfo;
    gretl_matrix *Pi;
    char tmp[NAMETRUNC];
    const char *vname;
    int firstlen = 10;
    int namelen;
    double x;
    int i, j;

    Pi = gretl_matrix_alloc(jv->Alpha->rows, jv->Beta->rows);
    if (Pi == NULL) {
        return E_ALLOC;
    }

    namelen = max_namelen_in_list(jvar->ylist, dset);
    if (firstlen <= namelen) {
        firstlen = namelen + 1;
    }
    if (namelen < 12) {
        namelen = 12;
    }

    gretl_matrix_multiply_mod(jv->Alpha, GRETL_MOD_NONE,
                              jv->Beta, GRETL_MOD_TRANSPOSE,
                              Pi, GRETL_MOD_NONE);

    pprintf(prn, "%s\n", _("long-run matrix (alpha * beta')"));

    vname = dset->varname[jvar->ylist[1]];
    maybe_trim_varname(tmp, vname);

    pprintf(prn, "%*s", firstlen + namelen, tmp);
    for (j=1; j<Pi->cols; j++) {
        vecm_beta_varname(tmp, jvar, dset, j);
        pprintf(prn, "%*s", namelen + 1, tmp);
    }

    pputc(prn, '\n');

    for (i=0; i<Pi->rows; i++) {
        vname = dset->varname[jvar->ylist[i+1]];
        maybe_trim_varname(tmp, vname);
        pprintf(prn, "%-*s", firstlen, tmp);
        for (j=0; j<Pi->cols; j++) {
            x = gretl_matrix_get(Pi, i, j);
            if (fabs(x) < 0.5e-14) {
                x = 0.0;
            }
            pprintf(prn, "%#*.5g ", namelen, x);
        }
        pputc(prn, '\n');
    }

    pputc(prn, '\n');

    gretl_matrix_free(Pi);

    return 0;
}

/* Compute Hamilton's Omega (Johansen 1991 calls it Lambda): the
   cross-equation variance matrix.
*/

static int compute_omega (GRETL_VAR *vecm)
{
    if (vecm->S == NULL) {
        vecm->S = gretl_matrix_alloc(vecm->neqns, vecm->neqns);
        if (vecm->S == NULL) {
            return E_ALLOC;
        }
    }

    gretl_matrix_multiply_mod(vecm->E, GRETL_MOD_TRANSPOSE,
                              vecm->E, GRETL_MOD_NONE,
                              vecm->S, GRETL_MOD_NONE);

    gretl_matrix_divide_by_scalar(vecm->S, vecm->T);

    return 0;
}

static void gretl_matrix_I (gretl_matrix *A, int n)
{
    int i;

    gretl_matrix_zero(A);
    for (i=0; i<n; i++) {
        gretl_matrix_set(A, i, i, 1.0);
    }
}

#define lag_wanted(v, i) (v->lags == NULL || in_gretl_list(v->lags, i))

/* After doing OLS estimation of the VECM conditional on \beta: copy
   the coefficients on the lagged differences (i.e. form the \Gamma
   matrices) so we can compute the VAR representation */

static void copy_coeffs_to_Gamma (GRETL_VAR *vecm, gretl_matrix **G)
{
    int nl = var_n_lags(vecm);
    int i, j, k, h;
    double x;

    for (i=0; i<vecm->neqns; i++) {
        for (k=0; k<vecm->order; k++) {
            if (!lag_wanted(vecm, k+1)) {
                gretl_matrix_zero(G[k]);
                continue;
            }
            h = k + vecm->ifc;
            /* successive lags (distinct \Gamma_i matrices) */
            for (j=0; j<vecm->neqns; j++) {
                /* successive \Delta x_j */
                x = gretl_matrix_get(vecm->B, h, i);
                gretl_matrix_set(G[k], i, j, x);
                h += nl;
            }
        }
    }

#if JDEBUG > 1
    for (k=0; k<vecm->order; k++) {
        char msg[32];

        sprintf(msg, "Gamma matrix, lag %d", k+1);
        gretl_matrix_print(G[k], msg);
    }
#endif
}

/* \Pi, as will be used in forming the VAR representation */

static void form_Pi (GRETL_VAR *v, gretl_matrix *Pi)
{
    gretl_matrix_multiply_mod(v->jinfo->Alpha, GRETL_MOD_NONE,
                              v->jinfo->Beta, GRETL_MOD_TRANSPOSE,
                              Pi, GRETL_MOD_NONE);
}

/* After doing OLS estimation of the VECM conditional on \beta:
   copy the coefficients on the EC terms (\beta' X) into the \alpha
   matrix.
*/

static int OLS_to_alpha (GRETL_VAR *v)
{
    int rank = v->jinfo->rank;
    int pos = v->ncoeff - rank;
    double x;
    int i, j;

    for (i=0; i<v->neqns; i++) {
        for (j=0; j<rank; j++) {
            x = gretl_matrix_get(v->B, pos + j, i);
            gretl_matrix_set(v->jinfo->Alpha, i, j, x);
        }
    }

    return 0;
}

/* flags for controlling "full" estimation of VECM */

enum {
    ALPHA_RESTRICTED = 1 << 0,
    BOOTSTRAPPING    = 1 << 1
};

#define bootstrap(f)        (f & BOOTSTRAPPING)
#define alpha_restricted(f) (f & ALPHA_RESTRICTED)
#define estimate_alpha(f)   (!(f & ALPHA_RESTRICTED))

/* VAR representation: transcribe the coefficient matrix A_i (for lag
   i) into its place in the full VAR coefficient matrix, A
*/

static void add_Ai_to_VAR_A (gretl_matrix *Ai, GRETL_VAR *vecm,
                             int k, int flags)
{
    int i, j, offset = k * vecm->neqns;
    int tr = bootstrap(flags);
    double x;

    for (i=0; i<vecm->neqns; i++) {
        for (j=0; j<vecm->neqns; j++) {
            x = gretl_matrix_get(Ai, i, j);
            if (tr) {
                gretl_matrix_set(vecm->A, j + offset, i, x);
            } else {
                gretl_matrix_set(vecm->A, i, j + offset, x);
            }
        }
    }
}

/* write pre-computed ML alpha into model structs */

static void transcribe_alpha (GRETL_VAR *v)
{
    MODEL *pmod;
    double aij, sij = NADBL;
    int r = jrank(v);
    int k = v->ncoeff - r;
    int i, j;

    for (i=0; i<v->neqns; i++) {
        pmod = v->models[i];
        for (j=0; j<r; j++) {
            aij = gretl_matrix_get(v->jinfo->Alpha, i, j);
            if (v->jinfo->Ase != NULL) {
                sij = gretl_matrix_get(v->jinfo->Ase, i, j);
            }
            pmod->coeff[k+j] = aij;
            pmod->sderr[k+j] = sij;
        }
    }
}

/* The X (data) and B (coefficient) matrices may need expanding
   to take account of the EC terms */

static int vecm_check_size (GRETL_VAR *v, int flags)
{
    int xc = (v->X != NULL)? v->X->cols : 0;
    int err = 0;

    if (bootstrap(flags)) {
        /* We're re-visiting this function: v->ncoeff will
           already be full-size, as will the v->X matrix.
        */
        int dim = v->ncoeff;

        if (alpha_restricted(flags)) {
            dim -= gretl_VECM_rank(v);
        }

        if (dim > 0 && (v->X == NULL || v->B == NULL)) {
            gretl_errmsg_set("vecm_check_size: X and/or B wrong");
            return E_DATA;
        } else if (dim == 0) {
            /* FIXME! */
            if (v->X != NULL) v->X->cols = 0;
            if (v->B != NULL) v->B->rows = 0;
#if 0
            fprintf(stderr, "vecm_check_size: weird case?\n");
#endif
        } else {
            v->X->cols = dim;
            v->B->rows = dim;
        }
        return 0;
    }

#if JDEBUG
    fprintf(stderr, "vecm_check_size: ncoeff: %d -> %d; xc = %d\n",
            v->ncoeff, v->ncoeff + jrank(v), xc);
    if (v->B != NULL) {
        fprintf(stderr, "v->B->rows = %d\n", v->B->rows);
    }
#endif

    v->ncoeff += jrank(v);

    if (estimate_alpha(flags)) {
        xc += jrank(v);
    }

#if JDEBUG
    fprintf(stderr, "after alpha: ncoeff = %d, xc = %d\n", v->ncoeff, xc);
#endif


    if (xc == 0) {
        /* nothing to be done */
        return 0;
    }

    if (v->X == NULL) {
        v->X = gretl_matrix_alloc(v->T, xc);
        if (v->X == NULL) {
            err = E_ALLOC;
        } else {
            /* record allocated size */
            v->xcols = xc;
        }
    } else if (v->X->cols < xc) {
        if (v->xcols == xc) {
            gretl_matrix_reuse(v->X, -1, xc);
        } else {
            err = gretl_matrix_realloc(v->X, v->T, xc);
            if (!err) {
                /* record revised full size */
                v->xcols = xc;
            }
        }
    }

    if (!err) {
        int rows = v->ncoeff; /* 2021-01-06: was xc */

        if (v->B == NULL) {
            v->B = gretl_matrix_alloc(rows, v->neqns);
            if (v->B == NULL) {
                err = E_ALLOC;
            }
        } else if (v->B->rows < rows) {
            err = gretl_matrix_realloc(v->B, rows, v->neqns);
        }
    }

    return err;
}

/* For estimating both alpha and Gamma: add the EC terms into
   the X data matrix */

static int add_EC_terms_to_X (GRETL_VAR *v, gretl_matrix *X,
                              const DATASET *dset)
{
    const gretl_matrix *B = v->jinfo->Beta;
    int rank = jrank(v);
    int k, k0 = v->ncoeff - rank;
    double xt, bxt, bij;
    int i, ii, j, s, t;
    int err = 0;

    for (j=0, k=k0; j<rank; j++, k++) {
        for (t=v->t1, s=0; t<=v->t2; t++, s++) {
            bxt = 0.0;
            ii = 0;

            /* beta * X(t-1) */
            for (i=0; i<v->neqns; i++) {
                xt = dset->Z[v->ylist[i+1]][t-1];
                bij = gretl_matrix_get(B, ii++, j);
                bxt += bij * xt;
            }

            /* restricted const or trend */
            if (auto_restr(v)) {
                bij = gretl_matrix_get(B, ii++, j);
                if (jcode(v) == J_REST_TREND) {
                    bij *= t;
                }
                bxt += bij;
            }

            /* restricted exog vars */
            if (v->rlist != NULL) {
                for (i=0; i<v->rlist[0]; i++) {
                    xt = dset->Z[v->rlist[i+1]][t]; /* was t-1 */
                    bij = gretl_matrix_get(B, ii++, j);
                    bxt += bij * xt;
                }
            }

            gretl_matrix_set(X, s, k, bxt);
        }
    }

    return err;
}

/* preparing for OLS conditional on beta: construct the
   appropriate dependent variable matrix, Y
*/

static int make_vecm_Y (GRETL_VAR *v, const DATASET *dset,
                        const gretl_matrix *Pi)
{
    int i, s, t, vi, vj;
    double pij, yti, xti;
    int err = 0;

    if (v->Y == NULL) {
        fprintf(stderr, "make_vecm_Y: v->Y is NULL\n");
        return E_DATA;
    }

    if (Pi == NULL) {
        /* "Y" is composed of plain DYt */
        for (i=0; i<v->neqns; i++) {
            vi = v->ylist[i+1];
            s = 0;
            for (t=v->t1; t<=v->t2; t++) {
                yti = dset->Z[vi][t] - dset->Z[vi][t-1];
                gretl_matrix_set(v->Y, s++, i, yti);
            }
        }
    } else {
        /* net out \alpha: "Y" = DY_t - \Pi Y*_t */
        int j, k, wexo, p1 = v->jinfo->Beta->rows;

        for (i=0; i<v->neqns; i++) {
            wexo = 1;
            vi = v->ylist[i+1];
            s = 0;
            for (t=v->t1; t<=v->t2; t++) {
                /* first difference */
                yti = dset->Z[vi][t] - dset->Z[vi][t-1];
                for (j=0; j<p1; j++) {
                    pij = gretl_matrix_get(Pi, i, j);
                    if (pij != 0.0) {
                        if (j < v->neqns) {
                            /* lagged Y level */
                            wexo = 0;
                            vj = v->ylist[j+1];
                            xti = dset->Z[vj][t-1];
                        } else if (j == v->neqns && auto_restr(v)) {
                            xti = (jcode(v) == J_REST_TREND)? t : 1;
                        } else {
                            k = j - v->ylist[0] - auto_restr(v) + 1;
                            vj = v->rlist[k];
                            xti = dset->Z[vj][t]; /* was t-1 */
                        }
                        yti -= pij * xti;
                    }
                }
                gretl_matrix_set(v->Y, s++, i, yti);
            }
            if (wexo) {
                fprintf(stderr, "make_vecm_Y: var %d is weakly exogenous\n", i);
            }
        }
    }

    return err;
}

/* As in PcGive: df = T - c, where c is the "average number of
   estimated parameters per equation, rounded towards zero".
   Note that this function is _not_ invoked in the case of
   "general" restriction on alpha/beta; in that case the df
   calculation is handled specially.
*/

static void vecm_set_df (GRETL_VAR *v, const gretl_matrix *H,
                         const gretl_matrix *R)
{
    int p = v->neqns;
    int r = v->jinfo->rank;
    int p1 = v->jinfo->Beta->rows;
    int Kpi, K = 0;
    double c;

    if (r == 0) {
        /* full-rank VAR */
        r = p;
    }

    /* lagged differences */
    K += var_n_lags(v) * p;

    /* deterministic stuff */
    K += v->jinfo->seasonals;
    if (jcode(v) > J_REST_CONST) {
        K++;
    }
    if (jcode(v) == J_UNREST_TREND) {
        K++;
    }

    /* unrestricted exogenous vars? */
    if (v->xlist != NULL) {
        K += v->xlist[0];
    }

    K *= p;

    /* df loss in \Pi = \alpha \beta': see Johansen (1995) chapter 7 */
    Kpi = (p + p1 - r) * r;
    if (H != NULL || R != NULL) {
        Kpi -= v->jinfo->lrdf;
    }
    K += Kpi;

    c = K / (double) p;
    v->df = v->T - floor(c);

#if JDEBUG
    fprintf(stderr, "vecm_set_df: T = %d, global K = %d, c = %g, df = %d\n",
            v->T, K, c, v->df);
#endif
}

/* In the case where the EC coefficients were not estimated via
   OLS but precomputed in v->jinfo->Alpha, augment the B matrix
   with the alpha values. We need to add a number of rows equal
   to v->ncoeff minus the current B->rows, then transcribe the
   transpose of v->jinfo->Alpha into the new rows.
*/

static int var_B_insert_alpha (GRETL_VAR *v)
{
    gretl_matrix *A = v->jinfo->Alpha;
    gretl_matrix *B = v->B;
    double aji, *dest, *src;
    int r_orig = B->rows;
    int r_full = v->ncoeff;
    size_t sz;
    int i, j;

    /* sanity check */
    if (r_full - r_orig <= 0 || A->rows != B->cols) {
        fprintf(stderr, "var_B_insert_alpha: matrix sizes wrong!\n");
	return E_DATA;
    }

    sz = r_orig * sizeof *B->val;

    for (j=B->cols-1; j>0; j--) {
	dest = B->val + r_full * j;
	src = B->val + r_orig * j;
	memmove(dest, src, sz);
    }

    gretl_matrix_reuse(B, r_full, B->cols);

    for (j=0; j<B->cols; j++) {
	for (i=0; i<A->cols; i++) {
	    aji = gretl_matrix_get(A, j, i);
	    gretl_matrix_set(B, r_orig + i, j, aji);
	}
    }

    return 0;
}

/* The following is designed to accommodate the case where alpha is
   restricted, in which case we can't just run OLS conditional on
   beta.

   DY_t = \Pi Y*_{t-1} + \sum_{i=1}^{k-1}\Gamma_i DY_{t-1} + ...

   Subtract \Pi Y*_t from both sides, call DY_t - \Pi Y*_t "Y",
   and call the lagged DYs "X".  Regress "Y" on "X" to find
   estimates of the \Gammas, etc.

   But the function also handles the case where \alpha will
   be estimated along with \Gamma.
*/

static int
VECM_estimate_full (GRETL_VAR *v, const gretl_restriction *rset,
                    const DATASET *dset, int flags)
{
    gretl_matrix *beta = v->jinfo->Beta;
    gretl_matrix *Pi = NULL;
    gretl_matrix *Ai = NULL;
    gretl_matrix **G = NULL;
    int order = v->order;
    int xc, n = v->neqns;
    int use_XTX = 0;
    int fix_Y = 0;
    int i, err;

#if JDEBUG
    fprintf(stderr, "VECM_estimate_full: %s\n",
            (estimate_alpha(flags))? "including alpha in estimation" :
            "netting out the EC terms");
#endif

    if (alpha_restricted(flags) && v->jinfo->Alpha == NULL) {
        /* error: alpha must be pre-computed */
        fprintf(stderr, " alpha must be pre-computed, but Alpha is NULL\n");
        return E_DATA;
    }

    err = vecm_check_size(v, flags);
    if (err) {
        fprintf(stderr, " error %d from vecm_check_size\n", err);
        return err;
    }

    xc = (v->X != NULL)? v->X->cols : 0;

    Pi = gretl_matrix_alloc(n, beta->rows);
    Ai = gretl_matrix_alloc(n, n);
    if (Pi == NULL || Ai == NULL) {
        err = E_ALLOC;
    }

    if (!err && order > 0) {
        G = gretl_matrix_array_new_with_size(order, n, n);
        if (G == NULL) {
            err = E_ALLOC;
        }
    }

    if (!err && alpha_restricted(flags)) {
        form_Pi(v, Pi);
    }

    if (!err && v->Y != NULL) {
        if (alpha_restricted(flags)) {
            fix_Y = 1;
            err = make_vecm_Y(v, dset, Pi);
        } else {
            err = make_vecm_Y(v, dset, NULL);
        }
        if (err) {
            fprintf(stderr, " error %d from make_vecm_Y\n", err);
        }
    }

    if (!err && estimate_alpha(flags)) {
        err = add_EC_terms_to_X(v, v->X, dset);
        if (err) {
            fprintf(stderr, " error %d from add_EC_terms_to_X\n", err);
        }
    }

    if (!err) {
        if (xc > 0) {
            /* run the regressions */
            if (v->B->rows > xc) {
                gretl_matrix_reuse(v->B, xc, -1);
            }
            if (bootstrap(flags)) {
                err = gretl_matrix_multi_ols(v->Y, v->X, v->B, v->E, NULL);
            } else {
                if (v->XTX != NULL) {
                    gretl_matrix_replace(&v->XTX, NULL);
                }
                err = gretl_matrix_multi_SVD_ols(v->Y, v->X, v->B, v->E, &v->XTX);
                use_XTX = 1;
            }
            if (err) {
                fprintf(stderr, " error %d from matrix_ols (xc = %d)\n", err, xc);
            }
        } else if (v->Y != NULL) {
            /* nothing to estimate, with alpha already in hand */
            gretl_matrix_copy_values(v->E, v->Y);
        } else {
            gretl_matrix_zero(v->E);
        }
    }

    if (!err && order > 0) {
        copy_coeffs_to_Gamma(v, G);
    }

    if (!err) {
        if (estimate_alpha(flags)) {
            err = OLS_to_alpha(v);
            if (!err) {
                form_Pi(v, Pi);
            }
        } else if (xc < v->ncoeff && v->B != NULL) {
            /* transcribe EC terms to v->B */
            var_B_insert_alpha(v);
        }
    }

    if (err) {
        goto bailout;
    }

    if (Pi->cols > n) {
        gretl_matrix_reuse(Pi, -1, n);
    }

    if (order == 0) {
        gretl_matrix_I(Ai, n);
        gretl_matrix_add_to(Ai, Pi);
        add_Ai_to_VAR_A(Ai, v, 0, flags);
    } else {
        for (i=0; i<=order; i++) {
            if (i == 0) {
                gretl_matrix_I(Ai, n);
                gretl_matrix_add_to(Ai, Pi);
                gretl_matrix_add_to(Ai, G[0]);
            } else if (i == order) {
                gretl_matrix_zero(Ai);
                gretl_matrix_subtract_from(Ai, G[i-1]);
            } else {
                gretl_matrix_copy_values(Ai, G[i]);
                gretl_matrix_subtract_from(Ai, G[i-1]);
            }
#if JDEBUG
            fprintf(stderr, "Ai matrix, lag %d\n\n", i+1);
            gretl_matrix_print(Ai, NULL);
#endif
            add_Ai_to_VAR_A(Ai, v, i, flags);
        }
    }

#if JDEBUG
    gretl_matrix_print(v->A, "vecm->A");
#endif

    if (!err && !bootstrap(flags)) {
        const gretl_matrix *XTX = NULL;

        if (use_XTX || estimate_alpha(flags)) {
            XTX = v->XTX;
        }
        if (fix_Y) {
            /* ensure "Y" holds plain differences */
            make_vecm_Y(v, dset, NULL);
        }
        transcribe_VAR_models(v, dset, XTX);
        if (alpha_restricted(flags)) {
            transcribe_alpha(v);
        }
    }

    if (!err && G != NULL) {
        /* see Johansen (1995, p. 45) */
        gretl_matrix *Tmp;

        Tmp = gretl_identity_matrix_new(n);
        for (i=0; i<order; i++) {
            gretl_matrix_subtract_from(Tmp, G[i]);
        }
        v->jinfo->Gamma = Tmp;
    }

 bailout:

    gretl_matrix_free(Pi);
    gretl_matrix_free(Ai);
    gretl_matrix_array_free(G, order);

    return err;
}

/* Print both "raw" and re-scaled versions of the beta and alpha
   matrices (cointegrating vectors and vectors of adjustments
   respectively).
*/

static int
print_beta_and_alpha (const GRETL_VAR *jvar, gretl_matrix *evals, int h,
                      const DATASET *dset, PRN *prn)
{
    int i, err = 0;

    pputs(prn, _("eigenvalue"));
    for (i=0; i<h; i++) {
        pprintf(prn, "%#12.5g ", evals->val[i]);
    }
    pputc(prn, '\n');

    /* "raw" vectors */
    print_beta_or_alpha(jvar, h, dset, prn, V_BETA, 0);
    print_beta_or_alpha(jvar, h, dset, prn, V_ALPHA, 0);

    /* re-scaled versions */
    print_beta_or_alpha(jvar, h, dset, prn, V_BETA, 1);
    print_beta_or_alpha(jvar, h, dset, prn, V_ALPHA, 1);

    pputc(prn, '\n');

    return err;
}

/*
   renormalize \beta such that its uppermost submatrix of
   size rank * rank is the identity matrix:

   \beta' = [ I | *free elements* ]
*/

static int phillips_normalize_beta (GRETL_VAR *vecm)
{
    gretl_matrix *c = NULL;
    gretl_matrix *beta_c = NULL;

    int r = jrank(vecm);
    int n = gretl_matrix_rows(vecm->jinfo->Beta);
    int i, j, err = 0;

    double x;

    c = gretl_matrix_alloc(r, r);
    beta_c = gretl_matrix_alloc(n, r);
    if (c == NULL || beta_c == NULL) {
        err = E_ALLOC;
        goto bailout;
    }

    for (i=0; i<r; i++) {
        for (j=0; j<r; j++) {
            x = gretl_matrix_get(vecm->jinfo->Beta, i, j);
            gretl_matrix_set(c, i, j, x);
        }
    }

    /* form \beta_c = \beta c^{-1} */
    err = gretl_invert_general_matrix(c);
    if (err) {
        fprintf(stderr, "phillips_normalize_beta: c is singular\n");
        goto bailout;
    }

    gretl_matrix_multiply(vecm->jinfo->Beta, c, beta_c);

    /* correct rounding error: set true zeros in \beta_c */
    for (i=0; i<n; i++) {
        for (j=0; j<r; j++) {
            if (i >= r) {
                if (gretl_matrix_get(beta_c, i, j) == -0) {
                    gretl_matrix_set(beta_c, i, j, 0);
                }
            } else if (i == j) {
                gretl_matrix_set(beta_c, i, j, 1.0);
            } else {
                gretl_matrix_set(beta_c, i, j, 0.0);
            }
        }
    }

#if JDEBUG
    gretl_matrix_print(vecm->jinfo->Beta, "original beta");
    gretl_matrix_print(beta_c, "beta_c = beta * c^{-1}");
#endif

    gretl_matrix_copy_values(vecm->jinfo->Beta, beta_c);

 bailout:

    gretl_matrix_free(c);
    gretl_matrix_free(beta_c);

    return err;
}

static int col_normalize_beta (GRETL_VAR *vecm, int vnorm)
{
    gretl_matrix *B = vecm->jinfo->Beta;
    double x, den;
    int i, j, row;

    for (j=0; j<B->cols; j++) {
        row = (vnorm == NORM_DIAG)? j : 0;
        den = gretl_matrix_get(B, row, j);
        if (den != 0.0) {
            for (i=0; i<B->rows; i++) {
                x = gretl_matrix_get(B, i, j);
                gretl_matrix_set(B, i, j, x / den);
            }
        }
    }

    return 0;
}

static int normalize_beta (GRETL_VAR *vecm, const gretl_matrix *H,
                           int *do_stderrs)
{
    int vnorm = libset_get_int(VECM_NORM);

    if (vnorm == NORM_NONE) {
        if (do_stderrs != NULL) {
            *do_stderrs = 0;
        }
        return 0;
    }

    if (H == NULL) {
        if (vnorm == NORM_PHILLIPS) {
            return phillips_normalize_beta(vecm);
        } else {
            if (do_stderrs != NULL) {
                *do_stderrs = 0;
            }
            return col_normalize_beta(vecm, vnorm);
        }
    } else {
        gretl_matrix *B = vecm->jinfo->Beta;

        if (B->cols == 1) {
            double den = B->val[0];
            int i;

            if (den != 0.0) {
                for (i=0; i<B->rows; i++) {
                    if (B->val[i] != 0) {
                        B->val[i] /= den;
                    }
                }
            }
        }
    }

    return 0;
}

static int restricted_beta_se (GRETL_VAR *v, int r, int p1)
{
    double x;
    int i;

    v->jinfo->Bse = gretl_matrix_alloc(p1, r);
    if (v->jinfo->Bse == NULL) {
        return E_ALLOC;
    }

    for (i=0; i<v->jinfo->Bvar->rows; i++) {
        x = gretl_matrix_get(v->jinfo->Bvar, i, i);
        v->jinfo->Bse->val[i] = sqrt(x);
    }

    return 0;
}

static int restricted_beta_variance (GRETL_VAR *vecm,
                                     gretl_matrix *Hin)
{
    gretl_matrix *H = NULL;
    gretl_matrix *O = NULL;
    gretl_matrix *aOa = NULL;
    gretl_matrix *K = NULL;
    gretl_matrix *Vphi = NULL;
    int r = jrank(vecm);
    int p = vecm->neqns;
    int p1 = p + n_restricted_terms(vecm);
    int nb = r * p1;
    int freeH = 0;
    int err = 0;

    if (r > 1) {
        H = gretl_matrix_I_kronecker_new(r, Hin, &err);
        if (err) {
            return err;
        }
        freeH = 1;
    } else {
        H = Hin;
    }

    clear_gretl_matrix_err();

    O = gretl_matrix_copy(vecm->S);
    aOa = gretl_matrix_alloc(r, r);
    K = gretl_matrix_alloc(nb, nb);
    Vphi = gretl_matrix_alloc(H->cols, H->cols);

    err = get_gretl_matrix_err();
    if (err) {
        goto bailout;
    }

    err = gretl_invert_symmetric_matrix(O);

    if (!err) {
        err = gretl_matrix_qform(vecm->jinfo->Alpha, GRETL_MOD_TRANSPOSE,
                                 O, aOa, GRETL_MOD_NONE);
    }

    if (!err) {
        err = gretl_matrix_kronecker_product(aOa, vecm->jinfo->S11, K);
    }

    if (!err) {
        gretl_matrix_qform(H, GRETL_MOD_TRANSPOSE, K,
                           Vphi, GRETL_MOD_NONE);
    }

    if (!err) {
        err = gretl_invert_symmetric_matrix(Vphi);
    }

    if (!err) {
        gretl_matrix_divide_by_scalar(Vphi, vecm->df);
    }

    if (!err) {
        vecm->jinfo->Bvar = gretl_matrix_alloc(nb, nb);
        if (vecm->jinfo->Bvar == NULL) {
            err = E_ALLOC;
        } else {
            gretl_matrix_qform(H, GRETL_MOD_NONE, Vphi,
                               vecm->jinfo->Bvar,
                               GRETL_MOD_NONE);
        }
    }

    if (!err) {
        err = restricted_beta_se(vecm, r, p1);
    }

 bailout:

    gretl_matrix_free(O);
    gretl_matrix_free(aOa);
    gretl_matrix_free(K);
    gretl_matrix_free(Vphi);

    if (freeH) {
        gretl_matrix_free(H);
    }

    return err;
}

/* VECM: compute the variance of the estimator of \beta, after doing
   Phillips normalization */

static int beta_variance (GRETL_VAR *vecm)
{
    gretl_matrix *O = NULL;
    gretl_matrix *aOa = NULL;
    gretl_matrix *HSH = NULL;
    double x;
    int r = jrank(vecm);
    int p1 = gretl_matrix_rows(vecm->jinfo->Beta);
    int nh = p1 - r;
    int i, j, k, err = 0;

    O = gretl_matrix_copy(vecm->S);
    aOa = gretl_matrix_alloc(r, r);
    HSH = gretl_matrix_alloc(nh, nh);

    if (O == NULL || aOa == NULL || HSH == NULL) {
        err = E_ALLOC;
        goto bailout;
    }

    /* compute \alpha' \Omega^{-1} \alpha */

    err = gretl_invert_symmetric_matrix(O);
    if (err) {
        goto bailout;
    }

    gretl_matrix_qform(vecm->jinfo->Alpha, GRETL_MOD_TRANSPOSE, O,
                       aOa, GRETL_MOD_NONE);

#if JDEBUG
    gretl_matrix_print(vecm->S, "vecm->S");
    gretl_matrix_print(O, "O = inverse(vecm->S)");
    gretl_matrix_print(vecm->jinfo->Alpha, "alpha_c");
    gretl_matrix_print(aOa, "aOa = alpha_c' * O * alpha_c");
#endif

    /* form H'*S11*H: just keep the south-east corner */

    for (i=r; i<p1; i++) {
        for (j=r; j<p1; j++) {
            x = gretl_matrix_get(vecm->jinfo->S11, i, j);
            gretl_matrix_set(HSH, i - r, j - r, x);
        }
    }

#if JDEBUG
    gretl_matrix_print(vecm->jinfo->S11, "full S11");
    gretl_matrix_print(HSH, "H'*S11*H");
#endif

    vecm->jinfo->Bvar = gretl_matrix_kronecker_product_new(aOa, HSH, &err);
    if (err) {
        goto bailout;
    }

    err = gretl_invert_symmetric_matrix(vecm->jinfo->Bvar);
    if (err) {
        goto bailout;
    }

    vecm->jinfo->Bse = gretl_zero_matrix_new(p1, r);
    if (vecm->jinfo->Bse == NULL) {
        err = E_ALLOC;
        goto bailout;
    }

    /* note use of df rather than T here */
    gretl_matrix_divide_by_scalar(vecm->jinfo->Bvar, vecm->df);

    k = 0;
    for (j=0; j<r; j++) {
        /* cointegrating vector j */
        for (i=r; i<p1; i++) {
            x = gretl_matrix_get(vecm->jinfo->Bvar, k, k);
            gretl_matrix_set(vecm->jinfo->Bse, i, j, sqrt(x));
            k++;
        }
    }

#if JDEBUG
    gretl_matrix_print(vecm->jinfo->Bvar, "var(beta)");
    gretl_matrix_print(vecm->jinfo->Bse, "se(beta)");
#endif

 bailout:

    gretl_matrix_free(O);
    gretl_matrix_free(aOa);
    gretl_matrix_free(HSH);

    return err;
}

int johansen_ll_calc (GRETL_VAR *jvar, const gretl_matrix *evals)
{
    gretl_matrix *S00;
    int n = jvar->neqns;
    int r = jrank(jvar);
    int h = (r > 0)? r : n;
    int i, err = 0;

    S00 = gretl_matrix_copy(jvar->jinfo->S00);

    if (S00 == NULL) {
        err = E_ALLOC;
    } else {
        double ldet = gretl_matrix_log_determinant(S00, &err);

        jvar->ll = n * (1.0 + LN_2_PI) + ldet;
        for (i=0; i<h; i++) {
            jvar->ll += log(1.0 - evals->val[i]);
        }
        jvar->ll *= -(jvar->T / 2.0);
        gretl_matrix_free(S00);
    }

    return err;
}

static int vecm_ll_stats (GRETL_VAR *vecm)
{
    int T = vecm->T;
    int g = vecm->neqns;
    int k = g * (vecm->order + 1); /* FIXME gappy */
    int err = 0;

    vecm->ldet = gretl_vcv_log_determinant(vecm->S, &err);
    if (err) {
        return err;
    }

    k += vecm->jinfo->seasonals;

    /* FIXME: is the following right for k? */

    if (jcode(vecm) >= J_UNREST_CONST) {
        k++;
    }
    if (jcode(vecm) == J_UNREST_TREND) {
        k++;
    }
    if (vecm->xlist != NULL) {
        k += vecm->xlist[0];
    }

    k *= g;

    vecm->AIC = (-2.0 * vecm->ll + 2.0 * k) / T;
    vecm->BIC = (-2.0 * vecm->ll + log(T) * k) / T;
    vecm->HQC = (-2.0 * vecm->ll + 2.0 * log(log(T)) * k) / T;

    VAR_portmanteau_test(vecm);

    return 0;
}

static void coint_test_print_exog (GRETL_VAR *jvar,
                                   const DATASET *dset,
                                   PRN *prn)
{
    int i, vi;

    if (jvar->xlist != NULL && jvar->xlist[0] > 0) {
        pprintf(prn, "\n%s: ", _("Exogenous regressor(s)"));
        for (i=1; i<=jvar->xlist[0]; i++) {
            vi = jvar->xlist[i];
            pprintf(prn, "%s ", dset->varname[vi]);
        }
    }
}

static int
compute_coint_test (GRETL_VAR *jvar, const DATASET *dset,
                    gretlopt opt, PRN *prn)
{
    gretl_matrix *evals = jvar->jinfo->evals;
    gretl_matrix *tests;
    gretl_matrix *pvals;
    int T = jvar->T;
    int n = jvar->neqns;
    int nexo = 0, nrexo = 0;
    int partial = 0;
    double llc, trace, lmax;
    double cumeig = 0.0;
    int jcase, i;

    tests = gretl_matrix_alloc(n, 2);
    pvals = gretl_matrix_alloc(n, 2);

    if (tests == NULL || pvals == NULL) {
        gretl_matrix_free(tests);
        gretl_matrix_free(pvals);
        return E_ALLOC;
    }

    for (i=n-1; i>=0; i--) {
        lmax = -T * log(1.0 - evals->val[i]);
        cumeig += lmax;
        trace = cumeig;
        gretl_matrix_set(tests, i, 0, trace);
        gretl_matrix_set(tests, i, 1, lmax);
    }

    if (jvar->xlist != NULL) {
        nexo = jvar->xlist[0];
    }

    if (jvar->rlist != NULL) {
        nrexo = jvar->rlist[0];
    }

    jcase = jcode(jvar);

    if (nrexo > 0 && jcase != J_UNREST_TREND) {
        /* use gamma-approx for partial system */
        partial = 1;
    }

    print_Johansen_test_case(jcase, prn);
    if (nexo > 0) {
        coint_test_print_exog(jvar, dset, prn);
    }
    pputc(prn, '\n');

    llc = (1.0 + LN_2_PI) * jvar->T;
    pprintf(prn, "\n%s = %g (%s: %g)\n", _("Log-likelihood"),
            jvar->ll + llc, _("including constant term"), jvar->ll);

    /* make sure the df is calculated */
    vecm_set_df(jvar, NULL, NULL);

    if (partial || nexo > 0) {
        pputc(prn, '\n');
        pputs(prn, _("Cointegration tests, ignoring exogenous variables"));
    }

    pprintf(prn, "\n%s %s %s  %s  %s  %s\n", _("Rank"), _("Eigenvalue"),
            _("Trace test"), _("p-value"),
            _("Lmax test"), _("p-value"));

    for (i=0; i<n; i++) {
        double pv[2] = {0};

        trace = gretl_matrix_get(tests, i, 0);
        lmax = gretl_matrix_get(tests, i, 1);
        gamma_LR_asy_pvals(trace, lmax, jcase, n - i, pv);
        pprintf(prn, "%4d%#11.5g%#11.5g [%6.4f]%#11.5g [%6.4f]\n",
                i, evals->val[i], trace, pv[0], lmax, pv[1]);
        gretl_matrix_set(pvals, i, 0, pv[0]);
        gretl_matrix_set(pvals, i, 1, pv[1]);
    }

    if (partial) {
        pputc(prn, '\n');
        pprintf(prn, _("Cointegration tests conditional on %d I(1) variable(s)"),
                       nrexo);
        pprintf(prn, "\n%s %s %s  %s  %s\n", _("Rank"), _("Eigenvalue"),
                _("Trace test"), _("p-value"), _("pval(T)"));

        for (i=0; i<n; i++) {
           trace = gretl_matrix_get(tests, i, 0);
            double pv;
            int det = jcase;

            /* bodge for comparability with CATS */
            if (jcase == J_UNREST_CONST) {
                det = J_REST_TREND;
            }
            pv = gamma_harbo_trace_pval(trace, det, n, nrexo, i, 0);
            pprintf(prn, "%4d%#11.5g%#11.5g [%6.4f]", i,
                    evals->val[i], trace, pv);
            gretl_matrix_set(pvals, i, 0, pv);
            gretl_matrix_set(pvals, i, 1, NADBL);
            pv = gamma_harbo_trace_pval(trace, det, n, nrexo, i, jvar->df);
            pprintf(prn, " [%6.4f]\n", pv);
        }
    } else {
        /* !partial */
        double pv;

        pputc(prn, '\n');
        pprintf(prn, _("Corrected for sample size (df = %d)"), jvar->df);
        pprintf(prn, "\n%s %s %s\n", _("Rank"), _("Trace test"),
                _("p-value"));
        for (i=0; i<n; i++) {
            trace = gretl_matrix_get(tests, i, 0);
            pv = gamma_LR_T_pval(trace, jcase, n - i, jvar->df);
            pprintf(prn, "%4d%#11.5g [%6.4f]\n", i, trace, pv);
            if (!(opt & OPT_Y)) {
                /* not asymptotic: prefer sample-corrected p-value */
                gretl_matrix_set(pvals, i, 0, pv);
            }
        }
    }

    pputc(prn, '\n');

    if (nexo > 0 || nrexo > 0) {
        if (partial && jcase == J_UNREST_CONST) {
            pputs(prn, _("Warning: the p-values shown are for the "
                  "case of\na restricted trend"));
            pputs(prn, "\n\n");
        } else if (!partial) {
            pputs(prn, _("Note: in general, the test statistics above "
                         "are valid only in the\nabsence of additional "
                         "regressors."));
            pputs(prn, "\n\n");
        }
    }

    record_matrix_test_result(tests, pvals);

    return 0;
}

static int johansen_get_eigenvalues (gretl_matrix *S00,
                                     const gretl_matrix *S01,
                                     const gretl_matrix *S11,
                                     gretl_matrix **M,
                                     gretl_matrix **evals,
                                     int rank)
{
    gretl_matrix *Tmp = NULL;
    int n = S11->cols;
    int err;

    err = gretl_invert_symmetric_matrix(S00);
    if (err) {
        return err;
    }

    Tmp = gretl_matrix_alloc(n, n);
    if (Tmp == NULL) {
        return E_ALLOC;
    }

    *M = gretl_matrix_alloc(n, n);
    if (*M == NULL) {
        gretl_matrix_free(Tmp);
        return E_ALLOC;
    }

    gretl_matrix_qform(S01, GRETL_MOD_TRANSPOSE,
                       S00, Tmp, GRETL_MOD_NONE);

    *evals = gretl_gensymm_eigenvals(Tmp, S11, *M, &err);

    if (!err) {
        err = gretl_symmetric_eigen_sort(*evals, *M, rank);
    }

    gretl_matrix_free(Tmp);

    return err;
}

/* Public entry point for cointegration test */

int johansen_coint_test (GRETL_VAR *jvar, const DATASET *dset,
                         gretlopt opt, PRN *prn)
{
    int p1 = jvar->jinfo->R1->cols;
    int p = jvar->neqns;
    int err = 0;

    jvar->jinfo->Beta = gretl_matrix_alloc(p1, p);
    jvar->jinfo->Alpha = gretl_matrix_alloc(p, p);
    jvar->jinfo->evals = gretl_column_vector_alloc(p);

    if (jvar->jinfo->Beta == NULL ||
        jvar->jinfo->Alpha == NULL ||
        jvar->jinfo->evals == NULL) {
        err = E_ALLOC;
    }

    if (!err) {
        err = gretl_matrix_SVD_johansen_solve(jvar->jinfo->R0,
                                              jvar->jinfo->R1,
                                              jvar->jinfo->evals,
                                              jvar->jinfo->Beta,
                                              jvar->jinfo->Alpha, 0);
    }

    if (err) {
        pputs(prn, _("Failed to find eigenvalues\n"));
    } else {
        johansen_ll_calc(jvar, jvar->jinfo->evals);
        compute_coint_test(jvar, dset, opt, prn);

        if (!(opt & OPT_Q)) {
            print_beta_and_alpha(jvar, jvar->jinfo->evals, p,
                                 dset, prn);
            print_long_run_matrix(jvar, dset, prn);
        }
    }

    return err;
}

static void set_beta_test_df (GRETL_VAR *jvar, const gretl_matrix *H)
{
    int r = jrank(jvar);
    int p1 = jvar->jinfo->Beta->rows;
    int s = H->cols;

    jvar->jinfo->lrdf = r * (p1 - s);
}

/* Likelihood ratio test calculation, for restriction on
   an existing VECM -- which is not modified */

int
johansen_LR_calc (const GRETL_VAR *jvar, const gretl_matrix *evals,
                  const gretl_matrix *H, gretl_restriction *rset,
                  int job, PRN *prn)
{
    gretl_matrix *S00;
    double llr = 0.0;
    double ldet = 0.0;
    double T_2 = (double) jvar->T / 2.0;
    int n = jvar->neqns;
    int r = jrank(jvar);
    int s = H->cols;
    int h = (r > 0)? r : n;
    int i, err = 0;

    S00 = gretl_matrix_copy(jvar->jinfo->S00);

    if (S00 == NULL) {
        err = E_ALLOC;
    } else {
        ldet = gretl_matrix_log_determinant(S00, &err);
    }

    if (!err) {
        llr = - T_2 * n * (1.0 + LN_2_PI) - T_2 * ldet;
        for (i=0; i<h; i++) {
            pprintf(prn, _("eigenvalue %d = %g\n"), i+1, evals->val[i]);
            llr -= T_2 * log(1.0 - evals->val[i]);
        }
        pputc(prn, '\n');
    }

    if (S00 != NULL) {
        gretl_matrix_free(S00);
    }

    if (!err) {
        double x = 2.0 * (jvar->ll - llr);
        int nb = gretl_matrix_rows(jvar->jinfo->Beta);
        int df;

        if (job == V_BETA) {
            df = h * (nb - s);
        } else {
            df = h * (n - s);
        }

        /* allow for possible prior restriction */
        df -= jvar->jinfo->lrdf;

        pprintf(prn, _("Unrestricted loglikelihood (lu) = %.8g\n"), jvar->ll);
        pprintf(prn, _("Restricted loglikelihood (lr) = %.8g\n"), llr);
        pprintf(prn, "2 * (lu - lr) = %g\n", x);

        if (df > 0) {
            double pv = chisq_cdf_comp(df, x);

            if (jvar->jinfo->lrdf > 0) {
                pprintf(prn, _("Allowing for prior restriction, df = %d\n"), df);
            }
            pprintf(prn, "P(%s(%d) > %g) = %g\n", _("Chi-square"), df, x, pv);
            rset_add_results(rset, x, pv, llr);
        }
    }

    return err;
}

#define NSMIN 1.0e-16

static gretl_matrix *johansen_nullspace (const gretl_matrix *R,
					 int *err)
{
    gretl_matrix *H = gretl_matrix_right_nullspace(R, err);

    if (!*err && H->cols == 1) {
	/* this normalization used to be applied in
	   gretl_matrix.c, prior to 2023-05-08
	*/
	double hij, x = 0;
	int i;

	for (i=0; i<H->rows; i++) {
	    if (fabs(H->val[i]) > x) {
		x = H->val[i];
	    }
	}
        for (i=0; i<H->rows; i++) {
            hij = H->val[i] / x;
            H->val[i] = fabs(hij) < NSMIN ? 0 : hij;
        }
    }

    return H;
}

/* see Johansen (1995), section 7.2: compute the p x s
   direct restriction matrix H = R_\perp and form the
   following two intermediate matrices:

   H' * S_{11} * H  in "S11"
   S_{01} * H       in "S01"
*/

static int prep_beta_restriction (const GRETL_VAR *jvar,
                                  const gretl_matrix *R,
                                  gretl_matrix **S01,
                                  gretl_matrix **S11,
                                  gretl_matrix **pH)
{
    gretl_matrix *H;
    int s, p = jvar->neqns;
    int err = 0;

    if (R == NULL) {
        return E_DATA;
    }

    H = johansen_nullspace(R, &err);
    if (err) {
        return err;
    }

    *pH = H;
    s = gretl_matrix_cols(H);

    *S11 = gretl_matrix_alloc(s, s);
    *S01 = gretl_matrix_alloc(p, s);
    if (*S11 == NULL || *S01 == NULL) {
        return E_ALLOC;
    }

    /* calculate S11 <- H' S11 H */
    err = gretl_matrix_qform(H, GRETL_MOD_TRANSPOSE,
                             jvar->jinfo->S11,
                             *S11, GRETL_MOD_NONE);

    if (!err) {
        /* S01 <- S01H */
        err = gretl_matrix_multiply(jvar->jinfo->S01, H, *S01);
    }

    return err;
}

/* check whether the restriction on @jvar in @rset is
   a simple restriction on \beta */

static int simple_beta_restriction (const GRETL_VAR *jvar,
                                    const gretl_restriction *rset)
{
    int ret = 0;

    if (rset_VECM_acols(rset) == 0) {
        const gretl_matrix *R = rset_get_R_matrix(rset);
        const gretl_matrix *q = rset_get_q_matrix(rset);
        int rcols = jvar->neqns + n_restricted_terms(jvar);

        ret = 1;

        if (!gretl_is_zero_matrix(q)) {
            /* non-homogeneous */
            ret = 0;
        } else if (R->cols > rcols) {
            /* not common to all columns */
            ret = 0;
        }
    }

    return ret;
}

/* check whether the restriction on @jvar in @rset is
   a simple restriction on \alpha */

static int simple_alpha_restriction (const GRETL_VAR *jvar,
                                     const gretl_restriction *rset)
{
    int ret = 0;

    if (rset_VECM_bcols(rset) == 0) {
        const gretl_matrix *Ra = rset_get_Ra_matrix(rset);
        const gretl_matrix *qa = rset_get_qa_matrix(rset);

        ret = 1;

        if (!gretl_is_zero_matrix(qa)) {
            /* non-homogeneous */
            ret = 0;
        } else if (Ra->cols > jvar->neqns) {
            /* not common to all columns */
            ret = 0;
        }
    }

    return ret;
}

static int
transcribe_restriction_matrices (const GRETL_VAR *jvar,
                                 const gretl_restriction *rset)
{
    int err = 0;

    if (rset_VECM_bcols(rset) > 0) {
        const gretl_matrix *R = rset_get_R_matrix(rset);
        const gretl_matrix *q = rset_get_q_matrix(rset);

        if (R != jvar->jinfo->R) {
            gretl_matrix_replace(&jvar->jinfo->R, gretl_matrix_copy(R));
        }

        if (q != jvar->jinfo->q) {
            gretl_matrix_replace(&jvar->jinfo->q, gretl_matrix_copy(q));
        }

        if (jvar->jinfo->R == NULL ||
            (q != NULL && jvar->jinfo->q == NULL)) {
            err = E_ALLOC;
        }
    }

    if (!err && rset_VECM_acols(rset) > 0) {
        const gretl_matrix *Ra = rset_get_Ra_matrix(rset);
        const gretl_matrix *qa = rset_get_qa_matrix(rset);

        if (Ra != jvar->jinfo->Ra) {
            gretl_matrix_replace(&jvar->jinfo->Ra, gretl_matrix_copy(Ra));
        }

        if (qa != jvar->jinfo->qa) {
            gretl_matrix_replace(&jvar->jinfo->qa, gretl_matrix_copy(qa));
        }

#if JDEBUG
        gretl_matrix_print(jvar->jinfo->Ra, "jinfo->Ra, transcribed");
        gretl_matrix_print(jvar->jinfo->qa, "jinfo->qa, transcribed");
#endif

        if (jvar->jinfo->Ra == NULL ||
            (qa != NULL && jvar->jinfo->qa == NULL)) {
            err = E_ALLOC;
        }
    }

    return err;
}

/* driver for VECM estimation subject to "general" restrictions
   on beta and/or alpha */

static int j_general_restrict (GRETL_VAR *jvar, gretl_restriction *rset,
                               const DATASET *dset, int flags,
                               PRN *prn)
{
    int err = 0;

#if JDEBUG
    fprintf(stderr, "j_general_restrict: starting\n");
#endif

    err = general_vecm_analysis(jvar, rset, dset, prn);

#if JDEBUG
    fprintf(stderr, " done general_vecm_analysis, err = %d\n", err);
#endif

    if (!err) {
        if (rset_VECM_acols(rset) > 0) {
            flags |= ALPHA_RESTRICTED;
        }
        err = VECM_estimate_full(jvar, rset, dset, flags);
#if JDEBUG
        fprintf(stderr, " VECM_estimate_full, err = %d\n", err);
#endif
    }

    if (!err) {
        err = gretl_VAR_do_error_decomp(jvar->S, jvar->C, NULL);
#if JDEBUG
        fprintf(stderr, " do_error_decomp, err = %d\n", err);
#endif
    }

    if (!err && !bootstrap(flags)) {
        /* FIXME 'k' for AIC etc? */
        err = vecm_ll_stats(jvar);
#if JDEBUG
        fprintf(stderr, " vecm_ll_stats, err = %d\n", err);
#endif
    }

    if (!err && !bootstrap(flags)) {
        err = transcribe_restriction_matrices(jvar, rset);
#if JDEBUG
        fprintf(stderr, " transcribe_restriction_matrices, err = %d\n", err);
#endif
    }

#if JDEBUG
    fprintf(stderr, "j_general_restrict: returning %d\n", err);
#endif

    return err;
}

/* Obtain the unrestricted log-likelihood for running the LR test.  We
   need do this (only) in the context where we're doing full
   estimation of a restricted system.  The prior system, relative to
   which the (new) restriction is defined, may have been restricted
   already, in which case the unrestricted ll is not available for
   comparison.

   This function is low-budget in that we don't bother with the
   eigenvectors, just the eigenvalues.
*/

static int get_unrestricted_ll (GRETL_VAR *jvar)
{
    gretl_matrix *S00 = NULL;
    gretl_matrix *Tmp = NULL;
    gretl_matrix *e = NULL;
    double ldet;
    int n1 = jvar->jinfo->S11->cols;
    int n = jvar->neqns;
    int r = jrank(jvar);
    int i, err = 0;

    S00 = gretl_matrix_copy(jvar->jinfo->S00);
    if (S00 == NULL) {
        return E_ALLOC;
    }

    Tmp = gretl_matrix_alloc(n1, n1);
    if (Tmp == NULL) {
        gretl_matrix_free(S00);
        return E_ALLOC;
    }

    err = gretl_invert_symmetric_matrix(S00);

    if (!err) {
        gretl_matrix_qform(jvar->jinfo->S01, GRETL_MOD_TRANSPOSE,
                           S00, Tmp, GRETL_MOD_NONE);
        e = gretl_gensymm_eigenvals(Tmp, jvar->jinfo->S11, NULL, &err);
    }

    if (!err) {
        gretl_matrix_copy_values(S00, jvar->jinfo->S00);
        ldet = gretl_matrix_log_determinant(S00, &err);
    }

    if (!err) {
        qsort(e->val, n1, sizeof *e->val, gretl_inverse_compare_doubles);
        jvar->jinfo->ll0 = n * (1.0 + LN_2_PI) + ldet;
        for (i=0; i<r; i++) {
            jvar->jinfo->ll0 += log(1.0 - e->val[i]);
        }
        jvar->jinfo->ll0 *= -(jvar->T / 2.0);
    }

    gretl_matrix_free(S00);
    gretl_matrix_free(Tmp);
    gretl_matrix_free(e);

    return err;
}

/* common finalization for estimation subject to simple beta
   restriction, simple alpha restriction, or no restriction.
*/

static int vecm_finalize (GRETL_VAR *jvar, gretl_matrix *H,
                          const gretl_matrix *Ra,
                          const DATASET *dset,
                          int flags)
{
    int do_stderrs = jrank(jvar) < jvar->neqns;
    int err = 0;

    if (Ra == NULL) {
        /* alpha unrestricted */
        err = normalize_beta(jvar, H, &do_stderrs);
    } else {
        do_stderrs = 0;
    }

    if (!err) {
        if (bootstrap(flags)) {
            do_stderrs = 0;
        } else {
            vecm_set_df(jvar, H, Ra);
        }
        err = VECM_estimate_full(jvar, NULL, dset, flags);
    }

    if (!err) {
        err = compute_omega(jvar);
    }

    if (!err && do_stderrs) {
        if (H != NULL) {
            err = restricted_beta_variance(jvar, H);
        } else {
            err = beta_variance(jvar);
        }
    }

    if (!err) {
        err = gretl_VAR_do_error_decomp(jvar->S, jvar->C, NULL);
    }

    if (!err && !bootstrap(flags)) {
        err = vecm_ll_stats(jvar);
    }

    return err;
}

/* estimation subject to "simple" restriction on alpha */

static int
est_simple_alpha_restr (GRETL_VAR *jvar, gretl_restriction *rset,
                        const DATASET *dset, int flags,
                        PRN *prn)
{
    const gretl_matrix *Ra = rset_get_Ra_matrix(rset);
    gretlopt opt = OPT_F;
    int err = 0;

#if JDEBUG
    fprintf(stderr, "\n*** starting est_simple_alpha_restr\n\n");
#endif

    if (bootstrap(flags)) {
        opt |= OPT_B;
    } else {
        err = get_unrestricted_ll(jvar);
    }

    if (!err) {
        err = vecm_alpha_test(jvar, rset, dset, opt, prn);
    }

    if (!err) {
        err = vecm_finalize(jvar, NULL, Ra, dset,
                            flags | ALPHA_RESTRICTED);
    }

    if (!err && !bootstrap(flags)) {
        err = transcribe_restriction_matrices(jvar, rset);
    }

    return err;
}

/* full estimation subject to "simple" restriction on beta */

static int
est_simple_beta_restr (GRETL_VAR *jvar, const gretl_restriction *rset,
                       const DATASET *dset, int flags)
{
    const gretl_matrix *R = NULL;
    gretl_matrix *H = NULL;
    gretl_matrix *M = NULL;
    gretl_matrix *S00 = NULL;
    gretl_matrix *S01 = NULL;
    gretl_matrix *S11 = NULL;
    gretl_matrix *evals = NULL;
    int r = jrank(jvar);
    int err = 0;

#if JDEBUG
    fprintf(stderr, "\n*** starting est_simple_beta_restr\n\n");
#endif

    if (!bootstrap(flags)) {
        err = get_unrestricted_ll(jvar);
    }

    if (!err) {
        R = rset_get_R_matrix(rset);
        err = prep_beta_restriction(jvar, R, &S01, &S11, &H);
    }

    if (!err) {
        S00 = gretl_matrix_copy(jvar->jinfo->S00);
        if (S00 == NULL) {
            err = E_ALLOC;
        }
    }

    if (!err) {
        err = johansen_get_eigenvalues(S00, S01, S11,
                                       &M, &evals, r);
    }

#if JDEBUG
    gretl_matrix_print(M, "raw eigenvector(s)");
#endif

    if (!err) {
        if (jvar->jinfo->Beta == NULL) {
            jvar->jinfo->Beta = gretl_matrix_multiply_new(H, M, &err);
        } else {
            err = gretl_matrix_multiply(H, M, jvar->jinfo->Beta);
        }
    }

    if (!err && !bootstrap(flags)) {
        set_beta_test_df(jvar, H);
        err = johansen_ll_calc(jvar, evals);
    }

    if (!err) {
        err = vecm_finalize(jvar, H, NULL, dset, flags);
    }

    if (!err && !bootstrap(flags)) {
        err = transcribe_restriction_matrices(jvar, rset);
    }

    gretl_matrix_free(H);
    gretl_matrix_free(M);
    gretl_matrix_free(evals);
    gretl_matrix_free(S00);
    gretl_matrix_free(S01);
    gretl_matrix_free(S11);

    return err;
}

/* "unrestricted" VECM estimation */

static int j_estimate_unrestr (GRETL_VAR *jvar,
                               const DATASET *dset)
{
    gretl_matrix *S00 = NULL;
    gretl_matrix *evals = NULL;
    int r = jrank(jvar);
    int err = 0;

#if JDEBUG
    fprintf(stderr, "\n*** starting j_estimate_unrestr\n\n");
#endif

    S00 = gretl_matrix_copy(jvar->jinfo->S00);
    if (S00 == NULL) {
        err = E_ALLOC;
    }

    if (!err) {
        err = johansen_get_eigenvalues(S00,
                                       jvar->jinfo->S01,
                                       jvar->jinfo->S11,
                                       &jvar->jinfo->Beta,
                                       &evals, r);
    }

#if JDEBUG
    gretl_matrix_print(jvar->jinfo->Beta, "raw Beta");
#endif

    if (!err) {
        err = johansen_ll_calc(jvar, evals);
    }

    if (!err) {
        err = vecm_finalize(jvar, NULL, NULL, dset, 0);
    }

    if (!err) {
        if (evals->rows > r) {
            /* include only the non-zero values */
            evals->rows = r;
        }
        jvar->jinfo->evals = evals;
    } else {
        gretl_matrix_free(evals);
    }

    gretl_matrix_free(S00);

    return err;
}

/* Here we prep the system with the initial eigen-analysis, then
   basically hand over to jrestrict.c
*/

static int j_estimate_general (GRETL_VAR *jvar, gretl_restriction *rset,
                               const DATASET *dset, int flags,
                               PRN *prn)
{
    gretl_matrix *S00 = NULL;
    gretl_matrix *evals = NULL;
    gretl_matrix *M = NULL;
    int r = jrank(jvar);
    int err = 0;

#if JDEBUG
    fprintf(stderr, "\n*** starting j_estimate_general\n\n");
#endif

    S00 = gretl_matrix_copy(jvar->jinfo->S00);
    if (S00 == NULL) {
        err = E_ALLOC;
    }

    if (!err) {
        err = johansen_get_eigenvalues(S00,
                                       jvar->jinfo->S01,
                                       jvar->jinfo->S11,
                                       &M, &evals, r);
    }

    if (bootstrap(flags)) {
        if (M != NULL) {
            gretl_matrix_copy_values(jvar->jinfo->Beta, M);
            gretl_matrix_free(M);
        }
    } else {
        jvar->jinfo->Beta = M;
    }

#if JDEBUG
    gretl_matrix_print(jvar->jinfo->Beta, "raw eigenvector(s)");
#endif

    if (!err && !bootstrap(flags)) {
        err = johansen_ll_calc(jvar, evals);
    }

    gretl_matrix_free(S00);
    gretl_matrix_free(evals);

    if (!err) {
#if JDEBUG
        fprintf(stderr, "j_estimate_general: calling j_general_restrict\n");
#endif
        err = j_general_restrict(jvar, rset, dset, flags, prn);
    }

#if JDEBUG
    fprintf(stderr, "j_estimate_general: returning %d\n", err);
#endif

    return err;
}

static int real_johansen_estimate (GRETL_VAR *jvar,
                                   gretl_restriction *rset,
                                   const DATASET *dset,
                                   int flags,
                                   PRN *prn)
{
    int err = 0;

#if JDEBUG
    fprintf(stderr, "johansen_estimate: starting\n");
#endif

    if (rset == NULL) {
        err = j_estimate_unrestr(jvar, dset);
    } else if (simple_beta_restriction(jvar, rset)) {
        err = est_simple_beta_restr(jvar, rset, dset, flags);
    } else if (simple_alpha_restriction(jvar, rset)) {
        err = est_simple_alpha_restr(jvar, rset, dset, flags, prn);
    } else {
        err = j_estimate_general(jvar, rset, dset, flags, prn);
    }

#if JDEBUG
    fprintf(stderr, "johansen_estimate: returning %d\n", err);
#endif

    return err;
}

/* Public entry point for VECM estimation, restricted or not */

int johansen_estimate (GRETL_VAR *jvar, gretl_restriction *rset,
                       const DATASET *dset, PRN *prn)
{
    return real_johansen_estimate(jvar, rset, dset, 0, prn);
}

/* bootstrap round for "unrestricted" VECM: take some
   shortcuts
*/

static int boot_round_simple (GRETL_VAR *jvar, const DATASET *dset)
{
    gretl_matrix *M = NULL;
    gretl_matrix *evals = NULL;
    int err;

    err = johansen_get_eigenvalues(jvar->jinfo->S00, jvar->jinfo->S01,
                                   jvar->jinfo->S11, &M, &evals,
                                   jrank(jvar));

#if JDEBUG
    gretl_matrix_print(M, "raw eigenvector(s)");
#endif

    if (!err) {
        gretl_matrix_copy_values(jvar->jinfo->Beta, M);
        err = normalize_beta(jvar, NULL, NULL);
        if (!err) {
            err = VECM_estimate_full(jvar, NULL, dset,
                                     BOOTSTRAPPING);
        }
        if (!err) {
            err = compute_omega(jvar);
        }
    }

    gretl_matrix_free(M);
    gretl_matrix_free(evals);

    return err;
}

/* Simplified version of the Johansen procedure, to be called in
   the process of computing bootstrap confidence intervals for
   impulse response functions.  We just have to do enough to
   generate the VAR representation.
*/

int johansen_boot_round (GRETL_VAR *jvar, const DATASET *dset)
{
    gretl_restriction *rset;
    int err = 0;

#if JDEBUG
    fprintf(stderr, "\n*** starting johansen_bootstrap_round()\n\n");
#endif

    rset = rset_from_VECM(jvar, &err);
    if (err) {
        return err;
    }

    if (rset != NULL) {
        /* experimental */
        err = real_johansen_estimate(jvar, rset, dset,
                                     BOOTSTRAPPING, NULL);
        free(rset);
    } else {
        err = boot_round_simple(jvar, dset);
    }

    return err;
}

void print_beta_alpha_Pi (const GRETL_VAR *jvar,
                          const DATASET *dset,
                          PRN *prn)
{
    int r = jrank(jvar);

    print_beta_or_alpha(jvar, r, dset, prn, V_BETA, 0);
    print_beta_or_alpha(jvar, r, dset, prn, V_ALPHA, 0);
    pputc(prn, '\n');
    print_long_run_matrix(jvar, dset, prn);
}

/* compute and print beta, alpha and alpha*beta', in the context where
   we've tested a (common, homogeneous) restriction on beta,
   represented by H, and verbose output has been requested.
*/

static int show_beta_alpha_etc (const GRETL_VAR *jvar,
                                const gretl_matrix *H,
                                const gretl_matrix *M,
                                const DATASET *dset,
                                PRN *prn)
{
    JohansenInfo *jv = jvar->jinfo;
    int err = 0;

    gretl_matrix_multiply_mod(H, GRETL_MOD_NONE,
                              M, GRETL_MOD_NONE,
                              jv->Beta, GRETL_MOD_NONE);

    if (jv->rank == 1) {
        /* and if r > 1? */
        double den = jv->Beta->val[0];

        if (!floateq(den, 0.0)) {
            gretl_matrix_divide_by_scalar(jv->Beta, den);
        }
    }

    if (!err) {
        err = compute_alpha(jv);
    }

    if (!err) {
        print_beta_alpha_Pi(jvar, dset, prn);
    }

    return err;
}

/* test for a common, homogeneous restriction on beta (only) */

static int vecm_beta_test (GRETL_VAR *jvar,
                           gretl_restriction *rset,
                           const DATASET *dset,
                           gretlopt opt,
                           PRN *prn)
{
    const gretl_matrix *R;
    gretl_matrix *H = NULL;
    gretl_matrix *M = NULL;
    gretl_matrix *S11 = NULL;
    gretl_matrix *S01 = NULL;
    gretl_matrix *S00 = NULL;
    gretl_matrix *evals = NULL;
    int verbose = (opt & OPT_V);
    int s, p, r;
    int err = 0;

    R = rset_get_R_matrix(rset);
    H = johansen_nullspace(R, &err);

    if (err) {
        return err;
    }

    p = jvar->neqns;
    r = jrank(jvar);
    s = gretl_matrix_cols(H);

    S11 = gretl_matrix_alloc(s, s);
    S01 = gretl_matrix_alloc(p, s);
    S00 = gretl_matrix_copy(jvar->jinfo->S00);

    if (S11 == NULL || S01 == NULL || S00 == NULL) {
        err = E_ALLOC;
        goto bailout;
    }

    pprintf(prn, "\n%s\n\n",
            _("Test of restrictions on cointegrating relations"));

    if (verbose) {
        gretl_matrix_print_to_prn(H, "Restriction matrix, H", prn);
    }

    /* calculate S11 <- H' S11 H */
    err = gretl_matrix_qform(H, GRETL_MOD_TRANSPOSE,
                             jvar->jinfo->S11, S11,
                             GRETL_MOD_NONE);

    if (verbose) {
        gretl_matrix_print_to_prn(S11, "H'*S11*H", prn);
    }

    if (!err) {
        /* S01 <- S01*H */
        err = gretl_matrix_multiply(jvar->jinfo->S01, H, S01);
    }

    if (verbose) {
        gretl_matrix_print_to_prn(S01, "S01*H", prn);
    }

    if (!err) {
        err = johansen_get_eigenvalues(S00, S01, S11,
                                       &M, &evals, r);
    }

    if (!err) {
        if (verbose) {
            gretl_matrix_print_to_prn(M, "M", prn);
        }
        johansen_LR_calc(jvar, evals, H, rset, V_BETA, prn);
    }

    if (!err && verbose) {
        show_beta_alpha_etc(jvar, H, M, dset, prn);
    }

 bailout:

    gretl_matrix_free(H);
    gretl_matrix_free(M);
    gretl_matrix_free(evals);
    gretl_matrix_free(S00);
    gretl_matrix_free(S11);
    gretl_matrix_free(S01);

    return err;
}

/* Test of linear restrictions on the cointegrating relations in a
   VECM.  If the restrictions are "simple" (homogeneous and in common)
   we do the test using the eigen-system approach.  If they are
   "general" restrictions we hand off to the specialized machinery in
   jrestrict.c.
*/

int vecm_test_restriction (GRETL_VAR *jvar,
                           gretl_restriction *rset,
                           const DATASET *dset,
                           gretlopt opt,
                           PRN *prn)
{
    gretl_matrix *B0 = NULL;
    gretl_matrix *A0 = NULL;
    PRN *vprn;
    int err = 0;

    B0 = gretl_matrix_copy(jvar->jinfo->Beta);
    A0 = gretl_matrix_copy(jvar->jinfo->Alpha);

    if (B0 == NULL || A0 == NULL) {
        return E_ALLOC;
    }

    vprn = (opt & OPT_S)? NULL : prn;

    if (simple_beta_restriction(jvar, rset)) {
        err = vecm_beta_test(jvar, rset, dset, opt, vprn);
    } else if (simple_alpha_restriction(jvar, rset)) {
        err = vecm_alpha_test(jvar, rset, dset, opt, vprn);
    } else {
        err = general_vecm_analysis(jvar, rset, dset, vprn);
    }

    if (!err) {
        rset_record_LR_result(rset);
    }

    /* restore orginal Beta, Alpha on exit */

    gretl_matrix_replace(&jvar->jinfo->Beta, B0);
    gretl_matrix_replace(&jvar->jinfo->Alpha, A0);

    return err;
}
