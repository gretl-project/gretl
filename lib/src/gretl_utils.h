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

#ifndef GRETL_UTILS_H
#define GRETL_UTILS_H

#include "libgretl.h"

#include <float.h>
#include <limits.h>

#define floateq(x, y)  (fabs((x) - (y)) < DBL_EPSILON)
#define floatneq(x, y) (fabs((x) - (y)) > DBL_EPSILON)
#define floatgt(x, y)  ((x) - (y) > DBL_EPSILON)
#define floatlt(x, y)  ((y) - (x) > DBL_EPSILON)

#define ok_int(x) (x <= (double) INT_MAX && x >= (double) INT_MIN)

enum {
    SESSION_CLEAR_ALL,
    SESSION_CLEAR_DATASET
};

enum {
    BLAS_UNKNOWN,
    BLAS_NETLIB,
    BLAS_ATLAS,
    BLAS_OPENBLAS,
    BLAS_MKL,
    BLAS_VECLIB,
    BLAS_BLIS
};

void libgretl_init (void);

#ifdef HAVE_MPI
int libgretl_mpi_init (int self, int np, int dcmt);
#else
/* dummy function */
int gretl_mpi_initialized (void);
#endif

int blas_is_threaded (void);

void blas_set_num_threads (int nt);

int blas_get_num_threads (void);

int auto_mpi_ok (void);

void libgretl_session_cleanup (int mode);

void libgretl_cleanup (void);

double date_as_double (int t, int pd, double sd0);

/* checks on variables */

int gretl_isdummy (int t1, int t2, const double *x);

int gretl_iszero (int t1, int t2, const double *x);

int gretl_isconst (int t1, int t2, const double *x);

int gretl_isstoch (int t1, int t2, const double *x);

int gretl_isunits (int t1, int t2, const double *x);

int gretl_iscount (int t1, int t2, const double *x);

int gretl_isint (int t1, int t2, const double *x);

int gretl_isdiscrete (int t1, int t2, const double *x);

int gretl_ispositive (int t1, int t2, const double *x, int strict);

int gretl_is_oprobit_ok (int t1, int t2, const double *x);

int accept_as_discrete (const DATASET *dset, int v, int strict);

int true_const (int v, const DATASET *dset);

/* setting observations */

char *format_obs (char *obs, int maj, int min, int pd);

int set_obs (const char *parm1, const char *parm2,
             DATASET *dset, gretlopt opt);

int simple_set_obs (DATASET *dset, int pd, const char *stobs,
                    gretlopt opt);

/* sorting and comparison */

int gretl_compare_doubles (const void *a, const void *b);

int gretl_inverse_compare_doubles (const void *a, const void *b);

int count_distinct_values (const double *x, int n);

int count_distinct_int_values (const int *x, int n);

int rearrange_id_array (double *x, int m, int n);

int gretl_compare_ints (const void *a, const void *b);

int gretl_compare_strings (const void *a, const void *b);

/* miscellaneous */

int libgretl_version (void);

int gretl_version_number (const char *version);

char *gretl_version_string (char *targ, int vnum);

void printlist (const int *list, const char *msg);

double gretl_double_from_string (const char *s, int *err);

int gretl_int_from_string (const char *s, int *err);

int positive_int_from_string (const char *s);

int varnum_from_string (const char *str, DATASET *dset);

int gretl_int_from_double (double x, int *err);

guint32 gretl_unsigned_from_double (double x, int *err);

guint64 gretl_uint53_from_double (double x, int *err);

GretlType gretl_type_from_name (const char *s, const DATASET *dset);

double *copyvec (const double *src, int n);

void doubles_array_free (double **X, int m);

double **doubles_array_new (int m, int n);

double **doubles_array_new0 (int m, int n);

int doubles_array_adjust_length (double **X, int m, int new_n);

double **data_array_from_model (const MODEL *pmod, double **Z,
                                int missv);

int ijton (int i, int j, int nrows);

int transcribe_array (double *targ, const double *src,
                      const DATASET *dset);

int gretl_copy_file (const char *src, const char *dest);

int gretl_delete_var_by_name (const char *s, PRN *prn);

int combine_categories (DATASET *dset, int v0, int v1, int v2);

#ifndef WIN32

int gretl_spawn (char *cmdline);

int gnuplot_make_image (const char *input_fname);

int gretl_pipe_output (gchar **argv, gchar **envp,
                       const char *currdir, PRN *prn,
                       gchar **errp);

int get_stack_size (void);

#endif

/* model selection criteria */

int gretl_calculate_criteria (double ess, int n, int k,
                              double *ll, double *aic, double *bic,
                              double *hqc);

int ls_criteria (MODEL *pmod);

/* hypothesis tests mechanism */

int get_last_test_type (void);

void record_test_result (double teststat, double pval);

void record_matrix_test_result (gretl_matrix *tests,
                                gretl_matrix *pvals);

void record_LR_test_result (double teststat, double pval, double lnl);

void record_QLR_test_result (double teststat, double pval, double brk);

double get_last_test_statistic (void);

double get_last_pvalue (void);

double get_last_lnl (void);

double get_last_break (void);

gretl_matrix *get_last_test_matrix (int *err);

gretl_matrix *get_last_pvals_matrix (int *err);

/* timer */

double gretl_stopwatch (void);

#if 0 /* see comments in gretl_utils.c */

/* aligned allocation */

void *gretl_aligned_malloc (size_t size, size_t alignment);

void gretl_aligned_free (void *mem);

#endif

/* search in path, etc. */

int check_for_program (const char *prog);

const char *blas_variant_string (void);

int blas_is_openblas (void);

int get_blas_details (char **s1, char **s2, char **s3);

char *get_cpu_details (void);

int avx_support (void);

gint64 gretl_monotonic_time (void);

/* gzip, gunzip */

int gretl_gzip (char *fname, char *zname);

int gretl_gunzip (char *zname, char *fname);

/* binary conversions */

gretl_matrix *dec2bin (double x, const gretl_matrix *v, int *err);

gretl_matrix *bin2dec (const gretl_matrix *m, int *err);

#endif /* GRETL_UTILS_H */

#ifdef __ARM_ARCH_ISA_A64

double _Complex __divdc3 (double a, double b, double c, double d);

#endif
