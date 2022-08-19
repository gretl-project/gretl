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

/* below: LOOP_PRINT, LOOP_MODEL and LOOP_STORE are
   used only in "progressive" loops.
*/

typedef struct {
    int lineno;    /* location: line number in loop */
    int n;         /* number of repetitions */
    int nvars;     /* number of variables */
    char **names;  /* names of vars to print */
    bigval *sum;   /* running sum of values */
    bigval *ssq;   /* running sum of squares */
    double *xbak;  /* previous values */
    int *diff;     /* indicator for difference */
    char *na;      /* indicator for NAs in calculation */
} LOOP_PRINT;

typedef struct {
    int lineno;             /* location: line number in loop */
    int n;                  /* number of repetitions */
    int nc;                 /* number of coefficients */
    MODEL *model0;          /* copy of initial model */
    bigval *bigarray;       /* global pointer array */
    bigval *sum_coeff;      /* sums of coefficient estimates */
    bigval *ssq_coeff;      /* sums of squares of coeff estimates */
    bigval *sum_sderr;      /* sums of estimated std. errors */
    bigval *ssq_sderr;      /* sums of squares of estd std. errs */
    double *cbak;           /* previous values of coeffs */
    double *sbak;           /* previous values of std. errs */
    int *cdiff;             /* indicator for difference in coeff */
    int *sdiff;             /* indicator for difference in s.e. */
} LOOP_MODEL;

typedef struct {
    int lineno;     /* location: line number in loop */
    int n;          /* number of observations */
    int nvars;      /* number of variables to store */
    char **names;   /* names of vars to print */
    char *fname;    /* filename for output */
    gretlopt opt;   /* formatting option */
    DATASET *dset;  /* temporary data storage */
} LOOP_STORE;

static int extend_loop_dataset (LOOP_STORE *lstore);
static void loop_model_free (LOOP_MODEL *lmod);
static void loop_print_free (LOOP_PRINT *lprn);
static void loop_store_free (LOOP_STORE *lstore);
static void loop_store_init (LOOP_STORE *lstore);
static int loop_print_save_model (MODEL *pmod, DATASET *dset,
                                  PRN *prn, ExecState *s);
