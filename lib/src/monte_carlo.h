/* gretl - The Gnu Regression, Econometrics and Time-series Library
 * Copyright (C) 1999-2000 Ramu Ramanathan and Allin Cottrell
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License 
 * along with this software; if not, write to the 
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

/* monte_carlo.h for gretl */

#if defined(ENABLE_GMP)
#include <gmp.h>
typedef mpf_t bigval;
#elif defined(HAVE_LONG_DOUBLE)
typedef long double bigval;
#else
typedef double bigval;
#endif

enum loop_types {
    COUNT_LOOP,
    WHILE_LOOP,
    FOR_LOOP
};

enum ifcodes {
    SET_FALSE,
    SET_TRUE,
    SET_ELSE,
    SET_ENDIF,
    IS_FALSE,
    RELAX
};

typedef struct {
    int ID;
    int *list;
    bigval *sum;
    bigval *ssq;
} LOOP_PRINT;   

typedef struct {
    int ID;                      /* ID number for model */
    int ci;                      /* command index for model */
    int t1, t2, nobs;            /* starting observation, ending
                                    observation, and number of obs */
    int ncoeff, dfn, dfd;        /* number of coefficents; degrees of
                                    freedom in numerator and denominator */
    int *list;                   /* list of variables by ID number */
    int ifc;                     /* = 1 if the equation includes a constant,
                                    else = 0 */
    bigval *sum_coeff;      /* sums of coefficient estimates */
    bigval *ssq_coeff;      /* sums of squares of coeff estimates */
    bigval *sum_sderr;      /* sums of estimated std. errors */
    bigval *ssq_sderr;      /* sums of squares of estd std. errs */
} LOOP_MODEL;

typedef struct {
    char type;
    int err;
    int ntimes;
    int lvar, rvar;
    double rval;
    int ineq;
    int ncmds;
    int nmod;
    int nprn;
    int nstore;
    int next_model;
    int next_print;
    char **lines;
    int *ci;
    MODEL **models;
    LOOP_MODEL *lmodels;
    LOOP_PRINT *prns;
    char **storename;
    char **storelbl;
    double *storeval;
} LOOPSET;

/* functions follow */

int ok_in_loop (int ci, const LOOPSET *ploop);

int parse_loopline (char *line, LOOPSET *ploop, 
		    DATAINFO *pdinfo, const double **Z);

int loop_condition (int k, LOOPSET *ploop, 
		    double **Z, DATAINFO *pdinfo); 

void monte_carlo_free (LOOPSET *ploop);

int loop_model_init (LOOP_MODEL *plmod, const MODEL *pmod,
		     int id);

int loop_print_init (LOOP_PRINT *pprn, const LIST list, int id);

int loop_store_init (LOOPSET *ploop, const LIST list, 
		     DATAINFO *pdinfo);

int update_loop_model (LOOPSET *ploop, int cmdnum, MODEL *pmod);

int update_loop_print (LOOPSET *ploop, int cmdnum, 
		       const LIST list, double ***pZ, 
		       const DATAINFO *pdinfo);

void print_loop_results (LOOPSET *ploop, 
			 const DATAINFO *pdinfo, 
			 PRN *prn, PATHS *ppaths,
			 char *loopstorefile);

int add_to_loop (LOOPSET *ploop, char *line, int ci, 
		 gretlopt oflags);

void get_cmd_ci (const char *line, CMD *command);

int get_modnum_by_cmdnum (LOOPSET *ploop, int cmdnum);

int if_eval (const char *line, double ***pZ, DATAINFO *pdinfo);

int ifstate (int code);
