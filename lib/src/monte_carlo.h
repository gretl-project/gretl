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

typedef enum {
    COUNT_LOOP,
    WHILE_LOOP,
    FOR_LOOP
} loop_types;

typedef struct {
    int ID;
    int *list;
    long double *sum;
    long double *ssq;
} LOOP_PRINT;   

typedef struct {
    int ID;                      /* ID number for model */
    int t1, t2, nobs;            /* starting observation, ending
                                    observation, and number of obs */
    int ncoeff, dfn, dfd;        /* number of coefficents; degrees of
                                    freedom in numerator and denominator */
    int *list;                   /* list of variables by ID number */
    int ifc;                     /* = 1 if the equation includes a constant,
                                    else = 0 */
    long double *sum_coeff;      /* sums of coefficient estimates */
    long double *ssq_coeff;      /* sums of squares of coeff estimates */
    long double *sum_sderr;      /* sums of estimated std. errors */
    long double *ssq_sderr;      /* sums of squares of estd std. errs */
} LOOP_MODEL;

typedef struct {
    char type;
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

int ok_in_loop (int ci);

int parse_loopline (char *line, LOOPSET *ploop, 
		    DATAINFO *pdinfo);

int loop_condition (int k, LOOPSET *ploop, 
		    double **Z, DATAINFO *pdinfo); 

void monte_carlo_free (LOOPSET *ploop);

int loop_model_init (LOOP_MODEL *plmod, const MODEL *pmod,
		     const int id);

int loop_print_init (LOOP_PRINT *pprn, const LIST list, const int id);

int loop_store_init (LOOPSET *ploop, const LIST list, 
		     DATAINFO *pdinfo);

int update_loop_model (LOOPSET *ploop, const int cmdnum, MODEL *pmod);

int update_loop_print (LOOPSET *ploop, const int cmdnum, 
		       const LIST list, double ***pZ, 
		       const DATAINFO *pdinfo);

void print_loop_results (LOOPSET *ploop, 
			 const DATAINFO *pdinfo, 
			 PRN *prn, PATHS *ppaths, int *model_count,
			 char *loopstorefile);

int add_to_loop (LOOPSET *ploop, char *line, const int ci,
		 const int opt);

void get_cmd_ci (const char *line, CMD *command);

int get_modnum_by_cmdnum (LOOPSET *ploop, const int cmdnum);

int if_eval (const char *line, double **Z, DATAINFO *pdinfo);
