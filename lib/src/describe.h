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

typedef struct {
    char varname[9];              /* for ID purposes */
    int numbins;                  /* number of bins or intervals */
    double xbar, sdx;             /* mean and std dev of variable */
    double *midpt, *endpt;        /* arrays of midpoints and endpoints
                                     of the intervals */
    int *f;                       /* frequencies in the intervals */
    double chisqu;                /* Chi-squared statistic for testing
                                     for a Gaussian distribution */
    int n;
    int t1, t2;
    int errcode;                  
    char errmsg[ERRLEN];
} FREQDIST;

/* functions follow */

void free_freq (FREQDIST *freq);

FREQDIST *freq_func (double **pZ, const DATAINFO *pdinfo, 
		     double *zz, const int nzz, 
		     const char *varname, const int params);

int corrgram (const int *list, const int order, 
	      double **pZ, DATAINFO *pdinfo, 
	      const PATHS *ppaths, const int batch, 
	      print_t *prn);

int periodogram (const int *list, 
		 double **pZ, const DATAINFO *pdinfo, 
		 const PATHS *ppaths, const int batch, 
		 const int opt, print_t *prn);

GRETLSUMMARY *summary (int *list, 
		       double **pZ, const DATAINFO *pdinfo,
		       print_t *prn);

void print_summary (GRETLSUMMARY *summ,
		    const DATAINFO *pdinfo,
		    print_t *prn, int batch); 

void free_summary (GRETLSUMMARY *summ);

CORRMAT *corrlist (int *list, 
		   double **pZ, const DATAINFO *pdinfo);

void free_corrmat (CORRMAT *corrmat);

int esl_corrmx (int *list, 
		double **pZ, const DATAINFO *pdinfo, 
		const int batch, print_t *prn);

int means_test (int *list, 
		double *Z, const DATAINFO *pdinfo, 
		const int vareq, print_t *prn);

int vars_test (int *list, 
	       double *Z, const DATAINFO *pdinfo, 
	       print_t *prn);

void matrix_print_corr (CORRMAT *corr, const DATAINFO *pdinfo,
			const int batch, print_t *prn);
