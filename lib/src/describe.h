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
} FREQDIST;

/* functions follow */

void free_freq (FREQDIST *freq);

FREQDIST *freqdist (double ***pZ, const DATAINFO *pdinfo, 
		    int varno, int params);

int corrgram (const int varno, const int order, 
	      double ***pZ, DATAINFO *pdinfo, 
	      PATHS *ppaths, const int batch, 
	      PRN *prn);

int periodogram (const int varno, 
		 double ***pZ, const DATAINFO *pdinfo, 
		 PATHS *ppaths, const int batch, 
		 const int opt, PRN *prn);

GRETLSUMMARY *summary (LIST list, 
		       double ***pZ, const DATAINFO *pdinfo,
		       PRN *prn);

void print_summary (GRETLSUMMARY *summ,
		    const DATAINFO *pdinfo,
		    const int pause, PRN *prn); 

void free_summary (GRETLSUMMARY *summ);

CORRMAT *corrlist (LIST list, 
		   double ***pZ, const DATAINFO *pdinfo);

void free_corrmat (CORRMAT *corrmat);

int esl_corrmx (LIST list, 
		double ***pZ, const DATAINFO *pdinfo, 
		const int pause, PRN *prn);

int means_test (LIST list, 
		double **Z, const DATAINFO *pdinfo, 
		const int vareq, PRN *prn);

int vars_test (LIST list, 
	       double **Z, const DATAINFO *pdinfo, 
	       PRN *prn);

void matrix_print_corr (CORRMAT *corr, const DATAINFO *pdinfo,
			const int pause, PRN *prn);
