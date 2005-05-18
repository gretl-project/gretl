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

typedef struct Summary_ Summary;
typedef struct FREQDIST_ FREQDIST;

struct Summary_ {
    int n;
    int *list;
    double *mean;
    double *median;
    double *sd;
    double *skew; 
    double *xkurt;
    double *low;
    double *high;
    double *cv;
};

struct FREQDIST_ {
    char varname[VNAMELEN];  /* for ID purposes */
    int dist;                /* code for theoretical distribution */
    int numbins;             /* number of bins or intervals */
    double xbar, sdx;        /* mean and std dev of variable */
    double *midpt;           /* array of midpoints of intervals */
    double *endpt;           /* array of endpoints of intervals */
    int *f;                  /* frequencies in the intervals */
    double test;             /* either Chi-squared statistic for testing
                                for a Gaussian distribution, or z statistic
			        for testing for Gamma dist.
			     */
    int n;
    int t1, t2;
};

/* functions follow */

int gretl_minmax (int t1, int t2, const double *x, 
		  double *min, double *max);

double gretl_mean (int t1, int t2, const double *x);

double gretl_median (int t1, int t2, const double *x);

double gretl_sst (int t1, int t2, const double *x);

double gretl_variance (int t1, int t2, const double *x);

double gretl_stddev (int t1, int t2, const double *x);

double gretl_covar (int t1, int t2, const double *x, const double *y);

double gretl_corr (int t1, int t2, const double *x, const double *y);

double gretl_corr_rsq (int t1, int t2, const double *x, const double *y);

int gretl_moments (int t1, int t2, const double *x, 
		   double *xbar, double *sd, 
		   double *skew, double *kurt, int k);

void free_freq (FREQDIST *freq);

FREQDIST *get_freq (int varno, const double **Z, const DATAINFO *pdinfo, 
		    int params, gretlopt opt);

int freqdist (int varno, const double **Z, const DATAINFO *pdinfo,
	      int graph, PRN *prn, gretlopt opt);

int model_error_dist (const MODEL *pmod, double ***pZ,
		      DATAINFO *pdinfo, PRN *prn);

int auto_acf_order (int pd, int nobs);

int corrgram (int varno, int order, 
	      double ***pZ, DATAINFO *pdinfo, 
	      int batch, PRN *prn);

int periodogram (int varno, 
		 double ***pZ, const DATAINFO *pdinfo, 
		 int batch, int opt, PRN *prn);

Summary *summary (const int *list, const double **Z, 
		       const DATAINFO *pdinfo,
		       PRN *prn);

void print_summary (const Summary *summ,
		    const DATAINFO *pdinfo,
		    PRN *prn); 

void free_summary (Summary *summ);

CorrMat *corrlist (const int *list, const double **Z, const DATAINFO *pdinfo);

void free_corrmat (CorrMat *corrmat);

int gretl_corrmx (int *list, const double **Z, const DATAINFO *pdinfo, 
		  PRN *prn);

int means_test (const int *list, const double **Z, 
		const DATAINFO *pdinfo, 
		gretlopt opt, PRN *prn);

int vars_test (const int *list, const double **Z, 
	       const DATAINFO *pdinfo, PRN *prn);

void matrix_print_corr (CorrMat *corr, const DATAINFO *pdinfo,
			PRN *prn);

double doornik_chisq (double skew, double kurt, int n);

int mahalanobis_distance (const int *list, double ***pZ,
			  DATAINFO *pdinfo, gretlopt opt,
			  PRN *prn);

char *unique_savename (char *vname, DATAINFO *pdinfo, int vmax);

