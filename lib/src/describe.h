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

#ifndef DESCRIBE_H
#define DESCRIBE_H

typedef struct Summary_ Summary;
typedef struct FreqDist_ FreqDist;
typedef struct Xtab_ Xtab;
typedef struct MahalDist_ MahalDist;

struct Summary_ {
    int n;
    int missing;
    int *list;
    double *stats;
    double *mean;
    double *median;
    double *sd;
    double *skew; 
    double *xkurt;
    double *low;
    double *high;
    double *cv;
};

struct FreqDist_ {
    char varname[VNAMELEN];  /* for ID purposes */
    int discrete;            /* 1 if variable contains integers */
    int dist;                /* code for theoretical distribution */
    int numbins;             /* number of bins or intervals */
    double xbar, sdx;        /* mean and std dev of variable */
    double *midpt;           /* array of midpoints of intervals */
    double *endpt;           /* array of endpoints of intervals */
    int *f;                  /* frequencies in the intervals */
    double test;             /* either Chi-squared statistic for testing
                                for a Gaussian distribution, or z statistic
			        for testing for Gamma dist. */
    int n;
    int t1, t2;
};

struct Xtab_ {
    char rvarname[VNAMELEN]; 
    char cvarname[VNAMELEN]; 
    int rows, cols;
    double *rval, *cval;
    double *rtotal, *ctotal;
    int **f;
    int n, missing;
    int t1, t2;
};

/* functions follow */

int eval_ytest (double y, GretlOp op, double test);

int gretl_minmax (int t1, int t2, const double *x, 
		  double *min, double *max);

double gretl_min (int t1, int t2, const double *x);

double gretl_max (int t1, int t2, const double *x);

double gretl_sum (int t1, int t2, const double *x);

double gretl_mean (int t1, int t2, const double *x);

double gretl_restricted_mean (int t1, int t2, const double *x,
			      const double *y, GretlOp yop, 
			      double yval);

double gretl_median (int t1, int t2, const double *x);

double gretl_sst (int t1, int t2, const double *x);

double gretl_variance (int t1, int t2, const double *x);

double gretl_restricted_variance (int t1, int t2, const double *x,
				  const double *y, GretlOp yop,
				  double yval);

double gretl_stddev (int t1, int t2, const double *x);

double gretl_restricted_stddev (int t1, int t2, const double *x,
				const double *y, GretlOp yop,
				double yval);

double gretl_long_run_variance (int t1, int t2, const double *x, int m);

double gretl_covar (int t1, int t2, const double *x, const double *y);

double gretl_corr (int t1, int t2, const double *x, const double *y,
		   int *missing);

double gretl_corr_rsq (int t1, int t2, const double *x, const double *y);

double gretl_acf (int k, int t1, int t2, const double *y);

double gretl_xcf (int k, int t1, int t2, const double *x, const double *y);

int gretl_moments (int t1, int t2, const double *x, 
		   double *xbar, double *sd, 
		   double *skew, double *kurt, int k);

void free_freq (FreqDist *freq);

int freq_setup (int v, const double **Z, const DATAINFO *pdinfo,
		int *pn, double *pxmax, double *pxmin, int *nbins, 
		double *binwidth);

FreqDist *get_freq (int varno, const double **Z, const DATAINFO *pdinfo, 
		    int nbins, int params, gretlopt opt, int *err);

int freqdist (int varno, const double **Z, const DATAINFO *pdinfo,
	      int graph, gretlopt opt, PRN *prn);

int crosstab (const int *list, const double **Z, 
	      const DATAINFO *pdinfo, gretlopt opt,
	      PRN *prn);

int model_error_dist (const MODEL *pmod, double ***pZ,
		      DATAINFO *pdinfo, PRN *prn);

int auto_acf_order (int pd, int nobs);

int auto_spectrum_order (int T, gretlopt opt);

int corrgram (int varno, int order, int nparam,
	      double ***pZ, DATAINFO *pdinfo, 
	      PRN *prn, gretlopt opt);

int xcorrgram (const int *list, int order, 
	      double ***pZ, DATAINFO *pdinfo, 
	      PRN *prn, gretlopt opt);

int periodogram (int varno, int width, 
		 double ***pZ, const DATAINFO *pdinfo, 
		 gretlopt opt, PRN *prn);

Summary *summary (const int *list, const double **Z, 
		  const DATAINFO *pdinfo,
		  PRN *prn);

void print_summary (const Summary *summ,
		    const DATAINFO *pdinfo,
		    PRN *prn); 

void free_summary (Summary *summ);

VMatrix *corrlist (int *list, const double **Z, const DATAINFO *pdinfo);

VMatrix *vmatrix_new (void);

void free_vmatrix (VMatrix *vmat);

int gretl_corrmx (int *list, const double **Z, const DATAINFO *pdinfo, 
		  PRN *prn);

int means_test (const int *list, const double **Z, 
		const DATAINFO *pdinfo, 
		gretlopt opt, PRN *prn);

int vars_test (const int *list, const double **Z, 
	       const DATAINFO *pdinfo, PRN *prn);

void matrix_print_corr (VMatrix *corr, const DATAINFO *pdinfo, PRN *prn);

double dh_root_b1_to_z1 (double rb1, double n);

double dh_b2_to_z2 (double b1, double b2, double n);

double doornik_chisq (double skew, double xkurt, int n);

int gretl_system_normality_test (const gretl_matrix *E, 
				 const gretl_matrix *Sigma, 
				 PRN *prn);

int mahalanobis_distance (const int *list, double ***pZ,
			  DATAINFO *pdinfo, gretlopt opt,
			  PRN *prn);

MahalDist *get_mahal_distances (const int *list, double ***pZ,
				DATAINFO *pdinfo, gretlopt opt,
				PRN *prn);

void free_mahal_dist (MahalDist *md);

const double *mahal_dist_get_distances (const MahalDist *md);

int mahal_dist_get_n (const MahalDist *md);

const int *mahal_dist_get_varlist(const MahalDist *md);

double gretl_gini (int t1, int t2, const double *x);

int gini (int vnum, const double **Z, DATAINFO *pdinfo, 
	  gretlopt opt, PRN *prn);

#endif /* DESCRIBE_H */
