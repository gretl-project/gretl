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

#ifndef DESCRIBE_H
#define DESCRIBE_H

typedef struct MahalDist_ MahalDist;

typedef struct Summary_ {
    gretlopt opt;
    int n;
    int weight_var;
    int *misscount;
    int *list;
    int *ival;
    double *stats;
    double *mean;
    double *median;
    double *sd;
    double *skew;
    double *xkurt;
    double *low;
    double *high;
    double *cv;
    double *perc05;
    double *perc95;
    double *iqr;
    double sw;
    double sb;
} Summary;

typedef struct FreqDist_ {
    char varname[VNAMELEN];  /* for ID purposes */
    char gname[MAXDISP];     /* name to use in plots */
    int discrete;            /* 1 if series contains integers */
    int strvals;             /* 1 if series is string-valued */
    int dist;                /* code for theoretical distribution */
    int numbins;             /* number of bins or intervals */
    double xbar, sdx;        /* mean and std dev of variable */
    double *midpt;           /* array of midpoints of intervals */
    double *endpt;           /* array of endpoints of intervals */
    char **S;                /* array of string values */
    int *f;                  /* frequencies */
    double test;             /* Chi-squared statistic for testing for a
                                Gaussian distribution, or z statistic
			        for testing for Gamma dist. */
    int n;                   /* included observation */
    int t1, t2;              /* sample limits */
} FreqDist;

typedef struct Xtab_ {
    char rvarname[VNAMELEN]; /* name of rows series */
    char cvarname[VNAMELEN]; /* name of cols series */
    char **Sr;               /* array of row string values */
    char **Sc;               /* array of column string values */
    int rows, cols;          /* dimensions of table */
    double *rval, *cval;     /* row and column numeric values */
    int *rtotal, *ctotal;    /* marginal totals */
    int **f;                 /* array of frequencies */
    int n, missing;          /* observation counts */
    int t1, t2;              /* sample limits */
    int rstrs;               /* row string-valued flag */
    int cstrs;               /* col string-valued flag */
} Xtab;

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

double gretl_quantile (int t1, int t2, const double *x,
		       double p, gretlopt opt, int *err);

int gretl_array_quantiles (double *a, int n, double *p, int k);

double gretl_array_quantile (double *a, int n, double p);

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

double gretl_long_run_variance (int t1, int t2, const double *x,
				int m, double mu);

double gretl_covar (int t1, int t2, const double *x, const double *y,
		    int *missing);

double gretl_corr (int t1, int t2, const double *x, const double *y,
		   int *missing);

double gretl_corr_rsq (int t1, int t2, const double *x, const double *y);

double gretl_skewness (int t1, int t2, const double *x);

double gretl_kurtosis (int t1, int t2, const double *x);

int gretl_moments (int t1, int t2, const double *x,
		   const double *wts,
		   double *xbar, double *sd,
		   double *skew, double *kurt, int k);

double *gretl_sorted_series (int v, const DATASET *dset,
			     gretlopt opt, int *n,
			     int *err);

void free_freq (FreqDist *freq);

int freq_setup (int v, const DATASET *dset, int *pn,
		double *pxmax, double *pxmin, int *nbins,
		double *binwidth);

FreqDist *get_freq (int varno, const DATASET *dset,
		    double fmin, double fwid, int nbins, int params,
		    gretlopt opt, int *err);

FreqDist *get_discrete_freq (int v, const DATASET *dset,
			     gretlopt opt, int *err);

int freqdist (int varno, const DATASET *dset,
	      gretlopt opt, PRN *prn);

gretl_matrix *freqdist_matrix (const double *x, int t1, int t2,
			       int *err);

void *get_last_result_data (GretlType *type, int *err);

void set_last_result_data (void *ptr, GretlType type);

void last_result_cleanup (void);

int crosstab (const int *list, const DATASET *dset,
	      gretlopt opt, PRN *prn);

int crosstab_from_matrix (gretlopt opt, PRN *prn);

int compare_xtab_rows (const void *a, const void *b);

Xtab *single_crosstab (const int *list, const DATASET *dset,
		       gretlopt opt, PRN *prn, int *err);

gretl_matrix *xtab_to_matrix (const Xtab *tab);

void free_xtab (Xtab *tab);

int correspondence (const double *x, const double *y,
		    int n, int *err);

int model_error_dist (const MODEL *pmod, DATASET *dset,
		      gretlopt opt, PRN *prn);

int auto_acf_order (int T);

int auto_spectrum_order (int T, gretlopt opt);

int corrgram (int varno, int order, int nparam,
	      DATASET *dset, gretlopt opt, PRN *prn);

int xcorrgram (const int *list, int order,
	       DATASET *dset, gretlopt opt, PRN *prn);

int periodogram (int varno, int width,
		 const DATASET *dset,
		 gretlopt opt, PRN *prn);

int residual_periodogram (const double *x, int width,
			  const DATASET *dset,
			  gretlopt opt, PRN *prn);

gretl_matrix *periodogram_matrix (const double *x, int t1, int t2,
				  int width, int *err);

int fractint (int varno, int order,
	      const DATASET *dset,
	      gretlopt opt, PRN *prn);

Summary *get_summary (const int *list, const DATASET *dset,
		      gretlopt opt, PRN *prn,
		      int *err);

Summary *get_summary_weighted (const int *list, const DATASET *dset,
			       int var, gretlopt opt, PRN *prn,
			       int *err);

Summary *get_summary_restricted (const int *list,
				 const DATASET *dset,
				 const double *rv,
				 gretlopt opt, PRN *prn,
				 int *err);

int list_summary (const int *list, int wgtvar,
		  const DATASET *dset,
		  gretlopt opt, PRN *prn);

int summary_statistics_by (const int *list, DATASET *dset,
			   gretlopt opt, PRN *prn);

void print_summary (const Summary *summ,
		    const DATASET *dset,
		    PRN *prn);

void print_summary_single (const Summary *s,
			   int digits, int places,
			   const DATASET *dset,
			   PRN *prn);

int summary_has_missing_values (const Summary *summ);

void free_summary (Summary *summ);

VMatrix *corrlist (int ci, int *list, const DATASET *dset,
		   gretlopt opt, int *err);

VMatrix *vmatrix_new (void);

void free_vmatrix (VMatrix *vmat);

int gretl_corrmx (int *list, const DATASET *dset,
		  gretlopt opt, PRN *prn);

int satterthwaite_df (double v1, int n1,
		      double v2, int n2);

int means_test (const int *list, const DATASET *dset,
		gretlopt opt, PRN *prn);

int vars_test (const int *list, const DATASET *dset,
	       gretlopt opt, PRN *prn);

void print_corrmat (VMatrix *corr, const DATASET *dset, PRN *prn);

double doornik_chisq (double skew, double xkurt, int n);

int multivariate_normality_test (const gretl_matrix *E,
				 const gretl_matrix *Sigma,
				 gretlopt opt, PRN *prn);

int mahalanobis_distance (const int *list, DATASET *dset,
			  gretlopt opt, PRN *prn);

MahalDist *get_mahal_distances (const int *list, DATASET *dset,
				gretlopt opt, PRN *prn,
				int *err);

void free_mahal_dist (MahalDist *md);

const double *mahal_dist_get_distances (const MahalDist *md);

int mahal_dist_get_n (const MahalDist *md);

const int *mahal_dist_get_varlist(const MahalDist *md);

gretl_matrix *distance (const gretl_matrix *X,
			const gretl_matrix *Y,
			const char *type, int *err);

double gretl_gini (int t1, int t2, const double *x);

int gini (int varno, DATASET *dset, gretlopt opt, PRN *prn);

int shapiro_wilk (const double *x, int t1, int t2, double *W, double *pval);

int gretl_normality_test (int varno, const DATASET *dset,
			  gretlopt opt, PRN *prn);

gretl_matrix *gretl_normtest_matrix (const double *y,
				     int t1, int t2,
				     gretlopt opt,
				     int *err);

gretl_matrix *acf_matrix (const double *x, int order,
			  const DATASET *dset, int n,
			  int *err);

gretl_matrix *xcf_vec (const double *x, const double *y,
		       int p, const DATASET *dset,
		       int n, int *err);

double ljung_box (int m, int t1, int t2, const double *y, int *err);

#endif /* DESCRIBE_H */
