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

#ifndef GENFUNCS_H
#define GENFUNCS_H

int sort_series (const double *x, double *y, int f, 
		 const DATASET *dset);

int gretl_sort_by (const double *x, const double *y, 
		   double *z, const DATASET *dset);

int rank_series (const double *x, double *y, int f, 
		 const DATASET *dset);

gretl_matrix *rank_vector (const gretl_matrix *x, int f, int *err);

int diff_series (const double *x, double *y, int f, 
		 const DATASET *dset);

int orthdev_series (const double *x, double *y, const DATASET *dset);

int cum_series (const double *x, double *y, 
		const DATASET *dset);

int resample_series (const double *x, double *y, 
		     const DATASET *dset);

int block_resample_series (const double *x, double *y, int blocklen,
			   const DATASET *dset);

int fracdiff_series (const double *x, double *y, double d,
		     int diff, int obs, const DATASET *dset);

int boxcox_series (const double *x, double *y, double d,
		   const DATASET *dset);

int filter_series (const double *x, double *y, const DATASET *dset, 
		   gretl_matrix *A, gretl_matrix *C, double y0);

gretl_matrix *filter_matrix (gretl_matrix *X, gretl_vector *A, gretl_vector *C, 
			     double y0, int *err);

int exponential_movavg_series (const double *x, double *y, 
			       const DATASET *dset,
			       double d, int n);

int movavg_series (const double *x, double *y, const DATASET *dset,
		   int k, int center);

int seasonally_adjust_series (const double *x, double *y, 
			      DATASET *dset, int tramo,
			      int use_log);

int panel_statistic (const double *x, double *y, const DATASET *dset, 
		     int k, const double *mask);

gretl_matrix *panel_shrink (const double *x, const DATASET *dset,
			    int *err);

int hp_filter (const double *x, double *hp, const DATASET *dset,
	       double lambda, gretlopt opt);

int bkbp_filter (const double *x, double *bk, const DATASET *dset, 
		 int bkl, int bku, int k);

int butterworth_filter (const double *x, double *bw, const DATASET *dset,
			int n, double cutoff);

int poly_trend (const double *x, double *fx, const DATASET *dset, int order);

int weighted_poly_trend (const double *x, double *fx, const DATASET *dset,
			 int order, gretlopt opt, double wratio, 
			 double midfrac);

void poly_weights (double *w, int T, double wmax, 
		   double midfrac, gretlopt opt);

gretl_matrix *hp_gain (double lambda, int hipass);

gretl_matrix *butterworth_gain (int n, double cutoff, int hipass);

int dummy (DATASET *dset, int center);

int panel_dummies (DATASET *dset, gretlopt opt, PRN *prn);

int gen_unit (DATASET *dset);

int panel_unit_first_obs (int t, const DATASET *dset);

int gen_time (DATASET *dset, int tm);

int gen_wkday (DATASET *dset);

const double *gretl_plotx (const DATASET *dset, gretlopt opt);

double *get_fit_or_resid (const MODEL *pmod, DATASET *dset, 
			  ModelDataIndex idx, char *vname, 
			  char *vlabel, int *err);

int get_observation_number (const char *s, const DATASET *dset);

int get_t_from_obs_string (const char *s, const DATASET *dset);

int list_linear_combo (double *y, const int *list, 
		       const gretl_vector *b, 
		       const DATASET *dset);

double imhof (const gretl_matrix *m, double arg, int *err);

double dw_pval (const gretl_matrix *u, const gretl_matrix *X, 
		double *pDW, int *err);

gretl_matrix *multi_acf (const gretl_matrix *m, 
			 const int *list, 
			 const DATASET *dset,
			 int p, int *err);

gretl_matrix *multi_xcf (const void *px, int xtype,
			 const void *py, int ytype,
			 const DATASET *dset,
			 int p, int *err);

gretl_matrix *forecast_stats (const double *y, const double *f,
			      int t1, int t2, gretlopt opt,
			      int *err);

double gretl_round (double x);

double gretl_bessel (char type, double v, double x, int *err);

double gretl_npv (int t1, int t2, const double *x, double r, 
		  int pd, int *err);

double gretl_irr (const double *x, int n, int pd, int *err);

double logistic_cdf (double x);

gretl_matrix *matrix_chowlin (const gretl_matrix *Y, 
			      const gretl_matrix *X,
			      int f, int *err);

int list_ok_dollar_vars (DATASET *dset, PRN *prn);

int nadaraya_watson (const double *y, const double *x, double h,
		     DATASET *dset, double *m);

int gretl_loess (const double *y, const double *x, int poly_order,
		 double bandwidth, gretlopt opt, DATASET *dset, 
		 double *m);

double series_get_nobs (int t1, int t2, const double *x);

double series_sum_all (int t1, int t2, const double *x);

gretl_matrix *aggregate_by (const double *x, 
			    const double *y,
			    const int *xlist,
			    const int *ylist,
			    const char *fncall,
			    const DATASET *dset,
			    int *err);

int fill_dataset_dates_series (const DATASET *dset, double *x);

#endif /* GENFUNCS_H */
