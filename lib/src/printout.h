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

/* printout.h for gretl */

#define GRETL_DIGITS 6
#define GRETL_MP_DIGITS 12

#define PAGELINES 21

/* functions follow */
 
void session_time (PRN *prn);

void logo (void);

void lib_logo (void);

void gui_script_logo (PRN *prn);

void gui_logo (PRN *prn);

void text_print_model_confints (const CoeffIntervals *cf, PRN *prn);

void print_freq (const FreqDist *freq, PRN *prn);

void print_xtab (const Xtab *tab, gretlopt opt, PRN *prn);

void print_smpl (const DATAINFO *pdinfo, 
		 int fulln, PRN *prn); 

void
print_contemp_covariance_matrix (const gretl_matrix *m, 
				 double ldet, PRN *prn);

int outcovmx (MODEL *pmod, const DATAINFO *pdinfo, PRN *prn);

void obs_marker_init (const DATAINFO *pdinfo);

void print_obs_marker (int t, const DATAINFO *pdinfo, PRN *prn);

void varlist (const DATAINFO *pdinfo, PRN *prn);

void maybe_list_vars (const DATAINFO *pdinfo, PRN *prn);

int get_printdata_blocks (void);

int printdata (const int *list, const char *mstr,
	       const double **Z, const DATAINFO *pdinfo, 
	       gretlopt opt, PRN *prn);

int print_data_sorted (const int *list, const int *obsvec, 
		       const double **Z, const DATAINFO *pdinfo, 
		       PRN *prn);

int text_print_fit_resid (const FITRESID *fr, 
			  const DATAINFO *pdinfo, 
			  PRN *prn);

int print_fit_resid (const MODEL *pmod, 
		     const double **Z, const DATAINFO *pdinfo, 
		     PRN *prn);

int text_print_forecast (const FITRESID *fr, DATAINFO *pdinfo, 
			 gretlopt opt, PRN *prn);

void print_iter_info (int iter, double crit, int type, int k, 
		      const double *b, const double *g, 
		      double sl, PRN *prn);

void text_print_vmatrix (VMatrix *vmat, PRN *prn);

void gretl_print_fullwidth_double (double x, int digits, PRN *prn);

void gretl_print_value (double x, PRN *prn);

char *gretl_fix_exponent (char *s);

void bufspace (int n, PRN *prn);

void print_centered (const char *s, int width, PRN *prn);

void gretl_printxn (double x, int n, PRN *prn);

int in_usa (void);

char *bufgets (char *s, size_t size, const char *buf);

void bufgets_init (const char *buf);

void bufgets_finalize (const char *buf);

void scroll_pause (void);

int scroll_pause_or_quit (void);
