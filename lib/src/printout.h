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

#define GRETL_MP_DIGITS 12
#define NAMETRUNC 18

void session_time (PRN *prn);

void logo (int quiet);

void lib_logo (void);

void gui_script_logo (PRN *prn);

void gui_logo (PRN *prn);

void print_freq (const FreqDist *freq, int varno, const DATASET *dset, 
		 PRN *prn);

void print_freq_test (const FreqDist *freq, PRN *prn);

void print_xtab (const Xtab *tab, const DATASET *dset,
		 gretlopt opt, PRN *prn);

void print_smpl (const DATASET *dset, int fulln,
		 gretlopt opt, PRN *prn); 

void
print_contemp_covariance_matrix (const gretl_matrix *m, 
				 double ldet, PRN *prn);

int outcovmx (MODEL *pmod, const DATASET *dset, PRN *prn);

int max_obs_marker_length (const DATASET *dset);

void print_obs_marker (int t, const DATASET *dset, int len, PRN *prn);

char *maybe_trim_varname (char *targ, const char *src);

int max_namelen_in_list (const int *list, const DATASET *dset);

void list_series (const DATASET *dset, gretlopt opt, PRN *prn);

void maybe_list_series (const DATASET *dset, PRN *prn);

int printdata (const int *list, const char *ostr,
	       DATASET *dset, gretlopt opt,
	       PRN *prn);

int print_data_in_columns (const int *list, const int *obsvec, 
			   const DATASET *dset, gretlopt opt,
			   PRN *prn);

int print_series_with_format (const int *list,
			      const DATASET *dset, 
			      char fmt, int digits, PRN *prn);

int text_print_fit_resid (const FITRESID *fr, 
			  const DATASET *dset, 
			  PRN *prn);

int print_fit_resid (const MODEL *pmod, 
		     const DATASET *dset, 
		     PRN *prn);

int text_print_forecast (const FITRESID *fr, DATASET *dset, 
			 gretlopt opt, PRN *prn);

int print_fcast_stats_matrix (const gretl_matrix *m, int T,
			      gretlopt opt, PRN *prn);

void print_iter_info (int iter, double crit, int type, int k, 
		      const double *b, const double *g, 
		      double sl, PRN *prn);

void text_print_vmatrix (VMatrix *vmat, PRN *prn);

void gretl_sprint_fullwidth_double (double x, int digits, char *targ, PRN *prn);

void gretl_print_fullwidth_double (double x, int digits, PRN *prn);

int set_gretl_digits (int d);

int get_gretl_digits (void);

void gretl_print_value (double x, PRN *prn);

char *gretl_fix_exponent (char *s);

void bufspace (int n, PRN *prn);

void dashline (int n, PRN *prn);

void print_centered (const char *s, int width, PRN *prn);

void gretl_printxn (double x, int n, PRN *prn);

int in_usa (void);

char *bufgets (char *s, size_t size, const char *buf);

char *safe_bufgets (char **pdest, size_t *psize, const char *buf);

size_t bufgets_peek_line_length (const char *buf);

int bufgets_init (const char *buf);

int query_bufgets_init (const char *buf);

void bufgets_finalize (const char *buf);

void buf_rewind (const char *buf);

int bufseek (const char *buf, long int offset);

int buf_back_lines (const char *buf, int n);

long buftell (const char *buf);

void bufgets_cleanup (void);
