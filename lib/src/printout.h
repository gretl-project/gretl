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

/* printout.h for gretl */

#define GRETL_DIGITS 6
#define GRETL_MP_DIGITS 12

/* functions follow */
 
void session_time (FILE *fp);

void logo (void);

void lib_logo (void);

void gui_logo (FILE *fp);

void text_print_model_confints (const CONFINT *cf,
				const DATAINFO *pdinfo, 
				PRN *prn);

void printfreq (FREQDIST *freq, 
		PRN *prn);

void printcorr (const CORRMAT *corrmat, 
		const DATAINFO *pdinfo, 
		PRN *prn);

void print_smpl (const DATAINFO *pdinfo, 
		 int fulln, PRN *prn); 

int outcovmx (MODEL *pmod, 
	      const DATAINFO *pdinfo, 
	      int pause, PRN *prn);

int page_break (int n, int *lineno, int quit_option);

void print_obs_marker (int t, const DATAINFO *pdinfo, PRN *prn);

int printdata (LIST list, 
	       double ***pZ, const DATAINFO *pdinfo, 
	       int pause, unsigned long oflag, PRN *prn);

int text_print_fit_resid (const FITRESID *fr, 
			  const DATAINFO *pdinfo, 
			  PRN *prn);

int print_fit_resid (const MODEL *pmod, 
		     double ***pZ, DATAINFO *pdinfo, 
		     PRN *prn);

int text_print_fcast_with_errs (const FITRESID *fr, 
				double ***pZ, DATAINFO *pdinfo, PRN *prn,
				PATHS *ppaths, int plot);

void text_print_matrix (const double *rr, const int *list, 
			MODEL *pmod, const DATAINFO *pdinfo, 
			int pause, PRN *prn);

void gretl_print_fullwidth_double (double x, int digits, PRN *prn);

void gretl_print_value (double x, PRN *prn);

char *gretl_fix_exponent (char *s);

void gretl_print_destroy (PRN *prn);

PRN *gretl_print_new (int prncode, const char *fname);

void gretl_print_attach_buffer (PRN *prn, char *buf);

void gretl_print_attach_file (PRN *prn, FILE *fp);

int pprintf (PRN *prn, const char *template, ...);

int pputs (PRN *prn, const char *s);

int pputc (PRN *prn, int c);

int do_printf (const char *line, double ***pZ, 
	       DATAINFO *pdinfo, MODEL *pmod,
	       PRN *prn);

int in_usa (void);

char *bufgets (char *s, size_t size, const char *buf);
