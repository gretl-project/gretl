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

/* functions follow */
 
void session_time (FILE *fp);

void logo (void);

void gui_logo (FILE *fp);

void printmodel (const MODEL *pmod, 
		 const DATAINFO *pdinfo, 
		 PRN *prn);

void print_model_confints (const MODEL *pmod, 
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
	      const int pause, PRN *prn);

void print_white_vcv (const MODEL *pmod, PRN *prn);

int page_break (const int n, int *lineno, const int quit_option);

int printdata (LIST list, 
	       double **pZ, const DATAINFO *pdinfo, 
	       int pause, int byobs, PRN *prn);

int print_fit_resid (const MODEL *pmod, 
		     double **pZ, DATAINFO *pdinfo, 
		     PRN *prn);

void printxx (const double xx, char *str, const int ci);

void gretl_print_destroy (PRN *prn);

PRN *gretl_print_new (int prncode, const char *fname);

