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
 
void session_time (void);

void logo (void);

void gui_logo (void);

void printmodel (const MODEL *pmod, 
		 const DATAINFO *pdinfo, 
		 print_t *prn);

void print_model_confints (const MODEL *pmod, 
			   const DATAINFO *pdinfo, 
			   print_t *prn);

void printfreq (FREQDIST *freq, 
		print_t *prn);

void printcorr (const CORRMAT *corrmat, 
		const DATAINFO *pdinfo, 
		print_t *prn);

void print_smpl (const DATAINFO *pdinfo, 
		 int fulln, print_t *prn); 

int outcovmx (MODEL *pmod, 
	      const DATAINFO *pdinfo, 
	      const int batch, print_t *prn);

void print_white_vcv (const MODEL *pmod, print_t *prn);

void takenotes (const int batch, const int runit);

int takenotes_quit (const int batch, const int runit);

void space (int n, print_t *prn);

int printdata (int *list, 
	       double **pZ, const DATAINFO *pdinfo, 
	       int batch, int byobs, print_t *prn);

int print_fit_resid (const MODEL *pmod, 
		     double **pZ, DATAINFO *pdinfo, 
		     print_t *prn);

void printxx (const double xx, char *str, const int ci);

void gretl_print_destroy (print_t *prn);

print_t *gretl_print_new (int prncode, const char *fname);

