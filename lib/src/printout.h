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

#define PAGELINES 21

enum prn_codes {
    GRETL_PRINT_STDOUT,
    GRETL_PRINT_STDERR,
    GRETL_PRINT_FILE,
    GRETL_PRINT_BUFFER,
    GRETL_PRINT_NULL
};

enum gretl_print_formats {
    GRETL_PRINT_FORMAT_PLAIN,
    GRETL_PRINT_FORMAT_TEX,
    GRETL_PRINT_FORMAT_TEX_DOC,
    GRETL_PRINT_FORMAT_RTF,
    GRETL_PRINT_FORMAT_FIXED
};

#define plain_format(p) (p->format == GRETL_PRINT_FORMAT_PLAIN)
#define rtf_format(p)   (p->format == GRETL_PRINT_FORMAT_RTF)
#define tex_format(p)   (p->format == GRETL_PRINT_FORMAT_TEX || \
                         p->format == GRETL_PRINT_FORMAT_TEX_DOC)
#define doc_format(p)   (p->format == GRETL_PRINT_FORMAT_TEX_DOC)
#define is_tex(f)       (f == GRETL_PRINT_FORMAT_TEX || \
                         f == GRETL_PRINT_FORMAT_TEX_DOC)
#define is_rtf(f)       (f == GRETL_PRINT_FORMAT_RTF)

/* functions follow */
 
void session_time (FILE *fp);

void logo (void);

void lib_logo (void);

void gui_script_logo (PRN *prn);

void gui_logo (FILE *fp);

void text_print_model_confints (const CONFINT *cf,
				const DATAINFO *pdinfo, 
				PRN *prn);

void print_freq (const FREQDIST *freq, PRN *prn);

void printcorr (const CorrMat *corrmat, 
		const DATAINFO *pdinfo, 
		PRN *prn);

void print_smpl (const DATAINFO *pdinfo, 
		 int fulln, PRN *prn); 

int outcovmx (MODEL *pmod, 
	      const DATAINFO *pdinfo, 
	      PRN *prn);

void print_obs_marker (int t, const DATAINFO *pdinfo, PRN *prn);

int printdata (int *list, const double **Z, const DATAINFO *pdinfo, 
	       gretlopt oflag, PRN *prn);

int text_print_fit_resid (const FITRESID *fr, 
			  const DATAINFO *pdinfo, 
			  PRN *prn);

int print_fit_resid (const MODEL *pmod, 
		     double ***pZ, DATAINFO *pdinfo, 
		     PRN *prn);

int text_print_fcast_with_errs (const FITRESID *fr, 
				double ***pZ, DATAINFO *pdinfo, PRN *prn,
				int plot);

void text_print_matrix (const double *rr, const int *list, 
			MODEL *pmod, const DATAINFO *pdinfo, 
			PRN *prn);

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

void bufspace (int n, PRN *prn);

void gretl_printxn (double x, int n, PRN *prn);

int do_printf (const char *line, double ***pZ, 
	       DATAINFO *pdinfo, MODEL *pmod,
	       PRN *prn);

int generate_obs_markers (double ***pZ, DATAINFO *pdinfo, char *s);

int in_usa (void);

char *bufgets (char *s, size_t size, const char *buf);

void scroll_pause (void);

int scroll_pause_or_quit (void);
