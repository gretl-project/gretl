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

/* graphing.h for gretl */

#ifndef GRAPHING_H
#define GRAPHING_H

#include <stdio.h>

#define GRAPH_NO_DATA -999

typedef enum {
    GPT_IMPULSES       = 1 << 0,  /* use impulses for plotting */
    GPT_LINES          = 1 << 1,  /* force use of lines for plotting */
    GPT_RESIDS         = 1 << 2,  /* doing residual plot */
    GPT_FA             = 1 << 3,  /* doing fitted/actual plot */
    GPT_DUMMY          = 1 << 4,  /* using a dummy for separation */
    GPT_XYZ            = 1 << 5,  /* X-Y, controlling for Z */
    GPT_BATCH          = 1 << 6,  /* working in batch mode */
    GPT_GUI            = 1 << 7,  /* called from GUI context */
    GPT_FIT_OMIT       = 1 << 8,  /* User said don't draw fitted line on graph */
    GPT_DATA_STYLE     = 1 << 9,  /* data style is set by user */
    GPT_FILE           = 1 << 10, /* send output to named file */
    GPT_IDX            = 1 << 11, /* plot against time or obs index */
    GPT_TS             = 1 << 12, /* doing time series plot */
    GPT_Y2AXIS         = 1 << 13, /* plot has second y-axis */
    GPT_AUTO_FIT       = 1 << 14, /* automatic (OLS) fitted line was added */
    GPT_FIT_HIDDEN     = 1 << 15, /* autofit line calculated, but suppressed */
    GPT_PNG_OUTPUT     = 1 << 16, /* output is to PNG file */
    GPT_PRINT_MARKERS  = 1 << 17, /* print observation markers */
    GPT_LETTERBOX      = 1 << 18, /* special format for time series graphs */
    GPT_PARAMETRIC     = 1 << 19, /* gnuplot should be in parametric mode */
    GPT_XZEROAXIS      = 1 << 20, /* show x = 0 line */
    GPT_YZEROAXIS      = 1 << 21, /* show y = 0 line */
    GPT_FILL_SWITCH    = 1 << 22, /* switching from errorbars to fill */
    GPT_ERR_SWITCH     = 1 << 23, /* switching from fill to errorbars */
    GPT_MONO           = 1 << 24  /* monochrome output */
} GptFlags; 

typedef struct gretlRGB_ gretlRGB;

struct gretlRGB_ {
    unsigned char r;
    unsigned char g;
    unsigned char b;
};

#define MAXTITLE 128
#define N_GP_COLORS 7
#define BOXCOLOR (N_GP_COLORS - 1)

typedef enum {
    PLOT_REGULAR = 0,
    PLOT_H_TEST,
    PLOT_PROB_DIST,
    PLOT_FORECAST,
    PLOT_GARCH,
    PLOT_FREQ_SIMPLE,
    PLOT_FREQ_NORMAL,
    PLOT_FREQ_GAMMA,
    PLOT_PERIODOGRAM,
    PLOT_CORRELOGRAM,
    PLOT_CUSUM,
    PLOT_MULTI_SCATTER,
    PLOT_TRI_GRAPH,
    PLOT_RANGE_MEAN,
    PLOT_HURST,
    PLOT_LEVERAGE,
    PLOT_IRFBOOT,
    PLOT_KERNEL,
    PLOT_VAR_ROOTS,
    PLOT_ELLIPSE,
    PLOT_MULTI_IRF,
    PLOT_PANEL,
    PLOT_BI_GRAPH,
    PLOT_MANY_TS,
    PLOT_RQ_TAU,
    PLOT_BOXPLOTS,
    PLOT_CURVE,
    PLOT_TYPE_MAX
} PlotType;

typedef enum {
    PLOT_FIT_NONE,
    PLOT_FIT_OLS,
    PLOT_FIT_QUADRATIC,
    PLOT_FIT_INVERSE,
    PLOT_FIT_LOESS,
    PLOT_FIT_NA       /* fit option not applicable */
} FitType;

typedef enum {
    GP_TERM_NONE,
    GP_TERM_PNG,
    GP_TERM_EPS,
    GP_TERM_PDF,
    GP_TERM_FIG,
    GP_TERM_TEX,
    GP_TERM_EMF,
    GP_TERM_SVG,
    GP_TERM_PLT
} TermType;

typedef enum {
    GP_PNG_NONE,  /* no PNG terminal available */
    GP_PNG_OLD,   /* old driver from gnuplot 3.N */
    GP_PNG_GD1,   /* libgd driver, no truecolor option */
    GP_PNG_GD2,   /* libgd with truecolor option */
    GP_PNG_CAIRO  /* newer cairo-based driver */
} PNGTerm;

typedef enum {
    GP_PDF_NONE,
    GP_PDF_PDFLIB,
    GP_PDF_CAIRO
} PDFTerm;

#define frequency_plot_code(c) (c == PLOT_FREQ_SIMPLE || \
				c == PLOT_FREQ_NORMAL || \
				c == PLOT_FREQ_GAMMA)

#define set_png_output(p) (p->flags |= GPT_PNG_OUTPUT)
#define get_png_output(p) (p->flags & GPT_PNG_OUTPUT) 
    
const char *get_gretl_png_term_line (PlotType ptype, GptFlags flags);

const char *get_gretl_emf_term_line (PlotType ptype, int color);

const char *gp_justification_string (int j);

const char *gnuplot_label_front_string (void);

int split_graph_fontspec (const char *s, char *name, int *psz);

void gnuplot_missval_string (FILE *fp);

int gnuplot_init (PlotType ptype, FILE **fpp);

FILE *gnuplot_batch_init (int *err);

int specified_gp_output_format (void);

int write_plot_type_string (PlotType ptype, FILE *fp);

void write_plot_line_styles (int ptype, FILE *fp);

void print_plot_bounding_box_request (FILE *fp);

PlotType plot_type_from_string (const char *str);

int gnuplot_make_graph (void);

void reset_plot_count (void);

int matrix_plot (gretl_matrix *m, const int *list, const char *literal, 
		 gretlopt opt);

int gnuplot (const int *plotlist, const char *literal,
	     const double **Z, const DATAINFO *pdinfo, 
	     gretlopt opt);

int multi_scatters (const int *list, const double **Z,
		    const DATAINFO *pdinfo, gretlopt opt);

int gnuplot_3d (int *list, const char *literal,
		double ***pZ, DATAINFO *pdinfo, 
		gretlopt opt);

int plot_freq (FreqDist *freq, DistCode dist);

int garch_resid_plot (const MODEL *pmod, const DATAINFO *pdinfo); 

int rmplot (const int *list, const double **Z, DATAINFO *pdinfo, 
	    PRN *prn);

int hurstplot (const int *list, const double **Z, DATAINFO *pdinfo, 
	       PRN *prn);

int theil_forecast_plot (const int *plotlist, const double **Z, 
			 const DATAINFO *pdinfo, gretlopt opt);

int 
gretl_panel_ts_plot (const int *list, const double **Z, DATAINFO *pdinfo,
		     gretlopt opt);

int plot_fcast_errs (const FITRESID *fr, const double *maxerr,
		     const DATAINFO *pdinfo, gretlopt opt);

int plot_tau_sequence (const MODEL *pmod, const DATAINFO *pdinfo,
		       int k);

int 
gretl_VAR_plot_impulse_response (GRETL_VAR *var,
				 int targ, int shock, int periods,
				 double alpha,
				 const double **Z,
				 const DATAINFO *pdinfo);

int 
gretl_VAR_plot_multiple_irf (GRETL_VAR *var, int periods,
			     double alpha,
			     const double **Z,
			     const DATAINFO *pdinfo);

int gretl_system_residual_plot (void *p, int ci, const DATAINFO *pdinfo);

int gretl_system_residual_mplot (void *p, int ci, const DATAINFO *pdinfo); 

int gretl_VAR_roots_plot (GRETL_VAR *var);

int confidence_ellipse_plot (gretl_matrix *V, double *b, 
			     double tcrit, double Fcrit, double alpha,
			     const char *iname, const char *jname);

int xy_plot_with_control (const int *list, const char *literal,
			  const double **Z, const DATAINFO *pdinfo,
			  gretlopt opt);

int is_auto_fit_string (const char *s);

int gnuplot_has_ttf (int reset);

int gnuplot_pdf_terminal (void);

int gnuplot_png_terminal (void);

int gnuplot_has_rgb (void);

int gnuplot_has_style_fill (void);

int gnuplot_has_latin5 (void);

int gnuplot_has_cp1250 (void);

int gnuplot_has_cp1254 (void);

int gnuplot_has_bbox (void);

int gnuplot_has_utf8 (void);

int gnuplot_uses_datafile_missing (void);

void set_graph_palette (int i, gretlRGB color);

void set_graph_palette_from_string (int i, const char *cstr);

void graph_palette_reset (int i);

void print_rgb_hash (char *s, const gretlRGB *color);

void gretl_rgb_get (gretlRGB *color, const char *s);

void print_palette_string (char *s);

const gretlRGB *get_graph_color (int i);

int gnuplot_test_command (const char *cmd);

void gnuplot_png_set_use_aa (int s);

#ifndef WIN32
int gnuplot_has_wxt (void);
#endif

#endif /* GRAPHING_H */

