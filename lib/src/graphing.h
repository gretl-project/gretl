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

typedef enum {
    GPT_IMPULSES       = 1 << 0,  /* use impulses for plotting */
    GPT_LINES          = 1 << 1,  /* force use of lines for plotting */
    GPT_RESIDS         = 1 << 2,  /* doing residual plot */
    GPT_FA             = 1 << 3,  /* doing fitted/actual plot */
    GPT_DUMMY          = 1 << 4,  /* using a dummy for separation */
    GPT_XYZ            = 1 << 5,  /* X-Y, controlling for Z */
    GPT_FIT_OMIT       = 1 << 6,  /* user said don't draw fitted line on graph */
    GPT_DATA_STYLE     = 1 << 7,  /* data style is set by user */
    GPT_IDX            = 1 << 8,  /* plot against time or obs index */
    GPT_TS             = 1 << 9,  /* doing time series plot */
    GPT_Y2AXIS         = 1 << 10, /* plot has second y-axis */
    GPT_AUTO_FIT       = 1 << 11, /* automatic (OLS) fitted line was added */
    GPT_FIT_HIDDEN     = 1 << 12, /* autofit line calculated, but suppressed */
    GPT_PNG_OUTPUT     = 1 << 13, /* output is to PNG file */
    GPT_PRINT_MARKERS  = 1 << 14, /* print observation markers */
    GPT_LETTERBOX      = 1 << 15, /* special format for time series graphs */
    GPT_PARAMETRIC     = 1 << 16, /* gnuplot should be in parametric mode */
    GPT_XZEROAXIS      = 1 << 17, /* show x = 0 line */
    GPT_YZEROAXIS      = 1 << 18, /* show y = 0 line */
    GPT_FILL_SWITCH    = 1 << 19, /* switching from errorbars to fill */
    GPT_ERR_SWITCH     = 1 << 20, /* switching from fill to errorbars */
    GPT_MONO           = 1 << 21, /* monochrome output */
    GPT_GRID_Y         = 1 << 22, /* display horizontal grid lines */
    GPT_GRID_X         = 1 << 23, /* display vertical grid lines */
    GPT_POLAR          = 1 << 24, /* plot is in polar mode */
    GPT_XL             = 1 << 25, /* large */
    GPT_XXL            = 1 << 26, /* extra-large */
    GPT_TIMEFMT        = 1 << 27  /* using gnuplot "timefmt" */
} GptFlags; 

typedef struct gretlRGB_ gretlRGB;

struct gretlRGB_ {
    unsigned char r;
    unsigned char g;
    unsigned char b;
};

typedef struct GPT_SPEC_ GPT_SPEC;

#define MAXTITLE 128
#define N_GP_COLORS 8
#define BOXCOLOR (N_GP_COLORS - 2)
#define SHADECOLOR (N_GP_COLORS - 1)

#define GP_WIDTH      640
#define GP_HEIGHT     480
#define GP_LB_WIDTH   680
#define GP_LB_HEIGHT  400
#define GP_XL_WIDTH   680
#define GP_XL_HEIGHT  510
#define GP_XXL_WIDTH  680
#define GP_XXL_HEIGHT 680
#define GP_SQ_SIZE    480

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
    PLOT_ROOTS,
    PLOT_ELLIPSE,
    PLOT_MULTI_IRF,
    PLOT_PANEL,
    PLOT_BI_GRAPH,
    PLOT_MANY_TS,
    PLOT_RQ_TAU,
    PLOT_FACTORIZED,
    PLOT_BOXPLOTS,
    PLOT_CURVE,
    PLOT_QQ,
    PLOT_USER,
    PLOT_XCORRELOGRAM,
    PLOT_TYPE_MAX
} PlotType;

typedef enum {
    PLOT_FIT_NONE,
    PLOT_FIT_OLS,
    PLOT_FIT_QUADRATIC,
    PLOT_FIT_CUBIC,
    PLOT_FIT_INVERSE,
    PLOT_FIT_LOESS,
    PLOT_FIT_LOGLIN,
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
    GP_PNG_GD1,   /* libgd driver, no truecolor option */
    GP_PNG_GD2,   /* libgd with truecolor option */
    GP_PNG_CAIRO  /* newer cairo-based driver */
} PNGTerm;

typedef enum {
    GP_PDF_NONE,
    GP_PDF_PDFLIB,
    GP_PDF_CAIRO
} PDFTerm;

typedef enum {
    GP_EPS_NONE,
    GP_EPS_PS,
    GP_EPS_CAIRO
} EPSTerm;

#define maybe_big_multiplot(c) (c == PLOT_MULTI_IRF || \
				c == PLOT_MULTI_SCATTER || \
				c == PLOT_PANEL)

#define frequency_plot_code(c) (c == PLOT_FREQ_SIMPLE || \
				c == PLOT_FREQ_NORMAL || \
				c == PLOT_FREQ_GAMMA)

#define set_png_output(p) (p->flags |= GPT_PNG_OUTPUT)
#define get_png_output(p) (p->flags & GPT_PNG_OUTPUT)
#define unset_png_output(p) (p->flags &= ~GPT_PNG_OUTPUT) 
    
const char *get_gretl_png_term_line (PlotType ptype, GptFlags flags);

const char *get_png_line_for_plotspec (const GPT_SPEC *spec);

const char *get_gretl_emf_term_line (PlotType ptype, int color);

const char *get_gretl_pdf_term_line (PlotType ptype, GptFlags flags);

const char *get_gretl_eps_term_line (PlotType ptype, GptFlags flags);

const char *gp_justification_string (int j);

int split_graph_fontspec (const char *s, char *name, int *psz);

void gnuplot_missval_string (FILE *fp);

FILE *open_plot_input_file (PlotType ptype, int *err);

int finalize_plot_input_file (FILE *fp);

int gnuplot_graph_wanted (PlotType ptype, gretlopt opt);

void gnuplot_cleanup (void);

int specified_gp_output_format (void);

void write_plot_output_line (const char *path, FILE *fp);

int write_plot_type_string (PlotType ptype, GptFlags flags, FILE *fp);

void write_plot_line_styles (int ptype, FILE *fp);

void write_plot_bounding_box_request (FILE *fp);

PlotType plot_type_from_string (const char *str);

void plot_get_scaled_dimensions (int *width, int *height, double scale);

int graph_written_to_file (void);

void reset_plot_count (void);

int matrix_plot (gretl_matrix *m, const int *list, const char *literal, 
		 gretlopt opt);

int gnuplot (const int *plotlist, const char *literal,
	     const DATASET *dset, gretlopt opt);

int multi_scatters (const int *list, const DATASET *dset, 
		    gretlopt opt);

int matrix_scatters (const gretl_matrix *m, const int *list, 
		     const DATASET *dset, gretlopt opt);

int gnuplot_3d (int *list, const char *literal,
		DATASET *dset, gretlopt opt);

int plot_freq (FreqDist *freq, DistCode dist);

int garch_resid_plot (const MODEL *pmod, const DATASET *dset); 

int rmplot (const int *list, DATASET *dset, 
	    gretlopt opt, PRN *prn);

int hurstplot (const int *list, DATASET *dset, gretlopt opt,
	       PRN *prn);

int qq_plot (const int *list, const DATASET *dset, gretlopt opt);

int correlogram_plot (const char *vname,
		      const double *acf, 
		      const double *pacf,
		      int m, double pm, 
		      gretlopt opt);

int periodogram_plot (const char *vname,
		      int T, int L, const double *x,
		      gretlopt opt);

int theil_forecast_plot (const int *plotlist, const DATASET *dset, 
			 gretlopt opt);

int gretl_panel_ts_plot (int vnum, DATASET *dset, gretlopt opt);

int plot_fcast_errs (const FITRESID *fr, const double *maxerr,
		     const DATASET *dset, gretlopt opt);

int plot_simple_fcast_bands (const MODEL *pmod, 
			     const FITRESID *fr, 
			     const DATASET *dset, 
			     gretlopt opt);

int plot_tau_sequence (const MODEL *pmod, const DATASET *dset,
		       int k);

int 
gretl_VAR_plot_impulse_response (GRETL_VAR *var,
				 int targ, int shock, 
				 int periods, double alpha,
				 const DATASET *dset,
				 gretlopt opt);

int gretl_VAR_plot_FEVD (GRETL_VAR *var, int targ, int periods, 
			 const DATASET *dset, gretlopt opt);

int 
gretl_VAR_plot_multiple_irf (GRETL_VAR *var, 
			     int periods, double alpha,
			     const DATASET *dset,
			     gretlopt opt);

int gretl_VECM_combined_EC_plot (GRETL_VAR *var, 
				 const DATASET *dset);

int gretl_system_residual_plot (void *p, int ci, int eqn, const DATASET *dset);

int gretl_system_residual_mplot (void *p, int ci, const DATASET *dset); 

int gretl_VAR_roots_plot (GRETL_VAR *var);

int confidence_ellipse_plot (gretl_matrix *V, double *b, 
			     double tcrit, double Fcrit, double alpha,
			     const char *iname, const char *jname);

int xy_plot_with_control (const int *list, const char *literal,
			  const DATASET *dset, gretlopt opt);

int gnuplot_process_file (gretlopt opt, PRN *prn);

void print_gnuplot_literal_lines (const char *s, FILE *fp);

int is_auto_fit_string (const char *s);

int gnuplot_has_ttf (int reset);

int gnuplot_pdf_terminal (void);

int gnuplot_eps_terminal (void);

int gnuplot_png_terminal (void);

int gnuplot_has_cp950 (void);

void set_graph_palette (int i, gretlRGB color);

void set_graph_palette_from_string (int i, const char *cstr);

void graph_palette_reset (int i);

void print_rgb_hash (char *s, const gretlRGB *color);

void gretl_rgb_get (gretlRGB *color, const char *s);

void print_palette_string (char *s);

const gretlRGB *get_graph_color (int i);

int gnuplot_test_command (const char *cmd);

void gnuplot_png_set_use_aa (int s);

void gnuplot_png_set_default_scale (double s);

void date_from_gnuplot_time (char *targ, size_t tsize, 
			     const char *fmt, double x);

double gnuplot_time_from_date (const char *s, 
			       const char *fmt);

double gnuplot_version (void);

#ifndef WIN32
int gnuplot_has_wxt (void);
#endif

#endif /* GRAPHING_H */

