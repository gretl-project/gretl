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
    GPT_MONO           = 1 << 19, /* monochrome output */
    GPT_GRID_Y         = 1 << 20, /* display horizontal grid lines */
    GPT_GRID_X         = 1 << 21, /* display vertical grid lines */
    GPT_POLAR          = 1 << 22, /* plot is in polar mode */
    GPT_XL             = 1 << 23, /* large */
    GPT_XXL            = 1 << 24, /* extra-large */
    GPT_XW             = 1 << 25, /* extra-wide */
    GPT_TIMEFMT        = 1 << 26, /* using gnuplot "timefmt" */
    GPT_ICON           = 1 << 27, /* saving plot "as icon" */
    GPT_STEPS          = 1 << 28, /* force steps for plot */
    GPT_LOGY           = 1 << 29  /* log y axis */
} GptFlags;

/* an extra "command" for use with GUI callback */
#define GP_ASYNC (NC+1)

typedef guint32 gretlRGB;

typedef struct GPT_SPEC_ GPT_SPEC;

#define N_GP_LINETYPES 8

#define GP_WIDTH      640
#define GP_HEIGHT     480
#define GP_LB_WIDTH   680
#define GP_LB_HEIGHT  400
#define GP_XL_WIDTH   680
#define GP_XL_HEIGHT  510
#define GP_XXL_WIDTH  680
#define GP_XXL_HEIGHT 680
#define GP_XW_WIDTH   800
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
    PLOT_FREQ_DISCRETE,
    PLOT_PERIODOGRAM,
    PLOT_CORRELOGRAM,
    PLOT_CUSUM,
    PLOT_MULTI_BASIC,
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
    PLOT_BAR,
    PLOT_STACKED_BAR,
    PLOT_3D,
    PLOT_BAND,
    PLOT_HEATMAP,
    PLOT_GEOMAP,
    PLOT_GRIDPLOT,
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
    PLOT_FIT_LINLOG,
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
    GP_TERM_HTM,
    GP_TERM_PLT,
    GP_TERM_VAR
} TermType;

#define maybe_big_multiplot(c) (c == PLOT_MULTI_IRF || \
				c == PLOT_MULTI_BASIC || \
				c == PLOT_PANEL)

#define frequency_plot_code(c) (c == PLOT_FREQ_SIMPLE || \
				c == PLOT_FREQ_NORMAL || \
				c == PLOT_FREQ_GAMMA || \
				c == PLOT_FREQ_DISCRETE)

#define set_png_output(p) (p->flags |= GPT_PNG_OUTPUT)
#define get_png_output(p) (p->flags & GPT_PNG_OUTPUT)
#define unset_png_output(p) (p->flags &= ~GPT_PNG_OUTPUT)

const char *gretl_gnuplot_term_line (TermType ttype,
				     PlotType ptype,
				     GptFlags flags,
				     const char *font);

const char *get_png_line_for_plotspec (const GPT_SPEC *spec);

char *gretl_png_font_string (void);

const char *gp_justification_string (int j);

int split_graph_fontspec (const char *s, char *name, int *psz);

double gnuplot_version (void);

char *gnuplot_version_string (void);

void gnuplot_missval_string (FILE *fp);

void write_gp_dataval (double x, FILE *fp, int final);

FILE *open_plot_input_file (PlotType ptype, GptFlags flags, int *err);

FILE *open_3d_plot_input_file (int *iact);

int finalize_plot_input_file (FILE *fp);

int finalize_3d_plot_input_file (FILE *fp);

int gnuplot_graph_wanted (PlotType ptype, gretlopt opt);

void gnuplot_cleanup (void);

int specified_gp_output_format (void);

int write_plot_output_line (const char *path, FILE *fp);

int write_plot_type_string (PlotType ptype, GptFlags flags, FILE *fp);

void write_plot_line_styles (int ptype, FILE *fp);

int write_plot_bounding_box_request (FILE *fp);

void set_effective_plot_ci (int ci);

void set_special_plot_size (float width, float height);

void set_special_font_size (int fsize);

int set_plotstyle (const char *style);

const char *get_plotstyle (void);

PlotType plot_type_from_string (const char *str);

void plot_get_scaled_dimensions (int *width, int *height, double scale);

int graph_written_to_file (void);

int plot_output_to_buffer (void);

int graph_displayed (void);

void reset_plot_count (void);

int matrix_plot (gretl_matrix *m, const int *list, const char *literal,
		 gretlopt opt);

int gnuplot (const int *plotlist, const char *literal,
	     const DATASET *dset, gretlopt opt);

int multi_plots (const int *list, const DATASET *dset,
		 gretlopt opt);

int matrix_multi_plots (const gretl_matrix *m, const int *list,
			const DATASET *dset, gretlopt opt);

int gnuplot_3d (int *list, const char *literal,
		DATASET *dset, gretlopt *opt);

int plot_freq (FreqDist *freq, DistCode dist, gretlopt opt);

int plot_corrmat (VMatrix *corr, gretlopt opt);

int garch_resid_plot (const MODEL *pmod, const DATASET *dset);

int rmplot (const int *list, DATASET *dset,
	    gretlopt opt, PRN *prn);

int hurstplot (const int *list, DATASET *dset, gretlopt opt,
	       PRN *prn);

int qq_plot (const int *list, const DATASET *dset, gretlopt opt);

int kd_plot (const int *list, const DATASET *dset, gretlopt opt);

int hf_plot (const int *list, const char *literal,
	     const DATASET *dset, gretlopt opt);

int correlogram_plot (const char *vname,
		      const double *acf,
		      const double *pacf,
		      const gretl_matrix *PM,
		      int m, double pm,
		      gretlopt opt);

int periodogram_plot (const char *vname,
		      int T, int L, const double *x,
		      gretlopt opt);

int arma_spectrum_plot (MODEL *pmod, const DATASET *dset,
			gretlopt opt);

int theil_forecast_plot (const int *plotlist, const DATASET *dset,
			 gretlopt opt);

int gretl_panel_ts_plot (int vnum, DATASET *dset, gretlopt opt);

int panel_means_XY_scatter (const int *list, const char *literal,
                            const DATASET *dset, gretlopt opt);

int cli_panel_plot (const int *list, const char *literal,
		    const DATASET *dset, gretlopt opt);

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

int gnuplot_process_input (const char *literal, gretlopt opt, PRN *prn);

int print_gnuplot_literal_lines (const char *s, int ci,
				 gretlopt opt, FILE *fp);

int is_auto_fit_string (const char *s);

void set_graph_color_from_string (int i, const char *cstr);

void graph_palette_reset (int i);

void print_rgb_hash (char *s, gretlRGB color);

gretlRGB gretl_rgb_get (const char *s);

void print_palette_string (char *s);

gretlRGB get_graph_color (int i);

gretlRGB get_boxcolor (void);

gretlRGB get_shadecolor (void);

void set_boxcolor (gretlRGB color);

void set_shadecolor (gretlRGB color);

gretlRGB numeric_color_from_string (const char *s, int *err);

int gnuplot_test_command (const char *cmd);

void set_default_png_scale (double s);

double get_default_png_scale (void);

int gnuplot_has_wxt (void);

int write_map_gp_file (void *ptr);

int transcribe_geoplot_file (const char *src,
			     const char *dest,
			     const char *datname);

int write_tdisagg_plot (const gretl_matrix *YY, int mult,
			const char *title, DATASET *dset);

#endif /* GRAPHING_H */
